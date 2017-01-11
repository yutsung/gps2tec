'''
This is GPS TEC Calculate Program

Author      : Wei-Han, Chen
Since       : 2016/11/25
Update notes:

'''

import os
import sys
import datetime
import time
import glob
import warnings
import multiprocessing
import shutil

import numpy as np

from gps2tec_modules import IO
from gps2tec_modules import Getdata
from gps2tec_modules import DataProcess
from gps2tec_modules import StationBias


def main():
    warnings.filterwarnings("ignore")

    input_para = IO.Input_file(sys.argv[1])
    main_pwd = os.getcwd()

    for run in range(input_para.total_run):
        # deal with time
        target_st_time = datetime.datetime(input_para.year, 1, 1) + datetime.timedelta(days=input_para.st_doy - 1)
        target_time = target_st_time + datetime.timedelta(days=run)
        target_doy = (target_time - datetime.datetime(target_time.year, 1, 1)).days + 1

        # search gps ofile
        gpsfile = Getdata.GPSourceFile(target_time.year, target_doy)
        if input_para.case_type == 'local':
            shutil.copy(input_para.crx2rnx_pwd, 'crx2rnx')
        elif input_para.case_type in 'igsrt':
            DataProcess.preprocess(target_time.year, target_doy, input_para.save_pwd, input_para.download_list_fn, input_para.crx2rnx_pwd) # create directory and copy needed file () to save_pwd
            gpsfile.download_data(input_para.case_type, input_para.download_list_fn)
        else:
            print "case_type error"

        gpsfile.decompress_file()
        gpsfile.d2o()

        # check gps ofile
        if input_para.case_type == 'igsrt':
            pass#stn_list = sorted(glob.glob('*{0:03}{1}{2:02}.{3:02}o'.format(target_doy, chr(97 + target_hh), target_mm, target_year % 100)))
        else:
            stn_list = sorted(glob.glob("*{0:03}0.{1:02}o".format(target_doy, target_time.year % 100)))
        if stn_list == []:
            print "No GPS observation data on {0}.{1:03}".format(target_time.year, target_doy)
            break

        # get GIM data
        gim = Getdata.GimData(target_time.year, target_doy)
        gim.get_gimdata()
        gim.create_mapfile()
        gim.read_station_bias()

        # get Navigation data
        navidata = Getdata.NavigationData(target_time.year, target_doy)
        navidata.get_navidata()
        satdata = navidata.read_navidata()

        # get P1P2 code
        p1p2 = Getdata.P1P2codeData(target_time.year, target_doy)
        p1p2.get_p1p2data()

        # get P1C1 code
        p1c1 = Getdata.P1C1codeData(target_time.year, target_doy)
        p1c1.get_p1c1data()

        print "Start to process GPS data..."
        interval = len(stn_list) // input_para.mp_num
        for i in range(input_para.mp_num):
            st_idx = i * interval
            ed_idx = i * interval + interval
            if i == input_para.mp_num-1 : ed_idx += len(stn_list) % input_para.mp_num
            multi_args = (stn_list[st_idx:ed_idx], input_para, target_time, target_doy, satdata, gim, i)
            p = multiprocessing.Process(target=multi_gps2tec, args=multi_args)
            p.start()


def multi_gps2tec(stn_list, input_para, target_time, target_doy, satdata, gim, subprocess_num):

    for stn_num, station in enumerate(stn_list, start=1):
        target_stn = station[0:4]
        print "Process:{4:>2} ({1:4}/{2:4}), DOY:{3:>3}, Station: {0}, ".format(target_stn, stn_num, len(stn_list), target_doy, subprocess_num),
        st_time = time.time()
        ofile = Getdata.GPSofileData(target_stn, target_time.year, target_doy)
        stnx, stny, stnz = ofile.read_ofile_xyz()
        if stnx == 0.0 or stny == 0.0 or stnz == 0.0:
            if os.path.isfile('marker.crd'):
                stnifo = open('marker.crd', 'r')
                try:
                    rdline = stnifo.readline()
                    while not (rdline[0:4] == target_stn):
                        if rdline == "":
                            raise EOFError
                        rdline = stnifo.readline()
                    stnx = float(rdline[5:19])
                    stny = float(rdline[20:34])
                    stnz = float(rdline[35:49])
                except EOFError:
                    print "Please write on the {0}'s XYZ information in marker.crd".format(target_stn)
                    continue
            else:
                print "The marker.crd is not exist!!"
                sys.exit()
        stn_lon, stn_lat, stn_alt = DataProcess.xyz2g(stnx, stny, stnz)
        ####    load gps ofile P1, P2, L1, L2
        try:
            # st = time.time()
            P1, P2, L1, L2, P1yn, dt = ofile.read_ofile_PLvalue()
            # ed = time.time()
            # print "read ofile:  {0}".format(ed - st)
        except TypeError:
            print "the ofile has no enough data!"
            continue

        # st = time.time()
        e = DataProcess.elevation(satdata, stnx, stny, stnz)
        sate_bias = DataProcess.satellite_bias(target_time.year, target_doy, P1yn)
        # ed = time.time()
        # print "calculate e, satellite bias:  {0}".format(ed - st)

        # check station bias
        stnbs_file = Getdata.StationBiasData(target_time.year, target_doy)
        if stnbs_file.file_exist():
            stn_bias = stnbs_file.read_stnbias(target_stn)
            if stn_bias == 0:
                check_stnbs = False
                write_stnbs = False
            else:
                check_stnbs = True
                write_stnbs = True
        else:
            stn_bias = 0
            check_stnbs = False
            write_stnbs = False

        # start to calculate stec, vtec, stnbias
        check_calstec = False
        GIMminTEC, GIMminTime = gim.min_gimTEC(stn_lon, stn_lat)
        while (not check_stnbs) or (not check_calstec):
            if check_calstec: check_stnbs = True
            Pcode, Lcode = DataProcess.cal_PLcode(P1, P2, L1, L2, sate_bias, stn_bias, e, input_para.elevation_angle)
            stec = DataProcess.cal_stec(Pcode, Lcode, 30)
            vtec, lon, lat, B = DataProcess.cal_vtec(stec, stnx, stny, stnz, satdata, e, input_para.elevation_angle,
                                                     input_para.h1, input_para.h2, GIMminTEC)
            check_calstec = True
            if not check_stnbs:
                if input_para.stnbias_method == 'gim':
                    stn_bias = StationBias.cal_stnbias_GIM(B, e, stn_lon, GIMminTime)
                    if type(stn_bias) == bool: break
                    continue
                elif input_para.stnbias_method == 'fitting':
                    ep = np.linspace(0, 2879, 2880) * 30 / 3600
                    stn_bias = StationBias.cal_stnbias_FITTING(stec, e, ep, satdata, stnx, stny, stnz)

        if type(stn_bias) == bool:
            print "cannot calculate station bias!"
            continue
        minTEC = DataProcess.cal_mintec(vtec, stn_lon, GIMminTime)

        # st = time.time()
        # start to output
        op_file = IO.OutputData(input_para.case_type, target_stn, target_time.year, target_doy)
        op_file.output_obs(lon, lat, vtec, stec, e)
        op_file.output_log(stn_lon, stn_lat, vtec.shape[0], minTEC, stn_bias)
        stnbs_file.write_stnbias(target_stn, stn_bias, write_stnbs, minTEC, GIMminTEC, input_para.stnbias_method)
        # ed = time.time()
        # print "write data:  {0}".format(ed - st)

        ed_time = time.time()
        print "Scucess!  Minmum TEC: {4:8.4f}, GIMminTEC: {5:8.4f}, Bias: {6:8.4f}, RunTime: {7:5.2f}".format(
            target_doy, stn_num + 1, len(stn_list), target_stn, minTEC, GIMminTEC, stn_bias, ed_time - st_time)


main()
