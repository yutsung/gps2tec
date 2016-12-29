'''
This is GPS TEC Calculate Program

Author      : Wei-Han, Chen
Since       : 2016/11/25
Update notes:

'''

import os
import sys
import datetime
import glob

from gps2tec_modules import Tools
from gps2tec_modules import Getdata
from gps2tec_modules import DataProcess


input_para = Tools.Input_file(sys.argv[1])
main_pwd = os.getcwd()

for run in range(input_para.total_run):
    # deal with time
    target_st_time = datetime.datetime(input_para.year, 1, 1) + datetime.timedelta(days=input_para.st_doy - 1)
    target_time = target_st_time + datetime.timedelta(days=run)
    target_doy = (target_time - datetime.datetime(target_time.year, 1, 1)).days + 1

    # search gps ofile
    gpsfile = Getdata.GPSourceFile(target_time.year, target_doy)
    if input_para.case_type == 'local':
        pass
    elif input_para.case_type in 'igsrt':
        Tools.preprocess(target_time.year, target_doy, input_para.save_pwd, input_para.download_list_fn) # create directory and copy needed file () to save_pwd
        gpsfile.download_data(input_para.case_type, input_para.download_list_fn)
    else:
        print "case_type error"

    gpsfile.decompress_file()
    gpsfile.d2o()

    # check gps ofile
    if input_para.case_type == 'igsrt':
        pass#stn_list = sorted(glob.glob('*{0:03}{1}{2:02}.{3:02}o'.format(target_doy, chr(97 + target_hh), target_mm, target_year % 100)))
    elif input_para.case_type == 'igs':
        stn_list = sorted(glob.glob("*{0:03}0.{1:02}o".format(target_doy, target_time.year % 100)))
    if stn_list == []:
        print "No GPS observation data on {0}.{1:03}".format(target_time.year, target_doy)
        break

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
    #### start to process o file
    for stn_num, station in enumerate(stn_list, start=1):
        target_stn = station[0:4]
        print "DOY: {3:03}, Station: {0} ({1:4}/{2:4}) ".format(target_stn, stn_num, len(stn_list), target_doy),
        st_time = datetime.datetime.now()
        ofile = Getdata.GPSofileData(target_stn, target_time.year, target_doy)
        stnx, stny, stnz = ofile.read_ofile_xyz()
        if stnx == 0.0 or stny == 0.0 or stnz ==0.0:
            if os.path.isfile('marker.crd'):
                stnifo = open('marker.crd','r')
                try:
                    rdline = stnifo.readline()
                    while not (rdline[0:4] == target_stn):
                        if rdline == "":
                            raise EOFError
                        rdline = stnifo.readline()
                    stnx = float(rdline[ 5:19])
                    stny = float(rdline[20:34])
                    stnz = float(rdline[35:49])
                except EOFError:
                    print "Please write on the {0}'s XYZ information in marker.crd".format(target_stn)
                    continue
            else:
                print "The marker.crd is not exist!!"
                sys.exit()
        stn_lon, stn_lat, stn_alt = Tools.TransCondinate().xyz2g(stnx, stny, stnz)
        ####    load gps ofile P1, P2, L1, L2
        try:
            st = datetime.datetime.now()
            P1, P2, L1, L2, P1yn, dt = ofile.read_ofile_PLvalue()
            ed = datetime.datetime.now()
        except TypeError:
            print "the ofile has no enough data!"
            continue

        print ed - st
        print "pasue"

        e = DataProcess.elevation(satdata, stnx, stny, stnz)
        sate_bias = DataProcess.satellite_bias(target_time.year, target_doy, P1yn)

