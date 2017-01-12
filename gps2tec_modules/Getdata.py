"""
This is the module of gps2tec.py
This module is to get data 
[GIM, Navigation Data, P1C1, P1P2, GPS ofile, station Bias]

create by Wei-Han Chen @2016.02.16
"""
class GimData:
    def __init__(self, year, doy):
        self.year  = int(year)
        self.doy   = int(doy)
        self.mapfn = "MAP{0}{1:03}.gim".format(self.year, self.doy)
        self.sourcefn = "CODG{0:03}0.{1}I.Z".format(self.doy, str(self.year)[2:4])
 
  
    def file_exist(self):
        import os
        return os.path.isfile(self.mapfn)

    def get_gimdata(self):
        import urllib, os
        if self.file_exist():
            print "No need to download GIM data..."
            return
        print "Start to download GIM data..."
        weblink = "ftp://ftp.unibe.ch/aiub/CODE/{0}/".format(self.year)
        if not os.path.isfile(self.sourcefn[:-2]): 
            if not os.path.isfile(self.sourcefn):
                download = urllib.URLopener()
                download.retrieve(weblink+self.sourcefn, self.sourcefn)
            os.system("gzip -fd {0}".format(self.sourcefn)) 

    def create_mapfile(self):
        fid = open(self.sourcefn[:-2],'r')
        rdline = fid.readline()
        while not rdline[60:73] == "END OF HEADER":
           if rdline[60:78] == "EPOCH OF FIRST MAP":
               st_time = int(rdline[18:24])
               rdline  = fid.readline()
               ed_time = int(rdline[18:24])
               if ed_time == 0: ed_time = 24
           if rdline[60:77] == "# OF MAPS IN FILE": map_num  = int(rdline[0:6]) 
           if rdline[60:68] == "INTERVAL":          time_inv = int(rdline[0:6])/3600
           if rdline[60:78] == "LAT1 / LAT2 / DLAT":
               latst   = float(rdline[2:8])
               lated   = float(rdline[8:14])
               latstep = float(rdline[14:20])
               latnum  = int(abs((lated - latst)/latstep + 1))
           if rdline[60:78] == "LON1 / LON2 / DLON":
               lonst   = float(rdline[2:8])
               loned   = float(rdline[8:14])
               lonstep = float(rdline[14:20])
               lonnum  = int(abs((loned - lonst)/lonstep + 1))          
           rdline = fid.readline()
        fmap = open(self.mapfn,'w')   
        fmap.write("{0:1}{1:1}{2:02}{3:03}".format(time_inv, st_time, ed_time, latnum))
        for oc in range(lonnum):
            fmap.write("{0:>7}".format(int(lonst+(oc*lonstep))))
        rdline = fid.readline()
        while not (rdline[60:74] == "END OF TEC MAP" and int(rdline[0:6]) == map_num):
            if rdline[60:80] == "LAT/LON1/LON2/DLON/H":
                fmap.write("\n")
                fmap.write("{0:>7}".format(rdline[1:8]))
                for k in range(int(lonnum/16)+1):
                    rdline = fid.readline()
                    if k == int(lonnum/16):
                        for j in range(int(lonnum%16)):
                            fmap.write("{0:>7}".format(int(rdline[(j)*5:(j+1)*5])))
                    else:
                        for j in range(16):
                            fmap.write("{0:>7}".format(int(rdline[(j)*5:(j+1)*5])))    
            rdline = fid.readline()
        fid.close()
        fmap.close()

    def min_gimTEC(self, target_lon, target_lat):
        minTEC  = 8888.0
        minTime =  888.0
        fid = open(self.mapfn,'r')
        rdline = fid.readline()
        time_inv  = int(rdline[0:1])
        st_time   = int(rdline[1:2])
        ed_time   = int(rdline[2:4])
        latnum    = int(rdline[4:7])
        lonst = float(rdline[7:14])
        loned = float(rdline[-8:-1])
        lonstep = float(rdline[14:21]) - lonst
        rdline = fid.readline()
        latst  = float(rdline[0:7])
        rdline = fid.readline()
        latstep = float(rdline[0:7]) - latst
        lated = latst + (latnum-1)*latstep
        fid.close()

        lon_left_idx   = int((target_lon - lonst)/lonstep) + 1
        lon_right_idx  = lon_left_idx + 1
        x_ratio = (target_lon - (lonst+(lon_left_idx-1)*lonstep) )/lonstep
        if target_lat >= latst:
            lat_buttom_idx = 1
            lat_top_idx    = 1
            y_ratio        = 0
        elif target_lat <= lated:
            lat_buttom_idx = int((lated - latst)/latstep) + 1
            lat_top_idx    = lat_buttom_idx
            y_ratio        = 1
        else:
            lat_top_idx    = int((target_lat - latst)/latstep) + 1
            lat_buttom_idx = lat_top_idx + 1
            y_ratio        = (target_lat - (latst+(lat_top_idx-1)*latstep) )/latstep
        
        fid = open(self.mapfn,'r')
        rdline = fid.readline()
        lat_idx = 1

        time_num = ((ed_time - st_time)/time_inv)+1
        time_idx = 1
        target_latidx = 1 + lat_top_idx + latnum*(time_idx-1)
        while True:
            rdline = fid.readline();lat_idx += 1
            if lat_idx == target_latidx:
                a   = float(rdline[lon_left_idx*7 : (lon_left_idx+1)*7])/10.0 
                b   = float(rdline[lon_right_idx*7 : (lon_right_idx+1)*7])/10.0
                rdline = fid.readline();lat_idx += 1
                if rdline:
                    c   = float(rdline[lon_left_idx*7 : (lon_left_idx+1)*7])/10.0
                    d   = float(rdline[lon_right_idx*7 : (lon_right_idx+1)*7])/10.0
                else:
                    c   = a
                    d   = b
                tec = d*x_ratio*y_ratio + c*y_ratio - c*x_ratio*y_ratio - b*x_ratio*y_ratio - a*y_ratio + a*x_ratio*y_ratio + b*x_ratio + a - a*x_ratio
                if tec < minTEC:
                    minTEC  = tec
                    minTime = (st_time + (time_idx-1)*time_inv) + target_lon/15.0
                    if minTime >= 24:
                        minTime = minTime - 24
                    elif minTime < 0:
                        minTime = minTime + 24
                time_idx += 1
                target_latidx = 1 + lat_top_idx + latnum*(time_idx-1)
            if time_idx > time_num: break 

        fid.close()
        return minTEC, minTime

    def read_station_bias(self):
        import os
        stnbias_fn = 'bias{0}{1:03}.dat'.format(self.year, self.doy)
        if os.path.isfile(stnbias_fn): return
        bias_fid = open (stnbias_fn, 'w')
        with open(self.sourcefn[:-2],'r') as fid:
            rdline = fid.readline() 
            while not rdline[60:73] == "END OF HEADER":
                if rdline[3:4]=="G" and rdline[60:80] == "STATION / BIAS / RMS":
                    stn  = rdline[6:10].lower()
                    bias = float(rdline[28:36])*0.299792458
                    bias_fid.write('{0}{1:>11.4f}{2:>10}\n'.format(stn, bias, 'CODG'))
                rdline = fid.readline()
            bias_fid.close()
        
        
class NavigationData: 
    def __init__(self, year, doy, hh=0, mm=0, download_type='igs'):
        import datetime, time
        self.year  = int(year)
        self.doy   = int(doy)
        self.mon   = (datetime.datetime(self.year, 1, 1)+datetime.timedelta(days =self.doy-1)).month
        self.day   = (datetime.datetime(self.year, 1, 1)+datetime.timedelta(days =self.doy-1)).day
        self.hh    = int(hh)
        self.mm    = int(mm)
        self.types = download_type
        if self.types=='igs':
            target_time1 = datetime.datetime(self.year, 1, 1) + datetime.timedelta(days =self.doy-1)
            target_time2 = datetime.datetime(self.year, 1, 1) + datetime.timedelta(days =self.doy)
            self.dweeks1 = int((target_time1 - datetime.datetime(1980,1,6)).days/7)
            self.ddays1  = int((target_time1 - datetime.datetime(1980,1,6)).days%7)
            self.dweeks2 = int((target_time2 - datetime.datetime(1980,1,6)).days/7)
            self.ddays2  = int((target_time2 - datetime.datetime(1980,1,6)).days%7)
            self.sourcefn_igs1 = "igs{0:04}{1}.sp3.Z".format(self.dweeks1, self.ddays1)
            self.sourcefn_igr1 = "igr{0:04}{1}.sp3.Z".format(self.dweeks1, self.ddays1)
            self.sourcefn_igs2 = "igs{0:04}{1}.sp3.Z".format(self.dweeks2, self.ddays2)
            self.sourcefn_igr2 = "igr{0:04}{1}.sp3.Z".format(self.dweeks2, self.ddays2)
        # for IGSRT 
        elif self.types=='igsrt':
            target_time = datetime.datetime(self.year,1,1,self.hh)+datetime.timedelta(days=self.doy-1)-datetime.timedelta(hours=8)
            dweeks = int((target_time - datetime.datetime(1980,1,6)).days/7)
            ddays  = int((target_time - datetime.datetime(1980,1,6)).days%7)
            self.sourcefn_igu = "igu{0:04}{1}_{2:02}.sp3.Z".format(dweeks, ddays, int(target_time.hour/6)*6)

    def file_exist(self):
        import os
        exist = True
        if self.types=='igs':
            if os.path.isfile(self.sourcefn_igs1[:-2]):
                self.sourcefn1 = self.sourcefn_igs1[:-2]    
            elif os.path.isfile(self.sourcefn_igr1[:-2]):
                self.sourcefn1 = self.sourcefn_igr1[:-2]
            else:
                exist = False  
            if os.path.isfile(self.sourcefn_igs2[:-2]):
                self.sourcefn2 = self.sourcefn_igs2[:-2]
            elif os.path.isfile(self.sourcefn_igr2[:-2]):
                self.sourcefn2 = self.sourcefn_igr2[:-2]
            else:
                exist = False
        elif self.types=='igsrt':
            if os.path.isfile(self.sourcefn_igu[:-2]):
                self.sourcefn = self.sourcefn_igu[:-2]
            else:
                exist = False 
        return exist

    def get_navidata(self):
        import urllib, os
        if self.file_exist():
            print "No need to download Navigation data..."
            return

        print "Start to download Navigation data..."
        if self.types=='igs':
            weblink = "ftp://igscb.jpl.nasa.gov/pub/product/"
            if not (os.path.isfile(self.sourcefn_igs1) or os.path.isfile(self.sourcefn_igr1)):
                try:
                    download = urllib.URLopener()
                    download.retrieve("{0}{1:04}/{2}".format(weblink, self.dweeks1, self.sourcefn_igs1), self.sourcefn_igs1)
                    self.sourcefn1 = self.sourcefn_igs1[:-2]
                except IOError:
                    download = urllib.URLopener()
                    download.retrieve("{0}{1:04}/{2}".format(weblink, self.dweeks1, self.sourcefn_igr1), self.sourcefn_igr1)
                    self.sourcefn1 = self.sourcefn_igr1[:-2]
            if not (os.path.isfile(self.sourcefn_igs2) or os.path.isfile(self.sourcefn_igr2)):
                try:
                    download = urllib.URLopener()
                    download.retrieve("{0}{1:04}/{2}".format(weblink, self.dweeks2, self.sourcefn_igs2), self.sourcefn_igs2)
                    self.sourcefn2 = self.sourcefn_igs2[:-2]
                except IOError:
                    download = urllib.URLopener()
                    download.retrieve("{0}{1:04}/{2}".format(weblink, self.dweeks2, self.sourcefn_igr2), self.sourcefn_igr2)
                    self.sourcefn2 = self.sourcefn_igr2[:-2]
        elif self.types=='igsrt':
            weblink = "ftp://cddis.gsfc.nasa.gov/pub/gps/products/{0}/".format(self.sourcefn_igu[3:7])
            download = urllib.URLopener()
            download.retrieve(weblink+self.sourcefn_igu, self.sourcefn_igu)
            self.sourcefn = self.sourcefn_igu[:-2]
        os.system("gzip -fd *sp3.Z")

    def read_navidata(self):
        import numpy as np
        import datetime
        from scipy import interpolate
        #from DataProcess import cubspl
        sp3     = np.zeros((97,32,3)) 
        satdata = np.zeros((2880,32,3))
        x       = np.zeros(97)
        y       = np.zeros(97)

        if self.types=='igs':
            fid = open(self.sourcefn1,'r')
            rdline = fid.readline()
            while not rdline == "":
                if rdline[0:3] == "*  ":
                    time = int(rdline[14:16])*3600 + int(rdline[17:19])*60 + int(rdline[20:22])
                    time = time /(15*60) + 1
                if rdline[0:1] == "P" and int(rdline[2:4]) > 0:    
                    snv = int(rdline[2:4])
                    sp3[time-1, snv-1, 0] = float(rdline[4:18])  * 1000
                    sp3[time-1, snv-1, 1] = float(rdline[18:32]) * 1000
                    sp3[time-1, snv-1, 2] = float(rdline[32:46]) * 1000
                rdline = fid.readline()
            fid.close()

            fid = open(self.sourcefn2,'r')
            rdline = fid.readline()
            i = 0
            while not rdline == "":
                if rdline[0:3] == "*  ":
                    time = int(rdline[14:16])*3600 + int(rdline[17:19])*60 + int(rdline[20:22])
                    time = time /(15*60) + 1
                    i = i + 1
                    if i == 2: break
                if rdline[0:1] == "P" and int(rdline[2:4]) > 0: 
                    snv = int(rdline[2:4])
                    sp3[96, snv-1, 0] = float(rdline[4:18])  * 1000
                    sp3[96, snv-1, 1] = float(rdline[18:32]) * 1000
                    sp3[96, snv-1, 2] = float(rdline[32:46]) * 1000
                rdline = fid.readline()
            fid.close()
        elif self.types=='igsrt':
            target_next = datetime.datetime(self.year, self.mon, self.day) + datetime.timedelta(days=1)
            load_P      = False
            fid = open(self.sourcefn,'r')
            rdline = fid.readline()
            while rdline:
                if rdline[0:3] == "*  ":
                    load_P      = False
                    time = int(rdline[14:16])*3600 + int(rdline[17:19])*60 + int(rdline[20:22])
                    if int(rdline[8:10])!=self.mon or int(rdline[11:13])!=self.day:
                        if int(rdline[8:10])==target_next.month and int(rdline[11:13])==target_next.day: break
                        rdline = fid.readline()
                        continue
                    load_P      = True
                    time = time /(15*60) + 1
                if rdline[0:1] == "P" and int(rdline[2:4]) > 0 and load_P:
                    snv = int(rdline[2:4])
                    sp3[time-1, snv-1, 0] = float(rdline[4:18])  * 1000
                    sp3[time-1, snv-1, 1] = float(rdline[18:32]) * 1000
                    sp3[time-1, snv-1, 2] = float(rdline[32:46]) * 1000
                rdline = fid.readline()

            rdline = fid.readline()
            while rdline[0:3] !="*  " and rdline:
                if rdline[0:1] == "P" and int(rdline[2:4]) > 0:
                   snv = int(rdline[2:4])
                   sp3[96, snv-1, 0] = float(rdline[4:18])  * 1000
                   sp3[96, snv-1, 1] = float(rdline[18:32]) * 1000
                   sp3[96, snv-1, 2] = float(rdline[32:46]) * 1000
                rdline = fid.readline()
            
            fid.close()
        
        for k in range(32):
            for j in range(3):
                for i in range(97):
                    x[i] = 0 + (i)*15*60
                    y[i] = sp3[i, k, j]
                tck = interpolate.splrep(x, y, s=0)
                xnew = np.arange(0, 2880*30,30)
                temp = interpolate.splev(xnew, tck, der=0)
                #temp = cubspl(x, y, 97, 2880)
                for i in range(2880):
                    satdata[i, k, j] = temp[i]
        return satdata

class P1P2codeData:
    def __init__(self, year, doy):
        import datetime, time
        from dateutil.relativedelta import relativedelta
        self.year   = int(year)
        self.doy    = int(doy)
        target_time = datetime.datetime(self.year, 1, 1) + datetime.timedelta(days =self.doy-1)
        delta = relativedelta(months=+1)
        self.mon    = int(target_time.month)
        self.year_bom  = (target_time-delta).year
        self.mon_bom   = (target_time-delta).month
        self.sourcefn = "P1P2{0}{1:02}.DCB.Z".format(str(self.year)[2:4], self.mon)
        self.sourcefn_bom = "P1P2{0}{1:02}.DCB.Z".format(str(self.year_bom)[2:4], self.mon_bom)

    def file_exist(self):
        import os
        return os.path.isfile(self.sourcefn[:-2])
    
    def get_p1p2data(self):
        import urllib, os
        if self.file_exist():
            print "No need to download P1P2 DCB data..."
            return
        print "Start to download P1P2 DCB data..."
        weblink = "ftp://ftp.unibe.ch/aiub/CODE/{0}/".format(self.year)
        if not os.path.isfile(self.sourcefn):
            try:
                download = urllib.URLopener()
                download.retrieve(weblink+self.sourcefn, self.sourcefn)
            except IOError:
                weblink = "ftp://ftp.unibe.ch/aiub/CODE/{0}/".format(self.year_bom)
                download = urllib.URLopener()
                download.retrieve(weblink+self.sourcefn_bom, self.sourcefn)
                
            
        os.system("gzip -fd {0}".format(self.sourcefn)) 

class P1C1codeData:
    def __init__(self, year, doy):
        import datetime, time
        from dateutil.relativedelta import relativedelta
        self.year   = int(year)
        self.doy    = int(doy)
        target_time = datetime.datetime(self.year, 1, 1) + datetime.timedelta(days =self.doy-1)
        delta = relativedelta(months=+1)
        self.mon    = int(target_time.month)
        self.year_bom  = (target_time-delta).year
        self.mon_bom   = (target_time-delta).month
        self.sourcefn = "P1C1{0}{1:02}.DCB.Z".format(str(self.year)[2:4], self.mon)
        self.sourcefn_bom = "P1C1{0}{1:02}.DCB.Z".format(str(self.year_bom)[2:4], self.mon_bom)

    def file_exist(self):
        import os
        return os.path.isfile(self.sourcefn[:-2])

    def get_p1c1data(self):
        import urllib, os
        if self.file_exist():
            print "No need to download P1C1 DCB data..."
            return
        print "Start to download P1C1 DCB data..."
        weblink = "ftp://ftp.unibe.ch/aiub/CODE/{0}/".format(self.year)
        if not os.path.isfile(self.sourcefn):
            try:
                download = urllib.URLopener()
                download.retrieve(weblink+self.sourcefn, self.sourcefn)
            except IOError:
                weblink = "ftp://ftp.unibe.ch/aiub/CODE/{0}/".format(self.year_bom)
                download = urllib.URLopener()
                download.retrieve(weblink+self.sourcefn_bom, self.sourcefn)

        os.system("gzip -fd {0}".format(self.sourcefn))

class GPSourceFile:
    def __init__(self, year, doy, hh=0, mm=0):
        self.year = year
        self.doy  = doy
        self.hh   =  hh
        self.mm   =  mm

    def download_data(self, download_type, download_file):
        from ftplib import FTP
        import sys
        import multiprocessing

        ####    different type to different ftp host, path, and filename list
        if download_type == 'igs':
            ftp_host = 'garner.ucsd.edu'
            ftp_pwd  = '/archive/garner/rinex/{0}/{1:03}'.format(self.year, self.doy)
            ftp_fn_filter = '*{1:03}0.{0:02}d.Z'.format(self.year%100, self.doy) 
        elif download_type == 'igsrt':
            ftp_host = 'cddis.gsfc.nasa.gov'
            ftp_pwd  = '/pub/gps/data/highrate/{0}/{2:03}/{1:02}d/{3:02}/'.format(self.year, self.year%100, self.doy, self.hh)
            ftp_fn_filter = '*{0:03}{1}{2:02}.{3:02}d.Z'.format(self.doy, chr(97+self.hh), self.mm, self.year%100)
            download_file = 'bias{0}{1:03}.dat'.format(self.year, self.doy)
        ####    connect to ftp
        try:
            ftp = FTP(ftp_host)
            ftp.login()
            ftp.cwd(ftp_pwd)
        except:
            print 'FTP connect to "', ftp_host, '" get some error !'
            sys.exit()

        filenames = []
        if download_file == '':
            filenames=ftp.nlst(ftp_fn_filter)
        else:
            with open(download_file, 'r') as f:
                file_temp = f.readlines()
                for stn in file_temp:
                    filenames.append('{0}{1}'.format(stn[:4], ftp_fn_filter[1:]))

        ftp.quit()
        total_file = len(filenames)
        interval = total_file // 2
        process_list =[]
        print 'Start to download GPS file from "{0}"'.format(ftp_host)
        for i in range(2):
            st = i*interval
            ed = i*interval+interval
            if i == 1 : ed += total_file % 2
            ftp = FTP(ftp_host)
            ftp.login()
            ftp.cwd(ftp_pwd)
            p = multiprocessing.Process(target=self.multi_download, args=(ftp, filenames[st:ed], ftp_fn_filter))
            p.start()
            process_list.append(p)

        for job in process_list: job.join()

    def multi_download(self, ftp, filenames, ftp_fn_filter):
        import os
        total_file = len(filenames)
        file_count = 1
        for filename in filenames:
            print "Processing:  {1} -- ".format(float(file_count)*100/total_file, filename[:4]),
            if os.path.isfile(filename) or os.path.isfile("{0}{1}".format(filename[:4], ftp_fn_filter[1:-2])) or os.path.isfile("{0}{1}o".format(filename[:4], ftp_fn_filter[1:-3])):
                print 'Already have the file !'
                file_count+=1
                continue
            save = open(filename, 'wb')
            try:
                ftp.retrbinary('RETR '+ filename, save.write)
                file_count+=1
                print "Download sucessful !"
            except:
                print "This station does not exist in FTP !"
                file_count+=1
                os.remove(filename)
            save.close()
        ftp.quit()


    def decompress_file(self):
        import os
        import glob
        print "Start to decompress file ..."
        filelist = glob.glob("*{1:03}*.{0:02}d.Z".format(self.year%100, self.doy))
        for filename in filelist:
            os.system("gzip -fd {0}".format(filename))

    def d2o(self):
        import os
        import glob
        print "Start to do d2o ..."
        filelist = glob.glob("*{1:03}*.{0:02}d".format(self.year%100, self.doy))
        if os.path.isfile('crx2rnx'):
            for filename in filelist:
                os.system("./crx2rnx {0}".format(filename))
                os.remove(filename)
        else:
            print 'The file "crx2rnx" not exist !'


class GPSofileData:
    def __init__(self, stn, year, doy, hh=0, mm=0, ofile_type='igs'):
        import datetime, time
        self.year   = int(year)
        self.doy    = int(doy)
        target_time = datetime.datetime(self.year, 1, 1) + datetime.timedelta(days =self.doy-1)
        self.mon    = int(target_time.month)
        self.day    = int(target_time.day)
        self.hh     = hh
        self.mm     = mm
        if ofile_type == 'igs':
            self.ofile  = "{0}{2:03}0.{1}o".format(stn, str(year)[2:4], self.doy)
        elif ofile_type == 'igsrt':
            self.ofile  = '{0}{1:03}{2}{3:02}.{4:02}o'.format(stn, self.doy, chr(97+self.hh), self.mm, self.year%100)

    def read_ofile_xyz(self):
        fid = open(self.ofile, 'r')
        while True:
            rdline = fid.readline()
            if rdline[60:79] == "APPROX POSITION XYZ":    
                self.stn_x = float(rdline[ 0:14])
                self.stn_y = float(rdline[14:28])
                self.stn_z = float(rdline[28:42])
                break
        fid.close() 
        return self.stn_x, self.stn_y, self.stn_z

    def read_ofile_PLvalue(self, dt=30):
        import numpy as np

        total_step = 86400 / dt
        P1 = np.zeros((total_step, 32))
        P2 = np.zeros((total_step, 32))
        L1 = np.zeros((total_step, 32))
        L2 = np.zeros((total_step, 32))

        with open(self.ofile, 'r') as fid:
            ofile_data = fid.readlines()

        data_line = 0
        obs_dt = 0
        obs_num = 0
        obs_list = []
        while not ofile_data[data_line][60:73] == "END OF HEADER":
            if ofile_data[data_line][60:68] == "INTERVAL":   # read GPS data time interval
                obs_dt = float(ofile_data[data_line][0:6])
            if ofile_data[data_line][60:79] == "# / TYPES OF OBSERV": # read numbers of GPS observation and kind
                if not obs_num: obs_num = int(ofile_data[data_line][4:6])
                if obs_num < 4:
                    return
                else:
                    obs_list.extend(ofile_data[data_line][10:60].split())
            data_line += 1

        each_obs_lines = ((obs_num-1)//5)+1
        data_line += 1

        if not obs_dt:
            # calculate time interval
            sate_num = int(ofile_data[data_line][30:32])
            first_sec = int(ofile_data[data_line][9:12]) * 3600 + int(ofile_data[data_line][12:15]) * 60 + int(ofile_data[data_line][15:18])
            second_idx = data_line + each_obs_lines * sate_num + (sate_num - 1) // 12 + 1
            second_sec = int(ofile_data[second_idx][9:12]) * 3600 + int(ofile_data[second_idx][12:15]) * 60 + int(ofile_data[second_idx][15:18])
            obs_dt = second_sec - first_sec

        if not ('L1' in obs_list and 'L2' in obs_list and 'P2' in obs_list and ('P1' in obs_list or 'C1' in obs_list)):  # if L1 L2 P2 and (C1 or P1) exist, start to read ofile
            return
        else:
            P1_exist = True if ('P1' in obs_list) else False
            P1_err_nums = 0
            while data_line < len(ofile_data):
                sate_list = []
                sate_dataline_list = []
                if int(ofile_data[data_line][28:29]) != 0:
                    if ofile_data[data_line][33:35] == '':
                        data_line += int(ofile_data[data_line][29:32]) + 1
                    else:
                        data_line += each_obs_lines * int(ofile_data[data_line][30:32]) + (int(ofile_data[data_line][30:32])-1)//12 + 1
                    continue
                current_second = int(ofile_data[data_line][9:12]) * 3600 + int(ofile_data[data_line][12:15]) * 60 + int(ofile_data[data_line][15:18])
                sate_num = int(ofile_data[data_line][30:32])
                if current_second % dt == 0:
                    time_idx = current_second//dt
                    # load satellites number [start]
                    for i in range(sate_num//12):
                        for j in range(12):
                            if ofile_data[data_line][32+j*3:33+j*3] == "G" or ofile_data[data_line][32+j*3:33+j*3] == " ":
                                sate_list.append(int(ofile_data[data_line][33+j*3:35+j*3]))
                                sate_dataline_list.append(i*12+j)
                        data_line += 1
                    for i in range(sate_num%12):
                        if ofile_data[data_line][32+i*3:33+i*3] == "G" or ofile_data[data_line][32+i*3:33+i*3] == " ":
                            sate_list.append(int(ofile_data[data_line][33+i*3:35+i*3]))
                            sate_dataline_list.append((sate_num//12)*12+i)
                    # load satellites number [end]
                    if sate_num%12 != 0: data_line += 1

                    # load P1 P2 L1 L2 C1 code data [start]
                    code_data_line_st = data_line
                    for i in range(len(sate_list)):

                        # L1
                        target_idx = obs_list.index('L1') + 1
                        target_data_line = code_data_line_st + sate_dataline_list[i] * each_obs_lines + (target_idx - 1) // 5
                        if target_idx > 5: target_idx -= 5
                        try:
                            L1[time_idx, sate_list[i] - 1] = float(ofile_data[target_data_line][1 + 16 * (target_idx - 1):14 + 16 * (target_idx - 1)])
                        except ValueError:
                            L1[time_idx, sate_list[i] - 1] = 0

                        # L2
                        target_idx = obs_list.index('L2') + 1
                        target_data_line = code_data_line_st + sate_dataline_list[i] * each_obs_lines + (target_idx - 1) // 5
                        if target_idx > 5 : target_idx -= 5
                        try:
                            L2[time_idx, sate_list[i] - 1] = float(ofile_data[target_data_line][1 + 16 * (target_idx - 1):14 + 16 * (target_idx - 1)])
                        except ValueError:
                            L2[time_idx, sate_list[i] - 1] = 0

                        # P2
                        target_idx = obs_list.index('P2') + 1
                        target_data_line = code_data_line_st + sate_dataline_list[i] * each_obs_lines + (target_idx - 1) // 5
                        if target_idx > 5: target_idx -= 5
                        try:
                            P2[time_idx, sate_list[i] - 1] = float(ofile_data[target_data_line][1 + 16 * (target_idx - 1):14 + 16 * (target_idx - 1)])
                        except ValueError:
                            P2[time_idx, sate_list[i] - 1] = 0

                        # P1
                        try:
                            target_idx = obs_list.index('P1') + 1
                            target_data_line = code_data_line_st + sate_dataline_list[i] * each_obs_lines + (target_idx - 1) // 5
                            if target_idx > 5: target_idx -= 5
                            P1[time_idx, sate_list[i] - 1] = float(ofile_data[target_data_line][1 + 16 * (target_idx - 1):14 + 16 * (target_idx - 1)])
                        except ValueError:
                            P1_err_nums += 1
                            target_idx = obs_list.index('C1') + 1
                            target_data_line = code_data_line_st + sate_dataline_list[i] * each_obs_lines + (target_idx - 1) // 5
                            if target_idx > 5: target_idx -= 5
                            try:
                                P1[time_idx, sate_list[i] - 1] = float(ofile_data[target_data_line][1 + 16 * (target_idx - 1):14 + 16 * (target_idx - 1)])
                            except ValueError:
                                P1[time_idx, sate_list[i] - 1] = 0

                    data_line += each_obs_lines * sate_num
                else:
                    data_line += each_obs_lines * sate_num + (sate_num-1)//12 + 1
                #print data_line
        if P1_err_nums > total_step:
            P1_exist = False
        #if obs_dt==-999:
        #    obs_dt = check_dt
        if obs_dt < dt: obs_dt = dt

        return P1, P2, L1, L2, P1_exist, obs_dt

class StationBiasData:
    def __init__(self, year, doy):
        self.biasfn = 'bias{0}{1:03}.dat'.format(year, doy)

    def file_exist(self):
        import os
        return os.path.isfile(self.biasfn)
        
    def read_stnbias(self, stn):
        with open(self.biasfn,'r') as f:
            rdlist=f.readlines()

        for i in range(len(rdlist)):
            if rdlist[i].split()[0] == stn:
                return float(rdlist[i].split()[1])

        return 0.0

    def write_stnbias(self, stn, stn_bias, write_stnbs, minTEC, GIMminTEC, method):
        
        if (not write_stnbs) and (abs(minTEC-GIMminTEC)<1 or method =='fitting'):
            with open(self.biasfn,'a') as f:
                f.write('{0}{1:>11.4f}{2:>10}\n'.format(stn, stn_bias, method.upper()))

         


    


