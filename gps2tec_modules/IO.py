class Input_file:
    def __init__(self, filename):
        import datetime
        import sys
        self.fn = filename
        with open(self.fn, 'r') as fid:
            rdline = fid.readline()
            while rdline:
                if len(rdline.split()) == 0:
                    rdline = fid.readline()
                    continue
                if rdline.split('=')[0].strip() == 'case_type': self.case_type = rdline.split('=')[1].strip().lower()

                if rdline.split('=')[0].strip() == 'year':      self.year      = int(rdline.split('=')[1].strip())
                if rdline.split('=')[0].strip() == 'st_doy':    self.st_doy    = int(rdline.split('=')[1].strip())
                if rdline.split('=')[0].strip() == 'ed_doy':    self.ed_doy    = int(rdline.split('=')[1].strip())

                if rdline.split('=')[0].strip() == 'download_list_fn':
                    if len(rdline.split('=')) == 1:
                        self.download_list_fn = None
                    else:
                        self.download_list_fn = rdline.split('=')[1].strip()
                if rdline.split('=')[0].strip() == 'save_pwd':
                    if len(rdline.split('=')) == 1:
                        self.save_pwd = None
                    else:
                        self.save_pwd = rdline.split('=')[1].strip()

                if rdline.split('=')[0].strip() == 'stnbias_method':       self.stnbias_method  = rdline.split('=')[1].strip().lower()
                if rdline.split('=')[0].strip() == 'elevation_angle':      self.elevation_angle = int(rdline.split('=')[1].strip())
                if rdline.split('=')[0].strip() == 'h1_km':                self.h1 = float(rdline.split('=')[1].strip())
                if rdline.split('=')[0].strip() == 'h2_km':                self.h2 = float(rdline.split('=')[1].strip())
                if rdline.split('=')[0].strip() == 'multiprocess_number':  self.mp_num = int(rdline.split('=')[1].strip())
                rdline = fid.readline()

        self.total_run = self.ed_doy - self.st_doy + 1


class OutputData:
    def __init__(self, download_type, stn, year, doy, hh=0, mm=0):
        self.year = year
        self.doy = doy
        self.hh = hh
        self.mm = mm
        self.stn = stn
        self.types = download_type

        self.ofn = 'o{0}{1}{2:03}.dat'.format(stn, year, doy)
        self.afn = 'a{0}{1}{2:03}.dat'.format(stn, year, doy)
        self.vfn = 'v{0}{1}{2:03}.dat'.format(stn, year, doy)
        self.sfn = 's{0}{1}{2:03}.dat'.format(stn, year, doy)
        self.efn = 'e{0}{1}{2:03}.dat'.format(stn, year, doy)
        self.log = '{0}{1}{2:03}.log'.format(stn, year, doy)

    def output_obs(self, lon, lat, vtec, stec, e):
        import numpy as np
        import os
        time_num = vtec.shape[0]
        sat_num = vtec.shape[1]

        if self.types == 'igsrt' and os.path.isfile(self.vfn):
            st_idx = self.hh * 60 * 2 + self.mm * 2
            ed_idx = st_idx + 30
            o_before = np.loadtxt(self.ofn)
            a_before = np.loadtxt(self.afn)
            v_before = np.loadtxt(self.vfn)
            s_before = np.loadtxt(self.sfn)
            e_before = np.loadtxt(self.efn)

            o_before[st_idx:ed_idx, :] = lon[st_idx:ed_idx, :].copy()
            a_before[st_idx:ed_idx, :] = lat[st_idx:ed_idx, :].copy()
            v_before[st_idx:ed_idx, :] = vtec[st_idx:ed_idx, :].copy()
            s_before[st_idx:ed_idx, :] = stec[st_idx:ed_idx, :].copy()
            e_before[st_idx:ed_idx, :] = e[st_idx:ed_idx, :].copy()

            lon = o_before.copy()
            lat = a_before.copy()
            vtec = v_before.copy()
            stec = s_before.copy()
            e = e_before.copy()

        ofn = open(self.ofn, 'w')
        afn = open(self.afn, 'w')
        vfn = open(self.vfn, 'w')
        sfn = open(self.sfn, 'w')
        efn = open(self.efn, 'w')

        for ttt in range(time_num):
            for snv in range(sat_num):
                if stec[ttt, snv] < 0 or vtec[ttt, snv] > 500:
                    stec[ttt, snv] = 0
                    vtec[ttt, snv] = 0
                    lon[ttt, snv] = 0
                    lat[ttt, snv] = 0
                    e[ttt, snv] = 0
                ofn.write('{0:8.2f}'.format(lon[ttt, snv]))
                afn.write('{0:8.2f}'.format(lat[ttt, snv]))
                vfn.write('{0:8.2f}'.format(vtec[ttt, snv]))
                sfn.write('{0:8.2f}'.format(stec[ttt, snv]))
                efn.write('{0:8.2f}'.format(e[ttt, snv]))
            ofn.write('\n')
            afn.write('\n')
            vfn.write('\n')
            sfn.write('\n')
            efn.write('\n')

        ofn.close()
        afn.close()
        vfn.close()
        sfn.close()
        efn.close()

    def output_log(self, stn_lon, stn_lat, time_num, minTEC, stn_bias):
        with open(self.log, 'w') as fn:
            fn.write('Station name: {0}\n'.format(self.stn))
            fn.write('Year: {0}\n'.format(self.year))
            fn.write('Julian date: {0}\n'.format(self.doy))
            fn.write('Station location: {0:.2f} N/ {1:.2f} E\n'.format(stn_lat, stn_lon))
            fn.write('Epoch: {0}\n'.format(time_num))
            fn.write('Minimum of VTEC: {0:.3f} E16 (electrons/m2)\n'.format(minTEC))
            fn.write('Station BIAS: {0:.6f} (ns)\n'.format(stn_bias))

