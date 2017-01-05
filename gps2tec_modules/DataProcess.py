def preprocess(year, doy, save_pwd, download_list_fn):
    import os, shutil
    if save_pwd != '':
        if save_pwd[-1]!= '/': save_pwd += '/'
    if not os.path.isdir('{0}{1}/{2:03}/'.format(save_pwd, year, doy)):
        os.makedirs('{0}{1}/{2:03}/'.format(save_pwd, year, doy))
    shutil.copy('crx2rnx','{0}{1}/{2:03}/'.format(save_pwd, year, doy))
    shutil.copy('marker.crd','{0}{1}/{2:03}/'.format(save_pwd, year, doy))
    if download_list_fn !='':
        shutil.copy(download_list_fn,'{0}{1}/{2:03}/'.format(save_pwd, year, doy))
    os.chdir('{0}{1}/{2:03}/'.format(save_pwd, year, doy))


def xyz2g(x, y, z):
    import numpy as np
    a = 6378.137
    f = 1 / 298.257223563
    Rad = 180 / np.pi

    e2 = 2 * f - f * f
    b = (1 - f) * a
    r = np.sqrt(x * x + y * y)
    Np = z
    for i in range(7):
        phi = np.arctan((z + e2 * Np) / r)
        N = a / np.sqrt(1 - e2 * np.sin(phi) * np.sin(phi))
        Np = N * np.sin(phi)
    alt = r / np.cos(phi) - N
    lon = np.arctan2(y, x) * Rad
    lat = phi * Rad
    return lon, lat, alt


def elevation(satpoi, stnx, stny, stnz):
    import numpy as np
    stnpoi2 = np.sqrt(stnx*stnx + stny*stny + stnz*stnz)
    abar = np.zeros(3)
    e    = np.zeros((2880, 32))
    rad = np.pi/180
    for snv in range(32):
        for i in range(2880):
            abar[0] = satpoi[i, snv, 0] - stnx
            abar[1] = satpoi[i, snv, 1] - stny
            abar[2] = satpoi[i, snv, 2] - stnz
            abar2 = np.sqrt(abar[0]*abar[0] + abar[1]*abar[1] + abar[2]*abar[2])
            e[i, snv] = np.arccos((stnx*abar[0] + stny*abar[1] + stnz*abar[2])/(stnpoi2*abar2))/rad
            e[i, snv] = 90 - e[i, snv]
    return e


def satellite_bias(year, doy, P1yn):
    import numpy as np
    import datetime
    import os
    bias = np.zeros(32)
    cn = 299792458. / 1000000000.
    target_time = datetime.datetime(year, 1, 1) + datetime.timedelta(days =doy-1)
    mon = int(target_time.month)
    fn = "P1P2{0}{1:02}.DCB".format(str(year)[2:4], mon)
    if not os.path.isfile(fn):
        fn = "P1P2{0}{1:02}.DCB".format(str(year)[2:4], mon-1)
    with open(fn, 'r') as fid:
        rdline = fid.readlines()
    for i in range(len(rdline)):
        if rdline[i][0:1] == 'G':
            snv = int(rdline[i][1:3])
            if snv < 33:
                bias[snv-1] = float(rdline[i][26:36])*cn

    if not P1yn:
        fn = "P1C1{0}{1:02}.DCB".format(str(year)[2:4], mon)
        if not os.path.isfile(fn):
            fn = "P1C1{0}{1:02}.DCB".format(str(year)[2:4], mon-1)
        with open(fn, 'r') as fid:
            rdline = fid.readlines()
        for i in range(len(rdline)):
            if rdline[i][0:1] == 'G':
                snv = int(rdline[i][1:3])
                if snv < 33:
                    bias[snv-1] = bias[snv-1] - float(rdline[i][26:36])*cn
    return bias


def cal_PLcode(P1, P2, L1, L2, sat_bias, stn_bias, e, target_elve):
    import numpy as np
    sat_num = P1.shape[1]
    time_num = P1.shape[0]
    Pcode = np.zeros((time_num, sat_num))
    Lcode = np.zeros((time_num, sat_num))
    f1 = 1575.42e6
    f2 = 1227.60e6
    c = 2.99792458e8
    fat = ((f1 ** 2) * (f2 ** 2) / (f1 ** 2 - f2 ** 2)) / (40.3e16)

    for j in range(32):
        for i in range(2880):
            if e[i, j] >= target_elve and P1[i, j] != 0 and P2[i, j] != 0 and L1[i, j] != 0 and L2[i, j] != 0:
                Pcode[i, j] = (P2[i, j] - P1[i, j] + sat_bias[j] + stn_bias) * fat
                # if P1[i,j] == P2[i,j]: Pcode[i,j]=0
                Lcode[i, j] = (c / f1 * L1[i, j] - c / f2 * L2[i, j]) * fat
        Pcode_sort = np.sort(Pcode[np.where(Pcode != 0)])
        if len(Pcode_sort) == 0: continue
        Pcode_mid = Pcode_sort[np.ceil(len(Pcode_sort) / 2) - 1]
        for i in range(1, 2880):
            if (Pcode[i, j] != 0) and (abs(Pcode[i, j] - Pcode[i - 1, j]) > 500 or abs(Pcode[i, j] - Pcode_mid) > 500):
                Pcode[i, j] = 0

    return Pcode, Lcode


def cal_stec(Pcode, Lcode, dt):
    import numpy as np
    time_num = Pcode.shape[0]
    sat_num = Pcode.shape[1]
    stec = np.zeros((time_num, sat_num))

    length = 2880 * 30 / dt

    Pdata = np.zeros((length + 1, sat_num))
    Ldata = np.zeros((length + 1, sat_num))
    temp2 = np.zeros(length)
    index = np.zeros(4)
    x = np.zeros(length)
    Pdata[0:-1, :] = Pcode.copy()
    Ldata[0:-1, :] = Lcode.copy()
    for snv in range(32):
        # Lcode
        temp2 = Ldata[:-1, snv].copy()
        temp = 0.0
        for i in range(3, length):
            count = 0
            check = False
            temp3 = 0.0
            if Ldata[i, snv] != 0:
                count2 = 0
                for j in range(i - 3, i):
                    if temp2[j] != 0:
                        check = True
                        temp3 = temp2[j]
                        count2 += 1
                        index[count2 - 1] = j
                if check and temp3 != 0 and abs(temp2[i] - temp3) > 1.5:
                    if snv == 2:
                        index[count2 - 1] = j
                    if count2 >= 2:
                        temp4 = temp2[int(index[int(count2 - 1 - 1)])]
                        if abs(temp4 - temp3) < 1.5:
                            temp = temp + (temp2[i] - temp3 + (temp4 - temp3) * (i - int(index[count2 - 1])))
                        else:
                            temp = temp + (temp2[i] - temp3)
                    else:
                        temp = temp + (temp2[i] - temp3)
                Ldata[i, snv] = Ldata[i, snv] - temp
        # Pcode
        temp2[0:3] = Pdata[0:3, snv].copy()
        temp2[-3:] = Pdata[-3:, snv].copy()
        for i in range(3, length - 3):
            if len(np.where(Pdata[i - 3:i + 3 + 1, snv] != 0)[0]) > 2:
                temp2[i] = np.mean(Pdata[np.where(Pdata[i - 3:i + 3 + 1, snv] != 0)[0] + i - 3, snv])
            else:
                temp2[i] = 0.0
        Pdata[3:-4, snv] = temp2[3:-3].copy()
        Ldata[np.where(Pdata[3:-4, snv] == 0)[0] + 3, snv] = 0
        Pdata[0:3, snv] = 0.0

    # STEC Offset by Pdata and Ldata
    for snv in range(32):
        temp2[0] = 0.0
        count = 0
        temp = 0.0
        first_run = True
        # temp_tmp=0.0
        for i in range(1, length):
            temp2[i] = Ldata[i, snv]
            if Pdata[i, snv] != 0 and Ldata[i, snv] != 0:
                if Pdata[i - 1, snv] != 0 and Ldata[i - 1, snv] != 0 and abs(
                                Ldata[i, snv] - Ldata[i - 1, snv]) <= 5 or count == 0:
                    count += 1
                    x[count - 1] = i
                    temp = temp + (Pdata[i, snv] - Ldata[i, snv])
                if (Pdata[i + 1, snv] == 0 or Ldata[i + 1, snv] == 0) and count != 0:
                    temp = temp / count
                    for j in range(count):
                        temp2[int(x[j])] = Ldata[int(x[j]), snv] + temp
                    if first_run == True and temp2[1] != 0 and Ldata[0, snv] != 0:
                        temp2[0] = Ldata[0, snv] + temp
                        temp2[1] = Ldata[1, snv] + temp
                        temp2[2] = Ldata[2, snv] + temp
                        first_run = False
                    if i == length - 1 and Ldata[i, snv] != 0:
                        temp2[i] = Ldata[i, snv] + temp
                    # temp_tmp = temp
                    count = 0
                    temp = 0.0
                    # else:
                    #    if Ldata[i,snv]!=0:
                    #        temp2[i]=Ldata[i,snv]+temp_tmp
        for i in range(length):
            stec[i, snv] = temp2[i]

    return stec


def cal_vtec(stec, stnx, stny, stnz, satpoi, e, target_elve, h1, h2, gim_min=0):
    import numpy as np
    f1 = 1575.42e6
    f2 = 1227.60e6
    c = 2.99792458e8
    fat = ((f1 ** 2) * (f2 ** 2) / (f1 ** 2 - f2 ** 2)) / (40.3e16)
    h1 *= 1000
    h2 *= 1000
    Re = 6378.137e3
    rad = 180.0 / np.pi

    time_num = stec.shape[0]
    sat_num = stec.shape[1]
    vtec = np.zeros((time_num, sat_num))
    lon = np.zeros((time_num, sat_num))
    lat = np.zeros((time_num, sat_num))
    B = np.zeros((time_num, sat_num))
    abar = np.zeros(3)
    stnpoi2 = np.sqrt(stnx * stnx + stny * stny + stnz * stnz)
    for j in range(sat_num):
        for i in range(time_num):
            B[i, j] = 0
            if e[i, j] >= target_elve:
                abar[0] = satpoi[i, j, 0] - stnx
                abar[1] = satpoi[i, j, 1] - stny
                abar[2] = satpoi[i, j, 2] - stnz
                abar2 = np.sqrt(abar[0] * abar[0] + abar[1] * abar[1] + abar[2] * abar[2])
                efat = 1 / (h1 - h2) * ( np.sqrt((Re * Re) * np.sin(e[i, j] / rad) ** 2 - Re * Re + (Re + h1) ** 2) - np.sqrt( (Re * Re) * np.sin(e[i, j] / rad) ** 2 - Re * Re + (Re + h2) ** 2))
                theta = np.arccos(stnpoi2 * np.cos(e[i, j] / rad) / (stnpoi2 + (h1 + h2) / 2)) - e[i, j] / rad
                rin = np.sqrt(stnx * stnx + stny * stny + stnz * stnz + (stnpoi2 + (h1 + h2) / 2) ** 2 - 2 * stnpoi2 * (
                stnpoi2 + (h1 + h2) / 2) * np.cos(theta))
                lon[i, j] = np.arctan2(abar2 * stny + abar[1] * rin, abar2 * stnx + abar[0] * rin) * rad
                lat[i, j] = np.arcsin((stnz + abar[2] * rin / abar2) / (stnpoi2 + (h1 + h2) / 2)) * rad
                vtec[i, j] = stec[i, j] / efat
                if stec[i, j] != 0: B[i, j] = (gim_min - vtec[i, j]) / (fat / efat)
            else:
                e[i, j] = 0

    return vtec, lon, lat, B


def cal_mintec(vtec, slon, gim_min_time):
    time_num = vtec.shape[0]
    sat_num = vtec.shape[1]

    vdata = 0.0
    count = 0.0
    minTEC = -999
    for ttt in range(time_num):
        for snv in range(sat_num):
            if vtec[ttt, snv] != 0:
                tim = (ttt + 1) / 120.0 + slon / 15.0
                if tim >= 24: tim = tim - 24
                if tim < 0: tim = tim + 24
                if tim >= (gim_min_time - 0.25) and tim < (gim_min_time + 0.25):
                    count += 1
                    vdata = vdata + vtec[ttt, snv]

    if count > 0:
        minTEC = vdata / count

    return minTEC


