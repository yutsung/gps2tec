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