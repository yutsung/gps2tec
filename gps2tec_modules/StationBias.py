def cal_stnbias_GIM(B, e, slon, gim_min_time):
    import numpy as np

    sat_num  = B.shape[1]
    time_num = B.shape[0]
    data  = 0
    e2    = 0
    count = 0

    for ttt in range(time_num):
        for snv in range(sat_num):
            if B[ttt, snv] !=0 and e[ttt, snv] !=0:
                tim = (ttt+1)/120.0+slon/15.0
                if tim >=24 : tim = tim - 24
                if tim < 0  : tim = tim + 24
                if tim >= (gim_min_time - 0.25) and tim < (gim_min_time + 0.25):
                    count += 1
                    data = data + B[ttt, snv]*(np.sqrt(e[ttt, snv]))
                    e2 = e2 + np.sqrt(e[ttt, snv])
    if count > 0:
        stnbias = data/e2
    else:
        return True

    return stnbias

def cal_stnbias_FITTING(stec, el, ep, satdata,  gpsx, gpsy, gpsz):
    import numpy as np
    from gps2tec_modules import DataProcess
    time_idx, satnum = stec.shape
    clight = 2.99792458e8
    f1 = 1.57542e9
    f2 = 1.2276e9
    ionfac = 40.3082e16
    m2tecu = (f1*f1)*(f2*f2)/(f1*f1-f2*f2)/ionfac
    tecu2ns = 1e9/clight/m2tecu
    rad = 180. /np.pi
    km2m = 1e3
    order = 4

    rx = gpsx
    ry = gpsy
    rz = gpsz

    rlon, rlat, ralt= DataProcess.xyz2g(rx, ry, rz)
    glon, glat, galt= DataProcess.xyz2g(satdata[:, :, 0], satdata[:, :, 1], satdata[:, :, 2])
    if rlon > 180: rlon -= 360
    glon[np.where(glon > 180)] -= 360
    A = azimuth(rlat, rlon, glat, glon)

    t_r = ep*2*np.pi/24
    H, L = Hmatrix(stec*tecu2ns, A, el, t_r, rlon, rlat, order)
    R, resid, rank, s=np.linalg.lstsq(H, L)
    return R[0, 0] * 0.299792458


def azimuth(lat1, lon1, lat2, lon2):
    import  numpy as np
    lat1 = lat1 * np.pi / 180
    lon1 = lon1 * np.pi / 180
    lat2 = lat2 * np.pi / 180
    lon2 = lon2 * np.pi / 180

    az = np.arctan2(np.cos(lat2) * np.sin(lon2-lon1), (np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(lon2 - lon1)))
    az[np.where(lat1 <= -np.pi / 2)] = 0
    az[np.where(lat2 >=  np.pi / 2)] = 0
    az[np.where(lat2 <= -np.pi / 2)] = np.pi
    az[np.where(lat1 >=  np.pi / 2)] = np.pi

    az = np.mod(az, 2 * np.pi)
    az = az * 180 / np.pi

    return az

def Hmatrix(stec, A, el, t_r, rlon, rlat, order):
    import numpy as np
    time_idx, satnum = stec.shape
    dhour = 12
    rearth = 6371.
    h = 450.
    rad = 180. /np.pi
    nsat = satnum

    datanumber = len(np.where(stec!=0)[0])
    A[np.where(stec==0)] = np.nan
    el[np.where(stec==0)] = np.nan
    stec[np.where(stec==0)] = np.nan
    zp = np.arcsin(rearth/(rearth+h)*np.sin((0.9782*(90 - el))*np.pi/180))
    MAPF = np.cos(zp)
    ltr=time_idx
    dltr=ltr/dhour
    t_r = t_r[:,None]*np.ones((time_idx, satnum))
    H = np.zeros((datanumber, (order+1)*(order+1)*12+1))
    L = np.zeros((datanumber,1))
    sbm = rlat*np.ones((time_idx, satnum))/rad
    slm = rlon*np.ones((time_idx, satnum))/rad
    b, s = Get_IPP(el/rad, A/rad, sbm, slm, zp, t_r)
    c = 0
    for i in range(dhour):
        k = np.arange(i*dltr,dltr*(i+1))
        zp1   = zp[k,:].copy()
        MAPF1 = MAPF[k,:].copy()
        stec1 = stec[k,:].copy()
        b1    = b[k,:].copy()
        s1    = s[k,:].copy()
        idx   = ~np.isnan(b1)
        b1    = b1[idx].copy()
        s1    = s1[idx].copy()
        stec1 = stec1[idx].copy()
        zp1   = zp1[idx].copy()
        MAPF1 = MAPF1[idx].copy()
        H[c:c+len(zp1),0]   = -1*MAPF1
        st = (order+1)*(order+1)*i+1
    	ed = (order+1)*(order+1)*(i+1)
        H[c:c+len(zp1),st:ed+1]  = Get_coef(b1, s1, order)
        L[c:c+len(zp1),0] = stec1*MAPF1
        c = c+len(zp1) 
    return H, L

def Get_IPP(E, A, B, L, z, t_r):
    import numpy as np

    t = np.pi/2 - E - z
    b = np.arcsin(np.sin(B)*np.cos(t)+np.cos(B)*np.sin(t)*np.cos(A))
    s = L+np.arcsin(np.sin(t)*np.sin(A)/np.cos(t))
    s = s+t_r-np.pi
    return b, s

def Get_coef(b,s,order):
    import numpy as np
    from scipy import special

    cof_P = np.zeros((len(b), (order+1)*(order+1)))
    ms = np.zeros((len(s),4))

    for i in range(len(s)):
        ms[i,:] = np.linspace(s[i], order*s[i], order)

    i = 0;
    x = np.sin(b)
    for n in range(order+1):
        for m in range(n+1):
            P = special.lpmv(m,n,x)
            if m==0:
                cof_P[:, i] = P*normP(n,m)
            else:
                cof_P[:, i] = P*normP(n,m)*np.cos(ms[:,m-1])
                i=i+1
                cof_P[:, i] = P*normP(n,m)*np.sin(ms[:,m-1]) 
            i=i+1
    return cof_P

def normP(n,m):
    from math import factorial
    import numpy as np

    if m==0:
        return np.sqrt(factorial(n-m)*(2*n+1)/factorial(n+m))
    else:
        return np.sqrt(factorial(n-m)*(4*n+2)/factorial(n+m))
    

