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

                rdline = fid.readline()

        self.total_run = self.ed_doy - self.st_doy + 1


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


class TransCondinate:
    def __init__(self):
        pass

    def xyz2g(self, x, y, z):
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