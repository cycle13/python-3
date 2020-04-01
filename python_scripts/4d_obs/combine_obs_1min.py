import sys
sys.path.append('../../common_python/scale2gradsv180101/')

import numpy as np
import datetime as dt
from letkfobsio import *
import sys

obsdir = '/home/ra001011/a03471/data/input_data/obs/QCED_1KM_v4_v500M_attn0.01/'
obsoutdir = '/home/ra001011/a03471/data/input_data/obs/QCED_1KM_v4_v500M_attn0.01_4D_1min/'

stime = dt.datetime(2013, 7, 13,  4,  0,  0)
etime = dt.datetime(2013, 7, 13,  6,  0,  0)
tintout = dt.timedelta(seconds=60)
tint = dt.timedelta(seconds=30)

obs4d_slot_s = -1
obs4d_slot_e = 0

time = stime
while time <= etime:
    for islot in range(obs4d_slot_s, obs4d_slot_e+1):
        istime = time + tint * islot
        print(istime.strftime('%Y:%m:%d %H-%M-%S'), '->', time.strftime('%Y:%m:%d %H-%M-%S'))
        sys.stdout.flush()

        filename = "{:s}/radar_{:s}.dat".format(obsdir, istime.strftime('%Y%m%d%H%M%S'))

        f = open(filename, 'rb')
        radar_lon, radar_lat, radar_z, obs = readobs_radar_all_fast(f, endian='>')
        f.close()

        obs_tmp = np.array([tuple(i) + ((tint * islot).total_seconds(),) for i in obs], dtype=np.dtype(obsrecords['scale']))

        if islot == obs4d_slot_s:
           obsout = obs_tmp
        else:
           obsout = np.concatenate((obsout, obs_tmp))

    fo = open("{:s}/radar_{:s}.dat".format(obsoutdir, time.strftime('%Y%m%d%H%M%S')), 'wb')
    writeobs_radar_all(fo, radar_lon, radar_lat, radar_z, obsout, endian='>')
    fo.close()

    time += tintout
