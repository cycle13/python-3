import numpy as np
import datetime as dt
from letkfobsio import *
import sys

obsdir = '/data15/gylien/PAWR/letkf_v4/20170716/QCED_1KM_v4_v500M_attn0.09'
obsoutdir = '/data15/gylien/PAWR/letkf_v4/20170716/QCED_1KM_v4_v500M_attn0.09_4d_-4_5'

stime = dt.datetime(2017, 7, 16,  3, 10, 0)
etime = dt.datetime(2017, 7, 16, 10, 55, 0)
tintout = dt.timedelta(seconds=300)
tint = dt.timedelta(seconds=30)

obs4d_slot_s = -4
obs4d_slot_e = 5

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
