import numpy as np
import numpy.ma as ma
from .fortranio import *
import datetime as dt
import time
import pyart

def pawr_read(filename, endian=''):
    dtype_real = endian + 'f4'
    dtype_int = endian + 'i4'

    t0 = time.time()
    data = {}

    f = open(filename, 'rb')

    buf = np.zeros(6, dtype=dtype_real)
    fort_seq_read(f, buf)
    try:
        data['time'] = dt.datetime(*buf)
    except ValueError:
        data['time'] = None

    buf = np.zeros(8, dtype=dtype_real)
    fort_seq_read(f, buf)
    radar_lon  = buf[0]
    radar_lat  = buf[1]
    radar_alt  = buf[2]
    beam_wid_h = buf[3]
    beam_wid_v = buf[4]
    beam_wid_r = buf[5]
    lambda     = buf[6]
    undef      = buf[7]

    buf = np.zeros(4, dtype=dtype_int)
    fort_seq_read(f, buf)
    na = buf[0]
    nr = buf[1]
    ne = buf[2]
    nvar = buf[3]

    radar = pyart.testing.make_empty_ppi_radar( nr , na , ne )
    radar.latitude['data'] = radar_lon
    radar.longitude['data'] = radar_lat
    radar.altitude['data'] = radar_alt
    radar.sweep_number['data'] = ne
    
    radar.ngates=nr
    radar.nrays=na*ne

    azim = np.zeros(na, dtype=dtype_real)
    fort_seq_read(f, azim)

    radi = np.zeros(nr, dtype=dtype_real)
    fort_seq_read(f, radi)
    radar.range['data'] = radi

    elev = np.zeros(ne, dtype=dtype_real)
    fort_seq_read(f, elev)
    radar.fixed_angle['data']=elev

    buf = np.zeros(1, dtype=dtype_real)
    fort_seq_read(f, buf)
    attn_fac = buf[0]

    radar.azimuth['data']=np.zeros( na * ne )
    radar.elevation['data']=np.zeros( na * ne )
    radar.time['data']=np.zeros( na * ne )


    buf = np.zeros((ne, nr, na), dtype=dtype_real)
    for ie in range(ne):
        fort_seq_read(f, buf[ie])
    ref = buf

    for ie in range(ne):
        fort_seq_read(f, ie)
    wind = buf

    f.close()

    #Save variables into radar object
    ref_dict = pyart.config.get_metadata('reflectivity')
    ref_dict['data'] = np.zeros( na * ne , nr )
    rv_dict = pyart.config.get_metadata('velocity')
    rv_dict['data'] = np.zeros( na * ne , nr )

    iray=0
    position=0
    for jj in range( ne ) :
        for ii in range( na ) :
            ref_dict['data'][iray,:]=ref[jj,:,ii]
            rv_dict['data'][iray,:]=rv[jj,:,ii]
            radar.azimuth['data'][iray] = azim[ii]
            radar.elevation['data'][iray] = elev[ii]
            iray=iray+1

    radar.fields = {'rv': rv_dict , 'ref':ref_dict}

    radar.ray_angle_res=dict()
    radar.ray_angle_res['data'] = 360.0 / na.astype(float)


    print("Radar data '{:s}' was read in {:.3f} seconds".format(filename, time.time() - t0))
    return radar


