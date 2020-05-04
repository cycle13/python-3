import sys
import calc
sys.stderr.write("hello from cython\n")

import numpy as np
import numpy.ma as ma
import numpy.testing as npt
from mpl_toolkits.basemap import Basemap

           
__all__ = ['set_bmap', 'calc_rotate_winds',
           'calc_height', 'interp_z', 'interp_p', 'calc_destagger', 'calc_destagger_uvw', 'calc_qhydro', 'calc_pt',
           'calc_ref', 'extrap_z_t0', 'extrap_z_pt', 'extrap_p_zt', 'calc_slp', 'calc_rhosfc_psfc']


missingv = 1.e-33
rsphere = 6.37122e6


def radar_int( sio , proj , topo , radar , t=None ) :

    bmap = calc.set_bmap(sio, proj, resolution=None, rtol=1.e-6)

    [dbz_,rv_,max_dbz_]=calc.calc_ref_rv(sio, bmap , radar , topo ,  min_dbz=-20. , t=t )

    if sio.bufsize == 0:
        lon_ = np.copy(sio.lon)
        lat_ = np.copy(sio.lat)
    else :
        lon_ = np.copy(sio.lon[sio.bufsize:-sio.bufsize, sio.bufsize:-sio.bufsize])
        lat_ = np.copy(sio.lat[sio.bufsize:-sio.bufsize, sio.bufsize:-sio.bufsize])

    [ z_ , z_h_ ] = calc_height( sio , topo=topo )

    na=radar.nrays
    nr=radar.ngates
    target_lat_ = np.reshape( radar.gate_latitude  , na * nr )
    target_lon_ = np.reshape( radar.gate_longitude , na * nr )
    target_z_   = np.reshape( radar.gate_altitude  , na * nr )

    radar_ref = vol_int( dbz_ , lon_ , lat_ , z_ , target_lon_ , target_lat_ , target_z_ )
    radar_rv  = vol_int( dbz_ , lon_ , lat_ , z_ , target_lon_ , target_lat_ , target_z_ )

    return radar_ref , radar_vr 


def general_int( var , lon , lat , z , tlon , tlat , tz ) :

    var_int = np.zeros( tlon.shape )


    return var_int




