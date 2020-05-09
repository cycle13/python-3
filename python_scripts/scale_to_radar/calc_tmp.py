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
    from scipy.interpolate import interpn

    bmap = calc.set_bmap(sio, proj, resolution=None, rtol=1.e-6)

    [ref_,rv_,max_ref_]=calc.calc_ref_rv(sio, bmap , radar , topo ,  min_dbz=-20. , t=t )

    if sio.bufsize == 0:
        lon_ = np.copy(sio.lon)
        lat_ = np.copy(sio.lat)
    else :
        lon_ = np.copy(sio.lon[sio.bufsize:-sio.bufsize, sio.bufsize:-sio.bufsize])
        lat_ = np.copy(sio.lat[sio.bufsize:-sio.bufsize, sio.bufsize:-sio.bufsize])
    dx_ = sio.dimdef['coor_g']['x'][sio.bufsize+1] - sio.dimdef['coor_g']['x'][sio.bufsize]
    dy_ = sio.dimdef['coor_g']['y'][sio.bufsize+1] - sio.dimdef['coor_g']['y'][sio.bufsize]
    nz_ , nx_ , ny_ = ref_.shape

    [ z_ , z_h_ ] = calc.calc_height( sio , topo=topo )

    na_=radar.nrays
    nr_=radar.ngates
    tlat_ = np.reshape( radar.gate_latitude['data']  , na_ * nr_ )
    tlon_ = np.reshape( radar.gate_longitude['data'] , na_ * nr_ )
    tz_   = np.reshape( radar.gate_altitude['data']  , na_ * nr_ )

    ri_, rj_ = bmap(tlon_,tlat_)
    ri_=ri_/dx_
    rj_=rj_/dy_
    rk_  = get_k( ri_ , rj_ , tz_ , z_ )

    rijk_interp=np.zeros((na_*nr_,3))
    rijk_interp[:,0]=rk_
    rijk_interp[:,1]=ri_
    rijk_interp[:,2]=rj_

    xc_ = np.arange(nx_) 
    yc_ = np.arange(ny_) 
    zc_ = np.arange(nz_) 

    model_ref = interpn((zc_,xc_,yc_), ref_ , rijk_interp, method='linear', bounds_error=False).reshape(na_,nr_)
    model_rv  = interpn((zc_,xc_,yc_), rv_  , rijk_interp, method='linear', bounds_error=False).reshape(na_,nr_)

    return model_ref , model_rv 

def get_k( np.ndarray[double, ndim=1] ri , np.ndarray[double, ndim=1] rj , np.ndarray[double, ndim=1] tz , np.ndarray[f_real, ndim=1] z ) :

    cdef int nk,nx,ny,ii,kk
    nz,nx,ny=np.shape(z)
    cdef np.ndarray[int, ndim=1] i=np.round(ri).astype(int)
    cdef np.ndarray[int, ndim=1] j=np.round(rj).astype(int)
    i[i<0]=0
    i[i>=nx]=nx-1
    j[j<0]=0
    j[j>=ny]=ny-1

    cdef int nk=ri.size
    cdef np.ndarray[double, ndim=1] k=np.zeros(nk)
    cdef np.ndarray[double, ndim=1] tmp_z

    for ii in range(nk) :
        tmp_z = z[:,i[ii],j[ii]]
        if tz[ii] < tmp_z[0] :
            k[ii] = 0
        if tz[ii] > tmp_z[nz-1] :
            k[ii] = nz-1
        if tz[ii] > tmp_z[0] and tz[ii] < tmp_z[nz-1] :
           for kk in range(nz-1) :
               if tz[ii] >= tmp_z[kk] and tz[ii] <= tmp_z[kk+1] :
                   k[ii]=kk + (tz[ii] - tmp_z[kk])/(tmp_z[kk+1]-tmp_z[kk])
                   break

    return k






