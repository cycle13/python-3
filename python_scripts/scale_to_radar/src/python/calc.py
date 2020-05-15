import sys
import numpy as np
import numpy.ma as ma
import numpy.testing as npt
from .calc_for import calc_for
from mpl_toolkits.basemap import Basemap


           
__all__ = ['set_bmap', 'calc_rotate_winds',
           'calc_height', 'interp_z', 'interp_p', 'calc_destagger', 'calc_destagger_uvw', 'calc_qhydro', 'calc_pt',
           'calc_ref', 'extrap_z_t0', 'extrap_z_pt', 'extrap_p_zt', 'calc_slp', 'calc_rhosfc_psfc']


missingv = 1.e-33
rsphere = 6.37122e6


def set_bmap(sio, proj, resolution=None, rtol=1.e-6):
    """
    Set map projection

    XXXXXX
    """
    llcrnrlon = sio.lon[0, 0]
    llcrnrlat = sio.lat[0, 0]
    urcrnrlon = sio.lon[-1, -1]
    urcrnrlat = sio.lat[-1, -1]
    if proj['type'] == 'LC':
        bmap = Basemap(projection='lcc', lat_1=proj['LC_lat1'], lat_2=proj['LC_lat2'], lon_0=proj['basepoint_lon'],
                       llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                       rsphere=rsphere, resolution=resolution)
    elif proj['type'] == 'MER':
        bmap = Basemap(projection='merc', lon_0=proj['basepoint_lon'],lat_0=proj['basepoint_lat'],lat_ts=proj['basepoint_lat'],
                       llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                       rsphere=rsphere, resolution=resolution)
    else:
        raise ValueError('[Error] Unsupport map projection.')

    # verify projection setting
    x, y = bmap(sio.lon, sio.lat)
    x += sio.dimdef['coor_g']['x'][0]
    y += sio.dimdef['coor_g']['y'][0]
    for iy in range(sio.dimdef['len_g']['y']):
        npt.assert_allclose(x[iy,:], sio.dimdef['coor_g']['x'], rtol=rtol, err_msg='[Error] Incorrect projection settings.')
    for ix in range(sio.dimdef['len_g']['x']):
        npt.assert_allclose(y[:,ix], sio.dimdef['coor_g']['y'], rtol=rtol, err_msg='[Error] Incorrect projection settings.')

    return bmap


def calc_rotate_winds(sio, bmap, u=None, v=None, t=None, dryrun=False):
    """
    Calculate the rotation of u, v winds

    Parameters
    ----------
    sio : <scale.io.ScaleIO> class
        Split SCALE I/O class
    bmap : <mpl_toolkits.basemap.Basemap> class
        Basemap class
    u : 3-D ndarray, optional
        Grid wind along i direction (m/s). Read from files if not given
    v : 3-D ndarray, optional
        Grid wind along j direction (m/s). Read from files if not given
    t : int or <datetime.datetime> class or None, optional
        Time to read
        * None -- all times (defalut)
    dryrun : bool, optional
        * True -- dry run mode, only reading necessary data using 'sio' class
        * False -- do real computation (default)

    Returns
    -------
    u_rot : 3-D ndarray
        Rotated u-wind (m/s)
    v_rot : 3-D ndarray
        Rotated v-wind (m/s)
    """
    u_ = (sio.readvar('U') if u is None else u)
    v_ = (sio.readvar('V') if v is None else v)
    if dryrun:
        return None, None

    if sio.bufsize == 0:
        lon = np.copy(sio.lon)
        lat = np.copy(sio.lat)
    else:
        lon = np.copy(sio.lon[sio.bufsize:-sio.bufsize, sio.bufsize:-sio.bufsize])
        lat = np.copy(sio.lat[sio.bufsize:-sio.bufsize, sio.bufsize:-sio.bufsize])
    ny, nx = lon.shape

    if len(u_.shape) == 3:
        nz = u_.shape[0]
        lon = np.repeat(lon[np.newaxis,:,:], nz, axis=0)
        lat = np.repeat(lat[np.newaxis,:,:], nz, axis=0)

    u_rot, v_rot = bmap.rotate_vector(u_, v_, lon, lat)

    mag_square = u_ * u_ + v_ * v_
    tmpcos = (u_rot * u_ + v_rot * v_) / mag_square
    tmpsin = (u_rot * v_ - v_rot * u_) / mag_square
    u_rot = (tmpcos * u_ - tmpsin * v_).astype(u_.dtype)
    v_rot = (tmpsin * u_ + tmpcos * v_).astype(u_.dtype)

    return u_rot, v_rot


def calc_height(sio, topo=None, dryrun=False):
    """
    Calculate the 3-D height

    Parameters
    ----------
    sio : <scale.io.ScaleIO> class
        Split SCALE I/O class
    topo : 2-D ndarray, optional
        Surface height (m). Read from files if not given
    dryrun : bool, optional
        * True -- dry run mode, only reading necessary data using 'sio' class
        * False -- do real computation (default)

    Returns
    -------
    height : 3-D ndarray
        Height in full levels (m)
    height_h : 3-D ndarray
        Height in half levels (m)
    """
    topo_ = (sio.readvar('TOPO') if topo is None else topo)
    if dryrun:
        return None, None

    height = np.zeros((sio.dimdef['len']['z'][0], topo_.shape[0], topo_.shape[1]), dtype=topo_.dtype)
    height_h = np.zeros((sio.dimdef['len']['zh'][0], topo_.shape[0], topo_.shape[1]), dtype=topo_.dtype)

    for k in range(len(sio.z)):
        height[k,:,:] = topo_ + (sio.zh[-1] - topo_) / sio.zh[-1] * sio.z[k]
    for k in range(len(sio.zh)):
        height_h[k,:,:] = topo_ + (sio.zh[-1] - topo_) / sio.zh[-1] * sio.zh[k]

    return height, height_h

def calc_destagger(var, axis=0, first_grd=False):
    """
    Calculate full-level values from the half-level (staggered-grid) values

    Parameters
    ----------
    var : ndarray
        Input variable
    axis : int, optional
        Staggered axis. Default: 0
    first_grd : bool, optional
        * True -- Addtional first-row grids are provided for interpolation
        * False -- No additional first-row grid (default)

    Returns
    -------
    varout : ndarray
        Destaggered variable
    """
    if axis < 0 or axis >= var.ndim:
        raise ValueError("'axis' is invalid. It must be within [0, var.ndim).")
    varshape = list(var.shape)
    if first_grd:
        varshape[axis] -= 1
    if type(var) == ma.MaskedArray:
        varout = ma.masked_all(varshape, dtype=var.dtype)
        varout.fill_value = var.fill_value
    else:
        varout = np.empty(varshape, dtype=var.dtype)

    slice_obj_1 = [slice(None)] * var.ndim
    slice_obj_2 = [slice(None)] * var.ndim
    if first_grd:
        slice_obj_1[axis] = slice(0, varshape[axis])
        slice_obj_2[axis] = slice(1, varshape[axis]+1)
        varout = 0.5 * (var[slice_obj_1] + var[slice_obj_2])
    else:
        slice_obj_1[axis] = 0
        varout[slice_obj_1] = var[slice_obj_1]
        slice_obj_1[axis] = slice(0, varshape[axis]-1)
        slice_obj_2[axis] = slice(1, varshape[axis])
        varout[slice_obj_2] = 0.5 * (var[slice_obj_1] + var[slice_obj_2])

    return varout


def calc_destagger_uvw(sio, rho=None, momx=None, momy=None, momz=None, destagger=True, first_grd=True, t=None, dryrun=False):
    """
    Calculate 3-D u, v, w winds

    Parameters
    ----------
    sio : <scale.io.ScaleIO> class
        Split SCALE I/O class
    rho : 3-D ndarray, optional
        Density (kg/m3). Read from files if not given
    momx : 3-D ndarray, optional
        x-momentum (kg/m2/s). Read from files if not given
    momy : 3-D ndarray, optional
        y-momentum (kg/m2/s). Read from files if not given
    momz : 3-D ndarray, optional
        z-momentum (kg/m2/s). Read from files if not given
    destaggered : bool, optional
        * True -- Destagger momx, momy, momz before calculation (default)
        * False -- Do not need to destagger momx, momy, momz
    first_grd : bool, optional
        * True -- Addtional first-row grids are provided for interpolation (default)
        * False -- No additional first-row grid
    t : int or <datetime.datetime> class or None, optional
        Time to read
        * None -- all times (defalut)
    dryrun : bool, optional
        * True -- dry run mode, only reading necessary data using 'sio' class
        * False -- do real computation (default)

    Returns
    -------
    u : 3-D ndarray
        Destaggered u-wind (m/s)
    v : 3-D ndarray
        Destaggered v-wind (m/s)
    w : 3-D ndarray
        Destaggered w-wind (m/s)
    momx : 3-D ndarray
        Destaggered x-momentum (kg/m2/s)
    momy : 3-D ndarray
        Destaggered y-momentum (kg/m2/s)
    momz : 3-D ndarray
        Destaggered z-momentum (kg/m2/s)
    """
    rho_ = (sio.readvar('DENS', t=t) if rho is None else rho)
    momz_ = (sio.readvar('MOMZ', t=t) if momz is None else momz)
    if momx is None:
        if first_grd:
            momx_ = sio.readvar('MOMX', t=t, bufsize=0)[:,sio.bufsize:-sio.bufsize,sio.bufsize-1:-sio.bufsize]
        else:
            momx_ = sio.readvar('MOMX', t=t)
    else:
        momx_ = momx
    if momy is None:
        if first_grd:
            momy_ = sio.readvar('MOMY', t=t, bufsize=0)[:,sio.bufsize-1:-sio.bufsize,sio.bufsize:-sio.bufsize]
        else:
            momy_ = sio.readvar('MOMY', t=t)
    else:
        momy_ = momy
    if dryrun:
        return None, None, None, None, None, None

    if destagger:
#        print(' --- destagger momx')
        momx_ = calc_destagger(momx_, axis=2, first_grd=first_grd)
#        print(' --- destagger momy')
        momy_ = calc_destagger(momy_, axis=1, first_grd=first_grd)
#        print(' --- destagger momz')
        momz_ = calc_destagger(momz_, axis=0, first_grd=False)
    momz_[0,:,:] = 0.

#    print(2)

    u = momx_ / rho_
    v = momy_ / rho_
    w = momz_ / rho_

#    print(3)

    return u, v, w, momx_, momy_, momz_


def calc_ref_rv(sio, bmap , radar , topo , min_dbz=-20., z=None ,  u=None , v=None , w=None , rho=None, qr=None, qs=None, qg=None, t=None, dryrun=False , use_wt=True ):
    """
    Calculate radar reflectivity

    Parameters
    ----------
    sio : <scale.io.ScaleIO> class
        Split SCALE I/O class
    min_dbz : float
        Minimum value of the reflectivity
    rho : 3-D ndarray, optional
        Density (kg/m3). Read from files if not given
    qr : 3-D ndarray, optional
        Rain water mixing ratio (kg/kg). Read from files if not given
    qs : 3-D ndarray, optional
        Snow mixing ratio (kg/kg). Read from files if not given
    qg : 3-D ndarray, optional
        Graupel mixing ratio (kg/kg). Read from files if not given
    t : int or <datetime.datetime> class or None, optional
        Time to read
        * None -- all times (defalut)
    dryrun : bool, optional
        * True -- dry run mode, only reading necessary data using 'sio' class
        * False -- do real computation (default)

    Returns
    -------
    dbz : 3-D ndarray
        Radar reflectivity (dBZ)
    max_dbz : 2-D ndarray
        Maximum radar reflectivity in vertical (dBZ)
    """
    import math


    #Define constants 
    min_q_ = 1.0e-10
    d2r = np.pi / 180.0
    r2d = 1.0 / d2r

    nor_=8.0e-2      #[cm^-4]
    nos_=3.0e-2      #[cm^-4]
    nog_=4.0e-4      #[cm^-4]
    ror_=1.0         #[g/cm3]
    ros_=0.1         #[g/cm3]
    rog_=0.917       #[g/cm3]
    gg_=9.81         #Gravity
    re_=6371.3e3     #Earth radious


    a_=2115e0   #[cm**1-b / s]
    b_=0.8e0
    c_=152.93e0 #[cm**1-b / s]
    d_=0.25e0
    Cd_=0.6e0

    #Read variables
    rho_ = (sio.readvar('DENS', t=t) if rho is None else rho)
    nz , ny, nx = rho_.shape
    qr_ = (sio.readvar('QR', t=t) if qr is None else qr)
    qs_ = (sio.readvar('QS', t=t) if qs is None else qs)
    qg_ = (sio.readvar('QG', t=t) if qg is None else qg)

    if ( u is None or v is None or w is None ) :
       [u_,v_,w_,momx_,momy_,momz_]=calc_destagger_uvw(sio, t=t , first_grd=False )
       [ u_ , v_ ]=calc_rotate_winds(sio, bmap, u=u_, v=v_, t=t)
    else :
       u_ = np.copy(u)
       v_ = np.copy(v)
       w_ = np.copy(w)

    if ( z is None ) :
       [ z_ , z_h_ ] = calc_height( sio , topo=topo )

    if sio.bufsize == 0:
        lon_ = np.copy(sio.lon)
        lat_ = np.copy(sio.lat)
    else:
        lon_ = np.copy(sio.lon[sio.bufsize:-sio.bufsize, sio.bufsize:-sio.bufsize])
        lat_ = np.copy(sio.lat[sio.bufsize:-sio.bufsize, sio.bufsize:-sio.bufsize])

    if dryrun:
        return None, None

    qr_[qr_ < 1.e-10] = min_q_
    qs_[qs_ < 1.e-10] = min_q_
    qg_[qg_ < 1.e-10] = min_q_

    rofactor_= ( 1.0 / rho_ ) ** 0.5


    #Compute reflectivity 
    refr_= 2.53e4 * (rho_ * qr_ * 1.0e3) ** 1.84
    refs_= 3.48e3 * (rho_ * qs_ * 1.0e3) ** 1.66
    refg_= 8.18e4 * (rho_ * qg_ * 1.0e3) ** 1.50

    dbz = 10. * np.log10(refr_ + refs_ + refg_ )
    dbz[dbz < min_dbz] = min_dbz
    max_dbz = ma.max(dbz, axis=0)

    #Compute terminal velocity
    if use_wt :
       #Rain terminal velocity
       gamma_=math.gamma( 4.0 + b_ )
       mask_= qr_ > min_q_
       wr_ = np.zeros( qr_.shape )
       lr_ = np.zeros( qr_.shape )
       lr_[mask_]= np.power( np.pi * ror_ * nor_ /  ( 1.0e-3 * rho_[mask_] * qr_[mask_] ) , 0.25) 

       wr_[mask_]= a_ * gamma_ / ( 6.0 * ( np.power( lr_[mask_] , b_ ) ) )
       wr_= 1.0e-2 * wr_ * rofactor_

       #Snow terminal velocity
       gamma_=math.gamma( 4.0 + d_ )
       mask_= qs_ > min_q_
       ws_ = np.zeros( qs_.shape )
       ls_ = np.zeros( qs_.shape )
       ls_[mask_]= ( np.pi * ros_ * nos_ / ( 1.0e-3 * rho_[mask_] * qs_[mask_] ) ) ** 0.25
       ws_[mask_]= c_ * gamma_ / ( 6.0 * ( ls_[mask_] ** d_ ) )
       ws_= 1.0e-2 * ws_ * rofactor_

       #Graupel terminal velocity
       gamma_=math.gamma( 4.5 )
       mask_= qg_ > min_q_
       wg_ = np.zeros( qg_.shape )
       lg_ = np.zeros( qg_.shape )
       lg_[mask_]= ( np.pi * rog_ * nog_ / ( 1.0e-3 * rho_[mask_] * qg_[mask_] ) ) ** 0.25
       wg_[mask_]= gamma_ *  np.power( ( 4.0 * gg_ * 100.0 * rog_ )/( 3.0 * Cd_ * 1.0e-3 * rho_[mask_] ) , 0.5 )
       wg_[mask_]= 1.0e-2 * wg_[mask_] / ( 6.0 * np.power( lg_[mask_] , 0.5 ) )

       wt_ = ( wr_ * refr_ + ws_ *  refs_ +  wg_ *  refg_ ) / ( refr_ + refs_ + refg_ )

    #Compute azimuth and elevation angle with respect to radar
    dlon_=lon_-radar['radar_lon']
    dlat_=lat_-radar['radar_lat']

    mask_ = np.logical_not( np.logical_and( dlon_ == 0.0 , dlat_ == 0.0 ) )
    az_ = np.zeros( dlon_.shape )
    az_[mask_] = np.arctan2(-dlon_[mask_]*np.cos(radar['radar_lat']*d2r)*d2r,-dlat_[mask_]*d2r) + np.pi
    #az_[ az_ < 0.0 ] = 2.0*np.pi + az_[az_ < 0.0]
    sin_az_ = np.sin( az_ )
    cos_az_ = np.cos( az_ )
    sin_az_ = np.repeat(sin_az_[np.newaxis,:,:], nz, axis=0)
    cos_az_ = np.repeat(cos_az_[np.newaxis,:,:], nz, axis=0)

    cosd_ =np.sin( lat_ * d2r )*np.sin( radar['radar_lat'] * d2r ) + np.cos( lat_*d2r )*np.cos( radar['radar_lat'] * d2r )*np.cos( ( lon_ - radar['radar_lon'] )*d2r )
    cosd_[cosd_ > 1.0 ]=1.0
    cosd_[cosd_ <-1.0 ]=-1.0
    dist_ = np.arccos( cosd_ ) * re_

    dist_ = np.repeat(dist_[np.newaxis,:,:], nz, axis=0)
    elev_ = np.arctan2( z_ - radar['radar_alt'] , dist_ )

    #Compute radial velocity
    rv = u_ * np.cos(elev_) * np.sin(az_)
    rv = rv + v_ * np.cos(elev_) * np.cos(az_)
    if( use_wt ) :
       rv = rv + (w_ - wt_)*np.sin(elev_)
    else  :
       rv = rv + (w_)*np.sin(elev_)

    #import matplotlib.pyplot as plt
    #plt.pcolor( lon_ , lat_ , rv[10,:,:]) ; plt.colorbar() ; plt.show()
    return dbz,rv,max_dbz

def radar_int( sio , proj , topo , radar , t=None ) :
    from scipy.interpolate import interpn

    bmap = set_bmap(sio, proj, resolution=None, rtol=1.e-6)

    [ref_,rv_,max_ref_]=calc_ref_rv(sio, bmap , radar , topo ,  min_dbz=-20.0 , t=t )

    if sio.bufsize == 0:
        lon_ = np.copy(sio.lon)
        lat_ = np.copy(sio.lat)
    else :
        lon_ = np.copy(sio.lon[sio.bufsize:-sio.bufsize, sio.bufsize:-sio.bufsize])
        lat_ = np.copy(sio.lat[sio.bufsize:-sio.bufsize, sio.bufsize:-sio.bufsize])
    dx_ = sio.dimdef['coor_g']['x'][sio.bufsize+1] - sio.dimdef['coor_g']['x'][sio.bufsize]
    dy_ = sio.dimdef['coor_g']['y'][sio.bufsize+1] - sio.dimdef['coor_g']['y'][sio.bufsize]
    [nz_ , ny_ , nx_] = ref_.shape


    [ z_ , z_h_ ] = calc_height( sio , topo=topo )

    na_=radar['na']
    nr_=radar['nr']
    ne_=radar['ne']
    tlat_ = np.reshape( radar['lat_gate']  , ne_ * nr_ * na_ ).astype('float64')
    tlon_ = np.reshape( radar['lon_gate'] , ne_ * nr_ * na_ ).astype('float64')
    tz_   = np.reshape( radar['z']  , ne_ * nr_ * na_ ).astype('float64')

    ri_, rj_ = bmap(tlon_,tlat_)
    ri_=ri_/dx_ + 1.0
    rj_=rj_/dy_ + 1.0
    rk_ = calc_for.get_k(z_full=np.transpose(z_,axes=[0,2,1]),nlev=nz_,nlon=nx_,nlat=ny_,ri=ri_+1.0,rj=rj_+1.0,rlev=tz_,nin=np.size(tz_) )
    rk_ = rk_ 

    rijk_interp=np.zeros((na_*nr_*ne_,3))
    rijk_interp[:,0]=rk_
    rijk_interp[:,1]=rj_
    rijk_interp[:,2]=ri_

    xc_ = np.arange(nx_)
    yc_ = np.arange(ny_)
    zc_ = np.arange(nz_)

    undef = radar['ref'].fill_value
    nin=na_*nr_*ne_

    ref_ = 10.0 ** ( ref_ / 10.0 )

    radar['model_ref'] = ( calc_for.itpl_3d_vec(var=np.transpose(ref_,axes=[0,2,1]),nlev=nz_,nlon=nx_,nlat=ny_,ri=ri_,rj=rj_,rk=rk_,fill_value=undef,nin=nin) ).reshape( ne_,nr_,na_ )
    radar['model_rv'] = ( calc_for.itpl_3d_vec(var=np.transpose(rv_,axes=[0,2,1]),nlev=nz_,nlon=nx_,nlat=ny_,ri=ri_,rj=rj_,rk=rk_,fill_value=undef,nin=nin) ).reshape( ne_,nr_,na_ )

    radar['model_ref'] = 10.0 * np.log10( radar['model_ref'] )

    #radar['model_ref'] = interpn((zc_,yc_,xc_), ref_ , rijk_interp, method='linear', bounds_error=False, fill_value = radar['ref'].fill_value ).reshape(ne_,nr_,na_)
    #radar['model_rv']  = interpn((zc_,yc_,xc_), rv_  , rijk_interp, method='linear', bounds_error=False, fill_value = radar['rv'].fill_value  ).reshape(ne_,nr_,na_)

    return radar

def radar_regrid( radar , grid )  :

    nvar = 4 #We will process only reflectivity and wind for the model and for the radar.

    radar_grid=dict()
    radar_grid['grid']=grid

    [ne,nr,na]=np.shape( radar['ref'] )
    nin = ne*nr*na
    x=radar['x'].reshape( nin )
    y=radar['y'].reshape( nin )
    z=radar['z'].reshape( nin )

    data=np.zeros( ( nin , nvar ) )
    data[:,0] = (10.0**( radar['ref'].data / 10.0 ) ).reshape( nin )
    data[:,1] = (radar['rv'].data).reshape( nin )
    data[:,2] = (10.0**( radar['model_ref'].data / 10.0 ) ).reshape( nin )
    data[:,3] = (radar['model_rv'].data).reshape( nin )

    undef = radar['ref'].fill_value
    
    [data_ave , data_max , data_min , data_std , data_n ] = calc_for.com_interp_boxavereg(xini=grid['xini'],dx=grid['dx'],nx=grid['nx']
                                                                                                  ,yini=grid['yini'],dy=grid['dy'],ny=grid['ny']
                                                                                                  ,zini=grid['zini'],dz=grid['dz'],nz=grid['nz']
                                                                                                  ,nvar=nvar,xin=x,yin=y,zin=z,datain=data
                                                                                                  ,undef=undef,nin=nin )

    radar_grid['data_ave']=np.ma.masked_values(data_ave, undef )
    radar_grid['data_max']=np.ma.masked_values(data_max, undef )
    radar_grid['data_min']=np.ma.masked_values(data_min, undef )
    radar_grid['data_n']  =data_n
    radar_grid['var_list']=['ref','rv','model_ref','model_rv']

    radar_grid['data_ave'][:,:,:,0]=10.0 * np.log10( radar_grid['data_ave'][:,:,:,0] )
    radar_grid['data_ave'][:,:,:,2]=10.0 * np.log10( radar_grid['data_ave'][:,:,:,2] )
    radar_grid['data_max'][:,:,:,0]=10.0 * np.log10( radar_grid['data_max'][:,:,:,0] )
    radar_grid['data_max'][:,:,:,2]=10.0 * np.log10( radar_grid['data_max'][:,:,:,2] )
    radar_grid['data_min'][:,:,:,0]=10.0 * np.log10( radar_grid['data_min'][:,:,:,0] )
    radar_grid['data_min'][:,:,:,2]=10.0 * np.log10( radar_grid['data_min'][:,:,:,2] )


    return radar_grid




