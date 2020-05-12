import numpy as np
import numpy.ma as ma
from .fortranio import *
from .radar_tools import radar_tools as rt
import datetime as dt
import time
import struct




Re = 6370.e3 # Earth radius in (m)


def pawr_read(filename, endian='',minref=0.0):
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
    data['radar_lon'] = buf[0]
    data['radar_lat'] = buf[1]
    data['radar_alt'] = buf[2]
    data['beam_wid_h'] = buf[3]
    data['beam_wid_v'] = buf[4]
    data['beam_wid_r'] = buf[5]
    data['lambda'] = buf[6]
    data['undef'] = buf[7]
    data['minref']= minref

    buf = np.zeros(4, dtype=dtype_int)
    fort_seq_read(f, buf)
    data['na'] = buf[0]
    data['nr'] = buf[1]
    data['ne'] = buf[2]
    data['nvar'] = buf[3]

    data['azim'] = np.zeros(data['na'], dtype=dtype_real)
    fort_seq_read(f, data['azim'])

    data['radi'] = np.zeros(data['nr'], dtype=dtype_real)
    fort_seq_read(f, data['radi'])

    data['elev'] = np.zeros(data['ne'], dtype=dtype_real)
    fort_seq_read(f, data['elev'])

    buf = np.zeros(1, dtype=dtype_real)
    fort_seq_read(f, buf)
    data['attn_fac'] = buf[0]

    buf = np.zeros((data['ne'], data['nr'], data['na']), dtype=dtype_real)
    for ie in range(data['ne']):
        fort_seq_read(f, buf[ie])
    data['ref'] = ma.masked_values(buf, data['undef'])

    for ie in range(data['ne']):
        fort_seq_read(f, buf[ie])
    data['rv'] = ma.masked_values(buf, data['undef'])

    for ie in range(data['ne']):
        fort_seq_read(f, buf[ie])
    data['qc'] = ma.masked_values(buf, data['undef'])

    for ie in range(data['ne']):
        fort_seq_read(f, buf[ie])
    data['attn'] = ma.masked_values(buf, data['undef'])

    f.close()

    mask_ = np.logical_and( data['qc'] == 0.0 , data['ref'] < minref )
    data['ref'][mask_]=minref



    print("Radar data '{:s}' was read in {:.3f} seconds".format(filename, time.time() - t0))

    t0 = time.time()

    [data['z'] , data['lon_gate'] , data['lat_gate'] , data['x'] , data['y'] , data['local_elevation'], data['distance_to_radar'] ]=rt.radar_georeference( 
                rrange=data['radi']       ,
                relev=data['elev']        ,
                razim=data['azim']        ,
                rlon0=data['radar_lon']   ,
                rlat0=data['radar_lat']   ,
                rz0=data['radar_alt']     ,
                ne=data['ne']             ,
                nr=data['nr']             ,
                na=data['na']               )

    #data = antenna_to_cartesian( data )

    #data = cartesian_to_geographic_aeqd( data, R=6370997.)

    print("Radar data '{:s}' was georeferenced in {:.3f} seconds".format(filename, time.time() - t0))

    return data


def dist_ll(lon1, lat1, lon2, lat2):
    return Re * np.arccos(np.sin(np.deg2rad(lat1)) * np.sin(np.deg2rad(lat2)) + 
                          np.cos(np.deg2rad(lat1)) * np.cos(np.deg2rad(lat2)) * np.cos(np.deg2rad(lon2 - lon1)))


def az_ll(lon1, lat1, lon2, lat2):
    return np.rad2deg(np.arctan2(np.sin(np.deg2rad(lon2 - lon1)) * np.cos(np.deg2rad(lat2)),
                                 np.cos(np.deg2rad(lat1)) * np.sin(np.deg2rad(lat2)) - np.sin(np.deg2rad(lat1)) * np.cos(np.deg2rad(lat2)) * np.cos(np.deg2rad(lon2 - lon1))))


def ll_arc_distance(lon0, lat0, arc_dist, az):
    if arc_dist == 0.:
        lon = lon0
        lat = lat0
    else:
        cdist = np.cos(arc_dist / Re)
        sdist = np.sin(arc_dist / Re)
        sinll1 = np.sin(np.deg2rad(lat0))
        cosll1 = np.cos(np.deg2rad(lat0))
        lat = np.rad2deg(np.arcsin(sinll1 * cdist + cosll1 * sdist * np.cos(np.deg2rad(az))))
        lon = lon0 + np.rad2deg(np.arctan2(sdist * np.sin(np.deg2rad(az)), 
                                           cosll1 * cdist - sinll1 * sdist * np.cos(np.deg2rad(az))))
    return lon, lat


def radar_georeference(data, lon=None, lat=None, radi_h=None):
    Ns = 1.21
    ke = 4. / 3.

    data['r_earth'] = Re

    data['symhgt'] = np.zeros((data['ne'], data['nr']), dtype='f4')

    for ie in range(data['ne']):
        for ir in range(data['nr']):
            # Two possibilities the second probably more accurate than the
            # first one, but both assume a standard constant value for the
            # refractivity index.
            #data['symhgt'][ie,ir] = data['radar_alt'] \
            #                   + data['radi'][ir] * np.sin(np.deg2rad(data['elev'][ie])) \
            #                   + (data['radi'][ir] ** 2) / (2. * Ns * Re)
            data['symhgt'][ie,ir] = data['radar_alt'] \
                               + np.sqrt(data['radi'][ir] ** 2 + (ke * Re) ** 2 +
                                         2. * data['radi'][ir] * ke * Re * np.sin(np.deg2rad(data['elev'][ie]))) \
                               - ke * Re

    data['radi_h'] = np.zeros((data['ne'], data['nr']), dtype='f4')
    data['lon'] = np.zeros((data['ne'], data['nr'], data['na']), dtype='f4')
    data['lat'] = np.zeros((data['ne'], data['nr'], data['na']), dtype='f4')
    data['hgt'] = np.zeros((data['ne'], data['nr'], data['na']), dtype='f4')

    #if (lon is None) or (lat is None) or (radi_h is None):
    for ie in range(data['ne']):

         print('ie =', ie)

         # The curvature of the radar beam is not taken into account.
         data['radi_h'][ie,:] = ke * Re * np.arcsin(data['radi'] * np.cos(np.deg2rad(data['elev'][ie])) / (ke * Re)) # Radar horizontal range

         for ir in range(data['nr']):
             for ia in range(data['na']):
                 data['lon'][ie,ir,ia], data['lat'][ie,ir,ia] = \
                     ll_arc_distance(data['radar_lon'], data['radar_lat'], data['radi_h'][ie,ir], data['azim'][ia])

#        np.save('radi_h.npy', data['radi_h'])
#        np.save('lon.npy', data['lon'])
#        np.save('lat.npy', data['lat'])

    #else:
    #    data['radi_h'] = np.load(radi_h)
    #    data['lon'] = np.load(lon)
    #    data['lat'] = np.load(lat)


    for ia in range(data['na']):
        data['hgt'][:,:,ia] = data['symhgt']

    return True


def antenna_to_cartesian(data):
    """
    Return Cartesian coordinates from antenna coordinates.
    Parameters
    ----------
    ranges : array
        Distances to the center of the radar gates (bins) in kilometers.
    azimuths : array
        Azimuth angle of the radar in degrees.
    elevations : array
        Elevation angle of the radar in degrees.
    Returns
    -------
    x, y, z : array
        Cartesian coordinates in meters from the radar.
    Notes
    -----
    The calculation for Cartesian coordinate is adapted from equations
    2.28(b) and 2.28(c) of Doviak and Zrnic [1]_ assuming a
    standard atmosphere (4/3 Earth's radius model).
    .. math::
        z = \\sqrt{r^2+R^2+2*r*R*sin(\\theta_e)} - R
        s = R * arcsin(\\frac{r*cos(\\theta_e)}{R+z})
        x = s * sin(\\theta_a)
        y = s * cos(\\theta_a)
    Where r is the distance from the radar to the center of the gate,
    :math:`\\theta_a` is the azimuth angle, :math:`\\theta_e` is the
    elevation angle, s is the arc length, and R is the effective radius
    of the earth, taken to be 4/3 the mean radius of earth (6371 km).
    References
    ----------
    .. [1] Doviak and Zrnic, Doppler Radar and Weather Observations, Second
        Edition, 1993, p. 21.
    """
    [ elevations , r , azimuths ] = np.meshgrid( data['elev'] , data['radi'] , data['azim'] , indexing='ij' )

    theta_e = elevations * np.pi / 180.0    # elevation angle in radians.
    theta_a = azimuths * np.pi / 180.0      # azimuth angle in radians.
    R = 6371.0 * 1000.0 * 4.0 / 3.0         # effective radius of earth in meters.

    data['z'] = (r ** 2 + R ** 2 + 2.0 * r * R * np.sin(theta_e)) ** 0.5 - R
    s = R * np.arcsin(r * np.cos(theta_e) / (R + data['z']))  # arc length in m.
    data['x'] = s * np.sin(theta_a)
    data['y'] = s * np.cos(theta_a)

    return data


def cartesian_to_geographic_aeqd(data, R=6370997.) :
    """
    Azimuthal equidistant Cartesian to geographic coordinate transform.
    Transform a set of Cartesian/Cartographic coordinates (x, y) to
    geographic coordinate system (lat, lon) using a azimuthal equidistant
    map projection [1]_.
    .. math::
        lat = \\arcsin(\\cos(c) * \\sin(lat_0) +
                       (y * \\sin(c) * \\cos(lat_0) / \\rho))
        lon = lon_0 + \\arctan2(
            x * \\sin(c),
            \\rho * \\cos(lat_0) * \\cos(c) - y * \\sin(lat_0) * \\sin(c))
        \\rho = \\sqrt(x^2 + y^2)
        c = \\rho / R
    Where x, y are the Cartesian position from the center of projection;
    lat, lon the corresponding latitude and longitude; lat_0, lon_0 are the
    latitude and longitude of the center of the projection; R is the radius of
    the earth (defaults to ~6371 km). lon is adjusted to be between -180 and
    180.
    Parameters
    ----------
    x, y : array-like
        Cartesian coordinates in the same units as R, typically meters.
    lon_0, lat_0 : float
        Longitude and latitude, in degrees, of the center of the projection.
    R : float, optional
        Earth radius in the same units as x and y. The default value is in
        units of meters.
    Returns
    -------
    lon, lat : array
        Longitude and latitude of Cartesian coordinates in degrees.
    References
    ----------
    .. [1] Snyder, J. P. Map Projections--A Working Manual. U. S. Geological
        Survey Professional Paper 1395, 1987, pp. 191-202.
    """
    x = np.atleast_1d(np.asarray(data['x']))
    y = np.atleast_1d(np.asarray(data['y']))

    lat_0_rad = np.deg2rad(data['radar_lat'])
    lon_0_rad = np.deg2rad(data['radar_lon'])

    rho = np.sqrt(x*x + y*y)
    c = rho / R

    #with warnings.catch_warnings():
    # division by zero may occur here but is properly addressed below so
    # the warnings can be ignored
    #warnings.simplefilter("ignore", RuntimeWarning)
    lat_rad = np.arcsin(np.cos(c) * np.sin(lat_0_rad) +
                            y * np.sin(c) * np.cos(lat_0_rad) / rho)
    data['lat_gate'] = np.rad2deg(lat_rad)
    # fix cases where the distance from the center of the projection is zero
    data['lat_gate'][rho == 0] = data['radar_lat']

    x1 = x * np.sin(c)
    x2 = rho*np.cos(lat_0_rad)*np.cos(c) - y*np.sin(lat_0_rad)*np.sin(c)
    lon_rad = lon_0_rad + np.arctan2(x1, x2)
    data['lon_gate'] = np.rad2deg(lon_rad)
    # Longitudes should be from -180 to 180 degrees
    data['lon_gate'][data['lon_gate'] > 180] -= 360.
    data['lon_gate'][data['lon_gate'] < -180] += 360.

    return data




