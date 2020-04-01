import numpy as np
import numpy.ma as ma
from fortranio import *
import datetime as dt
import time
from scipy.interpolate import interp2d


import sys
sys.path.append('../../common_python/common_functions/')
from common_functions import common_functions as cf


Re = cf.re #6370.e3 # Earth radius in (m)


def radarobs_read(filename, endian=''):

    dtype_real = endian + 'f4'
    dtype_int = endian + 'i4'

    t0 = time.time()
    data = {}

    f = open(filename, 'rb')

    buf = np.zeros(6, dtype=dtype_real)
    fort_read(f, buf)
    try:
        data['time'] = dt.datetime(*buf)
    except ValueError:
        data['time'] = None

    buf = np.zeros(8, dtype=dtype_real)
    fort_read(f, buf)
    data['radar_lon'] = buf[0]
    data['radar_lat'] = buf[1]
    data['radar_alt'] = buf[2]
    data['beam_wid_h'] = buf[3]
    data['beam_wid_v'] = buf[4]
    data['beam_wid_r'] = buf[5]
    data['lambda'] = buf[6]
    data['undef'] = buf[7]

    buf = np.zeros(4, dtype=dtype_int)
    fort_read(f, buf)
    data['na'] = buf[0]
    data['nr'] = buf[1]
    data['ne'] = buf[2]
    data['nvar'] = buf[3]

    data['azim'] = np.zeros(data['na'], dtype=dtype_real)
    fort_read(f, data['azim'])

    data['radi'] = np.zeros(data['nr'], dtype=dtype_real)
    fort_read(f, data['radi'])

    data['elev'] = np.zeros(data['ne'], dtype=dtype_real)
    fort_read(f, data['elev'])

    buf = np.zeros(1, dtype=dtype_real)
    fort_read(f, buf)
    data['attn_fac'] = buf[0]

    buf = np.zeros((data['ne'], data['nr'], data['na']), dtype=dtype_real)
    for ie in range(data['ne']):
        fort_read(f, buf[ie])
    data['ref'] = ma.masked_values(buf, data['undef'])

    for ie in range(data['ne']):
        fort_read(f, buf[ie])
    data['wind'] = ma.masked_values(buf, data['undef'])

    for ie in range(data['ne']):
        fort_read(f, buf[ie])
    data['qc'] = ma.masked_values(buf, data['undef'])

    for ie in range(data['ne']):
        fort_read(f, buf[ie])
    data['attn'] = ma.masked_values(buf, data['undef'])

    f.close()

    print("Radar data '{:s}' was read in {:.3f} seconds".format(filename, time.time() - t0))
    return data

def radarobs_read_preqc(filename, endian=''):

    dtype_real = endian + 'f4'
    dtype_int = endian + 'i4'

    t0 = time.time()
    data = {}

    f = open(filename, 'rb')

    buf = np.zeros(5, dtype=dtype_real)
    fort_read(f, buf,seq=False)
    data['radar_lon'] = buf[0]
    data['radar_lat'] = buf[1]
    data['radar_alt'] = buf[2]
    data['beam_wid_h'] = buf[3]
    data['beam_wid_v'] = buf[4]

    #print( data['radar_lat'] , data['radar_lon'] )

    buf = np.zeros(3, dtype=dtype_int)
    fort_read(f, buf , seq=False)
    data['na'] = buf[0]
    data['nr'] = buf[1]
    data['ne'] = buf[2]
  
    #print( data['na'] , data['nr'] , data['ne'] )

    buf = np.zeros((data['ne'], data['nr'], data['na']), dtype=dtype_real)
    for ie in range(data['ne']):
        fort_read(f, buf[ie],seq=False)
    #Warning assuming that the minimum value is the one corresponding to the 
    #missing values, since the missing values is not the defined in the pre-qc 
    #data format.
    data['data'] = ma.masked_values(buf, np.min( buf ) )

    buf = np.zeros(2, dtype=dtype_int)
    fort_read(f, buf,seq=False)

    tmp = np.zeros( (data['na'],data['ne']) , dtype=dtype_real)
    fort_read(f, tmp , seq=False)
    data['azim']=tmp[:,0]

    tmp = np.zeros( (data['na'],data['ne']) , dtype=dtype_real)
    fort_read(f, tmp ,seq=False)
    data['elev']=tmp[0,:]

    #Warning assuming that the range resolution is 100 meters!!
    #Pre qc files do not contain information about the resolution in range.    
    data['radi']=np.zeros(data['nr'])
    for ir in range(data['nr'])  :
      data['radi'][ir]=(ir)*100


    f.close()

    print("Radar data '{:s}' was read in {:.3f} seconds".format(filename, time.time() - t0))
    return data



def radar_georeference(data, lon=None, lat=None, radi_h=None):

    Ns = 1.21
    ke = 4. / 3.

    data['r_earth'] = Re

    data['symhgt'] = np.zeros((data['ne'], data['nr']), dtype='f4')

    for ie in range(data['ne']):
        for ir in range(data['nr']):
            data['symhgt'][ie,ir] = data['radar_alt'] \
                               + np.sqrt(data['radi'][ir] ** 2 + (ke * Re) ** 2 +
                                         2. * data['radi'][ir] * ke * Re * np.sin(np.deg2rad(data['elev'][ie]))) \
                               - ke * Re

    data['radi_h'] = np.zeros((data['ne'], data['nr']), dtype='f4')
    data['lon'] = np.zeros((data['ne'], data['nr'], data['na']), dtype='f4')
    data['lat'] = np.zeros((data['ne'], data['nr'], data['na']), dtype='f4')
    data['hgt'] = np.zeros((data['ne'], data['nr'], data['na']), dtype='f4')

    if (lon is None) or (lat is None) or (radi_h is None):
        for ie in range(data['ne']):

            # The curvature of the radar beam is not taken into account.
            data['radi_h'][ie] = ke * Re * np.arcsin(data['radi'] * np.cos(np.deg2rad(data['elev'][ie])) / (ke * Re)) # Radar horizontal range
            [tmp_a,tmp_r] = np.meshgrid(data['azim'],data['radi_h'][ie,:])
            [data['lon'][ie,:,:],data['lat'][ie,:,:]]=cf.com_ra_to_ll(data['radar_lon'],data['radar_lat'],tmp_r,tmp_a,data['nr'],data['na'])
            #for ir in range(data['nr']):
            #    for ia in range(data['na']):
            #        data['lon'][ie,ir,ia], data['lat'][ie,ir,ia] = \
            #            cf.com_ll_arc_distance_sngl(data['radar_lon'], data['radar_lat'], data['radi_h'][ie,ir], data['azim'][ia])

#        np.save('radi_h.npy', data['radi_h'])
#        np.save('lon.npy', data['lon'])
#        np.save('lat.npy', data['lat'])

    else:
        data['radi_h'] = np.load(radi_h)
        data['lon'] = np.load(lon)
        data['lat'] = np.load(lat)


    for ia in range(data['na']):
        data['hgt'][:,:,ia] = data['symhgt']

    return True


def  radar_topography( data , topo )    :
#Interpolate topography data to each ppi scan. 
   if not( 'radi_h' in data )   :
      #Data is not georeferenced... georeference data first.
      radar_georeference( data ) 

   topo_r = topo['range'][:,0]
   topo_a = topo['azimuth'][0,:]

   ne=data['ne']
   nr=data['nr']
   na=data['na']

   #Define a "interpolator"
   interpolator = interp2d(topo_a,topo_r,topo['mean'], kind='linear')

   data['topo'] = np.zeros((ne,nr,na), dtype='f4')

   #Loop over the vertical levels.
   for ie in range( 0 , ne )   :

       ppi_r = data['radi_h'][ie,:]
       ppi_a = data['azim']
        
       data['topo'][ie,:,:] = interpolator( ppi_a , ppi_r )

   return True


