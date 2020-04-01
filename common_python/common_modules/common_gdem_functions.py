import gdal
import os.path
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('../../common_python/common_functions/')

#The data that is converted by these functions comes from:
#https://search.earthdata.nasa.gov/search/granules?p=C197265171-LPDAAC_ECS&m=-0.140625!0!2!1!0!0%2C2&q=aster&ok=aster
#Aster Global Digital Elevation Model V0002


from common_functions import common_functions as cf


undef = -999 #Is there an undef value in raster gdem data?

def raster_gdem_v2_to_radar( cen_lon , cen_lat , rmin , dr , nr , amin , da , na , datapath , prefix='ASTGTM2_' , sufix='_dem.tif') :

    #Given a certain region this function will read the required tiles and merge them
    #into a single numpy array. Lat and lon will also be computed.

    #Inputs:
    #cen_lon : the longitud of the 0 range.
    #cen_lat : the latitude of the 0 range.
    #rmin: minimum range of the grid.
    #dr: radial resolution of the grid.
    #nr: number of grid points in the radial direction
    #amin: minimum azimuth
    #da: azimuth resolution
    #na: number of azimuth
    #datapath: is the path containing raster files.
    #prefix and sufix are part of the file naming convention. Default values are according to the
    #names originally given to the files. 
  
    #Define the polar coordinte grid

    output_r=np.arange(rmin,rmin+nr*dr,dr)
    output_a=np.arange(amin,amin+na*da,da)

    [ output_a ,  output_r ] = np.meshgrid( output_a , output_r )

    #print( output_r.shape  )
    #print( output_a.shape  )
    #print( output_r[:,0] )
    #print( output_a[0,:] )

    #Get lat and lon corresponding to the grid points in the polar coordinate.

    [output_lon,output_lat]=cf.com_ra_to_ll(cen_lon=cen_lon,cen_lat=cen_lat
                                          ,r=output_r,a=output_a,nr=nr,na=na)
 
    output_nx = output_lon.shape[0]
    output_ny = output_lon.shape[1]

    #print( output_r[:,0] )
    #print( output_a[0,:] )

    #print( output_lon[:,0] )
    #print( output_lat[0,:] )

    #DEBUG DEBUG Test if the ra-to-ll and ll-to-ra transformations are consistent.
    #[output_r_inv , output_a_inv ] = cf.com_ll_to_ra(cen_lon=cen_lon,cen_lat=cen_lat
    #                                                 ,lon=output_lon,lat=output_lat
    #                                                 ,nx=output_nx,ny=output_ny) 
    #print( np.max( np.abs( output_r_inv - output_r ) ) )
    #print( np.max( np.abs( output_a_inv - output_a ) ) )
    #
    #
    #print( np.mean( np.abs( output_r_inv - output_r ) ) )
    #print( np.mean( np.abs( output_a_inv - output_a ) ) )
    # According to these test mean absolute error for range is 0.2 meters
    #                         mean absolute error for azimuth is 0.01 degrees
    #                         larger errors are expected closer to the radar.


    #Get the intiger limits (for tile searching)

    maxlon_i = int( np.ceil( np.max(output_lon) )  )
    minlon_i = int( np.floor( np.min(output_lon) ) )
    maxlat_i = int( np.ceil( np.max(output_lat) )  )
    minlat_i = int( np.floor( np.min(output_lat) ) )

    #Allocate the output arrays.
    #Order flag is required to produce fortran contiguos arrays.
    output_mean_topo = np.zeros( ( nr , na ) , order='F' , dtype=np.float32 )
    output_max_topo  = np.zeros( ( nr , na ) , order='F' , dtype=np.float32 )
    output_min_topo  = np.zeros( ( nr , na ) , order='F' , dtype=np.float32 )
    output_n_topo    = np.zeros( ( nr , na ) , order='F' , dtype=np.int32 )

    #Double loops over raster files.

    for ilon in range( minlon_i , maxlon_i + 1 )    :
        for ilat in range( minlat_i , maxlat_i + 1 )   :

             #Form the name of the current tile.
             latstr='%0*d' % (2,ilat)
             lonstr='%0*d' % (3,ilon)
             sn='N'
             if ilat < 0 :
               sn='S'
             we='W'
             if ilon > 0 :
               we='E'

             my_tile= sn + latstr + we + lonstr

             #Read the data
             raster_file=datapath + '/' + prefix + my_tile + sufix
             print('Reading and interpolating data from: ' + raster_file )

             if not os.path.isfile(raster_file)            :
                print('[Warning]:Path ' + raster_file + ' does not exist: SKIP IT')
             else                                          :

                [raster_lon,raster_lat,raster_data]=read_raster_data(raster_file)
                raster_nx=raster_data.shape[0]
                raster_ny=raster_data.shape[1]

                if (np.min(raster_lon) > np.max(output_lon) ) |   \
                   (np.min(raster_lat) > np.max(output_lat) ) |   \
                   (np.max(raster_lat) < np.min(output_lat) ) |   \
                   (np.max(raster_lon) < np.min(output_lon) )     :
             
                   print('File ' + raster_file + ' does not have useful data: SKIP IT') 
             
                else                                             :

                   #plt.pcolor(raster_data[0:1000,0:1000,0])
                   #plt.colorbar()
                   #plt.show()

                   #Convert raster lat lon to range and azimuth.
                   [raster_r,raster_a]=cf.com_ll_to_ra(cen_lon=cen_lon,cen_lat=cen_lat
                                                         ,lon=raster_lon,lat=raster_lat
                                                         ,nx=raster_nx,ny=raster_ny)

                   #Interpolate using box averaging.
                   cf.com_interp_boxavereg(xini=rmin,dx=dr,nx=nr
                                    ,yini=amin,dy=da,ny=na
                                    ,xin=raster_r.reshape( raster_nx * raster_ny )
                                    ,yin=raster_a.reshape( raster_nx * raster_ny )
                                    ,datain=raster_data.reshape( raster_nx * raster_ny ) 
                                    ,nin=raster_nx*raster_ny 
                                    ,data_sum=output_mean_topo
                                    ,data_max=output_max_topo,data_min=output_min_topo
                                    ,data_n=output_n_topo,undef=undef)

    #Compute the mean topography.
    mask = output_n_topo > 0 
    output_mean_topo[ mask ] = output_mean_topo[ mask ] / output_n_topo[ mask ]

    #Complete missing values using neighborhood values.
    mask = output_n_topo == 0  #Identify missing data points
    output_mean_topo=cf.com_complete_missing_2d(field=output_mean_topo,missing_mask=mask,
                                             nx=output_nx,ny=output_ny,npass=4)

    output_max_topo=cf.com_complete_missing_2d(field=output_max_topo,missing_mask=mask,
                                             nx=output_nx,ny=output_ny,npass=4)

    output_min_topo=cf.com_complete_missing_2d(field=output_min_topo,missing_mask=mask,
                                             nx=output_nx,ny=output_ny,npass=4)



    #Put the output into a dictionary.
    output=dict()
    output['longitude']=output_lon
    output['latitude']=output_lat
    output['range']=output_r
    output['azimuth']=output_a
    output['mean']=output_mean_topo
    output['max']=output_max_topo
    output['min']=output_min_topo
    output['number']=output_n_topo

    #TODO : output should be a masked array

    return output

def read_raster_data(inputfile)  :

    my_raster = gdal.Open(inputfile)

    nx=my_raster.RasterXSize
    ny=my_raster.RasterYSize

    nb=my_raster.RasterCount

    my_raster_data=np.zeros((nx,ny,nb))

    #Read the data and store it into a numpy array

    for ib in range( 0 , nb )  :

        my_raster_data[:,:,ib] = my_raster.ReadAsArray(ib)


    #Get the lat and lons.
    [lon,lat]=get_lat_lon(inputfile)

    #print(lon[0,:])
    #print(lat[:,0])

    return lon , lat , my_raster_data



def get_lat_lon(inputfile)  :

   my_raster = gdal.Open(inputfile)

   gt = my_raster.GetGeoTransform()
   proj = my_raster.GetProjection()

   xres = gt[1]
   yres = gt[5]

   # get the edge coordinates and add half the resolution 
   # to go to center coordinates
   xmin = gt[0] + xres * 0.5
   xmax = gt[0] + (xres * my_raster.RasterXSize) - xres * 0.5
   ymin = gt[3] + (yres * my_raster.RasterYSize) + yres * 0.5
   ymax = gt[3] - yres * 0.5

   lon=np.zeros( my_raster.RasterXSize )
   lat=np.zeros( my_raster.RasterYSize )
   for ii in range( 0 , my_raster.RasterXSize )  :
     lon[ii] = xmin + ii * xres 
   for ii in range( 0 , my_raster.RasterYSize )  :
     lat[ii]= ymax + ii*yres 

   my_raster = None

   # create a grid of xy coordinates in the original projection
   [lon,lat] = np.meshgrid(lon,lat)

   return lon , lat


