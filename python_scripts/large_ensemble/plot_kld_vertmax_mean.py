# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016



@author:
"""
import sys
sys.path.append('../../common_python/common_functions/')
sys.path.append('../../common_python/common_modules/')

import numpy as np
import datetime as dt
import ctl_reader as ctlr
import os

import common_plot_functions as cpf
import common_mask_functions as cmf

basedir='/home/ra001011/a03471/data/output_data/'

expname = '/LE_D1_1km_5min/'

filetypes=['guesgp']   #analgp , analgz , guesgp , guesgz

max_dbz_levels=[20,40,50]   #Maximum reflectivity

plot_variables=['w','tk','qv','u','v','dbz']

lat_radar=34.823
lon_radar=135.523
radar_range=50.0e3   #Radar range in meters (to define the radar mask)

min_kld_value = 0.001 #Values below this wont be shown in the plots.

#Define initial and end times using datetime module.
itime = dt.datetime(2013,7,13,5, 5,0)  #Initial time.
etime = dt.datetime(2013,7,13,5,55,0)  #End time.

#Define the delta.
delta=dt.timedelta(seconds=300)

#Compute the total number of times
ntimes=int( 1 + np.around((etime-itime).seconds / delta.seconds) ) #Total number of times.

#=========================================================
#  LOOP OVER FILE TYPES
#=========================================================

for my_file_type in filetypes     :

  ctl_file = basedir + expname + '/ctl/' + my_file_type + '.ctl'

  plotbasedir=basedir + expname + '/plots/' + my_file_type + '/'

  #=========================================================
  #  READ CTL FILE
  #=========================================================

  ctl_dict = ctlr.read_ctl( ctl_file )

  nx=ctl_dict['nx']
  ny=ctl_dict['nx']
  nlev=ctl_dict['nz']
  nt=int(1)             #Force the number of times to be one.
  ctl_dict['nt']=int(1) #Force the number of times to be one.

  undef=np.float32( ctl_dict['undef'] )

  #=========================================================
  #  READ LAT LON AND GENERATE 2D RADAR MASK.
  #=========================================================

  latlon_file = basedir + expname + '/latlon/latlon.grd'

  tmp=ctlr.read_data_records(latlon_file,ctl=ctl_dict,records=np.array([0,1]))
  lat=tmp[:,:,1]
  lon=tmp[:,:,0]

  #Exclude areas outside the radar domain.
  radar_mask = cmf.distance_range_mask( lon_radar , lat_radar , radar_range , lon , lat )


  #=========================================================
  #  DEFINE REGIONS
  #=========================================================

  lati=np.array([34.75,34.6])
  late=np.array([35.25,34.9])
  loni=np.array([135.5,135.4])
  lone=np.array([136.25,135.7])

  zi=np.array([0,0,0])
  ze=np.array([nlev,nlev,nlev])

  reg_name='REG_1','REG_2','TOTAL'

  kld=dict()

  ctime=itime 

  #Add the global domain as a region.
  lati=np.append(lati,lat[0,0])
  late=np.append(late,lat[nx-1,ny-1])
  loni=np.append(loni,lon[0,0])
  lone=np.append(lone,lon[nx-1,ny-1])


  #Compute xi,xe,yi,ye corresponding to each region.
  xi , yi = cmf.lat_lon_to_i_j(lon,lat,loni,lati)
  xe , ye = cmf.lat_lon_to_i_j(lon,lat,lone,late)

  nregs=int( xi.shape[0] )

  #=========================================================
  #  DEFINE VARIABLES
  #=========================================================

  kld=dict()

  ens_mean=dict()

  kld_regional_mean=dict()
  kld_regional_max=dict()
  kld_regional_min=dict()

  kld_time_mean=dict()
  kld_time_std=dict()

  time_mean_max_dbz=np.zeros([nx,ny])

  my_kld=np.zeros([nx,ny,nlev])

  my_kld_sq=np.zeros([nx,ny,nlev])


  for var in plot_variables        :

    if var in ctl_dict['var_list']  :
 
       kld_regional_mean[var]=np.zeros([ntimes,nregs])
       kld_regional_max[var] =np.zeros([ntimes,nregs])
       kld_regional_min[var] =np.zeros([ntimes,nregs])

       kld_time_mean[var]=np.zeros([nx,ny])
       kld_time_std[var]=np.zeros([nx,ny])

  #=========================================================
  #  START TIME LOOP
  #=========================================================

  it=0
  while ( ctime <= etime ):

   print( ctime )

   print ( 'Reading the kld ')

   #=========================================================
   #  READ THE DATA
   #=========================================================

   my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ my_file_type + '/' + '/kldistance.grd'

   kld=ctlr.read_data_grads(my_file,ctl_dict,masked=False)

   mean_kld = dict()

   for my_var in kld  :
       kld[my_var][ kld[my_var] == undef ] = np.nan  #Filter undef bimodaility values.

       #Use the radar mask to define a masked array.
       tmp_mean_kld = np.squeeze( np.nanmean(kld[my_var],2) )

       tmp_mask = np.logical_not( radar_mask ) 
       mean_kld[my_var] = np.ma.masked_array( tmp_mean_kld  , mask = tmp_mask )

   #Read the ensemble mean to get the information from the storm location.
   my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ my_file_type + '/' + '/moment0001.grd'

   ens_mean=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
   
   #Compute max_dbz (we will use this to identify areas associated with clouds and convection)
   max_dbz = np.squeeze( np.nanmax(ens_mean['dbz'],2) )
   
   #=======================================================================================
   #Plot kld
   #=======================================================================================

   # for key in kld :

   for var in plot_variables        :

      my_kld = mean_kld[var]

      #========================================================================================
      #Generate regional averages of the moment.
      #========================================================================================

      kld_regional_mean[var][it,:],kld_regional_max[var][it,:],kld_regional_min[var][it,:]=cmf.get_regional_average_grid(my_kld.data,xi,xe,yi,ye,0,0,undef)

      #=========================================================================================
      #Acumulate the moment statistics in time. 
      #=========================================================================================

      kld_time_mean[var]=kld_time_mean[var] +  my_kld[:,:]                     #Accumulate the mean.
      kld_time_std[var]=kld_time_std[var]   + np.power(  my_kld[:,:]  , 2 )    #Accumulate the standard deviation.

   time_mean_max_dbz=time_mean_max_dbz + max_dbz  #To accumulate integrated liquid.

   #=========================================================================================
   #Advance time
   #=========================================================================================

   ctime = ctime + delta
 
   it = it + 1

  print ( "Finish time loop" )

  ntimes=it


  #=========================================================================================
  #Generate some time independent plots.
  #=========================================================================================


  #=========================================================================================
  #Plot the mean kld and its standard deviation.
  #=========================================================================================

  #for key in kld_time_mean  :
  for var in plot_variables        :

   if var in ctl_dict['var_list']  :

        
      my_kld_mean=kld_time_mean[var] / ntimes

      time_mean_max_dbz = time_mean_max_dbz / ntimes

      print('plotting the kld mean for variable ' + var)
      print('Min=' ,np.min(my_kld_mean),' Max=',np.max(my_kld_mean) )

      #Plot time mean of the moments.

      my_kld_sprd=np.sqrt( kld_time_std[var] / ntimes  - np.power( kld_time_mean[var], 2 ) / ntimes )

      cpf.set_default()  #Restore defaults
      my_map=cpf.cmap_discretize('Blues',10)
      cpf.figconf['figpath']=plotbasedir + '/time_independent_plots/'
      cpf.figconf['figsize']=(12,10)
      cpf.figconf['titlefontsize']=20
      cpf.figconf['labelfontsize']=20
      cpf.figconf['pcolor']=True
      cpf.figconf['shadedmin']=0.0
      if var == 'w'     :
         cpf.figconf['shadedmax']=0.1
      else              :
         cpf.figconf['shadedmax']=0.1
      cpf.figconf['shadedcolormap']='Blues'
      cpf.figconf['colorbar']=True
      cpf.figconf['colorbarfontsize']=15
      cpf.figconf['axessize']=[0.1,0.1,0.8,0.8]
      cpf.figconf['axesrange']=[134.97,136.09,34.36,35.30]
      cpf.figconf['contour']=True
      cpf.figconf['contourlevels']=[0.1e-6,0.2e-6,0.3e-6]
      cpf.figconf['contourcolormap']='Greys'
      cpf.figconf['gridline']=True
      cpf.figconf['xtick']=[134.5,135,135.5,136,136.5,137]
      cpf.figconf['ytick']=[34,34.5,35,35.5]
      cpf.figconf['ylabel']='Lat'
      cpf.figconf['xlabel']='Lon'

      cpf.figconf['title']='Mean KLD for variable ' + var 
      cpf.figconf['figname']='/Figure_KLD_time_mean_sprd_' + var 

      cpf.figconf['close'] = False
      cpf.plot_x_y_cartopy( lon , lat , my_kld_mean[:,:] , my_kld_sprd[:,:] , my_kld[:,:] ) 
      cpf.figconf['close']=True
      cpf.figconf['show']=True
      cpf.figconf['colorbar']=False
      cpf.figconf['shadedcolormap']='Greys'
      cpf.figconf['shadedmin']=0.0
      cpf.figconf['shadedmax']=1.0
      color_mask = 0.4*np.ones( np.shape( radar_mask ) )
      color_mask[ radar_mask ] = np.nan
      cpf.plot_x_y_cartopy( lon , lat , color_mask , my_kld_sprd[:,:] , color_mask )


  #=========================================================================================
  #Plot time evolution of kld over different regions.
  #=========================================================================================

  #for key in kld_time_mean  :
  for var in plot_variables        :

   if var in ctl_dict['var_list']  :


      #Plot time series of moments for different regions.

      for ireg in range(0,nregs) :

         iregstr="%04d" % ( ireg + 1 )
         cpf.set_default()  #Restore defaults
         cpf.figconf['figpath']=plotbasedir + '/time_independent_plots/'
         cpf.figconf['figname']='Figure_kld_reg_' + iregstr + '_var_' + var 
         cpf.figconf['figsize']=(12,10)
         cpf.figconf['title']='Mean KLD for var ' + var + ' for region ' + reg_name[ireg] 
         cpf.figconf['titlefontsize']=20
         cpf.figconf['labelfontsize']=20
         cpf.figconf['axessize']=[0.1,0.1,0.8,0.8]
         cpf.figconf['gridline']=True

         cpf.figconf['linestyle']=['-']
         cpf.figconf['linemarker']=['o']
         cpf.figconf['linecolor']=['b']
         cpf.figconf['ylabel']=['KLD']
         cpf.figconf['xlabel']=['Time (seconds)']
         #cpf.figconf['show']=True

         time = np.arange( 0 , np.shape(kld_regional_mean[var])[0] ) * delta.total_seconds()

         cpf.plot_lines( time , kld_regional_mean[var][:,ireg] )



