# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016

#Compute KLD profile over rainny and non-rainny grid points and store the result in a file.

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
from common_functions import common_functions as comm 

lat_radar=34.823
lon_radar=135.523
radar_range=50.0e3   #Radar range in meters (to define the radar mask)

file_name='moment0001'

basedir='/home/jruiz/share/LARGE_ENSEMBLE/output_data/'

expnames=['LE_D1_1km_5min_OFP_V2','LE_D1_1km_2min_OFP_V2','LE_D1_1km_1min','LE_D1_1km_30sec_OFP_V2','LE_D1_1km_5min_4D_OFP_V2','LE_D1_1km_1min_4D']
#basedir='/home/jruiz/share/LARGE_ENSEMBLE/output_data/home/ra001011/a03471/data/output_data/'
#expnames=['LE_D1_1km_5min','LE_D1_1km_2min','LE_D1_1km_1min','LE_D1_1km_30sec','LE_D1_1km_30sec_nospinup','LE_D1_1km_1min_4D']

cvars=['dbz','tk','qv','u','v','w']

deltat=[300,240,300,300,300,300]

init_date = ['20130713050500','20130713050400','20130713050500','20130713050500','20130713050500','20130713050500']

filetypes=['guesgp','analgp']   #analgp , analgz , guesgp , guesgz

rain_threshold = 30.0     #Maximum reflectivities above this threshold will be considered rainny points.   
norain_threshold = 0.0    #Maximum reflectivities below this threshold will be considered no rain points.

#Define initial and end times using datetime module.
etime = dt.datetime(2013,7,13,6,0,30)  #End time.

#=========================================================
#  LOOP OVER FILE TYPES
#=========================================================

for iexp , my_exp in enumerate( expnames ) :

   delta = dt.timedelta(seconds=deltat[iexp])

   itime = dt.datetime.strptime( init_date[iexp] , '%Y%m%d%H%M%S' )

   #Compute the total number of times
   ntimes=int( 1 + np.around((etime-itime).seconds / delta.seconds) ) #Total number of times.

   for my_file_type in filetypes     :

      ctl_file = basedir + '/' + my_exp + '/ctl/' + my_file_type + '.ctl'

      outputdir=basedir + '/' + my_exp + '/time_mean/' + my_file_type + '/'

      if not os.path.exists( outputdir)  :

         os.makedirs( outputdir )

      #=========================================================
      #  READ CTL FILE
      #=========================================================

      ctl_dict = ctlr.read_ctl( ctl_file )

      nx=ctl_dict['nx']
      ny=ctl_dict['nx']
      nlev=len( ctl_dict['full_lev_list'] )
      nt=int(1)             #Force the number of times to be one.
      ctl_dict['nt']=int(1) #Force the number of times to be one.
      nz=np.size( ctl_dict['vlevels'] )

      undef=np.float32( ctl_dict['undef'] )

      if  ctl_dict['big_endian']   :
         dtypein = '>f4'
         endian='big_endian'
      else                         :
         dtypein = 'f4'
         endian='little_endian'

      if  ctl_dict['sequential']   :
         access='sequential'
      else                         :
         access='direct'

      sequential=ctl_dict['sequential']

      #=========================================================
      #  READ LAT LON AND GENERATE 2D RADAR MASK.
      #=========================================================

      latlon_file = basedir + my_exp + '/latlon/latlon.grd'

      tmp=ctlr.read_data_records(latlon_file,ctl=ctl_dict,records=np.array([0,1]))
      lat=tmp[:,:,1]
      lon=tmp[:,:,0]

      #Exclude areas outside the radar domain.
      tmp_radar_mask = cmf.distance_range_mask( lon_radar , lat_radar , radar_range , lon , lat )

      #Transform 2D radar mask to 3D.
      radar_mask = np.zeros( (np.shape(lat)[0], np.shape(lat)[1],nz) ).astype(bool)

      for kk in range(0,nz)   :
         radar_mask[:,:,kk] = tmp_radar_mask 

      #=========================================================
      #  START TIME LOOP
      #=========================================================

      it=0
 
      ctime = itime 

      max_dbz = np.zeros( (np.shape(lat)[0], np.shape(lat)[1],nz) )

      parameter_rain = dict()
      parameter_norain = dict()

      while ( ctime <= etime )  :

         print( ctime )

         print ( 'Reading data')

         #=========================================================
         #  READ THE DATA
         #=========================================================

         my_file=basedir + my_exp + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/'+ my_file_type + '/' + '/' + file_name +'.grd'

         #Read all the variables and levels at once
         tmp_parameter = ctlr.read_data_grads(  my_file , ctl_dict ) 

         my_file=basedir + my_exp + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/'+ my_file_type + '/moment0001.grd'       
 
         ensemble_mean = ctlr.read_data_grads(  my_file , ctl_dict )
 
         ensemble_mean[ ensemble_mean == undef ] = np.nan

         tmp_max_dbz = np.nanmax( ensemble_mean['dbz'] , 2 )

         for kk in range(0,nz)   :
            max_dbz[:,:,kk] = np.squeeze( tmp_max_dbz[:] )

         #Several mask can be defined here.

         rain_mask = np.logical_and( max_dbz > rain_threshold , radar_mask )

         norain_mask = np.logical_and( max_dbz < norain_threshold , radar_mask )


         for my_var in cvars   :

            if it == 0 :
               parameter_rain[my_var] = np.zeros( (nz,ntimes) )
               parameter_norain[my_var] = np.zeros( (nz,ntimes) )
  
            tmp_rain = np.copy( np.squeeze( tmp_parameter[my_var] ) )
            tmp_rain[ tmp_rain == undef ] = np.nan
            tmp_rain[ np.logical_not( rain_mask ) ] = np.nan

            tmp_norain = np.copy( np.squeeze( tmp_parameter[my_var] ) )
            tmp_norain[ tmp_norain == undef ] = np.nan 
            tmp_norain[ np.logical_not( norain_mask ) ] = np.nan

            #import matplotlib.pyplot as plt 
 
            #plt.pcolor( rain_mask[:,:,5] )
            #plt.colorbar()

            #plt.show()

            #plt.pcolor( norain_mask[:,:,5] )
            #plt.colorbar()

            #plt.show()
            #plt.pcolor( tmp_rain[:,:,5] )
            #plt.colorbar()
            #plt.show()

            #plt.pcolor( tmp_norain[:,:,5] )
            #plt.colorbar()
            #plt.show()


            #Compute the vertical profile
            for kk in range(0,nz)   :

               if my_var == 'w' :
                  tmp_rain[:,:,kk][ tmp_rain[:,:,kk] < 0 ]=np.nan
                  tmp_norain[:,:,kk][ tmp_norain[:,:,kk] < 0 ]=np.nan
                  parameter_rain[my_var][kk,it] = np.nanmean( np.nanmean( tmp_rain[:,:,kk] , 1 ) , 0 )
                  parameter_norain[my_var][kk,it] = np.nanmean( np.nanmean( tmp_norain[:,:,kk] , 1 ) , 0 )
               else             :
                  parameter_rain[my_var][kk,it] = np.nanmean( np.nanmean( tmp_rain[:,:,kk] , 1 ) , 0 )
                  parameter_norain[my_var][kk,it] = np.nanmean( np.nanmean( tmp_norain[:,:,kk] , 1 ) , 0 ) 

         ctime = ctime + delta

         it = it + 1

      print ( "Finish time loop" )

      for my_var in cvars  :

         my_file=outputdir + '/' + file_name + '_' + my_var + '_rain_profile.npz'
         np.savez(my_file, parameter = parameter_rain[my_var] )

         my_file=outputdir + '/' + file_name +  '_' + my_var + '_norain_profile.npz'
         np.savez(my_file, parameter = parameter_norain[my_var] )

         #import matplotlib.pyplot as plt
         #plt.contourf( np.transpose( parameter_rain[my_var] ) )
         #plt.colorbar()
         #plt.show()
     

     




