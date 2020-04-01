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

basedir='/home/ra001011/a03471/data/output_data/'

expnames=['LE_D1_1km_5min','LE_D1_1km_2min','LE_D1_1km_1min','LE_D1_1km_30sec','LE_D1_1km_30sec_nospinup','LE_D1_1km_1min_4D']

cvars=['dbz','tk','qv','u','v','w']

deltat=[300,240,300,300,300,300]

init_date = ['20130713050500','20130713050400','20130713050500','20130713050500','20130713050500','20130713050500']

filetype='guesgp'   #analgp , analgz , guesgp , guesgz

rain_threshold = 30.0     #Maximum reflectivities above this threshold will be considered rainny points.   
norain_threshold = 0.0    #Maximum reflectivities below this threshold will be considered no rain points.

#Define initial and end times using datetime module.
etime = dt.datetime(2013,7,13,5,59,30)  #End time.

nbv=1000

variable_combination=[['dbz','dbz','dbz','dbz','dbz','v' ,'v' ,'v','v' ],
                      ['tk' ,'v'  ,'qv' ,'w'  ,'dbz','tk','qv','w','v' ]]

obs_increment = [5.0 , 5.0 , 5.0 , 5.0 , 5.0 , 2.0 , 2.0 , 2.0 , 2.0 ]
obs_error     = [5.0 , 5.0 , 5.0 , 5.0 , 5.0 , 2.0 , 2.0 , 2.0 , 2.0 ]



#=========================================================
#  LOOP OVER FILE TYPES
#=========================================================

for iexp , my_exp in enumerate( expnames ) :

   delta = dt.timedelta(seconds=deltat[iexp])

   itime = dt.datetime.strptime( init_date[iexp] , '%Y%m%d%H%M%S' )

   #Compute the total number of times
   ntimes=int( 1 + np.around((etime-itime).seconds / delta.seconds) ) #Total number of times.


   ctl_file_2 = basedir + '/' + my_exp + '/ctl/' + filetype + '.ctl'
   ctl_file   = basedir +  '/' + my_exp + '/ctl/update_mean_diff.ctl'

   outputdir=basedir + '/' + my_exp + '/time_mean/' + filetype + '/'

   if not os.path.exists( outputdir)  :

      os.makedirs( outputdir )

   #=========================================================
   #  READ CTL FILE
   #=========================================================

   ctl_dict = ctlr.read_ctl( ctl_file )
   ctl_dict_2 = ctlr.read_ctl( ctl_file_2 )

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

   if iexp == 0  :

      latlon_file = basedir + my_exp + '/latlon/latlon.grd'

      tmp=ctlr.read_data_records(latlon_file,ctl=ctl_dict,records=np.array([0,1]))
      lat=tmp[:,:,1]
      lon=tmp[:,:,0]

      #Exclude areas outside the radar domain.
      radar_mask = cmf.distance_range_mask( lon_radar , lat_radar , radar_range , lon , lat )

   #Loop over variables combination

   for iv in range( len( variable_combination[0] ) ) :   #Loop over selected variable combinations.

      var_obs = variable_combination[0][iv]
      var_upd = variable_combination[1][iv]


      #=========================================================
      #  START TIME LOOP
      #=========================================================

      it=0

      ctime = itime 

      rmse_rain = np.zeros((nz,ntimes))
      rmse_norain = np.zeros((nz,ntimes))
      bias_rain = np.zeros((nz,ntimes))
      bias_norain = np.zeros((nz,ntimes))
      kld_rain = np.zeros((nz,ntimes))
      kld_norain = np.zeros((nz,ntimes))

      while ( ctime <= etime )  :

         print( ctime )

         print ('Reading data')

         #=========================================================
         #  READ THE DATA
         #=========================================================

         my_file=basedir + '/' + my_exp + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + filetype + '/update_comp_meandiff_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '.grd'

         #Read all the variables and levels at once
         parameter = ctlr.read_data(  my_file , ctl_dict , undef2nan = True ) 

         my_file=basedir + '/' + my_exp + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + filetype + '/update_comp_kld_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '.grd'

         #Read all the variables and levels at once
         kld = ctlr.read_data(  my_file , ctl_dict , undef2nan = True )

         my_file=basedir + my_exp + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/moment0001.grd'       

         max_dbz =  np.squeeze( np.nanmax( ctlr.read_data_grads(  my_file , ctl_dict_2 , undef2nan = True )['dbz']  , 2 ) )

         #Several mask can be defined here.

         rain_mask = np.logical_and( max_dbz > rain_threshold , radar_mask )

         norain_mask = np.logical_and( max_dbz < norain_threshold , radar_mask )


         for kk in range( 0 , nz ) :

            nan_mask = np.logical_not( np.isnan( np.squeeze( parameter[:,:,kk] ) ) ) 
 
            tmp_rain_mask = np.logical_and( nan_mask , rain_mask )
            tmp_norain_mask = np.logical_and( nan_mask , norain_mask )

            tmp=np.squeeze( parameter[:,:,kk] )

            rmse_rain[kk,it] = np.sqrt( np.nanmean( np.power( tmp[tmp_rain_mask] , 2) ) )
            rmse_norain[kk,it] = np.sqrt( np.nanmean( np.power( tmp[ tmp_norain_mask ] , 2 ) ) )

            bias_rain[kk,it] = np.nanmean( tmp[ tmp_rain_mask ] ) 
            bias_norain[kk,it] = np.nanmean( tmp[ tmp_norain_mask ] ) 
 
            tmp=kld[:,:,kk]
 
            kld_rain[kk,it] = np.nanmean( tmp[tmp_rain_mask ] ) 
            kld_norain[kk,it] = np.nanmean( tmp[ tmp_norain_mask ] ) 


         ctime = ctime + delta

         it = it + 1

      print ( "Finish time loop" )


      my_file=outputdir+ '/update_comp_rmse_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '_rain_profile.npz'
      np.savez(my_file, parameter = rmse_rain )
      my_file=outputdir+ '/update_comp_rmse_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '_norain_profile.npz'
      np.savez(my_file, parameter = rmse_norain )

      my_file=outputdir+ '/update_comp_bias_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '_rain_profile.npz'
      np.savez(my_file, parameter = bias_rain )
      my_file=outputdir+ '/update_comp_bias_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '_norain_profile.npz'
      np.savez(my_file, parameter = bias_norain )
      my_file=outputdir+ '/update_comp_kld_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '_rain_profile.npz'
      np.savez(my_file, parameter = kld_rain )
      my_file=outputdir+ '/update_comp_kld_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '_norain_profile.npz'
      np.savez(my_file, parameter = kld_norain )

      print( my_exp )
      print('My var obs ',var_obs,' My vr upd ',var_upd)
      print('Print RMSE rain ',np.nanmax( rmse_rain ) , np.nanmin( rmse_rain ) )
      print('Print RMSE norain ',np.nanmax( rmse_norain ) , np.nanmin( rmse_norain ) )
      print('Print BIAS rain ',np.nanmax( bias_rain ) , np.nanmin( bias_rain ) )
      print('Print BIAS norain ',np.nanmax( bias_norain ) , np.nanmin( bias_norain ) )
      print('Print KLD rain ',np.nanmax( kld_rain ) , np.nanmin( kld_rain ) )
      print('Print KLD norain ',np.nanmax( kld_norain ) , np.nanmin( kld_norain ) )


