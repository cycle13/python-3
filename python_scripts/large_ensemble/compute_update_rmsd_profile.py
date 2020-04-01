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

basedir='/home/jruiz/share/LARGE_ENSEMBLE/output_data/home/ra001011/a03471/data/output_data/'

expnames=['LE_D1_1km_5min','LE_D1_1km_30sec']

deltat=[300,300]

init_date = ['20130713050500','20130713050500']

filetypes=['guesgp']   #analgp , analgz , guesgp , guesgz

variable_combination=[['dbz','dbz','dbz','dbz','dbz','v' ,'v' ,'v','v' ],
                      ['tk' ,'v'  ,'qv' ,'w'  ,'dbz','tk','qv','w','v' ]]

obs_increment = [5.0 , 5.0 , 5.0 , 5.0 , 5.0 , 2.0 , 2.0 , 2.0 , 2.0 ]
obs_error     = [5.0 , 5.0 , 5.0 , 5.0 , 5.0 , 2.0 , 2.0 , 2.0 , 2.0 ]

nbv=1000

rain_threshold = 30.0     #Maximum reflectivities above this threshold will be considered rainny points.   
norain_threshold = 0.0    #Maximum reflectivities below this threshold will be considered no rain points.

#Define initial and end times using datetime module.
etime = dt.datetime(2013,7,13,5,59,30)  #End time.

#=========================================================
#  LOOP OVER FILE TYPES
#=========================================================

for iexp , my_exp in enumerate( expnames ) :

   delta = dt.timedelta(seconds=deltat[iexp])

   itime = dt.datetime.strptime( init_date[iexp] , '%Y%m%d%H%M%S' )

   #Compute the total number of times
   ntimes=int( np.around((etime-itime).seconds / delta.seconds) ) #Total number of times.


   ctl_file = basedir + '/' + my_exp + '/ctl/update_mean_diff.ctl'

   ctl_file_2 = basedir + '/' + my_exp + '/ctl/moment0001_for.ctl'

   outputdir=basedir + '/' + my_exp + '/time_mean/guesgp/'

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

   #=========================================================
   #  START TIME LOOP
   #=========================================================

   for iv,my_increment in enumerate( obs_increment )   :

      var_obs = variable_combination[0][iv]
      var_upd = variable_combination[1][iv]

      profile_update_mean_rmsd_rain = np.zeros( (nz-1,ntimes) )
      profile_update_mean_kf_rain = np.zeros( (nz-1,ntimes) )
      profile_update_kld_rain = np.zeros( (nz-1,ntimes) )

      profile_update_mean_rmsd_norain = np.zeros( (nz-1,ntimes) )
      profile_update_mean_kf_norain = np.zeros( (nz-1,ntimes) )
      profile_update_kld_norain = np.zeros( (nz-1,ntimes) )
      
      profile_update_mean_bias_rain = np.zeros( (nz-1,ntimes) )
      profile_update_mean_bias_norain = np.zeros( (nz-1,ntimes) )

      it=0
      ctime = itime
      while ( ctime <= etime )  :

         print( ctime )

         print ( 'Reading data')

         #=========================================================
         #  OBTAIN THE RAIN MASK
         #=========================================================

         my_file=basedir + my_exp + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/moment0001.grd'

         #Read all the variables and levels at once
         ensemble_mean = ctlr.read_data_grads(  my_file , ctl_dict_2 , undef2nan = True ) 

         max_dbz = np.squeeze( np.nanmax( np.delete( ensemble_mean['dbz'] , 4 , 2) , 2 ) )

         rain_mask = np.logical_and( max_dbz > rain_threshold , radar_mask )

         norain_mask = np.logical_and( max_dbz < norain_threshold , radar_mask )

         #=========================================================
         #  READ THE DATA
         #=========================================================

         my_file=basedir + '/' + my_exp + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/update_comp_meandiff_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '.grd'
         update_mean_diff = np.squeeze( ctlr.read_data(  my_file , ctl_dict , undef2nan = True ) )
         my_file=basedir + '/' + my_exp + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/update_comp_updated_mean_kf_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '.grd'
         update_mean_kf = np.squeeze( ctlr.read_data(  my_file , ctl_dict , undef2nan = True ) )
         my_file=basedir + '/' + my_exp + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/update_comp_kld_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '.grd'
         update_kld = np.squeeze( ctlr.read_data(  my_file , ctl_dict , undef2nan = True ) )
         
         update_mean_bias = np.copy(update_mean_diff)
         
         update_mean_diff = np.power( np.delete( update_mean_diff , 4 ,2 )  , 2 )

         update_mean_kf   = np.power( np.delete( update_mean_kf - np.squeeze(ensemble_mean[ var_upd ]) , 4 ,2 ) , 2 )

         update_kld = np.delete( update_kld , 4 , 2 )

         update_mean_diff[ update_mean_diff == 0.0 ] = np.nan
         update_mean_kf[ update_mean_kf == 0.0 ] = np.nan
         update_kld[ update_kld == 0.0 ]=np.nan 

         #Compute the vertical profile
         for kk in range( update_mean_diff.shape[2])   :
             
            profile_update_mean_bias_rain[kk,it] = np.nanmean( update_mean_diff[:,:,kk][ rain_mask ]  ) 
            profile_update_mean_bias_norain[kk,it] = np.nanmean( update_mean_diff[:,:,kk][ norain_mask ]  ) 


            profile_update_mean_rmsd_rain[kk,it] = np.sqrt( np.nanmean( update_mean_diff[:,:,kk][ rain_mask ]  )  )
            profile_update_mean_rmsd_norain[kk,it] = np.sqrt( np.nanmean( update_mean_diff[:,:,kk][ norain_mask ]  )  )

            profile_update_mean_kf_rain[kk,it] = np.sqrt( np.nanmean( update_mean_kf[:,:,kk][ rain_mask ]  )  )
            profile_update_mean_kf_norain[kk,it] = np.sqrt( np.nanmean( update_mean_kf[:,:,kk][ norain_mask ]  )  )

            profile_update_kld_rain[kk,it] = np.nanmean( update_kld[:,:,kk][ rain_mask ]  ) 
            profile_update_kld_norain[kk,it] = np.nanmean( update_kld[:,:,kk][ norain_mask ]  ) 

         ctime = ctime + delta

         it = it + 1

      print ( "Finish time loop" )

      my_file=outputdir + '/update_mean_rmsd_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '_rain_profile.npz'
      np.savez(my_file, parameter = profile_update_mean_rmsd_rain )
      my_file=outputdir + '/update_mean_rmsd_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '_norain_profile.npz'
      np.savez(my_file, parameter = profile_update_mean_rmsd_norain )

      my_file=outputdir + '/update_mean_bias_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '_rain_profile.npz'
      np.savez(my_file, parameter = profile_update_mean_bias_rain )
      my_file=outputdir + '/update_mean_bias_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '_norain_profile.npz'
      np.savez(my_file, parameter = profile_update_mean_bias_norain )

      my_file=outputdir + '/update_mean_kf_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '_rain_profile.npz'
      np.savez(my_file, parameter = profile_update_mean_kf_rain )
      my_file=outputdir + '/update_mean_kf_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '_norain_profile.npz'
      np.savez(my_file, parameter = profile_update_mean_kf_norain )

      my_file=outputdir + '/update_kld_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '_rain_profile.npz'
      np.savez(my_file, parameter = profile_update_kld_rain )
      my_file=outputdir + '/update_kld_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '_norain_profile.npz'
      np.savez(my_file, parameter = profile_update_kld_norain )
      
     
      print( np.nanmin( profile_update_mean_rmsd_rain ) , np.nanmax( profile_update_mean_rmsd_rain ) ) 



