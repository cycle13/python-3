# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016

Porpouse compare linear and non-linear updates.
One variable will be assumed as the observed variable and a second variable (may be the same)
Will be the updated variable. The properties of the posterior distribution obtained using a PF and an ENKF 
will be compared.
The variables are updated at the same grid point where the observation is assumed. 

@author:
"""
import sys
sys.path.append('../../common_python/common_functions/')
sys.path.append('../../common_python/common_modules/')

import numpy as np
#import matplotlib.pyplot as plt
import datetime as dt
import ctl_reader as ctlr
import binary_io  as bio
import os

from common_functions import common_functions as comm


basedir='/home/ra001011/a03471/data/output_data/'

#expnames  = ['LE_D1_1km_30sec','LE_D1_1km_5min','LE_D1_1km_1min','LE_D1_1km_1min_4D','LE_D1_1km_30sec_nospinup']
#expnames  = ['LE_D1_1km_30sec','LE_D1_1km_5min','LE_D1_1km_30sec_nospinup','LE_D1_1km_1min','LE_D1_1km_1min_4D']
expnames  = ['LE_D1_1km_2min']

variable_combination=[['dbz','dbz','dbz','dbz','dbz','v' ,'v' ,'v','v' ],
                      ['tk' ,'v'  ,'qv' ,'w'  ,'dbz','tk','qv','w','v' ]]

#expdeltas = [300,300,300,300,300,300]
expdeltas = [120]

#delta=dt.timedelta(seconds=60)  #Original data is every 30 seconds

filetypes=['guesgp']   #Only guess will be processed.

smooth=False           #False- no smooth , True apply smooth
smooth_lambda=10       #Smooth length scale (in number of grid points)

nbv=1000              #Total number of ensemble members.

obs_increment = [5.0 , 5.0 , 5.0 , 5.0 , 5.0 , 2.0 , 2.0 , 2.0 , 2.0 ]
obs_error     = [5.0 , 5.0 , 5.0 , 5.0 , 5.0 , 2.0 , 2.0 , 2.0 , 2.0 ]

itime = dt.datetime(2013,7,13,5, 4,0)  #Initial time.
etime = dt.datetime(2013,7,13,5,58,0)  #End time.

#=========================================================
#  LOOP OVER FILE TYPES
#=========================================================

for iexp , my_exp_name in enumerate( expnames ) :

   delta=dt.timedelta( seconds = expdeltas[iexp] )

   for my_file_type in filetypes  :

      #=========================================================
      #  READ CTL FILE
      #=========================================================

      ctl_file = basedir + my_exp_name + '/ctl/' + my_file_type + '.ctl' 

      ctl_dict = ctlr.read_ctl( ctl_file )


      #Estimate number of bins in histogram. 
      nbins=int(np.round(np.sqrt(nbv))) #Number of bins to be used in histogram computation.


      nx=np.int(ctl_dict['nx'])
      ny=np.int(ctl_dict['ny'])
      nz=np.array(ctl_dict['end_record']).max() + 1 #Total number of records in binary file.
      nlev=ctl_dict['nz']                 #Number of vertical levels for 3D variables.

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

      #All the fields will be read so.
      n_selected_fields = nz
      selected_fields = np.arange(0,nz) + 1

      #=========================================================
      #  START COMPUTATION LOOP
      #=========================================================

      my_ensemble=np.zeros([nx,ny,nz,nbv]).astype('float32')

      ctime = itime

      while ( ctime <= etime ):

        print(ctime)

        print(' Reading the ensemble ')
 
        my_path=basedir + '/' + my_exp_name + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/' 

        if not os.path.exists(my_path)                       :
           print('[Warning]:Path ' + my_path + ' does not exist: SKIP IT')

        else                                                 :
  
          #Read the ensemble and the undef mask
          [my_ensemble , my_undefmask] = comm.read_ensemble(path=my_path,nx=nx,ny=ny,nbv=nbv
                                                        ,selected_fields=selected_fields
                                                        ,n_selected_fields=n_selected_fields
                                                        ,undef=undef,ie=endian,acc=access)

          for iv in range( len( variable_combination[0] ) ) :   #Loop over selected variable combinations.

             var_obs = variable_combination[0][iv] 
             var_upd = variable_combination[1][iv]
 
             print('Computing update comparisson for variables ',var_obs,' and ',var_upd) 

             #Find the indices corresponding to var_obs and var_upd in the ensemble array. 
             
             [var_list , lev_list , obs_record_list] =ctlr.record_subset( ctl_dict , [var_obs] , ctl_dict['vlevels'] )
             [var_list , lev_list , upd_record_list] =ctlr.record_subset( ctl_dict , [var_upd] , ctl_dict['vlevels'] )

             local_nz = len( obs_record_list )  #This is the number of levels for one single variable. 

             print('Obs_var_records ',obs_record_list)
             print('Update_var_records ',upd_record_list)

             print( obs_increment )

             [mean_diff , mode_diff , std_diff , kld , updated_mean_kf , updated_std_kf , updated_mode_kf ]=comm.compare_update( obs_var = my_ensemble[:,:,obs_record_list,:] , 
                                                       state_var = my_ensemble[:,:,upd_record_list,:] , undefmask = my_undefmask[:,:,obs_record_list] ,
                                                       nens = nbv , nx = nx , ny = ny , nz = local_nz , obs_increment = obs_increment[iv] , obs_error = obs_error[iv] )

             my_file=basedir + '/' + my_exp_name + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/update_comp_meandiff_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '.grd'
             comm.write_data(outfile=my_file,mydata=mean_diff,nx=nx,ny=ny,nz=local_nz,ie=endian,acc=access)

             my_file=basedir + '/' + my_exp_name + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/update_comp_modediff_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '.grd'
             comm.write_data(outfile=my_file,mydata=mode_diff,nx=nx,ny=ny,nz=local_nz,ie=endian,acc=access)

             my_file=basedir + '/' + my_exp_name + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/update_comp_vdiff_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '.grd'
             comm.write_data(outfile=my_file,mydata=std_diff,nx=nx,ny=ny,nz=local_nz,ie=endian,acc=access)

             my_file=basedir + '/' + my_exp_name + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/update_comp_kld_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '.grd'
             comm.write_data(outfile=my_file,mydata=kld,nx=nx,ny=ny,nz=local_nz,ie=endian,acc=access)

             my_file=basedir + '/' + my_exp_name + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/update_comp_updated_mean_kf_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '.grd'
             comm.write_data(outfile=my_file,mydata=updated_mean_kf,nx=nx,ny=ny,nz=local_nz,ie=endian,acc=access)

             my_file=basedir + '/' + my_exp_name + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/update_comp_updated_std_kf_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '.grd'
             comm.write_data(outfile=my_file,mydata=updated_std_kf,nx=nx,ny=ny,nz=local_nz,ie=endian,acc=access)

             my_file=basedir + '/' + my_exp_name + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/update_comp_updated_mode_kf_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '.grd'
             comm.write_data(outfile=my_file,mydata=updated_mode_kf,nx=nx,ny=ny,nz=local_nz,ie=endian,acc=access)

        ctime = ctime + delta

      print ( "Finish time loop" )



