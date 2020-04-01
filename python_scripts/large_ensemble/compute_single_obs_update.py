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
expnames  = ['LE_D1_1km_5min']

obs_variables=['dbz','u','v']    #These are the observed variables (all variables will be updated)

#variable_combination=[['dbz','dbz','dbz','dbz','dbz','v' ,'v' ,'v','v' ],
#                      ['tk' ,'v'  ,'qv' ,'w'  ,'dbz','tk','qv','w','v' ]]

#expdeltas = [300,300,300,300,300,300]
expdeltas = [120]

obs_location=[115,67,8]     #Max W KLD location 5MIN  (I took into account the missing vertical level)
#obs_location=[114,68,7]    #Max W KLD location 30SEC

#delta=dt.timedelta(seconds=60)  #Original data is every 30 seconds

filetypes=['guesgp']   #Only guess will be processed.

nbv=1000              #Total number of ensemble members.

obs_increment = [-5.0 , 2.0 , 2.0 ]
obs_error     = [5.0 , 2.0 , 2.0 ]

itime = dt.datetime(2013,7,13,5,30,0)  #Initial time.
etime = dt.datetime(2013,7,13,5,30,0)  #End time.

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

          for iobs , my_obs_var in enumerate( obs_variables ) :
 
             print('Computing single obs update for variable ',my_var 
             [i_ind,e_ind]=get_var_start_end( ctl , my_var )  #Get start and end index.

             hxf = my_ensemble[obs_location[0],obs_location[1],obs_location[2]+i_ind,:]
             yo  = np.mean( hxf ) + obs_increment[iobs]
             
             #EnKF update
             [umean,ustd,umode]compute_single_obs_update(hxf=hxf,ens=my_ensemble,mask=my_undefmask,nx=nx,ny=ny,nz=nz,nens=nbv,yo=yo,obs_error=obs_error[iobs]
                                                         ,undef=undef,update_type=1]
             prefix=basedir + '/' + my_exp_name + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/updateEnKF_' + my_obs_var + '_'
             sufix = '_obsloc_x'+str(obs_location[0])+'_y'+str(obs_location[1])+'_z'+str(obs_location[2])+'.grd'
             my_file= prefix + 'mean' + sufix 
             comm.write_data(outfile=my_file,mydata=umean,nx=nx,ny=ny,nz=nz,ie=endian,acc=access)
             my_file= prefix + 'std' + sufix
             comm.write_data(outfile=my_file,mydata=ustd,nx=nx,ny=ny,nz=nz,ie=endian,acc=access)
             my_file= prefix + 'mode' + sufix
             comm.write_data(outfile=my_file,mydata=umode,nx=nx,ny=ny,nz=nz,ie=endian,acc=access)

             #ISPF update
             [umean,ustd,umode]compute_single_obs_update(hxf=hxf,ens=my_ensemble,mask=my_undefmask,nx=nx,ny=ny,nz=nz,nens=nbv,yo=yo,obs_error=obs_error[iobs]
                                                         ,undef=undef,update_type=2]
             prefix=basedir + '/' + my_exp_name + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/updateISPF_' + my_obs_var + '_'
             sufix = '_obsloc_x'+str(obs_location[0])+'_y'+str(obs_location[1])+'_z'+str(obs_location[2])+'.grd'
             my_file= prefix + 'mean' + sufix
             comm.write_data(outfile=my_file,mydata=umean,nx=nx,ny=ny,nz=nz,ie=endian,acc=access)
             my_file= prefix + 'std' + sufix
             comm.write_data(outfile=my_file,mydata=ustd,nx=nx,ny=ny,nz=nz,ie=endian,acc=access)
             my_file= prefix + 'mode' + sufix
             comm.write_data(outfile=my_file,mydata=umode,nx=nx,ny=ny,nz=nz,ie=endian,acc=access)

        ctime = ctime + delta

      print ( "Finish time loop" )

