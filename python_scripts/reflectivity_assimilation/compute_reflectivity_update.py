# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016

Assimilate data one observation at a time using the LETKF.
The main purpose is to make simple fast tests of observation impact using realistic prior.

@author:
"""
import sys
sys.path.append('../../common_python/common_functions/')
sys.path.append('../../common_python/common_modules/')
sys.path.append('../../common_python/common_letkf/')

import numpy as np
#import matplotlib.pyplot as plt
import datetime as dt
import ctl_reader as ctlr
import binary_io  as bio
import os
import matplotlib.pyplot as plt

from common_functions import common_functions as comm
from cletkf           import common_da        as cda


basedir='/home/jruiz/share/LARGE_ENSEMBLE/output_data/volume64/data/ra001011/pmaldonado/LETKF_SCALE/exps/'

expname  = 'CORDOBA_2km_5min_60m_cumulus'


updated_variables = ['tk','qr','qs','qg','w','p']

expdelta = 300

file_type='guesgz'     #Only guess will be processed.

nbv=30                 #Total number of ensemble members.

obs_increment = -15.0
obs_error     = 5.0 ** 2
mindbz        = 0.0e0

itime = dt.datetime(2018,11,10,21,10,0)  #Initial time.
etime = dt.datetime(2018,11,10,21,10,0)  #End time.

#=========================================================
#  LOOP OVER FILE TYPES
#=========================================================

delta=dt.timedelta( seconds = expdelta )

#=========================================================
#  READ CTL FILE
#=========================================================

ctl_file = basedir + expname + '/ctl/guesgz_mean.ctl' 

ctl_dict = ctlr.read_ctl( ctl_file )

nx=np.int(ctl_dict['nx'])
ny=np.int(ctl_dict['ny'])
nz=np.array(ctl_dict['end_record']).max() + 1 #Total number of records in binary file.
nlev=ctl_dict['nz']                           #Number of vertical levels for 3D variables.

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

ctime = itime

while ( ctime <= etime ):

   print(ctime)
   print(' Reading the ensemble ')
   my_path=basedir + '/' + expname + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + file_type + '/' 

   if not os.path.exists(my_path)                       :
      print('[Warning]:Path ' + my_path + ' does not exist: SKIP IT')
   else                                                 :
      #Read the ensemble and the undef mask
      [my_ensemble , my_undefmask] = comm.read_ensemble(path=my_path,nx=nx,ny=ny,nbv=nbv
                                                        ,selected_fields=selected_fields
                                                        ,n_selected_fields=n_selected_fields
                                                        ,undef=undef,ie=endian,acc=access)
      my_ensemble = my_ensemble.astype('float32')
      nvar = len( updated_variables )
      xf = np.zeros((nx,ny,nlev,nbv,nvar)).astype('float32')
      
      for ivar,my_var in enumerate(updated_variables)  :
          [var_list , lev_list , record_list] =ctlr.record_subset( ctl_dict , [my_var] , ctl_dict['vlevels'] )
          xf[:,:,:,:,ivar] = my_ensemble[:,:,record_list,:]
          if my_var == 'qr'  :
              iqr=ivar
          if my_var == 'qg'  :
              iqg=ivar
          if my_var == 'qs'  :
              iqs=ivar
          if my_var == 'tk'  :
              itk=ivar
          if my_var == 'p'   :
              ip=ivar


      xf_mean = np.squeeze( np.nanmean( xf , 3 ) ).astype('float32')
      xf_std  = np.squeeze( np.nanstd( xf , 3 ) ).astype('float32')

              
      #Free some memory
      my_ensemble = np.zeros((1,1))

      print('Computing prior reflectivity')
      xfo = cda.calc_ref_ens(nx=nx,ny=ny,nz=nlev,nbv=nbv,qrens=xf[:,:,:,:,iqr],qsens=xf[:,:,:,:,iqs],
                                                qgens=xf[:,:,:,:,iqg],tens=xf[:,:,:,:,itk],
                                                pens=xf[:,:,:,:,ip])

      xfo_mean = cda.calc_ref_ens(nx=nx,ny=ny,nz=nlev,nbv=1,qrens=xf_mean[:,:,:,iqr],qsens=xf_mean[:,:,:,iqs],
                                                                 qgens=xf_mean[:,:,:,iqg],tens=xf_mean[:,:,:,itk],
                                                                 pens=xf_mean[:,:,:,ip])
      #Compute the simple analysis update
      print('Computing simple letkf update')

      dep = np.ones( np.shape( np.squeeze( xfo_mean ) ) ) * obs_increment 
      error = np.ones( np.shape( xfo_mean ) ) * obs_error

      #Aplicamos el limite inferior en los valores de reflectividad.
      xfo[ xfo < mindbz ]=mindbz                         #Aplicamos la restriccion.
      tmp_mean = np.squeeze( np.mean( xfo , 3 ) )        #Calculo la media en reflectividad como la calcula el letkf.
      tmp_mask = tmp_mean + dep < mindbz 
      dep[ tmp_mask ]= mindbz - tmp_mean[ tmp_mask ]  #La observacion no puede ser menor que mindbz

      xa=cda.simple_letkf(nx=nx,ny=ny,nz=nlev,nbv=nbv,nvar=nvar,xfo=xfo,xf=xf,dep=dep,error=error,min_o=mindbz).astype('float32')

      #Compute the reflectivity asociated to the ensemble mean.
      xa_mean = np.squeeze( np.nanmean( xa , 3 ) ).astype('float32')
      xa_std  = np.squeeze( np.nanstd( xa , 3 ) ).astype('float32')

      print('Computing posterior reflectivity')
      xao = cda.calc_ref_ens(nx=nx,ny=ny,nz=nlev,nbv=nbv,qrens=xa[:,:,:,:,iqr],qsens=xa[:,:,:,:,iqs],
                                                qgens=xa[:,:,:,:,iqg],tens=xa[:,:,:,:,itk],
                                                pens=xa[:,:,:,:,ip])

      xao_mean = cda.calc_ref_ens(nx=nx,ny=ny,nz=nlev,nbv=1,qrens=xa_mean[:,:,:,iqr],qsens=xa_mean[:,:,:,iqs],
                                                                 qgens=xa_mean[:,:,:,iqg],tens=xa_mean[:,:,:,itk],
                                                                 pens=xa_mean[:,:,:,ip])

      #Compute tot cond 
      tot_cond_f = xf[:,:,:,:,iqr] + xf[:,:,:,:,iqs] + xf[:,:,:,:,iqg] 
      tot_cond_a = xa[:,:,:,:,iqr] + xa[:,:,:,:,iqs] + xa[:,:,:,:,iqg]

      tot_cond_f[ tot_cond_f < np.min(tot_cond_f) ] = np.min(tot_cond_f)
      tot_cond_a[ tot_cond_a < np.min(tot_cond_f) ] = np.min(tot_cond_f) 
      
      tot_cond_fcero_freq = np.sum( ( tot_cond_f == np.min(tot_cond_f) ).astype('int'),3)/nbv
      tot_cond_acero_freq = np.sum( ( tot_cond_a == np.min(tot_cond_f) ).astype('int'),3)/nbv


      #Write the analysis for the update variables
      print('Writing data')
      for ivar,my_var in enumerate(updated_variables)  :
          my_file=basedir + '/' + expname + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + file_type + '/update_dbz_' + my_var + '_obs_inc_' + str(obs_increment) + '_obserr_' + str(obs_error) + '.grd'
          comm.write_data(outfile=my_file,mydata=xa_mean[:,:,:,ivar],nx=nx,ny=ny,nz=nlev,ie=endian,acc=access)

      #Write the analysis reflectivity (reflectivity of the ensemble mean)
      my_file=basedir + '/' + expname + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + file_type + '/update_dbz_dbz_obs_inc_' + str(obs_increment) + '_obserr_' + str(obs_error) + '.grd'
      comm.write_data(outfile=my_file,mydata=xao_mean,nx=nx,ny=ny,nz=nlev,ie=endian,acc=access)

     #Write the tot cond cero frec
      my_file=basedir + '/' + expname + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + file_type + '/posterior_tot_cond_' + str(obs_increment) + '_obserr_' + str(obs_error) + '.grd'
      comm.write_data(outfile=my_file,mydata=np.mean(tot_cond_a,3),nx=nx,ny=ny,nz=nlev,ie=endian,acc=access)

     #Write the analysis tot cond cero frec
      my_file=basedir + '/' + expname + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + file_type + '/posterior_tot_cond_cero_frec_' + str(obs_increment) + '_obserr_' + str(obs_error) + '.grd'
      comm.write_data(outfile=my_file,mydata=tot_cond_acero_freq,nx=nx,ny=ny,nz=nlev,ie=endian,acc=access)


      #Write the gues for the update variables
      print('Writing data')
      for ivar,my_var in enumerate(updated_variables)  :
          my_file=basedir + '/' + expname + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + file_type + '/prior_' + my_var + '.grd'
          comm.write_data(outfile=my_file,mydata=xf_mean[:,:,:,ivar],nx=nx,ny=ny,nz=nlev,ie=endian,acc=access)

      #Write the analysis reflectivity (reflectivity of the ensemble mean)
      my_file=basedir + '/' + expname + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + file_type + '/prior_dbz.grd'
      comm.write_data(outfile=my_file,mydata=xfo_mean,nx=nx,ny=ny,nz=nlev,ie=endian,acc=access)

     #Write the tot cond cero frec
      my_file=basedir + '/' + expname + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + file_type + '/prior_tot_cond_' + str(obs_increment) + '_obserr_' + str(obs_error) + '.grd'
      comm.write_data(outfile=my_file,mydata=np.mean(tot_cond_f,3),nx=nx,ny=ny,nz=nlev,ie=endian,acc=access)

     #Write the analysis tot cond cero frec
      my_file=basedir + '/' + expname + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + file_type + '/prior_tot_cond_cero_frec_' + str(obs_increment) + '_obserr_' + str(obs_error) + '.grd'
      comm.write_data(outfile=my_file,mydata=tot_cond_fcero_freq,nx=nx,ny=ny,nz=nlev,ie=endian,acc=access)


   ctime = ctime + delta
print ( "Finish time loop" )



plt.figure()
plt.pcolor( np.squeeze( xfo_mean[140:190,50:110,28] ) ,cmap='gist_ncar')
plt.colorbar()
plt.title('Forecast reflectivity')

plt.figure()
plt.pcolor( np.squeeze( xao_mean[140:190,50:110,28] ) ,cmap='gist_ncar')
plt.colorbar()
plt.title('Analisys reflectivity')

plt.figure()
plt.pcolor( np.squeeze( xao_mean[140:190,50:110,28] - xfo_mean[140:190,50:110,28] ) ,cmap='gist_ncar')
plt.colorbar()
plt.title('Analisys update')

plt.figure()
plt.pcolor( np.squeeze( dep[140:190,50:110,28] ) ,cmap='Blues')
plt.colorbar()
plt.title('Observation departure')

plt.figure()
plt.pcolor( np.squeeze( np.mean(tot_cond_f[140:190,50:110,28,:],2) ) ,cmap='Blues')
plt.colorbar()
plt.title('Forecast condensate')

plt.figure()
plt.pcolor( np.squeeze( np.mean(tot_cond_a[140:190,50:110,28,:],2) ) ,cmap='Blues')
plt.colorbar()
plt.title('Analysis condensate')

plt.figure()
plt.pcolor( np.squeeze( np.mean(tot_cond_a[140:190,50:110,28,:]-tot_cond_f[140:190,50:110,28,:],2) ) ,cmap='Blues')
plt.colorbar()
plt.title('Analysis condensate')


plt.figure()
plt.pcolor( np.squeeze( tot_cond_fcero_freq[140:190,50:110,28] ) ,cmap='Blues')
plt.colorbar()
plt.title('Forecast cero condensate probability')

plt.figure()
plt.pcolor( np.squeeze( tot_cond_acero_freq[140:190,50:110,28] ) ,cmap='Blues')
plt.colorbar()
plt.title('Analysis cero condensate probability')

plt.figure()
plt.plot( xfo[160,90,28,:] )
plt.plot( xao[160,90,28,:] )

plt.figure()
plt.plot( tot_cond_f[160,90,28,:] )
plt.plot( tot_cond_a[160,90,28,:] )

