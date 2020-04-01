# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016

@author:
"""
import sys
sys.path.append('../../common_python/common_covmatrix/')
sys.path.append('../../common_python/common_modules/')

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import ctl_reader as ctlr
import binary_io  as bio
import os

from covariance_matrix_tools import covariance_matrix_tools as cmt

#To compute some components of the ensemble covariance matrix.
config=dict()

config['datapath']='/home/ra001011/a03471/data/output_data/'      #path to the experiment data.
config['experiment']='LE_D1_5km'                                  #Experiment name
 
config['filetype']='guesgp'                                       #analgp,analgz,guesgp,guesgz... etc
config['np']=2                                                    #Number of processors for openmp
config['enssize']=1000                                            #Number of ensemble members
config['bootstrap_samples']=0                                     #Number of bootstrap samples (0 means no bootstrap)

config['skipx']=1                                                 #Skip values for covariance strength computation.
config['skipy']=1                                                 # 1 means no skip.
config['skipz']=1

config['delta']=40                                                #Maximum distance for covariance computation strhength.

config['smoothcov']=False                                         #Wheter fields will be smoothed before computing covariances.
config['smoothcovlength']=20.0e4                                  #Lanczos filter length for field smoothing.
config['dx']=1.0e3                                                #Horizontal resolution of the data.

config['basevars']=['dbz','u','v','tk']                           #The covariance between each basevar and each covvar will be computed.
config['covvars']=['dbz','u','v','tk']
config['levels']=[1000,900,500,200]                               #Select the vertical levels to be used in the computation.

config['itime'] = dt.datetime(2013,7,12,14,00,00)  #Initial time.
config['etime'] = dt.datetime(2013,7,12,14,00,00)  #End time.

#Define the delta.
config['delta']=dt.timedelta(seconds=3600)  #Original data is every 30 seconds

ctime = config['itime']

#=========================================================
#  PREPARE SOME ENVIROMENT VARIABLES
#=========================================================

os.system('ulimit -s unlimited')
os.system('export OMP_STACKSIZE=500M')
os.system('export KMP_STACKSIZE=500M')
os.system('export OMP_NUM_THREADS=' + str( config['np'] ) )

#=========================================================
#  READ CTL FILE
#=========================================================

ctl_file = config['datapath'] + config['experiment'] + '/ctl/' + config['filetype'] + '.ctl' 

ctl_dict = ctlr.read_ctl( ctl_file )

nx=np.int(ctl_dict['nx'])
ny=np.int(ctl_dict['ny'])
nz=np.array(ctl_dict['end_record']).max() + 1 #Total number of records in binary file.
nlev=ctl_dict['nz']                 #Number of vertical levels for 3D variables.
undef=ctl_dict['undef']
if ctl_dict['big_endian']  :
   endian='little_endian'
else                       :
   endian='big_endian'

#=========================================================
#  DETERMINE WICH RECORDS WILL BE RETAINED.
#=========================================================

[subset_var_list , subset_lev_list , subset_record_list ] =  ctlr.record_subset( ctl_dict , np.unique( config['basevars'] + config['covvars']) , config['levels'])

subset_record_list = subset_record_list + 1  #1 based instead of 0 based

subset_nrecord=np.size( subset_record_list )

#=========================================================
#  START THE LOOP IN TIME
#=========================================================

#while ( ctime <= etime ):

  #=========================================================
  #  READ THE ENSEMBLE
  #=========================================================
  my_path= config['datapath'] + config['experiment'] + '/' + datstr  + '/' + config['filetype'] + '/'
 
  [my_ensemble,undef_mask]=cmt.read_ensemble(path=my_path,nx=nx,ny=ny,nbv=enssize,selected_fields=subset_record_list,n_selected_fields=subset_nrecord,undef=undef,ie=endian)

  #if get_moments :
  #    print("Computing pdf moments")
  #    my_moments=comm.compute_moments(nx=nx,ny=ny,nz=nz,nbv=enssize,nmoments=nmoments,undef=undef) 

  #    for imoment in range (0,nmoments)  :
  #      momentstr="%04d" % ( imoment + 1 )
  #      my_file=basedir + '/' + expname + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + filetype + '/moment' + momentstr + '.grd'
  #      bio.write_data_direct_woundef(my_file,my_moments[:,:,:,imoment],'f4')

  #if get_kldistance    :
  #     print("Computeing kl distance") 
  #     kldist=comm.compute_kld(nx=nx,ny=ny,nz=nz,nbv=enssize,undef=undef,nbins=nbins)
  #     my_file=basedir + '/' + expname + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + filetype + '/kldistance.grd'
  #     bio.write_data_direct_woundef(my_file,kldist,'f4')

  #if get_histogram     :
  #     print("Computing histogram")
  # #Compute histogram explicitelly
  #     ensmin,ensmax,histogram=comm.compute_histogram(nx=nx,ny=ny,nz=nz,nbv=enssize,nbins=nbins,undef=undef)
  #     #Write ensemble minimum
  #     my_file=basedir + '/' + expname + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + filetype + '/ensmin.grd' 
  #     bio.write_data_direct_woundef(my_file,ensmin,'f4')
  #     #Write ensemble maximum
  #     my_file=basedir + '/' + expname + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + filetype + '/ensmax.grd'
  #     bio.write_data_direct_woundef(my_file,ensmax,'f4')
  #     #Write histogram (check!)
  #     my_file=basedir + '/' + expname + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + filetype + '/histogram.grd'
  #     histogram=np.reshape(histogram,[nx,ny,nz*nbins])
  #     bio.write_data_direct_woundef(my_file,ensmax,'i2') #Using low precission integer to store histogram.
  #
  #ctime = ctime + delta

print ( "Finish time loop" )



