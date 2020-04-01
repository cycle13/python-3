# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016

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
from scipy.stats import kurtosis
from scipy.stats import moment

import matplotlib.pyplot as plt

from common_functions import common_functions as comm


basedir='/home/ra001011/a03471/data/output_data/'

#expnames  = ['LE_D1_1km_30sec','LE_D1_1km_30sec_nospinup','LE_D1_1km_1min','LE_D1_1km_2min','LE_D1_1km_1min_4D','LE_D1_1km_5min']
expnames  = ['LE_D2_1km_5min']


#expdeltas = [30,30,60,120,60,300]
expdeltas = [300]

#delta=dt.timedelta(seconds=60)  #Original data is every 30 seconds

filetypes=['guesgp']   #analgp , analgz , guesgp , guesgz

smooth=False             #False- no smooth , True apply smooth
#smooth_lambda=10         #Smooth length scale (in number of grid points)

nbv=1000                 #Total number of ensemble members.

nmoments=4               #Cantidad total de momentos que vamos a calcular.

get_kldistance=True      #Wether we compute or not the Kullback-Leiber distance.

get_histogram=True       #Wether if histogram will be explicitelly calculated and stored.

get_moments=True         #Wether we will be computeing ensemble moments.

itime = dt.datetime(2013,7,13,5,5,0)  #Initial time.
etime = dt.datetime(2013,7,13,5,5,0)  #End time.

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
      #comm.allocate_ensemble(nx=nx,ny=ny,nz=nz,nbv=nbv) #Allocate the ensemble and the undefmask.

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

          print("Computing pdf moments fortran " + ctime.strftime("%Y%m%d%H%M%S") )
          #Get the moments of the PDF for different variables and levels.
          my_moments_fortran_100=comm.compute_moments(my_ensemble=my_ensemble[:,:,:,0:100],my_undefmask=my_undefmask 
                                        ,nx=nx,ny=ny,nz=nz,nbv=100,nmoments=nmoments,undef=undef) 

          my_moments_fortran_1000=comm.compute_moments(my_ensemble=my_ensemble[:,:,:,0:1000],my_undefmask=my_undefmask
                                        ,nx=nx,ny=ny,nz=nz,nbv=1000,nmoments=nmoments,undef=undef)


          #Gets moments scipy
          print("Computing pdf moments scipy " + ctime.strftime("%Y%m%d%H%M%S") )
          #my_moments_scipy = np.zeros( np.shape( my_moments_fortran ))
          #for imoment in range( nmoments )  :
          #   my_ensemble[my_ensemble == undef ] = np.nan
          #   my_moments_scipy[:,:,:,imoment]=moment( my_ensemble , moment=imoment+1, axis=3)

           

        ctime = ctime + delta

my_moments_fortran_100[ my_moments_fortran_100 == undef ] = np.nan
my_moments_fortran_1000[ my_moments_fortran_1000 == undef ] = np.nan

kurtosis_fortran_100= my_moments_fortran_100[:,:,:,3] / np.power( my_moments_fortran_100[:,:,:,1] , 2 ) -3

kurtosis_fortran_1000= my_moments_fortran_1000[:,:,:,3] / np.power( my_moments_fortran_1000[:,:,:,1] , 2 ) -3

my_ensemble_test = my_ensemble[99,99,29,:] 

print( kurtosis( my_ensemble_test ) , kurtosis_fortran_100[99,99,29] , kurtosis_fortran_1000[99,99,29] )


#kurtosis_scipy  = kurtosis( my_ensemble , axis=3 )



print ( "Finish time loop" )

for imoment in range( 1 ) :

   plt.figure()
   plt.subplot(1,2,1)
   p=plt.pcolor( np.squeeze( kurtosis_fortran_100[:,:,29]  ) )
   plt.colorbar(p)
   plt.subplot(1,2,2)
   p=plt.pcolor( np.squeeze( kurtosis_fortran_1000[:,:,29] ) )
   plt.colorbar(p)
   plt.show()


