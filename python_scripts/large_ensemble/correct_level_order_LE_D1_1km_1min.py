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

from common_functions import common_functions as comm


basedir='/home/ra001011/a03471/data/output_data/'

expname = 'LE_D1_1km_1min'

filetypes=['analgp']   #analgp , analgz , guesgp , guesgz

smooth=False             #False- no smooth , True apply smooth
smooth_lambda=10         #Smooth length scale (in number of grid points)

nbv=1000                 #Total number of ensemble members.

nmoments=4               #Cantidad total de momentos que vamos a calcular.

get_kldistance=True      #Wether we compute or not the Kullback-Leiber distance.

get_histogram=True       #Wether if histogram will be explicitelly calculated and stored.

get_moments=True         #Wether we will be computeing ensemble moments.

itime = dt.datetime(2013,7,13,5,0,0)  #Initial time.
etime = dt.datetime(2013,7,13,5,17,0)  #End time.

ctime = itime
#Define the delta.
delta=dt.timedelta(seconds=60)  #Original data is every 30 seconds

#=========================================================
#  LOOP OVER FILE TYPES
#=========================================================


for my_file_type in filetypes  :

  #=========================================================
  #  READ CTL FILE
  #=========================================================

  ctl_file = basedir + expname + '/ctl/' + my_file_type + '.ctl' 

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

  #All the fields will be read so..
  n_selected_fields = nz
  selected_fields = np.arange(0,nz) + 1


  #=========================================================
  #  READ LAT LON
  #=========================================================

  tmp=bio.read_data_direct(basedir + expname + '/latlon/latlon.grd',nx,ny,2,dtypein)
  lat=tmp[:,:,1]
  lon=tmp[:,:,0]

  #=========================================================
  #  START COMPUTATION LOOP
  #=========================================================

  my_ensemble=np.zeros([nx,ny,nz,nbv]).astype('float32')
  #comm.allocate_ensemble(nx=nx,ny=ny,nz=nz,nbv=nbv) #Allocate the ensemble and the undefmask.



  while ( ctime <= etime ):


    print(ctime)

    print(' Reading the ensemble ')
 
    my_path=basedir + '/' + expname + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/' 

    if not os.path.exists(my_path)                       :
       print('[Warning]:Path ' + my_path + ' does not exist: SKIP IT')

    else                                                 :
  
      #Read the ensemble and the undef mask
      [my_ensemble , my_undefmask] = comm.read_ensemble(path=my_path,nx=nx,ny=ny,nbv=nbv
                                                        ,selected_fields=selected_fields
                                                        ,n_selected_fields=n_selected_fields
                                                        ,undef=undef,ie=endian,acc=access)

      if get_moments :
          print("Computing pdf moments" + ctime.strftime("%Y%m%d%H%M%S") )
          #Get the moments of the PDF for different variables and levels.
          my_moments=comm.compute_moments(my_ensemble=my_ensemble,my_undefmask=my_undefmask 
                                        ,nx=nx,ny=ny,nz=nz,nbv=nbv,nmoments=nmoments,undef=undef) 

          for imoment in range (0,nmoments)  :
            momentstr="%04d" % ( imoment + 1 )
            my_file=basedir + '/' + expname + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/moment' + momentstr + '.grd'
            bio.write_data_direct_woundef(my_file,my_moments[:,:,:,imoment],'f4')

      if get_kldistance    :
           print("Computing kl divergence" + ctime.strftime("%Y%m%d%H%M%S") ) 
           kldist=comm.compute_kld(my_ensemble=my_ensemble,my_undefmask=my_undefmask
                                  ,nx=nx,ny=ny,nz=nz,nbv=nbv,undef=undef,nbins=nbins)
           my_file=basedir + '/' + expname + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/kldistance.grd'
           bio.write_data_direct_woundef(my_file,kldist,'f4')

      if get_histogram     :
           print("Computing histogram" + ctime.strftime("%Y%m%d%H%M%S") )
       #Compute histogram explicitelly
           [ensmin,ensmax,histogram]=comm.compute_histogram(my_ensemble=my_ensemble,my_undefmask=my_undefmask
                                                         ,nx=nx,ny=ny,nz=nz,nbv=nbv,nbins=nbins,undef=undef)
           #Write ensemble minimum
           my_file=basedir + '/' + expname + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/ensmin.grd' 
           bio.write_data_direct_woundef(my_file,ensmin,'f4')
           #Write ensemble maximum
           my_file=basedir + '/' + expname + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/ensmax.grd'
           bio.write_data_direct_woundef(my_file,ensmax,'f4')
           #Write histogram (check!)
           my_file=basedir + '/' + expname + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/histogram.grd'
           histogram=np.reshape(histogram,[nx,ny,nz*nbins])
           bio.write_data_direct_woundef(my_file,ensmax,'i2') #Using low precission integer to store histogram.

    ctime = ctime + delta

print ( "Finish time loop" )



