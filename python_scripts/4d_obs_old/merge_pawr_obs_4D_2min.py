# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016

This script merge several observation files
within a data assimilation window.

@author:
"""

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import os

import sys
sys.path.append('../../common_python/common_obs_scale/')
from common_obs_scale import common_obs_scale as oio

#####################################################################
# SETINGS 
#####################################################################

itime = dt.datetime(2013,7,13,0,0,0)  #Initial time.
etime = dt.datetime(2013,7,13,23,55,30)  #End time.

ori_obs_path='/home/ra001011/a03471/data/input_data/obs/QCED_1KM_v4_v500M_attn0.01/'  #Path to original observation files.

obs_file_prefix='radar_' #Prefix of observation file names.

out_obs_path='/home/ra001011/a03471/data/input_data/obs/QCED_1KM_v4_v500M_attn0.01_4D_2min/'  #Path were the new observation files will be stored.

ori_time_freq=30    #Original observation frequency in seconds.
out_time_freq=120   #Requested observation frequency in 4D output files.

nrecinput=7
nrecoutput=8

endian='b'          #Endianness of the input / output files.

#####################################################################

delta_t_da =dt.timedelta(seconds=out_time_freq)

delta_t_obs=dt.timedelta(seconds=ori_time_freq)

ctime = itime

while ( ctime <= etime ) :

  ctime = ctime + delta_t_da

  print('Generating obs for time: ' + ctime.strftime('%Y-%m-%d-%H-%M-%S'))

  #Search for all the files corresponding to the previous times.
  obs_time=ctime 

  while ( obs_time > ctime - delta_t_da )  :

    obs_dif=( obs_time - ctime ).total_seconds()

    print( 'Reading obs at time: ' + obs_time.strftime('%Y-%m-%d-%H-%M-%S') )

    filename =  ori_obs_path + '/' + obs_file_prefix + obs_time.strftime('%Y%m%d%H%M%S')  + '.dat' 

    [nobs,radarlon,radarlat,radarz] = oio.get_nobs_radar(cfile=filename,nrec=nrecinput,endian=endian)

    tmp_obs = oio.read_obs_radar(cfile=filename,nobs=nobs,nrec=nrecinput,endian=endian)

    tmp_ = np.ones( ( tmp_obs.shape[0] , 1 ) ) * obs_dif 

    tmp_obs = np.concatenate( (tmp_obs , tmp_) , axis = 1 )

    if obs_time == ctime    :
       obs = tmp_obs.copy()
    else                    :
       obs=np.concatenate( (obs , tmp_obs) , axis = 0 )


    obs_time = obs_time - delta_t_obs

  #End of the loop over the obs slots, now we are ready to write a new data file. 

  filename =  out_obs_path + '/' + obs_file_prefix + ctime.strftime('%Y%m%d%H%M%S')  + '.dat'

  nobs=obs.shape[0]

  if not os.path.exists(out_obs_path):
   os.mkdir(out_obs_path)

    
  oio.write_obs_radar(cfile=filename,nobs=nobs,rlon=radarlon,rlat=radarlat,rz=radarz,obs=obs,nrec=nrecoutput,append=False,endian=endian)

print('Finish!')



