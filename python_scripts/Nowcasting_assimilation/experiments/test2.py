# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 17:36:06 2017

@author: RISDA 2018
"""

#Run a LETKF experiment using the observations created by the script run_nature.py

import sys
sys.path.append('../model/')
sys.path.append('../data_assimilation/')

#from model  import lorenzn          as model          #Import the model (fortran routines)
#from obsope import common_obs       as hoperator      #Import the observation operator (fortran routines)
from da_2d import common_da_tools_2d as das            #Import the data assimilation routines (fortran routines)
import common_obs_lorenzN_2D as hoperator
import genadvback_radar as gen
import gen_pert as pert
import gen_start_field as start_f

import matplotlib.pyplot as plt
import numpy as np
import time
import conf_experiment_mv as conf         #Load the experiment configuration
from scipy import stats
import os
import generate_obs as ge
from datetime import datetime
from datetime import timedelta 

#=================================================================
# LOAD CONFIGURATION : 
#=================================================================

GeneralConf=conf.GeneralConf
DAConf     =conf.DAConf
ModelConf  =conf.ModelConf

#=================================================================
#  LOAD OBSERVATIONS AND NATURE RUN CONFIGURATION
#=================================================================

field_start_0 = np.load(GeneralConf['DataBase']+'/'+GeneralConf['radar']+'/'+GeneralConf['date']+'/'+GeneralConf['vol']+'/'+GeneralConf['path_ref']+'/cappi_2km_resh2km_'+GeneralConf['date_start_0']+'.npy')

field_start_1 = np.load(GeneralConf['DataBase']+'/'+GeneralConf['radar']+'/'+GeneralConf['date']+'/'+GeneralConf['vol']+'/'+GeneralConf['path_ref']+'/cappi_2km_resh2km_'+GeneralConf['date_start_1']+'.npy')

fileoutdir = GeneralConf['DataBase']+'/'+GeneralConf['radar']+'/'+GeneralConf['date']+'/'+GeneralConf['vol']+'/'+ GeneralConf['filedir_mv']

date_start = GeneralConf['date_start_1']

dt = ModelConf['dt']
dx = ModelConf['dx']
desp_max = ModelConf['desp_max']
sigma = ModelConf['sigma']
sigmatrend = ModelConf['sigmatrend']
nx = ModelConf['nx']
ny = ModelConf['ny']
box_size = ModelConf['box_size']
lat = np.load(GeneralConf['DataBase']+'/'+GeneralConf['radar']+'/'+GeneralConf['date']+'/'+GeneralConf['vol']+'/'+GeneralConf['path_ref']+'/lat_lon/latitud_anguil.npy')
lon = np.load(GeneralConf['DataBase']+'/'+GeneralConf['radar']+'/'+GeneralConf['date']+'/'+GeneralConf['vol']+'/'+GeneralConf['path_ref']+'/lat_lon/longitud_anguil.npy')

ObsLoc = np.zeros((lat.size, 3))
ObsLoc[:,0] = lat.flatten()
ObsLoc[:,1] = lon.flatten()
ObsLoc[:,2] = np.ones(lat.size)
ObsError = np.random.rand(ModelConf['Nx'])
XLoc = np.zeros((lat.size,2))
XLoc[:,0] = lat.flatten()
XLoc[:,1] = lon.flatten()
dt_entreimag = GeneralConf['dt_entreimag']

#U_x, V_x = start_f.gen_start_field(field_start_0, field_start_1, fileoutdir, date_start, dt, dx, desp_max, sigma, sigmatrend, nx, ny, box_size, UNDEF=-999)

U_x = np.zeros( nx * ny )
V_x = np.zeros( nx * ny )

#=================================================================
# INITIALIZATION : 
#=================================================================

Nx = ModelConf['Nx']    ####nx*ny#####
nvar = ModelConf['nvar']
nens = DAConf['nens']
DALength = 1
sigma_start = DAConf['sigma_start']
desv_gauss_start = DAConf['desv_gauss_start']

XF = np.zeros([Nx, nens, nvar, DALength])


#Generate a random initial conditions
XF[:,:,:,0] = pert.gen_pertfield_ini(nx, ny, U_x, V_x, sigma_start, desv_gauss_start, nens)


iref=np.round( nx/2 )*nx + np.round( ny/2)

cov_U=np.zeros( nx * ny )

for ii in range( Nx )  :

    cov_U[ii] = np.cov( [ XF[iref,:,0,0] , XF[ii,:,0,0] ]  )[0,1]


        

cov_U=np.reshape( cov_U , [nx,ny])


plt.pcolor( cov_U )
plt.colorbar()
plt.show()
plt.close()




