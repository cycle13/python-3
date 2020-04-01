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
from da_prueba_2d import common_da_tools_2d as das            #Import the data assimilation routines (fortran routines)
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

#field_start_0 = np.load(GeneralConf['DataBase']+'/'+GeneralConf['radar']+'/'+GeneralConf['date']+'/'+GeneralConf['vol']+'/'+GeneralConf['path_ref']+'/cappi_2km_resh2km_'+GeneralConf['date_start_0']+'.npy')

#field_start_1 = np.load(GeneralConf['DataBase']+'/'+GeneralConf['radar']+'/'+GeneralConf['date']+'/'+GeneralConf['vol']+'/'+GeneralConf['path_ref']+'/cappi_2km_resh2km_'+GeneralConf['date_start_1']+'.npy')

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
#lat = np.load(GeneralConf['DataBase']+'/'+GeneralConf['radar']+'/'+GeneralConf['date']+'/'+GeneralConf['vol']+'/'+GeneralConf['path_ref']+'/lat_lon/latitud_anguil.npy')
#lon = np.load(GeneralConf['DataBase']+'/'+GeneralConf['radar']+'/'+GeneralConf['date']+'/'+GeneralConf['vol']+'/'+GeneralConf['path_ref']+'/lat_lon/longitud_anguil.npy')

#Genero lats y lons simples 

lat=np.zeros((nx,ny))
lon=np.zeros((nx,ny))
for ii in range(nx):
    for jj in range(ny):
        lat[ii,jj]=jj/100
        lon[ii,jj]=ii/100
        
dx = 1.0/100.0
dy = 1.0/100.0
        


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

U_x = np.ones( nx * ny )
V_x = np.ones( nx * ny )*-1

#=================================================================
# INITIALIZATION : 
#=================================================================

Nx = ModelConf['Nx']    ####nx*ny#####
nvar = ModelConf['nvar']
nens = DAConf['nens']
DALength = DAConf['dalengh']
sigma_start = 0.5
desv_gauss_start = 10

XA = np.zeros([Nx, nens, nvar, DALength])
XF = np.zeros([Nx, nens, nvar, DALength])


#Generate a random initial conditions
XF[:,:,:,0] = pert.gen_pertfield_ini(nx, ny, U_x, V_x, sigma_start, desv_gauss_start, nens)


#UObs, VObs = ge.generate_obs(field_0, field_1, fileoutdir, idatename1, dt, dx, desp_max, sigma, sigmatrend, nx, ny, box_size, UNDEF=-999)

UObs= np.ones( nx * ny ) * -999.0
VObs= np.ones( nx * ny ) * -999.0

UObs[ int(np.round(ny/2)*nx + np.round(nx/2) ) ]=5
VObs[ int(np.round(ny/2)*nx + np.round(nx/2) ) ]=-5

print( np.max(UObs) , np.min(UObs) , np.max(VObs) , np.min(VObs) )

YObs = np.concatenate((UObs,VObs))

   
UF, VF, UObs, VObs, ObsLoc, ObsErrorW = hoperator.model_to_obs(nens,XF[:,:,0,0], XF[:,:,1,0], UObs, VObs, ObsLoc, ObsError)


#Pasar toodos los ensambles que estan en XF
YObs = np.concatenate((UObs, VObs))
YF = np.concatenate((UF,VF))
ObsLoc_da = np.vstack([ObsLoc,ObsLoc])
NObs = YObs.size   
ObsErrorW = np.concatenate((ObsErrorW, ObsErrorW))

print(ObsLoc)

plt.pcolor(  np.reshape( XF[:,3,0,0]  , [nx,ny] ) )
plt.colorbar()
plt.show()
plt.close()

XA[:,:,:,0] = das.da_letkf( nx=Nx , nt=1 , no=NObs , nens=nens ,  xloc=XLoc  ,
                              tloc=np.zeros(Nx)+1.0        , nvar=2                        , xfens=XF[:,:,:,0]           ,
                              obs=YObs             , obsloc=ObsLoc_da                , ofens=YF                       ,
                              rdiag=ObsErrorW       , loc_scale=np.array([30000.0,-1.0]) , inf_coefs=DAConf['InfCoefs']   ,
                              update_smooth_coef=0.0 , dx=dx , dy=dy )[:,:,:,0]  


xa_mean=np.mean( XA[:,:,:,0] , 1 )
xf_mean=np.mean( XF[:,:,:,0] , 1 )

xa_mean_u=np.reshape( xa_mean[:,0] , [nx,ny] )
xf_mean_u=np.reshape( xf_mean[:,0] , [nx,ny] )


plt.pcolor( xa_mean_u  )
plt.colorbar()
plt.show()
plt.close()

xa_mean_v=np.reshape( xa_mean[:,1] , [nx,ny] )
xf_mean_v=np.reshape( xf_mean[:,1] , [nx,ny] )


plt.pcolor( xa_mean_v )
plt.colorbar()
plt.show()
plt.close()



