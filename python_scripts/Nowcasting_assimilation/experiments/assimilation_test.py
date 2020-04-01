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
#import common_obs_lorenzN_2D as hoperator
import common_obs_mv as hoperator
import gen_pert as pert
import gen_start_field as start_f

import matplotlib.pyplot as plt
import numpy as np
import time
import conf_experiment_mv as conf         #Load the experiment configuration
import os
import generate_obs as ge
from datetime import datetime
from datetime import timedelta 
import common_functions as cf

#=================================================================
# LOAD CONFIGURATION : 
#=================================================================

GeneralConf=conf.GeneralConf
DAConf     =conf.DAConf
ModelConf  =conf.ModelConf

inidate =  datetime.strptime(GeneralConf['date_start_mv'],'%Y%m%d%H%M%S')
enddate =  datetime.strptime(GeneralConf['date_end_mv'  ],'%Y%m%d%H%M%S')

fileoutdir = '../tmp/'

#=================================================================
#  GET MOTION VECTOR FIELD TO PRODUCE THE FIRST FORECAST.
#=================================================================

dt = ModelConf['dt']
dx = ModelConf['dx']
desp_max = ModelConf['desp_max']
sigma = ModelConf['sigma']
sigmatrend = ModelConf['sigmatrend']
nx = ModelConf['nx']
ny = ModelConf['ny']
box_size = ModelConf['box_size']
lat = np.load('../datos/latitud_anguil.npy')
lon = np.load('../datos/longitud_anguil.npy')

format_date=GeneralConf['format_date']

#ObsError = np.random.rand(ModelConf['Nx'])
XLoc = np.zeros((lat.size,2))
XLoc[:,0] = lat.flatten()
XLoc[:,1] = lon.flatten()
dt_entreimag = GeneralConf['dt_entreimag']

file_0 = cf.Data_File_Name( GeneralConf , inidate )
field_0 = np.load(file_0)

file_1 = cf.Data_File_Name( GeneralConf , inidate + timedelta(seconds=GeneralConf['dt_entreimag']) )
field_1 = np.load(file_1)

if os.path.isfile(fileoutdir + '/UVStart.npz') :

   tmp = np.load( fileoutdir + '/UVStart.npz' ) 

   U_x = tmp['U_x']
   V_x = tmp['V_x']  

else                               :

   U_x, V_x = start_f.gen_start_field(field_0, field_1, fileoutdir, inidate , dt_entreimag, dx, desp_max, sigma, sigmatrend, nx, ny, box_size, UNDEF=-999)

   np.savez( fileoutdir + '/UVStart.npz' , U_x = U_x , V_x = V_x )


field_0[field_0<0]=np.nan
#cf.Plot_Data_Fast_2d( GeneralConf , field_0 )
field_1[field_1<0]=np.nan
#cf.Plot_Data_Fast_2d( GeneralConf , field_1 )

#cf.Plot_Data_Fast( GeneralConf , U_x , nx , ny )

#cf.Plot_Data_Fast( GeneralConf , V_x , nx , ny )
cf.Plot_Ref_Wind_Fast( GeneralConf , field_0 , U_x , V_x , figname= fileoutdir + '/WindRef_0.png')


#=================================================================
# INITIALIZATION : 
#=================================================================

Nx = ModelConf['Nx']    ####nx*ny#####
nvar = ModelConf['nvar']
nens = DAConf['nens']
DALength = DAConf['dalengh']
sigma_start = DAConf['sigma_start']
desv_gauss_start = DAConf['desv_gauss_start']

XA = np.zeros([Nx, nens, nvar, DALength])
XF = np.zeros([Nx, nens, nvar, DALength])

#Generate a random initial conditions
XA[:,:,:,0] = pert.gen_pertfield_ini(nx, ny, U_x, V_x, sigma_start, desv_gauss_start, nens)

cf.Plot_Data_Fast( GeneralConf , np.mean( XA[:,:,0,0] , 1 ) , nx , ny , figname= fileoutdir + '/XaMean_0.png')
cf.Plot_Data_Fast( GeneralConf , XA[:,0,0,0] - np.mean( XA[:,:,0,0] , 1 ) , nx , ny , figname= fileoutdir + '/XaPert_0.png')

#=================================================================
#  MAIN DATA ASSIMILATION LOOP : 
#=================================================================
sigma_da = DAConf['sigma_da']
sigma_gauss_da = DAConf['sigma_gauss_da']


cdate = inidate

it = 1 

#Adelanto el tiempo antes de iniciar el ciclo.
cdate = cdate + 2 * timedelta( seconds = dt_entreimag )

#for it in range(1, DALength):
while cdate <  enddate            :
 
   print('Data assimilation cycle # ', str(it) )

   #=================================================================
   #  ENSEMBLE FORECAST  : 
   #=================================================================   

   print('Generate observations')
   #Las observaciones del k-esimo ciclo son los motion vectors obtenidos en los
   #dt_entreimag segundos previos al tiempo del analisis. 

   date_0 = cdate - timedelta( seconds=dt_entreimag )
   date_1 = cdate 

   #El pronostico advecta los datos desde el tiempo k-esimo - 1 hasta el tiempo
   #k-esimo en donde se combinaran con las observaciones.

   if os.path.isfile( fileoutdir + '/UVObs' + str(it) + '.npz') :
  
      print('Observations already generated, load them from a file')

      tmp = np.load( fileoutdir + '/UVObs' + str(it) + '.npz' )

      UObs = tmp['UObs']
      VObs = tmp['VObs'] 

   else                                            :
 
      print('Observations not present, lets generate them ')

      date_0_str = date_0.strftime(format_date)
      date_1_str = date_1.strftime(format_date)

      if not (os.path.exists('../datos/cappi_2km_resh2km_'+ date_0_str + '.npy') and os.path.exists('../datos/cappi_2km_resh2km_'+ date_1_str + '.npy')):
         print("No se encuentran los archivos "+ date_0_str +" o "+ date_1_str )
         cdate = cdate+timedelta(minutes=dt_entreimag)
         continue

      field_0 = np.load('../datos/cappi_2km_resh2km_'+ date_0_str +'.npy')
      field_1 = np.load('../datos/cappi_2km_resh2km_'+ date_1_str +'.npy')

      UObs, VObs = ge.generate_obs(field_0, field_1, fileoutdir, date_0_str , dt_entreimag , dx, desp_max, sigma, sigmatrend, nx, ny, box_size, UNDEF=-999)

      np.savez(  fileoutdir + '/UVObs' + str(it) + '.npz' , UObs=UObs , VObs=VObs )

   print('Propagating the ensemble in time (persistance so far)')
   #######   UXF, VXF = gen.genadvback_radar(nx, ny, u_motion1, v_motion1, ddt, code)  ########       
   XF[:,:,:,it] = XA[:,:,:,it-1]  #-> Persistance
   XF[:,:,:,it] = pert.gen_pertfield(nx, ny, XF[:,:,:,it], sigma_da, sigma_gauss_da, nens)

   plot_UObs = np.copy( UObs )
   plot_VObs = np.copy( VObs )
   plot_UObs[ plot_UObs == -999 ] = np.nan
   plot_VObs[ plot_VObs == -999 ] = np.nan
   cf.Plot_Data_Fast( GeneralConf , plot_UObs , nx , ny , figname= fileoutdir + '/UObs_' + str(it) + '.png')
   cf.Plot_Data_Fast( GeneralConf , plot_VObs , nx , ny , figname= fileoutdir + '/VObs_' + str(it) + '.png')

   #=================================================================
   #  OBSERVATION OPERATOR  : 
   #================================================================= 

   #Apply h operator and transform from model space to observation space. 
   print("Run observation operator")

   ObsLoc = np.zeros((lat.size, 3))
   ObsLoc[:,0] = lat.flatten()
   ObsLoc[:,1] = lon.flatten()
   ObsLoc[:,2] = np.ones(lat.size) #*it 

   #El ObsError es la varianza del error aca lo corregi porque en su lugar habia un ruido random (eso es si queremos simular errores de observacion)
   #pero en este caso como las observaciones son reales los motion vectors estimados ya tienen un error (lo que ponemos aca es lo que suponemos es la 
   #varianza de ese error que no conocemos)

   ObsError = np.ones( ModelConf['Nx'] )

   UF, VF, UObs, VObs, ObsLoc, ObsErrorW = hoperator.model_to_obs(nens,XF[:,:,0,it], XF[:,:,1,it], UObs, VObs, ObsLoc, ObsError)

   #UF es el ensamble de U en el espacio de las observaciones.
   #VF es el ensamble de V en el espacio de las observaciones.
   #UObs es el vector que contiene solo las observaciones que no son undef.
   #VObs idem anteior pero para V.
   #ObsLoc es la ubicacion la ubicacion [lat,lon,tiempo]

   #Pasar toodos los ensambles que estan en XF
   YObs = np.concatenate((UObs, VObs))
   YF = np.concatenate((UF,VF))
   ObsLoc_da = np.vstack([ObsLoc,ObsLoc])
   NObs = YObs.size   
   ObsErrorW = np.concatenate((ObsErrorW, ObsErrorW))
   #print('Observation operator took ', time.time()-start, 'seconds.')

   #=================================================================
   #  LETKF DA  : 
   #================================================================= 

   print('Data assimilation')

   #STATE VARIABLES ESTIMATION:
 
   XA[:,:,:,it] = das.da_letkf( nx=Nx , nt=1 , no=NObs , nens=nens ,  xloc=XLoc  ,
                              tloc=np.ones(1)        , nvar=2                        , xfens=XF[:,:,:,it]           ,
                              obs=YObs             , obsloc=ObsLoc_da                , ofens=YF                       ,
                              rdiag=ObsErrorW       , loc_scale=DAConf['LocScales'] , inf_coefs=DAConf['InfCoefs']   ,
                              update_smooth_coef=0.0 )[:,:,:,0]

   #XA[:,:,:,it] = XF[:,:,:,it]

   cf.Plot_Data_Fast( GeneralConf , np.mean( XA[:,:,0,it] , 1 ) , nx , ny , figname= fileoutdir + '/XaMean_'+str(it)+'.png')
   cf.Plot_Data_Fast( GeneralConf , np.mean( XF[:,:,0,it] , 1 ) , nx , ny , figname= fileoutdir + '/XfMean_'+str(it)+'.png')

   cf.Plot_Data_Fast( GeneralConf , XA[:,0,0,0] - np.mean( XA[:,:,0,it] , 1 ) , nx , ny , figname= fileoutdir + '/XaPert_'+str(it)+'.png')
   cf.Plot_Data_Fast( GeneralConf , XF[:,0,0,0] - np.mean( XF[:,:,0,it] , 1 ) , nx , ny , figname= fileoutdir + '/XfPert_'+str(it)+'.png')


   cf.Plot_Data_Fast( GeneralConf , np.mean( XA[:,:,0,it] - XF[:,:,0,it] , 1 ) , nx , ny , figname= fileoutdir + '/XaUpdate_'+str(it)+'.png')

   it = it+1
   cdate = cdate + timedelta(seconds=dt_entreimag)

np.save(fileoutdir +'/asimilacion.npy', XA=XA,XF=XF)
