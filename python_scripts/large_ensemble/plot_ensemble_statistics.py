# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016

@author:
"""
import numpy as np
#import matplotlib as plt
import datetime as dt
#import binary_io as bio
#import bred_vector_functions as bvf
import os
import common_functions

print(dir(common_functions))

basedir='/data9/jruiz/EXPERIMENTS/'

expname = '/OsakaPAR_1km_control1000m_smallrandompert_new/'

filetype='analgp'   #analgp , analgz , guesgp , guesgz

compute_statistics=True

undef=1.0e33

plotbasedir=basedir + expname + '/plots/'

enssize=10       #Total number of bred vectors.

nmoments=4       #Cantidad total de momentos que vamos a calcular.

kldistance=True  #Wether we compute or not the Kullback-Leiber distance.

plotlevels=np.array([6,13,17])   #Which levels will be plotted.
plotvars=['U','V','W','T','QV','QHYD']    #Which variables will be plotted.

common_functions.com_mean(3,plotlevels)

#The following will be used to extract a particlar variable from the original data.
#This variables should be especified according to the data that we have in the binary files.
ctl_vars=['U','V','W','T','QV','QHYD']  #Complete list of variables in ctl file.
ctl_inirecord=[0,12,24,36,48,60]        #Starting record for each variable. From 0 to N
ctl_endrecord=[11,23,35,47,59,71]       #End record for each variable. From 0 to N.

#Create the plotbasedir
if not os.path.exists(plotbasedir):
   os.mkdir(plotbasedir)

#Defini initial and end times using datetime module.
itime = dt.datetime(2013,7,13,5,10,00)  #Initial time.
etime = dt.datetime(2013,7,13,5,11,00)  #End time.

#Define the delta.
delta=dt.timedelta(seconds=60)  #Original data is every 30 seconds 

nx=180
ny=180
nz=np.max(ctl_endrecord) + 1  #Total number of records in binary file.
nlev=12                       #Number of vertical levels for 3D variables.

ctime= itime + delta

#Get lat lon.

lat=bio.read_data_direct(basedir + expname + '/latlon/lat.grd',nx,ny,1,'>f4')[:,:,0]
lon=bio.read_data_direct(basedir + expname + '/latlon/lon.grd',nx,ny,1,'>f4')[:,:,0]

my_ensemble=np.zeros([nx,ny,nz,enssize])

while ( ctime <= etime ):

  #Ensemble reading loop.
  if compute_statistics  :

   for imem in range (0, enssize):

      memstr="%04d" % ( imem + 1 )

      print( ' Reading ensemble member ' + memstr )
 
      my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/' + filetype + '/' + memstr + '.grd'
      my_ensemble[:,:,:,imem-1]=bio.read_data_direct_woundef(my_file,nx,ny,nz,'f4')

   undefmask=np.any( abs(my_ensemble) > undef , axis = 3 )
    
   #Todo, sacar el parametro undefbin porque solo sirve para el ouput. Ver como enmascaro los nans.
   my_moments=comf.compute_moments(my_ensemble,nx,ny,nz,enssize,1,undefmask) 

   #for imoment in range (0,nmoments)  :
   #   momentstr="%04d" % ( imoment + 1 )
   #   my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/' + filetype + '/moment' + momentstr + '.grd'
   #   write_data_direct_woundef(my_file,my_moments[:,:,:,imoment],'f4')

   #kldist=comf.com_kld(ndim, var1, kld)

   plt.figure(1)
   plt.pcolor(my_ensemble[:,:,0,0])
   plt.show()

   plt.figure(1)
   plt.pcolor(my_moments[:,:,0,0])
   plt.show()

   #Escribo los estadisticos calculados en un archivo.
  
  #Fin de la lectura de los miembros del ensamble y del calculo de los estadisticos
  else   :
   
   #Leo los estadisticos ya calculados de un archivo.

  #Continuo con el ploteo.  

   #Todo ... in this section compute, plot and write distribution centered moments.

   #Todo compute K-L distance based on the ensemble.

   #Compute histogram and bimodality ?

   ctime = ctime + delta

print ( "Finish time loop" )





