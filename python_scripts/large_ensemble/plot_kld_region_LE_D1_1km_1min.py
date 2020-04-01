# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016

@author:
"""
import sys
sys.path.append('../../common_python/common_functions/')
sys.path.append('../../common_python/common_modules/')

import numpy as np
import datetime as dt
import ctl_reader as ctlr
import os

import common_plot_functions as cpf
import common_mask_functions as cmf

basedir='/home/ra001011/a03471/data/output_data/'

expname = '/LE_D1_1km_1min/'

filetypes=['guesgp']   #analgp , analgz , guesgp , guesgz

max_dbz_levels=[-55,-50,30]   #Levels for vertically averaged condensate.

plot_pressure_levels=np.array([600,500,300])   #Which levels will be plotted

plot_variables=['w','tk','qv','u','v']


#Define initial and end times using datetime module.
itime = dt.datetime(2013,7,13,5, 0,0)  #Initial time.
etime = dt.datetime(2013,7,13,5,25,0)  #End time.

#Define the delta.
delta=dt.timedelta(seconds=60)

#Compute the total number of times
ntimes=int( 1 + np.around((etime-itime).seconds / delta.seconds) ) #Total number of times.

#=========================================================
#  DEFINE REGION/S
#=========================================================

reg=dict()

reg['CONV']=dict()
reg['CONV']['ini_lat']=34.60
reg['CONV']['end_lat']=35.4
reg['CONV']['ini_lon']=135
reg['CONV']['end_lon']=136.3


#=========================================================
#  LOOP OVER FILE TYPES
#=========================================================

for my_file_type in filetypes     :

  ctl_file = basedir + expname + '/ctl/' + my_file_type + '.ctl'

  plotbasedir=basedir + expname + '/plots/' + my_file_type + '/'

  #=========================================================
  #  READ CTL FILE
  #=========================================================

  ctl_dict = ctlr.read_ctl( ctl_file )

  nx=ctl_dict['nx']
  ny=ctl_dict['nx']
  nlev=ctl_dict['nz']
  nt=int(1)             #Force the number of times to be one.
  ctl_dict['nt']=int(1) #Force the number of times to be one.

  undef=np.float32( ctl_dict['undef'] )

  plotlevels=list()
  for ilev in plot_pressure_levels  :
      tmp=np.nonzero( ilev == ctl_dict['vlevels'])  
      if np.size( tmp ) > 0  :
         plotlevels.append(tmp[0])
      else                   :
         print('Warning! Level ' + str(ilev) + ' not found in this dataset')
   
  plotlevels=np.array(plotlevels)  #From list to numpy array.

  #=========================================================
  #  READ LAT LON
  #=========================================================

  latlon_file = basedir + expname + '/latlon/latlon.grd'

  tmp=ctlr.read_data_records(latlon_file,ctl=ctl_dict,records=np.array([0,1]))
  lat=tmp[:,:,1]
  lon=tmp[:,:,0]

  #=========================================================
  #  DEFINE REGIONS
  #=========================================================

  lati=np.array([34.75,34.6])
  late=np.array([35.25,34.9])
  loni=np.array([135.5,135.4])
  lone=np.array([136.25,135.7])

  #=========================================================
  #  DEFINE VARIABLES
  #=========================================================

  kld=dict()

  ens_mean=dict()


  kld_time_mean=dict()
  kld_time_std=dict()

  max_dbz=np.zeros([nx,ny,nlev])

  time_mean_max_dbz=np.zeros([nx,ny,nlev])

  my_kld=np.zeros([nx,ny,nlev])

  my_kld_sq=np.zeros([nx,ny,nlev])


  for var in plot_variables        :

    if var in ctl_dict['var_list']  :

       varpos=[i for i,x in enumerate(ctl_dict['var_list']) if x == var][0]
       varsize=int( ctl_dict['var_size'][varpos] )
       if varsize == 0   :
          varsize = 1
       kld_time_mean[var]=np.zeros([nx,ny,varsize])
       kld_time_std[var]=np.zeros([nx,ny,varsize])

  #=========================================================
  #  START TIME LOOP
  #=========================================================

  it=0

  ctime=itime
 
  while ( ctime <= etime ):

   print( ctime )

   print ( 'Reading the kld ')

   #=========================================================
   #  READ THE DATA
   #=========================================================

   my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ my_file_type + '/' + '/kldistance.grd'

   kld=ctlr.read_data_grads(my_file,ctl_dict,masked=False)

   #Read the ensemble mean to get the information from the storm location.
   my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ my_file_type + '/' + '/moment0001.grd'

   ens_mean=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
   
   #Compute max_dbz (we will use this to identify areas associated with clouds and convection)
   tmp_max_dbz = np.squeeze( np.nanmax(ens_mean['dbz'],2) )
   
   for ilev in range(0,nlev)   :          #Create a fake 3D array for the max_dbz
                                        #This is because the plotting function expects a 3D array as input.
      max_dbz[:,:,ilev]=tmp_max_dbz


   #=======================================================================================
   #Plot KLD
   #=======================================================================================

   # for key in kld :

   for var in plot_variables        :

      if var in ctl_dict['var_list']  :

        for my_reg in reg               :

         #Plot moments.
         my_kld=kld[var]
         my_kld[ my_kld > 100 ] = np.nan #There are som inf values in the reflectivity field.
         my_kld[ my_kld == undef ] = np.nan

         print('Kld for Var ',var,' ',(np.nanmin(my_kld)),np.nanmax(my_kld))

         date=ctime.strftime("%Y%m%d%H%M%S")
         cpf.set_default()  #Restore defaults
         my_map=cpf.cmap_discretize('Blues',10)
         cpf.figconf['figpath']=plotbasedir
         cpf.figconf['figsize']=(12,6)
         cpf.figconf['titlefontsize']=20
         cpf.figconf['labelfontsize']=20
         cpf.figconf['pcolor']=True
         cpf.figconf['shadedmin']=0.0
         if var == 'w'     :
            cpf.figconf['shadedmax']=1.0
         else              :
            cpf.figconf['shadedmax']=0.15
         cpf.figconf['shadedcolormap']=my_map
         cpf.figconf['colorbar']=True
         cpf.figconf['colorbarfontsize']=15
         cpf.figconf['axessize']=[0.1,0.1,0.8,0.8]
         cpf.figconf['contour']=True
         cpf.figconf['contourlevels']=max_dbz_levels
         cpf.figconf['contourcolormap']='Reds'
         cpf.figconf['gridline']=True
         cpf.figconf['ylabel']='Lat'
         cpf.figconf['xlabel']='Lon'
         cpf.figconf['xtick']=[134.5,135,135.5,136,136.5,137]
         cpf.figconf['ytick']=[34,34.5,35,35.5]
         cpf.figconf['axesrange']=[reg[my_reg]['ini_lon'],reg[my_reg]['end_lon'],reg[my_reg]['ini_lat'],reg[my_reg]['end_lat']]
         #cpf.figconf['show']=True
 
         if my_kld.shape[2] > 1 :

             for il in plotlevels   : 
                cpf.figconf['title']='KLD for variable ' + var + ' at level ' + str( int(ctl_dict['vlevels'][np.asscalar(il)])) + ' at ' + date
                cpf.figconf['figname']='/Figure_KLD_' + my_reg + '_' + var + '_' + date + '_' + str(int(ctl_dict['vlevels'][np.asscalar(il)])) 

                cpf.plot_x_y_cartopy( lon , lat , my_kld[:,:,il] , max_dbz[:,:,il] , my_kld[:,:,il] )

         else                    :
                cpf.figconf['title']='KLD for variable ' + var + ' at ' + date 
                cpf.figconf['figname']='/Figure_KLD_' + var + '_' + date    
     
                cpf.plot_x_y_cartopy( lon , lat , my_kld[:,:] , max_dbz[:,:,0] , my_kld[:,:] )

   #=========================================================================================
   #Advance time
   #=========================================================================================

   ctime = ctime + delta
 
   it = it + 1

  print ( "Finish time loop" )

  ntimes=it




