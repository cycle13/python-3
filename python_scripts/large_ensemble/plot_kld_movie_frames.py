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
import matplotlib.pyplot as plt
import os

import common_plot_functions as cpf
import common_mask_functions as cmf

basedir='/home/ra001011/a03471/data/output_data/'

#expname = '/LE_D1_1km_5min/'

filetypes=['guesgp']   #analgp , analgz , guesgp , guesgz

max_dbz_levels=[-55,-50,30]   #Levels for vertically averaged condensate.

plot_pressure_levels=np.array([500,600])   #Which levels will be plotted

plot_variables=['qv','tk']


#Define initial and end times using datetime module.
itime = dt.datetime(2013,7,13,5, 0,0)  #Initial time.
etime = dt.datetime(2013,7,13,5,20,0)  #End time.

#Exp names
exp_names=['LE_D1_1km_5min','LE_D1_1km_2min','LE_D1_1km_1min','LE_D1_1km_30sec']

#Define the delta.
delta=[ dt.timedelta(seconds=300) , dt.timedelta(seconds=120) , dt.timedelta(seconds=60) , dt.timedelta(seconds=30) ]

nexp=np.size( exp_names )

#Compute the total number of times
ntimes=int( 1 + np.around((etime-itime).seconds / delta[3].seconds) ) #Total number of times.

#=========================================================
#  DEFINE REGION/S
#=========================================================

reg=dict()

reg['CONV']=dict()
reg['CONV']['ini_lat']=34.60
reg['CONV']['end_lat']=35.4
reg['CONV']['ini_lon']=135
reg['CONV']['end_lon']=136.3

#reg['TOTAL']=dict()
#reg['TOTAL']['ini_lat']=34.2
#reg['TOTAL']['end_lat']=35.5
#reg['TOTAL']['ini_lon']=134.75
#reg['TOTAL']['end_lon']=136.25


#=========================================================
#  LOOP OVER FILE TYPES
#=========================================================

for my_file_type in filetypes     :

  ctl_file = basedir + exp_names[0] + '/ctl/' + my_file_type + '.ctl'

  plotbasedir=basedir + '/plots/' + my_file_type + '/'

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

  latlon_file = basedir + exp_names[0] + '/latlon/latlon.grd'

  tmp=ctlr.read_data_records(latlon_file,ctl=ctl_dict,records=np.array([0,1]))
  lat=tmp[:,:,1]
  lon=tmp[:,:,0]

  #=========================================================
  #  DEFINE VARIABLES
  #=========================================================

  kld=dict()

  ens_mean=dict()

  max_dbz=np.zeros([nx,ny,nexp])

  for var in plot_variables        :

    if var in ctl_dict['var_list']  :
 
       varpos=[i for i,x in enumerate(ctl_dict['var_list']) if x == var][0]
       varsize=int( ctl_dict['var_size'][varpos] )
       if varsize == 0   :
          varsize = 1
       kld[var]=np.zeros([nx,ny,varsize,nexp])

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

   iexp = 0
   for my_exp in exp_names  : 

      #I need to round the time to the closest time for each experiment...
      time_diff= ( ctime - itime ).total_seconds()
     
      time_diff =  np.round( time_diff / delta[iexp].total_seconds() ) * delta[iexp].total_seconds() 

      my_time = itime + dt.timedelta( seconds = time_diff ) 

      print( my_exp , my_time.strftime("%Y%m%d%H%M%S") , ctime.strftime("%Y%m%d%H%M%S") )

      #Read KLD
      my_file=basedir + my_exp + '/' + my_time.strftime("%Y%m%d%H%M%S") + '/'+ my_file_type + '/' + '/kldistance.grd'
      print(my_file)

      my_kld=ctlr.read_data_grads(my_file,ctl_dict,masked=False)

      for var in plot_variables  :
         if var in my_kld           :
            kld[var][:,:,:,iexp] = np.squeeze( my_kld[var] )
 
      my_file=basedir + my_exp + '/' + my_time.strftime("%Y%m%d%H%M%S") + '/'+ my_file_type + '/' + '/moment0001.grd'
      my_mean=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
  

      #Compute max_dbz (we will use this to identify areas associated with clouds and convection)
      max_dbz[:,:,iexp] = np.squeeze( np.nanmax(my_mean['dbz'],2) )

      iexp = iexp + 1
   #=======================================================================================
   #Plot KLD
   #=======================================================================================

   # for key in kld :

   for var in plot_variables        :

      if var in ctl_dict['var_list']  :

         #Plot moments.
         my_kld=kld[var]
         my_kld[ my_kld > 100 ] = np.nan #There are som inf values in the reflectivity field.
         my_kld[ my_kld == undef ] = np.nan

         print('Kld for Var ',var,' ',(np.nanmin(my_kld)),np.nanmax(my_kld))

         date=ctime.strftime("%Y%m%d%H%M%S")
         cpf.set_default()  #Restore defaults
         my_map=cpf.cmap_discretize('Blues',10)
         cpf.figconf['figpath']=plotbasedir
         cpf.figconf['figsize']=(12,10)
         cpf.figconf['titlefontsize']=20
         cpf.figconf['labelfontsize']=12
         cpf.figconf['pcolor']=True
         cpf.figconf['shadedmin']=0.0
         if var == 'w'     :
            cpf.figconf['shadedmax']=2.0
         else              :
            cpf.figconf['shadedmax']=0.15
         cpf.figconf['shadedcolormap']=my_map
         cpf.figconf['colorbar']=False
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
         cpf.figconf['nsubplotrows']=2
         cpf.figconf['nsubplotcolumns']=2
         cpf.figconf['subplottitles']=['5 min.','2 min.','1 min.','30 sec.']
         cpf.figconf['text']=True
         cpf.figconf['textstring']=ctime.strftime("%H:%M:%S") + ' UTC'
         cpf.figconf['textlocx']=0.5
         cpf.figconf['textlocy']=0.485


         #cpf.figconf['show']=True

         for il in plotlevels   : 
             for my_reg in reg      :
               #Select the appropiate region
               cpf.figconf['axesrange']=[reg[my_reg]['ini_lon'],reg[my_reg]['end_lon'],reg[my_reg]['ini_lat'],reg[my_reg]['end_lat']]

               cpf.figconf['figname']='/Figure_KLD_' + my_reg + '_' + var + '_' + date + '_' + str(int(ctl_dict['vlevels'][np.asscalar(il)])) 
               cpf.plot_x_y_subplot( lon , lat , np.squeeze( my_kld[:,:,il,:] ) , max_dbz[:,:,:] , np.squeeze( my_kld[:,:,il,:] ) )

   #=========================================================================================
   #Advance time
   #=========================================================================================

   ctime = ctime + delta[3]
 
   it = it + 1

  print ( "Finish time loop" )

  ntimes=it





