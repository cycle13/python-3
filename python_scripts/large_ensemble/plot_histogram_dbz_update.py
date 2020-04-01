# -*- coding: utf-8 -*-

#Grafica el histograma de cada variable dentro de un determinado dominio para el 
#punto con el maximo indice de bimodalidad. E indica la posicion de dicho punto.
#Esto sirve para tener una idea visual de que significa un indice de bimodalidad grande para 
#diferentes variables y casos.

"""
Created on Tue Nov  1 18:45:15 2016

@author:
"""

import sys
sys.path.append('../../common_python/common_functions/')
sys.path.append('../../common_python/common_modules/')


import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import binary_io as bio
import ctl_reader as ctlr

import common_histogram_functions as chf
import common_plot_functions as cpf
import common_mask_functions as cmf

import scipy.ndimage.filters as spyf


import bred_vector_functions as bvf  #Candidato a salir.
import os

nbv=1000                                      #Number of ensemble members
nbins=int(np.round(np.sqrt(nbv)))             #Number of bins used to compute the histogram
thresholdmin=0.005   
smooth_range=2                                #To smooth the histogram

dbz_threshold = 30.0                          #Only values with max dbz above this will be used

basedir='/home/ra001011/a03471/data/output_data/'

expname = '/LE_D1_1km_30sec/'

plotbasedir=basedir + expname + '/plots/'

filetypes=['guesgp']   #analgp , analgz , guesgp , guesgz

plot_variables=['u','v','w','tk','qv','dbz']


#Define initial and end times using datetime module.
itime = dt.datetime(2013,7,13,5,5,0)  #Initial time.
etime = dt.datetime(2013,7,13,5,55,0)  #End time.

#Define the delta.
delta=dt.timedelta(seconds=300)

#Compute the total number of times
ntimes=int( 1 + np.around((etime-itime).seconds / delta.seconds) ) #Total number of times.
times=np.zeros(ntimes)

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

  #=========================================================
  #  READ LAT LON
  #=========================================================

  latlon_file = basedir + expname + '/latlon/latlon.grd'

  tmp=ctlr.read_data_records(latlon_file,ctl=ctl_dict,records=np.array([0,1]))
  lat=np.squeeze(tmp[:,:,1])
  lon=np.squeeze(tmp[:,:,0])

  #=========================================================
  #  DEFINE A RADAR MASK (GRID POINTS WITHIN THE RADAR RANGE)
  #=========================================================

  lat_radar=34.823
  lon_radar=135.523
  radar_range=50.0e3   #Radar range in meters (to define the radar mask)
                       #I intentionally using a smaller range to discard grid points
                       #near the edge of the radar coverage.
  #Exclude areas outside the radar domain.
  radar_mask = cmf.distance_range_mask( lon_radar , lat_radar , radar_range , lon , lat )

  #=========================================================
  #  LOOP OVER TIME
  #=========================================================
  #Define dictionaries.
  ctime=itime 

  it=0

  while ( ctime <= etime ):

    print( ctime )

    print ( 'Reading the histogram ')

    hist_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/' + '/histogram.grd'
    max_file =basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/' + '/ensmax.grd'
    min_file =basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/' + '/ensmin.grd'

    hist_forecast=chf.read_histogram(hist_file,max_file,min_file,nx,ny,nbins,ctl_dict,dtypein='i2')

    hist_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/analgp/' + '/histogram.grd'
    max_file =basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/analgp/' + '/ensmax.grd'
    min_file =basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/analgp/' + '/ensmin.grd'

    hist_analysis=chf.read_histogram(hist_file,max_file,min_file,nx,ny,nbins,ctl_dict,dtypein='i2')


    # hist_properties=analyze_histogram_fun( my_hist , thresholdmin )
   
    my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/' + '/moment0001.grd'
    ens_mean_forecast=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=True)
    my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/analgp/' + '/moment0001.grd'
    ens_mean_analysis=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=True)

    my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/' + '/moment0002.grd'
    ens_var_forecast=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=True)
    my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/analgp/' + '/moment0002.grd'
    ens_var_analysis=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=True)

    my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/' + '/kldistance.grd'
    ens_kld_forecast=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=True)
    my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/analgp/' + '/kldistance.grd'
    ens_kld_analysis=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=True)

    var='qv'
    #Apply the radar mask, we will only consider grid points within the radar domain.
    ref_update = ens_mean_analysis[var] - ens_mean_forecast[var]
    for iz in range(0,ref_update.shape[2])   :
       ref_update[:,:,iz][np.logical_not(radar_mask)]=np.nan
       if iz == 4 :
          ref_update[:,:,iz]=np.nan

    ref_update[ ens_mean_analysis[var] > 1.0e-4 ]=np.nan
   
    #Busco el punto donde mas se reduce la reflectividad. 
    max_loc=np.nanargmin( ref_update )

    [tmpnx , tmpny , tmpnz ]=np.shape( np.squeeze( ref_update ) )
    [xmax  , ymax  , zmax  ]=np.unravel_index( max_loc , (tmpnx,tmpny,tmpnz) ) 

    print(xmax,ymax,zmax)
    #Get the histogram limits for the corresponding maximum bimodality index location.
    hist_max=np.squeeze( hist_forecast[var]['maxval'][xmax,ymax,zmax] )
    hist_min=np.squeeze( hist_forecast[var]['minval'][xmax,ymax,zmax] )
    hist_delta=(hist_max-hist_min)/nbins
    hist_range=hist_min + hist_delta / 2 + hist_delta *  np.arange(0,nbins,1)
    hist_bars=hist_forecast[var]['hist'][xmax,ymax,zmax,:] / np.sum( hist_forecast[var]['hist'][xmax,ymax,zmax,:] )
    mean=ens_mean_forecast[var][xmax,ymax,zmax]
    vvar=ens_var_forecast[var][xmax,ymax,zmax] 
    kld=ens_kld_forecast[var][xmax,ymax,zmax]
    
    plt.figure()
    plt.subplot(1,2,1)     
    #Make the plot and save it 
    the_bars=plt.bar( hist_range , hist_bars , width = hist_delta ,color='r' ) 
    plt.title('Forecast: Var ' + str(vvar) + ' Mean ' + str(mean) + ' KLD ' + str(kld))
    #plt.show()

    #Get the histogram limits for the corresponding maximum bimodality index location.
    hist_max=np.squeeze( hist_analysis[var]['maxval'][xmax,ymax,zmax] )
    hist_min=np.squeeze( hist_analysis[var]['minval'][xmax,ymax,zmax] )
    hist_delta=(hist_max-hist_min)/nbins
    hist_range=hist_min + hist_delta / 2 + hist_delta *  np.arange(0,nbins,1)
    hist_bars=hist_analysis[var]['hist'][xmax,ymax,zmax,:] / np.sum( hist_analysis[var]['hist'][xmax,ymax,zmax,:] )
    mean=ens_mean_analysis[var][xmax,ymax,zmax]
    vvar=ens_var_analysis[var][xmax,ymax,zmax]
    kld=ens_kld_analysis[var][xmax,ymax,zmax]

    plt.subplot(1,2,2)
    #Make the plot and save it 
    the_bars=plt.bar( hist_range , hist_bars , width = hist_delta ,color='r' )
    plt.title('Analysis: Var ' + str(vvar) + ' Mean ' + str(mean) + ' KLD ' + str(kld))
    #plt.show()

    my_fig = basedir + expname  + '/plots/' + my_file_type + '/Figure_histogram_max_bimodality_' + var + '_' + ctime.strftime("%Y%m%d%H%M%S") + '.png'
    print( 'Generationg the following figure : ' + my_fig )
    plt.show()
    plt.savefig( my_fig )
    plt.close()

    ctime = ctime + delta
 
  it = it + 1

  print ( "Finish time loop" )


