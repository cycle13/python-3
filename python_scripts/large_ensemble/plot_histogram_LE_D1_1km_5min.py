# -*- coding: utf-8 -*-

#Grafica el histograma de cada variable dentro de un determinado dominio para el 
#punto con el maximo kld. E indica la posicion de dicho punto.
#Esto sirve para tener una idea visual de que significa un KLD grande para 
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

basedir='/home/ra001011/a03471/data/output_data/'

expname = '/LE_D1_1km_5min/'

plotbasedir=basedir + expname + '/plots/'

filetypes=['guesgp']   #analgp , analgz , guesgp , guesgz

plot_variables=['u','v','w','tk','qv','dbz']


#Define initial and end times using datetime module.
itime = dt.datetime(2013,7,13,5,15,0)  #Initial time.
etime = dt.datetime(2013,7,13,5,15,0)  #End time.

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

    hist_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ my_file_type + '/' + '/histogram.grd'
    max_file =basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ my_file_type + '/' + '/ensmax.grd'
    min_file =basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ my_file_type + '/' + '/ensmin.grd'

    hist=chf.read_histogram(hist_file,max_file,min_file,nx,ny,nbins,ctl_dict,dtypein='i2')

    # hist_properties=analyze_histogram_fun( my_hist , thresholdmin )
   
    my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ my_file_type + '/' + '/moment0001.grd'
    ens_mean=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
    my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ my_file_type + '/' + '/moment0002.grd'
    ens_var=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
    #my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ my_file_type + '/' + '/moment0003.grd'
    #ens_skew=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
    #my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ my_file_type + '/' + '/moment0004.grd'
    #ens_kurt=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
    my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ my_file_type + '/' + '/kldistance.grd'
    ens_kld=ctlr.read_data_grads(my_file,ctl_dict,masked=False)


    #Get the location of the maximum kld for the selected variables.
    for var in plot_variables    :
       if var in ens_kld     :
         if ( np.size( np.shape( ens_kld[var] ) ) >= 2 )  :
            nz = np.shape(ens_kld[var])[2] 
         else                                             :
            nz = 1
 
         #Apply the radar mask, we will only consider grid points within the radar domain.
         for iz in range(0,nz)   :
            ens_kld[var][:,:,iz][radar_mask]=np.nan
         ens_kld[var][ ens_kld[var] > 100 ]=np.nan
         ens_kld[var][ ens_kld[var] == undef ]=np.nan
 
         max_loc=np.nanargmax( ens_kld[var] )

         [tmpnx , tmpny , tmpnz ]=np.shape( np.squeeze( ens_kld[var] ) )
         [xmax  , ymax  , zmax  ]=np.unravel_index( max_loc , (tmpnx,tmpny,tmpnz) ) 

         #Get the histogram limits for the corresponding maximum kld location.
         hist_max=np.squeeze( hist[var]['maxval'][xmax,ymax,zmax] )
         hist_min=np.squeeze( hist[var]['minval'][xmax,ymax,zmax] )
         hist_delta=(hist_max-hist_min)/nbins
         hist_range=hist_min + hist_delta / 2 + hist_delta *  np.arange(0,nbins,1)
         hist_bars=hist[var]['hist'][xmax,ymax,zmax,:] / np.sum( hist[var]['hist'][xmax,ymax,zmax,:] )
         hist_mean=ens_mean[var][xmax,ymax,zmax]
         hist_std=np.sqrt( ens_var[var][xmax,ymax,zmax] )
         #hist_skew=ens_skew[var][xmax,ymax,zmax]/np.power(ens_var[var][xmax,ymax,zmax],3/2)
         #hist_kurt=ens_kurt[var][xmax,ymax,zmax]/np.power(ens_var[var][xmax,ymax,zmax],2) -3
         hist_kld=ens_kld[var][xmax,ymax,zmax]

         #Smooth the histogram
         hist_bars_s=spyf.uniform_filter(hist_bars, size=smooth_range, mode='reflect')
         #Get the Gaussian fit to the histogram
         hist_gaussfit=hist_delta*(1/(np.sqrt(2*np.pi*np.power(hist_std,2))))*np.exp( -0.5*np.power( (hist_range-hist_mean)/(hist_std) ,2) )
         
         #Make the plot and save it 
         hist_label= var + ' KLD=' +  str( hist_kld )
         the_bars=plt.bar( hist_range , hist_bars_s , width = hist_delta ,color='r',label=hist_label ) 
         the_lines=plt.plot( hist_range , hist_gaussfit , 'k--',linewidth=4) 
         plt.legend(handles=[the_bars])
         print( 'Generationg the following figure : ' + 'Figure_histogram_maxkld_' + var + '_' + ctime.strftime("%Y%m%d%H%M%S") + '.png' )
         plt.savefig( basedir + expname  + '/plots/' + my_file_type + '/Figure_histogram_maxkld_' + var + '_' + ctime.strftime("%Y%m%d%H%M%S") + '.png' )
         #plt.show()
         plt.close()



    ctime = ctime + delta
 
  it = it + 1

  print ( "Finish time loop" )


