# -*- coding: utf-8 -*-

#Grafica el histograma de cada variable dentro de un determinado dominio para el 
#punto con el maximo indice de kld. E indica la posicion de dicho punto.
#Esto sirve para tener una idea visual de que significa un indice de kld grande para 
#diferentes variables y casos.

"""
Created on Tue Nov  1 18:45:15 2016

@author:
"""

import sys
sys.path.append('../../../common_python/common_functions/')
sys.path.append('../../../common_python/common_modules/')


import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import ctl_reader as ctlr

import common_histogram_functions as chf
import common_mask_functions as cmf

import scipy.ndimage.filters as spyf
import os

nbv=1000                                      #Number of ensemble members
nbins=int(np.round(np.sqrt(nbv)))             #Number of bins used to compute the histogram
thresholdmin=0.005   
smooth_range=2                                #To smooth the histogram

dbz_threshold = 30.0                          #Only values with max dbz above this will be used

basedir='/home/jruiz/share/LARGE_ENSEMBLE/output_data/home/ra001011/a03471/data/output_data/'

expname = '/LE_D1_1km_5min/'


filetypes=['guesgp']   #analgp , analgz , guesgp , guesgz

plot_variables=['dbz','v','w','tk']

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

ctl_file = basedir + expname + '/ctl/guesgp.ctl'

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
radar_range=45.0e3   #Radar range in meters (to define the radar mask)
                     #I intentionally using a smaller range to discard grid points
                     #near the edge of the radar coverage.
#Exclude areas outside the radar domain.

radar_mask = cmf.distance_range_mask( lon_radar , lat_radar , radar_range , lon , lat )


fig, axs = plt.subplots( 2,2 , figsize=[10,10] )

fig.subplots_adjust(wspace=0.15,hspace=0.2,bottom=0.05,left=0.05,right=0.97,top=0.97)

icol=0
irow=0

delta = dt.timedelta(seconds=300)
#Search for the histogram associated to the maximum KLD

for ivar , my_var in enumerate( plot_variables ) :
    max_kld = 0.0
    ctime=itime

    while ( ctime <= etime )  :
       date=ctime.strftime("%Y%m%d%H%M%S")
       my_file=basedir + '/'+expname+'/' + date + '/guesgp/moment0001.grd'
       ens_mean=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
       my_file=basedir + '/'+expname+'/' + date + '/guesgp/moment0002.grd'
       ens_var=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
       my_file=basedir + '/'+expname+'/' + date + '/guesgp/kldistance.grd'
       ens_kld=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
       dbz_mask = ens_mean['dbz'] < dbz_threshold
       ens_kld[my_var][ dbz_mask ] = np.nan
       ens_kld[my_var]=np.delete( ens_kld[my_var] , 4 , 2 )
       nz = np.shape(ens_kld[my_var])[2]
       #Apply the radar mask, we will only consider grid points within the radar domain.
       for iz in range(0,nz)   :
          ens_kld[my_var][:,:,iz][np.logical_not(radar_mask)]=np.nan
          ens_kld[my_var][ ens_kld[my_var] == undef ]=np.nan
       print(np.nanmax( ens_kld[my_var] ) )
       if np.nanmax( ens_kld[my_var] ) > max_kld :
          max_kld = np.nanmax( ens_kld[my_var] )
          max_loc=np.nanargmax( ens_kld[my_var] )
          [tmpnx , tmpny , tmpnz ]=np.shape( np.squeeze( ens_kld[my_var] ) )
          [xmax  , ymax  , zmax  ]=np.unravel_index( max_loc , (tmpnx,tmpny,tmpnz) )
          hist_file = basedir + '/'+expname+'/' + date + '/guesgp/histogram.grd'
          max_file  = basedir + '/'+expname+'/' + date + '/guesgp/ensmax.grd'
          min_file  = basedir + '/'+expname+'/' + date + '/guesgp/ensmin.grd'
          hist=chf.read_histogram(hist_file,max_file,min_file,nx,ny,nbins,ctl_dict,dtypein='i2') 
          hist_max=np.squeeze( np.delete(hist[my_var]['maxval'],4,2)[xmax,ymax,zmax] )
          hist_min=np.squeeze( np.delete(hist[my_var]['minval'],4,2)[xmax,ymax,zmax] )
          hist_var=np.squeeze( np.delete(ens_var[my_var],4,2)[xmax,ymax,zmax] )
          hist_mean=np.squeeze( np.delete(ens_mean[my_var],4,2)[xmax,ymax,zmax] )
          hist_delta=(hist_max-hist_min)/nbins
          hist_range=hist_min + hist_delta / 2 + hist_delta *  np.arange(0,nbins,1)
          print('hist shape ',hist[my_var]['hist'].shape)
          hist_bars=np.delete(hist[my_var]['hist'],4,2)[xmax,ymax,zmax,:] / np.sum( np.delete(hist[my_var]['hist'],4,2)[xmax,ymax,zmax,:] )
          #Smooth the histogram
          hist_bars_s=spyf.uniform_filter(hist_bars, size=smooth_range, mode='reflect')
          #Get the Gaussian fit to the histogram
          hist_range_ext=hist_min + hist_delta / 2 + hist_delta *  np.arange(-10,nbins+10,1)
          hist_gaussfit=hist_delta*(1/(np.sqrt(2*np.pi*hist_var)))*np.exp( -0.5*np.power(hist_range_ext-hist_mean,2)/hist_var )
          #Make the plot and save it 
          hist_kld = 100 * max_kld
          hist_label= '(KLD=' + '%1.1f' % hist_kld + ')'

       ctime = ctime + delta

    if ivar == 0 :
      title = '(a) Ref. ' + hist_label
      ax=axs[0,0]
      xamin=45.0
      xamax=65.0
      deltaa=5.0
    if ivar == 1 :
      title = '(b) U ' + hist_label
      ax=axs[0,1]
      xamin=4.0
      xamax=15.0
      deltaa=4.0
    if ivar == 2 :
      title = '(c) W ' + hist_label
      ax=axs[1,0]
      xamin=-5.0
      xamax=15.0
      deltaa=5.0
    if ivar == 3 :
      title = '(d) T ' + hist_label
      ax=axs[1,1]
      xamin=255.0
      xamax=265.0
      deltaa=5.0

    the_bars=ax.bar( hist_range , hist_bars_s , width = hist_delta ,color='grey' )
    the_lines=ax.plot( hist_range_ext , hist_gaussfit , 'k--',linewidth=4)
    ax.set_xlim(xmin=xamin, xmax=xamax )
    ax.set_ylim(ymin=0.0, ymax=0.4 )
    ax.set_xticks( np.arange( xamin , xamax+deltaa, deltaa) )
    #ax.set_xticklabels( np.round(1.0e2*hist_range[::5])/1.0e2 ,fontsize=12)
    ax.set_title(title,fontsize=14)
    print('dBZ',ens_mean['dbz'][xmax,ymax,zmax])
    #ax.text( np.max(hist_range) - 30*hist_delta , 0.09 ,title ,fontsize=15,color='k',bbox={'facecolor':'white', 'alpha':0.8,'edgecolor':'white'})
    #ax.text( np.max(hist_range) - 23*hist_delta , 0.09 ,var_name + hist_label ,fontsize=15,color='k',bbox={'facecolor':'white', 'alpha':0.8,'edgecolor':'white'})



#plt.show()
my_fig = './Figure_examples_largekld_histogram.png'
plt.savefig( my_fig , format='png' , dpi=300)
plt.close()
        




