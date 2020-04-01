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

basedir='/home/ra001011/a03471/data/output_data/'

expname = '/LE_D1_1km_5min/'


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
radar_range=50.0e3   #Radar range in meters (to define the radar mask)
                     #I intentionally using a smaller range to discard grid points
                     #near the edge of the radar coverage.
#Exclude areas outside the radar domain.

radar_mask = cmf.distance_range_mask( lon_radar , lat_radar , radar_range , lon , lat )


nrows=2
ncols=2

fig, axs = plt.subplots( nrows,ncols , figsize=[10,10] )

fig.subplots_adjust(wspace=0.15,hspace=0.2,bottom=0.05,left=0.05,right=0.97,top=0.97)

icol=0
irow=0

for ii in range( 0 , 4 )  :

   

   if ii == 0 :
      print('Case 1')
      date = '20130713052000'
      var='qv'
      title='(a) - $q_v$'
      irow=0;icol=0
   if ii == 1 :
      print('Case 2')
      date = '20130713051500'
      var='tk'
      title='(b) - $T$'
      irow=0;icol=1
   if ii == 2 :
      print('Case 3')
      date = '20130713052000'
      var='w'
      title='(c) - $W$'
      irow=1;icol=0
   if ii == 3 :
      print('Case 4')
      date = '20130713053000'
      var='u'
      title='(d) - $U$'
      irow=1;icol=1

   #ax = axs[irow,icol]

   hist_file = basedir + '/LE_D1_1km_5min/' + date + '/guesgp/histogram.grd'
   max_file  = basedir + '/LE_D1_1km_5min/' + date + '/guesgp/ensmax.grd'
   min_file  = basedir + '/LE_D1_1km_5min/' + date + '/guesgp/ensmin.grd'

   hist=chf.read_histogram(hist_file,max_file,min_file,nx,ny,nbins,ctl_dict,dtypein='i2')

   my_file=basedir + '/LE_D1_1km_5min/' + date + '/guesgp/moment0001.grd'
   ens_mean=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
   my_file=basedir + '/LE_D1_1km_5min/' + date + '/guesgp/moment0002.grd'
   ens_var=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
   my_file=basedir + '/LE_D1_1km_5min/' + date + '/guesgp/bimodality_index.grd'
   ens_bimodality=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
 
   max_dbz = np.squeeze( np.nanmax( ens_mean['dbz'] , 2 ) )
   dbz_mask = max_dbz >= dbz_threshold
 
   nz = np.shape(ens_bimodality[var])[2] 
 
   #Apply the radar mask, we will only consider grid points within the radar domain.
   for iz in range(0,nz)   :
       ens_bimodality[var][:,:,iz][np.logical_not(radar_mask)]=np.nan
       ens_bimodality[var][:,:,iz][np.logical_not(dbz_mask)  ]=np.nan
       ens_bimodality[var][ ens_bimodality[var] == undef ]=np.nan
   max_loc=np.nanargmax( ens_bimodality[var] )
   [tmpnx , tmpny , tmpnz ]=np.shape( np.squeeze( ens_bimodality[var] ) )
   [xmax  , ymax  , zmax  ]=np.unravel_index( max_loc , (tmpnx,tmpny,tmpnz) ) 
   #Get the histogram limits for the corresponding maximum bimodality index location.
   hist_max=np.squeeze( hist[var]['maxval'][xmax,ymax,zmax] )
   hist_min=np.squeeze( hist[var]['minval'][xmax,ymax,zmax] )
 
   if ii == 0 :
      hist_max = hist_max * 1000
      hist_min = hist_min * 1000

   hist_delta=(hist_max-hist_min)/nbins
   hist_range=hist_min + hist_delta / 2 + hist_delta *  np.arange(0,nbins,1)
   hist_bars=hist[var]['hist'][xmax,ymax,zmax,:] / np.sum( hist[var]['hist'][xmax,ymax,zmax,:] )
   hist_mean=ens_mean[var][xmax,ymax,zmax]
   hist_var=ens_var[var][xmax,ymax,zmax] 

   if ii == 0 :
      hist_mean = hist_mean * 1000
      hist_var  = hist_var * 1e6
   #hist_skew=ens_skew[var][xmax,ymax,zmax]/np.power(ens_var[var][xmax,ymax,zmax],3/2)
   #hist_kurt=ens_kurt[var][xmax,ymax,zmax]/np.power(ens_var[var][xmax,ymax,zmax],2) -3
   hist_bimodality=ens_bimodality[var][xmax,ymax,zmax]
   #Smooth the histogram
   hist_bars_s=spyf.uniform_filter(hist_bars, size=smooth_range, mode='reflect')
   #Get the Gaussian fit to the histogram
   hist_gaussfit=hist_delta*(1/(np.sqrt(2*np.pi*hist_var)))*np.exp( -0.5*np.power(hist_range-hist_mean,2)/hist_var )
          
   #Make the plot and save it 
   hist_bimodality = 1000 * hist_bimodality
   hist_label= ' Bimodality Index=' + '%1.1f' % hist_bimodality  
   
   the_bars=axs[irow,icol].bar( hist_range , hist_bars_s , width = hist_delta ,color='grey' ) 
   the_lines=axs[irow,icol].plot( hist_range , hist_gaussfit , 'k--',linewidth=4) 
   #axs[irow,icol].legend(handles=[the_bars])
   #axs[irow,icol].set_title(title,fontsize = 15)
   axs[irow,icol].set_xlim(xmin=np.min( hist_range), xmax=np.max( hist_range ) )
   axs[irow,icol].set_ylim(ymin=0.0, ymax=0.1 )
   axs[irow,icol].set_xticklabels( np.round(1.0e2*hist_range[::5])/1.0e2 ,fontsize=12)
   print( hist_bars_s )
   #plt.title('Location x=' + str(xmax) + ' y=' + str(ymax) + ' z=' + str(zmax))
   #plt.show()

   axs[irow,icol].text( np.max(hist_range) - 30*hist_delta , 0.09 ,title ,fontsize=15,color='k',bbox={'facecolor':'white', 'alpha':0.8,'edgecolor':'white'})
   axs[irow,icol].text( np.max(hist_range) - 23*hist_delta , 0.09 ,hist_label ,fontsize=15,color='k',bbox={'facecolor':'white', 'alpha':0.8,'edgecolor':'white'})


#plt.show()
my_fig = './Figure_examples_bimodal_histogram.eps'
plt.savefig( my_fig , format='eps' , dpi=300)
plt.close()
        




