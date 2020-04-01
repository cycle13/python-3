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

basedir='/home/jruiz/share/LARGE_ENSEMBLE/output_data/home/ra001011/a03471/data/output_data/'

expname = '/LE_D1_1km_5min/'


filetypes=['guesgp']   #analgp , analgz , guesgp , guesgz



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

titles=['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)']

nrows=3
ncols=3

fig, axs = plt.subplots( nrows,ncols , figsize=[9,9] )

fig.subplots_adjust(wspace=0.1,hspace=0.1,bottom=0.05,left=0.04,right=0.97,top=0.97)

icol=0
irow=0

varlist=['w','qv','dbz']

for ii in range( 0 , 3 )  :

   
   if ii == 0 :
      iy = 152
      ix = 31
      #ix , iy  = cmf.lat_lon_to_i_j(lon,lat,np.array([136.2]),np.array([34.275])) 
      iz = 8
      file_path = basedir + '/LE_D1_1km_5min/20130713050500/guesgp/'

   if ii == 1 :
      iy = 68
      ix = 114
      #ix , iy  = cmf.lat_lon_to_i_j(lon,lat,np.array([135.28]),np.array([35.05]))
      iz = 8 
      file_path = basedir + '/LE_D1_1km_5min/20130713053000/guesgp/'

   if ii == 2 :
      iy = 68
      ix = 115
      #ix , iy  = cmf.lat_lon_to_i_j(lon,lat,np.array([135.28]),np.array([35.05]))
      iz = 7
      file_path = basedir + '/LE_D1_1km_30sec/20130713053000/guesgp/'


   ix=int(ix) ; iy=int(iy)

   hist_file = file_path + 'histogram.grd'
   max_file  = file_path + 'ensmax.grd'
   min_file  = file_path + 'ensmin.grd'
   mean_file = file_path + 'moment0001.grd'
   var_file = file_path + 'moment0002.grd'
   bim_file = file_path + 'bimodality_index.grd'
   kld_file = file_path + 'kldistance.grd'
   hist=chf.read_histogram(hist_file,max_file,min_file,nx,ny,nbins,ctl_dict,dtypein='i2')
   ens_mean=ctlr.read_data_grads(mean_file,ctl_dict,masked=False)
   ens_var=ctlr.read_data_grads(var_file,ctl_dict,masked=False)
   ens_bim=ctlr.read_data_grads(bim_file,ctl_dict,masked=False)
   ens_kld=ctlr.read_data_grads(kld_file,ctl_dict,masked=False)


   for ivar,my_var in enumerate(varlist)  :

       
      #Get the histogram limits for the corresponding maximum bimodality index location.
      hist_max=np.squeeze( hist[my_var]['maxval'][ix,iy,iz] )
      hist_min=np.squeeze( hist[my_var]['minval'][ix,iy,iz] )
      if ivar == 1 :
         hist_max = 1000 * hist_max
         hist_min = 1000 * hist_min

      hist_delta=(hist_max-hist_min)/nbins
      hist_range=hist_min + hist_delta / 2 + hist_delta *  np.arange(0,nbins,1)
      hist_bars=hist[my_var]['hist'][ix,iy,iz,:] / np.sum( hist[my_var]['hist'][ix,iy,iz,:] )
      #Smooth the histogram
      #hist_bars_s=spyf.uniform_filter(hist_bars, size=smooth_range, mode='reflect')
      hist_bars_s = np.convolve(hist_bars, np.ones((3,))/3, mode='same')
      #hist_bars_s = hist_bars 
      hist_mean = ens_mean[my_var][ix,iy,iz]
      hist_var  = ens_var[my_var][ix,iy,iz]
      hist_kld  = 100*ens_kld[my_var][ix,iy,iz] 
      hist_bim  = 1000*ens_bim[my_var][ix,iy,iz] 
      if ivar == 1 :
          hist_mean = 1000 * hist_mean
          hist_var  = 1e6 * hist_var
      #Get the Gaussian fit to the histogram
      hist_gaussfit=hist_delta*(1/(np.sqrt(2*np.pi*hist_var)))*np.exp( -0.5*np.power(hist_range-hist_mean,2)/hist_var )
      ax=axs[ivar,ii] 
      the_bars=ax.bar( hist_range , hist_bars_s , width = hist_delta ,color='grey' )
      the_lines=ax.plot( hist_range , hist_gaussfit ,'k--',linewidth=1.0) 
      ax.set_xlim(xmin=np.min( hist_range), xmax=np.max( hist_range ) )
      ax.set_ylim(ymin=0.0, ymax=0.2 )
      ax.set_xticklabels( np.round(1.0e2*hist_range[::5])/1.0e2 ,fontsize=12)
      ax.text( np.max(hist_range) - 30*hist_delta , 0.17 , titles[ ii + ivar * ncols ] , fontsize = 15,bbox={'facecolor':'white', 'alpha':0.8,'edgecolor':'white'} )
      texto = 'KLD=' + '%1.1f' % hist_kld + ' B=' + '%1.1f' % hist_bim 
      ax.text( np.max(hist_range) - 25*hist_delta , 0.17 , texto , fontsize = 12 ,bbox={'facecolor':'white', 'alpha':0.8,'edgecolor':'white'})

      ax.set_yticks([0.0,0.1,0.2])
      ax.set_xticks(hist_range[::7])
      if ii != 0  :
         plt.setp(ax.get_yticklabels(), visible=False)
      ax.tick_params(axis='both', which='both', length=0)





#plt.show()
my_fig = './Figure_histograms_CASE_CONVECTION.eps'
plt.savefig( my_fig , format='eps' , dpi=300)
plt.close()
        




