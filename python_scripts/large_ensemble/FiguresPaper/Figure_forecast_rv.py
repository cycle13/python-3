# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016

@author:
"""
import sys
sys.path.append('../../../common_python/common_functions/')
sys.path.append('../../../common_python/common_modules/')

import numpy as np
import datetime as dt
import ctl_reader as ctlr
import os
import matplotlib.pyplot as plt
from common_functions import common_functions as cf
import pickle as pkl

import matplotlib
#import matplotlib.cm     as cm
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
#from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
from scipy.stats.stats import pearsonr

import common_plot_functions as cpf
import common_mask_functions as cmf

basedir='/home/jruiz/share/LARGE_ENSEMBLE/output_data/'

figname='./Figure_forecast_rv'

exps=['LE_D1_1km_5min_OFP_V2','LE_D1_1km_2min_OFP_V2','LE_D1_1km_1min','LE_D1_1km_30sec_OFP_V2','LE_D1_1km_1min_4D']

lat_radar=34.823
lon_radar=135.523
radar_range=50.0e3   #Radar range in meters (to define the radar mask)

init_time = '20130713051000'  #Forecast initialization time. 

lead_time = '20130713051000'  #Lead time in seconds.

level = 2

#=========================================================================================
#Plot the observed and forecasted REF
#=========================================================================================

ncols=3
nrows=2
fig, axs = plt.subplots( nrows,ncols , figsize=[10,9.5] , sharex = 'col' , sharey = 'row')
fig.subplots_adjust(wspace=0.075,hspace=0.28,bottom=0.03,left=0.05,right=0.97,top=0.97)

icol = 0
irow = 0

filehandler=open( basedir + '/' + exps[0] +  '/' + init_time + '/fcstrad/mean/fcstrad_grid' + lead_time + '.pkl' , 'rb' )
radar_grid=pkl.load( filehandler )

ax=axs[irow,icol]
obs_rv= radar_grid['data_ave'][level,:,:,1] 
p=ax.pcolor( obs_rv , vmin=0.0 , vmax=60.0 )
print( np.shape( radar_grid['data_ave'] ) )

for iexp,my_exp in enumerate(exps)  :
   #Smooth output
   print( basedir + '/' + exps[iexp] +  '/' + init_time + '/fcstrad/mean/fcstrad_grid' + lead_time + '.pkl' )
   filehandler=open( basedir + '/' + exps[iexp] +  '/' + init_time + '/fcstrad/mean/fcstrad_grid' + lead_time + '.pkl' , 'rb' )
   radar_grid=pkl.load( filehandler )

   icol = icol + 1
   if icol > ncols-1   :
      icol = 0
      irow = irow + 1


   ax=axs[irow,icol]
   for_rv= radar_grid['data_ave'][level,:,:,3] 
   obs_rv= radar_grid['data_ave'][level,:,:,1] 
   p=ax.pcolor( for_rv , vmin=0.0 , vmax=60.0 )
   p2=ax.contour( obs_rv , np.array([20,40]) )
   for_dbzmax=for_rv[40:,60:100]
   obs_dbzmax=obs_rv[40:,60:100]
   print( 'RMSE=', np.sqrt( np.mean( np.power(  for_rv - obs_rv , 2 ) ) ) )
   tmp_1 = for_rv.reshape( np.size(for_rv) )
   tmp_2 = obs_rv.reshape( np.size(obs_rv) )
   tmp_1=tmp_1.data[np.logical_not( tmp_1.mask ) ]
   tmp_2=tmp_2.data[np.logical_not( tmp_2.mask ) ]
   print( 'CORR=', pearsonr(tmp_1,tmp_2) )


 
plt.show()
plt.savefig( figname + '.png' , format='png' , dpi=300  )
plt.close()


