# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016

@author:
"""
import sys
sys.path.append('../../../common_python/common_functions/')
sys.path.append('../../../common_python/common_modules/')

from matplotlib.colors import ListedColormap, LinearSegmentedColormap

import numpy as np
import datetime as dt
import ctl_reader as ctlr
import os
import matplotlib.pyplot as plt
from common_functions import common_functions as cf
import pickle as pkl

import matplotlib
import matplotlib.cm     as cm
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
#from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
from scipy.stats.stats import pearsonr

import common_plot_functions as cpf
import common_mask_functions as cmf

basedir='/home/jruiz/share/LARGE_ENSEMBLE/output_data/'

figname='./Figure_forecast_ref'

exps=['LE_D1_1km_5min_OFP_V2','LE_D1_1km_2min_OFP_V2','LE_D1_1km_1min','LE_D1_1km_30sec_OFP_V2','LE_D1_1km_1min_4D']

lat_radar=34.823
lon_radar=135.523
radar_range=50.0e3   #Radar range in meters (to define the radar mask)

#init_time = '20130713060000'  #Forecast initialization time. 
#lead_time = '20130713063000'  #Lead time in seconds.

for_len=1800
out_freq=30
init_freq=600

clevs=np.arange(0,70,1)

itime = dt.datetime(2013,7,13,5,10,0)  #
etime = dt.datetime(2013,7,13,6,0,0)

titles=['5MIN','2MIN','1MIN','30SEC','1MIN4D']

my_map = cm.get_cmap('gist_ncar',80)
tmp_colors = my_map( range( 80 ) )
for ii in range( 10 ) :
    tmp_colors[ii,:] =np.array([1.0,1.0,1.0,1.0])
    tmp_colors=np.delete( tmp_colors , np.arange( 70 , 80 ) , axis = 0 )
    my_map = ListedColormap( tmp_colors )


#=========================================================================================
#Plot the observed and forecasted REF
#=========================================================================================


#Time loop to read forecast and interpolate them to nearest radar data.
ctime = itime
ifor = 0
while ( ctime <= etime ) :

    fi_time = ctime
    fe_time = ctime + dt.timedelta( seconds = for_len )
    cftime = fi_time

    forecast_time = 0
    while ( cftime <= fe_time ) :

        ncols=3
        nrows=2
        fig, axs = plt.subplots( nrows,ncols , figsize=[10,9.5] , sharex = 'col' , sharey = 'row')
        fig.subplots_adjust(wspace=0.075,hspace=0.28,bottom=0.03,left=0.05,right=0.97,top=0.97)

        icol = 0
        irow = 0
        ax=axs[irow,icol]
        filehandler=open( basedir + '/' + exps[0] +  '/' + dt.datetime.strftime( fi_time ,'%Y%m%d%H%M%S' ) + '/fcstrad/mean/fcstrad_grid' + dt.datetime.strftime( cftime ,'%Y%m%d%H%M%S' ) + '.pkl' , 'rb' )
        radar_grid=pkl.load( filehandler )
        obs_dbzmax=np.nanmax( radar_grid['data_ave'][:,:,:,0] , 0 )
        p=ax.contourf( obs_dbzmax , clevs , cmap = my_map )
        fig.colorbar(p,ax=axs[irow,icol])
        ax.grid()
        ax.set_title('OBS')
        print( np.shape( radar_grid['data_ave'] ) )

        for iexp,my_exp in enumerate(exps)  :
           #Smooth output
           print( basedir + '/' + exps[iexp] +  '/' + dt.datetime.strftime( fi_time ,'%Y%m%d%H%M%S' ) + '/fcstrad/mean/fcstrad_grid' + dt.datetime.strftime( cftime ,'%Y%m%d%H%M%S' ) + '.pkl' )
           filehandler=open( basedir + '/' + exps[iexp] +  '/' + dt.datetime.strftime( fi_time ,'%Y%m%d%H%M%S' ) + '/fcstrad/mean/fcstrad_grid' + dt.datetime.strftime( cftime ,'%Y%m%d%H%M%S' ) + '.pkl' , 'rb' )
           radar_grid=pkl.load( filehandler )

           icol = icol + 1
           if icol > ncols-1   :
              icol = 0
              irow = irow + 1


           ax=axs[irow,icol]
           for_dbzmax= np.nanmax( radar_grid['data_ave'][:,:,:,2] , 0 )
           obs_dbzmax= np.nanmax( radar_grid['data_ave'][:,:,:,0] , 0 )
           p=ax.contourf( for_dbzmax , clevs , cmap=my_map )
           fig.colorbar(p,ax=axs[irow,icol])
           #p2=ax.contour( obs_dbzmax , np.array([20,40]) )
           for_dbzmax=for_dbzmax[40:,60:100]
           obs_dbzmax=obs_dbzmax[40:,60:100]
           rmse = np.round( np.sqrt( np.mean( np.power(  for_dbzmax - obs_dbzmax , 2 ) ) )*100.0) / 100.0
           print( 'RMSE=', rmse )
           tmp_1 = for_dbzmax.reshape( np.size(for_dbzmax) )
           tmp_2 = obs_dbzmax.reshape( np.size(obs_dbzmax) )
           tmp_1=tmp_1.data[np.logical_not( tmp_1.mask ) ]
           tmp_2=tmp_2.data[np.logical_not( tmp_2.mask ) ]
           rcoef =  np.round(pearsonr(tmp_1,tmp_2)[0] * 100.0 )/100.0
           print( 'CORR=', rcoef )
           ax.set_title(titles[iexp] + ' corr=' + str(rcoef) + ' RMSE=' + str(rmse) )
           ax.grid()
 
        #plt.show()
        plt.savefig( 'forecast_ref_i' + dt.datetime.strftime( fi_time ,'%Y%m%d%H%M%S' ) + '_' + dt.datetime.strftime( cftime ,'%Y%m%d%H%M%S' ) + '.png' , format='png' , dpi=300  )
        plt.close()

        forecast_time = forecast_time + 1
        cftime = cftime + dt.timedelta( seconds = out_freq )
    ifor  = ifor + 1
    ctime = ctime + dt.timedelta( seconds = init_freq )


   


