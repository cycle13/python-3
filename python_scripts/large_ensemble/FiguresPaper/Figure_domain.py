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


import common_plot_functions as cpf
import common_mask_functions as cmf

basedir='/home/ra001011/a03471/data/output_data/'

exps=['LE_D1_1km_5min']

lat_radar=34.823
lon_radar=135.523
radar_range=60.0e3   #Radar range in meters (to define the radar mask)

#=========================================================
#  READ LAT LON
#=========================================================
 
latlon_file = basedir + '/LE_D1_1km_5min/latlon/latlon.grd'
latlon_ctl  = basedir + '/LE_D1_1km_5min/latlon/latlon.ctl'

ctl_dict = ctlr.read_ctl( latlon_ctl )

my_data=ctlr.read_data_grads(latlon_file,ctl=ctl_dict)

lon = np.squeeze( my_data['glon'] )
lat = np.squeeze( my_data['glat'] )
topo= np.squeeze( my_data['topo'] )

#Exclude areas outside the radar domain.
radar_mask = cmf.distance_range_mask( lon_radar , lat_radar , radar_range , lon , lat )

#=========================================================================================
#Plot the mean topography and PAWR radar range
#=========================================================================================

xtick=[134.0,134.5,135,135.5,136,136.5,137,137.5]
ytick=[33.0,33.5,34,34.5,35,35.5,36.0,36.5]
axesrange=[134.4,136.7,33.75,35.8]

#Plot time mean of the moments.

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs


blues=['#001a33','#004d99','#0073e6','#3399ff','#66b3ff','#b3d9ff']
reds =['#ffb3b3','#ff6666','#ff0000','#b30000','#660000','#1a0000']


fig=plt.figure(1,figsize=[6.5,5])

ax = fig.add_subplot(111, projection=ccrs.PlateCarree())

#The pcolor
smin=np.nanmin( my_data['topo'] )
smax=np.nanmax( my_data['topo'] ) 

p=ax.pcolor( lon , lat ,  topo ,
     transform=ccrs.PlateCarree(),vmin=smin , vmax=smax ,cmap=cpf.cmap_discretize('copper_r',41) )
ax.set_extent( axesrange , ccrs.PlateCarree())

topo[ topo > 1.2 ] = np.nan 
#Colorbar
cb=plt.colorbar(p, ax=ax, orientation='vertical' , shrink=0.9 )
cb.ax.tick_params(labelsize=10)



p=ax.pcolor( lon , lat ,  topo ,
     transform=ccrs.PlateCarree(),vmin=-1.0 , vmax=10.0 ,cmap=cpf.cmap_discretize('Blues_r',10) )
ax.set_extent( axesrange , ccrs.PlateCarree())



gl=ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
             linewidth=1.0, color='k',alpha=0.5, linestyle='--')
#Grid and ticks
gl.xlocator=mticker.FixedLocator(xtick)
gl.ylocator=mticker.FixedLocator(ytick)
#Coastline plot
ax.coastlines('10m',linewidth=1.0)

ax.contour( lon , lat , np.logical_not( radar_mask ).astype(int) , levels = [0.9] , colors=reds[3] , linestyles='solid' , linewidths=3 )
clevelsneg = [-100,-40,-20]
clevelspos = [20,40,100]
plt.xlabel('Latitude')
plt.ylabel('Longitude')

gl.xlabel_style = {'size': 12, 'color': 'k' }
gl.ylabel_style = {'size': 12, 'color': 'k' }
gl.xlabels_top = False
gl.ylabels_right = False


#plt.show()
plt.savefig( 'Figure_domain.eps' , format='eps' , dpi=300  )
plt.close()


