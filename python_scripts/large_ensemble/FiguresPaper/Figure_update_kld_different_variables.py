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

figname='./Figure_update_kld_vertmean_diff_variables'

#exps=['LE_D1_1km_5min','LE_D1_1km_2min','LE_D1_1km_1min','LE_D1_1km_30sec','LE_D1_1km_30sec_nospinup','LE_D1_1km_1min_4D']

deltat=[300,120,60,30,30,60]

filetype='guesgp'   #analgp , analgz , guesgp , guesgz

obs_variables=['dbz','dbz','dbz','v','v' ,'v']
upd_variables=['dbz','tk' ,'w'  ,'v','tk','w']

obs_increment=[5.0  ,5.0  ,5.0  ,2.0,2.0 ,2.0]

nbv=1000

lat_radar=34.823
lon_radar=135.523
radar_range=50.0e3   #Radar range in meters (to define the radar mask)

init_times = ['20130713050500','20130713050400','20130713050500','20130713050500','20130713050500','20130713050500'] 

sigma_smooth=2.0

#Define initial and end times using datetime module.
#itime = dt.datetime(2013,7,13,5, 5,0)  #Initial time.
etime = dt.datetime(2013,7,13,5,55,0)  #End time.

#=========================================================
#  LOOP OVER FILE TYPES
#=========================================================

parameter = dict()   

ensemble_mean = dict()

my_exp='LE_D1_1km_30sec'

ctl_file = basedir + '/' + my_exp + '/ctl/update_mean_diff.ctl'

ctl_file_2 = basedir + '/' + my_exp + '/ctl/moment0001_for.ctl'

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

ctl_dict_2 = ctlr.read_ctl( ctl_file_2 ) 

#=========================================================
#  READ LAT LON
#=========================================================

latlon_file = basedir + '/' + my_exp + '/latlon/latlon.grd'

tmp=ctlr.read_data_records(latlon_file,ctl=ctl_dict,records=np.array([0,1]))
lat=tmp[:,:,1]
lon=tmp[:,:,0]

#Exclude areas outside the radar domain.
radar_mask = cmf.distance_range_mask( lon_radar , lat_radar , radar_range , lon , lat )

#=========================================================
#  READ THE DATA
#=========================================================

parameter_5MIN = []
parameter_30SEC = []

for iv,my_var in enumerate( obs_variables ) :

   var_obs = obs_variables[iv]
   var_upd = upd_variables[iv]

   my_file=basedir + '/LE_D1_1km_5min/time_mean/guesgp/update_comp_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '_kld' + '.grd'
   parameter_5MIN.append( np.nanmean( np.delete( ctlr.read_data(my_file,ctl_dict,undef2nan=False) , 4 , 2 ) , 2 ) )
   my_file=basedir + '/LE_D1_1km_30sec/time_mean/guesgp/update_comp_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '_kld' + '.grd'
   parameter_30SEC.append( np.nanmean( np.delete( ctlr.read_data(my_file,ctl_dict,undef2nan=False) , 4 , 2 ) , 2 ) )

   parameter_5MIN[iv]=cf.gaussian_filter(field0=parameter_5MIN[iv],dx=1.0,sigma=sigma_smooth,nx=nx,ny=ny,undef=undef)
   parameter_30SEC[iv]=cf.gaussian_filter(field0=parameter_30SEC[iv],dx=1.0,sigma=sigma_smooth,nx=nx,ny=ny,undef=undef)

   parameter_5MIN[iv][ parameter_5MIN[iv] == undef ] = np.nan
   parameter_30SEC[iv][ parameter_30SEC[iv] == undef ] = np.nan 


   print( my_var )
   print( np.nanmin( parameter_5MIN[iv] ) , np.nanmax( parameter_5MIN[iv] ) )
   print( np.nanmin( parameter_30SEC[iv] ) , np.nanmax( parameter_30SEC[iv] ) )

my_file=basedir + '/LE_D1_1km_30sec/time_mean/guesgp/moment0001_mean.grd'

max_dbz= np.squeeze( np.nanmax( np.delete( ctlr.read_data_grads(my_file,ctl_dict_2,masked=False,undef2nan=True)['dbz'] , 4 , 2 ) , 2 ) )
max_dbz[ np.logical_not( radar_mask )  ] = np.nan

print('Finish reading data')  

#=========================================================================================
#Plot the mean KLD and its standard deviation.
#=========================================================================================

#Start subplots
ncols = 3
nrows = 2
icoldelta = 1.0/ncols
irowdelta = 1.0/nrows
hmargin=0.0
vmargin=-0.5
hoffset=0.0
voffset=0.0

icol = 0
irow = 0

xtick=[134.5,135,135.5,136,136.5,137]
ytick=[34,34.5,35,35.5]
axesrange=[134.97,136.09,34.36,35.30]
titles = ['(a) Z - Z','(b) Z - T','(c) Z - W','(d) V - V','(e) V - T','(f) V - W']

#Plot time mean of the moments.
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs


ifig = 0 

fig, axs = plt.subplots( nrows , ncols , subplot_kw=dict(projection=ccrs.PlateCarree() ) , figsize=[14,7] , sharex = 'col' , sharey = 'row' )

fig.subplots_adjust(wspace=0.1,hspace=0.0,bottom=0.05,left=0.04,right=0.92,top=0.97)

blues=['#001a33','#004d99','#0073e6','#3399ff','#66b3ff','#b3d9ff']
reds =['#ffb3b3','#ff6666','#ff0000','#b30000','#660000','#1a0000']


for ivar,my_var in enumerate(obs_variables) :

   ifig = icol + ncols*(irow) 

   plot_kld_5min   = np.squeeze( parameter_5MIN[ivar] )
   plot_kld_30sec  = np.squeeze( parameter_30SEC[ivar] )

   plot_kld_diff = 100 * ( ( plot_kld_30sec - plot_kld_5min ) / plot_kld_5min )

   plot_kld_5min[ np.logical_not( radar_mask )  ] = np.nan
   plot_kld_diff[ np.logical_not( radar_mask )  ] = np.nan


   #Axes limits
   #my_axes = [icoldelta*(icol-1)+hmargin+hoffset,irowdelta*(irow-1)+vmargin+voffset,irowdelta-2*hmargin,icoldelta-2*vmargin]
   ax = axs[irow,icol]  #plt.axes( my_axes , facecolor=None , projection=ccrs.PlateCarree() )

   #The pcolor
   ncolors = 22
   smin=-110
   smax=110
   delta = (smax-smin)/ncolors
   clevels = np.arange( smin , smax + delta , delta )
   #p=ax.pcolor(lon , lat ,  np.transpose( np.squeeze( plot_bim_5min ) ) ,
   #     transform=ccrs.PlateCarree(),vmin=smin , vmax=smax ,cmap=cpf.cmap_discretize('YlGn',21) )
   my_map = cpf.cmap_discretize('coolwarm',ncolors)
   p=ax.contourf(lon , lat ,  np.transpose( np.squeeze( plot_kld_diff ) ) ,
        transform=ccrs.PlateCarree(),levels=clevels ,cmap=my_map )

   ax.set_extent( axesrange , ccrs.PlateCarree())
   gl=ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                linewidth=1.0, color='k',alpha=0.5, linestyle='--')
   #Grid and ticks
   gl.xlocator=mticker.FixedLocator(xtick)
   gl.ylocator=mticker.FixedLocator(ytick)
   #Coastline plot
   ax.coastlines('10m',linewidth=1.0)
   #Contour map
   if ivar == 0  :
      clevs = [50,100,150,200]
   if ivar == 1  :
      clevs = [3,6,9,12]
   if ivar == 2  :
      clevs = [3,6,9,12,15]
   if ivar == 3  :
      clevs = [1,2,3,4]
   if ivar == 4  :
      clevs = [1,2,3,4]
   if ivar == 5  :
      clevs = [0.5,1.0,1.5,2.0,2.5]

   cpos=ax.contour( lon , lat , 100*np.transpose( plot_kld_5min ) ,  colors='k',linestyles='solid')
   plt.clabel(cpos, inline=1, fontsize=10,fmt='%1.1f')

   #Axes label and title

   if irow == 1   :
      plt.xlabel('Latitude')
   if icol == 1   :
      plt.ylabel('Longitude')

   #ax.set_title(titles[ifig],fontsize = 15)
   ax.text(135.01,35.19,titles[ifig],fontsize=15,color='k',bbox={'facecolor':'white', 'alpha':0.8,'edgecolor':'white'})
   if icol == 0  :
      gl.ylabel_style = {'size': 12, 'color': 'k' }
   else          :
      gl.ylabel_style = {'visible': False }
   if irow == 1  :
      gl.xlabel_style = {'size': 12, 'color': 'k' }
   else          :
      gl.xlabel_style = {'visible': False }


   gl.xlabels_top = False
   gl.ylabels_right = False

   irow = irow + 1
   if irow >= nrows   :
      irow = 0
      icol = icol + 1

cbar_ax = fig.add_axes([0.93, 0.05, 0.02, 0.9])
m = plt.cm.ScalarMappable(cmap=my_map)
m.set_array(np.transpose(plot_kld_diff))
m.set_clim(smin,smax)
cb=plt.colorbar(m,cax=cbar_ax,boundaries=np.arange(smin,smax+delta,delta))


 
#plt.show()
plt.savefig( figname + '.eps' , format='eps' , dpi=300  )
plt.close()


