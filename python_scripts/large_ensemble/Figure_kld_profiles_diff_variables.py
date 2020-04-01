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
from scipy.ndimage.filters import gaussian_filter


import common_plot_functions as cpf
import common_mask_functions as cmf

basedir='/home/ra001011/a03471/data/output_data/'

figname='./Figure_kld_profiles_diff_variables.eps'

exps=['LE_D1_1km_5min','LE_D1_1km_30sec']

deltat=[300,120,60,30,30,60]

filetype='guesgp'   #analgp , analgz , guesgp , guesgz

plot_variables=['tk','u','w','dbz']

lat_radar=34.823
lon_radar=135.523
radar_range=50.0e3   #Radar range in meters (to define the radar mask)

sigma_smooth=2.0

#=========================================================
#  LOOP OVER FILE TYPES
#=========================================================

profile_rain = dict()   
profile_norain = dict()

profile_rain_w = dict()

for iexp , my_exp in enumerate(exps)   :

  ctl_file = basedir + '/' + my_exp + '/ctl/' + filetype + '.ctl'

  #=========================================================
  #  READ CTL FILE (To get vertical levels)
  #=========================================================

  if iexp == 0  :
 
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
 
  if iexp == 0 :

     latlon_file = basedir + '/' + my_exp + '/latlon/latlon.grd'

     tmp=ctlr.read_data_records(latlon_file,ctl=ctl_dict,records=np.array([0,1]))
     lat=tmp[:,:,1]
     lon=tmp[:,:,0]

     #Exclude areas outside the radar domain.
     radar_mask = cmf.distance_range_mask( lon_radar , lat_radar , radar_range , lon , lat )

  #=========================================================
  #  READ THE DATA
  #=========================================================

  profile_rain[ my_exp ] = dict()
  profile_norain[ my_exp ] = dict()

  for my_var in plot_variables  :
      my_file=basedir + '/' + my_exp + '/time_mean/'+ filetype + '/kldistance_' + my_var + '_rain_profile.npz' 
      #Remove the fourth element due to post-processing error.
      profile_rain[my_exp][my_var]=np.delete( np.load(my_file)['parameter'] , [0,4] , axis=0 )
      my_file=basedir + '/' + my_exp + '/time_mean/'+ filetype + '/kldistance_' + my_var + '_norain_profile.npz'
      #Remove the fourth element due to post-processing error.
      profile_norain[my_exp][my_var]=np.delete( np.load(my_file)['parameter'] , [0,4] , axis=0 )

  my_file=basedir + '/' + my_exp + '/time_mean/'+ filetype + '/moment0001_w_rain_profile.npz'
  #Remove the fourth element due to post-processing error.
  profile_rain_w[my_exp] = np.delete( np.load(my_file)['parameter'] , [0,4] , axis=0 )

 
print(' Finish the loop over experiments ' )    

#=========================================================================================
#Plot the mean KLD and its standard deviation.
#=========================================================================================

nexp = len( exps )

levels=ctl_dict['vlevels']
levels_str=list()

#To solve the postprocessing error.
levels=np.delete(levels,[0,4],axis=0)
levels[3]=850.0


import matplotlib as mpl
#import matplotlib.cm     as cm
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
#from mpl_toolkits.basemap import Basemap
#import cartopy.crs as ccrs

ncols=4
nrows=2

fig, axs = plt.subplots( nrows,ncols , figsize=[15,5] , )

fig.subplots_adjust(wspace=0.3,hspace=0.25,bottom=0.05,left=0.05,right=0.99,top=0.95)

titles = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)']

#my_map=cpf.cmap_discretize('Blues',20)
my_map=cpf.cmap_discretize('YlGn',20)
clevelsneg = [-100,-40,-20,-10]
clevelspos = [10,20,40,100]

time_index=np.arange( 0 , profile_rain['LE_D1_1km_5min']['tk'].shape[1] ) * 5.0 + 5.0 #Time index in minutes.
xtick= time_index[1::2]

cmapb = mpl.cm.Blues_r(np.linspace(0,1,11))
cmapb = mpl.colors.ListedColormap(cmapb[:-5,:-1])
cmapr = mpl.cm.Reds(np.linspace(0,1,11))
cmapr = mpl.colors.ListedColormap(cmapr[3:,:-1])

blues=['#001a33','#004d99','#0073e6','#3399ff','#66b3ff','#b3d9ff']
reds =['#ffb3b3','#ff6666','#ff0000','#b30000','#660000','#1a0000']

#for key in kld_time_mean  :
for ivar , my_var in enumerate( plot_variables ) :

   irow = 0
   icol = ivar
   isubplot = icol + ncols*irow
   #Plot rainny profiles.
   ax = axs[irow,icol]

   #Shaded
   profile_rain['LE_D1_1km_5min'][my_var][ profile_rain['LE_D1_1km_5min'][my_var] == 0 ] = np.nan
   smin=np.nanmin( profile_rain['LE_D1_1km_5min'][my_var] )
   smax=np.nanmax( profile_rain['LE_D1_1km_5min'][my_var] )
   smax=smax + (smax-smin)/3.0
   p=ax.pcolor( time_index , -np.log( levels ) ,  np.squeeze( profile_rain['LE_D1_1km_5min'][my_var] ) 
           ,vmin=smin , vmax=smax ,cmap=my_map )

   #Grid lines
   ax.grid(linewidth=1.0, color='k',alpha=0.5, linestyle='--')
   #Grid and ticks
   ax.set_xticks(xtick)
   my_levs = levels[[0,3,4,6,8]]
   ytick=-np.log(my_levs)
   ax.set_yticks(ytick)
   #Get the level string list.
   levels_str=[]
   for ilev in my_levs  :
      levels_str.append( str(int(ilev)) )
   ax.set_yticklabels(levels_str,fontsize=10,color='k')
   #Colorbar
   cb=plt.colorbar(p,ax=ax,shrink=0.9)
   cb.ax.tick_params(labelsize=10)
   #Contour map
   percent_diff = 100 * np.squeeze( ( profile_rain[ 'LE_D1_1km_30sec' ][my_var] - profile_rain['LE_D1_1km_5min'][my_var] ) / profile_rain['LE_D1_1km_5min'][my_var] )
   for ii in range( np.size( clevelspos ) )  :
       cpos=ax.contour( time_index , -np.log(levels) , percent_diff , levels = clevelspos[ii] ,colors=reds[3],linestyles='solid')
       plt.clabel(cpos, inline=1, fontsize=10,fmt='%1.0f')
   for ii in range( np.size( clevelsneg ) )  :
       cneg=ax.contour( time_index , -np.log(levels) , percent_diff , levels = clevelsneg[ii] ,colors=blues[1],linestyles='solid')
       plt.clabel(cneg, inline=1, fontsize=10,fmt='%1.0f')
   ax.spines['right'].set_visible(False)
   ax.spines['top'].set_visible(False)
   ax.set_title(titles[isubplot],fontsize = 15)
 
   #Plot no-rainny profiles
   irow=1
   isubplot = icol + ncols*irow
   #Plot rainny profiles.
   ax = axs[irow,icol]

   #Shaded
   profile_norain['LE_D1_1km_5min'][my_var][ profile_norain['LE_D1_1km_5min'][my_var] == 0 ] = np.nan
   smin=np.nanmin( profile_norain['LE_D1_1km_5min'][my_var] )
   smax=np.nanmax( profile_norain['LE_D1_1km_5min'][my_var] )
   smax=smax + (smax-smin)/3.0
   p=ax.pcolor( time_index , -np.log( levels ) ,  np.squeeze( profile_norain['LE_D1_1km_5min'][my_var] )
           ,vmin=smin , vmax=smax ,cmap=my_map )

   #Grid lines
   ax.grid(linewidth=1.0, color='k',alpha=0.5, linestyle='--')
   #Grid and ticks
   ax.set_xticks(xtick)
   my_levs = levels[[0,3,4,6,8]]
   ytick=-np.log(my_levs)
   ax.set_yticks(ytick)
   #Get the level string list.
   levels_str=[]
   for ilev in my_levs  :
      levels_str.append( str(int(ilev)) )
   ax.set_yticklabels(levels_str,fontsize=10,color='k')
   #Colorbar
   cb=plt.colorbar(p,ax=ax,shrink=0.9)
   cb.ax.tick_params(labelsize=10)
   #Contour map
   percent_diff = 100 * np.squeeze( ( profile_norain[ 'LE_D1_1km_30sec' ][my_var] - profile_norain['LE_D1_1km_5min'][my_var] ) / profile_norain['LE_D1_1km_5min'][my_var]  )
   for ii in range( np.size( clevelspos ) )  :
       cpos=ax.contour( time_index , -np.log(levels) , percent_diff , levels = clevelspos[ii] ,colors=reds[3],linestyles='solid')  
       plt.clabel(cpos, inline=1, fontsize=10,fmt='%1.0f')
   for ii in range( np.size( clevelsneg ) )  :
       cneg=ax.contour( time_index , -np.log(levels) , percent_diff , levels = clevelsneg[ii] ,colors=blues[1],linestyles='solid')
       plt.clabel(cneg, inline=1, fontsize=10,fmt='%1.0f')
   ax.spines['right'].set_visible(False)
   ax.spines['top'].set_visible(False)
   ax.set_title(titles[isubplot],fontsize = 15)


plt.show()
plt.savefig( figname , format='eps' , dpi=300)
plt.close()


