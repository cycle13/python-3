# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016

@author:
"""
import sys
sys.path.append('../../../common_python/common_functions/')
sys.path.append('../../../common_python/common_modules/')
sys.path.append('../../../python_scripts/gyl_scripts/')
sys.path.append('../../../python_scripts/radar_qc/src/fortran/')

import numpy as np
import datetime as dt
import ctl_reader as ctlr
import os
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from common_functions import common_functions as cf

#For radar
import radartools as rt 
from common_qc_tools  import qc  #Fortran code routines.

import common_mask_functions as cmf
import common_plot_functions as cpf



basedir='/home/jruiz/share/LARGE_ENSEMBLE/output_data/home/ra001011/a03471/data/output_data/'

figname='./Figure_reflectivity_dbzmax_domain'


filetype='guesgp'   #analgp , analgz , guesgp , guesgz

lat_radar=34.823
lon_radar=135.523
radar_range=60.0e3   #Radar range in meters (to define the radar mask)

scale_date = '20130713055000'
radar_date = '20130713-145010'


xtick=np.arange(134.0,138.0,0.5)
ytick=np.arange(34.0,36.0,0.5)

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib import patches
import cartopy.crs as ccrs

titles=['(b)','(c)']

ncols=2
nrows=2
fig, axs = plt.subplots( nrows,ncols , subplot_kw=dict(projection=ccrs.PlateCarree() ), figsize=[10,8.5] )

fig.subplots_adjust(wspace=0.13,hspace=0.05,bottom=-0.02,left=0.07,right=0.98,top=0.99)

#=========================================================
#  PLOT MODEL DOMAIN
#=========================================================

lat_radar=34.823
lon_radar=135.523
radar_range=60.0e3   #Radar range in meters (to define the radar mask)

#=========================================================
#  READ LAT LON AND TOPO
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
axesrange=[134.49,136.55,33.99,35.65]

blues=['#001a33','#004d99','#0073e6','#3399ff','#66b3ff','#b3d9ff']
reds =['#ffb3b3','#ff6666','#ff0000','#b30000','#660000','#1a0000']

#fig=plt.figure(1,figsize=[6.5,5])
#ax=plt.subplot(121,projection=ccrs.PlateCarree())
ax1 = axs[0,0]   #plt.axes([0.05, 0.075, 0.44, 1.0],projection=ccrs.PlateCarree())
#ax = fig.add_subplot(111, projection=ccrs.PlateCarree())

#The pcolor
smin=np.nanmin( my_data['topo'] )
smax=np.nanmax( my_data['topo'] ) 

p=ax1.pcolor( lon , lat ,  topo ,
             transform=ccrs.PlateCarree(),vmin=smin , vmax=smax ,cmap=cpf.cmap_discretize('copper_r',41) )

topo[ topo > 1.2 ] = np.nan 
#Colorbar
p=ax1.pcolor( lon , lat ,  topo ,
             transform=ccrs.PlateCarree(),vmin=-1.0 , vmax=10.0 ,cmap=cpf.cmap_discretize('Blues_r',10) )
gl=ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1.0, color='k',alpha=0.5, linestyle='--')

ax1.plot( lon_radar , lat_radar , "rX" , transform=ccrs.PlateCarree() , markersize=15 )
#Grid and ticks
gl.xlocator=mticker.FixedLocator(xtick)
gl.ylocator=mticker.FixedLocator(ytick)

#Coastline plot
ax1.coastlines('10m',linewidth=1.0)
ax1.contour( lon , lat , np.logical_not( radar_mask ).astype(int) , levels = [0.9] , colors=reds[3] , linestyles='solid' , linewidths=3 )
#ax1.text(134.53,35.496,'(a)',fontsize=20,color='k',bbox={'facecolor':'white', 'alpha':0.8,'edgecolor':'white'})

#Plot the box for the reflectivity panels
ax1.plot(np.array([134.97,134.97]),np.array([34.5,35.3]), "k-",linewidth=3, transform=ccrs.PlateCarree())
ax1.plot(np.array([136.09,136.09]),np.array([34.5,35.3]), "k-",linewidth=3, transform=ccrs.PlateCarree())
ax1.plot(np.array([134.97,136.09]),np.array([35.3,35.3]), "k-",linewidth=3, transform=ccrs.PlateCarree())
ax1.plot(np.array([134.97,136.09]),np.array([34.5,34.5]), "k-",linewidth=3, transform=ccrs.PlateCarree())
#axesrange=[134.97,136.09,34.65,35.30]

#ax.plot(np.array([135.28,135.28]),np.array([34.9,35.15]), "b-",linewidth=3, transform=ccrs.PlateCarree())
#ax.plot(np.array([136.2,136.2]),np.array([34.2,34.4]), "b-",linewidth=3, transform=ccrs.PlateCarree())

clevelsneg = [-100,-40,-20]
clevelspos = [20,40,100]
ax1.set_xlabel('Latitude')
ax1.set_ylabel('Longitude')

gl.xlabel_style = {'size': 12, 'color': 'k' }
gl.ylabel_style = {'size': 12, 'color': 'k' }
gl.xlabels_top = False
gl.ylabels_right = False
ax1.set_xlim( axesrange[0],axesrange[1] , ccrs.PlateCarree())
ax1.set_ylim( axesrange[2],axesrange[3] , ccrs.PlateCarree())
ax1.set_title('(a)',fontsize=14)

#COLORBAR FOR TOPOGRAPHY
cbar_ax = fig.add_axes([0.08, 0.49, 0.40, 0.02])
m = plt.cm.ScalarMappable(cmap=cpf.cmap_discretize('copper_r',41) )
m.set_array(np.transpose(topo))
m.set_clim(smin,smax)
cb=plt.colorbar(m,cax=cbar_ax,boundaries=np.arange(0,2000,250),orientation='horizontal')

xtick=np.arange(134.0,138.0,0.25)
ytick=np.arange(34.0,36.0,0.25)

axesrange2=[134.97,136.09,34.5,35.40]
#=========================================================
#  PLOT RADAR OBSERVATIONS
#=========================================================
ax=axs[0,1]

radar_file='/home/jruiz/share/LARGE_ENSEMBLE/radar_data/PAWR_LETKF_INPUTV4_' + radar_date + '.dat'

if  os.path.isfile('./radar_data'+radar_date+'.npz') :

   radar_data=np.load('./radar_data'+radar_date+'.npz')
   lon_rad=radar_data['lon_rad']
   lat_rad=radar_data['lat_rad']
   ref=radar_data['ref']
   na=radar_data['na']
   nr=radar_data['nr']
   ne=radar_data['ne']
   hgt=radar_data['hgt']
   rrange=radar_data['rrange']
   undef=np.min(ref)
   #ref[ref==undef]=np.nan
    
else                                   :
    
   radar_data=rt.radarobs_read(radar_file, endian='')
   rt.radar_georeference(radar_data)
   lon_rad=radar_data['lon']
   lat_rad=radar_data['lat']
   ref=radar_data['ref']
   na=radar_data['na']
   nr=radar_data['nr']
   ne=radar_data['ne']
   hgt=radar_data['symhgt']
   rrange=radar_data['radi_h']
   undef=np.min(ref)
   #ref[ref==undef]=np.nan
   
   np.savez('./radar_data'+radar_date+'.npz', lon_rad=lon_rad,lat_rad=lat_rad,ref=ref,na=na,nr=nr,ne=ne,hgt=hgt,rrange=rrange)
   
ref=np.transpose(ref, (2, 1, 0))
lon_rad=np.transpose( lon_rad,(2,1,0))[:,:,0]
lat_rad=np.transpose( lat_rad,(2,1,0))[:,:,0]
hgt=np.transpose(hgt)
rrange=np.transpose(rrange)

[var3d,var2d]=qc.echo_top(reflectivity=ref,heigth=hgt,rrange=rrange,na=na,nr=nr,ne=ne,undef=undef,nx=1,ny=1,nz=1)

max_dbz_obs=var2d[:,:,6]
#max_dbz_obs[max_dbz_obs < 0.1]=np.nan

smin=-5
smax=70
ncols=75
delta = ( smax - smin )/ncols
clevs=np.arange(smin,smax+delta,delta)

my_map = cm.get_cmap('gist_ncar',ncols+10)
tmp_colors = my_map( range( ncols +10 ) )
for ii in range( 5 ) :
   tmp_colors[ii,:] =np.array([1.0,1.0,1.0,1.0])
tmp_colors=np.delete( tmp_colors , np.arange( 80 , 85 ) , axis = 0 )
#for ii in range( 80 , 85 ) :
#   tmp_colors[ii,:] = tmp_colors[80,:]
my_map = ListedColormap( tmp_colors )

#my_map = cpf.cmap_discretize('gist_ncar',100)

p=ax.contourf(lon_rad , lat_rad ,  max_dbz_obs  , clevs ,
             transform=ccrs.PlateCarree() ,cmap=my_map)

ax.coastlines('10m',linewidth=1.0)

gl=ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
              linewidth=1.0, color='k',alpha=0.5, linestyle='--')
#Grid and ticks
gl.xlocator=mticker.FixedLocator(xtick)
gl.ylocator=mticker.FixedLocator(ytick)

gl.ylabel_style = {'size': 12, 'color': 'k' }
gl.xlabel_style = {'size': 12, 'color': 'k' }

gl.xlabels_top = False
gl.ylabels_right = False

ax.set_title('(b)',fontsize=14)
#ax.text(135.0,35.217,titles[iexp],fontsize=18,color='k',bbox={'facecolor':'white', 'alpha':0.8,'edgecolor':'white'})

ax.plot(np.array([135.28,135.28]),np.array([34.9,35.15]), "r-",linewidth=3, transform=ccrs.PlateCarree())
ax.set_xlim( axesrange2[0],axesrange2[1] , ccrs.PlateCarree())
ax.set_ylim( axesrange2[2],axesrange2[3] , ccrs.PlateCarree())

#=========================================================
#  PLOT 5 MIN EXPERIMENT
#=========================================================
ax=axs[1,0]

ctl_file = basedir + '/LE_D1_1km_5min/ctl/moment0001_for.ctl'
ctl_dict = ctlr.read_ctl( ctl_file )
nx=ctl_dict['nx']
ny=ctl_dict['nx']
nlev=ctl_dict['nz']
nt=int(1)             #Force the number of times to be one.
ctl_dict['nt']=int(1) #Force the number of times to be one.
undef=np.float32( ctl_dict['undef'] )

#  READ LAT LON
latlon_file = basedir + '/LE_D1_1km_5min/latlon/latlon.grd'

tmp=ctlr.read_data_records(latlon_file,ctl=ctl_dict,records=np.array([0,1]))
lat=tmp[:,:,1]
lon=tmp[:,:,0]

#Exclude areas outside the radar domain.
radar_mask = cmf.distance_range_mask( lon_radar , lat_radar , radar_range , lon , lat )

#  READ THE DATA
my_file=basedir + '/LE_D1_1km_5min/'+scale_date+'/analgp/moment0001.grd' 
ens_mean=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
dbz=np.delete( ens_mean['dbz'] , 4 , 2 )
dbz_max = np.nanmax( dbz , 2 )
dbz_max[np.logical_not(radar_mask)]=np.nan
#dbz_max[dbz_max < 0]=np.nan

p=ax.contourf(lon , lat ,  np.squeeze( dbz_max ) , clevs , 
               transform=ccrs.PlateCarree() ,cmap=my_map)

#c=ax.contour(lon_rad , lat_rad ,  max_dbz_obs ,colors='k',
#               transform=ccrs.PlateCarree(),levels=[5,40],linewidths=2,linestyles='-' )
ax.coastlines('10m',linewidth=1.0)
ax.set_title('(c)',fontsize=14)
ax.plot(np.array([135.28,135.28]),np.array([34.9,35.15]), "r-",linewidth=3, transform=ccrs.PlateCarree())
ax.set_xlim( axesrange2[0],axesrange2[1] , ccrs.PlateCarree())
ax.set_ylim( axesrange2[2],axesrange2[3] , ccrs.PlateCarree()) 

 
gl=ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                linewidth=1.0, color='k',alpha=0.5, linestyle='--')
#Grid and ticks
gl.xlocator=mticker.FixedLocator(xtick)
gl.ylocator=mticker.FixedLocator(ytick)    
      
gl.ylabel_style = {'size': 12, 'color': 'k' }
gl.xlabels_top = False
gl.ylabels_right = False

#=========================================================
#  PLOT 30SEC EXPERIMENT
#=========================================================
ax=axs[1,1]

#  READ THE DATA
my_file=basedir + '/LE_D1_1km_30sec/'+scale_date+'/analgp/moment0001.grd' 
ens_mean=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
dbz=np.delete( ens_mean['dbz'] , 4 , 2 )
dbz_max = np.nanmax( dbz , 2 )
dbz_max[np.logical_not(radar_mask)]=np.nan
#dbz_max[dbz_max < 0]=np.nan

p=ax.contourf(lon , lat ,  np.squeeze( dbz_max ) , clevs , 
               transform=ccrs.PlateCarree() ,cmap=my_map)

#c=ax.contour(lon_rad , lat_rad ,  max_dbz_obs ,colors='k',
#               transform=ccrs.PlateCarree(),levels=[5,40],linewidths=2,linestyles='-' )
ax.coastlines('10m',linewidth=1.0)
ax.set_title('(d)',fontsize=14)
ax.plot(np.array([135.28,135.28]),np.array([34.9,35.15]), "r-",linewidth=3, transform=ccrs.PlateCarree())
ax.set_xlim( axesrange2[0],axesrange2[1] , ccrs.PlateCarree())
ax.set_ylim( axesrange2[2],axesrange2[3] , ccrs.PlateCarree()) 

 
gl=ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                linewidth=1.0, color='k',alpha=0.5, linestyle='--')
#Grid and ticks
gl.xlocator=mticker.FixedLocator(xtick)
gl.ylocator=mticker.FixedLocator(ytick)    
      
gl.ylabel_style = {'size': 12, 'color': 'k' }
gl.xlabels_top = False
gl.ylabels_right = False 



#COLORBAR FOR REFLECTIVITY
cbar_ax = fig.add_axes([0.565, 0.49, 0.40, 0.02])
m = plt.cm.ScalarMappable(cmap=my_map)
m.set_array(np.transpose(dbz_max))
m.set_clim(smin,smax)
cb=plt.colorbar(m,cax=cbar_ax,orientation='horizontal',boundaries=np.arange(smin,smax+delta,delta))



#plt.show()
#plt.savefig( figname + '.eps', format='eps' , dpi=300)
plt.savefig( figname + '.png' , format='png' , dpi=300)
plt.close()

