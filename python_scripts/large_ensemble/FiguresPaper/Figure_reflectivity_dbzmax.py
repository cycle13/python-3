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
from common_functions import common_functions as cf

#For radar
import radartools as rt 
from common_qc_tools  import qc  #Fortran code routines.

import common_mask_functions as cmf
import common_plot_functions as cpf



basedir='/home/jruiz/share/LARGE_ENSEMBLE/output_data/home/ra001011/a03471/data/output_data/'

figname='./Figure_reflectivity_dbzmax.eps'

exps=['LE_D1_1km_5min','LE_D1_1km_2min','LE_D1_1km_1min','LE_D1_1km_30sec']

filetype='guesgp'   #analgp , analgz , guesgp , guesgz

lat_radar=34.823
lon_radar=135.523
radar_range=50.0e3   #Radar range in meters (to define the radar mask)

scale_date = '20130713055000'
radar_date = '20130713-145010'


xtick=[134.5,135,135.5,136,136.5,137]
ytick=[34,34.25,34.5,34.75,35,35.25,35.5]
axesrange=[134.97,136.09,34.65,35.30]

#=========================================================
#  READ'N PLOT
#=========================================================

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib import patches
import cartopy.crs as ccrs

titles=['(a)','(b)','(c)','(d)']

ncols=2
nrows=2
fig, axs = plt.subplots( nrows,ncols , subplot_kw=dict(projection=ccrs.PlateCarree() ), figsize=[10,6.0] , sharex = 'col' , sharey = 'row')

fig.subplots_adjust(wspace=0.01,hspace=-0.2,bottom=0.02,left=0.07,right=0.88,top=0.95)


#=========================================================
#  READ RADAR DATA
#=========================================================


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

#ref[ref==undef]=np.nan
max_dbz_obs=var2d[:,:,6]
max_dbz_obs[max_dbz_obs < 0.0]=0.0


for iexp , my_exp in enumerate(exps)   :

  ctl_file = basedir + '/' + my_exp + '/ctl/moment0001_for.ctl'

  #=========================================================
  #  READ CTL FILE
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

  my_file=basedir + '/' + my_exp + '/20130713055000/guesgp/moment0001.grd' 

  ens_mean=ctlr.read_data_grads(my_file,ctl_dict,masked=False)

  dbz=np.delete( ens_mean['dbz'] , 4 , 2 )
  dbz_max = np.nanmax( dbz , 2 )
  
  dbz_max[np.logical_not(radar_mask)]=np.nan
  dbz_max[dbz_max < 0]=np.nan

  if iexp == 0 :
     irow=0;icol=0
  if iexp == 1 :
     irow=0;icol=1
  if iexp == 2 :
     irow=1;icol=0
  if iexp == 3 :
     irow=1;icol=1
      
  ax=axs[irow,icol]
  smin=0
  smax=70
  ncols=70
  delta = ( smax - smin )/ncols
  clevs=np.arange(smin,smax+delta,delta)
  my_map = cpf.cmap_discretize('gist_ncar',100)
  p=ax.contourf(lon , lat ,  np.squeeze( dbz_max ) , clevs , 
               transform=ccrs.PlateCarree() ,cmap=my_map)

  c=ax.contour(lon_rad , lat_rad ,  max_dbz_obs ,colors='k',
               transform=ccrs.PlateCarree(),levels=[5,40],linewidths=2,linestyles='-' )
  #c=ax.pcolor( lon_rad , lat_rad ,  max_dbz_obs )
  
  #plt.clabel(c, inline=1, fontsize=10,fmt='%1.0f')
  
  ax.set_extent( axesrange , ccrs.PlateCarree())
  ax.coastlines('10m',linewidth=1.0)
   
  if irow == 1   :
      plt.xlabel('Latitude')
  if icol == 1   :
      plt.ylabel('Longitude')
      
  gl=ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                linewidth=1.0, color='k',alpha=0.5, linestyle='--')
  #Grid and ticks
  gl.xlocator=mticker.FixedLocator(xtick)
  gl.ylocator=mticker.FixedLocator(ytick)    
      
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
  
  ax.text(135.0,35.22,titles[iexp],fontsize=18,color='k',bbox={'facecolor':'white', 'alpha':0.8,'edgecolor':'white'})


cbar_ax = fig.add_axes([0.90, 0.09, 0.03, 0.8])
m = plt.cm.ScalarMappable(cmap=my_map)
m.set_array(np.transpose(dbz_max))
m.set_clim(smin,smax)
cb=plt.colorbar(m,cax=cbar_ax,boundaries=np.arange(smin,smax+delta,delta))

#plt.show()
plt.savefig( figname , format='eps' , dpi=300)
plt.close()

