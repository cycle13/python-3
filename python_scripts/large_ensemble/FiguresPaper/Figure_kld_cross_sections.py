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

import common_plot_functions as cpf
import common_mask_functions as cmf


basedir='/home/jruiz/share/LARGE_ENSEMBLE/output_data/home/ra001011/a03471/data/output_data/'

#cases=['CASE-A','CASE-B','CASE-C']

expnames = ['/LE_D1_1km_5min/','/LE_D1_1km_30sec/']

cen_lon = [135.28]
#cen_lon = [135.39]
#cen_lon = [135.732]
ini_lat = [34.9]
#ini_lat = [34.8]
end_lat = [35.15]

buffer_size=0.007  #Delta longtitude to perform zonal average.

filetypes=['guesgp']   #analgp , analgz , guesgp , guesgz

#Define initial and end times using datetime module.
ctimes = [dt.datetime(2013,7,13,5,30,0)]  #Initial time.
#ctimes = [dt.datetime(2013,7,13,5,50,0)]
#ctimes = [dt.datetime(2013,7,13,5,55,0)]

#Define the delta.
delta=dt.timedelta(seconds=60)

ctl_file = basedir + expnames[0] + '/ctl/moment0001_for.ctl'

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

latlon_file = basedir + expnames[0] + '/latlon/latlon.grd'

tmp=ctlr.read_data_records(latlon_file,ctl=ctl_dict,records=np.array([0,1]))
lat=np.squeeze(tmp[:,:,1])
lon=np.squeeze(tmp[:,:,0])

#=========================================================
#  READ'N PLOT
#=========================================================

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib import patches

ncolors=10
my_map=cpf.cmap_discretize('YlGnBu',ncolors)

tmp_lat = lat[:,90]
levels=ctl_dict['vlevels']

tick_levels=[1000,850,700,500,300]
levels_str=list()
levels=np.delete(levels,4,axis=0)
levels[3]=850.0

titles=['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)']

#Get the level string list.
levels_str=[]
for ilev in tick_levels  :
   levels_str.append( str(int(ilev)) )


ncols=3
nrows=2
fig, axs = plt.subplots( nrows,ncols , figsize=[15,10] , )
fig.subplots_adjust(wspace=0.05,hspace=0.1,bottom=0.1,left=0.05,right=0.97,top=0.95)

variables=['w','tk','dbz']

ititle = 0

for iexp,my_exp in enumerate( expnames ) :
 
   if iexp == 0 :
       my_title_2 = '5MIN'

   if iexp == 1 :
       my_title_2 = '30SEC'

   irow = iexp 

   mask = cmf.box_mask( cen_lon[0] - buffer_size , cen_lon[0] + buffer_size , ini_lat[0] , end_lat[0] , lon , lat ) 

   #=========================================================
   #  READ
   #=========================================================

   my_file=basedir + expnames[iexp] + ctimes[0].strftime("%Y%m%d%H%M%S") + '/'+filetypes[0]+'/moment0001.grd'
   m01=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=True)
   my_file=basedir + expnames[iexp] + ctimes[0].strftime("%Y%m%d%H%M%S") + '/'+filetypes[0]+'/moment0002.grd'
   m02=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=True)
   my_file=basedir + expnames[iexp] + ctimes[0].strftime("%Y%m%d%H%M%S") + '/'+filetypes[0]+'/kldistance.grd'
   kld=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=True)

   #Compute the variables that we want to plot over the cross section.
   for ivar,my_var in enumerate( variables ) :

      icol = ivar

      var_mean=np.squeeze(np.delete( m01[my_var],4,2))
      var_kld =np.squeeze(np.delete(kld[my_var],4,2))
      var_var=np.sqrt( np.squeeze(np.delete(m02[my_var],4,2)) )
      dbz=np.squeeze(np.delete(m01['dbz'],4,2))
      v=np.squeeze(np.delete(m01['v'],4,2))
      w=np.squeeze(np.delete(m01['w'],4,2))
      tk=np.squeeze(np.delete(m01['tk'],4,2))

      for ii in range( np.shape( var_mean )[2] ) :

         var_mean[:,:,ii][np.logical_not( mask )]=np.nan  
         var_kld[:,:,ii][np.logical_not( mask )]=np.nan
         var_var[:,:,ii][np.logical_not( mask )]=np.nan
         dbz[:,:,ii][np.logical_not( mask )]=np.nan
         v[:,:,ii][np.logical_not( mask )]=np.nan
         w[:,:,ii][np.logical_not( mask )]=np.nan
         tk[:,:,ii][np.logical_not( mask )]=np.nan

      var_mean=np.squeeze(np.nanmax( var_mean , 1 ))
      var_kld=np.squeeze(np.nanmax( var_kld , 1 ))
      var_var=np.squeeze(np.nanmax( var_var , 1 ))
      dbz=np.squeeze(np.nanmax( dbz , 1 ))
      v=np.squeeze(np.nanmean( v , 1 ))
      w=np.squeeze(np.nanmean( w , 1 ))
      tk=np.squeeze(np.nanmean( tk , 1 ))

      #PLOTTING W AND WKLD
      minlat = np.nanmin( lat[ mask ] )
      maxlat = np.nanmax( lat[ mask ] )
      ax=axs[irow,icol]
      if ivar == 0 :
        smin = 0.0
        smax = 60.0
        var_levs=np.array([1.0,2.0,3.0,4.0])
        colorbar_axes=[0.072, 0.04, 0.25, 0.025]
        my_title = 'W'
        print( np.nanmax( var_kld ) , my_title , my_title_2 )
      if ivar == 1 :
        smin = 0.0
        smax = 20.0
        var_levs=np.array([0.2,0.4,0.6,0.8])
        var_var = var_var #* 1000.0
        colorbar_axes=[0.38, 0.04, 0.25, 0.025]
        my_title='T'
      if ivar == 2 :
        smin = 0.0
        smax = 60.0
        var_levs=np.array([1.0,2.5,5.0])
        var_kld[var_kld > 0.59 ] = 0.59
        colorbar_axes=[0.70, 0.04, 0.25, 0.025]
        my_title='Ref.'


      delta = (smax-smin)/ncolors
      clevs=np.arange(smin,smax+delta,delta)
      p=ax.contourf( tmp_lat , -np.log( levels ) ,  100*np.transpose(var_kld) 
           , clevs  , cmap=my_map)
      c=ax.contour(tmp_lat , -np.log( levels ) , np.transpose(var_var),levels=var_levs,colors='r',linewidths=2.0,linestyles='solid')
      plt.clabel(c, inline=1, fontsize=15,fmt='%1.1f')
      cdbz=ax.contour(tmp_lat , -np.log( levels ) , np.transpose(dbz),levels=[30.0],colors='k',linewidths=3.0,linestyles='dashed')

      #if ivar == 2 :
      #  c2=ax.contour(tmp_lat , -np.log( levels ) , np.transpose(v),levels=np.arange(1,10,1),colors='k',linewidths=0.5,linestyles='solid')   
      #  c3=ax.contour(tmp_lat , -np.log( levels ) , np.transpose(v),levels=np.arange(-10,-1,1),colors='k',linewidths=0.5,linestyles='dashed')
      #if ivar == 1 :
      #  c2=ax.contour(tmp_lat , -np.log( levels ) , np.transpose(tk),levels=np.arange(250,310,2),colors='k',linewidths=0.5,linestyles='solid')

      ax.set_xlim([minlat,maxlat])
      ytick=-np.log(tick_levels)
      ax.set_yticks(ytick)
      ax.set_yticklabels(levels_str,fontsize=14,color='k')
      xtick=np.array([34.9,35.0,35.1])
      ax.set_xticks(xtick)
      ax.set_xticklabels(xtick,fontsize=14,color='k')

      if icol != 0  :
         plt.setp(ax.get_yticklabels(), visible=False)
      if irow != 1  :
         plt.setp(ax.get_xticklabels(), visible=False)
      ax.tick_params(axis='both', which='both', length=0)
      ax.grid(linewidth=0.5, color='k',alpha=0.2, linestyle='--')
      ax.set_title(titles[ititle] + ' ' + my_title + ' ' + my_title_2 ,fontsize=17)

      if irow == 1 :
         cbar_ax = fig.add_axes(colorbar_axes)
         m = plt.cm.ScalarMappable(cmap=my_map)
         m.set_array(np.transpose(100*var_kld))
         m.set_clim(smin,smax)
         cb=plt.colorbar(m,cax=cbar_ax,orientation='horizontal',shrink=0.9,boundaries=np.arange(smin,smax+delta,delta))
         cb.ax.tick_params(labelsize=14)

      ititle = ititle + 1

#plt.show()
plt.savefig('Figure_kld_cross_sections.png' , format='png' , dpi=300  )
plt.close()














