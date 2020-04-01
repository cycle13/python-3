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

cases=['CASE-A']

expnames = ['/LE_D1_1km_5min/','/LE_D1_1km_30sec/']

cen_lat = [34.3,34.3]
ini_lon = [134.8,134.8]
end_lon = [135.6,135.6]

buffer_size=0.007  #Delta longtitude to perform zonal average.

filetypes=['guesgp']   #analgp , analgz , guesgp , guesgz

#Define initial and end times using datetime module.
ctimes = [dt.datetime(2013,7,13,5,5,0),dt.datetime(2013,7,13,5,5,0)]  #Initial time.

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
my_map=cpf.cmap_discretize('Blues',ncolors)

tmp_lon = lon[90,:]
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


ncols=1
nrows=3
fig, axs = plt.subplots( nrows,ncols , figsize=[3,9] , )
fig.subplots_adjust(wspace=0.1,hspace=0.0,bottom=0.05,left=0.13,right=0.97,top=0.95)


for icase,my_case in enumerate( cases ) :

   mask = cmf.box_mask( ini_lon[icase] , end_lon[icase] , cen_lat[icase] - buffer_size , cen_lat[icase] + buffer_size , lon , lat ) 

   #=========================================================
   #  READ
   #=========================================================

   my_file=basedir + expnames[icase] + ctimes[icase].strftime("%Y%m%d%H%M%S") + '/guesgp/moment0001.grd'

   m01=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=True)

   my_file=basedir + expnames[icase] + ctimes[icase].strftime("%Y%m%d%H%M%S") + '/guesgp/moment0002.grd'

   m02=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=True)

   my_file=basedir + expnames[icase] + ctimes[icase].strftime("%Y%m%d%H%M%S") + '/guesgp/moment0003.grd'

   m03=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=True)

   my_file=basedir + expnames[icase] + ctimes[icase].strftime("%Y%m%d%H%M%S") + '/guesgp/moment0004.grd'

   m04=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=True)

   my_file=basedir + expnames[icase] + ctimes[icase].strftime("%Y%m%d%H%M%S") + '/guesgp/kldistance.grd'

   kld=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=True)

   my_file=basedir + expnames[icase] + ctimes[icase].strftime("%Y%m%d%H%M%S") + '/guesgp/bimodality_index.grd'

   bim=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=True)

   #Compute the variables that we want to plot over the cross section.

   w=np.squeeze(np.delete( m01['w'],4,2))
   wkld=np.squeeze(np.delete(kld['w'],4,2))
   tkld=np.squeeze(np.delete(kld['tk'],4,2))
   qkld=np.squeeze(np.delete(kld['qv'],4,2))
   vkld=np.squeeze(np.delete(kld['v'],4,2))
   ukld=np.squeeze(np.delete(kld['u'],4,2))
   dbz=np.squeeze(np.delete(m01['dbz'],4,2))
   wvar=np.sqrt( np.squeeze(np.delete(m02['w'],4,2)) )
   wbim=np.sqrt( np.squeeze(np.delete(bim['w'],4,2)) )


   for ii in range( np.shape( w )[2] ) :

       w[:,:,ii][np.logical_not( mask )]=np.nan  
       wkld[:,:,ii][np.logical_not( mask )]=np.nan
       tkld[:,:,ii][np.logical_not( mask )]=np.nan
       qkld[:,:,ii][np.logical_not( mask )]=np.nan
       vkld[:,:,ii][np.logical_not( mask )]=np.nan
       ukld[:,:,ii][np.logical_not( mask )]=np.nan
       dbz[:,:,ii][np.logical_not( mask )]=np.nan 
       wvar[:,:,ii][np.logical_not( mask )]=np.nan
       wbim[:,:,ii][np.logical_not( mask )]=np.nan

   w=np.squeeze(np.nanmean( w , 0 ))
   wkld=np.squeeze(np.nanmean( wkld , 0 ))
   tkld=np.squeeze(np.nanmean( tkld , 0 ))
   qkld=np.squeeze(np.nanmean( qkld , 0 ))
   vkld=np.squeeze(np.nanmean( vkld , 0 ))
   ukld=np.squeeze(np.nanmean( ukld , 0 ))
   dbz=np.squeeze(np.nanmean( dbz , 0 ))
   wvar=np.squeeze(np.nanmean( wvar , 0 ))
   wbim=np.squeeze(np.nanmean( wbim , 0 ))

   print('CASE ',cases[icase])
   print('Wkld',np.nanmin(wkld),np.nanmax(wkld),np.nanmean(wkld))
   print('tkld',np.nanmin(tkld),np.nanmax(tkld),np.nanmean(tkld))
   print('qkld',np.nanmin(qkld),np.nanmax(qkld),np.nanmean(qkld))
   print('vkld',np.nanmin(vkld),np.nanmax(vkld),np.nanmean(vkld))
   print('ukld',np.nanmin(ukld),np.nanmax(ukld),np.nanmean(ukld))
   print('wbim',np.nanmin(wbim),np.nanmax(wbim),np.nanmean(wbim))


   icol = icase
 
   irow = 0

   #PLOTTING W AND WKLD
   minlon = np.nanmin( lon[ mask ] )
   maxlon = np.nanmax( lon[ mask ] )
   ax=axs[irow]
   #if icase == 1 or icase == 2 :
   smin = 0.0
   smax = 1.5
   #else                        :
   #   smin = 0.0
   #   smax = 1.2
   delta = (smax-smin)/ncolors
   clevs=np.arange(smin,smax+delta,delta)
   p=ax.contourf( tmp_lon , -np.log( levels ) ,  np.transpose(wkld) 
           , clevs , cmap=my_map)
   cpos=ax.contour(tmp_lon , -np.log( levels ) , np.transpose(w),colors='r',levels=[0.1,0.4,0.8],linewidths=2.0,linestyles='solid')
   plt.clabel(cpos, inline=1, fontsize=10,fmt='%1.1f')
   cneg=ax.contour(tmp_lon , -np.log( levels ) , np.transpose(w), levels=[-0.8,-0.4,-0.1],colors='b',linewidths=2.0,linestyles='solid')
   plt.clabel(cneg, inline=1, fontsize=10,fmt='%1.1f')
   cdbz=ax.contour(tmp_lon , -np.log( levels ) , np.transpose(dbz),levels=[-18.0],colors='y',linewidths=3.0,linestyles='solid')
   ax.set_xlim([minlon,maxlon])
   ytick=-np.log(tick_levels)
   ax.set_yticks(ytick)
   ax.set_yticklabels(levels_str,fontsize=10,color='k')
   #if icol != 0  :
   #   plt.setp(ax.get_yticklabels(), visible=False)
   #   ax.tick_params(axis='both', which='both', length=0)
   if irow != 2  :
      plt.setp(ax.get_xticklabels(), visible=False)
   ax.tick_params(axis='both', which='both', length=0)

   m = plt.cm.ScalarMappable(cmap=my_map)
   m.set_array(np.transpose(wkld))
   m.set_clim(smin,smax)
   cb=plt.colorbar(m,ax=ax,shrink=0.9,boundaries=np.arange(smin,smax+delta,delta))
   ax.grid(linewidth=0.5, color='k',alpha=0.2, linestyle='--')
   ax.text(minlon+0.03,-5.49,titles[icol+ncols*irow],fontsize=15,color='k',bbox={'facecolor':'white', 'alpha':0.8,'edgecolor':'white'})
   if icol == 0 :
      #circle1 = plt.Circle((34.27,-5.99),0.01, color='k')
      #ax.add_artist(circle1)
      e1 = patches.Ellipse((135.275,-6.21), 0.055, 0.09,
         angle=0.0, linewidth=0.5, fill=True, zorder=2,facecolor='g',edgecolor='k')
      ax.add_patch(e1)



   #PLOTTING BIMODALITY AND W SPRD
   irow = 1 

   minlon = np.nanmin( lon[ mask ] )
   maxlon = np.nanmax( lon[ mask ] )
   ax=axs[irow]
   #if icase == 1 or icase == 2  :
   smin = 0.0
   smax = 1.0
   clevels = [3,6]
   #else                         :
   #   smin = 0.0
   #   smax = 10.0
   #   clevels = [5,10]
   delta = (smax-smin)/ncolors
   clevs=np.arange(smin,smax+delta,delta)
   p=ax.contourf( tmp_lon , -np.log( levels ) ,  np.transpose(wvar)
           , clevs , cmap=my_map)
   c=ax.contour(tmp_lon , -np.log( levels ) , 100*np.transpose(wbim),levels=clevels,colors='m',linewidths=2.0,linestyles='solid')
   plt.clabel(c, inline=1, fontsize=10,fmt='%1.0f')
   #c=ax.contour(tmp_lon , -np.log( levels ) , np.transpose(w),levels=[-1,-5,-10],colors='b',linewidths=1.0,linestyles='solid')
   cdbz=ax.contour(tmp_lon , -np.log( levels ) , np.transpose(dbz),levels=[-18.0],colors='y',linewidths=3.0,linestyles='solid')
   ax.set_xlim([minlon,maxlon])
   ytick=-np.log(tick_levels)
   ax.set_yticks(ytick)
   ax.set_yticklabels(levels_str,fontsize=10,color='k')
   #if icol != 0  :
   #   plt.setp(ax.get_yticklabels(), visible=False)
   #   ax.tick_params(axis='both', which='both', length=0)
   if irow != 2  :
      plt.setp(ax.get_xticklabels(), visible=False)
   ax.tick_params(axis='both', which='both', length=0)

   m = plt.cm.ScalarMappable(cmap=my_map)
   m.set_array(np.transpose(wvar))
   m.set_clim(smin,smax)
   cb=plt.colorbar(m,ax=ax,shrink=0.9,boundaries=np.arange(smin,smax+delta,delta))
   ax.grid(linewidth=0.5, color='k',alpha=0.2, linestyle='--')
   ax.text(minlon+0.03,-5.49,titles[icol+ncols*irow],fontsize=15,color='k',bbox={'facecolor':'white', 'alpha':0.8,'edgecolor':'white'})


   #PLOTTING QKLD AND VKLD
   irow = 2

   ax=axs[irow]
   #if icase == 1 or icase == 2 :
   smin = 0.0
   smax = 1.0
   clevels = [1,4,8]

   #else                        :
   #   smin = 0.0
   #   smax = 0.6
   #   clevels = [5,10,20]
   delta= (smax-smin)/ncolors
   clevs=np.arange(smin,smax+delta,delta)
   p=ax.contourf( tmp_lon , -np.log( levels ) ,  np.transpose(qkld)
           , clevs , cmap=my_map)
   c=ax.contour(tmp_lon , -np.log( levels ) , 100*np.transpose(vkld) , levels=clevels,colors='m',linewidths=2.0,linestyles='solid')
   plt.clabel(c, inline=1, fontsize=10,fmt='%1.0f')
   #c=ax.contour(tmp_lon , -np.log( levels ) , np.transpose(w),levels=[-1],colors='b',linewidths=1.0,linestyles='solid')
   cdbz=ax.contour(tmp_lon , -np.log( levels ) , np.transpose(dbz),levels=[-18.0],colors='y',linewidths=3.0,linestyles='solid')
   ax.set_xlim([minlon,maxlon])
   ytick=-np.log(tick_levels)
   ax.set_yticks(ytick)
   ax.set_yticklabels(levels_str,fontsize=10,color='k')
   #if icol != 0  :
   #plt.setp(ax.get_yticklabels(), visible=False)
   #ax.tick_params(axis='both', which='both', length=0)
   if irow != 2  :
      plt.setp(ax.get_xticklabels(), visible=False)
   ax.tick_params(axis='both', which='both', length=0)

   m = plt.cm.ScalarMappable(cmap=my_map)
   m.set_array(np.transpose(wvar))
   m.set_clim(smin,smax)
   cb=plt.colorbar(m,ax=ax,shrink=0.9,boundaries=np.arange(smin,smax+delta,delta))
   ax.grid(linewidth=0.5, color='k',alpha=0.2, linestyle='--')
   ax.text(minlon+0.03,-5.49,titles[icol+ncols*irow],fontsize=15,color='k',bbox={'facecolor':'white', 'alpha':0.8,'edgecolor':'white'})





#plt.show()
plt.savefig('Figure_CrossSection_CLEAR_AIR.eps' , format='eps' , dpi=300  )
plt.close()










