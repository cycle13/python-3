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

cases=['CASE-A','CASE-B','CASE-C']

expnames = ['/LE_D1_1km_5min/','/LE_D1_1km_5min/','/LE_D1_1km_30sec/']

cen_lon = [136.2,135.28,135.28]
ini_lat = [34.2,34.9,34.9]
end_lat = [34.4,35.15,35.15]

buffer_size=0.007  #Delta longtitude to perform zonal average.

nbv=1000

filetypes=['guesgp']   #analgp , analgz , guesgp , guesgz

#Define initial and end times using datetime module.
ctimes = [dt.datetime(2013,7,13,5,5,0),dt.datetime(2013,7,13,5,30,0),dt.datetime(2013,7,13,5,30,0)]  #Initial time.


obs_variables=['dbz','dbz','v']
upd_variables=['tk' ,'w'  ,'tk'  ]

obs_increment=[5.0 , 5.0 , 2.0]

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
nrows=3
fig, axs = plt.subplots( nrows,ncols , figsize=[10,10] , )
fig.subplots_adjust(wspace=0.15,hspace=0.0,bottom=0.05,left=0.05,right=0.97,top=0.95)


for icase,my_case in enumerate( cases ) :

   if icase == 0 :

      ctl_file = basedir + '/' + expnames[icase] + '/ctl/update_mean_diff.ctl'
      ctl_dict = ctlr.read_ctl( ctl_file )

      nx=ctl_dict['nx']
      ny=ctl_dict['nx']
      nlev=ctl_dict['nz']
      nt=int(1)             #Force the number of times to be one.
      ctl_dict['nt']=int(1) #Force the number of times to be one.

      undef=np.float32( ctl_dict['undef'] )

      ctl_file_2 = basedir + '/' + expnames[icase] + '/ctl/guesgp.ctl'
      
      ctl_dict_2 = ctlr.read_ctl( ctl_file_2 )

   #=========================================================
   #  READ
   #=========================================================
   mask = cmf.box_mask( cen_lon[icase] - buffer_size , cen_lon[icase] + buffer_size , ini_lat[icase] , end_lat[icase] , lon , lat ) 

   mean_diff = []
   kld       = []
   
   
   my_file=basedir + expnames[icase] + ctimes[icase].strftime("%Y%m%d%H%M%S") + '/guesgp/moment0001.grd'

   m01=ctlr.read_data_grads(my_file,ctl_dict_2,masked=False,undef2nan=True)

   for iv,my_var in enumerate( obs_variables ) :

      var_obs = obs_variables[iv]
      var_upd = upd_variables[iv]
                                                                                                                 
      my_file=basedir + '/' + expnames[icase] + '/' + ctimes[icase].strftime("%Y%m%d%H%M%S") +'/guesgp/update_comp_meandiff_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv])  + '.grd'
      mean_diff.append( np.delete( ctlr.read_data(my_file,ctl_dict,undef2nan=False) , 4 , 2 ) )
      my_file=basedir + '/' + expnames[icase] + '/' + ctimes[icase].strftime("%Y%m%d%H%M%S") +'/guesgp/update_comp_kld_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '.grd'
      kld.append( np.delete( ctlr.read_data(my_file,ctl_dict,undef2nan=False) , 4 , 2 ) )
      
      print(var_obs,var_upd,np.min(mean_diff[iv]),np.max(mean_diff[iv]))
      print(var_obs,var_upd,np.min(kld[iv]),np.max(kld[iv]))
      
      for ii in range( mean_diff[iv].shape[2] ) :

         mean_diff[iv][:,:,ii][np.logical_not( mask )]=np.nan
         kld[iv][:,:,ii][np.logical_not( mask )]=np.nan
         
      mean_diff[iv] = np.squeeze( np.nanmean( mean_diff[iv] , 1 ) )
      kld[iv]       = np.squeeze( np.nanmean( kld[iv] , 1 ) )
         
   dbz=np.squeeze(np.delete(m01['dbz'],4,2))
   for ii in range( dbz.shape[2] ) :
       dbz[:,:,ii][ np.logical_not(mask) ] = np.nan
   dbz = np.squeeze( np.nanmean( dbz , 1 ) )
      
   icol = icase
 
   irow = 0

   #UPDATE DIFF DBZ-W
   minlat = np.nanmin( lat[ mask ] )
   maxlat = np.nanmax( lat[ mask ] )
   
   ncolors=20
   my_map=cpf.cmap_discretize('coolwarm',ncolors)
   ax=axs[irow,icol]
   
   
   if icase == 1 or icase == 2 :
      smin = -0.6 
      smax = 0.6 
      levs = [0.05, 0.10, 0.15]
   else                        :
      smin = -7.0
      smax = 7.0
      levs = [0.25,1.0,2.0]
   
   #else                        :
   #   smin = 0.0
   #   smax = 120.0
   delta = (smax-smin)/ncolors
   clevs=np.arange(smin,smax+delta,delta)
   p=ax.contourf( tmp_lat , -np.log( levels ) ,  np.transpose(mean_diff[1]) 
           , clevs , cmap=my_map)
   cpos=ax.contour(tmp_lat , -np.log( levels ) , np.transpose(kld[1]),levels=levs,colors='g',linewidths=2.0,linestyles='solid')
   plt.clabel(cpos, inline=1, fontsize=10,fmt='%1.2f')

   cdbz=ax.contour(tmp_lat , -np.log( levels ) , np.transpose( dbz ),levels=[30.0],colors='y',linewidths=3.0,linestyles='solid')
   ax.set_xlim([minlat,maxlat])
   ytick=-np.log(tick_levels)
   ax.set_yticks(ytick)
   ax.set_yticklabels(levels_str,fontsize=10,color='k')
   if icol != 0  :
      plt.setp(ax.get_yticklabels(), visible=False)
   if irow != 2  :
      plt.setp(ax.get_xticklabels(), visible=False)
   ax.tick_params(axis='both', which='both', length=0)

   m = plt.cm.ScalarMappable(cmap=my_map)
   m.set_array(np.transpose(mean_diff))
   m.set_clim(smin,smax)
   cb=plt.colorbar(m,ax=ax,shrink=0.9,boundaries=np.arange(smin,smax+delta,delta))
   ax.grid(linewidth=0.5, color='k',alpha=0.2, linestyle='--')
   ax.text(minlat+0.01,-5.49,titles[icol+ncols*irow],fontsize=15,color='k',bbox={'facecolor':'white', 'alpha':0.8,'edgecolor':'white'})

   #UPDATE DIFF DBZ-TK
   irow=1
   ax=axs[irow,icol]
   
   if icase == 1 or icase == 2 :
      smin = -0.15 
      smax = 0.15 
      levs = [0.025,0.05, 0.10, 0.15]
   else                        :
      smin = -1.0
      smax = 1.0
      levs = [0.25,1.0,2.0]
   
   #else                        :
   #   smin = 0.0
   #   smax = 120.0
   delta = (smax-smin)/ncolors
   clevs=np.arange(smin,smax+delta,delta)
   p=ax.contourf( tmp_lat , -np.log( levels ) ,  np.transpose(mean_diff[0]) 
           , clevs , cmap=my_map)
   cpos=ax.contour(tmp_lat , -np.log( levels ) , np.transpose(kld[0]),levels=levs,colors='g',linewidths=2.0,linestyles='solid')
   plt.clabel(cpos, inline=1, fontsize=10,fmt='%1.2f')

   cdbz=ax.contour(tmp_lat , -np.log( levels ) , np.transpose( dbz ),levels=[30.0],colors='y',linewidths=3.0,linestyles='solid')
   ax.set_xlim([minlat,maxlat])
   ytick=-np.log(tick_levels)
   ax.set_yticks(ytick)
   ax.set_yticklabels(levels_str,fontsize=10,color='k')
   if icol != 0  :
      plt.setp(ax.get_yticklabels(), visible=False)
   if irow != 2  :
      plt.setp(ax.get_xticklabels(), visible=False)
   ax.tick_params(axis='both', which='both', length=0)

   m = plt.cm.ScalarMappable(cmap=my_map)
   m.set_array(np.transpose(mean_diff))
   m.set_clim(smin,smax)
   cb=plt.colorbar(m,ax=ax,shrink=0.9,boundaries=np.arange(smin,smax+delta,delta))
   ax.grid(linewidth=0.5, color='k',alpha=0.2, linestyle='--')
   ax.text(minlat+0.01,-5.49,titles[icol+ncols*irow],fontsize=15,color='k',bbox={'facecolor':'white', 'alpha':0.8,'edgecolor':'white'})


  #UPDATE DIFF DBZ-TK
   irow=2
   ax=axs[irow,icol]
   
   if icase == 1 or icase == 2 :
      smin = -0.025 
      smax = 0.025 
      levs = [0.02, 0.03, 0.04]
   else                        :
      smin = -0.25
      smax = 0.25
      levs = [0.05,0.1,0.15]
   
   #else                        :
   #   smin = 0.0
   #   smax = 120.0
   delta = (smax-smin)/ncolors
   clevs=np.arange(smin,smax+delta,delta)
   p=ax.contourf( tmp_lat , -np.log( levels ) ,  np.transpose(mean_diff[2])
           , clevs , cmap=my_map)
   cpos=ax.contour(tmp_lat , -np.log( levels ) , np.transpose(kld[2]),levels=levs,colors='g',linewidths=2.0,linestyles='solid')
   plt.clabel(cpos, inline=1, fontsize=10,fmt='%1.2f')

   cdbz=ax.contour(tmp_lat , -np.log( levels ) , np.transpose( dbz ),levels=[30.0],colors='y',linewidths=3.0,linestyles='solid')
   ax.set_xlim([minlat,maxlat])
   ytick=-np.log(tick_levels)
   ax.set_yticks(ytick)
   ax.set_yticklabels(levels_str,fontsize=10,color='k')
   if icol != 0  :
      plt.setp(ax.get_yticklabels(), visible=False)
   if irow != 2  :
      plt.setp(ax.get_xticklabels(), visible=False)
   ax.tick_params(axis='both', which='both', length=0)

   m = plt.cm.ScalarMappable(cmap=my_map)
   m.set_array(np.transpose(mean_diff))
   m.set_clim(smin,smax)
   cb=plt.colorbar(m,ax=ax,shrink=0.9,boundaries=np.arange(smin,smax+delta,delta))
   ax.grid(linewidth=0.5, color='k',alpha=0.2, linestyle='--')
   ax.text(minlat+0.01,-5.49,titles[icol+ncols*irow],fontsize=15,color='k',bbox={'facecolor':'white', 'alpha':0.8,'edgecolor':'white'})





#plt.show()
plt.savefig('Figure_Update_CrossSection_CASE_CONVECTION.eps' , format='eps' , dpi=300  )
plt.close()














