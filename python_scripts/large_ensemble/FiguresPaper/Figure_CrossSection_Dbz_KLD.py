# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016

@author:
"""

import sys
sys.path.append('../../../common_python/common_functions/')
sys.path.append('../../../common_python/common_modules/')

import scipy.ndimage.filters as spyf

import numpy as np
import datetime as dt
import ctl_reader as ctlr
import os

import common_plot_functions as cpf
import common_mask_functions as cmf

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib import patches
import cartopy.crs as ccrs
import common_histogram_functions as chf


nbv=1000                                      #Number of ensemble members
nbins=int(np.round(np.sqrt(nbv)))             #Number of bins used to compute the histogram
smooth_range=2                                #To smooth the histogram


basedir='/home/jruiz/share/LARGE_ENSEMBLE/output_data/home/ra001011/a03471/data/output_data/'

expnames = ['/LE_D1_1km_5min/','/LE_D1_1km_30sec/']

cen_lon = [135.28]
#cen_lon = [135.8]
ini_lat = [34.9]
end_lat = [35.25]

buffer_size=0.007  #Delta longtitude to perform zonal average.

filetype=['guesgp']   #analgp , analgz , guesgp , guesgz

#Define initial and end times using datetime module.
ctimes = [dt.datetime(2013,7,13,5,30,0)]  #Initial time.

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


ncolors=10
my_map=cpf.cmap_discretize('Blues',ncolors)

tmp_lat = lat[:,90]
levels=ctl_dict['vlevels']

tick_levels=[1000,850,700,500,300,200]
levels_str=list()
levels=np.delete(levels,4,axis=0)
levels[3]=850.0

#titles=['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)']

#Get the level string list.
levels_str=[]
for ilev in tick_levels  :
   levels_str.append( str(int(ilev)) )


#ncols=3
#nrows=3
fig= plt.figure( figsize=[15,8] , )
#fig.subplots_adjust(wspace=0.1,hspace=0.0,bottom=0.05,left=0.05,right=0.97,top=0.95)


#=========================================================
#  PLOT CROSS-SECTIONS OF REFLECTIVITY AND W.
#=========================================================

for iexp,my_exp in enumerate( expnames ) :

   mask = cmf.box_mask( cen_lon[0] - buffer_size , cen_lon[0] + buffer_size , ini_lat[0] , end_lat[0] , lon , lat ) 

   #=========================================================
   #  READ
   #=========================================================

   my_file=basedir + expnames[iexp] + ctimes[0].strftime("%Y%m%d%H%M%S") + '/' + filetype[0] + '/moment0001.grd'
   m01=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=True)
   w=np.squeeze(np.delete( m01['w'],4,2))
   dbz3d=np.squeeze(np.delete(m01['dbz'],4,2))

   for ii in range( np.shape( w )[2] ) :
       w[:,:,ii][np.logical_not( mask )]=np.nan  
       dbz3d[:,:,ii][np.logical_not( mask )]=np.nan 
   w=np.squeeze(np.nanmax( w , 1 ))
   dbz=np.squeeze(np.nanmax( dbz3d , 1 ))

   if iexp == 0 :
       ax = plt.axes([0.035, 0.12, 0.20, 0.84])
       my_title='(a) 5MIN'
       #For the 5MIN experiment find the location with the maximum bimodality in W
       my_file=basedir + expnames[iexp] + ctimes[0].strftime("%Y%m%d%H%M%S") + '/' + filetype[0] + '/bimodality_index.grd'
       ens_bimodality=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
       my_file=basedir + expnames[iexp] + ctimes[0].strftime("%Y%m%d%H%M%S") + '/' + filetype[0] + '/kldistance.grd'
       ens_kld=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
       wkld=np.squeeze(np.delete( ens_kld['w'],4,2))
       wbim=np.squeeze(np.delete( ens_bimodality['w'],4,2))
       wkld[ dbz3d < 30.0 ] = np.nan
       wbim[ dbz3d < 30.0 ] = np.nan
       for ii in range( np.shape( wbim  )[2] ) :
           wbim[:,:,ii][np.logical_not(mask)]=np.nan
           wkld[:,:,ii][np.logical_not(mask)]=np.nan
           
       max_loc=np.nanargmax( wkld )
       #print( iexp , np.nanmax( wkld ) )
       [tmpnx , tmpny , tmpnz ]=np.shape( np.squeeze( wkld ) )
       [hist_x_5MIN,hist_y_5MIN,hist_z_5MIN]=np.unravel_index( max_loc , (tmpnx,tmpny,tmpnz) )
       lon_hist=lon[hist_x_5MIN,hist_y_5MIN]
       lat_hist=lat[hist_x_5MIN,hist_y_5MIN]
       lev_hist=levels[hist_z_5MIN]
       print('xyz 5MIN', hist_x_5MIN , hist_y_5MIN , hist_z_5MIN )
       #print( iexp , wkld[hist_x_5MIN,hist_y_5MIN,hist_z_5MIN] )
       #print( lon_hist , lat_hist , lev_hist )
   if iexp == 1 :
       ax = plt.axes([0.245, 0.12, 0.20, 0.84])
       plt.setp(ax.get_yticklabels(), visible=False)
       my_title='(b) 30SEC'
       #For the 30SEC experiment find the location with the maximum bimodality in W
       my_file=basedir + expnames[iexp] + ctimes[0].strftime("%Y%m%d%H%M%S") + '/' + filetype[0] + '/bimodality_index.grd'
       ens_bimodality=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
       my_file=basedir + expnames[iexp] + ctimes[0].strftime("%Y%m%d%H%M%S") + '/' + filetype[0] + '/kldistance.grd'
       ens_kld=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
       wkld=np.squeeze(np.delete( ens_kld['w'],4,2))
       wbim=np.squeeze(np.delete( ens_bimodality['w'],4,2))
       wkld[ dbz3d < 30.0 ] = np.nan
       wbim[ dbz3d < 30.0 ] = np.nan
       for ii in range( np.shape( wbim  )[2] ) :
           wbim[:,:,ii][np.logical_not(mask)]=np.nan
           wkld[:,:,ii][np.logical_not(mask)]=np.nan

       max_loc=np.nanargmax( wkld )
       #print( iexp , np.nanmax( wkld ) )
       [tmpnx , tmpny , tmpnz ]=np.shape( np.squeeze( wkld ) )
       [hist_x_30SEC,hist_y_30SEC,hist_z_30SEC]=np.unravel_index( max_loc , (tmpnx,tmpny,tmpnz) )
       #print( iexp , wkld[hist_x_30SEC,hist_y_30SEC,hist_z_30SEC] )
       print('xyz 30SEC', hist_x_30SEC , hist_y_30SEC , hist_z_30SEC )
       lon_hist=lon[hist_x_30SEC,hist_y_30SEC]
       lat_hist=lat[hist_x_30SEC,hist_y_30SEC]
       lev_hist=levels[hist_z_30SEC]
       #print( lon_hist , lat_hist , lev_hist )

   #PLOTTING W AND DBZ
   minlat = np.nanmin( lat[ mask ] )
   maxlat = np.nanmax( lat[ mask ] )
   smin=0
   smax=70
   ncols=70
   delta = ( smax - smin )/ncols
   clevs=np.arange(smin,smax+delta,delta)
   my_map = cpf.cmap_discretize('gist_ncar',100)
   p=ax.contourf( tmp_lat , -np.log( levels )  ,  np.transpose( dbz ) , clevs  ,cmap=my_map)
   #p=ax.pcolor( tmp_lat , -np.log( levels )  ,  np.transpose( dbz )  ,cmap=my_map)
   cpos=ax.contour(tmp_lat , -np.log( levels ) , np.transpose(w),levels=[1,2.5,5,10],colors='k',linewidths=2.0,linestyles='solid')
   cpos=ax.contour(tmp_lat , -np.log( levels ) , np.transpose(w),levels=[-2.5,-1.0],colors='k',linewidths=2.0,linestyles='dashed')
   plt.clabel(cpos, inline=1, fontsize=14,fmt='%1.0f')

   ax.set_xlim([minlat,maxlat])
   ax.set_ylim([-np.log(levels[0]),-np.log(levels[-1])+0.01])
   ytick=-np.log(tick_levels)
   ax.set_yticks(ytick)
   ax.set_yticklabels(levels_str,fontsize=14,color='k')
   ax.set_xticks(np.array([34.9,35.0,35.1,35.2]))
   ax.set_xticklabels(np.array([34.9,35.0,35.1,35.2]),fontsize=14,color='k')
   ax.tick_params(axis='both', which='both', length=0)

   ax.grid(linewidth=0.5, color='k',alpha=0.2, linestyle='--')
   #ax.text(minlat+0.0,-5.35,my_title,fontsize=17,color='k',bbox={'facecolor':'white', 'alpha':0.0,'edgecolor':'white'})
   ax.set_title(my_title,fontsize=17)

   e1 = patches.Ellipse((lat_hist,-np.log(lev_hist)), 0.015, 0.042,
   angle=0.0, linewidth=0.5, fill=True, zorder=2,facecolor='r',edgecolor='k')
   ax.add_patch(e1)


cbar_ax = fig.add_axes([0.01, 0.05, 0.44, 0.04])
m = plt.cm.ScalarMappable(cmap=my_map )
m.set_array(np.transpose(dbz))
m.set_clim(smin,smax)
cb=plt.colorbar(m,cax=cbar_ax,boundaries=np.arange(smin,smax+delta,delta),orientation='horizontal')
cb.ax.tick_params(labelsize=14)

#=========================================================
#  PLOT HISTOGRAMS AT THE LOCATION OF MAXIMUM BIMODALITY
#=========================================================

for iexp,my_exp in enumerate( expnames ) :

    my_file=basedir + expnames[iexp] + ctimes[0].strftime("%Y%m%d%H%M%S") + '/' + filetype[0]  + '/moment0001.grd'
    my_file=basedir + expnames[iexp] + ctimes[0].strftime("%Y%m%d%H%M%S") + '/' + filetype[0]  + '/moment0001.grd'

    hist_file = basedir + expnames[iexp] + ctimes[0].strftime("%Y%m%d%H%M%S") + '/' + filetype[0]  + '/histogram.grd'
    max_file  = basedir + expnames[iexp] + ctimes[0].strftime("%Y%m%d%H%M%S") + '/' + filetype[0]  + '/ensmax.grd'
    min_file  = basedir + expnames[iexp] + ctimes[0].strftime("%Y%m%d%H%M%S") + '/' + filetype[0]  + '/ensmin.grd'

    hist=chf.read_histogram(hist_file,max_file,min_file,nx,ny,nbins,ctl_dict,dtypein='i2')
    my_file=basedir + expnames[iexp] + ctimes[0].strftime("%Y%m%d%H%M%S") + '/' + filetype[0]  + '/moment0001.grd'
    mean=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
    my_file=basedir + expnames[iexp] + ctimes[0].strftime("%Y%m%d%H%M%S") + '/' + filetype[0]  + '/moment0002.grd'
    variance=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
    my_file=basedir + expnames[iexp] + ctimes[0].strftime("%Y%m%d%H%M%S") + '/' + filetype[0] + '/kldistance.grd'
    kld=ctlr.read_data_grads(my_file,ctl_dict,masked=False)

    if iexp == 0 :
        ax1= plt.axes([0.485, 0.54, 0.16, 0.42])
        ax2= plt.axes([0.66, 0.54, 0.16, 0.42])
        ax3= plt.axes([0.835, 0.54, 0.16, 0.42])
        titles=['(c)','(d)','(e)']
        titles2='5MIN'
        #plt.setp(ax1.get_xticklabels(), visible=False)
        #plt.setp(ax2.get_xticklabels(), visible=False)
        #plt.setp(ax3.get_xticklabels(), visible=False)
        hist_x = hist_x_5MIN
        hist_y = hist_y_5MIN
        hist_z = hist_z_5MIN

    if iexp == 1 :
        ax1= plt.axes([0.485, 0.035, 0.16, 0.42])
        ax2= plt.axes([0.66, 0.035, 0.16, 0.42])
        ax3= plt.axes([0.835, 0.035, 0.16, 0.42])
        titles=['(f)','(g)','(h)']
        titles2='30SEC'
        hist_x = hist_x_30SEC
        hist_y = hist_y_30SEC
        hist_z = hist_z_30SEC

    #W
    ax=ax1
    var='w'
    #title=titles[0]
    xmin=-3.0
    xmax=15.0
    my_title=titles[0]+' W - ' + titles2
    hist_min=np.squeeze(np.delete( hist[var]['minval'],4,2))[hist_x,hist_y,hist_z]
    hist_max=np.squeeze(np.delete( hist[var]['maxval'],4,2))[hist_x,hist_y,hist_z]
    hist_mean=np.squeeze(np.delete( mean[var],4,2))[hist_x,hist_y,hist_z]
    hist_var=np.squeeze(np.delete( variance[var],4,2))[hist_x,hist_y,hist_z]
    hist_kld=np.squeeze(np.delete( kld[var],4,2))[hist_x,hist_y,hist_z]

    hist_delta=(hist_max-hist_min)/nbins
    hist_range=hist_min + hist_delta / 2 + hist_delta *  np.arange(0,nbins,1)
    hist_bars=np.squeeze(np.delete( hist[var]['hist'],4,2))[hist_x,hist_y,hist_z,:] / np.sum( np.squeeze(np.delete( hist[var]['hist'],4,2))[hist_x,hist_y,hist_z,:] )
    hist_bars_s=spyf.uniform_filter(hist_bars, size=smooth_range, mode='reflect')
    #Get the Gaussian fit to the histogram
    hist_gaussfit=hist_delta*(1/(np.sqrt(2*np.pi*hist_var)))*np.exp( -0.5*np.power(hist_range-hist_mean,2)/hist_var )
    hist_label= ' KLD=' + '%1.3f' % hist_kld
    the_bars=ax.bar( hist_range , hist_bars_s , width = hist_delta ,color='grey' )
    the_lines=ax.plot( hist_range , hist_gaussfit , 'k--',linewidth=4)

    xmin=np.floor(np.min(hist_range)) -1.0
    xmax=xmax=np.ceil(np.max(hist_range)) +1.0
    print(xmin,xmax)
    ax.set_xlim(xmin=xmin , xmax=xmax )
    ax.set_ylim(ymin=0.0, ymax=0.15 )
    ax.set_xticks( 2.5 + np.arange( xmin , xmax , 4.0 )[:-1] )
    ax.set_xticklabels( 2.5 + np.arange( xmin , xmax , 4.0 )[:-1] ,fontsize=14)
    ax.set_yticks( np.array([0,0.05,0.1,0.15]))
    ax.set_yticklabels( np.array([0,0.05,0.1,0.15]),fontsize=14)
    #ax.text( xmin + 0.5, 0.135 ,title ,fontsize=17,color='k',bbox={'facecolor':'white', 'alpha':0.0,'edgecolor':'white'})
    ax.text( xmin + 0.25 , 0.14 ,hist_label ,fontsize=17,color='k',bbox={'facecolor':'white', 'alpha':0.0,'edgecolor':'white'})
    ax.set_title( my_title , fontsize=17)
    #Temperature
    #print( mean['tk'][hist_x,hist_y,hist_z] )
    ax=ax2
    var='tk'
    title=titles[1]
    my_title=titles[1]+' T - ' + titles2
    plt.setp(ax.get_yticklabels(), visible=False)
    hist_min=np.squeeze(np.delete( hist[var]['minval'],4,2))[hist_x,hist_y,hist_z]
    hist_max=np.squeeze(np.delete( hist[var]['maxval'],4,2))[hist_x,hist_y,hist_z]
    hist_mean=np.squeeze(np.delete( mean[var],4,2))[hist_x,hist_y,hist_z]
    hist_var=np.squeeze(np.delete( variance[var],4,2))[hist_x,hist_y,hist_z]
    hist_kld=np.squeeze(np.delete( kld[var],4,2))[hist_x,hist_y,hist_z]
    hist_delta=(hist_max-hist_min)/nbins
    hist_range=hist_min + hist_delta / 2 + hist_delta *  np.arange(0,nbins,1)
    hist_bars=np.squeeze(np.delete( hist[var]['hist'],4,2))[hist_x,hist_y,hist_z,:] / np.sum( np.squeeze(np.delete( hist[var]['hist'],4,2))[hist_x,hist_y,hist_z,:] )
    hist_bars_s=spyf.uniform_filter(hist_bars, size=smooth_range, mode='reflect')

    #Get the Gaussian fit to the histogram
    hist_gaussfit=hist_delta*(1/(np.sqrt(2*np.pi*hist_var)))*np.exp( -0.5*np.power(hist_range-hist_mean,2)/hist_var )
    hist_label= ' KLD=' + '%1.3f' % hist_kld
    the_bars=ax.bar( hist_range , hist_bars_s , width = hist_delta ,color='grey' )
    the_lines=ax.plot( hist_range , hist_gaussfit , 'k--',linewidth=4)
    xmin=np.floor(np.min(hist_range)) -1.0
    xmax=xmax=np.ceil(np.max(hist_range)) +1.0
    print(xmin,xmax)
    ax.set_xlim(xmin=xmin , xmax=xmax )
    ax.set_ylim(ymin=0.0, ymax=0.15 )
    ax.set_xticks( 0.5 + np.arange( xmin , xmax , 2.0 ) )
    ax.set_xticklabels( 0.5 + np.arange( xmin , xmax , 2.0 ) ,fontsize=14)
    ax.set_ylim(ymin=0.0, ymax=0.15 )
    ax.set_yticks( np.array([0,0.05,0.1,0.15]))
    ax.set_yticklabels( np.array([0,0.05,0.1,0.15]),fontsize=14)
    #ax.text( xmin + 0.2, 0.135 ,title ,fontsize=17,color='k',bbox={'facecolor':'white', 'alpha':0.0,'edgecolor':'white'})
    ax.text( xmin + 0.1 , 0.14 ,hist_label ,fontsize=17,color='k',bbox={'facecolor':'white', 'alpha':0.0,'edgecolor':'white'})
    ax.set_title( my_title , fontsize=17)
    #Dbz
    #print( mean['dbz'][hist_x,hist_y,hist_z] )
    ax=ax3
    var='dbz'
    title=titles[2]
    xmin=30
    xmax=60
    my_title=titles[2]+' Ref. - ' + titles2
    plt.setp(ax.get_yticklabels(), visible=False)
    hist_min=np.squeeze(np.delete( hist[var]['minval'],4,2))[hist_x,hist_y,hist_z]
    hist_max=np.squeeze(np.delete( hist[var]['maxval'],4,2))[hist_x,hist_y,hist_z]
    hist_mean=np.squeeze(np.delete( mean[var],4,2))[hist_x,hist_y,hist_z]
    hist_var=np.squeeze(np.delete( variance[var],4,2))[hist_x,hist_y,hist_z]
    hist_kld=np.squeeze(np.delete( kld[var],4,2))[hist_x,hist_y,hist_z]

    hist_delta=(hist_max-hist_min)/nbins
    hist_range=hist_min + hist_delta / 2 + hist_delta *  np.arange(0,nbins,1)
    hist_bars=np.squeeze(np.delete( hist[var]['hist'],4,2))[hist_x,hist_y,hist_z,:] / np.sum( np.squeeze(np.delete( hist[var]['hist'],4,2))[hist_x,hist_y,hist_z,:] )
    hist_bars_s=spyf.uniform_filter(hist_bars, size=smooth_range, mode='reflect')
    #Get the Gaussian fit to the histogram
    hist_gaussfit=hist_delta*(1/(np.sqrt(2*np.pi*hist_var)))*np.exp( -0.5*np.power(hist_range-hist_mean,2)/hist_var )
    hist_label= ' KLD=' + '%1.3f' % hist_kld
    the_bars=ax.bar( hist_range , hist_bars_s , width = hist_delta ,color='grey' )
    the_lines=ax.plot( hist_range , hist_gaussfit , 'k--',linewidth=4)
    ax.set_xlim(xmin=xmin, xmax=xmax+0.01 )
    ax.set_ylim(ymin=0.0, ymax=0.15 )
    ax.set_xticks( 5.0 + np.arange( xmin , xmax , (xmax-xmin)/3.0 ) )
    ax.set_xticklabels( 5.0 + np.arange( xmin , xmax , (xmax-xmin)/3.0 ),fontsize=14)
    ax.set_yticks( np.array([0,0.05,0.1,0.15]))
    ax.set_yticklabels( np.array([0,0.05,0.1,0.15]),fontsize=14)
    #ax.text( xmin + 1.5, 0.135 ,title ,fontsize=17,color='k',bbox={'facecolor':'white', 'alpha':0.0,'edgecolor':'white'})
    ax.text( xmin + 0.7 , 0.14 ,hist_label ,fontsize=17,color='k',bbox={'facecolor':'white', 'alpha':0.0,'edgecolor':'white'})
    ax.set_title( my_title , fontsize=17)

plt.show()
plt.savefig('Figure_CrossSection_dBZ_KLD.png' , format='png' , dpi=300  )
plt.close()














