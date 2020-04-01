# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016

#Este script calcula:
-El aumento en KLD integrado verticalmente para cada variable. (durante el paso de pronostico)
-El aumento en el spread en la energia estatica humeda.

Objetivo: Explorar la relacion entre la no-gaussianidad y otras parametros como el crecimiento de las perturbaciones.

@author:
"""
import sys
sys.path.append('../../common_python/common_functions/')
sys.path.append('../../common_python/common_modules/')

import numpy as np
import datetime as dt
import ctl_reader as ctlr
import os

import common_plot_functions   as cpf
import common_mask_functions   as cmf
import common_smooth_functions as csf

basedir='/home/ra001011/a03471/data/output_data/'

expname = '/LE_D1_1km_30sec/'
delta_data =dt.timedelta(seconds=30)   #Experiment data delta t (to compute growth rate)

figname='Kld_growth'

plot_variables      =['u','v','w','tk','qv','dbz'] 

#Define initial and end times using datetime module.
itime = dt.datetime(2013,7,13,5,5, 0)  #Initial time.
etime = dt.datetime(2013,7,13,5,55, 0)  #End time.

#Define the delta.
delta_plots=dt.timedelta(seconds=300)  #Plot delta t, which times will be ploted.

#Compute the total number of times

#=========================================================
#  LOOP OVER FILE TYPES
#=========================================================

ctl_file = basedir + expname + '/ctl/guesgp.ctl'

plotbasedir=basedir + expname + '/plots/guesgp/'

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

latlon_file = basedir + expname + '/latlon/latlon.grd'

tmp=ctlr.read_data_records(latlon_file,ctl=ctl_dict,records=np.array([0,1]))
lat=tmp[:,:,1]
lon=tmp[:,:,0]


#=========================================================
#  START TIME LOOP
#=========================================================

it=0

ctime=itime

while ( ctime <= etime ):

   print( ctime )

   print ( 'Reading the moments and computing growth rate ')

   #=========================================================
   #  READ THE DATA
   #=========================================================

   ptime = ctime - delta_data

   #Analysis kld at time T - 1
   my_file=basedir + expname + ptime.strftime("%Y%m%d%H%M%S") + '/analgp/kldistance.grd'

   kld_p=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=True)

   #Gues kld at time T
   my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/kldistance.grd'

   kld_c=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=True)

   kldgr=dict()
   sprdgr=dict()
   #Compute KLD growth
   for my_var in plot_variables   :
       kldgr[my_var] = np.squeeze( np.nanmean( np.delete(kld_c[my_var] - kld_p[my_var] ,4,2) ,2)  ) / delta_data.total_seconds()

       #Smooth the growing rate.
       #kldgr[my_var]=csf.gaussian_smooth_2d( kldgr[my_var] , 2.0 , 5.0 )

   #Gues spread at time t
   my_file=basedir + expname + ptime.strftime("%Y%m%d%H%M%S") + '/analgp/moment0002.grd'

   sprd_p=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=True)

   #Analysis spread at time t
   my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/moment0002.grd'

   sprd_c=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=True)

   for my_var in plot_variables :
       sprdgr[my_var] = np.squeeze( np.nanmean( np.delete( sprd_c[my_var] - sprd_p[my_var] ,4,2) ,2)  ) / delta_data.total_seconds()
       
       #Smooth the growing rate.
       #sprdgr[my_var]=csf.gaussian_smooth_2d( sprdgr[my_var] , 2.0 , 5.0 )



   #Compute the total moist energy associated with the spread.
   #cp = 1004 ; tr = 280 ; L = 2264.705e3 ; Rd=287 ; pr=1e5 
   #moist_total_energy_c = np.squeeze( np.nanmean( sprd_c['u'] + sprd_c['v'] + sprd_c['tk']*(cp/tr) + sprd_c['qv']*(np.power(L,2)/(cp*tr)) ,2) ) + np.squeeze( ( Rd * tr / np.power( pr , 2 ) )*sprd_c['slp'] )

   #Read the ensemble mean to get the information from the storm location.
   my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/moment0001.grd'

   ens_mean=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
   
   #Compute max_dbz (we will use this to identify areas associated with clouds and convection)
   max_dbz = np.squeeze( np.nanmax(ens_mean['dbz'],2) )

   #sprdgr = ( moist_total_energy_c - moist_total_energy_p )
   #sprdgr=csf.gaussian_smooth_2d( sprdgr , 2.0 , 5.0 )

   #So far we have kld growth integrated in the vertical and for each variable.
   #sprd growth in the moist total energy. 

   #=======================================================================================
   #Plot G.R.  one panel with all the variables.
   #=======================================================================================

   #Plot time freq of the moments.
   import matplotlib
   import matplotlib.pyplot as plt
   import matplotlib.ticker as mticker
   import cartopy.crs as ccrs


   fig=plt.figure(1,figsize=[10,6.5])
   #plt.figure(1)

   #Start subplots
   ncols = 3
   nrows = 2
   icoldelta = 1.0/ncols
   irowdelta = 1.0/nrows
   hmargin=0.0
   vmargin=0.0
   hoffset=-0.06
   voffset=0.1

   icol = 1
   irow = nrows

   xtick=[134.5,135,135.5,136,136.5,137]
   ytick=[34,34.5,35,35.5]
   axesrange=[134.97,136.09,34.36,35.30]
   titles = ['(a)','(b)','(c)','(d)','(e)','(f)']



   for ivar,my_var in enumerate(plot_variables)  :

      #kld growth rate
      print('KLD analysis growth rate for var ',my_var)
      print(np.nanmin( kldgr[my_var]),np.nanmax(kldgr[my_var])) 

      print('SPRD analysis growth rate for var ')
      print(np.nanmin( sprdgr[my_var] ),np.nanmax( sprdgr[my_var] ))


      #smin = np.min([np.min(kldgr[my_var]),-np.max(kldgr[my_var])] )
      #smax = np.max([np.max(kldgr[my_var]),-np.min(kldgr[my_var])] )
      #my_range = smax - smin 
      my_range = np.std( np.reshape(  kldgr[my_var] , (nx*ny) )  )
      smin = -10.0  
      smax = 10.0 
      #Axes limits
      my_axes = [icoldelta*(icol-1)+hmargin+hoffset,irowdelta*(irow-1)+vmargin+voffset,irowdelta-2*hmargin,icoldelta-2*vmargin]

      ax = plt.axes( my_axes , facecolor=None , projection=ccrs.PlateCarree() )

      #The pcolor
      my_map=cpf.cmap_discretize('seismic',21)
      p=ax.pcolor(lon , lat ,   np.squeeze( kldgr[my_var] / my_range ) ,
        transform=ccrs.PlateCarree(),vmin=smin, vmax=smax ,cmap=my_map )

      ax.set_extent( axesrange , ccrs.PlateCarree())

      #Grid lines
      gl=ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                linewidth=1.0, color='k',alpha=0.5, linestyle='--')

      #Grid and ticks
      gl.xlocator=mticker.FixedLocator(xtick)
      gl.ylocator=mticker.FixedLocator(ytick)

      #Coastline plot
      ax.coastlines('10m',linewidth=1.0)

      #Contour map
      #if my_var == 'u' or my_var == 'v' or my_var == 'w' :
      #   smax = 20
      #if my_var == 'tk'  :
      #   smax = 2
      #if my_var == 'dbz' :
      #   smax = 1000
      #if my_var == 'qv'  :
      #   smax = 2.0e-6
      #smin = -np.max([np.max(sprdgr[my_var]),-np.min(sprdgr[my_var])] )
      #smax = np.max([np.max(sprdgr[my_var]),-np.min(sprdgr[my_var])] )
      my_range = np.std( np.reshape(  sprdgr[my_var] , (nx*ny) )  )

      #sdelta = (smax-smin)/2
      #my_range = smax - smin

      c=ax.contour( lon , lat , sprdgr[my_var] / my_range , levels = [(0.5 , 3.0)] ,colors='k' ,linestyles='solid')

      #smin = -np.max([np.max(sprdgr[my_var]),-np.min(sprdgr[my_var])] )
      #smax = 0
      #sdelta = (smax-smin)/2
      
      c=ax.contour( lon , lat , sprdgr[my_var] / my_range , levels = [(-3.0,-0.5)] ,colors='k' ,linestyles='dashed')



      c2=ax.contour( lon , lat ,  max_dbz ,levels = [-100,20] ,colors='#FFA500' , linewidths=3.0)

      #Texts
      sp=1   #Number of significant digits to be retained.
      my_max = np.nanmax( kldgr[my_var] )
      my_max=np.around( my_max / np.power(10, np.floor( np.log10( np.abs(my_max))) - sp ) ) * np.power(10, np.floor( np.log10( np.abs(my_max))) - sp )
      my_min = np.nanmin( kldgr[my_var] )
      #print( my_min , np.log10( my_min ) )
      my_min=np.around( my_min / np.power(10, np.floor( np.log10( np.abs(my_min))) - sp ) ) * np.power(10, np.floor( np.log10(np.abs(my_min))) - sp )
      
      print( my_min , my_max )
      #plt.text(135.2, 34.6,'Max=' + str(my_max) + ' Min=' + str(my_min),horizontalalignment='right', verticalalignment='top', transform=ccrs.PlateCarree() )
      plt.text(135.0 ,34.45 ,'Max=' + str(my_max) + ' Min=' + str(my_min),withdash=False ,ha="left", va="top",bbox=dict(boxstyle="round",ec=(1., 0.5, 0.5),fc=(1., 0.8, 0.8 )))


      #Axes label and title
      if irow == 1   :
         plt.xlabel('Latitude')

      if icol == 1   :
         plt.ylabel('Longitude')

      plt.title(titles[ivar] + ' ' + my_var,fontsize = 15)

      gl.xlabel_style = {'size': 12, 'color': 'k' }
      gl.ylabel_style = {'size': 12, 'color': 'k' }
      gl.xlabels_top = False
      gl.ylabels_right = False

      icol = icol + 1
      if icol > ncols   :
         icol = 1
         irow = irow - 1

      #Colorbar
      if ivar == 0  :
         cbar_ax = fig.add_axes([0.17, 0.52, 0.7, 0.03])
         cb=fig.colorbar(p, cax=cbar_ax , orientation='horizontal')
         cb.ax.tick_params(labelsize=10) 
      #cb=plt.colorbar(p)
      #cb.ax.tick_params(labelsize=10)


   plt.show()
   plt.savefig( plotbasedir + '/' + figname + '_' + my_var + '_' + ctime.strftime("%Y%m%d%H%M%S") + '.png' , format='png' , dpi=300)
   plt.close()


   ctime = ctime + delta_plots



