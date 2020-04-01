# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016

#Este script genera una figura que muestra tiempo a tiempo el coeficiente 
#de correlacion entre la dispersion y la kld para los puntos convectivos y no convectivos.

@author:
"""
import sys
sys.path.append('../../../common_python/common_functions/')
sys.path.append('../../../common_python/common_modules/')

import numpy as np
import datetime as dt
import ctl_reader as ctlr
import os

import common_plot_functions   as cpf
import common_mask_functions   as cmf
import common_smooth_functions as csf

basedir='/home/ra001011/a03471/data/output_data/'

expnames = ['LE_D1_1km_5min','LE_D1_1km_30sec']

deltas_data = [300,300,30,60]


plot_variables      =['u','w','tk','qv'] 

#Define initial and end times using datetime module.
itime = dt.datetime(2013,7,13,5,5, 0)  #Initial time.
etime = dt.datetime(2013,7,13,5,55, 0)  #End time.

#Define the delta.
delta_plots=dt.timedelta(seconds=300)  #Plot delta t, which times will be ploted.


lat_radar=34.823
lon_radar=135.523
radar_range=50.0e3   #Radar range in meters (to define the radar mask)

#Compute the total number of times

corrcoef_rain = dict()
corrcoef_norain = dict()


ntimes=int( 1 + np.around((etime-itime).seconds / delta_plots.seconds) ) #Total number of times.

for iexp , my_exp in enumerate( expnames )  :

   if iexp == 0 :

      #=========================================================
      #  LOOP OVER FILE TYPES
      #=========================================================

      ctl_file = basedir + '/' + my_exp + '/' + '/ctl/guesgp.ctl'

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

      latlon_file = basedir + '/' + my_exp + '/' + '/latlon/latlon.grd'

      tmp=ctlr.read_data_records(latlon_file,ctl=ctl_dict,records=np.array([0,1]))
      lat=tmp[:,:,1]
      lon=tmp[:,:,0]

      radar_mask = cmf.distance_range_mask( lon_radar , lat_radar , radar_range , lon , lat )

      radar_mask = np.reshape( radar_mask , nx * ny ) 

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

      #ptime = ctime - delta_data

      #Analysis kld at time T - 1
      #m y_file=basedir + expname + ptime.strftime("%Y%m%d%H%M%S") + '/analgp/kldistance.grd'

      #kld_p=ctlr.read_data_grads(my_file,ctl_dict,masked=False)

      #Gues kld at time T
      my_file=basedir + '/' + my_exp + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/kldistance.grd'

      kld_c=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
 
      #Initialize variables.
      if it == 0 :
         corrcoef_rain[my_exp] = dict()
         corrcoef_norain[my_exp] = dict()
  
         for my_var in plot_variables :
            corrcoef_rain[my_exp][my_var]=np.zeros( ntimes )
            corrcoef_norain[my_exp][my_var]=np.zeros( ntimes )

      #kldgr=dict()
      #sprdgr=dict()
      #Compute KLD growth
      for my_var in plot_variables   :
          kld_c[my_var][ kld_c[my_var] == undef ] = np.nan 
          kld_c[my_var] = np.squeeze( np.nanmean( np.delete( kld_c[my_var] , 4 , 2 ) , 2 ) )
          kld_c[ my_var ] = csf.gaussian_smooth_2d( kld_c[my_var] , 2.0 , 5.0 )

      my_file=basedir + '/' + my_exp + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/moment0002.grd'

      sprd_c=ctlr.read_data_grads(my_file,ctl_dict,masked=False)

      for my_var in plot_variables :
         sprd_c[my_var][ sprd_c[my_var] == undef ] = np.nan
         sprd_c[my_var] = np.squeeze( np.nanmean( np.delete( sprd_c[my_var] , 4 ,2 ) , 2 ) )
         sprd_c[ my_var ] = csf.gaussian_smooth_2d( sprd_c[my_var] , 2.0 , 5.0 )

      my_file=basedir + '/' + my_exp + '/' + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/moment0001.grd'

      ens_mean=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
   
      #Compute max_dbz (we will use this to identify areas associated with clouds and convection)
      dbz = np.squeeze( ens_mean['dbz'] )
      max_dbz = np.nanmax( dbz , 2 )

      #Compute correlation coefficients:
      zbar = np.reshape( max_dbz , ( 1 , np.size( max_dbz ) ) )

      for my_var in plot_variables :

         ybar = np.reshape( kld_c[my_var] , ( 1 , np.size( kld_c[my_var] ) ) )
         xbar = np.reshape( sprd_c[my_var] , ( 1 , np.size( kld_c[my_var] ) ) )

         mask = np.logical_and( np.logical_not( np.isnan( xbar) ) , np.logical_not( np.isnan( ybar) ) )
         mask = np.logical_and( mask , zbar <= 0.0 ) 
         mask = np.logical_and( mask , radar_mask )

         corrcoef_norain[my_exp][my_var][it]=np.corrcoef(xbar[mask],ybar[mask])[0,1]

         mask = np.logical_and( np.logical_not( np.isnan( xbar) ) , np.logical_not( np.isnan( ybar) ) )
         mask = np.logical_and( mask , zbar > 30.0 ) 
         mask = np.logical_and( mask , radar_mask )

         corrcoef_rain[my_exp][my_var][it]=np.corrcoef(xbar[mask],ybar[mask])[0,1] 

         print( corrcoef_rain[my_exp][my_var][it] , corrcoef_norain[my_exp][my_var][it] )

      #print(corrcoef_norain)

      ctime = ctime + delta_plots
   
      it = it + 1

#=======================================================================================
#Plot G.R.  one panel with all the variables.
#=======================================================================================

#Plot time freq of the moments.
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


#plt.figure(1)

#Start subplots
ncols = 2
nrows = 2

irow=0
icol=0

fig, axs = plt.subplots( nrows,ncols , figsize=[10,10] )
fig.subplots_adjust(wspace=0.4,hspace=0.2,bottom=0.07,left=0.05,right=0.97,top=0.97)

titles = ['(a) - $U$','(b) - $W$','(c) - $T$','(d) - $q_v$']

colors=['b','r']

time_index=np.arange( 0 , ntimes ) * 5.0 + 5.0 #Time index in minutes.

xtick= time_index[1::2]

for ivar,my_var in enumerate(plot_variables)  :

   #Axes limits
   ax = axs[icol,irow]

   #for iexp , my_exp in enumerate( expnames ) :
   p=ax.plot( time_index , corrcoef_rain['LE_D1_1km_5min'][my_var] , color='b' ,linestyle='solid' , linewidth=2,label='5MIN - Max. Z > 30.0 dBz')
   p=ax.plot( time_index , corrcoef_rain['LE_D1_1km_30sec'][my_var] , color='r' ,linestyle='solid', linewidth=2 ,label='30SEC - Max. Z > 30.0 dBz')
   p=ax.plot( time_index , corrcoef_norain['LE_D1_1km_5min'][my_var] , color='b' ,linestyle='dashed' ,linewidth=2 ,label='5MIN - Max. Z < 0.0 dBz')
   p=ax.plot( time_index , corrcoef_norain['LE_D1_1km_30sec'][my_var] , color='r' ,linestyle='dashed' ,linewidth=2 ,label='30SEC - Max. Z < 0.0 dBz')

   if ivar == 0 :
     ax.legend()


   #Grid lines
   ax.grid(linewidth=1.0, color='k',alpha=0.5, linestyle='--')
   ax.set_xticks(xtick)
   ax.set_yticks(np.arange(-0.6,1.0,0.2))
   ax.set_ylim([-0.8,1.0])

   #Grid lines
   if icol == 1   :
      ax.set_xlabel('Time (minutes)')

   if irow == 0   :
      ax.set_ylabel('Correlation coefficient')

   ax.set_title(titles[ivar] , fontsize = 15)

   icol = icol + 1
   if icol >= ncols   :
      icol = 0
      irow = irow + 1

#plt.show()
plt.savefig( './Figure_kld_sprd_correlation.eps' , format='eps' , dpi=300)
plt.close()




