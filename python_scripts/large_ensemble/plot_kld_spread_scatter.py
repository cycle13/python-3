# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016

#Este script calcula:
Un scatter plot con la relacion entre el crecimiento del ensemble spread y la no gaussianidad.

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

expname = '/LE_D1_1km_30sec_nospinup/'
delta_data =dt.timedelta(seconds=30)   #Experiment data delta t (to compute growth rate)


figname='Kld_spread_scatter'

plot_variables      =['u','v','w','tk','qv','dbz'] 

#Define initial and end times using datetime module.
itime = dt.datetime(2013,7,13,5,5,  0)  #Initial time.
etime = dt.datetime(2013,7,13,5,55, 0)  #End time.

#Define the delta.
delta_plots=dt.timedelta(seconds=300)  #Plot delta t, which times will be ploted.

#Compute the total number of times

lat_radar=34.823
lon_radar=135.523
radar_range=50.0e3   #Radar range in meters (to define the radar mask)


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

radar_mask = cmf.distance_range_mask( lon_radar , lat_radar , radar_range , lon , lat )

radar_mask = np.reshape( radar_mask , ( nx * ny ) )

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
   #my_file=basedir + expname + ptime.strftime("%Y%m%d%H%M%S") + '/analgp/kldistance.grd'

   #kld_p=ctlr.read_data_grads(my_file,ctl_dict,masked=False)

   #Gues kld at time T
   my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/kldistance.grd'

   kld_c=ctlr.read_data_grads(my_file,ctl_dict,masked=False)

   #kldgr=dict()
   #sprdgr=dict()
   #Compute KLD growth
   for my_var in plot_variables   :
       kld_c[my_var][ kld_c[my_var] == undef ] = np.nan 
       kld_c[my_var] = np.squeeze( np.nanmean( np.delete( kld_c[my_var],4,2) , 2 ) )
       kld_c[ my_var ] = csf.gaussian_smooth_2d( kld_c[my_var] , 2.0 , 5.0 )

   #Analysis kld at time T - 1
   #my_file=basedir + expname + ptime.strftime("%Y%m%d%H%M%S") + '/analgp/moment0002.grd'

   #sprd_p=ctlr.read_data_grads(my_file,ctl_dict,masked=False)

   #Compute the total moist energy associated with the spread.
   #cp = 1004 ; tr = 280 ; L = 2264.705e3 ; Rd=287 ; pr=1e5 
   #moist_total_energy_p = np.squeeze( np.nanmean( sprd_p['u'] + sprd_p['v'] + sprd_p['tk']*(cp/tr) + sprd_p['qv']*(np.power(L,2)/(cp*tr)) ,2) ) + np.squeeze( ( Rd * tr / np.power( pr , 2 ) )*sprd_p['slp'] )

   #Gues kld at time T
   my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/moment0002.grd'

   sprd_c=ctlr.read_data_grads(my_file,ctl_dict,masked=False)

   for my_var in plot_variables :
   #    sprd_p[my_var][ sprd_p[my_var] == undef ] = np.nan
       sprd_c[my_var][ sprd_c[my_var] == undef ] = np.nan
       sprd_c[my_var] = np.squeeze( np.nanmean( np.delete( sprd_c[my_var] ,4,2) , 2 ) )
       sprd_c[ my_var ] = csf.gaussian_smooth_2d( sprd_c[my_var] , 2.0 , 5.0 )

   #Compute the total moist energy associated with the spread.
   #cp = 1004 ; tr = 280 ; L = 2264.705e3 ; Rd=287 ; pr=1e5 
   #moist_total_energy_c = np.squeeze( np.nanmean( sprd_c['u'] + sprd_c['v'] + sprd_c['tk']*(cp/tr) + sprd_c['qv']*(np.power(L,2)/(cp*tr)) ,2) ) + np.squeeze( ( Rd * tr / np.power( pr , 2 ) )*sprd_c['slp'] )

   #Read the ensemble mean to get the information from the storm location.
   my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/moment0001.grd'

   ens_mean=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
   
   #Compute max_dbz (we will use this to identify areas associated with clouds and convection)
   dbz = np.squeeze( ens_mean['dbz'] )
   max_dbz = np.nanmax( dbz , 2 )

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


   plt.figure(1,figsize=[11.5,5])
   #plt.figure(1)

   #Start subplots
   ncols = 3
   nrows = 2
   icoldelta = 1.0/ncols
   irowdelta = 1.0/nrows
   hmargin=0.125
   vmargin=0.0
   hoffset=-0.052
   voffset=0.1

   icol = 1
   irow = nrows

   xtick=[134.5,135,135.5,136,136.5,137]
   ytick=[34,34.5,35,35.5]
   axesrange=[134.97,136.09,34.36,35.30]
   titles = ['(a)','(b)','(c)','(d)','(e)','(f)']

   for ivar,my_var in enumerate(plot_variables)  :

      #kld growth rate
      print('KLD growth rate for var ',my_var)
      print(np.nanmin( kld_c[my_var]),np.nanmax(kld_c[my_var])) 

      print('SPRD growth rate for var ')
      print(np.nanmin( sprd_c[my_var] ),np.nanmax( sprd_c[my_var] ))

      #Axes limits
      my_axes = [icoldelta*(icol-1)+hmargin+hoffset,irowdelta*(irow-1)+vmargin+voffset,irowdelta-2*hmargin,icoldelta-2*vmargin]

      ax = plt.axes( my_axes , facecolor=None )

      #The plot
      ybar = np.reshape( kld_c[my_var] , ( 1 , np.size( kld_c[my_var] ) ) )
      xbar = np.reshape( sprd_c[my_var] , ( 1 , np.size( kld_c[my_var] ) ) )
      zbar = np.reshape( max_dbz , ( 1 , np.size( max_dbz ) ) )

      mask = np.logical_and( np.logical_not( np.isnan( xbar) ) , np.logical_not( np.isnan( ybar) ) )
      mask = np.logical_and( mask , np.logical_not( np.isnan( zbar ) ) )
      mask = np.logical_and( mask , radar_mask )

      xbar=xbar[ mask ] 
      ybar=ybar[ mask ]
      zbar=zbar[ mask ]

      xbar = xbar / np.std( xbar ) 
      ybar = ybar / np.std( ybar )

      r=str( np.round( 1000*( np.corrcoef(xbar,ybar)[0,1] ) )/1000 )

      p=ax.plot( xbar , ybar , 'bo' ,linestyle='None' , label='All (r=' + r + ')')

      xbar=xbar[zbar > 30.0 ]
      ybar=ybar[zbar > 30.0 ]
      
      r=str( np.round( 1000*( np.corrcoef(xbar,ybar)[0,1] ) )/1000 )

      p=ax.plot( xbar , ybar , 'ro' , linestyle='None' , label='Z > 30.0 dBz (r=' + r + ')')

      ax.legend(loc='upper center')
      print( 'R for var ' + my_var + ' ' + str( np.corrcoef(xbar,ybar)[0,1] ) )

      #Grid lines
      ax.grid(linewidth=1.0, color='k',alpha=0.5, linestyle='--')


      #Grid lines
      if irow == 1   :
         plt.xlabel('Spread')

      if icol == 1   :
         plt.ylabel('KLD')

      plt.title(titles[ivar] + ' ' + my_var,fontsize = 15)

      icol = icol + 1
      if icol > ncols   :
         icol = 1
         irow = irow - 1

   #plt.show()
   plt.savefig( plotbasedir + '/' + figname + '_' + my_var + '_' + ctime.strftime("%Y%m%d%H%M%S") + '.png' , format='png' , dpi=300)
   plt.close()


   ctime = ctime + delta_plots


