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

import common_plot_functions   as cpf
import common_mask_functions   as cmf
import common_smooth_functions as csf

basedir='/home/ra001011/a03471/data/output_data/'

expname = '/LE_D1_1km_30sec_nospinup/'

plot_variables = ['u','w','tk','qv','dbz'] #['u','v','w','tk','qv','dbz','slp']

plot_colors    = ['b','r','k','y','m']

#Define initial and end times using datetime module.
itime = dt.datetime(2013,7,13,5,0, 0)  #Initial time.
etime = dt.datetime(2013,7,13,5,59,30)  #End time.

#Define the delta.
delta=dt.timedelta(seconds=30)

#Compute the total number of times
ntimes=int( 1 + np.around((etime-itime).seconds / delta.seconds) ) #Total number of times.

#=========================================================
#  LOOP OVER FILE TYPES
#=========================================================

ctl_file = basedir + expname + '/ctl/guesgp.ctl'

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
#  DEFINE REGIONS
#=========================================================

lati=np.array([34.8,34.4])
late=np.array([35.15,34.6])
loni=np.array([135.5,135.4])
lone=np.array([136.0,135.8])

zi=np.array([5,5])
ze=np.array([9,9])

reg_names=['HI-Z','LOW-Z']

kldgr=dict()

ctime=itime 

#Compute xi,xe,yi,ye corresponding to each region.
xi , yi = cmf.lat_lon_to_i_j(lon,lat,loni,lati)
xe , ye = cmf.lat_lon_to_i_j(lon,lat,lone,late)

nregs=int( xi.shape[0] )

#=========================================================
#  DEFINE VARIABLES
#=========================================================

kldgr=dict()

ens_mean=dict()

kld_mean=dict()
bim_mean=dict()
m01_mean=dict()
m02_mean=dict()

for my_reg in reg_names  :
   kld_mean[my_reg] = dict()
   bim_mean[my_reg] = dict()
   m01_mean[my_reg] = dict()
   m02_mean[my_reg] = dict()
     
   for var in plot_variables        :
      kld_mean[my_reg][var] =np.zeros([ntimes])
      bim_mean[my_reg][var] =np.zeros([ntimes])
      m01_mean[my_reg][var] =np.zeros([ntimes])
      m02_mean[my_reg][var] =np.zeros([ntimes]) 


#=========================================================
#  START TIME LOOP
#=========================================================

it=0

time_index=np.zeros([ntimes])

while ( ctime <= etime ):

   print( ctime )

   print ( 'Reading the moments and computing growth rate ')

   #=========================================================
   #  READ THE DATA
   #=========================================================

   #Read the first guess
   my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/kldistance.grd'
   kld=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=False)

   my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/bimodality_index.grd'
   bim=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=False)

   my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/moment0001.grd'
   m01=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=False)

   my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/guesgp/moment0002.grd'
   m02=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=False)


   print( kld.keys() )

   for ireg,my_reg in enumerate( reg_names )  : 
      for ivar,my_var in enumerate( plot_variables ) :
         kld_mean[my_reg][my_var][it] , _ , _ = cmf.get_regional_average_grid(kld[my_var],xi[ireg],xe[ireg],yi[ireg],ye[ireg],zi[ireg],ze[ireg],undef)
         bim_mean[my_reg][my_var][it] , _ , _ = cmf.get_regional_average_grid(bim[my_var],xi[ireg],xe[ireg],yi[ireg],ye[ireg],zi[ireg],ze[ireg],undef)
         m01_mean[my_reg][my_var][it] , _ , _ = cmf.get_regional_average_grid(m01[my_var],xi[ireg],xe[ireg],yi[ireg],ye[ireg],zi[ireg],ze[ireg],undef)
         m02_mean[my_reg][my_var][it] , _ , _ = cmf.get_regional_average_grid(m02[my_var],xi[ireg],xe[ireg],yi[ireg],ye[ireg],zi[ireg],ze[ireg],undef)


   #time_index[it] = np.floor( it/2 )

   time_index[it] = it * 30.0 / 60.0

   it = it + 1


   ##Read the analysis
   #my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/analgp/kldistance.grd'
   #kld=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=False)

   #my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/analgp/bimodality_index.grd'
   #bim=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=False)

   #my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/analgp/moment0001.grd'
   #m01=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=False)

   #my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/analgp/moment0002.grd'
   #m02=ctlr.read_data_grads(my_file,ctl_dict,masked=False,undef2nan=False)


   #for ireg,my_reg in enumerate( reg_names )  :
   #   for ivar,my_var in enumerate( plot_variables ) :
   #      kld_mean[my_reg][my_var][it] , _ , _ = cmf.get_regional_average_grid(kld[my_var],xi[ireg],xe[ireg],yi[ireg],ye[ireg],zi[ireg],ze[ireg],undef)
   #      bim_mean[my_reg][my_var][it] , _ , _ = cmf.get_regional_average_grid(bim[my_var],xi[ireg],xe[ireg],yi[ireg],ye[ireg],zi[ireg],ze[ireg],undef)
   #      m01_mean[my_reg][my_var][it] , _ , _ = cmf.get_regional_average_grid(m01[my_var],xi[ireg],xe[ireg],yi[ireg],ye[ireg],zi[ireg],ze[ireg],undef)
   #      m02_mean[my_reg][my_var][it] , _ , _ = cmf.get_regional_average_grid(m02[my_var],xi[ireg],xe[ireg],yi[ireg],ye[ireg],zi[ireg],ze[ireg],undef)


   ctime = ctime + delta

   #time_index[it] = np.floor( it/2 )

   #it = it + 1

print ( "Finish time loop" )

ntimes=it


for my_reg in reg_names :
   for my_var in plot_variables :

       kld_mean[my_reg][my_var][ np.logical_or( np.isinf( kld_mean[my_reg][my_var] ) , kld_mean[my_reg][my_var] == 0 )] = np.nan
       bim_mean[my_reg][my_var][ np.logical_or( np.isinf( bim_mean[my_reg][my_var] ) , bim_mean[my_reg][my_var] == 0 )] = np.nan
       m01_mean[my_reg][my_var][ np.logical_or( np.isinf( m01_mean[my_reg][my_var] ) , m01_mean[my_reg][my_var] == 0 )] = np.nan
       m02_mean[my_reg][my_var][ np.logical_or( np.isinf( m01_mean[my_reg][my_var] ) , m01_mean[my_reg][my_var] == 0 )] = np.nan

#=========================================================================================
#Plot time evolution of G.R. over different regions.
#=========================================================================================

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

#Start subplots
ncols = 2
nrows = 2

titles = ['(a)','(b)','(c)','(d)','(e)','(f)']
labels = ['$U$','$W$','$T$','$q_v$','$Z$']

fig, axs = plt.subplots( nrows,ncols ,  figsize=[6,5] )
fig.subplots_adjust(wspace=0.2,hspace=0.0,bottom=0.05,left=0.05,right=0.97,top=0.93)
icol=0
irow=0
for ireg,my_reg in enumerate( reg_names ) :
   ax = axs[irow,icol]
   for ivar,my_var in enumerate( plot_variables ) :
       ax.plot( time_index , ( kld_mean[my_reg][my_var] - np.nanmean( kld_mean[my_reg][my_var] ) )/ np.nanstd( kld_mean[my_reg][my_var] ) , linewidth = 2.0 , color=plot_colors[ivar],label=labels[ivar])
       ax.set_ylim([-5,5])
   ax.grid(linewidth=0.5, color='k',alpha=0.5, linestyle='--')
   ax.set_title(titles[icol+ncols*irow])
   if icol == 0 :
      ax.legend(fontsize=8,loc='lower right')
   icol = icol + 1
icol=0
irow=1
#for ireg,my_reg in enumerate( reg_names ) :
#   ax = axs[irow,icol]
#   for ivar,my_var in enumerate( plot_variables ) :
#       ax.plot( time_index , ( bim_mean[my_reg][my_var] - np.nanmean( bim_mean[my_reg][my_var] ) )/ np.nanstd( bim_mean[my_reg][my_var] ) , linewidth = 2.0 , color=plot_colors[ivar])
#       ax.set_ylim([-5,5])
#   icol = icol + 1
#   ax.grid(linewidth=0.5, color='k',alpha=0.5, linestyle='--')
#icol=0
#irow=2
for ireg,my_reg in enumerate( reg_names ) :
   ax = axs[irow,icol]
   for ivar,my_var in enumerate( plot_variables ) :
       ax.plot( time_index , ( m02_mean[my_reg][my_var] - np.nanmean( m02_mean[my_reg][my_var] ) )/ np.nanstd( m02_mean[my_reg][my_var] ) , linewidth = 2.0 , color=plot_colors[ivar])
       ax.set_ylim([-5,5])
   icol = icol + 1
   ax.grid(linewidth=0.5, color='k',alpha=0.5, linestyle='--')
   ax.set_title(titles[icol+ncols*irow])
#plt.show()
plt.savefig( 'Figure_spinup_timeseries.eps' , format='eps' , dpi=300)


