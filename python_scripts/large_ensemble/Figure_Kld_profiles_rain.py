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
import matplotlib.pyplot as plt
from common_functions import common_functions as cf


import common_plot_functions as cpf
import common_mask_functions as cmf

basedir='/home/ra001011/a03471/data/output_data/'

figname='./Figure_Kld_profile_rain_Exps'

filename='kldistance_mean.grd'

exps=['LE_D1_1km_5min','LE_D1_1km_2min','LE_D1_1km_1min','LE_D1_1km_30sec','LE_D1_1km_30sec_nospinup','LE_D1_1km_1min_4D']

deltat=[300,120,60,30,30,60]

filetype='guesgp'   #analgp , analgz , guesgp , guesgz

plot_variables=['tk','u','v','w','dbz','qv']

lat_radar=34.823
lon_radar=135.523
radar_range=50.0e3   #Radar range in meters (to define the radar mask)

init_times = ['20130713050500','20130713050400','20130713050500','20130713050500','20130713050500','20130713050500'] 

sigma_smooth=2.0

#Define initial and end times using datetime module.
#itime = dt.datetime(2013,7,13,5, 5,0)  #Initial time.
etime = dt.datetime(2013,7,13,5,55,0)  #End time.

#=========================================================
#  LOOP OVER FILE TYPES
#=========================================================

profile = dict()   

profile_w = dict()

for iexp , my_exp in enumerate(exps)   :

  itime = dt.datetime.strptime(init_times[iexp] , '%Y%m%d%H%M%S' )

  delta=dt.timedelta(seconds=deltat[iexp])

  #Compute the total number of times
  ntimes=int( 1 + np.around((etime-itime).seconds / delta.seconds) ) #Total number of times.

  ctl_file = basedir + '/' + my_exp + '/ctl/' + filetype + '.ctl'

  #=========================================================
  #  READ CTL FILE (To get vertical levels)
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

  profile[ my_exp ] = dict()

  for my_var in plot_variables  :

      my_file=basedir + '/' + my_exp + '/time_mean/'+ filetype + '/kldistance_' + my_var + '_rain_profile.npz' 

      tmp = np.load(my_file)

      #Remove the fourth element due to post-processing error.
      profile[my_exp][my_var]=np.delete( tmp['parameter'] , 4 , axis=0 )


  my_file=basedir + '/' + my_exp + '/time_mean/'+ filetype + '/moment0001_w_rain_profile.npz'

  tmp = np.load(my_file)

  #Remove the fourth element due to post-processing error.
  profile_w[my_exp] = np.delete( tmp['parameter'] , 4 , axis=0 )
 
print(' Finish the loop over experiments ' )    

#=========================================================================================
#Plot the mean KLD and its standard deviation.
#=========================================================================================

nexp = len( exps )

levels=ctl_dict['vlevels']
levels_str=list()

#To solve the postprocessing error.
levels=np.delete(levels,4,axis=0)
levels[3]=850.0


import matplotlib
#import matplotlib.cm     as cm
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
#from mpl_toolkits.basemap import Basemap
#import cartopy.crs as ccrs



#for key in kld_time_mean  :
for my_var in plot_variables :

   #Plot time mean of the moments.

   plt.figure(1,figsize=[11,5])
   #plt.figure(1)

   #Start subplots
   ncols = 3
   nrows = 2
   icoldelta = 1.0/ncols
   irowdelta = 1.0/nrows
   hmargin=0.15
   vmargin=0.01
   hoffset=-0.1
   voffset=0.1

   icol = 1 
   irow = nrows

   titles = ['(a)','(b)','(c)','(d)','(e)','(f)']

   my_map=cpf.cmap_discretize('YlGnBu',20)

   for iexp,my_exp in enumerate(exps)  : 

      time_index=np.arange( 0 , profile[my_exp][my_var].shape[1] ) * 5.0 + 5.0 #Time index in minutes.

      xtick= time_index[1::2]

      if my_var == 'w'     :
        smax=0.1
      elif my_var == 'qv'  :
        smax=0.1
      elif my_var == 'dbz' :
        smax=0.4
      else              :
        smax=0.06

      smin=0.0

      #Axes limits
      my_axes = [icoldelta*(icol-1)+hmargin+hoffset,irowdelta*(irow-1)+vmargin+voffset,irowdelta-2*hmargin,icoldelta-2*vmargin]
 
      ax = plt.axes( my_axes , facecolor=None )

      #The pcolor
      print(np.shape(time_index),np.shape(levels),np.shape( np.transpose( profile[my_exp][my_var] ) ) )
      p=ax.pcolor( time_index , -np.log( levels ) ,  np.squeeze( profile[my_exp][my_var] ) 
           ,vmin=smin , vmax=smax ,cmap=my_map )

      #ax.set_extent( axesrange  ccrs.PlateCarree())

      #Grid lines
      ax.grid(linewidth=1.0, color='k',alpha=0.5, linestyle='--')

      #Grid and ticks
      ax.set_xticks(xtick)
      my_levs = levels[[0,3,4,6,8]]
      ytick=-np.log(my_levs)
      ax.set_yticks(ytick)
      #Get the level string list.
      levels_str=[]
      for ilev in my_levs  :
         levels_str.append( str(int(ilev)) )
      ax.set_yticklabels(levels_str,fontsize=10,color='k')

      #Colorbar
      cb=plt.colorbar(p)
      cb.ax.tick_params(labelsize=12)

      #Contour map
      c=ax.contour( time_index[:-1] , -np.log(levels) , np.squeeze( profile_w[my_exp][:,:-1] ) ,
        levels = [0.5,1.0,1.5] ,colors='k')
 
      #Axes label and title
      if irow == 1   :
         plt.xlabel('Time (minutes)')
      if icol == 1   :
         plt.ylabel('Pressure (hPa) ')

      ax.spines['right'].set_visible(False)
      ax.spines['top'].set_visible(False)

      plt.title(titles[iexp],fontsize = 15)

      icol = icol + 1
      if icol > ncols   :
         icol = 1
         irow = irow - 1
 
   #plt.show()
   plt.savefig( figname + '_' + my_var + '.png' , format='png' , dpi=300)
   plt.close()


