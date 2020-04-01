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

figname='./Figure_bimodality_Frequency_Exps'

filename='bimodality_index_freq.grd'

exps=['LE_D1_1km_5min','LE_D1_1km_2min','LE_D1_1km_1min','LE_D1_1km_30sec','LE_D1_1km_30sec_nospinup','LE_D1_1km_1min_4D']

deltat=[300,120,60,30,30,60]

filetype='guesgp'   #analgp , analgz , guesgp , guesgz

plot_pressure_levels=np.array([600])   #Which levels will be plotted

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

parameter = dict()   

ensemble_freq = dict()

for iexp , my_exp in enumerate(exps)   :

  itime = dt.datetime.strptime(init_times[iexp] , '%Y%m%d%H%M%S' )

  delta=dt.timedelta(seconds=deltat[iexp])

  #Compute the total number of times
  ntimes=int( 1 + np.around((etime-itime).seconds / delta.seconds) ) #Total number of times.

  ctl_file = basedir + '/' + my_exp + '/ctl/' + filetype + '.ctl'

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

     plotlevels=list()
     for ilev in plot_pressure_levels  :
        tmp=np.nonzero( ilev == ctl_dict['vlevels'])  
        if np.size( tmp ) > 0  :
           plotlevels.append(tmp[0])
        else                   :
           print('Warning! Level ' + str(ilev) + ' not found in this dataset')
   
     plotlevels=np.array(plotlevels)  #From list to numpy array.

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

  my_file=basedir + '/' + my_exp + '/time_mean/'+ filetype + '/' + '/' + filename

  parameter[my_exp]=ctlr.read_data_grads(my_file,ctl_dict,masked=False)

  print('KLD values')
  print( my_exp )
  for var in parameter[my_exp] :
     print( var )
     print( np.min( parameter[my_exp][var] ) , np.max( parameter[my_exp][var] ) )

  my_file=basedir + '/' + my_exp + '/time_mean/'+ filetype + '/' + '/moment0001_mean.grd'

  ensemble_freq[my_exp]=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
  
  print('State variable values')
  print( my_exp )
  for var in ensemble_freq[my_exp] :
     print( var )
     print( np.min( ensemble_freq[my_exp][var] ) , np.max( ensemble_freq[my_exp][var] ) )

 
 
print(' Finish the loop over experiments ' )    

#=========================================================================================
#Plot the freq KLD and its standard deviation.
#=========================================================================================

nexp = len( exps )

plot_bim_freq = np.zeros( [ nx , ny , nexp ] )
plot_dbz_mean = np.zeros( [ nx , ny , nexp ] )

#for key in bim_time_freq  :
for var in plot_variables :

   for iexp,my_exp in enumerate(exps)  :

      #Prepare the data for plotting. 

      #Smooth output

      parameter[my_exp][var][ parameter[my_exp][var] == undef  ] = np.nan

      plot_bim_freq[:,:,iexp]  = np.squeeze( np.nanmean( parameter[my_exp][var] , 2) )
      plot_dbz_mean[:,:,iexp]  = np.squeeze( np.nanmean( ensemble_freq[my_exp]['dbz'] , 2) )

      print('Smoothing parameter for var=',var)
      plot_bim_freq[:,:,iexp]=cf.gaussian_filter(field0=plot_bim_freq[:,:,iexp],dx=1.0,sigma=sigma_smooth,nx=nx,ny=ny,undef=undef)

      if iexp == 0 :
         plot_bim_freq[:,:,iexp][ plot_bim_freq[:,:,iexp]  == 0 ] = np.nan

      if iexp > 0  :
         plot_bim_freq[:,:,iexp] = 100 * ( ( plot_bim_freq[:,:,iexp] - plot_bim_freq[:,:,0] ) / plot_bim_freq[:,:,0] )

      tmp = np.copy(plot_bim_freq[:,:,iexp])
      tmp[ np.logical_not( radar_mask ) ] = np.nan 
      plot_bim_freq[:,:,iexp] = tmp[:]
      print('Parameter Min=',np.nanmin( tmp ),' Parameter Max=',np.nanmax( tmp ) )
      tmp = np.copy( plot_dbz_mean[:,:,iexp] )
      tmp[ np.logical_not( radar_mask ) ] = np.nan 
      plot_dbz_mean[:,:,iexp] = tmp[:]
      print('Ref Min=',np.nanmin( tmp ),' Ref Max=',np.nanmax( tmp ) )

   #Plot time freq of the moments.

   plt.figure(1,figsize=[10,5])
   #plt.figure(1)

   #Start subplots
   ncols = 3
   nrows = 2
   icoldelta = 1.0/ncols
   irowdelta = 1.0/nrows
   hmargin=0.0
   vmargin=0.0
   hoffset=-0.15
   voffset=0.1

   icol = 1 
   irow = nrows

   xtick=[134.5,135,135.5,136,136.5,137]
   ytick=[34,34.5,35,35.5] 
   axesrange=[134.97,136.09,34.36,35.30]
   titles = ['(a)','(b)','(c)','(d)','(e)','(f)']

   import matplotlib
   #import matplotlib.cm     as cm
   import matplotlib.pyplot as plt
   import matplotlib.ticker as mticker
   #from mpl_toolkits.basemap import Basemap
   import cartopy.crs as ccrs



   for iexp,my_exp in enumerate(exps)  : 
 
      varsh = plot_bim_freq[:,:,iexp]
      varc  = plot_dbz_mean[:,:,iexp]

      if iexp == 0 :
         my_map =cpf.cmap_discretize('Blues',10)
         smin = 0.0
         smax = 0.1
      else         :
         my_map =cpf.cmap_discretize('jet',11)
         smin = -50.0
         smax =  50.0

      #Axes limits
      my_axes = [icoldelta*(icol-1)+hmargin+hoffset,irowdelta*(irow-1)+vmargin+voffset,irowdelta-2*hmargin,icoldelta-2*vmargin]
 
      ax = plt.axes( my_axes , facecolor=None , projection=ccrs.PlateCarree() )

      #The pcolor
      p=ax.pcolor(lon , lat ,  np.transpose( np.squeeze( varsh ) ) ,
        transform=ccrs.PlateCarree(),vmin=smin , vmax=smax ,cmap=my_map )

      ax.set_extent( axesrange , ccrs.PlateCarree())

      #Grid lines
      gl=ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                linewidth=1.0, color='k',alpha=0.5, linestyle='--')

      #Grid and ticks
      gl.xlocator=mticker.FixedLocator(xtick)
      gl.ylocator=mticker.FixedLocator(ytick)

      #Coastline plot
      ax.coastlines('10m',linewidth=1.0)

      #Colorbar
      cb=plt.colorbar(p)
      cb.ax.tick_params(labelsize=12)

      #Contour map
      c=ax.contour( lon , lat , np.transpose( np.squeeze( varc ) ) ,
        levels = [1,20] ,colors='k')
 
      #Axes label and title
      if irow == 1   :
         plt.xlabel('Latitude')
      if icol == 1   :
         plt.ylabel('Longitude')

      plt.title(titles[iexp],fontsize = 15)

      gl.xlabel_style = {'size': 12, 'color': 'k' }
      gl.ylabel_style = {'size': 12, 'color': 'k' }
      gl.xlabels_top = False
      gl.ylabels_right = False

      icol = icol + 1
      if icol > ncols   :
         icol = 1
         irow = irow - 1
 
   #plt.show()
   plt.savefig( figname + '_' + var + '.png' , format='png' , dpi=300)
   plt.close()


