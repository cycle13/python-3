# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016

@author:
"""
import sys
sys.path.append('../../common_python/common_functions/')
sys.path.append('../../common_python/common_modules/')

import numpy as np
import datetime as dt
import ctl_reader as ctlr
import os
import matplotlib.pyplot as plt
from common_functions import common_functions as cf


import common_plot_functions as cpf
import common_mask_functions as cmf

basedir='/home/ra001011/a03471/data/output_data/'


exps=['LE_D1_1km_5min','LE_D1_1km_2min','LE_D1_1km_1min','LE_D1_1km_30sec','LE_D1_1km_30sec_nospinup','LE_D1_1km_1min_4D']

deltat=[300,120,60,30,30,60]

filetype='guesgp'   #analgp , analgz , guesgp , guesgz


lat_radar=34.823
lon_radar=135.523
radar_range=50.0e3   #Radar range in meters (to define the radar mask)

init_times = ['20130713050500','20130713050400','20130713050500','20130713050500','20130713050500','20130713050500'] 

sigma_smooth=2.0

#Define initial and end times using datetime module.
#itime = dt.datetime(2013,7,13,5, 5,0)  #Initial time.
etime = dt.datetime(2013,7,13,5,55,0)  #End time.

variable_combination=[['dbz','dbz','dbz','dbz','dbz','v' ,'v' ,'v','v' ],
                      ['tk' ,'v'  ,'qv' ,'w'  ,'dbz','tk','qv','w','v' ]]

obs_increment = [5.0 , 5.0 , 5.0 , 5.0 , 5.0 , 2.0 , 2.0 , 2.0 , 2.0 ]
obs_error     = [5.0 , 5.0 , 5.0 , 5.0 , 5.0 , 2.0 , 2.0 , 2.0 , 2.0 ]


nbv=1000


#=========================================================
#  LOOP OVER FILE TYPES
#=========================================================

for iv in range( len( variable_combination[0] ) ) :   #Loop over selected variable combinations.

   var_obs = variable_combination[0][iv]
   var_upd = variable_combination[1][iv]

   rmse=dict()

   ensemble_mean =dict()
   
   for iexp , my_exp in enumerate(exps)   :

      itime = dt.datetime.strptime(init_times[iexp] , '%Y%m%d%H%M%S' )

      delta=dt.timedelta(seconds=deltat[iexp])

      #Compute the total number of times
      ntimes=int( 1 + np.around((etime-itime).seconds / delta.seconds) ) #Total number of times.

      #=========================================================
      #  READ CTL FILE
      #=========================================================

      if iexp == 0 and iv == 0  :

         ctl_file = basedir + '/' + my_exp + '/ctl/update_mean_diff.ctl'
 
         ctl_dict = ctlr.read_ctl( ctl_file )

         nx=ctl_dict['nx']
         ny=ctl_dict['nx']
         nlev=ctl_dict['nz']
         nt=int(1)             #Force the number of times to be one.
         ctl_dict['nt']=int(1) #Force the number of times to be one.

         undef=np.float32( ctl_dict['undef'] )

         ctl_file_2 = basedir + '/' + my_exp + '/ctl/guesgp.ctl'
         
         ctl_dict_2 = ctlr.read_ctl( ctl_file_2 )

      #=========================================================
      #  READ LAT LON
      #=========================================================
 
      if iexp == 0 and iv == 0 :

         latlon_file = basedir + '/' + my_exp + '/latlon/latlon.grd'

         tmp=ctlr.read_data_records(latlon_file,ctl=ctl_dict,records=np.array([0,1]))
         lat=tmp[:,:,1]
         lon=tmp[:,:,0]

         #Exclude areas outside the radar domain.
         radar_mask = cmf.distance_range_mask( lon_radar , lat_radar , radar_range , lon , lat )

      #=========================================================
      #  READ THE DATA
      #=========================================================

      my_file=basedir + '/' + my_exp + '/time_mean/'+ filetype + '/update_comp_meandiff_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '_rmse' + '.grd'

      rmse[my_exp]= np.nanmean( np.delete( ctlr.read_data(my_file,ctl_dict,undef2nan=True) , 4 , 2 ) , 2 )

      rmse[my_exp][ np.isnan( rmse[my_exp] ) ] = undef
      rmse[my_exp]=cf.gaussian_filter(field0=rmse[my_exp],dx=1.0,sigma=sigma_smooth,nx=nx,ny=ny,undef=undef)
      rmse[my_exp][ rmse[my_exp] == undef ] = np.nan


      print('RMSE values')
      print( my_exp )
      print( np.min( rmse[my_exp] ) , np.max( rmse[my_exp] ) )

      my_file=basedir + '/' + my_exp + '/time_mean/'+ filetype + '/' + '/moment0001_mean.grd'

      ensemble_mean[my_exp]= np.delete( ctlr.read_data_grads(my_file,ctl_dict_2)['dbz'] , 4 ,2 )
      ensemble_mean[my_exp]=np.max( ensemble_mean[my_exp] , 2 ) 
      ensemble_mean[my_exp]=cf.gaussian_filter(field0=ensemble_mean[my_exp],dx=1.0,sigma=sigma_smooth,nx=nx,ny=ny,undef=undef)
      ensemble_mean[my_exp][ ensemble_mean[my_exp] == undef ] = np.nan
 
 
      print('State variable values')
      print( my_exp )
      print( np.min( ensemble_mean[my_exp]) , np.max( ensemble_mean[my_exp]) )

 
   print(' Finish the loop over experiments ' )    

   #=========================================================================================
   #Plot the mean KLD and its standard deviation.
   #=========================================================================================

   nexp = len( exps )

   for iexp,my_exp in enumerate(exps)  :

      if iexp > 0  :
 
         rmse[my_exp] = 100 * ( ( rmse[my_exp] - rmse[exps[0]] ) / rmse[exps[0]] )

      #Apply radar mask
      print('Var Obs ',var_obs,' Var Upd ',var_upd)
      tmp = np.copy(rmse[my_exp])
      tmp[ np.logical_not( radar_mask ) ] = np.nan 
      rmse[my_exp] = tmp[:]
      print('Parameter Min=',np.nanmin( tmp ),' Parameter Max=',np.nanmax( tmp ) )
      tmp = np.copy( ensemble_mean[my_exp] )
      tmp[ np.logical_not( radar_mask ) ] = np.nan 
      ensemble_mean[my_exp] = tmp[:]
      print('Ref Min=',np.nanmin( tmp ),' Ref Max=',np.nanmax( tmp ) )

   #Plot time mean of the moments.

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
 
      varsh = rmse[my_exp]
      varc  = ensemble_mean[my_exp]

      if iexp == 0 :
         my_map =cpf.cmap_discretize('Blues',10)
         smin = 0.0
         smax = np.nanmax( varsh )
      else         :
         my_map =cpf.cmap_discretize('seismic',11)
         smin = -100.0
         smax =  100.0

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
   plt.savefig('Figure_update_rmse_VertMean_Exps' + '_vobs_' + var_obs + '_upd_obs_' + var_upd + '_ens_size_' + str(nbv) + '_obs_inc_' + str(obs_increment[iv]) + '.png'  , format='png' , dpi=300)
   plt.close()


