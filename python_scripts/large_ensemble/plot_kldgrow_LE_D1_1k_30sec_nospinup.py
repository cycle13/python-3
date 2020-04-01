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

import common_plot_functions   as cpf
import common_mask_functions   as cmf
import common_smooth_functions as csf

basedir='/home/ra001011/a03471/data/output_data/'

expname = '/LE_D1_1km_30sec_nospinup/'

filetypes=['guesgp']   #guesgp , analgp 

max_dbz_levels=[-55,-50,10,30,60]   #Levels for vertically averaged condensate.

plot_pressure_levels=np.array([950,500,300]) #np.array([950,500,300])   #Which levels will be plotted

plot_variables      =['u','v','w','tk','qv','dbz'] #['u','v','w','tk','qv','dbz','slp']


#Define initial and end times using datetime module.
itime = dt.datetime(2013,7,13,5,0, 30)  #Initial time.
etime = dt.datetime(2013,7,13,5,14,30)  #End time.

#Define the delta.
delta=dt.timedelta(seconds=30)

#Compute the total number of times
ntimes=int( 1 + np.around((etime-itime).seconds / delta.seconds) ) #Total number of times.

#=========================================================
#  LOOP OVER FILE TYPES
#=========================================================

for my_file_type in filetypes     :

  ctl_file = basedir + expname + '/ctl/' + my_file_type + '.ctl'

  plotbasedir=basedir + expname + '/plots/' + my_file_type + '/'

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

  latlon_file = basedir + expname + '/latlon/latlon.grd'

  tmp=ctlr.read_data_records(latlon_file,ctl=ctl_dict,records=np.array([0,1]))
  lat=tmp[:,:,1]
  lon=tmp[:,:,0]

  #=========================================================
  #  DEFINE REGIONS
  #=========================================================

  lati=np.array([34.75,34.6])
  late=np.array([35.25,34.9])
  loni=np.array([135.5,135.4])
  lone=np.array([136.25,135.7])

  zi=np.array([0,0,0])
  ze=np.array([nlev,nlev,nlev])

  reg_name='REG_1','REG_2','TOTAL'

  kldgr=dict()

  ctime=itime 

  #Add the global domain as a region.
  lati=np.append(lati,lat[0,0])
  late=np.append(late,lat[nx-1,ny-1])
  loni=np.append(loni,lon[0,0])
  lone=np.append(lone,lon[nx-1,ny-1])


  #Compute xi,xe,yi,ye corresponding to each region.
  xi , yi = cmf.lat_lon_to_i_j(lon,lat,loni,lati)
  xe , ye = cmf.lat_lon_to_i_j(lon,lat,lone,late)

  nregs=int( xi.shape[0] )

  #=========================================================
  #  DEFINE VARIABLES
  #=========================================================

  kldgr=dict()

  ens_mean=dict()

  kldgr_regional_mean=dict()
  kldgr_regional_max=dict()
  kldgr_regional_min=dict()

  kldgr_time_mean=dict()
  kldgr_time_std=dict()

  max_dbz=np.zeros([nx,ny,nlev])

  time_mean_max_dbz=np.zeros([nx,ny,nlev])

  my_kldgr=np.zeros([nx,ny,nlev])

  my_kldgr_sq=np.zeros([nx,ny,nlev])


  for var in plot_variables        :

    if var in ctl_dict['var_list']  :
 
       kldgr_regional_mean[var]=np.zeros([ntimes,nregs])
       kldgr_regional_max[var] =np.zeros([ntimes,nregs])
       kldgr_regional_min[var] =np.zeros([ntimes,nregs])

       varpos=[i for i,x in enumerate(ctl_dict['var_list']) if x == var][0]
       varsize=int( ctl_dict['var_size'][varpos] )
       if varsize == 0   :
          varsize = 1
       kldgr_time_mean[var]=np.zeros([nx,ny,varsize])
       kldgr_time_std[var]=np.zeros([nx,ny,varsize])

  #=========================================================
  #  START TIME LOOP
  #=========================================================

  it=0
  while ( ctime <= etime ):

   print( ctime )

   print ( 'Reading the moments and computing growth rate ')

   #=========================================================
   #  READ THE DATA
   #=========================================================

   ptime = ctime - delta 

   #Analysis spread at time T
   my_file=basedir + expname + ptime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/' + '/kldistance.grd'

   kld_p=ctlr.read_data_grads(my_file,ctl_dict,masked=False)

   #Gues spread at time T+1
   my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/' + '/kldistance.grd'

   kld_c=ctlr.read_data_grads(my_file,ctl_dict,masked=False)

   #Compute remove undef values, compute kld growth 
   #and smooth the result.
   for my_var in kld_c   :
       if my_var in plot_variables    :
          kld_c[my_var][ kld_c[my_var] == undef ] = 0 #np.nan 
          kld_p[my_var][ kld_p[my_var] == undef ] = 0 #np.nan
          kld_c[my_var][ kld_c[my_var] >  100   ] = 0 #np.nan
          kld_p[my_var][ kld_p[my_var] >  100   ] = 0 #np.nan
          kldgr[my_var] = ( kld_c[my_var] - kld_p[my_var] ) / delta.total_seconds()

          #Smooth the growing rate.
          kldgr[my_var][:,:,:,0]=csf.gaussian_smooth_2d( kldgr[my_var] , 2.0 , 5.0 )

   #Read the ensemble mean to get the information from the storm location.
   my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/' + my_file_type + '/' + '/moment0001.grd'

   ens_mean=ctlr.read_data_grads(my_file,ctl_dict,masked=False)
   
   #Compute max_dbz (we will use this to identify areas associated with clouds and convection)
   tmp_max_dbz = np.squeeze( np.nanmax(ens_mean['dbz'],2) )
   
   for ilev in range(0,nlev)   :          #Create a fake 3D array for the max_dbz
                                          #This is because the plotting function expects a 3D array as input.
      max_dbz[:,:,ilev]=tmp_max_dbz


   #=======================================================================================
   #Plot G.R. 
   #=======================================================================================

   # for key in kldgr :

   for var in plot_variables        :

      if var in ctl_dict['var_list']  :

         #Plot moments.
         my_kldgr=kldgr[var]
         my_kldgr[ my_kldgr > 100 ] = np.nan #There are som inf values in the reflectivity field.
         my_kldgr[ my_kldgr == undef ] = np.nan

         print('G.R. for Var ',var,' ',(np.nanmin(my_kldgr)),np.nanmax(my_kldgr))

         date=ctime.strftime("%Y%m%d%H%M%S")
         cpf.set_default()  #Restore defaults
         my_map='RdBu_r'   #cpf.cmap_discretize('Blues',10)
         cpf.figconf['figpath']=plotbasedir
         cpf.figconf['figsize']=(12,10)
         cpf.figconf['titlefontsize']=20
         cpf.figconf['labelfontsize']=20
         cpf.figconf['pcolor']=True
         cpf.figconf['shadedmin']=-0.0002
         if var == 'w'     :
            cpf.figconf['shadedmax']=0.0002
         else              :
            cpf.figconf['shadedmax']=0.0002
         cpf.figconf['shadedcolormap']=my_map
         cpf.figconf['colorbar']=True
         cpf.figconf['colorbarfontsize']=15
         cpf.figconf['axessize']=[0.1,0.1,0.8,0.8]
         cpf.figconf['contour']=True
         cpf.figconf['contourlevels']=max_dbz_levels
         cpf.figconf['contourcolormap']='Reds'
         cpf.figconf['gridline']=True
         cpf.figconf['ylabel']='Lat'
         cpf.figconf['xlabel']='Lon'
         cpf.figconf['xtick']=[134.5,135,135.5,136,136.5,137]
         cpf.figconf['ytick']=[34,34.5,35,35.5]
         #cpf.figconf['show']=True
 
         if my_kldgr.shape[2] > 1 :

             for il in plotlevels   : 
                cpf.figconf['title']='KLD G.R. for variable ' + var + ' at level ' + str( int(ctl_dict['vlevels'][np.asscalar(il)])) + ' at ' + date
                cpf.figconf['figname']='/Figure_KLDGR_' + var + '_' + date + '_' + str(int(ctl_dict['vlevels'][np.asscalar(il)])) 

                cpf.plot_x_y_cartopy( lon , lat , my_kldgr[:,:,il] , max_dbz[:,:,il] , my_kldgr[:,:,il] )

         else                    :
                cpf.figconf['title']='KLD G.R. for variable ' + var + ' at ' + date 
                cpf.figconf['figname']='/Figure_KLDGR_' + var + '_' + date    
     
                cpf.plot_x_y_cartopy( lon , lat , my_kldgr[:,:] , max_dbz[:,:,0] , my_kldgr[:,:] )

      #========================================================================================
      #Generate regional averages of the moment.
      #========================================================================================

      kldgr_regional_mean[var][it,:],kldgr_regional_max[var][it,:],kldgr_regional_min[var][it,:]=cmf.get_regional_average_grid(my_kldgr,xi,xe,yi,ye,zi,ze,undef)

      #=========================================================================================
      #Acumulate the moment statistics in time. 
      #=========================================================================================

      kldgr_time_mean[var]=kldgr_time_mean[var] +  my_kldgr[:,:,:,0]                     #Accumulate the mean.
      kldgr_time_std[var]=kldgr_time_std[var]   + np.power(  my_kldgr[:,:,:,0]  , 2 )    #Accumulate the standard deviation.

   time_mean_max_dbz=time_mean_max_dbz + max_dbz  #To accumulate integrated liquid.

   #=========================================================================================
   #Advance time
   #=========================================================================================

   ctime = ctime + delta
 
   it = it + 1

  print ( "Finish time loop" )

  ntimes=it


  #=========================================================================================
  #Generate some time independent plots.
  #=========================================================================================


  #=========================================================================================
  #Plot the mean GR and its standard deviation.
  #=========================================================================================

  #for key in kldgr_time_mean  :
  for var in plot_variables        :

   if var in ctl_dict['var_list']  :

      print("plotting the kld gr mean")

      #Plot time mean of the moments.

      cpf.set_default()  #Restore defaults
      my_map='RdBu_r'
      cpf.figconf['figpath']=plotbasedir + '/time_independent_plots/'
      cpf.figconf['figsize']=(12,10)
      cpf.figconf['titlefontsize']=20
      cpf.figconf['labelfontsize']=20
      cpf.figconf['pcolor']=True
      cpf.figconf['shadedmin']=-0.00008
      if var == 'w'     :
         cpf.figconf['shadedmax']=0.0003
         cpf.figconf['shadedmin']=-0.0003
      else              :
         cpf.figconf['shadedmax']=0.00008
      cpf.figconf['shadedcolormap']=my_map
      cpf.figconf['colorbar']=True
      cpf.figconf['colorbarfontsize']=15
      cpf.figconf['axessize']=[0.1,0.1,0.8,0.8]
      cpf.figconf['contour']=True
      cpf.figconf['contourlevels']=np.arange(0,1,0.2)
      cpf.figconf['contourcolormap']='Reds'
      cpf.figconf['gridline']=True
      cpf.figconf['xtick']=[134.5,135,135.5,136,136.5,137]
      cpf.figconf['ytick']=[34,34.5,35,35.5]
      cpf.figconf['ylabel']='Lat'
      cpf.figconf['xlabel']='Lon'
      cpf.figconf['show']='True'

      time_mean_max_dbz = time_mean_max_dbz / ntimes

      if my_kldgr.shape[2] > 1 :

          for il in plotlevels   :

             my_kldgr_mean=np.squeeze(kldgr_time_mean[var][:,:,il]) / ntimes
             my_kldgr_sprd=np.sqrt( np.squeeze(kldgr_time_std[var][:,:,il]) / ntimes  - np.power( my_kldgr_mean,2) / ntimes )

             cpf.figconf['title']='Mean KLD G.R. for variable ' + var + ' at level ' + str( ctl_dict['vlevels'][np.asscalar(il)])
             cpf.figconf['figname']='/Figure_KLDGR_mean_sprd_' + var + '_' + date + '_' + str(ctl_dict['vlevels'][il])
             cpf.plot_x_y_cartopy( lon , lat , my_kldgr_mean , my_kldgr_sprd , my_kldgr )

      else                    :
          cpf.figconf['title']='Mean KLD G.R. for variable ' + var 
          cpf.figconf['figname']='/Figure_KLDGR_mean_sprd_' + var + '_' + date

          cpf.plot_x_y_cartopy( lon , lat , my_kldgr_mean , my_kldgr_sprd , my_kldgr ) 

  #=========================================================================================
  #Plot time evolution of G.R. over different regions.
  #=========================================================================================

  #for key in kldgr_time_mean  :
  for var in plot_variables        :

   if var in ctl_dict['var_list']  :


      #Plot time series of moments for different regions.

      for ireg in range(0,nregs) :

         iregstr="%04d" % ( ireg + 1 )
         cpf.set_default()  #Restore defaults
         cpf.figconf['figpath']=plotbasedir + '/time_independent_plots/'
         cpf.figconf['figname']='Figure_kldgr_reg_' + iregstr + '_var_' + var 
         cpf.figconf['figsize']=(12,10)
         cpf.figconf['title']='Mean KLD G.R. for var ' + var + ' for region ' + reg_name[ireg] 
         cpf.figconf['titlefontsize']=20
         cpf.figconf['labelfontsize']=20
         cpf.figconf['axessize']=[0.1,0.1,0.8,0.8]
         cpf.figconf['gridline']=True

         cpf.figconf['linestyle']=['-']
         cpf.figconf['linemarker']=['o']
         cpf.figconf['linecolor']=['b']
         cpf.figconf['ylabel']=['Growing Rate']
         cpf.figconf['xlabel']=['Time (seconds)']
         #cpf.figconf['show']=True

         time = np.arange( 0 , np.shape(kldgr_regional_mean[var])[0] ) * delta.total_seconds()

         cpf.plot_lines( time , kldgr_regional_mean[var][:,ireg] )



