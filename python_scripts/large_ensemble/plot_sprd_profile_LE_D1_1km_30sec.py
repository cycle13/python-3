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

import common_plot_functions as cpf
import common_mask_functions as cmf


basedir='/home/ra001011/a03471/data/output_data/'

expname = '/LE_D1_1km_30sec/'

plotbasedir=basedir + expname + '/plots/'

buffer_zone_size=20  #Data close to the domain borders will be ignored.

tr_rain=30          #Ref above this tr means rainy grid point.
tr_norain=0.0       #Ref below this tr means no rain

filetypes=['guesgp']   #analgp , analgz , guesgp , guesgz

plot_variables=['u','v','w','tk','qv','dbz']   


#Define initial and end times using datetime module.
itime = dt.datetime(2013,7,13,5, 0,30)  #Initial time.
etime = dt.datetime(2013,7,13,5,20,0)  #End time.

lat_radar=34.823
lon_radar=135.523
radar_range=50.0e3   #Radar range in meters (to define the radar mask)

#Define the delta.
delta=dt.timedelta(seconds=30)

#Compute the total number of times
ntimes=int( 1 + np.around((etime-itime).seconds / delta.seconds) ) #Total number of times.
times=np.zeros(ntimes)

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

  #=========================================================
  #  READ LAT LON
  #=========================================================

  latlon_file = basedir + expname + '/latlon/latlon.grd'

  tmp=ctlr.read_data_records(latlon_file,ctl=ctl_dict,records=np.array([0,1]))
  lat=np.squeeze(tmp[:,:,1])
  lon=np.squeeze(tmp[:,:,0])

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

  sprd=dict()

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

  sprd=dict()

  ens_mean=dict()

  sprd_regional_mean=dict()
  sprd_regional_max=dict()
  sprd_regional_min=dict()

  sprd_time_mean=dict()
  sprd_time_std=dict()

  max_dbz=np.zeros([nx,ny,nlev])

  time_mean_max_dbz=np.zeros([nx,ny,nlev])

  my_sprd=np.zeros([nx,ny,nlev])

  my_sprd_sq=np.zeros([nx,ny,nlev])


  for var in plot_variables        :

    if var in ctl_dict['var_list']  :

       sprd_regional_mean[var]=np.zeros([ntimes,nregs])
       sprd_regional_max[var] =np.zeros([ntimes,nregs])
       sprd_regional_min[var] =np.zeros([ntimes,nregs])

       varpos=[i for i,x in enumerate(ctl_dict['var_list']) if x == var][0]
       varsize=int( ctl_dict['var_size'][varpos] )
       if varsize == 0   :
          varsize = 1
       sprd_time_mean[var]=np.zeros([nx,ny,varsize])
       sprd_time_std[var]=np.zeros([nx,ny,varsize])



  #=========================================================
  #  START TIME LOOP
  #=========================================================

  it=0
  while ( ctime <= etime ):

    times[it]=it

    print( ctime )

    print ( 'Reading the sprd ')

    #=========================================================
    #  READ THE DATA
    #=========================================================

    my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ my_file_type + '/' + '/moment0002.grd'
 
    sprd=ctlr.read_data_grads(my_file,ctl_dict,masked=False)

    #Read the ensemble mean to get the information from the storm location.
    my_file=basedir + expname + ctime.strftime("%Y%m%d%H%M%S") + '/'+ my_file_type + '/' + '/moment0001.grd'
 
    ens_mean=ctlr.read_data_grads(my_file,ctl_dict,masked=False)

    #Compute max_dbz (we will use this to identify areas associated with clouds and convection)
    max_dbz = np.squeeze( np.nanmax(ens_mean['dbz'],2) )

    for my_var in sprd  :

        sprd[my_var][sprd[my_var]==undef] = np.nan
        sprd[my_var][sprd[my_var]> 100  ] = np.nan
 
        sprd[my_var] = np.sqrt( sprd[my_var] )

    #==========================================================
    # COMPUTE THE MASK 
    #==========================================================

    mask=dict()
    #Several mask can be defined here.
    mask['rain'] = np.logical_and( max_dbz > tr_rain  , np.logical_not( max_dbz == undef ) )    #Rain mask
    mask['no_rain'] = np.logical_and( max_dbz < tr_norain , np.logical_not( max_dbz == undef ) ) #No rain mask 
    mask['total'] = np.ones((nx,ny), dtype=bool)

    addvar=dict()
    #Contour variable that will be ploted for each mask.
    addvar['rain'] = np.abs( ens_mean['w'] )
    addvar['no_rain'] = np.abs( ens_mean['w'] )
    addvar['total']   = np.abs( ens_mean['w'] )

    addvar['rain'][ ens_mean['w'] == undef ] = undef
    addvar['no_rain'][ ens_mean['w'] == undef ] = undef
    addvar['total'][ ens_mean['w'] == undef ] = undef

    for my_mask in addvar   :
        addvar[my_mask][addvar[my_mask] == undef ] = np.nan

    #Exclude areas outside the radar domain.
    #tmp_mask=np.zeros((nx,ny), dtype=bool)
    #tmp_mask[buffer_zone_size:-buffer_zone_size,buffer_zone_size:-buffer_zone_size]=True
    #for my_mask in mask   :
    #    mask[my_mask]=np.logical_and( mask[my_mask] , tmp_mask )
    radar_mask = cmf.distance_range_mask( lon_radar , lat_radar , radar_range , lon , lat )
    #Combine with the previously defined masks.
    for my_mask in mask   :
        mask[my_mask]=np.logical_and( mask[my_mask] , radar_mask )


    #==========================================================
    # COMPUTE THE VERTICAL PROFILE
    #==========================================================
    if( ctime == itime ): #Initialize profiles

       vprofile=dict()
       for var in plot_variables  :
          if var in sprd    :
             vprofile[var]=dict()
             for my_mask in mask  :
                 vprofile[var][my_mask]=np.zeros([nlev,ntimes])

       addvprofile=dict()
       ndata      =dict()
       for my_mask in mask  :
          addvprofile[my_mask]=np.zeros([nlev,ntimes])
 
    #For each horizontal grid point accumulate the vertical profile discriminating between rain and no-rain
    for my_mask in mask :
       for ilev in range(0,nlev) :
          for var in plot_variables        :
             if var in sprd    :

                vprofile[var][my_mask][ilev,it]=np.nanmean( sprd[var][:,:,ilev][ mask[my_mask] ])

          addvprofile[my_mask][ilev,it]=np.nanmean( addvar[my_mask][:,:,ilev][ mask[my_mask] ] )
        

    ctime = ctime + delta
 
    it = it + 1

print ( "Finish time loop" )


#==========================================================
#  PLOTTING THE RESULTS
#==========================================================

#Compute the time averaged moments and plot them.
my_plotdir= plotbasedir + '/time_independent_plots/'

#Create output directory

levels=ctl_dict['vlevels']
levels_str=list()

#To solve the postprocessing error.
levels=np.delete(levels,4,axis=0)
levels[3]=850.0

#Get the level string list.
for ilev in levels  :
   levels_str.append( str(int(ilev)) )


for my_mask in mask  :

   for var in vprofile  :


      #Define my profile and sort it (to correct error in postprocessing)        
      my_profile=vprofile[var][my_mask]
      my_addvar =addvprofile[my_mask]

      #To solve the postprocessing error.
      my_profile=np.delete(my_profile,4,axis=0)
      my_addvar =np.delete(my_addvar ,4,axis=0)

      time_index=np.arange( 0 , my_profile.shape[1] ) * delta.total_seconds()

      cpf.set_default()  #Restore defaults
      my_map=cpf.cmap_discretize('YlGnBu',20)
      cpf.figconf['figpath']=my_plotdir
      cpf.figconf['figsize']=(12,10)
      cpf.figconf['titlefontsize']=20
      cpf.figconf['labelfontsize']=15
      cpf.figconf['contourf']=True
      cpf.figconf['shadedmin']=0.5
      if var == 'w'     :
         cpf.figconf['shadedmax']=2.0
      elif var == 'qv'  :
         cpf.figconf['shadedmax']=1.0e-3
         cpf.figconf['shadedmin'] =0.0
      elif var == 'tk'  :
         cpf.figconf['shadedmax']=0.5
         cpf.figconf['shadedmin'] =0.2
      elif var == 'dbz' :
         cpf.figconf['shadedmax']=10.0
         cpf.figconf['shadedmin']=0.5
      else              :
         cpf.figconf['shadedmax']=2.0
         
      cpf.figconf['shadedcolormap']=my_map
      cpf.figconf['colorbar']=True
      cpf.figconf['colorbarfontsize']=15
      cpf.figconf['axessize']=[0.1,0.1,0.8,0.8]
      cpf.figconf['contour']=True
      cpf.figconf['contourlevels']=np.arange(0,5.0,0.25)
      #cpf.figconf['contourflevels']=np.arange(0,0.04,0.005)
      cpf.figconf['contourcolormap']='Reds'
      cpf.figconf['gridline']=False
      cpf.figconf['ylabel']='Pressure (hPa)'
      cpf.figconf['xlabel']='Time (seconds) '
      cpf.figconf['xtick']=np.arange(0,np.max(time_index),300)
      cpf.figconf['ytick']=-np.log(levels)
      cpf.figconf['yticklabel']=levels_str

      cpf.figconf['title']='SPRD Profile ' + my_mask + ' for var ' + var
      cpf.figconf['figname']='Figure_sprd_profile_rain_var_' + var
      #cpf.figconf['show']=True

      my_profile[ my_profile > 100 ] = np.nan #There are som inf values in the reflectivity field.
      my_profile[ my_profile == undef ] = np.nan

      time_index=np.arange( 0 , my_profile.shape[1] ) * delta.total_seconds() 
      cpf.plot_x_y(  time_index  , -np.log(levels) , my_profile , my_addvar , my_addvar )

    
