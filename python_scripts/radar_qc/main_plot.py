#!/home/qcradar/.conda/envs/da/bin/python
import datetime as dt                #Datetime operations
import numpy as np                   #Numpy
import numpy.ma as ma
import matplotlib  
matplotlib.use('agg')
import matplotlib.pyplot as plt   
import pyart
import os
import pickle
#=========================================================================================================
# CONFIGURATION SECTION
#=========================================================================================================

#General section
qc_path = "/home/qcradar/scripts/"

datapath = '/home/qcradar/data/'  #Main data path.
datapath_out = '/home/qcradar/data/'    #Out data path
deltat = dt.timedelta( seconds=600 )   #Time window (seconds)
deltat_archive = dt.timedelta( seconds=86400 ) #Time window that will be kept in the remote ftp server.
time_offset = 0.0                      #Time offset (from current time)
instrument_list = ['RMA1','RMA2','RMA4','RMA6','RMA8','RMA11','PAR','PER','ANG']  #Instrument list.


#Ftp section
ftp_host='mate.cima.fcen.uba.ar'
ftp_user='ftp_alertar'
ftp_pass='Dra6h&b3wUDr'
ftp_path='qc_figures'
compress=False

#=========================================================================================================
# END OF CONFIGURATION SECTION
#=========================================================================================================

import sys
sys.path.append( qc_path + '/radar_qc/' )
sys.path.append( qc_path + '/radar_qc/src/python/' )
sys.path.append( qc_path + '/radar_qc/src/fortran/')
sys.path.append( qc_path + '/radar_so/fortran/')
sys.path.append( qc_path + '/radar_so/')

import operational_tools as ot       #Operational tools.
import conf_defaults as conf         #Radar qc default configuration

options = conf.options #This is the default configuration.

options['name_ref'] ='ZH'    #'dBZ'              #Reflectivity
options['name_cref']='CZH'   #Corrected reflectivity (qc output)
options['name_v']   ='VRAD'  #'V'                #Dopper velocity
options['name_cv']  ='CVRAD' #Corrected wind (qc ouput)
options['name_rho'] ='RHOHV' #'RhoHV'            #Rho HV

#Set the dates that will be processed. 
current_date = dt.datetime.utcnow()

ref_date=dt.datetime(current_date.year, current_date.month, 1, 0, 0, 0)
freqtimes = deltat.total_seconds()*np.floor((current_date-ref_date).total_seconds()/deltat.total_seconds())
c_end_date=( ref_date + dt.timedelta( seconds=freqtimes ) - deltat*2   ).strftime('%Y%m%d%H%M%S')
c_ini_date=( ref_date + dt.timedelta( seconds=freqtimes ) - deltat*3 ).strftime('%Y%m%d%H%M%S')

#Set the dates that will be archived. 
a_end_date=( ref_date + dt.timedelta( seconds=freqtimes ) + deltat         ).strftime('%Y%m%d%H%M%S')
a_ini_date=( ref_date + dt.timedelta( seconds=freqtimes ) - deltat_archive ).strftime('%Y%m%d%H%M%S')



print('')
print('=============================================================================')
print('We will plot all the files within the following dates:' )
print( c_ini_date )
print( c_end_date )
print('=============================================================================')
print('')

#Set the dates that will be archived in the remote server.
a_end_date=( ref_date + dt.timedelta( seconds=freqtimes ) ).strftime('%Y%m%d%H%M%S')
a_ini_date=( ref_date + dt.timedelta( seconds=freqtimes ) - deltat_archive ).strftime('%Y%m%d%H%M%S')

print('')
print('=============================================================================')
print(' GETTING FILE LIST ')
print('=============================================================================')
print('')


#Obtenemos la lista de archivos.
file_list = ot.get_file_list( datapath , c_ini_date , c_end_date , time_search_type='filename' , file_type_list = ['.nc'] )

#file_list = ['/home/qcradar/data/cfradial/cfrad.20181001_130839.RMA4.nc']
#file_list = []

fig_list=[]

for my_file in file_list :

      print('=============================================================================')
      print(' PLOT FILE : ',my_file)
      print('=============================================================================')
      print('')
 
      #Call QC routine


      radar=pyart.io.read( my_file )

      file_time   = dt.datetime.strftime( ot.get_time_from_filename( my_file ) , '%Y%m%d_%H%M%S' )     

      file_instrument = ot.get_instrument_type_from_filename( my_file ) 

      os.makedirs( datapath + '/figures/' ,exist_ok=True)

      if options['name_ref'] in radar.fields   :

         figname= datapath + '/figures/' + 'CFRAD_REF_' + file_instrument + '_' + file_time + '.png'

         plt.figure()

         plt.subplot(1,2,1)
         plt.pcolormesh( radar.fields[options['name_ref']]['data'] ,vmin=options['plot']['DbzMin'],vmax=options['plot']['DbzMax'],cmap=options['plot']['CmapDbz'])
         plt.subplot(1,2,2)
         plt.pcolormesh( radar.fields[options['name_cref']]['data'] ,vmin=options['plot']['DbzMin'],vmax=options['plot']['DbzMax'],cmap=options['plot']['CmapDbz'])
      
         plt.colorbar()
         #plt.show()
         plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',
            orientation='portrait', papertype=None, format=None,
            transparent=False, bbox_inches=None, pad_inches=0.1,
            frameon=None)

         fig_list.append( figname )

         figname= datapath + '/figures/' + 'CFRAD_REF_SWEEPS_' + file_instrument + '_' + file_time + '.png'

         display = pyart.graph.RadarDisplay(radar)

         figsize_x = 2.0 * radar.nsweeps 
         figsize_y = 4.0
         fig = plt.figure(figsize=(figsize_x,figsize_y))

         for isweept in range( radar.nsweeps )  :


            ax = fig.add_subplot(2,radar.nsweeps,isweept+1)    

            start_index = radar.sweep_start_ray_index['data'][isweept]
            end_index   = radar.sweep_end_ray_index['data'][isweept]

            # Datos de lat/lon
            lat = radar.gate_latitude['data'][start_index:end_index]
            lon = radar.gate_longitude['data'][start_index:end_index]

            # Variables a graficar.
            tmp_var = radar.fields[options['name_ref']]['data'][start_index:end_index]
            tmp_var = ma.masked_invalid(tmp_var)

            #Grafico las variables
            plt.pcolormesh(lon, lat, tmp_var,vmin=options['plot']['DbzMin'],vmax=options['plot']['DbzMax'],cmap=options['plot']['CmapDbz'])
            #plt.title('Elevation ' + str(elev) + '° - ' + hour + ':' + min + 'Z' + day + 'JAN' + year + '\n OBS', fontsize=fs_title)
            #cbar=plt.colorbar(shrink=1, extend = 'both')
            #cbar.cmap.set_under('0.75')
            #cbar.cmap.set_over('black')
            #cbar.set_label('[dBZ]', fontsize= fs_bar)

            #colorbar_flag = False         

            #display.plot_ppi(options['name_ref'], isweept, vmin=options['plot']['DbzMin'],vmax=options['plot']['DbzMax'],cmap=options['plot']['CmapDbz'],
            #                 colorbar_flag = colorbar_flag , title_flag = False , axislabels_flag = False )
            #display.plot_range_rings([60, 120, 240, 300])

            ax = fig.add_subplot(2,radar.nsweeps,radar.nsweeps+isweept+1)

            start_index = radar.sweep_start_ray_index['data'][isweept]
            end_index   = radar.sweep_end_ray_index['data'][isweept]

            # Datos de lat/lon
            lat = radar.gate_latitude['data'][start_index:end_index]
            lon = radar.gate_longitude['data'][start_index:end_index]
                             
            # Variables a graficar.
            tmp_var = radar.fields[options['name_cref']]['data'][start_index:end_index]
            tmp_var = ma.masked_invalid(tmp_var)

            #Grafico las variables
            plt.pcolormesh(lon, lat, tmp_var,vmin=options['plot']['DbzMin'],vmax=options['plot']['DbzMax'],cmap=options['plot']['CmapDbz'])

            #display.plot_ppi(options['name_cref'], isweept, vmin=options['plot']['DbzMin'],vmax=options['plot']['DbzMax'],cmap=options['plot']['CmapDbz'],
            #                 colorbar_flag = colorbar_flag , title_flag = False , axislabels_flag = False )

            #display.plot_range_rings([60, 120, 240, 300])


         plt.tight_layout()
         #plt.show()
         plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',
            orientation='portrait', papertype=None, format=None,
            transparent=False, bbox_inches=None, pad_inches=0.1,
            frameon=None)

         fig_list.append( figname )


      if options['name_v'] in radar.fields   :
      
         figname= datapath + '/figures/' + 'CFRAD_V_' + file_instrument + '_' + file_time + '.png'
 
         plt.figure()
         plt.subplot(1,2,1)
         plt.pcolormesh( radar.fields[options['name_v']]['data'] ,vmin=options['plot']['VrMin'],vmax=options['plot']['VrMax'],cmap=options['plot']['CmapWind'])
         plt.subplot(1,2,2)
         plt.pcolormesh( radar.fields[options['name_cv']]['data'] ,vmin=options['plot']['VrMin'],vmax=options['plot']['VrMax'],cmap=options['plot']['CmapWind'])
      
         plt.colorbar()
         plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',
            orientation='portrait', papertype=None, format=None,
            transparent=False, bbox_inches=None, pad_inches=0.1,
            frameon=None)
     
         fig_list.append( figname )

         figname= datapath + '/figures/' + 'CFRAD_V_SWEEPS_' + file_instrument + '_' + file_time + '.png'

         display = pyart.graph.RadarDisplay(radar)

         figsize_x = 2.0 * radar.nsweeps
         figsize_y = 4.0
         fig = plt.figure(figsize=(figsize_x,figsize_y))

         for isweept in range( radar.nsweeps )  :

            ax = fig.add_subplot(2,radar.nsweeps,isweept+1)


            start_index = radar.sweep_start_ray_index['data'][isweept]
            end_index   = radar.sweep_end_ray_index['data'][isweept]

            # Datos de lat/lon
            lat = radar.gate_latitude['data'][start_index:end_index]
            lon = radar.gate_longitude['data'][start_index:end_index]

            # Variables a graficar.
            tmp_var = radar.fields[options['name_v']]['data'][start_index:end_index]
            tmp_var = ma.masked_invalid(tmp_var)

            #Grafico las variables
            plt.pcolormesh(lon, lat, tmp_var,vmin=options['plot']['VrMin'],vmax=options['plot']['VrMax'],cmap=options['plot']['CmapWind'])

            #colorbar_flag = False

            #display.plot_ppi(options['name_v'], isweept,vmin=options['plot']['VrMin'],vmax=options['plot']['VrMax'],cmap=options['plot']['CmapWind'], 
            #                 colorbar_flag = colorbar_flag , title_flag = False , axislabels_flag = False )
            #display.plot_range_rings([60, 120, 240, 300])

            ax = fig.add_subplot(2,radar.nsweeps,radar.nsweeps+isweept+1)

            start_index = radar.sweep_start_ray_index['data'][isweept]
            end_index   = radar.sweep_end_ray_index['data'][isweept]

            # Datos de lat/lon
            lat = radar.gate_latitude['data'][start_index:end_index]
            lon = radar.gate_longitude['data'][start_index:end_index]

            # Variables a graficar.
            tmp_var = radar.fields[options['name_cv']]['data'][start_index:end_index]
            tmp_var = ma.masked_invalid(tmp_var)

            #Grafico las variables
            plt.pcolormesh(lon, lat, tmp_var,vmin=options['plot']['VrMin'],vmax=options['plot']['VrMax'],cmap=options['plot']['CmapWind'])


            #display.plot_ppi(options['name_cv'], isweept,vmin=options['plot']['VrMin'],vmax=options['plot']['VrMax'],cmap=options['plot']['CmapWind'],
            #                 colorbar_flag = colorbar_flag , title_flag = False , axislabels_flag = False )

            #display.plot_range_rings([60, 120, 240, 300])


         plt.tight_layout()
         #plt.show()
         plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',
            orientation='portrait', papertype=None, format=None,
            transparent=False, bbox_inches=None, pad_inches=0.1,
            frameon=None)

         fig_list.append( figname )

print('')
print('=============================================================================')
print(' GETTING PICKLE FILE LIST ')
print('=============================================================================')
print('')


#Obtenemos la lista de archivos.
file_list = ot.get_file_list( datapath , c_ini_date , c_end_date , time_search_type='filename' , file_type_list = ['.pkl'] )

#file_list = ['/home/qcradar/data/grid/PER_20181001160000.pkl']

print('')
print('=============================================================================')
print(' READING PICKLE FILES  ')
print('=============================================================================')
print('')


for my_file in file_list :

    print('Ploting : ' , my_file )


    file_time   = dt.datetime.strftime( ot.get_time_from_filename( my_file ) , '%Y%m%d_%H%M%S' )

    file_instrument = ot.get_instrument_type_from_filename( my_file )


    with open( my_file , 'rb') as filein  :
         so=pickle.load(filein)

    if 'grid_' + options['name_cref'] in so : 

       figname= datapath + '/figures/' + 'GRID_' + options['name_cref'] + '_LEVELS_' + file_instrument + '_' + file_time + '.png'

  
       figsize_x = 10
       figsize_y = 6
       fig = plt.figure(figsize=(figsize_x,figsize_y))

       tmp_var = np.ma.masked_array( so['grid_' + options['name_cref']]['data'] , mask=  so['grid_' + options['name_cref']]['nobs'] == 0 )

       for ilevel in range( (so['grid_' + options['name_cref']]['data']).shape[0] )  :

            ax = fig.add_subplot(3,5,ilevel+1)

            plt.pcolormesh( 10.0*np.log10(tmp_var[ilevel,:,:]) ,vmin=options['plot']['DbzMin'],vmax=options['plot']['DbzMax'],cmap=options['plot']['CmapDbz'])

       plt.tight_layout()
       plt.colorbar()
       #plt.show()
       plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',
          orientation='portrait', papertype=None, format=None,
          transparent=False, bbox_inches=None, pad_inches=0.1,
          frameon=None)

       fig_list.append( figname )

       plt.close()

#       for key in so['grid_' + options['name_cref']]  :

#          if key == 'az' or key == 'ra' or key == 'el'  :

#               max_var=np.max( so['grid_' + options['name_cref']][key] )
#               min_var= 0.0

#               figname= datapath + '/figures/' + 'GRID_' + options['name_cref'] + '_' + key + '_LEVELS_' + file_instrument + '_' + file_time + '.png'

#               figsize_x = 10
#               figsize_y = 6
#               fig = plt.figure(figsize=(figsize_x,figsize_y))

#               tmp_var = np.ma.masked_array( so['grid_' + options['name_cref']][key] , mask=  so['grid_' + options['name_cref']]['nobs'] == 0 )

#               for ilevel in range( (so['grid_' + options['name_cref']][key]).shape[0] )  :

#                   ax = fig.add_subplot(3,5,ilevel+1)
#                   plt.pcolormesh( 10.0*np.log10( tmp_var[ilevel,:,:] ) , vmin = min_var , vmax = max_var  )

#               plt.tight_layout()
#               plt.colorbar()
#               plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',
#                   orientation='portrait', papertype=None, format=None,
#                   transparent=False, bbox_inches=None, pad_inches=0.1,
#                   frameon=None)

#               fig_list.append( figname )
#               plt.close()



    file_time   = dt.datetime.strftime( ot.get_time_from_filename( my_file ) , '%Y%m%d_%H%M%S' )

    file_instrument = ot.get_instrument_type_from_filename( my_file )


    if 'grid_' + options['name_cv'] in so :

       figname= datapath + '/figures/' + 'GRID_' + options['name_cv'] + '_LEVELS_' + file_instrument + '_' + file_time + '.png'


       figsize_x = 10
       figsize_y = 6
       fig = plt.figure(figsize=(figsize_x,figsize_y))

       tmp_var = np.ma.masked_array( so['grid_' + options['name_cv']]['data'] , mask=  so['grid_' + options['name_cv']]['nobs'] == 0 )

       for ilevel in range( (so['grid_' + options['name_cv']]['data']).shape[0] )  :

            ax = fig.add_subplot(3,5,ilevel+1)

            plt.pcolormesh( tmp_var[ilevel,:,:] ,vmin=options['plot']['VrMin'],vmax=options['plot']['VrMax'],cmap=options['plot']['CmapWind'] )

       plt.tight_layout()
       plt.colorbar()
       #plt.show()
       plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',
          orientation='portrait', papertype=None, format=None,
          transparent=False, bbox_inches=None, pad_inches=0.1,
          frameon=None)
       plt.close()

       fig_list.append( figname )

       for key in so['grid_' + options['name_cv']]  :

           if key == 'az' or key == 'ra' or key == 'el'  :

               max_var=np.max( so['grid_' + options['name_cref']][key] )
               min_var= 0.0  

               figname= datapath + '/figures/' + 'GRID_' + options['name_cv'] + '_' + key + '_LEVELS_' + file_instrument + '_' + file_time + '.png'

               figsize_x = 10
               figsize_y = 6
               fig = plt.figure(figsize=(figsize_x,figsize_y))

               tmp_var = np.ma.masked_array( so['grid_' + options['name_cv']][key] , mask=  so['grid_' + options['name_cv']]['nobs'] == 0 )

               for ilevel in range( (so['grid_' + options['name_cv']][key]).shape[0] )  :

                   ax = fig.add_subplot(3,5,ilevel+1)
                   plt.pcolormesh( tmp_var[ilevel,:,:] , vmin=min_var,vmax=max_var)

               plt.tight_layout()
               plt.colorbar()
               plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',
                   orientation='portrait', papertype=None, format=None,
                   transparent=False, bbox_inches=None, pad_inches=0.1,
                   frameon=None)

               fig_list.append( figname )
               plt.close()


#fig_list=[]

#for file in file_list :

print('')
print('=============================================================================')
print(' UPLOADING DATA  ')
print('=============================================================================')
print('')

#Call remote server uploading routine

ot.upload_to_ftp( fig_list , ftp_host, ftp_user, ftp_pass , ftp_path , compress=compress )




