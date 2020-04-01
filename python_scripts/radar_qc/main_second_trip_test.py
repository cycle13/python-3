#!/home/qcradar/.conda/envs/da/bin/python
import datetime as dt                #Datetime operations
import numpy as np                   #Numpy
import os
import scipy
import matplotlib.pyplot as plt
import pyart

PlotData=True
WriteData=False

#=========================================================================================================
# CONFIGURATION SECTION
#=========================================================================================================

#General section
qc_path = "/home/jruiz/Dropbox/DATA/"

datapath = '/home/jruiz/Dropbox/DATA/TMP_DATOS_RADAR/cfradial/'  #Main data path.
datapath_out = '/home/jruiz/Dropbox/DATA/TMP_DATOS_RADAR/cfradial/'    #Out data path
#deltat = dt.timedelta( seconds=600 )   #Time window (seconds)
#deltat_archive = dt.timedelta( seconds=86400 ) #Time window that will be kept in the remote ftp server.
#time_offset = 0.0                      #Time offset (from current time)
instrument_list = ['RMA1','RMA2','RMA4','RMA6','RMA8','RMA11']  #Instrument list.
#instrument_list = [ 'ANG' ]

time_threshold= 600  #Maximun difference time between volumens in seconds.

c_ini_date = '20181110000000'
c_end_date = '20181111000000'

file_type_list = ['.H5','.vol','.nc']

write_cfradial=True                    #Write a cfradial file with the qced data

#Second trip correction section.
name_ref='ZH'                                               #Reflectivity name
name_v  ='VRAD'                                             #Doppler velocity name

#=========================================================================================================
# END OF CONFIGURATION SECTION
#=========================================================================================================

import sys
sys.path.append( qc_path + '/radar_qc/' )
sys.path.append( qc_path + '/radar_qc/src/python/' )
sys.path.append( qc_path + '/radar_qc/src/fortran/')
sys.path.append( qc_path + '/radar_so/fortran/')

import radar_qc_module as rqc        #Radar qc python modules
import conf_defaults   as conf         #Radar qc default configuration
import operational_tools as ot       #Operational tools.
from   common_qc_tools  import qc  #Fortran code routines.

print('')
print('=============================================================================')
print('We will process all the files within the following dates:' )
print( c_ini_date )
print( c_end_date )
print('=============================================================================')
print('')

print('')
print('=============================================================================')
print(' GETTING FILE LIST ')
print('=============================================================================')
print('')


#Obtenemos la lista de archivos.

for my_instrument in instrument_list  :

   my_instrument_list = [ my_instrument ]

   #Get the file list corresponding to the current instrument.
   file_list = ot.get_file_list( datapath , c_ini_date , c_end_date , time_search_type='filename' , file_type_list = file_type_list , instrument_type_list = my_instrument_list )


   print('')
   print('=============================================================================')
   print(' PROCESING ' , my_instrument)
   print('=============================================================================')
   print('')


   #1)For each file in the list identify those with the strategy 301_02
   #2)Get the volume after and before this one (301_01)
   #3)Average these two and compare them to the 301_02.

   for my_file in file_list :

       if '0301_02' in my_file  :  #301-2 STRATEGY (this is the one we need to correct)

          current_file = my_file

          current_time = ot.get_time_from_filename( my_file )

          #Find the volume 301_01 after.

          file_before = None
          file_after  = None

          seconds_before = time_threshold
 
          seconds_after  = time_threshold

          for my_file_2 in file_list :

             if '0301_01' in my_file_2  : 

                tmp_time = ot.get_time_from_filename( my_file_2 )

                time_delta = ( tmp_time - current_time ).total_seconds() 
            
                if time_delta > 0 and time_delta <= seconds_after  :

                   seconds_after = time_delta 

                   file_after = my_file_2

                if time_delta < 0 and np.abs( time_delta ) >= seconds_before  :

                   seconds_before = time_delta 

                   file_before = my_file_2
                   
          print('Current file ', current_file)
          print('File after ',file_after)
          print('File before ',file_before)


          my_format = ot.get_format_from_filename( current_file )

          current_radar = ot.read_file( current_file , my_format )

          current_radar = ot.rename_fields( current_radar )

          current_ref = dict()

          undef = current_radar.fields['ZH']['_FillValue']

          current_range = current_radar.range['data'][:]

          [ current_ref['ref'] , current_ref['az'] , current_ref['level'] , current_ref['time'] , current_ref['az_exact'] ]=rqc.order_variable( current_radar , 'ZH' , undef )

          [ current_ref['dv'] , current_ref['az'] , current_ref['level'] , current_ref['time'] , current_ref['az_exact'] ]=rqc.order_variable( current_radar , 'VRAD' , undef )


          print( np.max( current_ref['ref'] ) , np.min( current_ref['ref'] ) ) 

          if file_before != None  :
             my_format = ot.get_format_from_filename( file_before )
             radar_bf = ot.read_file( file_before , my_format )
             radar_bf = ot.rename_fields( radar_bf )
           
             bf_ref = dict()
             bf_range = radar_bf.range['data'][:]

             [ bf_ref['ref'] , bf_ref['az'] , bf_ref['level'] , bf_ref['time'] , bf_ref['az_exact'] ]=rqc.order_variable( radar_bf , 'ZH' , undef )
 
             bf_ref['ref'] = bf_ref['ref'][:,:,[0,2,4]]

             
             #print( current_ref['level'] , bf_ref['level'][[0,2,4]] )

          if file_after  != None  :
             my_format = ot.get_format_from_filename( file_after )
             radar_af = ot.read_file( file_after , my_format )
             radar_af = ot.rename_fields( radar_af )

             af_ref = dict()
             af_range = radar_af.range['data'][:]

             [ af_ref['ref'] , af_ref['az'] , af_ref['level'] , af_ref['time'] , af_ref['az_exact'] ]=rqc.order_variable( radar_af , 'ZH' , undef )

             af_ref['ref'] = af_ref['ref'][:,:,[0,2,4]]

             #print( current_ref['level'] , af_ref['level'][[0,2,4]] )

          proceed = False
          if file_after != None and file_before != None :

             ctrl_ref = 0.5*( af_ref['ref'] + bf_ref['ref'] )   
             ctrl_range = af_range
             proceed = True 

          if file_after != None and file_before == None :

             ctrl_ref = af_ref['ref']
             ctrl_range = af_range
             proceed = True
             print( 'Warning : no file before ' )

          if file_after == None and file_before != None :

             ctrl_ref = bf_ref['ref']
             ctrl_range = bf_range
             proceed = True

             print( 'Warning : no file after ' )

          if file_after == None and file_before == None :

             proceed = False 

             print( 'Warning : no file after or before ... doing nothing ')

          if proceed            :

              print( np.shape( ctrl_ref ) )

              #Interpolate ctrl_ref to the reference    

              azimuths = current_ref['az']
              ctrl_ref_120 = np.zeros( np.shape( current_ref['ref'] ) )

              for ilev in range( 0 , 3 )   :
                 interpolator = scipy.interpolate.interp2d( ctrl_range , azimuths , ctrl_ref[:,:,ilev] , kind='linear', copy=True)  
                 ctrl_ref_120[:,:,ilev] = interpolator( current_range , azimuths )
            
              #Filter current reflectivity based on ctrl_ref_120 data.

              mask = np.zeros( np.shape( current_ref['ref'] ) ).astype(bool)

              #Take the local maximum of ctrl_reflectivity to take account cell movement between different times.
              tmp_ref = qc.box_functions_2d(datain=ctrl_ref_120 ,na=np.shape( ctrl_ref_120)[0],nr=np.shape( ctrl_ref_120)[1],ne=np.shape( ctrl_ref_120)[2],undef=undef
                                                   ,boxx=2,boxy=2,boxz=0,operation='MAXN',threshold=0.0)

              diff_ref = current_ref['ref'] - tmp_ref

              #Where current ctrl reflectivity is very small and difference between both ref is large then eliminate the echos in the current volum
              mask[ np.logical_and( diff_ref > 10 , tmp_ref < 10.0 ) ]=True


              current_ref_old = np.copy( current_ref['ref'][:] )
              current_dv_old  = np.copy( current_ref['dv'][:] )


              current_ref['ref'][ mask ] = undef 
              current_ref['dv'][ mask ]  = undef
              

              #Replace the original field by the corrected field.             
              tmp=rqc.order_variable_inv( current_radar , current_ref['ref'] , undef )
              current_radar.fields['ZH']['data']=np.ma.masked_array(tmp , mask = (tmp==undef) )

              tmp=rqc.order_variable_inv( current_radar , current_ref['dv'] , undef )
              current_radar.fields['VRAD']['data']=np.ma.masked_array(tmp , mask = (tmp==undef) )

              #Overwrite data clev in options['plot']['Elevs']  :
              if PlotData :

                 for ilev in range( 0 , 3 ) :
                    plt.figure(figsize=(15, 5))
                    plt.subplot(1,3,1)

                    plt.pcolor(current_ref_old[:,:,ilev],vmin=-10,vmax=70,cmap='pyart_NWSRef')
                    plt.title('Original Reflectivity')
                    plt.colorbar()

                    plt.subplot(1,3,2)
                    plt.pcolor(current_ref['ref'][:,:,ilev],vmin=-10,vmax=70,cmap='pyart_NWSRef')
                    plt.title('Corrected Reflectivity')
                    plt.colorbar()

                    plt.subplot(1,3,3)
                    plt.pcolor(tmp_ref[:,:,ilev],vmin=-10,vmax=70,cmap='pyart_NWSRef')
                    plt.title('Reference reflectivity')
                    plt.colorbar()

                    plt.figure(figsize=(10, 5))
                    plt.subplot(1,2,1)

                    plt.pcolor(current_dv_old[:,:,ilev],vmin=-30,vmax=30,cmap='pyart_NWSVel')
                    plt.title('Original Doppler')
                    plt.colorbar()

                    plt.subplot(1,2,2)
                    plt.pcolor(current_ref['dv'][:,:,ilev],vmin=-30,vmax=30,cmap='pyart_NWSVel')
                    plt.title('Corrected Doppler')
                    plt.colorbar()


                    plt.show()

                 #Overwrite the data corresponding to the current volume.
                 if WriteData   :
                    #ot.save_cfradial( current_file , current_radar )
                    pyart.io.cfradial.write_cfradial(current_file , current_radar, format='NETCDF4', time_reference=None, arm_time_variables=False)

             





