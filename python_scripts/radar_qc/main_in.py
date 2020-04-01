#!/home/qcradar/.conda/envs/da/bin/python
import datetime as dt                #Datetime operations
import numpy as np                   #Numpy
import os

#=========================================================================================================
# CONFIGURATION SECTION
#=========================================================================================================

#General section
qc_path = "/home/qcradar/scripts/"

datapath = '/ms-36/mrugna/RMA/datos/'  #Main data path.
datapath_out = '/home/qcradar/data/'    #Out data path
deltat = dt.timedelta( seconds=600 )   #Time window (seconds)
deltat_archive = dt.timedelta( seconds=86400 ) #Time window that will be kept in the remote ftp server.
time_offset = 0.0                      #Time offset (from current time)
instrument_list = ['RMA1','RMA2','RMA4','RMA6','RMA8','RMA11','PAR','PER','ANG']  #Instrument list.
#instrument_list = [ 'ANG' ]

file_type_list = ['.H5','.vol','.nc']

remove_local_pkl = True                #Remove intermediate gridded data in pkl format.
remove_local_dat = False               #Remove gridded data in letkf format.
remove_remote_dat = True               #Remove remote letkf files.

write_cfradial=True                    #Write a cfradial file with the qced data
#Qc section
name_ref='ZH'                                               #Reflectivity name
name_v  ='VRAD'                                             #Doppler velocity name
name_rho='RHOHV'                                            #RhoHV name
toporawdatapath=qc_path + "/data/terrain_data/raw/"         #Raw topography data path (tif format)
toporadardatapath=qc_path + "/data/terrain_data/radar/"     #Interpolated topography data path (unformatted)

#Superobbing section 
output_freq = 600
#        dx    dz   zmax  rmax
grid = [10000, 1000, 15e3, 240e3]
#opts = {'CZH': [4001, 5, 0], 'CVRAD': [4002, 2]}
opts = {'CZH': [4001, 5, 0]}
outputpath = datapath_out

#Ftp section
ftp_host='ftp.smn.gob.ar'
ftp_user='rra'
ftp_pass='TwV/27gm7YQsw'
ftp_path='radar'
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

import radar_qc_module as rqc        #Radar qc python modules
import conf_defaults   as conf         #Radar qc default configuration
import operational_tools as ot       #Operational tools.
import radar_so as so                #Superobbing module

options = conf.options #This is the default configuration.
options['name_ref'] = name_ref         #Reflectivity
options['name_v']   = name_v           #Dopper velocity
options['name_rho'] = name_rho         #Rho HV
options['toporawdatapath']= toporawdatapath
options['toporadardatapath']= toporadardatapath 

#Set the dates that will be processed. 
current_date = dt.datetime.utcnow()

ref_date=dt.datetime(current_date.year, current_date.month, 1, 0, 0, 0)
freqtimes = deltat.total_seconds()*np.floor((current_date-ref_date).total_seconds()/deltat.total_seconds())
c_end_date=( ref_date + dt.timedelta( seconds=freqtimes ) - deltat   ).strftime('%Y%m%d%H%M%S')
c_ini_date=( ref_date + dt.timedelta( seconds=freqtimes ) - deltat*2 ).strftime('%Y%m%d%H%M%S')

#Set the dates that will be archived. 
a_end_date=( ref_date + dt.timedelta( seconds=freqtimes ) + deltat         ).strftime('%Y%m%d%H%M%S')
a_ini_date=( ref_date + dt.timedelta( seconds=freqtimes ) - deltat_archive ).strftime('%Y%m%d%H%M%S')


print('')
print('=============================================================================')
print('We will process all the files within the following dates:' )
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
file_list = ot.get_file_list( datapath , c_ini_date , c_end_date , time_search_type='filename' , file_type_list = file_type_list )

#file_list = ['/home/qcradar/data//cfradial/cfrad.20181028_033655.RMA2.nc']
#file_list = ['/home/qcradar/data/cfradial/cfrad.20181028_032427.ANG.nc']
#file_list = ['/home/qcradar/data/cfradial/cfrad.20181027_184006.PAR.nc']

print(file_list)

print('')
print('=============================================================================')
print(' READING FILE LIST ')
print('=============================================================================')
print('')

#Obtenemos la lista de objetos radares.
radar_list = ot.read_multiple_files(  file_list , instrument_list )
my_updated_dirs =  []

my_updated_tars =  []

for radar in radar_list :

      print('=============================================================================')
      print(' QC')
      print('=============================================================================')
      print('')
 
      #Call QC routine

      [ radar , qc_output ] = rqc.main_qc( options , radar )

      #Save data in cfradial format

      if write_cfradial  :

         print('=============================================================================')
         print(' WRITING QC OUTPUT IN CFRADIAL FORMAT')
         print('=============================================================================')
         print('')


         ot.save_cfradial( datapath_out + '/cfradial/' , radar )

      print('')
      print('=============================================================================')
      print(' SUPEROBBING')
      print('=============================================================================')
      print('')

      #Call SO routine 

      letkf_filelist = so.main_radar_so(radar, output_freq, grid, opts, datapath_out  )

      for my_file in letkf_filelist :

          my_time_datetime = ot.get_time_from_filename( my_file )

          my_time = dt.datetime.strftime( my_time_datetime  + dt.timedelta( seconds=3600) , '%Y%m%d_%H' )

          my_minute = dt.datetime.strftime( my_time_datetime , '%M' )
 
          complete_path = datapath_out + '/radar/' + my_time

          if not complete_path in my_updated_dirs :

             my_updated_dirs.append( complete_path ) 
  
          if not os.path.isdir( complete_path )  :

             os.makedirs( complete_path )

          os.system('ln -sf ' + my_file + ' ' + complete_path + '/' + os.path.basename(my_file) )

          if my_minute == '00'   :  #Copy the file in the next hour folder as well.

             my_time = dt.datetime.strftime( my_time_datetime , '%Y%m%d_%H' )

             complete_path = datapath_out + '/radar/' + my_time

             if not complete_path in my_updated_dirs :

                my_updated_dirs.append( complete_path )

             if not os.path.isdir( complete_path )  :

                os.makedirs( complete_path )

             os.system('ln -sf ' + my_file + ' ' + complete_path + '/' + os.path.basename(my_file) )
             

print('')
print('=============================================================================')
print(' UPLOADING FILES TO REMOTE FTP SERVER ' + ftp_host )
print('=============================================================================')
print('')

for my_dir in my_updated_dirs  :

   my_tar_file = my_dir + '.tar.gz'
 
   if not my_tar_file in my_updated_tars :

      my_updated_tars.append( my_tar_file )

   os.system('rm -f ' + my_tar_file )
   os.system('tar --dereference -czvf ' + my_tar_file + ' -C ' + my_dir + ' .') 
   #os.system('gzip -f ' + my_tar_file )


ot.upload_to_ftp( my_updated_tars , ftp_host, ftp_user, ftp_pass , ftp_path , compress=compress ) 

for my_tar_file in my_updated_tars :

    os.system('rm -f ' + my_tar_file)


print('')
print('=============================================================================')
print('We will keep all the files within the following dates:' )
print( a_ini_date )
print( a_end_date )
print('=============================================================================')
print('')

print('')
print('=============================================================================')
print(' REMOVE OLD FILES FROM REMOTE SERVER ' + ftp_host )
print('=============================================================================')
print('')

#Call remote server deleting routine

if remove_remote_dat :

   ot.remove_from_ftp_timebased( ftp_host, ftp_user, ftp_pass , ftp_path , a_ini_date , a_end_date , file_format_list = ['tgz'] ) 

print('')
print('=============================================================================')
print(' REMOVE OLD FILES FROM LOCAL SERVER ' + datapath_out )
print('=============================================================================')
print('')

#Call local server deleting routine 
tmp_format_list=[]
if remove_local_pkl : 
   tmp_format_list.append('pickle')
if remove_local_dat :
   tmp_format_list.append('letkf')

ot.remove_from_localpath_timebased( datapath_out , a_ini_date , a_end_date , file_format_list = tmp_format_list )




