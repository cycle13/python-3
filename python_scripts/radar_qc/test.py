#!/home/qcradar/.conda/envs/da/bin/python
qc_path = "/home/qcradar/scripts/"

import sys
sys.path.append( qc_path + '/radar_qc/' )
sys.path.append( qc_path + '/radar_qc/src/python/' )
sys.path.append( qc_path + '/radar_qc/src/fortran/')
sys.path.append( qc_path + '/radar_so/fortran/')
sys.path.append( qc_path + '/radar_so/')

import radar_qc_module as rqc        #Radar qc python modules
import conf_defaults as conf         #Radar qc default configuration
import operational_tools as ot       #Operational tools.
import datetime as dt                #Datetime operations
import numpy as np                   #Numpy
import radar_so as so                #Superobbing module

import conf_defaults as conf

datapath = '/ms-36/mrugna/RMA/datos/'  #Main data path.
datapath_out = '/home/qcradar/data'    #Out data path
deltat = dt.timedelta( seconds=1800 )   #Time window (seconds)
instrument_list = ['ANG','PAR','PER']

file_type_list = ['.h5','.vol','.nc']

current_date = dt.datetime.utcnow()

ref_date=dt.datetime(current_date.year, current_date.month, 1, 0, 0, 0)
freqtimes = deltat.total_seconds()*np.floor((current_date-ref_date).total_seconds()/deltat.total_seconds())
c_end_date=( ref_date + dt.timedelta( seconds=freqtimes ) ).strftime('%Y%m%d%H%M%S')
c_ini_date=( ref_date + dt.timedelta( seconds=freqtimes ) - deltat ).strftime('%Y%m%d%H%M%S')


print('')
print('=============================================================================')
print('We will process all the files within the following dates:' )
print( c_ini_date )
print( c_end_date )
print('=============================================================================')
print('')

#Obtenemos la lista de archivos.
file_list = ot.get_file_list( datapath , c_ini_date , c_end_date , time_search_type='filename' , file_type_list = file_type_list )

for i in file_list :
   print (i)
#Obtenemos la lista de objetos radares.
radar_list = ot.read_multiple_files(  file_list , instrument_list )


