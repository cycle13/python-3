from . import io
#from . import pawr_read_pyart as pawr_read 
from . import pawr_read 
#import calc
#import pyart
from datetime import datetime as dt
from datetime import timedelta 
import numpy as np
import glob
import os


def get_radar_files(  data_path , data_type='pawr' ) :


   if data_type == 'pawr' :
        file_prefix='PAWR'
        file_sufix ='.dat'

   radar_files_dict=dict()

   file_list=glob.glob( data_path + '/' + file_prefix + '*' + file_sufix )
   file_list.sort()

   radar_files_dict['times'] = list()
   for cfile in file_list :
       times = list()

       if( data_type == 'pawr' ) :

           filename = os.path.basename( cfile ) 
           date_str = filename[19:27] + filename[28:34]
           #Read the date and convert into datetime format. Substract nine hours to convert to UTC time. 
           radar_files_dict['times'].append( dt.strptime( date_str , '%Y%m%d%H%M%S' ) + timedelta( hours = -9.0 ) )

   radar_files_dict['files'] = file_list

   return radar_files_dict

def get_radar_data( requested_time , data_path , radar_files = None , time_tol = 15 , data_type = 'pawr' ) :

    #If file list is not present, then generate the list of files.
    if radar_files is None :
        radar_files = get_radar_files( data_path , data_type=data_type )

    #Compute time distance
    time_dist = np.zeros( len( radar_files['times'] ) )
    for ii,my_time in  enumerate( radar_files['times'] ) :
        time_dist[ii] = ( my_time - requested_time ).total_seconds()
    if np.abs( time_dist ).min() <= time_tol  :
       print('Found radar data at ', dt.strftime( radar_files['times'][np.abs(time_dist).argmin()] , '%Y%m%d%H%M%S' ),' to be compared with model data valid at ',dt.strftime( requested_time , '%Y%m%d%H%M%S' ) )
       #radar = pawr_read.pawr_read( radar_files['files'][np.abs(time_dist).argmin()]  )
       radar = pawr_read.pawr_read( radar_files['files'][np.abs(time_dist).argmin()]  )
    else                            :
       radar = None

    return radar , radar_files



