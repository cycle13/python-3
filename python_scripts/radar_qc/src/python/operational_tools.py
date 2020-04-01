#print __doc__
# Author: Rapid Refresh Argentina Team
# License: BSD 3 clause
import datetime as dt
import os
import numpy as np

def read_multiple_files(  file_list , instrument_list = None )          :

   #import numpy as np
   #import datetime as dt
   #file_list is a list of filenames corresponding to radar files in different formats.
   #instrument_list (if present) a list of instrument names. Only names in instrument_list will be incorporated to the 
   #radar object list.

   #Generate a list of radar objects
   dtype = 'float64'  #On output radar object fields will have this precission.

   radar_list   = []                                       #Final list of radar objects.

   used_file = np.zeros( len( file_list ) ).astype(bool)    #Keep record of which files has been used

   for ifile , filename in enumerate( file_list )    :

      file_format = get_format_from_filename( filename )
 
      file_time   = get_time_from_filename( filename )

      file_instrument = get_instrument_type_from_filename( filename )

      #print('Main',file_format,file_time,file_instrument)

      if ( ( file_instrument in instrument_list ) or ( instrument_list == None ) ) and ( not used_file[ifile] ) and ( file_format != None )  and ( file_time != None ) :
         #Read the radar 
         my_radar = read_file( filename , file_format )  

         if my_radar != None  :
            #Add the current radar to the list. 
            #radar_list.append = my_radar
            my_radar.files = [ filename ]
            used_file[ ifile ] = True

            #Check if we can associate other radars in the list to the current radar.
  
            for tmp_ifile , tmp_filename  in enumerate( file_list )   :

		       #print('Second loop',get_format_from_filename( tmp_filename ),get_time_from_filename( tmp_filename ),get_instrument_type_from_filename( tmp_filename ) , used_file[tmp_ifile])
  
               if ( not used_file[tmp_ifile] ) and ( ifile != tmp_ifile )   :
                  if ( ( file_format     ==  get_format_from_filename( tmp_filename ) ) and 
                       ( file_time       ==  get_time_from_filename( tmp_filename ) ) and
                       ( file_instrument ==  get_instrument_type_from_filename( tmp_filename ) ) )  :

                     #Read the data
                     tmp_radar = read_file( tmp_filename , file_format )
                     if tmp_radar !=  None  :
                        
                        #Check if we can merge tmp_radar and my_radar objects.
                        [ my_radar , merged ] = merge_radar_object( my_radar , tmp_radar )
                        if merged   :
                           my_radar.files.append( tmp_filename )
                           used_file[ tmp_ifile ] = True
                           
                        else        :       
                           print('Warning: Inconsistent shapes found for ' + file_instrument + ' ' + file_time.strftime('%Y-%m-%d-%H-%M-%S') )
                           print('This volume will be processed separately')  
                     if tmp_radar == None  :
                           used_file[ tmp_ifile ] = True

            #So far my_radar contains all the variables corresponding to this instrument and initial time.
            my_radar = get_strat( filename , my_radar )  #Additional metadata that will be required by QC
            my_radar = rename_fields( my_radar )         #Rename fields different instruments share the same name convention.
          
            radar_list.append( my_radar )       

            print('RADAR : ' + my_radar.metadata['instrument_name'] + ' ' + file_time.strftime('%Y-%m-%d-%H-%M-%S') )
            for ifile in my_radar.files   :
                print('   FILE: ' + ifile )
         else       :
             used_file[ifile]=True 
   return radar_list

def get_format_from_filename( filename )            :

   file_format = None 

   if ('.h5' in filename ) or ( '.H5' in filename )    :
      file_format = 'h5'
   if ( '.vol' in filename ) or ( '.VOL' in filename ) :
      file_format = 'vol' 
   if ( '.nc'  in filename ) or ( '.NC' in filename )  or ( 'cfrad' in filename ) :
      file_format = 'cfrad' 
   if ( '.dat' in filename ) or ( '.DAT' in filename)  :
      file_format = 'letkf'
   if ( '.pkl' in filename ) or ( '.PKL' in filename)  :
      file_format = 'pickle' 
   if ( '.tar.gz' in filename ) or ( '.TAR.GZ' in filename ) :
      file_format = 'tgz'

   return file_format

def read_file( filename , format_file )    :
   import pyart
   from pyart.aux_io.sinarame_h5 import read_sinarame_h5
   from pyart.aux_io.rainbow_wrl import read_rainbow_wrl

   radar = None


   try   :
      if format_file == 'h5'   :
         radar = read_sinarame_h5(filename, file_field_names=True)
      if format_file == 'vol'  :
         radar = read_rainbow_wrl(filename, file_field_names=True)
      if format_file == 'cfrad'   :
         radar = pyart.io.read(filename)

      print( ' ' )
      print( '=============================================================================')
      print( 'Reading file:' + filename )
      print( '=============================================================================')
      print( ' ' )
      for my_key in radar.fields   :    
         print( 'Found the following variable: ' + my_key )
      print( ' ' )
   except  :
      print('Warning: Could not read file ' + filename )
      radar = None

   return radar

def rename_fields ( radar )  :

  #Unify different names for different variables depending on the data source.
  #If the intended name is already present, then we will use this variable 
  #and ignore any other variable that may contain a similar field.

  if not 'ZH' in radar.fields  :

     if 'TH' in radar.fields    :
        radar.fields['ZH'] = radar.fields.pop('TH')

     if 'DBZH' in radar.fields  :
        radar.fields['ZH'] = radar.fields.pop('DBZH')

  if not 'VRAD' in radar.fields  :
     if 'V' in radar.fields     :
        radar.fields['VRAD'] = radar.fields.pop('V')

     if 'dBZ' in radar.fields    :
        radar.fields['ZH'] = radar.fields.pop('dBZ')

  if not 'WRAD' in radar.fields  :
     if 'W' in radar.fields     :
        radar.fields['WRAD'] = radar.fields.pop('W')

  if not 'RHOHV' in radar.fields  :
     if 'RhoHV' in radar.fields :
        radar.fields['RHOHV'] = radar.fields.pop('RhoHV') 

  return radar   

def get_strat ( filename , radar )  :

#Include some parameters that are required by the QC module.

    #import numpy as np
    import numpy.ma as ma
    #import os

    local_fill_value = -9999.0
    levels=np.unique(radar.elevation['data'])

    #Set the instrument name 
    radar.metadata['instrument_name'] = get_instrument_type_from_filename( filename )

    #Add missing structures to the radar object.
    if radar.altitude_agl == None :
       radar.altitude_agl = dict()
    if radar.metadata == None :
       radar.metadata = dict()
    if radar.instrument_parameters == None :
       radar.instrument_parameters = dict()

    #Get the corresponding radar strategy depending on the filename.
    radar.altitude_agl['data'] = 0.0


    #TODO: Ver si podemos calcular esto en funcion de lo que esta
    #en la estructura radar en lugar de tener que sacarlo del nombre del archivo.
    if 'RMA' in filename  :


       if '9005_01' in filename  :  #9005-1 STRATEGY
          nyquist_velocity     = 6.63 * np.ones( np.shape( np.unique( radar.elevation['data'] ) ) )
       if '9005_02' in filename  :  #9005-2 STRATEGY
          nyquist_velocity     = 33.04 * np.ones( np.shape( np.unique( radar.elevation['data'] ) ) )
       if '9005_03' in filename  :  #9005-3 STRATEGY
          nyquist_velocity     = 3.98 * np.ones( np.shape( np.unique( radar.elevation['data'] ) ) )

       if '0117_01' in filename  :  #122-1 STRATEGY
          nyquist_velocity     = 6.63 * np.ones( np.shape( np.unique( radar.elevation['data'] ) ) )
       if '0117_02' in filename  :  #122-2 STRATEGY
          nyquist_velocity     = 13.25 * np.ones( np.shape( np.unique( radar.elevation['data'] ) ) )

       if '0121_01' in filename  :  #122-1 STRATEGY
          nyquist_velocity     = 6.63 * np.ones( np.shape( np.unique( radar.elevation['data'] ) ) )
       if '0121_02' in filename  :  #122-2 STRATEGY
          nyquist_velocity     = 13.25 * np.ones( np.shape( np.unique( radar.elevation['data'] ) ) )

       if '0122_01' in filename  :  #122-1 STRATEGY
          nyquist_velocity     = 8.28 * np.ones( np.shape( np.unique( radar.elevation['data'] ) ) )
       if '0122_02' in filename  :  #122-2 STRATEGY
          nyquist_velocity     = 39.79 * np.ones( np.shape( np.unique( radar.elevation['data'] ) ) )
       if '0122_03' in filename  :  #122-3 STRATEGY
          nyquist_velocity     = 13.35 * np.ones( np.shape( np.unique( radar.elevation['data'] ) ) )

       if '0123_01' in filename  :  #123-1 STRATEGY
          nyquist_velocity     = 8.28 * np.ones( np.shape( np.unique( radar.elevation['data'] ) ) )
       if '0123_02' in filename  :  #123-2 STRATEGY
          nyquist_velocity     = 39.79 * np.ones( np.shape( np.unique( radar.elevation['data'] ) ) )
       if '0123_03' in filename  :  #123-3 STRATEGY
          nyquist_velocity     = 13.25 * np.ones( np.shape( np.unique( radar.elevation['data'] ) ) )
       if '0123_04' in filename  :  #123-4 STRATEGY
          nyquist_velocity     = 8.28 * np.ones( np.shape( np.unique( radar.elevation['data'] ) ) )

       if '0200_01' in filename  :  #200-1 STRATEGY
          nyquist_velocity     = 4.42 * np.ones( np.shape( np.unique( radar.elevation['data'] ) ) )
       if '0200_02' in filename  :  #200-2 STRATEGY
          nyquist_velocity     = 13.25 * np.ones( np.shape( np.unique( radar.elevation['data'] ) ) )

       if '0300_01' in filename  :  #300-1 STRATEGY
          nyquist_velocity     = 4.42 * np.ones( np.shape( np.unique( radar.elevation['data'] ) ) )
       if '0300_02' in filename  :  #300-2 STRATEGY
          nyquist_velocity     = 16.56 * np.ones( np.shape( np.unique( radar.elevation['data'] ) ) )

       if '0301_01' in filename  :  #300-1 STRATEGY

       #Esta estrategia tiene una velocidad nyquist que varia con el angulo de elevacion.
          nyquist_velocity     = np.array([4.42,4.42,4.42,6.63,6.63,6.63,6.63,8.28,8.28,8.28,8.28,8.28,8.28,8.28,8.28])
       if '0301_02' in filename  :  #300-2 STRATEGY
          nyquist_velocity     = 13.25 * np.ones( np.shape( np.unique( radar.elevation['data'] ) ) )


       if '0201_01' in filename  :  #201-1 STRATEGY
          nyquist_velocity     = 4.42 * np.ones( np.shape( np.unique( radar.elevation['data'] ) ) )
       if '0201_02' in filename  :  #201-2 STRATEGY
          nyquist_velocity     = 13.25 * np.ones( np.shape( np.unique( radar.elevation['data'] ) ) )
       if '0201_03' in filename  :  #201-3 STRATEGY
          nyquist_velocity     = 8.28 * np.ones( np.shape( np.unique( radar.elevation['data'] ) ) )

       if '0202_01' in filename  :  #200-1 STRATEGY
          nyquist_velocity     = 4.42 * np.ones( np.shape( np.unique( radar.elevation['data'] ) ) )
       if '0202_02' in filename  :  #200-2 STRATEGY
          nyquist_velocity     = 13.25 * np.ones( np.shape( np.unique( radar.elevation['data'] ) ) )
       


    if ( 'PAR' in filename ) or ( 'ANG' in filename ) or ( 'PER' in filename )  :
          if np.max( radar.range['data']  ) == 119875.0 :
             nyquist_velocity = 39.8 * np.ones( np.shape( np.unique( radar.elevation['data'] ) ) )  #120
          if np.max( radar.range['data']  ) == 239750.0 :
             nyquist_velocity = 6.63 * np.ones( np.shape( np.unique( radar.elevation['data'] ) ) ) #240



    #Correct instrument altitude.
    if 'PAR' in filename   :

       radar.altitude_agl['data']= np.array( 30.0)
       radar.altitude['data']    = np.array(122.0)

    if 'ANG' in filename   :

       radar.altitude_agl['data']= np.array( 30.0)
       radar.altitude['data']    = np.array(190.0)

    if 'PER' in filename   :

       radar.altitude_agl['data']= np.array( 30.0)
       radar.altitude['data']    = np.array(100.0)

    if 'RMA1' in filename  :

       radar.altitude_agl['data']= np.array( 35.0)
       radar.altitude['data']    = np.array(476.0)

    if 'RMA2' in filename   :

       radar.altitude_agl['data']= np.array( 35.0)
       radar.altitude['data']    = np.array( 47.0)

    if 'RMA3' in filename  :

       radar.altitude_agl['data']= np.array( 35.0)
       radar.altitude['data']    = np.array(197.0)

    if 'RMA4' in filename   :

       radar.altitude_agl['data']= np.array( 35.0)
       radar.altitude['data']    = np.array(119.0)

    if 'RMA5' in filename   :

       radar.altitude_agl['data']= np.array( 35.0)
       radar.altitude['data']    = np.array(841.0)

    if 'RMA6' in filename   :

       radar.altitude_agl['data']= np.array( 35.0)
       radar.altitude['data']    = np.array( 80.0)

    if 'RMA7' in filename   :

       radar.altitude_agl['data']= np.array( 35.0)
       radar.altitude['data']    = np.array(388.0)

    if 'RMA8' in filename   :

       radar.altitude_agl['data']= np.array( 35.0)
       radar.altitude['data']    = np.array(131.0)

    #Some common parameters
    ray_angle_res        = 1.0
    radar_beam_width_h   = 1.0
    radar_beam_width_v   = 1.0
    meters_between_gates  = radar.range['data'][1]-radar.range['data'][0]

    #Apply the missing parameters to the radar structure.

    if (not  'nyquist_velocity' in radar.instrument_parameters ) or ( radar.instrument_parameters['nyquist_velocity'] == None )  :

       radar.instrument_parameters['nyquist_velocity']=dict()
       radar.instrument_parameters['nyquist_velocity']['long_name']='unambiguous_doppler_velocity'
       radar.instrument_parameters['nyquist_velocity']['units']='meters per second'
       radar.instrument_parameters['nyquist_velocity']['_FillValue']= local_fill_value
       radar.instrument_parameters['nyquist_velocity']['meta_group']='instrument_parameters'

       tmp_nyquist = np.ones( np.shape( radar.azimuth['data'] ) )
       for ielev , elev in enumerate( np.unique( radar.elevation['data'] ) ) :

          tmp_nyquist[ radar.elevation['data'] == elev ] = nyquist_velocity[ ielev ]

       radar.instrument_parameters['nyquist_velocity']['data'] = ma.array( tmp_nyquist , mask = np.zeros( np.shape( tmp_nyquist ) , dtype=bool ) , fill_value = local_fill_value )

    if (not 'radar_beam_width_v' in radar.instrument_parameters ) or ( radar.instrument_parameters['radar_beam_width_v'] == None ) :
       radar.instrument_parameters['radar_beam_width_v']=dict()
       radar.instrument_parameters['radar_beam_width_v']['long_name']='half_power_radar_beam_width_v_channel'
       radar.instrument_parameters['radar_beam_width_v']['units']='degrees'
       radar.instrument_parameters['radar_beam_width_v']['_FillValue']= local_fill_value
       radar.instrument_parameters['radar_beam_width_v']['meta_group']='instrument_parameters'

       radar.instrument_parameters['radar_beam_width_v']['data'] = ma.array( radar_beam_width_v , mask =False , fill_value = local_fill_value )

    if (not 'radar_beam_width_h' in radar.instrument_parameters ) or ( radar.instrument_parameters['radar_beam_width_h'] == None ) :
       radar.instrument_parameters['radar_beam_width_h']=dict()
       radar.instrument_parameters['radar_beam_width_h']['long_name']='half_power_radar_beam_width_h_channel'
       radar.instrument_parameters['radar_beam_width_h']['units']='degrees'
       radar.instrument_parameters['radar_beam_width_h']['_FillValue']= local_fill_value
       radar.instrument_parameters['radar_beam_width_h']['meta_group']='instrument_parameters'
 
       radar.instrument_parameters['radar_beam_width_h']['data'] = ma.array( radar_beam_width_h , mask =False , fill_value = local_fill_value )

    if radar.ray_angle_res == None :
       radar.ray_angle_res = dict()
       radar.ray_angle_res['long_name']='angular_resolution_between_rays'
       radar.ray_angle_res['units']='degrees'
       radar.ray_angle_res['_FillValue']= local_fill_value
    
       radar.ray_angle_res['data'] = ma.array( np.ones( np.shape( levels ) )*ray_angle_res , mask = np.zeros( np.shape( levels ) , dtype=bool ) , fill_value = local_fill_value )

    if radar.range == None :
       radar.range = dict()

    if not 'meters_between_gates' in radar.range :
       radar.range['meters_between_gates']= meters_between_gates

    if not 'meters_to_center_of_first_gate' in radar.range :
       radar.range['meters_to_center_of_first_gate']= radar.range['data'][0]  #meters_between_gates / 2.0

    return radar


def get_file_list( datapath , init_time , end_time , time_search_type = None , file_type_list = None , instrument_type_list = None )     :

   #datapath : base path of radar data
   #init time: [yyyymmddhhMMss] beginning of the time window
   #end time : [yyyymmddhhMMss] end of the time window
   #time_search_type : [filename] or [timestamp]
   #file_types_list  : a list with file extensions that will be included in the file_list

   #import os
   #import datetime as dt 
   #import numpy as np

   if time_search_type == None :
      time_search_type = 'timestamp'

   
   date_min = dt.datetime.strptime( init_time , '%Y%m%d%H%M%S')
   date_max = dt.datetime.strptime( end_time  , '%Y%m%d%H%M%S')

   file_list=[]

   for (dirpath, dirnames, filenames) in os.walk( datapath ):

      for filename in filenames            :
         current_filename = '/'.join([dirpath,filename])

         if time_search_type == 'filename'   :
            date_c = get_time_from_filename( current_filename )
         if time_search_type == 'timestamp'  :
            date_c = dt.fromtimestamp( os.stat(current_filename).st_ctime )
         if date_c != None  :
            if date_c >= date_min and date_c <= date_max  :
               file_list.append( current_filename )
   
   #Keep only some file names and some paths.

   tmp_file_list = []

   if file_type_list != None :

      for my_file in file_list  :
     
         filename = os.path.basename( my_file )

         if any(ft in filename for ft in file_type_list ):
 
            tmp_file_list.append( my_file )

      file_list = tmp_file_list[:]

   tmp_file_list = []

   if instrument_type_list != None :
      
      for my_file in file_list :

         instrument_type = get_instrument_type_from_filename( my_file )

         for it in instrument_type_list :
             if it == instrument_type  : 
                tmp_file_list.append( my_file ) 

      file_list = tmp_file_list[:]


   

   return file_list



def get_time_from_filename( file_complete_path )    :

   #import datetime as dt
   #import os 

   filename = os.path.basename( file_complete_path )
   file_time = None

   format = get_format_from_filename( file_complete_path )

   if format == 'h5'    :

      file_time  = dt.datetime.strptime(filename.split('_')[-1][:15], '%Y%m%dT%H%M%S')  

   if format == 'vol'   :

      file_time  = dt.datetime.strptime(filename[:14], '%Y%m%d%H%M%S')

   if format == 'cfrad'   :
      
      file_time  = dt.datetime.strptime( filename.split('.')[1] , '%Y%m%d_%H%M%S')

   if format == 'letkf'   :
 
      file_time  = dt.datetime.strptime(filename.split('_')[-1][:14], '%Y%m%d%H%M%S')

   if format == 'pickle'  :

      file_time  = dt.datetime.strptime(filename.split('_')[-1][:14], '%Y%m%d%H%M%S')

   if format == 'letkf'   :

      file_time  = dt.datetime.strptime(filename.split('_')[-1][:14], '%Y%m%d%H%M%S')

   if format == 'tgz'     :

      file_time  = dt.datetime.strptime(filename[:10], '%Y%m%d_%H')

   return file_time


def get_instrument_type_from_filename( file_complete_path ) :
   #import os

   filename = os.path.basename( file_complete_path )
   instrument_name = None
   if 'RMA' in filename    :
      index = filename.find('RMA')
      instrument_name = filename[index:index+4]
   if 'ANG' in file_complete_path    :
      instrument_name = 'ANG'
   if 'PAR' in file_complete_path    :
      instrument_name = 'PAR'
   if 'PER' in file_complete_path    :
      instrument_name = 'PER'

   return instrument_name    



def merge_radar_object( radar_1 , radar_2 )    :

   #import numpy as np

   na_1 = np.shape( radar_1.azimuth['data'] )[0]
   na_2 = np.shape( radar_2.azimuth['data'] )[0]

   nr_1 = np.shape( radar_1.range['data'] )[0]
   nr_2 = np.shape( radar_2.range['data'] )[0]


   azimuth_1 = radar_1.azimuth['data'] 
   azimuth_2 = radar_2.azimuth['data']

   range_1   = radar_1.range['data']
   range_2   = radar_2.range['data']

   elev_1    = radar_1.elevation['data']
   elev_2    = radar_2.elevation['data']

   merged = False  #Wether the two radars has been successfully merged.

   if na_1 == na_2  :
      diff_a = np.sum( azimuth_1 - azimuth_2 )
   else             :
      diff_a = 0.0
   if nr_1 == nr_2  :
      diff_r = np.sum( range_1   - range_2   )
   else             :
      diff_r = 0.0

   #Check the ideal case
   if ( na_1 == na_2 ) and ( nr_1 == nr_2 ) and ( diff_a == 0.0 ) and ( diff_r == 0.0 )  :
      #Dimensions of radar_1 and radar_2 are the same.
      merged = True 
      #Check if the field is not already present. If it is present do not merge.
      for my_key in radar_2.fields   :
         if not my_key in radar_1.fields   :
            radar_1.fields[my_key] = radar_2.fields[my_key]
   else                                                                                  :
      print('Warning: Inconsistent shapes',na_1,na_2,nr_1,nr_2,diff_a,diff_r)
      #Test if radial shapes conform.
      if  nr_1 != nr_2  :
         print('Different sizes in radial dimension ',nr_1,' ',nr_2)
         #These objects have different ranges. We will try to solve this issue.
         if nr_1 < nr_2 :
            small_range = range_1
            big_range   = range_2
            small_radar = radar_1
            big_radar   = radar_2
            small_na    = na_1
            big_na      = na_2
            small_nr    = nr_1
            big_nr      = nr_2
         else           :
            small_range = range_2
            big_range   = range_1
            big_radar   = radar_1
            small_radar = radar_2
            small_na    = na_2
            big_na      = na_1
            small_nr    = nr_2
            big_nr      = nr_1

         if np.shape( np.intersect1d( small_range , big_range ) )[0] == np.shape( small_range )[0]  :
            #We take the filds corresponding to the small_radar and we extend them so they have the same
            #number of ranges as the large_radar.
            for my_key in small_radar.fields    :
                if not my_key in big_radar.fields   :
                    tmp_field = small_radar.fields[ my_key ] 
                    undef = tmp_field['_FillValue']
                    tmp_data = np.ones( ( small_na , big_nr ) ) * undef 
                    tmp_data[ : , 0:small_nr ] = tmp_field['data'].data
                    small_radar.fields[ my_key ]['data'] = np.ma.masked_array( tmp_data , tmp_data == undef )
                    small_nr = big_nr
         else                                                                                        :
             merged = False
             return radar_1 , merged 

      #At this point nr is the same for both radars. Before combining the radars we need to check if 
      #azimuths can be combined. 
      if ( na_1 == na_2 ) and ( diff_a == 0.0 )   :
         for my_key in radar_2.fields    :
            if not my_key in radar_1.fields   :
               radar_1.fields[my_key] = radar_2.fields[my_key]
               merged = True

      if ( na_1 != na_2 ) or ( diff_a != 0.0 )  :
         print('Warning: Inconsistency in azimuth ',na_1,' ',na_2,' ',diff_a)
         azimuth_1 = np.round( azimuth_1 * 10.0 )/10.0 
         azimuth_2 = np.round( azimuth_2 * 10.0 )/10.0
         ##print( np.size( np.intersect1d( azimuth_1 , azimuth_2 ) ) , na_1  , np.max( azimuth_1[1:2888]-azimuth_2[0:2887]) )
         ##Azimuths differ. We will take azimuths_1 as a reference and try to make them compatible.
         ##if ( np.shape( np.intersect1d( azimuth_1 , azimuth_2 ) )[0] / na_1 ) >= 0.95    :
         ##Both objects are similar we will try to merge them.
         for my_key in radar_2.fields    :
            if not my_key in radar_1.fields   :
                undef = radar_2.fields[ my_key ][ '_FillValue' ] 
                tmp_data = np.ones( ( na_1 , nr_1 ) ) * undef 
                for i1 in range( 0 , na_1 )   :

                     my_ind=np.nonzero( np.logical_and( azimuth_2 == azimuth_1[i1] , elev_2 == elev_1[i1] ) )
                     if np.size( my_ind ) == 1 :
                        tmp_data[i1,:] = radar_2.fields[ my_key ]['data'][my_ind,:]
                radar_1.fields[my_key] = dict()
                radar_1.fields[my_key] = radar_2.fields[my_key] 
                radar_1.fields[my_key]['data'] = np.ma.masked_array( tmp_data , tmp_data == undef )
         merged = True

   return radar_1  , merged

def upload_to_ftp( filename_list , ftp_host, ftp_user, ftp_pass , ftp_path , ftp_passive=False , compress=False ) :
    from ftplib import FTP
    #import os
    
    ftp = FTP(ftp_host, ftp_user, ftp_pass)
    ftp.cwd(ftp_path)
    ftp.set_pasv(ftp_passive)
    for my_file in filename_list  :
       my_path=os.path.dirname( my_file )
       my_file_name=os.path.basename( my_file )
       my_current_path=os.getcwd()
       os.chdir(my_path)

       if compress  :
          os.system('gzip -f ' + my_file )
          my_file_name = my_file_name + '.gz'
       ftp.storbinary('STOR ' + my_file_name , open( my_file_name ,'rb') )
          
       os.chdir(my_current_path)

def remove_from_ftp_timebased( ftp_host, ftp_user, ftp_pass , ftp_path , ini_time , end_time , file_format_list = None ) :
    from ftplib import FTP
    #import os
    #import datetime as dt

    if file_format_list == None :
       file_format_list = []

    ftp = FTP(ftp_host, ftp_user, ftp_pass)
    ftp.cwd(ftp_path)

    #Get the file list in the remote server at the current directory
    files = []

    try:
       files = ftp.nlst()
    except :
       print ('No file list available')

    #Keep only the files which are within the time period and remove the rest.
    date_min = dt.datetime.strptime( ini_time  , '%Y%m%d%H%M%S')
    date_max = dt.datetime.strptime( end_time  , '%Y%m%d%H%M%S')

    for my_file in files   :

       date_c = get_time_from_filename( my_file )
       file_format = get_format_from_filename( my_file )

       if date_c != None  :
          if date_c <= date_min or date_c >= date_max  :
             if file_format in file_format_list  :
                 print('Deleting ' + my_file + ' from remote server' )
                 ftp.delete( my_file )
          
   
def remove_from_localpath_timebased( local_path , ini_time , end_time , file_format_list = None , time_search_type = None )  :

   if time_search_type == None :
      time_search_type = 'filename'

   if file_format_list == None :
      file_format_list = []

   date_min = dt.datetime.strptime( ini_time  , '%Y%m%d%H%M%S')
   date_max = dt.datetime.strptime( end_time  , '%Y%m%d%H%M%S')

   file_list=[]

   for (dirpath, dirnames, filenames) in os.walk( local_path ):

      for filename in filenames            :
         current_filename = '/'.join([dirpath,filename])
 
         file_format = get_format_from_filename( filename )

         if time_search_type == 'filename'   :
            date_c = get_time_from_filename( current_filename )
         if time_search_type == 'timestamp'  :
            date_c = dt.fromtimestamp( os.stat(current_filename).st_ctime )
         if date_c != None  :
            if date_c <= date_min or date_c >= date_max  :
               if file_format in file_format_list :
                  #We will remove this file.
                  print('Deleting ' + current_filename )
                  os.system('rm ' + current_filename )


def save_cfradial( local_path , radar , fileformat='NETCDF4' )  :

    import pyart

    os.makedirs( local_path ,exist_ok=True)

    input_file = radar.files[0] 

    file_instrument = get_instrument_type_from_filename( input_file )

    file_time   = dt.datetime.strftime( get_time_from_filename( input_file ) , '%Y%m%d_%H%M%S' )
 
    filename = local_path + '/cfrad.' + file_time + '.' + file_instrument + '.nc'

    print('Writing file : ',filename)

    pyart.io.cfradial.write_cfradial(filename , radar, format=fileformat  ) #,  arm_time_variables=False)

  




