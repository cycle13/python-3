#print __doc__
# Author: Rapid Refresh Argentina Team
# License: BSD 3 clause

def main_qc( options , radar = None ) :

   import sys
   import time
   sys.path.append('../fortran')

   import numpy as np
   import numpy.ma as ma
   import pyart
   from common_qc_tools  import qc  #Fortran code routines.
   import os


   start=time.time()

   output=dict() #Initialize output dictionary.

   computed_etfilter=False   #Flag to indicate if echo top has been computed already.

   #Read the data

   if radar == None  :

      print('')
      print('-------------------------------------------')
      print('Reading the data')
      print('-------------------------------------------')
      print('')

      radar = pyart.io.read(options['filename'])

   output['maxw_ref']=0.0
   output['maxw_v']  =0.0

   #===================================================
   # CHECK IF WE HAVE ENOUGH VARIABLES TO PROCEED.
   #===================================================

   if ( not options['name_v'] in radar.fields) and ( not options['name_ref'] in radar.fields) :

      print('')
      print('-------------------------------------------')
      print('ERROR: Not enough variables to proceed with QC for this file')
      print('-------------------------------------------')
      print('')
      return radar , output

   #===================================================
   # CHECK IF WE HAVE INFINITE OR NAN VALUES IN DATA 
   #===================================================

   #Check different aspects of the input data, numerical preccission
   #presence of invalid data like Inf or NaN.

   print('')
   print('-------------------------------------------')
   print('Checking if input data is valid')
   print('-------------------------------------------')
   print('')

   radar = check_valid_data( radar , options )

   #===================================================
   # SOME PRELIMINARY PROCEDURES.
   #===================================================

   #From cfradial format to azimuth,range,elevation order

   print('')
   print('-------------------------------------------')
   print('Reshaping the variables')
   print('-------------------------------------------')
   print('')

   [ radar , output ] = reshape_variables( radar , output , options )
   #Compute X,Y,lat and lon for each pixel

   print('')
   print('-------------------------------------------')
   print('Georeferencing the data')
   print('-------------------------------------------')
   print('')

   [ radar , output ] = georeference_data( radar , output , options )

   print('')
   print('-------------------------------------------')
   print('Interpolating the terrain')
   print('-------------------------------------------')
   print('')

   #Compute the terrain height corresponding to each pixel
   [ radar , output ] = get_topography( radar , output , options )

   #===================================================
   # FILTER APPLICATION
   #===================================================

   #This functions runs the selected filters in the selected
   #order.

   print('')
   print('-------------------------------------------')
   print('Quality control routines')
   print('-------------------------------------------')
   print('')

   [ radar , output ] = run_filters( radar , output , options )

   print('')
   print('-------------------------------------------')
   print('Updating radar object')
   print('-------------------------------------------')
   print('')


   update_radar_object( radar , output , options )

   if options['output_to_file']   :

      print('')
      print('-------------------------------------------')
      print('Writing the data')
      print('-------------------------------------------')
      print('')

      pyart.io.cfradial.write_cfradial(options['filename_out'],radar,format=options['file_out_format'],time_reference=True, arm_time_variables=True)


   print('')
   print('-------------------------------------------')
   print('End of QC routine')
   print('-------------------------------------------')
   print('')

   end=time.time()

   print("The elapsed time in {:s} is {:2f}".format("the entire QC",end-start) )


   return radar , output


#===================================================
# FUNCIONES GENERALES ------------------------------
#===================================================

#===================================================
# Check if input data is valid.
#===================================================

def check_valid_data( radar , options )   :

    import numpy as np
    import time

    start=time.time()

    #=========================
    # Check numeric precission
    #=========================

    #Sometimes radar filds have precissions other than numpy default. This breakes the code since the code relies on
    #default numpy dtype.
    default_type = np.ones( (1) ).dtype
    for my_key in radar.fields  :
       radar.fields[my_key]['data'] = ( radar.fields[my_key]['data'] ).astype( default_type )

    #=========================
    # Check if variables are stored as masked arrays 
    # If not then convert them to masked arrays.
    #=========================

    for my_key in radar.fields    :

       if not np.ma.is_masked( radar.fields[ my_key ]['data'] ) :

          radar.fields[ my_key ]['data'] = np.ma.masked_array( radar.fields[ my_key ]['data'] , mask = radar.fields[ my_key ]['data'] == ['_FillValue'] )


    #=========================
    # Check Nan and Inf
    #=========================

    for my_key in radar.fields    :
       print( my_key  )
       if my_key == options['name_ref']   :

          if np.sum( radar.fields[ my_key ]['data'].data != radar.fields[ my_key ]['_FillValue'] ) > 0.005*np.size( radar.fields[ my_key ]['data'].data )   :
             
             radar.fields[ my_key ]['data'].data[ np.isinf( radar.fields[ my_key ]['data'].data )  ] = options['norainrefval']
             radar.fields[ my_key ]['data'].data[ np.isnan( radar.fields[ my_key ]['data'].data )  ] = options['norainrefval']
             radar.fields[ my_key ]['data'].data[ radar.fields[ my_key ]['data'].data < options['norainrefval']  ] = options['norainrefval']
             radar.fields[ my_key ]['data'].data[ radar.fields[ my_key ]['data'].data == radar.fields[ my_key ]['_FillValue']  ] = options['norainrefval']
             radar.fields[ my_key ]['data'].mask= radar.fields[ my_key ]['data'].data == radar.fields[ my_key ]['_FillValue']

       else                               :
          radar.fields[ my_key ]['data'].data[ np.isinf( radar.fields[ my_key ]['data'].data )  ] = radar.fields[ my_key ]['_FillValue']
          radar.fields[ my_key ]['data'].data[ np.isnan( radar.fields[ my_key ]['data'].data )  ] = radar.fields[ my_key ]['_FillValue']
          radar.fields[ my_key ]['data'].mask= radar.fields[ my_key ]['data'].data == radar.fields[ my_key ]['_FillValue']

    end=time.time()
    print("The elapsed time in {:s} is {:2f}".format("Input data check",end-start) )

    return radar

#==================================================a
#' RESHAPE VARIABLES
#===================================================

def reshape_variables( radar , output , options )    :

   import numpy as np
   import time
   #From time,range -> azimuth,range,elevation

   start=time.time()

   if options['name_ref'] in radar.fields :

        output['undef_ref']=radar.fields[options['name_ref']]['_FillValue']

        [ output['ref'] , output['az'] , output['level'] , output['time'] , output['az_exact'] ]=order_variable( radar , options['name_ref'] , output['undef_ref'] )
        na=output['ref'].shape[0]
        nr=output['ref'].shape[1]
        ne=output['ref'].shape[2]
        output['na']=na
        output['nr']=nr
        output['ne']=ne

        output['cref'] = np.ones( (na,nr,ne))
        output['cref'] = np.copy(output['ref'])           #Initialize the corrected reflectivity array.
        output['ref_input'] = np.copy(output['ref'])      #Initialize the input reflectivity array (for ploting only)
 
        output['cref'][ output['cref'] == output['undef_ref'] ]=options['norainrefval']
        output['ref'] [ output['ref']  == output['undef_ref'] ]=options['norainrefval']

        output['qcref'] = np.zeros(output['cref'].shape)  #Set the qc flag array to 0.
        output['wref']  = np.zeros(output['cref'].shape)  #Set the weigths to 0.

        output['elevations']=np.unique(radar.elevation['data'])


   if options['name_v'] in radar.fields  : 

        output['undef_v']=radar.fields[ options['name_v'] ]['_FillValue']

        [ output['v'] , output['az'] , output['level'] , output['time'] , output['az_exact']  ]=order_variable( radar , options['name_v'] , output['undef_v']  )
 
        na=output['v'].shape[0]
        nr=output['v'].shape[1]
        ne=output['v'].shape[2]
        output['na']=na
        output['nr']=nr
        output['ne']=ne

        output['cv'] = np.ones((na,nr,ne))
        output['cv'] = np.copy(output['v'])              #Initialize the corrected velocity array
        output['input_v'] = np.copy(output['v'])         #Initialize the input velocity array (for ploting only)

        output['qcv'] = np.zeros(output['v'].shape)      #Set the qc flag array to 0.

        output['wv'] = np.zeros(output['v'].shape)       #Set the weigths to 0.

        output['elevations']=np.unique(radar.elevation['data'])

   end=time.time()
   print("The elapsed time in {:s} is {:2f}".format("Data reshape",end-start) )

   return radar , output

#===================================================
# GEOREFERENCE RADAR DATA
#===================================================

def georeference_data( radar , output , options )    :

   import numpy as np
   import time
   #Use pyart rutines to compute x, y and z at each grid point
   #Reorder data
   start=time.time()

   #dm is a dummy variable

   [ output['altitude'] , dm , dm , dm , dm ] = order_variable( radar , 'altitude' ,  options['undef'] )
   [ output['x']        , dm , dm , dm , dm ] = order_variable( radar , 'x' ,  options['undef'] ) 
   [ output['y']        , dm , dm , dm , dm ] = order_variable( radar , 'y' ,  options['undef'] )


   #Compute distance to radar for the first azimuth (we assume that distance to radar will be the same.
   #for all other azimuths.
   output['distance']=np.power( np.power(output['x'][0,:,:],2)+np.power(output['y'][0,:,:],2) , 0.5 )

   end=time.time()

   print("The elapsed time in {:s} is {:2f}".format("Data georeferentiation",end-start) )

   return radar , output 

#===================================================
# READ AND INTERPOLATE TOPOGRAPHY DATA
#===================================================

def get_topography( radar , output , options )  :

   import time
   import numpy as np
   import os
   from scipy.interpolate import interp2d

   na=output['na']
   nr=output['nr']
   ne=output['ne']

   start=time.time()

   #The idea is to store the topography interpolated to the earth surface in a polar geometry 
   #sourrounding each radar into a binary file. If this binary file is already created we read it
   #and interpolate the topography to the 3D radar geometry. If not we read the raw raster data.

   #Name of the file containing the topography interpolated the earth surface surrounding the radar.
   polar_coord_topo_file=options['toporadardatapath'] + '/' + radar.metadata['instrument_name'] + '_' +  \
             str( radar.range['meters_between_gates'] ) + '_' + str( np.size( radar.range['data'] ) * radar.range['meters_between_gates'] / 1000 ) \
             + '.dat'
   raster_topo_file=options['toporadardatapath'] + '/' + radar.metadata['instrument_name'] + '_' +  \
             str( radar.range['meters_between_gates'] ) + '_' + str( np.size( radar.range['data'] ) * radar.range['meters_between_gates'] / 1000 ) \
             + '.tif'

   try    :
      print('Trying to get a binary file with the interpolated topography')
      my_topo=read_topo( polar_coord_topo_file )

   except :
      print('I could not find an adequate binary file. We will get the data')
      os.makedirs(options['toporawdatapath'],exist_ok=True)
      os.makedirs(options['toporadardatapath'],exist_ok=True)

      #Generate topo file ( inputs are , radar lat and lon, range and azimuth )
      my_topo=generate_topo_file( radar.longitude['data'][0] , radar.latitude['data'][0] , radar.range['data'] , output['az'] ,
                                  options['toporawdatapath'] , polar_coord_topo_file )
                        

   #if ( os.path.isfile( polar_coord_topo_file ) )  :
   #   read_raster = False
   #   print('Using a previously generated topography file')
   #else                                :
   #   read_raster = True
   #   print('Topography file not found. We will generate a new file from raw raster data')

   #if read_raster    :   #We read the original data and interpolate it to a polar grid centered at the radar.

   #   os.makedirs(options['toporawdatapath'],exist_ok=True)
   #   os.makedirs(options['toporadardatapath'],exist_ok=True)

   #   #Generate topo file ( inputs are , radar lat and lon, range and azimuth )
   #   my_topo=generate_topo_file( radar.longitude['data'][0] , radar.latitude['data'][0] , radar.range['data'] , output['az'] , 
   #                               options['toporawdatapath'] , polar_coord_topo_file )

   #else   :

   #   print('Reading polar coordinate topography from a file')

   #   my_topo=read_topo( polar_coord_topo_file )

   #Interpolate topography data to the 3D radar structure.

   topo_r = my_topo['range'][0,:]
   topo_a = my_topo['azimuth'][:,0]

   #Define a "interpolator"
   interpolator = interp2d(topo_r,topo_a,my_topo['mean'], kind='linear')

   output['topo'] = np.zeros((output['na'],output['nr'],output['ne']))

   #Loop over the vertical levels.
   for ie in range( 0 , ne )   :

      ppi_r = output['distance'][:,ie]
      ppi_a = output['az']

      output['topo'][:,:,ie] =  interpolator( ppi_r, ppi_a )

   end=time.time()

   #Compute AGL height for each pixel.
   output['altitude_agl']=output['altitude']-output['topo'] 

   print("The elapsed time in {:s} is {:2f}".format("topography interpolation",end-start) )

   return radar , output

#===================================================
# COMPUTE THE FINAL WEIGHT AND APPLY FUZZY LOGIC QC
#===================================================
 
def weight_computation( output , options )   :

   import numpy as np

   if output['maxw_ref'] == 0    :
      output['wref'][:]=0.0
   else                          :
      output['wref']=output['wref'] / output['maxw_ref']

   if output['maxw_v'] == 0      :
      output['wv'][:]=0.0
   else                          :
      output['wv']=output['wv'] / output['maxw_v']

   #All the grid points where the probability of non-meteorological echoes
   #is greather than 0.5 will be flagged as missing.
   output['cref'][ output['wref'] > options['w_tr'] ] = output['undef_ref']
   output['cv'][output['wv'] > options['w_tr'] ] = output['undef_v']

   return output

#===================================================
# ADD CORRECTED DATA TO RADAR OBJECT
#===================================================

def  update_radar_object( radar , output , options )   :

   import numpy as np

   if options['keep_original_fields']   : 
      #We will add variables to the file corresponding to the corrected REF and VRAD.

      if options['name_ref'] in radar.fields :

         tmp=order_variable_inv( radar , output['cref'] , output['undef_ref'] )
         radar.add_field_like( options['name_ref'] , options['name_cref'] , radar.fields[ options['name_ref'] ]['data'] , True)
         radar.fields[ options['name_cref'] ]['data']=np.ma.masked_array(tmp , mask = ( tmp==output['undef_ref'] ) )

         #print( np.min( radar.fields[options['name_cref']]['data'] ) ) 

      if options['name_v'] in radar.fields :

         tmp=order_variable_inv( radar , output['cv'] , output['undef_v'] )
         radar.add_field_like( options['name_v'] , options['name_cv'] , radar.fields[ options['name_v'] ]['data'] , True)
         radar.fields[ options['name_cv'] ]['data']=np.ma.masked_array(tmp , mask = ( tmp==output['undef_v'] ) )

   else                                 :
      #We will overwrite the original variable names in the radar structure.
      if options['name_ref'] in radar.fields   : 
         tmp=order_variable_inv( radar , output['cref'] , output['undef_ref'] )
         radar.fields[ options['name_ref'] ]['data']=np.ma.masked_array(tmp , mask = ( tmp==output['undef_ref'] ) )
      if options['name_v'] in radar.fields     :
         tmp=order_variable_inv( radar , output['cv'] , output['undef_v'] )
         radar.fields[ options['name_v'] ]['data']=np.ma.masked_array(tmp , mask = ( tmp==output['undef_v'] ) )

   return radar , output 


#====================================================
# UPDATE OUTPUT BASED ON FILTER DATA
#====================================================

def output_update( output , qc_index , options , filter_name )  :

   from common_qc_tools  import qc  #Fortran code routines.
   import numpy as np

   na=output['na']
   nr=output['nr']
   ne=output['ne']

   #Compute the corresponding weigth.

   weigth = qc.multiple_1d_interpolation( field=qc_index , nx=na , ny=nr , nz=ne
                                         , undef=options['undef'] , xx=options[filter_name]['ifx']
                                         , yy=options[filter_name]['ify'] , nxx=np.size(options[filter_name]['ifx']) )

   weigth[ weigth == options['undef'] ] = 0.0


   if ( 'ref' in options[filter_name]['var_update_list'] ) and ( 'ref' in output )  : 
      #Update reflectivity
      if options[filter_name]['sequential']   :
         #Preserve the data used as input for the filter for ploting pourposes.
         output['input_ref']=np.copy(output['ref'])

      if not options[filter_name]['force']   :
         output['wref']=output['wref'] + weigth * options[filter_name]['w']
         output['qcref'][ weigth > 0.5 ] = options[filter_name]['code']
         output['maxw_ref']=output['maxw_ref'] + options[filter_name]['w']
      else                                   :
         if options[filter_name]['fill_value']  == 'undef'       :
            tmp_mask = np.logical_and( weigth > options[filter_name]['force_value'] , output['cref'] != output['undef_ref'] )
            output['cref'][ tmp_mask ]=output['undef_ref']
            if options[filter_name]['sequential']   :
               tmp_mask = np.logical_and( weigth > options[filter_name]['force_value'] , output['ref'] != output['undef_ref'] )
               output['ref'][ tmp_mask ]=output['undef_ref']
         elif options[filter_name]['fill_value']  == 'min_ref'   :
            tmp_mask = np.logical_and( weigth > options[filter_name]['force_value'] , output['cref'] != output['undef_ref'] )
            output['cref'][ tmp_mask ]=options['norainrefval'] 
            if options[filter_name]['sequential']   :
               tmp_mask = np.logical_and( weigth > options[filter_name]['force_value'] , output['cref'] != output['undef_ref'] )
               output['ref'][ tmp_mask ]=options['norainrefval']
         else                                                    :
            tmp_mask = np.logical_and( weigth > options[filter_name]['force_value'] , output['cref'] != output['undef_ref'] )
            output['cref'][ tmp_mask ]=options[filter_name]['fill_value']
            if options[filter_name]['sequential']   :
               tmp_mask = np.logical_and( weigth > options[filter_name]['force_value'] , output['ref'] != output['undef_ref'] )
               output['ref'][ tmp_mask ]=options[filter_name]['fill_value']


         output['qcref'][ weigth > options[filter_name]['force_value'] ] = options[filter_name]['code']
   
   if ( 'v' in options[filter_name]['var_update_list'] ) and ( 'v' in output )   :
      #Update radial velocity
      if options[filter_name]['sequential']   :
         #Preserve the data used as input for the filter for ploting pourposes.
         output['input_v']=np.copy(output['v'])

      if not options[filter_name]['force']   :
         output['wv']=output['wv'] + weigth * options[filter_name]['w']
         output['qcv'][ weigth > 0.5 ] = options[filter_name]['code']
         output['maxw_v']=output['maxw_v'] + options[filter_name]['w']
      else                                   :
         if (options[filter_name]['fill_value']  == 'undef') or  (options[filter_name]['fill_value']  == 'min_ref')  :
            output['cv'][ weigth > options[filter_name]['force_value'] ]=output['undef_v']
            if options[filter_name]['sequential']   :
               output['v'][ weigth > options[filter_name]['force_value'] ]=output['undef_v']
         else                                                    :
            output['cv'][ weigth > options[filter_name]['force_value'] ]= options[filter_name]['fill_value']
            if options[filter_name]['sequential']   :
               output['v'][ weigth > options[filter_name]['force_value'] ]=options[filter_name]['fill_value']


         output['qcv'][ weigth > options[filter_name]['force_value'] ] = options[filter_name]['code']


   #If "save" is true, then store the qc_index and the weigth computed for this particular 
   #filter. This will be stored in a dictionary named after the filter and which is part of the 
   #output dictionary. 
   if  options[filter_name]['save']      :
      output[filter_name]=dict()
      output[filter_name]['qc_index']=qc_index
      output[filter_name]['weight']=weigth
      output[filter_name]['undef']=options['undef']

   #Plot filter diagnostics if this option is available.
   if  options['plot']['Enable']    :
      plot_filter( output , qc_index , weigth , options , filter_name )

   return output

#====================================================
# PLOT THE FILTER
#====================================================

def plot_filter( output , qc_index , weigth , options , filter_name )  :

   import os.path 
   import numpy as np
   import matplotlib.pyplot as plt

   if ( 'ref' in options[filter_name]['var_update_list'] ) and ( 'ref' in output )  :

    
    #tmp_ref=np.ma.masked_array( output['input_ref'] , np.logical_or( output['input_ref'] == output['undef_ref'] , output['input_ref'] == options['norainrefval'] ) )
    #tmp_cref=np.ma.masked_array( output['cref'] , np.logical_or( output['cref'] == output['undef_ref'] ,  output['input_ref'] == options['norainrefval'] ) )

    tmp_ref=np.ma.masked_array( output['input_ref'] , output['input_ref'] == output['undef_ref'] )
    tmp_cref=np.ma.masked_array( output['cref'] , output['cref'] == output['undef_ref'] )

    tmp_qc_index=np.ma.masked_array( qc_index , qc_index == options['undef'] )

    tmp_x = np.ma.masked_array( output['x'] , output['x'] == options['undef'] )
    tmp_y = np.ma.masked_array( output['y'] , output['y'] == options['undef'] )

    for ilev in options['plot']['Elevs']  :

       plt.figure(figsize=(8, 8))
       plt.subplot(2,2,1)

       plt.pcolor(tmp_x[:,:,ilev]/1e3,tmp_y[:,:,ilev]/1e3, tmp_ref[:,:,ilev],vmin=options['plot']['DbzMin'],vmax=options['plot']['DbzMax'],cmap=options['plot']['CmapDbz'])
       plt.title('Original Reflectivity')
       plt.colorbar()

       plt.subplot(2,2,2)
       plt.pcolor(tmp_x[:,:,ilev]/1e3,tmp_y[:,:,ilev]/1e3, tmp_cref[:,:,ilev],vmin=options['plot']['DbzMin'],vmax=options['plot']['DbzMax'],cmap=options['plot']['CmapDbz'])
       plt.title('Corrected Reflectivity')
       plt.colorbar()

       plt.subplot(2,2,3)
       plt.pcolor(tmp_x[:,:,ilev]/1e3,tmp_y[:,:,ilev]/1e3, ( output['qcref'][:,:,ilev]==options[filter_name]['code'] ).astype(float) )
       plt.title('Pixels corrected by ' + filter_name)
       plt.colorbar()

       plt.subplot(2,2,4)
       plt.pcolor(tmp_x[:,:,ilev]/1e3,tmp_y[:,:,ilev]/1e3, ( tmp_qc_index[:,:,ilev] ) )
       plt.title('QC Index')
       plt.colorbar()
                                                                                                    
       if options['plot']['Show']  :

           plt.show()

       figname_prefix= options['plot']['Path'] + os.path.basename( options['filename'] )       
       figname=figname_prefix + '_elev_' + str(ilev) + '_' + filter_name + '_ref_' + options['plot']['FigNameSufix']
       plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',
                   orientation='portrait', papertype=None, format=None,
                   transparent=False, bbox_inches=None, pad_inches=0.1,
                   frameon=None)
        
       plt.close()

   if ( 'v' in options[filter_name]['var_update_list'] ) and ( 'v' in output )  :

    tmp_v=np.ma.masked_array( output['input_v'] , output['input_v'] == output['undef_v'] )
    tmp_cv=np.ma.masked_array( output['cv'] , output['cv'] == output['undef_v'] )
    tmp_qc_index=np.ma.masked_array( qc_index , qc_index == options['undef'] )

    for ilev in options['plot']['Elevs']   :

       plt.figure(figsize=(8, 8))
       plt.subplot(2,2,1)

       plt.pcolor(output['x'][:,:,ilev]/1e3,output['y'][:,:,ilev]/1e3,tmp_v[:,:,ilev],vmin=options['plot']['VrMin'],vmax=options['plot']['VrMax'],cmap=options['plot']['CmapWind'])
       plt.title('Original Doppler Velocity')
       plt.colorbar()

       plt.subplot(2,2,2)
       plt.pcolor(output['x'][:,:,ilev]/1e3,output['y'][:,:,ilev]/1e3,tmp_cv[:,:,ilev],vmin=options['plot']['VrMin'],vmax=options['plot']['VrMax'],cmap=options['plot']['CmapWind'])
       plt.title('Corrected Doppler Velocity')
       plt.colorbar()

       plt.subplot(2,2,3)
       plt.pcolor(output['x'][:,:,ilev]/1e3,output['y'][:,:,ilev]/1e3, ( output['qcv'][:,:,ilev]==options[filter_name]['code'] ).astype(float) )
       plt.title('Pixels corrected by ' + filter_name)
       plt.colorbar()

       plt.subplot(2,2,4)
       plt.pcolor(output['x'][:,:,ilev]/1e3,output['y'][:,:,ilev]/1e3, tmp_qc_index[:,:,ilev] )
       plt.title('QC index')
       plt.colorbar()

       if options['plot']['Show']  :

           plt.show()

       figname_prefix= options['plot']['Path'] + os.path.basename( options['filename'] )
       figname= figname_prefix + '_elev_' + str(ilev) + '_' + filter_name + '_v_' + options['plot']['FigNameSufix']
       plt.savefig(figname, dpi=None, facecolor='w', edgecolor='w',
                   orientation='portrait', papertype=None, format=None,
                   transparent=False, bbox_inches=None, pad_inches=0.1,
                   frameon=None)

       plt.close()

    return

#====================================================
# RUN THE SELECTED FILTERS IN THE CORRESPONDING ORDER
#====================================================

def run_filters( radar , output , options )   :

   import time
   import numpy as np

   my_filters       = list()

   my_filters_order = list()

   #First get the list of enabled filters from the options structure.
   #Note that 'order' can include several instances so filters can be 
   #Applied several times to the same data set
  
   for key in options  :
       
       if type( options[key] ) == dict  :

          if 'flag' in options[key]  :
           
              if options[key]['flag']  :
                 
                 for ii in range(len(options[key]['order'])) :
 
                     #Filters can be applied several times to the same data.
                     my_filters.append( key )
                     my_filters_order.append( options[key]['order'][ii] )

   #Sort the filters according to the order associated with each filter.   
 
   my_sorted_filters=[x for _,x in sorted(zip(my_filters_order,my_filters))]
   print('-------------------------------------------')
   print('The following filters will be applied')
   print('In the following order')
   print('-------------------------------------------')
   print('')
   for ifilter in my_sorted_filters  :
       print( ifilter )

   print('')
   print('-------------------------------------------')
              
   print('')
   print('-------------------------------------------')
   print('Running the filters')
   print('-------------------------------------------')
   print('')

   #Run the filtesr in the selected order.
   #If ploting is enabled then plot diagnostics.

   for ifilter in my_sorted_filters :

       start=time.time()

       if globals().get(ifilter) is None  :
          #If I can not find the corresponding function raise an error.
          print("Error: No filter called {}".format(ifilter))
       else                                      :
          #Run the filter (call the function named after the filter name.
          print('')
          print('-------------------------------------------')
          print('Running filter ' + ifilter )
          print('-------------------------------------------')
          print('')
          exec("[ radar , output ] = " + ifilter + "( radar , output , options )")

          #if 'cref' in output.keys()  :
          #    print( np.min( output['cref'] ) , output['undef_ref'] )

          end=time.time()
          print('')
          print("The elapsed time in {:s} is {:2f}".format(ifilter,end-start) )
          print('')

   print('')
   print('-------------------------------------------')
   print('Finish running the filters')
   print('-------------------------------------------')
   print('')


   return radar , output 

#===================================================
# FUNCIONES PARA LA IMPLEMENTACION DE DIVERSOS 
# FILTROS QUE COMPONEN EL QC.
#===================================================
 
#===================================================
# DEALIASING 
#===================================================

def Dealiasing( radar , output , options )   :

   import numpy as np
   import pyart

   filter_name='Dealiasing'
   if options[filter_name]['flag']  and ( options['name_v'] in radar.fields ) :

      #Define a new instance of the radar strcture containing the wind (potentially affected by previous filters).
      radar.add_field_like( options['name_v'], options['name_cv'] , radar.fields[ options['name_v'] ]['data'] , True)
      tmp = order_variable_inv(  radar , output['v'] , output['undef_v'] )
      radar.fields[ options['name_cv'] ]['data']= np.ma.masked_array( tmp , mask = tmp==output['undef_v'] )

      #Uso una de las funciones de dealiasing de pyart con los parametros por defecto
#      winddealias=pyart.correct.region_dealias.dealias_region_based(radar,interval_splits=20,interval_limits=None, 
#                 skip_between_rays=0,skip_along_ray=0,centered=True,nyquist_vel=None,
#                 check_nyquist_uniform=True,gatefilter=None,rays_wrap_around=True,keep_original=True,set_limits=True,
#                 vel_field=options['name_cv'] ,corr_vel_field=None)

      winddealias=pyart.correct.region_dealias.dealias_region_based(radar,interval_splits=3,interval_limits=None, 
                 skip_between_rays=30,skip_along_ray=200,centered=True,nyquist_vel=None,
                 check_nyquist_uniform=True,gatefilter=None,rays_wrap_around=True,keep_original=True,set_limits=True,
                 vel_field=options['name_cv'] ,corr_vel_field=None)

      #Replace cv wind by dealiased winds.
      radar.fields[ options['name_cv'] ]['data'] = winddealias['data']

      #Re-order dealiased wind data.
      [ output['cv'] , output['az'] , output['level'] , output['time'] , output['az_exact']  ]=order_variable( radar , options['name_cv'] , output['undef_v']  )
    
      mask=np.logical_and( output['cv'] != output['v'] , output['cv'] != output['undef_v'] )
      output['qcv'][ mask ]=options[filter_name]['code']

      radar.fields.pop(options['name_cv'], None)

      #Plot the data (this filter do not call update_output so we need to plot it here).
      tmp_index  = output['cv'] - output['v']
      tmp_weigth = np.zeros( np.shape( tmp_index ) )
      if options['plot']['Enable'] :
         plot_filter( output , tmp_index , tmp_weigth , options , filter_name )   

      #Store the difference between original and dealiased velocity.
      output['Dealiasing']=dict()
      output['Dealiasing']['vdiff']=( output['cv'] - output['v'] )
      output['Dealiasing']['vda']=np.copy( output['cv'] )
      #output['Dealiasing']['v']=np.copy( output['v'] )


      if options[filter_name]['sequential']     :
         #Following filters will be computed using the corrected velocity.
         output['v'] = output['cv'] 
   
   return radar , output 


#===================================================
# LOW ELEVATION ANGLE REFLECTIVITY FILTER
#===================================================

def DopplerRefFilter( radar , output , options)  :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.
   #Remove doppler gates associated with an undef refelctivity or with reflectivity under a certain threshold.

   filter_name='DopplerRefFilter'

   na=output['na']
   nr=output['nr']
   ne=output['ne']

   if   options['name_v'] in radar.fields   :

     if options['name_ref'] in radar.fields   :

        nx=options[filter_name]['nx']
        ny=options[filter_name]['ny']
        nz=options[filter_name]['nz']

        output['smooth_ref']=qc.box_functions_2d(datain=output['ref'],na=na,nr=nr,ne=ne,undef=output['undef_ref']
                                                  ,boxx=nx,boxy=ny,boxz=nz,operation='MEAN',threshold=0.0)

        tmp_index=np.zeros([na,nr,ne])

        if options[filter_name]['filter_undef']  :
           tmp_index[ output['smooth_ref'] == output['undef_ref']  ] = 1.0

        tmp_index[ np.logical_and( output['smooth_ref'] <= options[filter_name]['threshold'] , output['smooth_ref'] != output['undef_ref'] )] = 1.0

     else                                      :

        tmp_index=np.ones([na,nr,ne])
     
     output = output_update( output , tmp_index , options , filter_name )

   return radar , output

#===================================================
# MODEL FILTER 
#===================================================

def ModelFilter( radar , output , options)   :

   import numpy as np

   #TODO
   #TODO
   print('This filter is not ready yet')
     
   return radar , output 

#===================================================
# DEALIASING BORDER FILTER
#===================================================

def DealiasingEdgeFilter( radar , output , options )    :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.

   na=output['na']
   nr=output['nr']
   ne=output['ne']

   filter_name='DealiasingEdgeFilter'  

   #Test if Dealiasing has been applied to this data.
   if 'Dealiasing' in output :

      nyquistv=np.max( radar.get_nyquist_vel(0,check_uniform=True) )

      #Use a simple edge detection routine to detect the borders between dealiased and non-dealiased regions.
      tmp_index =qc.doppler_edge_filter( vdiff=output['Dealiasing']['vdiff'] , v=output['Dealiasing']['vda'], nx=na,ny=nr,nz=ne,undef=output['undef_v'],
                                                                               nboxx=options[filter_name]['nx'],nboxy=options[filter_name]['ny'],
                                                                               nboxz=options[filter_name]['nz'],edge_tr=( nyquistv/3.0 ) )

      tmp_index[ tmp_index == output['undef_v'] ] = options['undef']

      output = output_update( output , tmp_index , options , filter_name ) 

   return radar , output 

#===================================================
# ECHO TOP FILTER  
#===================================================

def EchoTopFilter( radar , output , options )   :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.

   filter_name='EchoTopFilter'
   if options['name_ref'] in radar.fields  :

      na=output['na']
      nr=output['nr']
      ne=output['ne']
      nx=options[filter_name]['nx']
      ny=options[filter_name]['ny']
      nz=options[filter_name]['nz']
                                   
      tmp_max_z=np.zeros((na,nr,ne))

      for ii in range(0,output['ne'])     :       #Estimate the maximum radar data height assoicated with each gate.
         tmp_max_z[:,:,ii]=output['altitude'][:,:,output['ne']-1]


      if  options[filter_name]['fast_computation']          :

         [tmp_index,tmp_data_2d]=qc.echo_top_fast(reflectivity=output['ref'],heigth=output['altitude'][0,:,:]
                                                ,rrange=output['distance'],na=na,nr=nr,ne=ne
                                                ,undef=output['undef_ref'],nx=nx,ny=ny,nz=nz)

      else                                     : 
         [tmp_data_3d,tmp_data_2d]=qc.echo_top(reflectivity=output['ref'],heigth=output['altitude'][0,:,:]
                                                ,rrange=output['distance'],na=na,nr=nr,ne=ne
                                                ,undef=output['undef_ref'],nx=nx,ny=ny,nz=nz)
         tmp_index=tmp_data_3d[:,:,:,0]

      tmp_index[ tmp_index != output['undef_ref'] ] = tmp_index[ tmp_index != output['undef_ref'] ] - output['topo'][ tmp_index != output['undef_ref'] ]  #Compute the echo top height over the terrain.
        
      computed_etfilter = True  #In case we need any of the other variables computed in this routine.

      #JUAN: Comente esta linea para poder filtrar tambien los ecos de terreno cerca del radar.
      #tmp_index[ tmp_max_z < options[filter_name]['heigthtr'] ] = output['undef_ref']  #Do not consider this filter when the volume maximum heigth is below
                                                                  #the specified threshold (i.e. pixels close to the radar)
      #If the pixel is already a "norain" pixel then leave it as it is.                           
      tmp_index[ output['ref'] == options['norainrefval'] ] = options['undef']
      #If the pixel is already undef then leave it as it is.
      tmp_index[ tmp_index == output['undef_ref'] ] = options['undef']

      output = output_update( output , tmp_index , options , filter_name ) 
    
   return radar , output 

#===================================================
# ECHO DEPTH FILTER 
#===================================================

def EchoDepthFilter( radar , output , options )   :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.

   filter_name='EchoDepthFilter'
   if options['name_ref'] in radar.fields  :

     na=output['na']
     nr=output['nr']
     ne=output['ne']
     nx=options[filter_name]['nx']
     ny=options[filter_name]['ny']
     nz=options[filter_name]['nz']

     tmp_max_z=np.zeros((na,nr,ne))

     for ii in range(0,ne)     :       #Estimate the maximum radar data height assoicated with each gate.
        tmp_max_z[:,:,ii]=output['altitude'][:,:,output['ne']-1]

     [tmp_data_3d,tmp_data_2d]=qc.echo_top(reflectivity=output['ref'],heigth=output['altitude'][0,:,:]
                                                ,rrange=output['distance'],na=na,nr=nr,ne=ne
                                                ,undef=output['undef_ref'],nx=nx,ny=ny,nz=nz)

     tmp_index=tmp_data_3d[:,:,:,2]

     tmp_index[ tmp_index == output['undef_ref'] ] = options['undef']

     tmp_index[ tmp_max_z < options[filter_name]['heigthtr'] ] = 1.0e6  #Do not consider this filter when the volume maximum heigth is below
                                                                  #the specified threshold (i.e. pixels close to the radar)

     output = output_update( output , tmp_index , options , filter_name ) 


   return radar , output 

#===================================================
# LOW ELEVATION ANGLE REFLECTIVITY FILTER
#===================================================

def LowElevFilter( radar , output , options)  :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.
   #Remove echos which are present at low elevation angles but not at higher elevation 
   #angles. This will help to remove clutter produced by anomalous propagation. 
   #This can also help to remove second trip echoes which are usually observed only at low elevation
   #angles.
   #This can also eliminate distant convective cells that are only seen in low 
   #elevation angles near the radar range edge. However these cells are of little interest
   #for data assimilation since that information is usually not enough to adequatelly constrain
   #the evolution of the convective cells.

   filter_name='LowElevFilter'

   if  options['name_ref'] in radar.fields :

     na=output['na']
     nr=output['nr']
     ne=output['ne']
     nx=options[filter_name]['nx']
     ny=options[filter_name]['ny']
     nz=options[filter_name]['nz']

     #Get the angles that will be used based on the selected threshold.
     tmp_angles= output['elevations'][ output['elevations'] < options[filter_name]['min_angle']]
     tmp_n_angles = np.size( tmp_angles )

     output['smooth_ref']=qc.box_functions_2d(datain=output['ref'],na=na,nr=nr,ne=ne,undef=output['undef_ref']
                                               ,boxx=nx,boxy=ny,boxz=nz,operation='MEAN',threshold=0.0)
     
     tmp_index=np.zeros([na,nr,ne]) 
     for ie in range( 0 , tmp_n_angles )  :
        tmp_index[:,:,ie]=np.logical_and( output['ref'][:,:,ie] > options['norainrefval'] , output['smooth_ref'][:,:,tmp_n_angles] <= options['norainrefval'] )
        tmp_index[:,:,ie][ output['altitude'][:,:,tmp_n_angles] > options[filter_name]['height_thr'] ] = 0.0

     tmp_index[ tmp_index == output['undef_ref'] ] = options['undef']

     output = output_update( output , tmp_index , options , filter_name ) 

   return radar , output 

#===================================================
# RHO HV FILTER
#===================================================

def RhoFilter( radar , output , options )  :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.

   filter_name = 'RhoFilter'

   if options['name_rho'] in radar.fields  :

      na=output['na']
      nr=output['nr']
      ne=output['ne']
      nx=options[filter_name]['nx']
      ny=options[filter_name]['ny']
      nz=options[filter_name]['nz']

      output['undef_rho']=radar.fields[ options['name_rho'] ]['_FillValue']

      [ rhohv , dm , dm , dm , dm  ]=order_variable( radar , options['name_rho'] , output['undef_rho']  )

      #Compute the filter parameter
      tmp_index=qc.box_functions_2d(datain=rhohv,na=na,nr=nr,ne=ne,undef=output['undef_rho']
                                              ,boxx=nx,boxy=ny,boxz=nz,operation='MEAN',threshold=0.0)
      #Compute smoothed reflectivity (areas with large reflectivity wont be affected by the filter)
      smooth_ref=qc.box_functions_2d(datain=output['ref'],na=na,nr=nr,ne=ne,undef=output['undef_ref']
                                              ,boxx=nx,boxy=ny,boxz=nz,operation='MEAN',threshold=0.0)

      tmp_index[ tmp_index == output['undef_rho'] ] = options['undef']

      tmp_mask = np.logical_or( smooth_ref > options[filter_name]['ref_threshold'] , output['cref'] == options['norainrefval'] )

      tmp_index[tmp_mask] = options['undef']

      #import matplotlib.pyplot as plt
      #tmp_index[tmp_index==options['undef']]=np.nan
      #plt.pcolor( tmp_index[:,:,0] )
      #plt.colorbar()
      #plt.show()

      output = output_update( output , tmp_index , options , filter_name ) 

   else   :
 
      print('Could not found variable ' + options['name_rho'] )
      print(radar.fields.keys() )
      

   return radar , output
 
#===================================================
# REFLECTIVITY SPECKLE FILTER
#===================================================

def RefSpeckleFilter( radar , output , options )    :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.

   filter_name='RefSpeckleFilter'

   if options['name_ref'] in radar.fields :

      #Compute the number pixels with reflectivities over spfiltertr sourrounding each pixels in the box defined by nx,ny,nz.
      na=output['na']
      nr=output['nr']
      ne=output['ne']

      nx=options[filter_name]['nx']
      ny=options[filter_name]['ny']
      nz=options[filter_name]['nz']
      tr=options[filter_name]['reftr']


      tmp_index=qc.box_functions_2d(datain=output['ref'].data,na=na,nr=nr,ne=ne,undef=output['undef_ref']
                                                ,boxx=nx,boxy=ny,boxz=nz,operation='COU2',threshold=tr) 

      tmp_index[ tmp_index == output['undef_ref'] ] = options['undef']

      tmp_index[ output['ref'] == options['norainrefval'] ] = options['undef']

      output = output_update( output , tmp_index , options , filter_name )

   return radar , output 

#===================================================
# DOPPLER SPECKLE FILTER
#===================================================

def DopplerSpeckleFilter( radar , output , options )   :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.

   filter_name='DopplerSpeckleFilter'

   if  options['name_v'] in radar.fields :

      na=output['na']
      nr=output['nr']
      ne=output['ne']
      nx=options[filter_name]['nx']
      ny=options[filter_name]['ny']
      nz=options[filter_name]['nz']
      tr=options[filter_name]['dvtr']

      tmp=np.abs(output['v'].data)

      tmp[ output['v'] == output['undef_v'] ] = output['undef_v']

      tmp_index=qc.box_functions_2d(datain=tmp,na=na,nr=nr,ne=ne,undef=output['undef_v']
                                                ,boxx=nx,boxy=ny,boxz=nz,operation='COU2',threshold=tr)

      tmp_index[ tmp_index == output['undef_v'] ] = options['undef']

      output = output_update( output , tmp_index , options , filter_name )

   return radar , output 

#===================================================
# ATTENUATION FILTER
#===================================================

def AttenuationFilter( radar , output , options )   :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.

   filter_name='AttenuationFilter'

   if options['name_ref'] in radar.fields :

      na=output['na']
      nr=output['nr']
      ne=output['ne']

      beaml=radar.range['data'][1]-radar.range['data'][0] #Get beam length

      [tmp_index,corrected_ref]=qc.get_attenuation( var=output['cref'],na=na,nr=nr,ne=ne,undef=output['undef_ref']
                                                ,beaml=beaml,cal_error=options[filter_name]['attcalerror']
                                                ,is_power=options[filter_name]['is_power']
                                                ,coefs=options[filter_name]['att_coefs'],mindbz=options['norainrefval'] )

      tmp_index[ tmp_index == output['undef_ref'] ] = options['undef']

      if options[filter_name]['attenuation_correction'] & ( options['name_ref'] in radar.fields )  :

          output['cref'] = corrected_ref 

      output = output_update( output , tmp_index , options , filter_name )

   return radar , output 

#===================================================
# TOPOGRAPHY BLOCKING FILTER
#===================================================

def BlockingFilter( radar , output , options )   :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.

   na=output['na']
   nr=output['nr']
   ne=output['ne']

   filter_name='BlockingFilter'  #This is not included in the fuzzy logic algorithm

   tmp_index=qc.compute_blocking( radarz=output['altitude'] , topo=output['topo'] , na=na , nr=nr , ne=ne      ,
                                              undef=options['undef'],
                                              radar_beam_width_v=radar.instrument_parameters['radar_beam_width_v']['data'] , 
                                              beam_length=radar.range['meters_between_gates']                              , 
                                              radarrange=radar.range['data'] , radarelev=output['elevations'] )  

   #Compute correction 
   if options[filter_name]['blocking_correction'] & ( options['name_ref'] in radar.fields )  :
      #Correct partially blocked precipitation echoes.
      mask=np.logical_and( tmp_index > 0.1 , tmp_index <= 0.3 )
      mask=np.logical_and( mask , output['cref'] > options['norainrefval'] )

      output['cref'][mask] = output['cref'][mask] + 1.0

      mask=np.logical_and( tmp_index > 0.3 , tmp_index <= 0.4 )
      mask=np.logical_and( mask , output['cref'] > options['norainrefval'] )

      output['cref'][mask] = output['cref'][mask] + 2.0

      mask= tmp_index > 0.4 
      mask=np.logical_and( mask , output['cref'] > options['norainrefval'] )

      output['cref'][mask] = output['cref'][mask] + 3.0

   output = output_update( output , tmp_index , options , filter_name ) 

   return radar , output 

#===================================================
# DOPPLER NOISE FILTER
#===================================================

def DopplerNoiseFilter( radar , output , options )   :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.

   na=output['na']
   nr=output['nr']
   ne=output['ne']

   filter_name='DopplerNoiseFilter'

   if options['name_v'] in radar.fields  :


     v=np.copy(output['v'])
     tmp_v=np.copy(output['v'])

     for ii in range(len(options[filter_name]['n_filter_pass']))  :

      nx=options[filter_name]['nx'][ii]
      ny=options[filter_name]['ny'][ii]
      nz=options[filter_name]['nz'][ii]
      tr=options[filter_name]['threshold'][ii]
      nf=options[filter_name]['n_filter_pass'][ii]

      #First round with small local domains.
      for ip in range(0,nf) :

         #Smooth v is the smoothed velocity without the outliers.
         smooth_v=qc.box_functions_2d(datain=tmp_v,na=na,nr=nr,ne=ne,undef=output['undef_v']
                                               ,boxx=nx,boxy=ny,boxz=nz,operation='MEAN',threshold=0.0)

         #Compute the distance between the original wind field and the smoothed wind field without the outliers.
         distance = np.abs( smooth_v - v )

         #Remove those pixels that are far from the local mean.
         tmp_v=np.copy( v )
         tmp_v[ distance > tr ] = output['undef_v']
         tmp_v[ distance == output['undef_v'] ] = output['undef_v']

      #At the end of the loop (multiple filter passes)
      #Remove all the pixels that are far from the local mean.
      v[ distance > tr ] = output['undef_v']
      tmp_v = np.copy( v ) 

     #Get the pixels in which v is undef but the original v is not (those are the ones removed by the filter)
     tmp_index = np.logical_and( output['v'] != output['undef_v'] , v == output['undef_v'] ).astype(int)   

     tmp_index[ tmp_index == output['undef_v'] ] = options['undef']


     output = output_update( output , tmp_index , options , filter_name )

   return radar , output 

#===================================================
#  DOPPLER LOCAL STD FILTER
#===================================================

def DopplerLocalStdFilter( radar , output , options )   :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.

   filter_name = 'DopplerLocalStdFilter'

   if options['name_v'] in radar.fields :

      na=output['na']
      nr=output['nr']
      ne=output['ne']
      nx=options[filter_name]['nx']
      ny=options[filter_name]['ny']
      nz=options[filter_name]['nz']

      #Compute the filter parameter
      tmp_index=qc.box_functions_2d(datain=output['v'],na=na,nr=nr,ne=ne,undef=output['undef_v']
                                               ,boxx=nx,boxy=ny,boxz=nz,operation='SIGM',threshold=0.0)

      tmp_index[ tmp_index == output['undef_v'] ] = options['undef']

      output = output_update( output , tmp_index , options , filter_name ) 

   return radar , output 

#===================================================
# INTERFERENCE FILTER
#===================================================

def InterferenceFilter( radar , output , options )   :

   import numpy as np
   #Clean rays affected by interference.
   filter_name = 'InterferenceFilter'

   if options['name_ref'] in radar.fields :

      tmp_ref=np.copy( output['ref'] )

      #This function returns 1 if the pixel is sopossed to be contaminated by interference and
      # 0 otherwise.
      tmp_index = interference_filter ( tmp_ref , output['undef_ref'] , options['norainrefval'] 
                                            , radar.range['data'] , options[filter_name] ) 

      tmp_index[ tmp_index == output['undef_ref'] ] = options['undef']

      output = output_update( output , tmp_index , options , filter_name ) 

   return radar , output 

#===================================================
# DOPPLER SPATIAL COHERENCE FILTER
#===================================================

def DopplerSpatialCoherenceFilter( radar , output , options )   :

   import numpy as np
   #Clean pixels which shows no coherence with their neighbors

   filter_name = 'DopplerSpatialCoherenceFilter'

   if options['name_v'] in radar.fields  :

      tmp_v=np.copy( output['v'] )

      tmp_index = dopplerspatialcoherence_filter( tmp_v , output['undef_v'] , options[filter_name] )

      tmp_index[ tmp_index == output['undef_v'] ] = options['undef']

      output = output_update( output , tmp_index , options , filter_name )

   return radar , output 

#===================================================
# POWER FILTER
#===================================================

def PowerFilter( radar , output , options )   :

   import numpy as np

   na=output['na']
   nr=output['nr']
   ne=output['ne']

   filter_name = 'PowerFilter' 

   if  options['name_ref'] in radar.fields  :

      #TODO TODO TODO
      #a = 0.01
      #C = 0
      #EXCLUDE_BELOW = -10  # 3.16

      #rango = radar.range['data'] / 1000
      #ref = np.copy(radar.fields[var]['data'])
      #pot = ref - 20*np.log10(rango) - 2*a*rango - C
      #pot = np.ma.masked_where(pot<-200, pot)
      #ref_cor = np.copy(ref)
      #umbral = np.floor(np.nanmin(pot[np.isfinite(pot)])) + umbral_desde_piso - 5
      #print(umbral)
      #ref_cor = np.ma.masked_where(pot<umbral, ref_cor)
      #radar.add_field_like(var, 'TH_cor', ref_cor, True)

      tmp_index[ tmp_index == output['undef_ref'] ] = options['undef']

      output = output_update( output , tmp_index , options , filter_name )

   #TODO

   return radar , output 

#===================================================
# DETECT MISSING REFLECTIVITY VALUES
#===================================================

def MissingRefFilter( radar , output , options )    :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.

   filter_name='MissingRefFilter'

   if options['name_ref'] in radar.fields  :

      na=output['na']
      nr=output['nr']
      ne=output['ne']

      tmp_index = qc.detect_missing(  output['ref'],na=na,nr=nr,ne=ne,undef=output['undef_ref']
                                                   ,min_ref=options['norainrefval'],threshold=options[filter_name]['threshold'] 
                                                   ,nmissing_max=options[filter_name]['nmissing_max'] )

      tmp_index = tmp_index.astype(int)

      tmp_index[ tmp_index == output['undef_ref'] ] = options['undef']

      output = output_update( output , tmp_index , options , filter_name )

   return radar , output 

#===================================================
# REFLECTIIVTY TEXTURE FILTER
#===================================================

def ReflectivityTextureFilter( radar , output , options )   :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.

   filter_name='ReflectivityTextureFilter'

   if  options['name_ref'] in radar.fields  : 

     na=output['na']
     nr=output['nr']
     ne=output['ne']

     nx=options[filter_name]['nx']
     ny=options[filter_name]['ny']
     nz=options[filter_name]['nz']
     tmp_index = qc.compute_texture(var=output['ref'],na=na,nr=nr,ne=ne,undef=output['undef_ref'],nx=nx,ny=ny,nz=nz)


     if options[filter_name]['use_smooth_ref'] :
        smooth_ref =qc.box_functions_2d(datain=output['ref'],na=na,nr=nr,ne=ne,undef=output['undef_ref'],
                                                 boxx=0,boxy=0,boxz=0,operation='MEAN',threshold=0.0)

        #High reflectivity cores will not be affected by texture filter.
        tmp_index[ smooth_ref >= options[filter_name]['smooth_ref_tr']  ]= 0.0

     tmp_index[output['ref']==output['undef_ref']] = 0.0

     tmp_index[ tmp_index == output['undef_ref'] ] = options['undef']

     output = output_update( output , tmp_index , options , filter_name )

   return radar , output 

#===================================================
# DOPPER VELOCITY TEXTURE FILTER
#===================================================

def DopplerTextureFilter( radar , output , options )   :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.

   filter_name='DopplerTextureFilter'

   if  options['name_v'] in radar.fields  :

     na=output['na']
     nr=output['nr']
     ne=output['ne']

     nx=options[filter_name]['nx']
     ny=options[filter_name]['ny']
     nz=options[filter_name]['nz']


     tmp_index = qc.compute_texture(var=output['v'],na=na,nr=nr,ne=ne,undef=output['undef_v'],nx=nx,ny=ny,nz=nz)

     tmp_index[ tmp_index == output['undef_v'] ] = options['undef']

     output = output_update( output , tmp_index , options , filter_name )

   return radar , output

#===================================================
# LOW DOPPLER VOLOCITY FILTER
#===================================================

def LowDopplerFilter( radar , output , options )   :
 
   #This filter removes reflectivity values which are close to the surface and 
   #have a low associated doppler velocity. 
   import numpy as np

   filter_name='LowDopplerFilter'

   if  options['name_v'] in radar.fields :

      na=output['na']
      nr=output['nr']
      ne=output['ne']

      tmp_index=np.abs(output['v'])
      tmp_index[ (output['altitude'] - output['topo']) > options[filter_name]['height_thr'] ]=10.0

      tmp_index[ tmp_index == output['undef_v'] ] = options['undef']
 
      output = output_update( output , tmp_index , options , filter_name )
 
   return radar , output 

#===================================================
# REF RANGE FILTER
#===================================================

def RefRangeFilter( radar , output , options )  :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.

   filter_name = 'RefRangeFilter'

   if options['name_ref'] in radar.fields  :

      tmp_index = np.zeros( np.shape( output['cref'] ) )

      tmp_index[ np.logical_or( output['cref'] > options[filter_name]['max'] , output['cref'] < options[filter_name]['min'] ) ] = 1.0

      output = output_update( output , tmp_index , options , filter_name )

   return radar , output

#===================================================
# DOPPLER RANGE FILTER
#===================================================

def DopplerRangeFilter( radar , output , options )  :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.

   filter_name = 'DopplerRangeFilter'

   if options['name_v'] in radar.fields  :

      tmp_index = np.zeros( np.shape( output['cv'] ) )

      tmp_index[ np.logical_or( output['cv'] > options[filter_name]['max'] , output['cv'] < options[filter_name]['min'] ) ] = 1.0
 
      output = output_update( output , tmp_index , options , filter_name )

   return radar , output

#===========================================================================================================
# OTRAS FUNCIONES CONTENIDAS EN ESTE MODULO
#===========================================================================================================   

#From Pyart order to ARE (Azmimuth , range , elevation )

def order_variable ( radar , var_name , undef )  :  

   import numpy as np
   import numpy.ma as ma
   #import warnings 
   #import matplotlib.pyplot as plt

   #From azimuth , range -> azimuth , range , elevation 

   if radar.ray_angle_res != None   :
      #print( radar.ray_angle_res , radar.ray_angle_res == None )
      ray_angle_res = np.unique( radar.ray_angle_res['data'] )
   else                             :
      print('Warning: ray_angle_res no esta definido, estimo la resolucion en radio como la diferencia entre los primeros angulos')
      ray_angle_res = np.min( np.abs( radar.azimuth['data'][1:] - radar.azimuth['data'][0:-1] ) )
      print('La resolucion en rango estimada es: ',ray_angle_res)


   if( np.size( ray_angle_res ) >= 2 )  :
      print('Warning: La resolucion en azimuth no es uniforme en los diferentes angulos de elevacion ')
      print('Warning: El codigo no esta preparado para considerar este caso y puede producir efectos indeseados ')
   ray_angle_res=np.nanmean( ray_angle_res )

   levels=np.sort( np.unique(radar.elevation['data']) )
   nb=radar.azimuth['data'].shape[0]

   order_azimuth=np.arange(0.0,360.0,ray_angle_res) #Asuming a regular azimuth grid

   na=np.size(order_azimuth)
   ne=np.size(levels)
   nr=np.size(radar.range['data'].data) 


   var = np.ones( (nb,nr) )

   if ( var_name == 'altitude' ) :
       var[:]=radar.gate_altitude['data']  
   elif( var_name == 'longitude' ) :
       var[:]=radar.gate_longitude['data']
   elif( var_name == 'latitude'  ) :
       var[:]=radar.gate_latitude['data']
   elif( var_name == 'x' )         :
       var[:]=radar.gate_x['data']
   elif( var_name == 'y' )         : 
       var[:]=radar.gate_y['data']
   else  :
       var[:]=radar.fields[var_name]['data'].data


   #Allocate arrays
   order_var    =np.zeros((na,nr,ne))
   order_time   =np.zeros((na,ne)) 
   azimuth_exact=np.zeros((na,ne))
   order_n      =np.zeros((na,nr,ne),dtype='int')
   
   current_lev = radar.elevation['data'][0]
   ilev = np.where( levels == current_lev )[0]

   for iray in range( 0 , nb )  :   #Loop over all the rays
 
     #Check if we are in the same elevation.
     if  radar.elevation['data'][iray] != current_lev  :
         current_lev = radar.elevation['data'][iray]
         ilev=np.where( levels == current_lev  )[0]

     #Compute the corresponding azimuth index.
     az_index = np.round( radar.azimuth['data'][iray] / ray_angle_res ).astype(int)
     #Consider the case when azimuth is larger than na*ray_angle_res-(ray_angle_res/2)
     if az_index >= na   :  
        az_index = 0

     tmp_var = var[iray,:]
     undef_mask = tmp_var == undef 
     tmp_var[ undef_mask ] = 0.0
    
     order_var [ az_index , : , ilev ] = order_var [ az_index , : , ilev ] + tmp_var
     order_n   [ az_index , : , ilev ] = order_n   [ az_index , : , ilev ] + np.logical_not(undef_mask).astype(int)

     order_time[ az_index , ilev ] = order_time[ az_index , ilev ] + radar.time['data'][iray]
     azimuth_exact[ az_index , ilev ] = azimuth_exact[ az_index , ilev ] + radar.azimuth['data'][ iray ]

   order_var[ order_n > 0 ] = order_var[ order_n > 0 ] / order_n[ order_n > 0 ]
   order_var[ order_n == 0] = undef

   return order_var , order_azimuth , levels , order_time , azimuth_exact

def order_variable_inv (  radar , var , undef )  :

   import numpy as np
   
   #From azimuth , range , elevation -> azimuth , range

   na=var.shape[0]
   nr=var.shape[1]
   ne=var.shape[2]

   nb=radar.azimuth['data'].shape[0]

   levels=np.sort( np.unique(radar.elevation['data']) )

   if radar.ray_angle_res != None   :
      #print( radar.ray_angle_res , radar.ray_angle_res == None )
      ray_angle_res = np.unique( radar.ray_angle_res['data'] )
   else                             :
      print('Warning: ray_angle_res no esta definido, estimo la resolucion en radio como la diferencia entre los primeros angulos')
      ray_angle_res = np.min( np.abs( radar.azimuth['data'][1:] - radar.azimuth['data'][0:-1] ) )
      print('La resolucion en rango estimada es: ',ray_angle_res)

   if( np.size( ray_angle_res ) >= 2 )  :
      print('Warning: La resolucion en azimuth no es uniforme en los diferentes angulos de elevacion ')
      print('Warning: El codigo no esta preparado para considerar este caso y puede producir efectos indesaedos ')
   ray_angle_res=np.nanmean( ray_angle_res )

   current_lev = radar.elevation['data'][0]
   ilev = np.where( levels == current_lev  )[0]

   output_var = np.zeros((nb,nr) )
   output_var[:] = undef

   for iray in range( 0 , nb )  :   #Loop over all the rays

      #Check if we are in the same elevation.
      if  radar.elevation['data'][iray] != current_lev  :
          current_lev = radar.elevation['data'][iray]
          ilev=np.where( levels == current_lev  )[0]

      #Compute the corresponding azimuth index.
      az_index = np.round( radar.azimuth['data'][iray] / ray_angle_res ).astype(int)
      #Consider the case when azimuth is larger than na*ray_angle_res-(ray_angle_res/2)
      if az_index >= na   :
         az_index = 0

      output_var[ iray , : ] = var[ az_index , : , ilev ]

   return output_var

def dopplerspatialcoherence_filter( v , undef , my_conf ) : 

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.
   from sklearn import linear_model, datasets
   import matplotlib.pyplot as plt

   #Compute a robust correlation between a ray and its neighbors.
   #Detect outliers based on the correlation. 

   na=v.shape[0]
   nr=v.shape[1]
   ne=v.shape[2]

   #Coherence mask starts flagging the undef values
   available_mask = ( v != undef )

   coherence_index = np.zeros( np.shape( v ) )

   ransac = linear_model.RANSACRegressor()

   if my_conf['compute_horizontal_coherence']   :
   #CONSIDER THE COHERENCE WITH THE NEIGHBOR BEAMS IN AZIMUTH DIRECTION

      for k in range(ne)  :

         for i in range(na)  :

              iprev = i - 1
              if i == 0    :
                   iprev = na-1


              #First get the number of undef values and generate a common undef mask.
              local_valid_mask = np.logical_and( available_mask[i,:,k] , available_mask[iprev,:,k] )

              n_valid = np.sum( local_valid_mask.astype(int) ) 

              coherence_index[i,:,k][ np.logical_not( local_valid_mask ) ]=coherence_index[i,:,k][ np.logical_not( local_valid_mask ) ] + 1.0

              if ( n_valid / nr ) > my_conf['threshold_undef']   :

                 tmp_vi = np.copy( v[i,:,k][local_valid_mask] )
                 tmp_viprev = np.copy( v[iprev,:,k][local_valid_mask] )

                 if my_conf['consistency_metric'] == 'ransac'   :
                    #Use a robust correlation method to decide which
                    #elements are spatially coherent.

                    ransac.fit( tmp_vi.reshape(-1, 1) , tmp_viprev.reshape(-1, 1) )
                    inlier_mask = ransac.inlier_mask_
                    outlier_mask = np.logical_not( inlier_mask )

                 if my_conf['consistency_metric'] == 'constant'  :
                     #Use a constant threshold to define spatial coherence.

                     inlier_mask = np.abs( tmp_vi - tmp_viprev ) < my_conf['constant_consistency_threshold'] 
                     outlier_mask = np.logical_not( inlier_mask )


                 corrcoef=np.corrcoef( tmp_vi[inlier_mask] ,  tmp_viprev[inlier_mask] )[0,1]

                 tmp=coherence_index[i,:,k][local_valid_mask]
                 tmp[outlier_mask]=tmp[outlier_mask] + 3.0
                 coherence_index[i,:,k][local_valid_mask]=tmp

                 tmp=coherence_index[iprev,:,k][local_valid_mask]
                 tmp[outlier_mask]=tmp[outlier_mask] + 3.0
                 coherence_index[iprev,:,k][local_valid_mask]=tmp

                 if corrcoef < my_conf['threshold_corr']   :

                    coherence_index[i,:,k] = coherence_index[i,:,k] + 3.0
                    coherence_index[iprev,:,k] = coherence_index[iprev,:,k] + 3.0
              else  :

                 coherence_index[i,:,k] = coherence_index[i,:,k] + 3.0

   #import matplotlib.pyplot
   #plt.pcolor( coherence_index[:,:,2] )
   #plt.colorbar()

   if my_conf['compute_vertical_coherence']   :
   #CONSIDER THE COHERENCE WITH THE NEIGHBOR BEAMS IN VERTICAL DIRECTION

      for k in range(1,ne)  :
 
         kprev = k - 1

         for i in range(na)  :


              #First get the number of undef values and generate a common undef mask.
              local_valid_mask = np.logical_and( available_mask[i,:,k] , available_mask[i,:,kprev] )

              n_valid = np.sum( local_valid_mask.astype(int) )

              coherence_index[i,:,k][ np.logical_not( local_valid_mask ) ]=coherence_index[i,:,k][ np.logical_not( local_valid_mask ) ] + 1.0

              if ( n_valid / nr ) > my_conf['threshold_undef']   :

                 tmp_vi = np.copy( v[i,:,k][local_valid_mask] )
                 tmp_viprev = np.copy( v[i,:,kprev][local_valid_mask] )

                 if my_conf['consistency_metric'] == 'ransac'   :
                    #Use a robust correlation method to decide which
                    #elements are spatially coherent.

                    ransac.fit( tmp_vi.reshape(-1, 1) , tmp_viprev.reshape(-1, 1) )
                    inlier_mask = ransac.inlier_mask_
                    outlier_mask = np.logical_not( inlier_mask )

                 if my_conf['consistency_metric'] == 'constant'  :
                     #Use a constant threshold to define spatial coherence.

                     inlier_mask = np.abs( tmp_vi - tmp_viprev ) < my_conf['constant_consistency_threshold']
                     outlier_mask = np.logical_not( inlier_mask )

                 #Compute the correlation between the adjacent rays for the inlier mask. 
                 corrcoef=np.corrcoef( tmp_vi[inlier_mask] ,  tmp_viprev[inlier_mask] )[0,1]

                 tmp=coherence_index[i,:,k][local_valid_mask]
                 tmp[outlier_mask]=tmp[outlier_mask] + 3.0
                 coherence_index[i,:,k][local_valid_mask]=tmp

                 tmp=coherence_index[i,:,kprev][local_valid_mask]
                 tmp[outlier_mask]=tmp[outlier_mask] + 3.0
                 coherence_index[i,:,kprev][local_valid_mask]=tmp

                 #coherence_index[i,:,k][local_valid_mask][outlier_mask] = coherence_index[i,:,k][local_valid_mask][outlier_mask] + 3.0 
                 #coherence_index[iprev,:,k][local_valid_mask][outlier_mask] = coherence_index[iprev,:,k][local_valid_mask][outlier_mask] + 3.0 

                 if corrcoef < my_conf['threshold_corr']   :

                    coherence_index[i,:,k] = coherence_index[i,:,k] + 3.0
                    coherence_index[i,:,kprev] = coherence_index[iprev,:,kprev] + 3.0

   for ifilter in range( my_conf['npass_filter'] )  :

      v[ coherence_index > my_conf['threshold_coherence_index'] ] = undef

      for k in range(ne) :

         for i in range(na) :
            if  my_conf['azimuthfilter']         :  #DETECT ISOLATED PIXELS IN AZIMUTH

               if ( i > 1 ) & ( i < na-2 ) :

                  #If we have in only one ray but not in the neighbors this suggest an interference pattern.
                  tmp_mask = np.logical_and( v[i-1,:,k] == undef , v[i+1,:,k] == undef )
                  coherence_index[i,:,k][ np.logical_and( tmp_mask , v[i,:,k] != undef  ) ] = 10.0

                  tmp_mask = np.logical_and( v[i-2,:,k] == undef , v[i+2,:,k] == undef )
                  coherence_index[i,:,k][ np.logical_and( tmp_mask , v[i,:,k] != undef ) ] = 10.0

               elif  i==na-1   :
                  tmp_mask = np.logical_and( v[i-1,:,k] == undef , v[0,:,k] == undef )
                  coherence_index[i,:,k][ np.logical_and( tmp_mask , v[i,:,k] != undef ) ] = 10.0

                  tmp_mask = np.logical_and( v[i-2,:,k] == undef , v[1,:,k] == undef )
                  coherence_index[i,:,k][ np.logical_and( tmp_mask , v[i,:,k] != undef ) ] = 10.0

               elif  i==na-2   :
                  tmp_mask = np.logical_and( v[i-1,:,k] == undef , v[i,:,k] == undef )
                  coherence_index[i,:,k][ np.logical_and( tmp_mask , v[i,:,k] != undef ) ] = 10.0

                  tmp_mask = np.logical_and( v[i-2,:,k] == undef , v[0,:,k] == undef )
                  coherence_index[i,:,k][ np.logical_and( tmp_mask , v[i,:,k] != undef ) ] = 10.0

               elif  i==0      :

                  tmp_mask = np.logical_and( v[na-1,:,k] == undef , v[i+1,:,k] == undef )
                  coherence_index[i,:,k][ np.logical_and( tmp_mask , v[i,:,k] != undef ) ] = 10.0

                  tmp_mask = np.logical_and( v[na-2,:,k] == undef , v[i+2,:,k] == undef )
                  coherence_index[i,:,k][ np.logical_and( tmp_mask , v[i,:,k] != undef ) ] = 10.0

               elif  i==1      :

                  tmp_mask = np.logical_and( v[i-1,:,k] == undef  , v[i+1,:,k] == undef )
                  coherence_index[i,:,k][ np.logical_and( tmp_mask , v[i,:,k] != undef ) ] = 10.0

                  tmp_mask = np.logical_and( v[na-1,:,k] == undef , v[i+2,:,k] == undef )
                  coherence_index[i,:,k][ np.logical_and( tmp_mask , v[i,:,k] != undef ) ] = 10.0

            if  my_conf['elevationfilter']       :   #DETECT ISOLATED PIXELS IN ELEVATION

               if ( k > 0 ) & ( k < ne-1 ) :

                  tmp_mask = np.logical_and( v[i,:,k-1] == undef , v[i,:,k+1] == undef )
                  coherence_index[i,:,k][ np.logical_and( tmp_mask , v[i,:,k] != undef ) ] = 10.0

               if ( k == 0 )                :

                  tmp_mask = np.logical_and( v[i,:,k+2] == undef , v[i,:,k+1] == undef )
                  coherence_index[i,:,k][ np.logical_and( tmp_mask , v[i,:,k] != undef ) ] = 10.0

               if ( k == ne-1 )            :

                  tmp_mask = np.logical_and( v[i,:,k-2] == undef  , v[i,:,k-1] == undef )
                  coherence_index[i,:,k][ np.logical_and( tmp_mask , v[i,:,k] != undef ) ] = 10.0

         for i in range(nr) :

            if  my_conf['rangefilter']           :   #DETECT ISOLATED PIXELS IN RANGE

               if ( i  > 0 ) & ( i < nr-1 ) :

                  tmp_mask = np.logical_and( v[:,i-1,k] == undef , v[:,i+1,k] == undef )
                  coherence_index[:,i,k][ np.logical_and( tmp_mask , v[:,i,k] != undef ) ] = 10.0

               if ( i == 0 )                :

                  tmp_mask = np.logical_and( v[:,i+2,k] == undef , v[:,i+1,k] == undef )
                  coherence_index[:,i,k][ np.logical_and( tmp_mask , v[:,i,k] != undef ) ] = 10.0

               if ( i == nr-1 )            :

                  tmp_mask = np.logical_and( v[:,i-2,k] == undef  , v[:,i-1,k] == undef )
                  coherence_index[:,i,k][ np.logical_and( tmp_mask , v[:,i,k] != undef ) ] = 10.0

   if my_conf['enable_speckle']   :

      tr=0.0
      nx=my_conf['nx']
      ny=my_conf['ny']
      nz=my_conf['nz']
      tmp=np.abs( v )
      tmp[ v == undef ] = undef

      speckle_v=qc.box_functions_2d(datain=tmp,na=na,nr=nr,ne=ne,undef=undef
                                 ,boxx=nx,boxy=ny,boxz=nz,operation='COU2',threshold=tr)

      coherence_index[ speckle_v < my_conf['speckle_threshold'] ] = 10.0

   #v[ coherence_index > my_conf['threshold_coherence_index'] ] = undef

   return coherence_index

def interference_filter ( ref , undef , min_ref , r , my_conf )  :

   import numpy as np
   from common_qc_tools  import qc  #Fortran code routines.
   from sklearn import linear_model, datasets
   import matplotlib.pyplot as plt

   #ref a 3D array with the reflectivity in the na,nr,ne grid.
   #undef: value corresponding to the missing value in the reflectivity field.
   #min_ref: value corresponding to the no-rain value in the reflectivity field.
   #r : radar range (vector)
   #my_conf : configuration for this filter.

   ref[ ref == min_ref ] = undef

   offset= my_conf['offset']

   na=ref.shape[0]
   nr=ref.shape[1]
   ne=ref.shape[2]

   nx=my_conf['nx']
   ny=my_conf['ny']
   nz=my_conf['nz']

   if my_conf['Smooth_Ref']  :

      tmp_ref = qc.box_functions_2d(datain=ref,na=na,nr=nr,ne=ne,undef=undef
                                 ,boxx=nx,boxy=ny,boxz=nz,operation='MEAN',threshold=0.0)
   else                      :
      tmp_ref = np.copy( ref )

   tmp_index=np.zeros( np.shape( ref ) ) 

   #tmp_z = np.power( 10.0, tmp_ref / 10.0  )

   #tmp_z[tmp_ref == undef ] = undef

   
   #offset=my_conf['offset']
   att=my_conf['att']
   AzimuthFilter=my_conf['AzimuthFilter']
   ElevationFilter=my_conf['ElevationFilter']
   npass_filter=my_conf['npass_filter']
   percent_valid_threshold=my_conf['percent_valid_threshold']
   corr_threshold=my_conf['corr_threshold']
   ref_threshold=my_conf['ref_threshold']
   percent_ref_threshold=my_conf['percent_ref_threshold']
   azimuth_ref_diff_threshold=my_conf['azimuth_ref_diff_threshold']
   
   Power_Regression = my_conf['Power_Regression']

   #Main loops

   for k in range(ne)  :

      for i in range(na)  :

           local_sref = np.copy( tmp_ref[i,:,k] )
           local_sref[0:offset]=undef

           undef_mask = local_sref != undef 

           #Get the following and previous rays.

           local_sref_m1 = np.copy( tmp_ref[i-1,:,k] )
 
           if i < na-1   :
              local_sref_p1 = np.copy( tmp_ref[i+1,:,k] )
           else          :
              local_sref_p1 = np.copy( tmp_ref[na-1,:,k] )
              

           local_ref = np.copy( ref[i,:,k] )
           local_ref[0:offset] =undef

           tmp_count = np.sum( (undef_mask).astype(int) )/nr


           #Local smooth power
           local_spower = np.power( 10.0 , ( local_sref  - 20.0 * np.log10( r ) - 2.0 * att * r ) / 10.0 )
           local_spower[ local_sref == undef  ] = undef
           #Local power
           local_power  = np.power( 10.0 , ( ref[i,:,k] - 20.0 * np.log10( r ) - 2.0 * att * r ) / 10.0 )
           local_power[ ref[i,:,k] == undef ] = undef 

           local_inlier_mask = undef_mask 

           if tmp_count > percent_valid_threshold   :
              ransac = linear_model.RANSACRegressor()     

              ransac.fit( r[undef_mask].reshape(-1, 1) , local_spower[undef_mask].reshape(-1, 1) )

              local_inlier_mask [ undef_mask ] = ransac.inlier_mask_
              #inlier_mask = ransac.inlier_mask_

              local_fit_power = ransac.predict( r.reshape(-1,1) )[:,0]
              local_fit_power[ local_fit_power < 0.0 ]=10e-10

              local_fit_ref = 10.0*np.log10( local_fit_power ) + 20.0 * np.log10( r ) + 2.0 * att * r


           else:

              local_fit_power = np.zeros( nr ) 
              local_inlier_mask = np.zeros( nr ).astype(bool) 

           if ( np.sum( local_inlier_mask.astype(int) ) >= 10 )  :
               if ( np.std(local_sref[local_inlier_mask]) > 0 ) :
                  corrcoef=np.corrcoef( local_fit_ref[ local_inlier_mask ],local_sref[ local_inlier_mask ] )[0,1]
                  corrcoef_1=np.corrcoef( r[local_inlier_mask] , local_sref[ local_inlier_mask ] )[0,1]

                  tmp_mask = np.logical_and( local_inlier_mask , local_sref_p1 != undef ) 
                  if np.sum(tmp_mask) >= 10   :
                     azimuth_ref_diff = np.power( np.mean( local_sref[ tmp_mask ] - local_sref_p1[ tmp_mask ] ) , 2) 
                  else                        :
                     azimuth_ref_diff = 0.0
                  tmp_mask = np.logical_and( local_inlier_mask , local_sref_m1 != undef )
                  if np.sum(tmp_mask) >= 10   :
                     azimuth_ref_diff = azimuth_ref_diff + np.power( np.mean( local_sref[ tmp_mask ] - local_sref_m1[ tmp_mask ] ) , 2) 
                  else                        :
                     azimuth_ref_diff = 0.0

                  azimuth_ref_diff = np.sqrt( azimuth_ref_diff / 2.0 )  / np.mean( local_sref[ local_inlier_mask ] ) 

               else                                                 :
                  corrcoef = np.array(0.0)
                  corrcoef_1 = 0.0
                  azimuth_ref_diff = 0.0

           else                                             :

              corrcoef=np.array(0.0)
              corrcoef_1 = 0.0
              azimuth_ref_diff = 0.0


           #if k == 0 :
           #   print( i , corrcoef , np.sum( local_inlier_mask.astype(int) )  )
           filter_ray = False
           
           if np.sum( local_inlier_mask )/(nr-offset) > percent_ref_threshold  :

              for it in range(np.size(corr_threshold) ) :

                  if ( corrcoef > corr_threshold[it] ) & ( azimuth_ref_diff > azimuth_ref_diff_threshold[it] )  :
                     filter_ray = True



           if filter_ray  :

              #This means that this ray is likely to be contaminated by interference.

              undef_mask = ( local_ref != undef )

              #If the reflectivity is far from the fitted interference, and is greather than the fitted
              #Interference, then correct the power substracting the interference power.         
              tmp_mask = np.logical_and( local_sref - local_fit_ref  > ref_threshold ,  undef_mask  )  
              tmp_mask = np.logical_and( local_ref  - local_fit_ref  > 0 , tmp_mask )

              ref[i, tmp_mask ,k] = 10.0*np.log10( local_power[tmp_mask] - local_fit_power[tmp_mask] ) + 20.0 * np.log10( r[tmp_mask] ) + 2.0 * att * r[tmp_mask]
 
              #If the reflectivity is far from the fitted interference, and is smaller than the fitted interference
              #then set that pixel as an undef pixel.
              tmp_mask = np.logical_and( local_sref - local_fit_ref <= ref_threshold , undef_mask )

              ref[i, tmp_mask ,k] = undef

              tmp_index[i , tmp_mask , k] = 1.0


              #if (k == 0) and (i == 61) :
              #   tmp_mask = local_spower != undef
              #   local_spower[tmp_mask]
              #   plt.figure()
              #   plt.plot( r[tmp_mask] , local_spower[tmp_mask] )
              #   plt.plot( r[tmp_mask] , local_fit_power[tmp_mask])
              #   plt.show()

              #   local_fit_ref = 10.0*np.log10( local_fit_power ) + 20.0 * np.log10( r ) + 2.0 * att * r
              #   plt.figure()
                
              #   plt.plot( r[tmp_mask] , local_sref[tmp_mask]      ,'ok' )
              #   plt.plot( r[local_inlier_mask] , local_sref[local_inlier_mask]      ,'ob' )
              #   plt.plot( r[tmp_mask] , local_fit_ref[tmp_mask] )
              #   plt.plot( r[tmp_mask] , local_fit_ref[tmp_mask] - ref_threshold , '--')
              #   plt.plot( r[tmp_mask] , local_fit_ref[tmp_mask] + ref_threshold , '--')
              #   plt.show()


   #Additional filter for the remaining echoes
   #consider cyclic boundary conditions in azimuth.

   for ifilter in range( npass_filter )  :

      for k in range(ne) :

         for i in range(na) :
            if  AzimuthFilter         :  #DETECT ISOLATED PIXELS IN AZIMUTH

               if ( i > 1 ) & ( i < na-2 ) :
 
                  #If we have reflectivity in only one ray but not in the neighbors this suggest an interference pattern.
                  tmp_mask = np.logical_and( ref[i-1,:,k] == undef , ref[i+1,:,k] == undef )
                  tmp_mask = np.logical_and( ref[i,:,k]   != undef , tmp_mask )
                  ref[i,:,k][ tmp_mask ] = undef
                  tmp_index[i,:,k][ tmp_mask ] = 1.0

                  tmp_mask = np.logical_and( ref[i-2,:,k] == undef , ref[i+2,:,k] == undef )
                  tmp_mask = np.logical_and( ref[i,:,k]   != undef , tmp_mask )
                  ref[i,:,k][ tmp_mask ] = undef
                  tmp_index[i,:,k][ tmp_mask ] = 1.0
 
               elif  i==na-1   :
                  tmp_mask = np.logical_and( ref[i-1,:,k] == undef, ref[0,:,k] == undef )
                  tmp_mask = np.logical_and( ref[i,:,k]   != undef , tmp_mask )
                  ref[i,:,k][ tmp_mask ] = undef
                  tmp_index[i,:,k][ tmp_mask ] = 1.0

                  tmp_mask = np.logical_and( ref[i-2,:,k] == undef , ref[1,:,k] == undef )
                  tmp_mask = np.logical_and( ref[i,:,k]   != undef , tmp_mask )
                  ref[i,:,k][ tmp_mask ] = undef
                  tmp_index[i,:,k][ tmp_mask ] = 1.0

               elif  i==na-2   :
                  tmp_mask = np.logical_and( ref[i-1,:,k] == undef , ref[i,:,k] == undef )
                  tmp_mask = np.logical_and( ref[i,:,k]   != undef , tmp_mask )
                  ref[i,:,k][ tmp_mask ] = undef
                  tmp_index[i,:,k][ tmp_mask ] = 1.0

                  tmp_mask = np.logical_and( ref[i-2,:,k] == undef , ref[0,:,k] == undef )
                  tmp_mask = np.logical_and( ref[i,:,k]   != undef , tmp_mask )
                  ref[i,:,k][ tmp_mask ] = undef
                  tmp_index[i,:,k][ tmp_mask ] = 1.0

               elif  i==0      :

                  tmp_mask = np.logical_and( ref[na-1,:,k] == undef , ref[i+1,:,k] == undef )
                  tmp_mask = np.logical_and( ref[i,:,k]    != undef , tmp_mask )
                  ref[i,:,k][ tmp_mask ] = undef
                  tmp_index[i,:,k][ tmp_mask ] = 1.0

                  tmp_mask = np.logical_and( ref[na-2,:,k] == undef , ref[i+2,:,k] == undef )
                  tmp_mask = np.logical_and( ref[i,:,k]    != undef , tmp_mask )
                  ref[i,:,k][ tmp_mask ] = undef
                  tmp_index[i,:,k][ tmp_mask ] = 1.0

               elif  i==1      :

                  tmp_mask = np.logical_and( ref[i-1,:,k] == undef , ref[i+1,:,k] == undef )
                  tmp_mask = np.logical_and( ref[i,:,k]    != undef , tmp_mask )
                  ref[i,:,k][ tmp_mask ] = undef
                  tmp_index[i,:,k][ tmp_mask ] = 1.0

                  tmp_mask = np.logical_and( ref[na-1,:,k] == undef , ref[i+2,:,k] == undef )
                  tmp_mask = np.logical_and( ref[i,:,k]    != undef , tmp_mask )
                  ref[i,:,k][ tmp_mask ] = undef
                  tmp_index[i,:,k][ tmp_mask ] = 1.0


            if ElevationFilter         :   #DETECT ISOLATED PIXELS IN ELEVATION

               if ( k > 0 ) & ( k < ne-1 ) :

                  tmp_mask = np.logical_and( ref[i,:,k-1] == undef , ref[i,:,k+1] == undef )
                  tmp_mask = np.logical_and( ref[i,:,k]    != undef , tmp_mask )
                  ref[i,:,k][ tmp_mask ] = undef
                  tmp_index[i,:,k][ tmp_mask ] = 1.0

               if ( k == 0 )                :

                  tmp_mask = np.logical_and( ref[i,:,k+2] == undef , ref[i,:,k+1] == undef )
                  tmp_mask = np.logical_and( ref[i,:,k]    != undef , tmp_mask )
                  ref[i,:,k][ tmp_mask ] = undef
                  tmp_index[i,:,k][ tmp_mask ] = 1.0

               if ( k == ne-1 )            :

                  tmp_mask = np.logical_and( ref[i,:,k-2] == undef , ref[i,:,k-1] == undef )
                  tmp_mask = np.logical_and( ref[i,:,k]    != undef , tmp_mask )
                  ref[i,:,k][ tmp_mask ] = undef
                  tmp_index[i,:,k][ tmp_mask ] = 1.0

   return tmp_index 


def smooth_fft( x )   :


    import numpy as np
    import scipy.fftpack

    N = np.size( x )

    w = scipy.fftpack.rfft( x )
    #f = scipy.fftpack.rfftfreq( N , 1.0 )
    spectrum = w**2

    cutoff_idx = spectrum > (spectrum.max()/10.0)
    #w2 = w.copy()
    w[cutoff_idx] = 0

    x_smooth = scipy.fftpack.irfft(w)


    return x_smooth


#Read topography data
def read_topo( my_file )           :
    import numpy as np

    f=open(my_file,'r')

    #First read the header containing nr and na
    tmpdata=np.fromfile(f,dtype='f4',count=2)
    na=int(tmpdata[0])
    nr=int(tmpdata[1])

    my_topo=dict()
    my_topo['mean']=np.reshape( np.fromfile(f,dtype='f4',count=nr*na) , (na,nr) )
    my_topo['max']=np.reshape( np.fromfile(f,dtype='f4',count=nr*na) , (na,nr) )
    my_topo['min']=np.reshape( np.fromfile(f,dtype='f4',count=nr*na) , (na,nr) )
    my_topo['number']=np.reshape( np.fromfile(f,dtype='f4',count=nr*na) , (na,nr) )
    my_topo['range']=np.reshape( np.fromfile(f,dtype='f4',count=nr*na) , (na,nr) )
    my_topo['azimuth']=np.reshape( np.fromfile(f,dtype='f4',count=nr*na) , (na,nr) )
    my_topo['latitude']=np.reshape( np.fromfile(f,dtype='f4',count=nr*na) , (na,nr) )
    my_topo['longitude']=np.reshape( np.fromfile(f,dtype='f4',count=nr*na) , (na,nr) )


    return my_topo

#Write topography data
def write_topo( my_topo , my_file )           :
    import numpy as np

    f=open(my_file,'w')

    na=my_topo['mean'].shape[0]
    nr=my_topo['mean'].shape[1]

    #First write nr and na
    np.array(na).astype('f4').tofile(f)
    np.array(nr).astype('f4').tofile(f)

    #Write the components of the my_topo dictionary.
    tmp=my_topo['range'].data
    np.reshape( my_topo['mean'].astype('f4') , (nr*na) ).tofile(f)
    np.reshape( my_topo['max'].astype('f4') , (nr*na) ).tofile(f)
    np.reshape( my_topo['min'].astype('f4') , (nr*na) ).tofile(f)
    np.reshape( my_topo['number'].astype('f4') , (nr*na) ).tofile(f)
    try     : 
       tmp=my_topo['range'].data
       np.reshape( tmp.astype('f4') , (nr*na)  ).tofile(f)
    except  :
       tmp=my_topo['range']
       np.reshape( tmp.astype('f4') , (nr*na)  ).tofile(f)
       

    #np.reshape( (my_topo['range'].data).astype('f4') , (nr*na) ).tofile(f)
    np.reshape( my_topo['azimuth'].astype('f4') , (nr*na) ).tofile(f)
    np.reshape( my_topo['latitude'].astype('f4') , (nr*na) ).tofile(f)
    np.reshape( my_topo['longitude'].astype('f4') , (nr*na) ).tofile(f)

    return my_topo

#Read a raster file an return lon,lat and data.
def read_raster_data(inputfile)  :
    import osgeo.gdal as gdal
    import numpy as np

    my_raster = gdal.Open(inputfile)

    nx=my_raster.RasterXSize
    ny=my_raster.RasterYSize

    nb=my_raster.RasterCount

    my_raster_data=np.zeros((nx,ny,nb))

    #Read the data and store it into a numpy array

    for ib in range( 0 , nb )  :

        my_raster_data[:,:,ib] = my_raster.ReadAsArray(ib)


    #Get the lat and lons.
    [lon,lat]=get_lat_lon(inputfile)

    return lon , lat , my_raster_data

#Read a raster file and get lat lon structure.
def get_lat_lon(inputfile)  :
    import osgeo.gdal as gdal
    import numpy as np

    my_raster = gdal.Open(inputfile)

    gt = my_raster.GetGeoTransform()
    proj = my_raster.GetProjection()

    xres = gt[1]
    yres = gt[5]

    # get the edge coordinates and add half the resolution 
    # to go to center coordinates
    xmin = gt[0] + xres * 0.5
    xmax = gt[0] + (xres * my_raster.RasterXSize) - xres * 0.5
    ymin = gt[3] + (yres * my_raster.RasterYSize) + yres * 0.5
    ymax = gt[3] - yres * 0.5

    lon=np.zeros( my_raster.RasterXSize )
    lat=np.zeros( my_raster.RasterYSize )
    for ii in range( 0 , my_raster.RasterXSize )  :
       lon[ii] = xmin + ii * xres
    for ii in range( 0 , my_raster.RasterYSize )  :
       lat[ii]= ymax + ii*yres

    my_raster = None

    # create a grid of xy coordinates in the original projection
    [lon,lat] = np.meshgrid(lon,lat)

    return lon , lat

#Download topography data and generate binary topography file.
def generate_topo_file( rlon , rlat , rrange , razimuth , raster_path , topo_file )    :
    import elevation
    import numpy as np
    import os
    from common_functions import common_functions as cf

    product='SRTM1'  #SRTM1 - 30 m res , SRTM3 - 90 m res.
         
    my_topo=dict()

    [ my_topo['range'] ,  my_topo['azimuth'] ] = np.meshgrid( rrange , razimuth )
    nr=np.size( rrange )
    na=np.size( razimuth )

    #Get lat and lon corresponding to the grid points in the polar coordinate.

    [my_topo['longitude'],my_topo['latitude']]=cf.com_ra_to_ll(cen_lon=rlon,cen_lat=rlat
                                                              ,r=my_topo['range'],a=my_topo['azimuth'],nr=nr,na=na)

    #Define grid limits
    maxlon = int( np.ceil( np.max(my_topo['longitude']) )  )
    minlon = int( np.floor( np.min(my_topo['longitude']) ) )
    maxlat = int( np.ceil( np.max(my_topo['latitude']) )  )
    minlat = int( np.floor( np.min(my_topo['latitude']) ) )

    array_size = np.shape( my_topo['range'] )

    #Download the data to a raster file.
    #west, south, east, north 

    my_topo['mean'] = np.zeros( array_size , order='F' , dtype=np.float32 )  
    my_topo['min']  = np.zeros( array_size , order='F' , dtype=np.float32 )
    my_topo['max']  = np.zeros( array_size , order='F' , dtype=np.float32 )
    my_topo['number']  = np.zeros( array_size , order='F' , dtype=np.int32 ) 

    #Data is downloaded and processed into 1 deg patches.
    for ilon in range( minlon , maxlon )  :
      for ilat in range( minlat , maxlat )  : 

         print('Downloading data for tile Lon=',ilon,' Lat=',ilat) 

         #If raster file is not present, then download it from the internet.
         raster_file = raster_path + '/' + product + str(ilon) + '_' + str(ilat) + '.tif'
         print('Downloading ' + raster_file )

         if ( not os.path.isfile( raster_file ) )  :
             elevation.clip(bounds=(ilon,ilat,ilon+1,ilat+1),output=raster_file,product=product)

         #Read data from raster file.
         print('Reading ' + raster_file ) 
         [raster_lon,raster_lat,raster_data]=read_raster_data(raster_file)
         raster_nx=raster_data.shape[0]
         raster_ny=raster_data.shape[1]

         #Convert raster lat lon to range and azimuth.
         [raster_r,raster_a]=cf.com_ll_to_ra(cen_lon=rlon,cen_lat=rlat
                                            ,lon=raster_lon,lat=raster_lat
                                            ,nx=raster_nx,ny=raster_ny)

         #Interpolate raster data to polar coordinates surrounding the radar.
         rmin=np.min( rrange )
         dr  =rrange[1] - rrange[0]
         amin=np.min( razimuth )
         da  =razimuth[1]  - razimuth[0]

         print('Interpolating patch Lon=',ilon,' Lat=',ilat) 
         cf.com_interp_boxavereg(xini=amin,dx=da,nx=na
                                ,yini=rmin,dy=dr,ny=nr
                                ,xin=raster_a.reshape( raster_nx * raster_ny )
                                ,yin=raster_r.reshape( raster_nx * raster_ny )
                                ,datain=raster_data.reshape( raster_nx * raster_ny )
                                ,nin=raster_nx*raster_ny
                                ,data_sum=my_topo['mean']
                                ,data_max=my_topo['max']
                                ,data_min=my_topo['min']
                                ,data_n=my_topo['number'],undef=-999)

          
    #Compute the mean topography.
    mask = my_topo['number'] > 0
    my_topo['mean'][ mask ] = my_topo['mean'][ mask ] / my_topo['number'][ mask ]

    #Complete missing values using neighborhood values.
    mask = my_topo['number'] == 0  #Identify missing data points
    my_topo['mean']=cf.com_complete_missing_2d(field=my_topo['mean'],missing_mask=mask,
                                               nx=na,ny=nr,npass=4)
    my_topo['max']=cf.com_complete_missing_2d(field=my_topo['max'],missing_mask=mask,
                                               nx=na,ny=nr,npass=4)
    my_topo['min']=cf.com_complete_missing_2d(field=my_topo['min'],missing_mask=mask,
                                               nx=na,ny=nr,npass=4)


    #Write topo file in polar coordintes into binary format.
    #This data will be used for fast interpolation to the corresponding volume elevation angles.
    write_topo( my_topo , topo_file )

    return my_topo








