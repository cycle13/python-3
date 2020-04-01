#print __doc__
# Author: Rapid Refresh Argentina Team
# License: BSD 3 clause

def main_qc( filename , options ) :

   import sys
   import time
   sys.path.append('../fortran')

   import numpy as np
   import numpy.ma as ma
   import matplotlib.pyplot as plt
   import pyart
   from common_qc_tools  import qc  #Fortran code routines.
   from scipy.interpolate import interp2d
   import netCDF4
   import os
   from netCDF4 import Dataset

   import matplotlib.pyplot as plt

   undef          = options['undef']   #This set the undef value for the rest of the script.

   output=dict() #Initialize output dictionary.

   computed_etfilter=False   #Flag to indicate if echo top has been computed already.

#  Constant parameters

   #Shortcut to variable names
   name_v   =options['name_v']
   name_ref =options['name_ref']
   name_rho =options['name_rho']

   name_cref=options['name_cref']
   name_cv  =options['name_cv']

   radar = pyart.io.read(filename)

   if options['is_rma']  :

      radar = get_rma_strat( filename , radar ) 

   if radar.altitude_agl['data'] == 0.0    :
      #Radar is assumed to be at ground level. Add 20 meters to account for the height of the tower.
      radar.altitude['data']=radar.altitude['data']+options['radar_altitude_agl']


   #Get the nyquist velocity
   if not options['is_rma']  :
      nyquistv=radar.get_nyquist_vel(0,check_uniform=True)
      

   output['maxw_ref']=0.0
   output['maxw_v']  =0.0

   #===================================================
   # RESHAPE VARIABLES
   #===================================================
   #From time,range -> azimuth,range,elevation

   startt=time.time()


   if name_ref in radar.fields :
        start=time.time()

        output['undef_ref']=radar.fields[name_ref]['_FillValue']

        [ output['ref'] , output['az'] , output['level'] , output['time'] , output['index'] , output['az_exact'] ]=order_variable( radar , name_ref , output['undef_ref'] )
        na=output['ref'].shape[0]
        nr=output['ref'].shape[1]
        ne=output['ref'].shape[2]
        output['na']=na
        output['nr']=nr
        output['ne']=ne

        output['cref'] = np.zeros(output['ref'].shape) 

        output['cref'][:] = output['ref']                 #Initialize the corrected reflectivity array.

        output['cref'][ output['cref'] == output['undef_ref'] ]=options['norainrefval']
        output['ref'] [ output['ref']  == output['undef_ref'] ]=options['norainrefval']

        output['qcref'] = np.zeros(output['cref'].shape)  #Set the qc flag array to 0.
        output['wref']  = np.zeros(output['cref'].shape)  #Set the weigths to 0.

        output['elevations']=np.unique(radar.elevation['data'])

        end=time.time()

        print("The elapsed time in {:s} is {:2f}".format("ref -> az,r,el",end-start) )
 
   if name_v in radar.fields  : 

        start=time.time()

        output['undef_v']=radar.fields[ name_v ]['_FillValue']

        [ output['v'] , output['az'] , output['level'] , output['time'] , output['index'] , output['az_exact']  ]=order_variable( radar , name_v , output['undef_v']  )
 
        na=output['ref'].shape[0]
        nr=output['ref'].shape[1]
        ne=output['ref'].shape[2]
        output['na']=na
        output['nr']=nr
        output['ne']=ne
 
        output['cv'] = np.zeros(output['v'].shape) 

        output['cv'][:] = output['v']                     #Initialize the corrected doppler velocity array

        output['qcv'] = np.zeros(output['v'].shape)       #Set the qc flag array to 0.

        output['wv'] = np.zeros(output['v'].shape)       #Set the weigths to 0.

        output['elevations']=np.unique(radar.elevation['data'])

        end=time.time()

        print("The elapsed time in {:s} is {:2f}".format("v -> az,r,el",end-start) )

   #===================================================
   # GEOREFERENCE RADAR DATA
   #===================================================

   #Use pyart rutines to compute x, y and z at each grid point
   #Reorder data
   start=time.time()

   #dm is a dummy variable

   [ output['altitude'] , dm , dm , dm , dm , dm ] = order_variable( radar , 'altitude' , undef )
   [ output['x']        , dm , dm , dm , dm , dm ] = order_variable( radar , 'x' , undef ) 
   [ output['y']        , dm , dm , dm , dm , dm ] = order_variable( radar , 'y' , undef )


   #Compute distance to radar for the first azimuth (we assume that distance to radar will be the same.
   #for all other azimuths.
   output['distance']=np.power( np.power(output['x'][0,:,:],2)+np.power(output['y'][0,:,:],2) , 0.5 )

   end=time.time()

   print("The elapsed time in {:s} is {:2f}".format("x,y,z -> az,r,el",end-start) )

   #===================================================
   # READ AND INTERPOLATE TOPOGRAPHY DATA
   #===================================================

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

   if ( os.path.isfile( polar_coord_topo_file ) )  :
      read_raster = False
      print('Using a previously generated topography file')
   else                                :
      read_raster = True
      print('Topography file not found. We will generate a new file from raw raster data')

   if read_raster    :   #We read the original data and interpolate it to a polar grid centered at the radar.

      #Generate topo file ( inputs are , radar lat and lon, range and azimuth )
      my_topo=generate_topo_file( radar.longitude['data'][0] , radar.latitude['data'][0] , radar.range['data'] , output['az'] , 
                                  options['toporawdatapath'] , polar_coord_topo_file )

   else   :

      print('Reading polar coordinate topography from a file')

      my_topo=read_topo( polar_coord_topo_file )

   #Interpolate topography data to the 3D radar structure.

   topo_r = my_topo['range'][0,:]
   topo_a = my_topo['azimuth'][:,0]

   #Define a "interpolator"
   interpolator = interp2d(topo_r,topo_a,my_topo['mean'], kind='linear')

   output['topo'] = np.zeros((na,nr,ne))

   #Loop over the vertical levels.
   for ie in range( 0 , ne )   :

      ppi_r = output['distance'][:,ie]
      ppi_a = output['az']

      output['topo'][:,:,ie] =  interpolator( ppi_r, ppi_a )

   end=time.time()

   #Compute AGL height for each pixel.
   output['altitude_agl']=output['altitude']-output['topo'] 

   print("The elapsed time in {:s} is {:2f}".format("topography interpolation",end-start) )

 
   #===================================================
   # DEALIASING 
   #===================================================

   filter_name='Dealiasing'
   if options[filter_name]['flag']  and ( name_v in radar.fields ) :
      start=time.time()

      if options[filter_name]['texture_filter']  :

         #Perform a texture filter before applying the dealiasing to the velocity field.
         #This is required to avoid applying the dealiasing algorithm to noisy data. 

         output['dv_texture']=qc.compute_texture(var=output['v'],na=na,nr=nr,ne=ne,
                                                 undef=output['undef_v'],
                                                 nx=options[filter_name]['nx'],
                                                 ny=options[filter_name]['ny'],
                                                 nz=options[filter_name]['nz'])

         output['cv'][ output['dv_texture'] > options[filter_name]['texture_thr'] ] = output['undef_v']
         output['qcv'][ output['dv_texture'] > options[filter_name]['texture_thr'] ] = options[filter_name]['texture_code']

         tmp = order_variable_inv(  radar , output['cv'] , output['index'] , output['undef_v'] )

         radar.fields[ name_cv ] = dict()

         radar.fields[ name_cv ] = radar.fields[ name_v ]

         radar.fields[name_cv]['data']= np.ma.masked_array( tmp , tmp==output['undef_v'] ) 

      else   :

          radar.fields[ name_cv ] = dict()

          radar.fields[ name_cv ] = radar.fields[ name_v ]

          tmp = order_variable_inv(  radar , output['cv'] , output['index'] , output['undef_v'] )

          radar.fields[name_cv]['data']= np.ma.masked_array( tmp , tmp==output['undef_v'] )

      #Uso una de las funciones de dealiasing de pyart con los parametros por defecto
      winddealias=pyart.correct.region_dealias.dealias_region_based(radar,interval_splits=options[filter_name]['interval_split'],interval_limits=None, 
                 skip_between_rays=options[filter_name]['skip_between_ray'],skip_along_ray=options[filter_name]['skip_along_ray'],centered=True,nyquist_vel=None,
                 check_nyquist_uniform=True,gatefilter=False,rays_wrap_around=True,keep_original=False,set_limits=True,
                 vel_field=name_cv,corr_vel_field=None)

      #Replace cv wind by dealiased winds.
      radar.fields[ name_cv ]['data'] = winddealias['data']

      #Re-order dealiased wind data.
      [ output['cv'] , output['az'] , output['level'] , output['time'] , output['index'] , output['az_exact']  ]=order_variable( radar , name_cv , output['undef_v']  )
    
      mask=np.logical_and( output['cv'] != output['v'] , output['cv'] != output['undef_v'] )
      output['qcv'][ mask ]=options[filter_name]['code']

      end=time.time()

      print("The elapsed time in {:s} is {:2f}".format("dealiasing",end-start) )


   #===================================================
   # MODEL FILTER 
   #===================================================

     #TODO
     #TODO
      

   #===================================================
   # DEALIASING BORDER FILTER
   #===================================================

   filter_name='DealiasingBorderFilter'  

   if options[filter_name]['flag']  and options['ifdealias'] :

     #Use a simple edge detection routine to detect the borders between dealiased and non-dealiased regions.
     #Flag this grid points and eliminate all the non-dealiased corrected pixels near by.

     tmp_diff=output['cv'] - output['v']

     mask = np.logical_or( output['cv'] == output['undef_v'] , output['v'] == output['undef_v'] )
     tmp_diff[ mask ] = output['undef_v']
    
     v_nyquist = radar.instrument_parameters['nyquist_velocity']['data'].max()
 
     if v_nyquist < 0 :
        v_nyquist = 40.0 #Put some value that can detect strong jumps in the doppler velocity field.

     [ output['edge_intensity'] , output['edge_mask'] ]=qc.simple_edge_filter( field=tmp_diff , nx=na,ny=nr,nz=ne,undef=output['undef_v'],
                                                                               nboxx=options[filter_name]['nx'],nboxy=options[filter_name]['ny'],
                                                                               nboxz=options[filter_name]['nz'],edge_tr=v_nyquist )

     #Find the pixels which are close to a dealiased region but which has not been corrected by the dealiasing.
     tmp_w = np.logical_and( output['edge_mask'] , tmp_diff == 0 ).astype(int) 

     if not options[filter_name]['force']   :   
       output['wv']=output['wv'] + tmp_w * options[filter_name]['w']
       output['qcv'][ tmp_w > 0.5 ] = options[filter_name]['code']
       output['maxw_v']=output['maxw_v'] + options[filter_name]['w']
     else                                   :
       output['cv'][ tmp_w > options[filter_name]['force_value'] ]=output['undef_v']
       output['qcv'][ tmp_w > options[filter_name]['force_value'] ] = options[filter_name]['code']

     #If requested store the auxiliary fields and data in the output dictionary.
     if  ( not options[filter_name]['save'] )     :
        output.pop('edge_intensity')
        output.pop('edge_mask')

     end=time.time()

   #===================================================
   # ECHO TOP FILTER  
   #===================================================
   filter_name='EchoTopFilter'
   if options[filter_name]['flag'] & ( name_ref in radar.fields ) :
     start=time.time()

     if ( not computed_etfilter )     :
     #Compute the filter field.

        nx=options[filter_name]['nx']
        ny=options[filter_name]['ny']
        nz=options[filter_name]['nz']
                                   
        tmp_max_z=np.zeros((na,nr,ne))

        for ii in range(0,output['ne'])     :       #Estimate the maximum radar data height assoicated with each gate.
           tmp_max_z[:,:,ii]=output['altitude'][:,:,output['ne']-1]
 
        [tmp_data_3d,tmp_data_2d]=qc.echo_top(reflectivity=output['ref'],heigth=output['altitude'][0,:,:]
                                                ,rrange=output['distance'],na=na,nr=nr,ne=ne
                                                ,undef=output['undef_ref'],nx=nx,ny=ny,nz=nz)

        output['echo_top']=tmp_data_3d[:,:,:,0]
        output['echo_depth']=tmp_data_3d[:,:,:,2]

        computed_etfilter = True  #In case we need any of the other variables computed in this routine.

     #Compute the corresponding weigth.
     tmp_w = qc.multiple_1d_interpolation( field=output['echo_top'] , nx=na , ny=nr , nz=ne
                                           , undef=output['undef_ref'] , xx=options[filter_name]['ifx']
                                           , yy=options[filter_name]['ify'] , nxx=np.size(options[filter_name]['ifx']) )

     tmp_w[ tmp_max_z < options[filter_name]['heigthtr'] ] = 0.0  #Do not consider this filter when the volume maximum heigth is below
                                                                  #the specified threshold (i.e. pixels close to the radar)

     if not options[filter_name]['force']   :
       output['wref']=output['wref'] + tmp_w * options[filter_name]['w']
       output['qcref'][ tmp_w > 0.5 ] = options[filter_name]['code']
       output['maxw_ref']=output['maxw_ref'] + options[filter_name]['w']
     else                                   :
       output['cref'][ tmp_w > options[filter_name]['force_value'] ]=output['undef_ref']
       output['qcref'][ tmp_w > options[filter_name]['force_value'] ] = options[filter_name]['code']


     #If requested store the auxiliary fields and data in the output dictionary.
     if  ( not options[filter_name]['save'] )     : 
        output.pop('echo_top') 
        computed_etfilter = False

     end=time.time()

     print("The elapsed time in {:s} is {:2f}".format(filter_name,end-start) )

   #===================================================
   # ECHO DEPTH FILTER 
   #===================================================
   filter_name='EchoDepthFilter'
   if options[filter_name]['flag'] & ( name_ref in radar.fields ) :
     start=time.time()

     if ( not computed_etfilter )     :
     #Compute the filter field.

        nx=options[filter_name]['nx']
        ny=options[filter_name]['ny']
        nz=options[filter_name]['nz']

        tmp_max_z=np.zeros((na,nr,ne))

        for ii in range(0,output['ne'])     :       #Estimate the maximum radar data height assoicated with each gate.
           tmp_max_z[:,:,ii]=output['altitude'][:,:,output['ne']-1]

        [tmp_data_3d,tmp_data_2d]=qc.echo_top(reflectivity=output['ref'],heigth=output['altitude'][0,:,:]
                                                ,rrange=output['distance'],na=na,nr=nr,ne=ne
                                                ,undef=output['undef_ref'],nx=nx,ny=ny,nz=nz)

        output['echo_top']=tmp_data_3d[:,:,:,0]
        output['echo_depth']=tmp_data_3d[:,:,:,2]

        computed_etfilter = True  #In case we need any of the other variables computed in this routine.


     #Compute the corresponding weigth.
     tmp_w = qc.multiple_1d_interpolation( field=output['echo_depth'] , nx=na , ny=nr , nz=ne
                                           , undef=output['undef_ref'] , xx=options[filter_name]['ifx']
                                           , yy=options[filter_name]['ify'] , nxx=np.size(options[filter_name]['ifx']) )

     tmp_w[ tmp_max_z < options[filter_name]['heigthtr'] ] = 0.0  #Do not consider this filter when the volume maximum heigth is below
                                                                  #the specified threshold (i.e. pixels close to the radar)

     if not options[filter_name]['force']   :
       output['wref']=output['wref'] + tmp_w * options[filter_name]['w']
       output['qcref'][ tmp_w > 0.5 ] = options[filter_name]['code']
       output['maxw_ref']=output['maxw_ref'] + options[filter_name]['w']
     else                                   :
       output['cref'][ tmp_w > options[filter_name]['force_value'] ]=output['undef_ref']
       output['qcref'][ tmp_w > options[filter_name]['force_value'] ] = options[filter_name]['code']


     #If requested store the auxiliary fields and data in the output dictionary.
     if  ( not options[filter_name]['save'] )     :
        output.pop('echo_top')
        computed_etfilter = False

     end=time.time()

     print("The elapsed time in {:s} is {:2f}".format(filter_name,end-start) )

   #===================================================
   # LOW ELEVATION ANGLE REFLECTIVITY FILTER
   #===================================================

   #Remove echos which are present at low elevation angles but not at higher elevation
   #angles. This will help to remove clutter produced by anomalous propagation.
   #This can also help to remove second trip echoes which are usually observed only at low elevation
   #angles.
   #This can also eliminate distant convective cells that are only seen in low 
   #elevation angles near the radar range edge. However these cells are of little interest
   #for data assimilation since that information is usually not enough to adequatelly constrain
   #the evolution of the convective cells.

   filter_name='LowElevFilter'

   if options[filter_name]['flag']  & ( name_ref in radar.fields ) :

      start=time.time()

      nx=options[filter_name]['nx']
      ny=options[filter_name]['ny']
      nz=options[filter_name]['nz']

      #Get the angles that will be used based on the selected threshold.
      tmp_angles= output['elevations'][ output['elevations'] < options[filter_name]['min_angle']]
      tmp_n_angles = np.size( tmp_angles )

      output['smooth_ref']=qc.box_functions_2d(datain=output['ref'],na=na,nr=nr,ne=ne,undef=output['undef_ref']
                                               ,boxx=nx,boxy=ny,boxz=nz,operation='MEAN',threshold=0.0)
     
      tmp_w=np.zeros([na,nr,ne]) 
      for ie in range( 0 , tmp_n_angles )  :
         tmp_w[:,:,ie]=np.logical_and( output['ref'][:,:,ie] > options['norainrefval'] , output['smooth_ref'][:,:,tmp_n_angles] <= options['norainrefval'] )
         tmp_w[:,:,ie][ output['altitude'][:,:,tmp_n_angles] > options[filter_name]['height_thr'] ] = 0.0

      if not options[filter_name]['force']   :
         output['wref']=output['wref'] + tmp_w * options[filter_name]['w']
         output['qcref'][ tmp_w > 0.5 ] = options[filter_name]['code']
         output['maxw_ref']=output['maxw_ref'] + options[filter_name]['w']
      else                                   :
         output['cref'][ tmp_w > options[filter_name]['force_value'] ]=output['undef_ref']
         output['qcref'][ tmp_w > options[filter_name]['force_value'] ] = options[filter_name]['code']

      #If requested store the auxiliary fields and data in the output dictionary.
      if  ( not options[filter_name]['save'] )     :
          output.pop('smooth_ref')
          computed_etfilter = False


      end=time.time()
      print("The elapsed time in {:s} is {:2f}".format(filter_name,end-start) )

   #===================================================
   # RHO HV FILTER
   #===================================================

   filter_name = 'RhoFilter'

   if options[filter_name]['flag']  :

      start=time.time()

      nx=options[filter_name]['nx']
      ny=options[filter_name]['ny']
      nz=options[filter_name]['nz']

      if name_rho in radar.fields  :

         output['undef_rho']=radar.fields[ name_rho ]['_FillValue']

         [ output['rho'] , dm , dm , dm , dm , dm  ]=order_variable( radar , name_rho , output['undef_rho']  )

         #Compute the filter parameter
         output['rho_smooth']=qc.box_functions_2d(datain=output['rho'],na=na,nr=nr,ne=ne,undef=output['undef_rho']
                                                 ,boxx=nx,boxy=ny,boxz=nz,operation='MEAN',threshold=0.0)


         #Compute the corresponding weigth.
         tmp_w = qc.multiple_1d_interpolation( field=output['rho_smooth'] , nx=na , ny=nr , nz=ne 
                                              , undef=output['undef_rho'] , xx=options[filter_name]['ifx'] 
                                              , yy=options[filter_name]['ify'] , nxx=np.size(options[filter_name]['ifx']) )

         output['wref']=output['wref'] + tmp_w * options[filter_name]['w']

         if name_ref in radar.fields :
            if not options[filter_name]['force']   :
               output['wref']=output['wref'] + tmp_w * options[filter_name]['w']
               output['qcref'][ tmp_w > 0.5 ] = options[filter_name]['code']
               output['maxw_ref']=output['maxw_ref'] + options[filter_name]['w']
            else                                   :
               output['cref'][ tmp_w > options[filter_name]['force_value'] ]=output['undef_ref']
               output['qcref'][ tmp_w > options[filter_name]['force_value'] ] = options[filter_name]['code']

      else   :
         display('Warning: could not perform RHO-HV filter because rho was not found on this file')

      print( options[filter_name]['save'] )
      if ( not options[filter_name]['save'] ) :
          output.pop('rho_smooth')
          output.pop('rho')

      end=time.time()

      print("The elapsed time in {:s} is {:2f}".format(filter_name,end-start) )

 
   #===================================================
   # REFLECTIVITY SPECKLE FILTER
   #===================================================

   filter_name='RefSpeckleFilter'

   if options[filter_name]['flag']  & ( name_ref in radar.fields ) :

       start=time.time()


       #Compute the number pixels with reflectivities over spfiltertr sourrounding each pixels in the box defined by nx,ny,nz.
       nx=options[filter_name]['nx']
       ny=options[filter_name]['ny']
       nz=options[filter_name]['nz']
       tr=options[filter_name]['reftr']

       output['speckle_ref']=qc.box_functions_2d(datain=output['ref'].data,na=na,nr=nr,ne=ne,undef=output['undef_ref']
                                                ,boxx=nx,boxy=ny,boxz=nz,operation='COU2',threshold=tr) 


       #Compute the corresponding weigth.
       tmp_w = qc.multiple_1d_interpolation( field=output['speckle_ref'] , nx=na , ny=nr , nz=ne
                                            , undef=output['undef_ref'] , xx=options[filter_name]['ifx']
                                            , yy=options[filter_name]['ify'] , nxx=np.size(options[filter_name]['ifx']) )

       tmp_w[ output['ref'] == output['undef_ref'] ] = 0.0 #If a grid point is already undef leave it as it is.


       if not options[filter_name]['force']   :
          output['wref']=output['wref'] + tmp_w * options[filter_name]['w']
          output['qcref'][ tmp_w > 0.5 ] = options[filter_name]['code']
          output['maxw_ref']=output['maxw_ref'] + options[filter_name]['w']
       else                                   :
          output['cref'][ tmp_w > options[filter_name]['force_value'] ]=output['undef_ref']
          output['qcref'][ tmp_w > options[filter_name]['force_value'] ] = options[filter_name]['code']



       if ( not options[filter_name]['save'] ) :
          output.pop('speckle_ref')

       end=time.time()

       print("The elapsed time in {:s} is {:2f}".format(filter_name,end-start) )


   #===================================================
   # DOPPLER SPECKLE FILTER
   #===================================================

   filter_name='DopplerSpeckleFilter'

   if options[filter_name]['flag'] & ( name_v in radar.fields ) :

       start=time.time()


       #Compute the number pixels with reflectivities over spfiltertr sourrounding each pixels in the box defined by nx,ny,nz.
       nx=options[filter_name]['nx']
       ny=options[filter_name]['ny']
       nz=options[filter_name]['nz']
       tr=options[filter_name]['dvtr']


       tmp=np.abs(output['v'].data)

       tmp[ output['v'] == output['undef_v'] ] = output['undef_v']

       output['speckle_v']=qc.box_functions_2d(datain=tmp,na=na,nr=nr,ne=ne,undef=output['undef_v']
                                                ,boxx=nx,boxy=ny,boxz=nz,operation='COU2',threshold=tr)


       #Compute the corresponding weigth.
       tmp_w = qc.multiple_1d_interpolation( field=output['speckle_v'] , nx=na , ny=nr , nz=ne
                                            , undef=output['undef_v'] , xx=options[filter_name]['ifx']
                                            , yy=options[filter_name]['ify'] , nxx=np.size(options[filter_name]['ifx']) )

       tmp_w[ output['v'] == output['undef_v'] ]=0.0 #Pixels which are already flagged as undef should remain undef.

       if not options[filter_name]['force']   :
          output['wv']=output['wv'] + tmp_w * options[filter_name]['w']
          output['qcv'][ tmp_w > 0.5 ] = options[filter_name]['code']
          output['maxw_v']=output['maxw_v'] + options[filter_name]['w']
       else                                   :
          output['cv'][ tmp_w > options[filter_name]['force_value'] ]=output['undef_v']
          output['qcv'][ tmp_w > options[filter_name]['force_value'] ] = options[filter_name]['code']

       if ( not options[filter_name]['save'] ) :
          output.pop('speckle_v')

       end=time.time()

       print("The elapsed time in {:s} is {:2f}".format(filter_name,end-start) )

   #===================================================
   # ATTENUATION FILTER
   #===================================================

   filter_name='AttenuationFilter'

   if options[filter_name]['flag'] & ( name_ref in radar.fields ) :
      
      start=time.time()

      beaml=radar.range['data'][1]-radar.range['data'][0] #Get beam length

      output['attenuation']=qc.get_attenuation( var=output['cref'],na=na,nr=nr,ne=ne,undef=output['undef_ref']
                                                ,beaml=beaml,cal_error=options[filter_name]['attcalerror']
                                                ,is_power=options[filter_name]['is_power']
                                                ,coefs=options[filter_name]['att_coefs'] )

      #Compute the corresponding weigth.
      tmp_w = qc.multiple_1d_interpolation( field=output['attenuation'] , nx=na , ny=nr , nz=ne
                                            , undef=output['undef_ref'] , xx=options[filter_name]['ifx']
                                            , yy=options[filter_name]['ify'] , nxx=np.size(options[filter_name]['ifx']) )

      if not options[filter_name]['force']   :
         output['wref']=output['wref'] + tmp_w * options[filter_name]['w']
         output['qcref'][ tmp_w > 0.5 ] = options[filter_name]['code']
         output['maxw_ref']=output['maxw_ref'] + options[filter_name]['w']
      else                                   :
         output['cref'][ tmp_w > options[filter_name]['force_value'] ]=output['undef_ref']
         output['qcref'][ tmp_w > options[filter_name]['force_value'] ] = options[filter_name]['code']

      if  ( not options[filter_name]['save'] )  :
          output.pop('attenuation')

      end=time.time()

      print("The elapsed time in {:s} is {:2f}".format(filter_name,end-start) )

   #===================================================
   # TOPOGRAPHY BLOCKING FILTER
   #===================================================

   filter_name='BlockingFilter'  #This is not included in the fuzzy logic algorithm

   if options[filter_name]['flag']   :

      start=time.time()

      output['blocking']=qc.compute_blocking( radarz=output['altitude'] , topo=output['topo'] , na=na , nr=nr , ne=ne      ,
                                              undef=output['undef_ref']                                                    ,
                                              radar_beam_width_v=radar.instrument_parameters['radar_beam_width_v']['data'] , 
                                              beam_length=radar.range['meters_between_gates']                              , 
                                              radarrange=radar.range['data'] , radarelev=output['elevations'] )  

      #Compute correction 
      if options[filter_name]['blocking_correction'] & ( name_ref in radar.fields )  :
         #Correct partially blocked precipitation echoes.
         mask=np.logical_and( output['blocking'] > 0.1 , output['blocking'] <= 0.3 )
         mask=np.logical_and( mask , output['cref'] > options['norainrefval'] )

         output['cref'][mask] = output['cref'][mask] + 1.0

         mask=np.logical_and( output['blocking'] > 0.3 , output['blocking'] <= 0.4 )
         mask=np.logical_and( mask , output['cref'] > options['norainrefval'] )

         output['cref'][mask] = output['cref'][mask] + 2.0

         mask= output['blocking'] > 0.4 
         mask=np.logical_and( mask , output['cref'] > options['norainrefval'] )

         output['cref'][mask] = output['cref'][mask] + 3.0

      #Set the pixels with values below the threshold as undef. 
      if name_ref in radar.fields :
          output['cref'][ output['blocking'] > options[filter_name]['blocking_threshold']  ] = output['undef_ref']
          output['qcref'][ output['blocking'] > options[filter_name]['blocking_threshold'] ] = options[filter_name]['code']

      if name_v in radar.fields :
          output['cv'][ output['blocking']   > options[filter_name]['blocking_threshold']  ] = output['undef_v']
          output['qcv']  [ output['blocking'] > options[filter_name]['blocking_threshold'] ] = options[filter_name]['code']

      if  ( not options[filter_name]['save'] ) :
          output.pop('blocking')  

      end=time.time()

      print("The elapsed time in {:s} is {:2f}".format(filter_name,end-start) )

   #===================================================
   # DOPPLER NOISE FILTER
   #===================================================

   filter_name='DopplerNoiseFilter'

   if options[filter_name]['flag'] & ( name_v in radar.fields ) :

     start=time.time()

     nx=options[filter_name]['nx']
     ny=options[filter_name]['ny']
     nz=options[filter_name]['nz']

     nx2=options[filter_name]['nx2']
     ny2=options[filter_name]['ny2']
     nz2=options[filter_name]['nz2']

     tr_1=options[filter_name]['threshold_1']
     tr_2=options[filter_name]['threshold_2']

     tmp_dv_1=np.copy(output['v'])
     tmp_dv_2=np.copy(output['v'])

     for ip in range(0,options[filter_name]['n_filter_pass']) :

       output['distance_1']=qc.compute_distance(tmp_dv_1,tmp_dv_2,nx=na,ny=nr,nz=ne,
                                                undef=output['undef_v'],nx_box=nx,ny_box=ny,nz_box=nz)

       tmp_dv_2=np.copy(output['v']) 

       tmp_dv_2[ output['distance_1'] > tr_1 ]=output['undef_v']

     #Compute the corresponding weigth.
     tmp_w = qc.multiple_1d_interpolation( field=output['distance_1'], nx=na , ny=nr , nz=ne
                                           , undef=output['undef_v'] , xx=options[filter_name]['ifx_1']
                                           , yy=options[filter_name]['ify_1'] , nxx=np.size(options[filter_name]['ifx_1']) )


     if not options[filter_name]['force']   :
       output['wv']=output['wv'] + tmp_w * options[filter_name]['w']
       output['qcv'][ tmp_w > 0.5 ] = options[filter_name]['code']
     else                                   :
       output['cv'][ tmp_w > options[filter_name]['force_value'] ]=output['undef_v']
       output['qcv'][ tmp_w > options[filter_name]['force_value'] ] = options[filter_name]['code']


     tmp_dv_1[ np.logical_or( output['distance_1'] > tr_1 , output['distance_1']==output['undef_v'] ) ]=output['undef_v']

     tmp_dv_2=np.copy(tmp_dv_1)

     for ip in range(0,options[filter_name]['n_filter_pass']) :

       output['distance_2']=qc.compute_distance(tmp_dv_1,tmp_dv_2,nx=na,ny=nr,nz=ne,
                                                undef=output['undef_v'],nx_box=nx2,ny_box=ny2,nz_box=nz2)

       tmp_dv_2=np.copy(output['v'])

       tmp_dv_2[ output['distance_2'] > tr_2 ]=output['undef_v']


     #Compute the corresponding weigth.
     tmp_w = qc.multiple_1d_interpolation( field=output['distance_2'] , nx=na , ny=nr , nz=ne
                                           , undef=output['undef_v'] , xx=options[filter_name]['ifx_2']
                                           , yy=options[filter_name]['ify_2'] , nxx=np.size(options[filter_name]['ifx_2']) )

     if not options[filter_name]['force']   :
       output['wv']=output['wv'] + tmp_w * options[filter_name]['w']
       output['qcv'][ tmp_w > 0.5 ] = options[filter_name]['code']
       output['maxw_ref']=output['maxw_ref'] + options[filter_name]['w']
     else                                   :
       output['cv'][ tmp_w > options[filter_name]['force_value'] ]=output['undef_v']
       output['qcv'][ tmp_w > options[filter_name]['force_value'] ] = options[filter_name]['code']

     if  not options[filter_name]['save']  :
          output.pop('distance_1')
          output.pop('distance_2')

     end=time.time()


     print("The elapsed time in {:s} is {:2f}".format(filter_name,end-start) )

   #===================================================
   # INTERFERENCE FILTER
   #===================================================

   #Clean rays affected by interference.

   filter_name = 'InterferenceFilter'

   if options[filter_name]['flag'] & ( name_ref in radar.fields ) :

      start=time.time()

      tmp_ref=np.copy( output['cref'] )

      nx=options[filter_name]['nx']
      ny=options[filter_name]['ny']
      nz=options[filter_name]['nz']


      #output['smooth_ref']=qc.box_functions_2d(datain=output['ref'],na=na,nr=nr,ne=ne,undef=output['undef_ref']
      #                                                             ,boxx=nx,boxy=ny,boxz=nz,operation='MEAN',threshold=0.0)


      output['cref'] = interference_filter ( output['ref'] , output['undef_ref'] , options['norainrefval'] 
                                            , radar.range['data'] , options[filter_name] ) 

      output['qcref'][ np.abs( tmp_ref - output['cref'] ) > 0 ]=options[filter_name]['code']



      end=time.time()
 
      print("The elapsed time in {:s} is {:2f}".format(filter_name,end-start) )  


   #===================================================
   # POWER FILTER
   #===================================================

   #TODO

   #===================================================
   # DETECT MISSING REFLECTIVITY VALUES
   #===================================================

   filter_name='MissingRefFilter'

   if options[filter_name]['flag'] & ( name_ref in radar.fields ) :

      start=time.time()

      options['missing_mask'] = qc.detect_missing(  output['ref'],na=na,nr=nr,ne=ne,undef=output['undef_ref']
                                                   ,min_ref=options['norainrefval'],threshold=options[filter_name]['threshold'] 
                                                   ,nmissing_max=options[filter_name]['nmissing_max'] )

      tmp_w = options['missing_mask'].astype(int)

      if not options[filter_name]['force']   :
         output['wref']=output['wref'] + tmp_w * options[filter_name]['w']
         output['qcref'][ tmp_w > 0.5 ] = options[filter_name]['code']
         output['maxw_ref']=output['maxw_ref'] + options[filter_name]['w']
      else                                   :
         output['cref'][ tmp_w > options[filter_name]['force_value'] ]=output['undef_ref']
         output['qcref'][ tmp_w > options[filter_name]['force_value'] ] = options[filter_name]['code']

      if  not options[filter_name]['save']  :
           output.pop('missing_mask')

      end=time.time()

      print("The elapsed time in {:s} is {:2f}".format(filter_name,end-start) )

   #===================================================
   # REFLECTIIVTY TEXTURE FILTER
   #===================================================

   filter_name='ReflectivityTextureFilter'

   if  options[filter_name]['flag'] & ( name_ref in radar.fields ) : 

     nx=options[filter_name]['nx']
     ny=options[filter_name]['ny']
     nz=options[filter_name]['nz']
     output['texture_ref']=qc.compute_texture(var=output['ref'],na=na,nr=nr,ne=ne,undef=output['undef_ref'],nx=nx,ny=ny,nz=nz)


     if options[filter_name]['use_smooth_ref'] :
        output['smooth_ref']=qc.box_functions_2d(datain=output['ref'],na=na,nr=nr,ne=ne,undef=output['undef_ref'],
                                                 boxx=0,boxy=0,boxz=0,operation='MEAN',threshold=0.0)

        #High reflectivity cores will not be affected by texture filter.
        output['texture_ref'][ output['smooth_ref'] >= options[filter_name]['smooth_ref_tr']  ]= 0.0

     output['texture_ref'][output['ref']==output['undef_ref']] = 0.0

     #Compute the corresponding weigth.
     tmp_w = qc.multiple_1d_interpolation( field=output['texture_ref'] , nx=na , ny=nr , nz=ne
                                           , undef=output['undef_ref'] , xx=options[filter_name]['ifx']
                                           , yy=options[filter_name]['ify'] , nxx=np.size(options[filter_name]['ifx']) )


     if not options[filter_name]['force']   :
        output['wref']=output['wref'] + tmp_w * options[filter_name]['w']
        output['qcref'][ tmp_w > 0.5 ] = options[filter_name]['code']
        output['maxw_ref']=output['maxw_ref'] + options[filter_name]['w']
     else                                   :
        output['cref'][ tmp_w > options[filter_name]['force_value'] ]=output['undef_ref']
        output['qcref'][ tmp_w > options[filter_name]['force_value'] ] = options[filter_name]['code']

     if  not options[filter_name]['save']  :
          output.pop('texture_v')

     end=time.time()


     print("The elapsed time in {:s} is {:2f}".format(filter_name,end-start) )

   #===================================================
   # DOPPER VELOCITY TEXTURE FILTER
   #===================================================

   filter_name='DopplerTextureFilter'

   if  options[filter_name]['flag'] & ( name_ref in radar.fields ) :

     nx=options[filter_name]['nx']
     ny=options[filter_name]['ny']
     nz=options[filter_name]['nz']
     output['texture_v']=qc.compute_texture(var=output['v'],na=na,nr=nr,ne=ne,undef=output['undef_v'],nx=nx,ny=ny,nz=nz)

     #Compute the corresponding weigth.
     tmp_w = qc.multiple_1d_interpolation( field=output['texture_v'] , nx=na , ny=nr , nz=ne
                                           , undef=output['undef_v'] , xx=options[filter_name]['ifx']
                                           , yy=options[filter_name]['ify'] , nxx=np.size(options[filter_name]['ifx']) )

     if not options[filter_name]['force']   :
        output['wv']=output['wv'] + tmp_w * options[filter_name]['w']
        output['qcv'][ tmp_w > 0.5 ] = options[filter_name]['code']
        output['maxw_v']=output['maxw_v'] + options[filter_name]['w']
     else                                   :
        output['cv'][ tmp_w > options[filter_name]['force_value'] ]=output['undef_v']
        output['qcv'][ tmp_w > options[filter_name]['force_value'] ] = options[filter_name]['code']


     if  not options[filter_name]['save']  :
          output.pop('texture_v')

     end=time.time()

     print("The elapsed time in {:s} is {:2f}".format(filter_name,end-start) )

   #===================================================
   # LOW DOPPLER VOLOCITY FILTER
   #===================================================
 
   #This filter removes reflectivity values which are close to the surface and 
   #have a low associated doppler velocity. 

   filter_name='LowDopplerFilter'

   if options[filter_name]['flag'] & ( name_v in radar.fields ) :


      #Compute the corresponding weigth.
      tmp_w = qc.multiple_1d_interpolation( field=output['v'] , nx=na , ny=nr , nz=ne
                                            , undef=output['undef_v'] , xx=options[filter_name]['ifx']
                                            , yy=options[filter_name]['ify'] , nxx=np.size(options[filter_name]['ifx']) )

      #The filter will not produce any effect below a certain height
      if  options[filter_name]['use_terrain']    :
          tmp_w[ output['altitude_agl'] >  options[filter_name]['height_thr'] ] = 0.0
      else                                    :
          tmp_w[ output['altitude'] > options[filter_name]['height_thr'] ] = 0.0

      if name_ref in radar.fields   :
         if not options[filter_name]['force']   :
         
            output['wref']=output['wref'] + tmp_w * options[filter_name]['w']
            output['qcref'][ tmp_w > 0.5 ] = options[filter_name]['code']
            output['maxw_ref']=output['maxw_ref'] + options[filter_name]['w']
         else                                   :
            output['cref'][ tmp_w > options[filter_name]['force_value'] ]=output['undef_ref']
            output['qcref'][ tmp_w > options[filter_name]['force_value'] ] = options[filter_name]['code']

      if name_v in radar.fields   :
         if not options[filter_name]['force']   :

            output['wv']=output['wv'] + tmp_w * options[filter_name]['w']
            output['qcv'][ tmp_w > 0.5 ] = options[filter_name]['code']
            output['maxw_v']=output['maxw_v'] + options[filter_name]['w']
         else                                   :
            output['cv'][ tmp_w > options[filter_name]['force_value'] ]=output['undef_v']
            output['qcv'][ tmp_w > options[filter_name]['force_value'] ] = options[filter_name]['code']


      end=time.time()

      print("The elapsed time in {:s} is {:2f}".format(filter_name,end-start) )

   #===================================================
   # COMPUTE THE FINAL WEIGHT AND APPLY FUZZY LOGIC QC
   #===================================================
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

   #===================================================
   # ADD CORRECTED DATA TO RADAR OBJECT
   #===================================================

   if name_ref in radar.fields :

      tmp=order_variable_inv( radar , output['cref'] , output['index'] , output['undef_ref'] )

      radar.fields[ name_cref ] = dict()

      radar.fields[ name_cref ] = radar.fields[ name_ref ] 

      radar.fields[ name_cref ]['data']=np.ma.masked_array(tmp , tmp==output['undef_ref'] )

   if name_v in radar.fields :

      tmp=order_variable_inv( radar , output['cv'] , output['index'] , output['undef_v'] )

      radar.fields[ name_cv ]= dict()

      radar.fields[ name_cv ]= radar.fields[ name_v ] 

      radar.fields[ name_cv ]['data']=np.ma.masked_array(tmp , tmp==output['undef_v'] )


   #===================================================
   # END
   #===================================================


   endt=time.time()

   print("The elapsed time in {:s} is {:2f}".format("the entire QC",endt-startt) )

   return radar , output

#===========================================================================================================
# OTRAS FUNCIONES CONTENIDAS EN ESTE MODULO
#===========================================================================================================   

#From Pyart order to ARE (Azmimuth , range , elevation )

def order_variable ( radar , var_name , undef )  :  

   import numpy as np
   import numpy.ma as ma

   #order_var es la variable ordenada con los azimuths entre 0 y 360 (si hay rayos repetidos se promedian).
   #order_azimuth es el azimuth "aproximado" utilizando 0 como azimuth inicial y avanzando en intervalos regulares e iguales a la resolucion
   #del azimuth en grados.
   #levels son los angulos de elevacion.
   #azimuth_exact es un array que contiene para cada nivel el valor exacto del azimuth que corresponde a cada rayo. 

   ray_angle_res = np.unique( radar.ray_angle_res['data'] )
   if( np.size( ray_angle_res ) >= 2 )  :
      print('Warning: La resolucion en azimuth no es uniforme en los diferentes angulos de elevacion ')
      print('Warning: El codigo no esta preparado para considerar este caso y puede producir efectos indesaedos ')
   ray_angle_res=np.nanmean( ray_angle_res )
   #print ('The resolution in azimuth is: %5.3f' % ( ray_angle_res ) )


   levels=np.unique(radar.elevation['data'])
   azimuth=radar.azimuth['data']
   time=radar.time['data']

   order_azimuth=np.arange(0.0,360.0-ray_angle_res,ray_angle_res) #Asuming quasi regular azimuth location
   naz=np.size(order_azimuth)
   nel=np.size(levels)

   if ( var_name == 'altitude' ) :
      var=radar.gate_altitude['data']
   elif( var_name == 'longitude' ) :
      var=radar.gate_longitude['data'] 
   elif( var_name == 'latitude'  ) :
      var=radar.gate_latitude['data'] 
   elif( var_name == 'x' )         :
      var=radar.gate_x['data']
   elif( var_name == 'y' )         : 
      var=radar.gate_y['data']
   else  :
      var=radar.fields[var_name]['data'].data


      var[ var == undef ] = np.nan

   nr=var.shape[1]

   #Allocate arrays
   order_var    = undef + np.zeros((naz,nr,nel))
   order_time   =np.zeros((naz,nel)) 
   order_index  =undef + np.zeros((naz,nel))   #This variable can be used to convert back to the azimuth - range array
   azimuth_exact=undef + np.zeros((naz,nel))

   order_var[:]     = undef 
   order_time[:]    = undef 
   azimuth_exact[:] = undef

   

   for ilev in range(0, nel) :

      levmask= radar.elevation['data'] == levels[ilev] 

      min_index = np.min( np.where( levmask ) )


      #Find the azimuths corresponding to the current elevation.
      azlev=azimuth[ levmask ]
      timelev=time[ levmask ]
      #Get variabile values corresponding to the current elevation
      varlev=var[ levmask , : ]

      #For the first azimuth which is a special case because it contains zero.
      az_index=np.logical_or( azlev <= ray_angle_res/2.0 , azlev >= 360 - ray_angle_res/2.0 )
     
      if ( np.sum(az_index) > 0 ) : 
         order_var[0,:,ilev] = np.nanmean( varlev[az_index,:] , 0 )
         order_time[0,ilev] = np.nanmean( timelev[ az_index ] )
         azimuth_exact[0,ilev] = np.nanmean( azlev[ az_index ] )
         order_index[0,ilev]   = np.where( az_index )[0][0] + min_index

      #Para los que vienen despues.
      for iaz in range(1,naz) :
         #Search for all the rays that are close to the current azimuth and level.
         az_index=np.logical_and( azlev <= order_azimuth[iaz] + ray_angle_res/2.0 , azlev >= order_azimuth[iaz] - ray_angle_res/2.0 )
         if( np.sum( az_index ) > 0 ) :
            order_var[iaz,:,ilev] = np.nanmean( varlev[az_index,:] , 0 )
            order_time[iaz,ilev] = np.nanmean( timelev[ az_index ] )
            azimuth_exact[iaz,ilev] = np.nanmean( azlev[ az_index ] )
            order_index[iaz,ilev] = np.where(az_index)[0][0] + min_index #If multiple levels corresponds to a single azimuth / elevation chose the first one.

   order_var[ np.isnan( order_var ) ]= undef
   order_index[ np.isnan( order_index ) ]=undef

   return order_var , order_azimuth , levels , order_time , order_index , azimuth_exact

#From ARE order to Pyart order 
def order_variable_inv (  radar , var , order_index , undef )  :
    import numpy as np
   
    #Esta funcion es la inversa de la funcion order variable. Es decir que toma un array ordenado como azimuth , range y elevation y lo vuelve
    #a ordenar como azimuth-elevation y range. Es decir el orden original que se encuentra en los archivos con formato cfradial y que heredan los objetos radar de pyart.

    #var es la variable ordenada como var(azimuth,range,elevation)
    #order_index (azimuth,elevation) contiene la posicion original de los haces que fueron asignados a cada azimuth y elevacion por la funcion order_variable.
    #nr numero de puntos en la direccion del rango.
    #nb numero de beams. 

    na=var.shape[0]
    nr=var.shape[1]
    ne=var.shape[2]

    nb=radar.azimuth['data'].shape[0]

    output_var=np.ones((nb,nr)) * undef 
    
    for ia in range(0,na)  :

       for ie in range(0,ne)  :

          if ( not order_index[ia,ie] == undef )  :
     
              output_var[int(order_index[ia,ie]),:]=var[ia,:,ie] 

    return output_var


def get_rma_strat ( filename , radar )  :

    #Get the corresponding radar strategy depending on the filename.
    if radar.altitude_agl == None :
       radar.altitude_agl = dict()
    if radar.metadata == None :
       radar.metadata = dict()
    if radar.instrument_parameters == None :
       radar.instrument_parameters = dict()
    if radar.range == None :
       radar.range = dict()
    if radar.ray_angle_res == None :
       radar.ray_angle_res = dict()

    radar.altitude_agl['data'] = 0.0   

    radar.metadata['instrument_name'] = 'RMA' + filename[ filename.find('RMA') + 3 ]

    if '0122_03' in filename  :  #122-3 STRATEGY

        radar.instrument_parameters['nyquist_velocity']=13.25
        radar.ray_angle_res['data']=1.0
        radar.range['meters_between_gates']=300.0


    return radar


def interference_filter ( ref , undef , min_ref , r , my_conf )  :

   import numpy as np
   import matplotlib.pyplot as plt

   #ref a 3D array with the reflectivity in the na,nr,ne grid.
   #undef: value corresponding to the missing value in the reflectivity field.
   #min_ref: value corresponding to the no-rain value in the reflectivity field.
   #r : radar range (vector)
   #my_conf : configuration for this filter.

   na=ref.shape[0]
   nr=ref.shape[1]
   ne=ref.shape[2]

   tmp_ref=np.copy( ref ) 

   offset=my_conf['offset']
   att=my_conf['att']
   npass_filter=my_conf['npass_filter']
   percent_valid_threshold=my_conf['percent_valid_threshold']
   corr_threshold=my_conf['corr_threshold']
   ref_threshold=my_conf['ref_threshold']
   percent_ref_threshold=my_conf['percent_ref_threshold']

   tmp_ref[:,0:offset,:]=undef

   #Main loops

   for k in range(ne)  :

      for i in range(na)  :

         for iter in range(2) :
           #Perform an iterative fitting.
           if iter == 0 :
               tmp_mask = np.logical_and( tmp_ref[i,:,k] != undef , tmp_ref[i,:,k] < 10.0 )

               tmp_count = np.sum( (tmp_mask).astype(int) )/nr

           raux=r[tmp_mask]
           dbm= tmp_ref[i,:,k][tmp_mask]  - 20.0*np.log10( raux ) - 2.0*att*raux

           if tmp_count > percent_valid_threshold   :
              #TODO reemplazar por una regresion robusta que sea menos sensible
              #a los outliers. SCIKIT-LEARN tiene implementaciones en python de diversos
              #metodos para regresiones robustas. Seria cuestion de probar cual es el mas adecuado
              #para este problema.
              p=np.polyfit(raux,dbm,1)
              #tmpcorr=np.corrcoef(raux,dbm)
              #corrcoef=tmpcorr[0,1]
              a=p[0]
              b=p[1]


           else:
              a=0.0
              b=1.0

           dbzrayo=a*r+b + 20.0*np.log10( r ) + 2.0*att*r

           if (i==257) and (k==0) :
              print( dbzrayo[tmp_mask] , tmp_ref[i,:,k][tmp_mask] )
              print(tmp_count)
           corrcoef=np.corrcoef( dbzrayo[tmp_mask],tmp_ref[i,:,k][tmp_mask] )[0,1]
           #if i == 203 :
           #    print(iter,corrcoef,np.sum(tmp_mask.astype(int))/nr)

           #print(a,b,corrcoef )
           if i == 257 :
             plt.plot( dbzrayo[tmp_mask] ) 
             plt.plot( tmp_ref[i,:,k][tmp_mask] )
             #print( tmp_ref[i,:,k][tmp_mask][0:10])
             #print( dbzrayo[tmp_mask][0:10])
             #plt.show()

           #Consider the grid points that are close to the fitted interference reflectivity.
           #tmp_mask = tmp_ref[i,:,k] != undef
           tmp_mask = np.logical_and( np.abs( dbzrayo - tmp_ref[i,:,k] ) < ref_threshold , tmp_ref[i,:,k] != undef )
         plt.show()
         if k == 0 :
            print(i,corrcoef,np.sum(tmp_mask.astype(int))/nr)
         if( ( corrcoef > corr_threshold ) & ( np.sum(tmp_mask.astype(int))/nr > percent_ref_threshold ) )  :
         #This means that this ray is likely to be contaminated by interference.

            ref[i,:,k][ tmp_mask ] = undef

            zrayo = np.power(10.0,dbzrayo/10.0) 
            z     = np.power(10.0,ref[i,:,k]/10.0)

            #If the reflectivity is far from the fitted interference, and is greather than the fitted
            #Interference, then correct the power substracting the interference power.         
            tmp_mask2 = np.logical_and( np.logical_not( tmp_mask ), z - zrayo > 0.0  ) 
            ref[i,:,k][ tmp_mask2 ] = ( 10.0*np.log10( z - zrayo ) )[ tmp_mask2 ]
 
            #If the reflectivity is far from the fitted interference, and is smaller than the fitted interference
            #then set that pixel as an undef pixel.
           
            ref[i,:,k][ np.logical_and( np.logical_not( tmp_mask ) , z - zrayo <= 0.0 ) ] = undef

   #Additional filter for the remaining echoes
   #consider cyclic boundary conditions in azimuth.
   for ifilter in range( npass_filter )  :

      for k in range(ne) :

         for i in range(na) :
            if ( i > 1 ) & ( i < na-2 ) :
 
               #If we have reflectivity in only one ray but not in the neighbors this suggest an interference pattern.
               tmp_mask = np.logical_and( ref[i-1,:,k] <= min_ref , ref[i+1,:,k] <= min_ref ) 
               ref[i,:,k][ np.logical_and( tmp_mask , ref[i,:,k] > min_ref ) ] = min_ref

               tmp_mask = np.logical_and( ref[i-2,:,k] <= min_ref , ref[i+2,:,k] <= min_ref )
               ref[i,:,k][ np.logical_and( tmp_mask , ref[i,:,k] > min_ref ) ] = min_ref

            elif  i==na-1   :
               tmp_mask = np.logical_and( ref[i-1,:,k] <= min_ref , ref[0,:,k] <= min_ref )
               ref[i,:,k][ np.logical_and( tmp_mask , ref[i,:,k] > min_ref ) ] = min_ref

               tmp_mask = np.logical_and( ref[i-2,:,k] <= min_ref , ref[1,:,k] <= min_ref )
               ref[i,:,k][ np.logical_and( tmp_mask , ref[i,:,k] > min_ref ) ] = min_ref
           
            elif  i==na-3   :
               tmp_mask = np.logical_and( ref[i-1,:,k] <= min_ref , ref[i,:,k] <= min_ref )
               ref[i,:,k][ np.logical_and( tmp_mask , ref[i,:,k] > min_ref ) ] = min_ref

               tmp_mask = np.logical_and( ref[i-2,:,k] <= min_ref , ref[0,:,k] <= min_ref )
               ref[i,:,k][ np.logical_and( tmp_mask , ref[i,:,k] > min_ref ) ] = min_ref

            elif  i==0      :
               tmp_mask = np.logical_and( ref[na-1,:,k] <= min_ref , ref[i+1,:,k] <= min_ref )
               ref[i,:,k][ np.logical_and( tmp_mask , ref[i,:,k] > min_ref ) ] = min_ref

               tmp_mask = np.logical_and( ref[na-2,:,k] <= min_ref , ref[i+2,:,k] <= min_ref )
               ref[i,:,k][ np.logical_and( tmp_mask , ref[i,:,k] > min_ref ) ] = min_ref

            elif  i==1      :
               tmp_mask = np.logical_and( ref[i-1,:,k] <= min_ref , ref[i+1,:,k] <= min_ref )
               ref[i,:,k][ np.logical_and( tmp_mask , ref[i,:,k] > min_ref ) ] = min_ref

               tmp_mask = np.logical_and( ref[na-1,:,k] <= min_ref , ref[i+2,:,k] <= min_ref )
               ref[i,:,k][ np.logical_and( tmp_mask , ref[i,:,k] > min_ref ) ] = min_ref


   return ref


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
    np.reshape( my_topo['mean'].astype('f4') , (nr*na) ).tofile(f)
    np.reshape( my_topo['max'].astype('f4') , (nr*na) ).tofile(f)
    np.reshape( my_topo['min'].astype('f4') , (nr*na) ).tofile(f)
    np.reshape( my_topo['number'].astype('f4') , (nr*na) ).tofile(f)
    np.reshape( (my_topo['range'].data).astype('f4') , (nr*na) ).tofile(f)
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

    #Download the data to a raster file.
    #west, south, east, north 

    my_topo['mean'] = np.zeros( np.shape( my_topo['range'] ) , order='F' , dtype=np.float32 )  
    my_topo['min']  = np.zeros( np.shape( my_topo['range'] ) , order='F' , dtype=np.float32 )
    my_topo['max']  = np.zeros( np.shape( my_topo['range'] ) , order='F' , dtype=np.float32 )
    my_topo['number']  = np.zeros( np.shape( my_topo['range'] ) , order='F' , dtype=np.int32 ) 

    #Data is downloaded and processed into 1 deg patches.
    for ilon in range( minlon , maxlon )  :
      for ilat in range( minlat , maxlat )  : 

         print('Downloading data for tile Lon=',ilon,' Lat=',ilat) 

         #If raster file is not present, then download it from the internet.
         raster_file = raster_path + '/' + product + str(ilon) + '_' + str(ilat) + '.tif'

         if ( not os.path.isfile( raster_file ) )  :
             elevation.clip(bounds=(ilon,ilat,ilon+1,ilat+1),output=raster_file,product=product)

         #Read data from raster file.
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




