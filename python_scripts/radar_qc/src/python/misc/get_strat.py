
def get_strat ( filename , radar )  :

    import numpy as np
    import numpy.ma as ma
    import os

    local_fill_value = -9999.0
    levels=np.unique(radar.elevation['data'])

    #Add missing structures to the radar object.
    if radar.altitude_agl == None :
       radar.altitude_agl = dict()
    if radar.metadata == None :
       radar.metadata = dict()
    if radar.instrument_parameters == None :
       radar.instrument_parameters = dict()

    if (not  'nyquist_velocity' in radar.instrument_parameters ) or ( radar.instrument_parameters['nyquist_velocity'] == None )  :

       radar.instrument_parameters['nyquist_velocity']=dict()
       radar.instrument_parameters['nyquist_velocity']['long_name']='unambiguous_doppler_velocity'
       radar.instrument_parameters['nyquist_velocity']['units']='meters per second'
       radar.instrument_parameters['nyquist_velocity']['_FillValue']= local_fill_value
       radar.instrument_parameters['nyquist_velocity']['meta_group']='instrument_parameters'

    if (not 'radar_beam_width_v' in radar.instrument_parameters ) or ( radar.instrument_parameters['radar_beam_width_v'] == None ) :
       radar.instrument_parameters['radar_beam_width_v']=dict()
       radar.instrument_parameters['radar_beam_width_v']['long_name']='half_power_radar_beam_width_v_channel'
       radar.instrument_parameters['radar_beam_width_v']['units']='degrees'
       radar.instrument_parameters['radar_beam_width_v']['_FillValue']= local_fill_value
       radar.instrument_parameters['radar_beam_width_v']['meta_group']='instrument_parameters'

    if (not 'radar_beam_width_h' in radar.instrument_parameters ) or ( radar.instrument_parameters['radar_beam_width_h'] == None ) :
       radar.instrument_parameters['radar_beam_width_h']=dict()
       radar.instrument_parameters['radar_beam_width_h']['long_name']='half_power_radar_beam_width_h_channel'
       radar.instrument_parameters['radar_beam_width_h']['units']='degrees'
       radar.instrument_parameters['radar_beam_width_h']['_FillValue']= local_fill_value
       radar.instrument_parameters['radar_beam_width_h']['meta_group']='instrument_parameters'

    if radar.range == None :
       radar.range = dict()
    if radar.ray_angle_res == None :
       radar.ray_angle_res = dict()
       radar.ray_angle_res['long_name']='angular_resolution_between_rays'
       radar.ray_angle_res['units']='degrees'
       radar.ray_angle_res['_FillValue']= local_fill_value

    #Get the corresponding radar strategy depending on the filename.
    radar.altitude_agl['data'] = 0.0


    if 'RMA' in filename  :

       radar.metadata['instrument_name'] = ( os.path.basename( filename ) ).split('_')[0]

       if '0117_01' in filename  :  #9005-1 STRATEGY
          nyquist_velocity     = 6.63
       if '0117_02' in filename  :  #9005-2 STRATEGY
          nyquist_velocity     = 33.04
       if '0117_03' in filename  :  #9005-3 STRATEGY
          nyquist_velocity     = 3.98

       if '0117_01' in filename  :  #122-1 STRATEGY
          nyquist_velocity     = 6.63
       if '0117_02' in filename  :  #122-2 STRATEGY
          nyquist_velocity     = 13.25

       if '0121_01' in filename  :  #122-1 STRATEGY
          nyquist_velocity     = 6.63
       if '0121_02' in filename  :  #122-2 STRATEGY
          nyquist_velocity     = 13.25

       if '0122_01' in filename  :  #122-1 STRATEGY
          nyquist_velocity     = 8.28
       if '0122_02' in filename  :  #122-2 STRATEGY
          nyquist_velocity     = 39.79
       if '0122_03' in filename  :  #122-3 STRATEGY
          nyquist_velocity     = 13.35

       if '0123_01' in filename  :  #123-1 STRATEGY
          nyquist_velocity     = 8.28
       if '0123_02' in filename  :  #123-2 STRATEGY
          nyquist_velocity     = 39.79
       if '0123_03' in filename  :  #123-3 STRATEGY
          nyquist_velocity     = 13.25
       if '0123_04' in filename  :  #123-4 STRATEGY
          nyquist_velocity     = 8.28

       if '0200_01' in filename  :  #200-1 STRATEGY
          nyquist_velocity     = 4.42
       if '0200_02' in filename  :  #200-2 STRATEGY
          nyquist_velocity     = 13.25

       if '0201_01' in filename  :  #201-1 STRATEGY
          nyquist_velocity     = 4.42
       if '0201_02' in filename  :  #201-2 STRATEGY
          nyquist_velocity     = 13.25
       if '0201_03' in filename  :  #201-3 STRATEGY
          nyquist_velocity     = 8.28

       if '0202_01' in filename  :  #200-1 STRATEGY
          nyquist_velocity     = 4.42
       if '0202_02' in filename  :  #200-2 STRATEGY
          nyquist_velocity     = 13.25


    if ( 'PAR' in filename ) or ( 'ANG' in filename ) or ( 'PER' in filename )  :

       if '/120/' in filename   : 
          nyquist_velocity = 39.8  #120

       if '/240/' in filename   :
          nyquist_velocity = 6.63  #240

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
    radar.instrument_parameters['nyquist_velocity']['data'] = ma.array(np.ones( np.shape(radar.azimuth['data']) )*nyquist_velocity , mask = np.zeros( np.shape(radar.azimuth['data']) , dtype=bool ) , fill_value = local_fill_value )
    radar.instrument_parameters['radar_beam_width_h']['data'] = ma.array( radar_beam_width_h , mask =False , fill_value = local_fill_value )
    radar.instrument_parameters['radar_beam_width_v']['data'] = ma.array( radar_beam_width_v , mask =False , fill_value = local_fill_value )
    radar.ray_angle_res['data'] = ma.array( np.ones( np.shape( levels ) )*ray_angle_res , mask = np.zeros( np.shape( levels ) , dtype=bool ) , fill_value = local_fill_value )
    radar.range['meters_between_gates']= meters_between_gates
    radar.range['meters_to_center_of_first_gate']= radar.range['data'][0]  #meters_between_gates / 2.0

    return radar


