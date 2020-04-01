# Author: Rapid Refresh Argentina Team
# License: BSD 3 clause

# History:  Created 10-2017

#Add path to additional modules.

#======================================
# ENVIRONMENTAL VARIABLES 
#======================================

import os
import numpy as np

os.environ['OMP_STACKSIZE'] = "1G" # visible in this process + all children

options = {}

#======================================
# I/O
#======================================

options['filename'] = ''          #The name of the input netcdf file
options['output_to_file'] = False #Wether the output will be written to a netcdf file.
options['filename_out'] = ''      #The name of the netcdf file where output will be written.
options['file_out_format']='NETCDF4'
options['keep_original_fields']=True
#Topography database paths
options['toporawdatapath']=""
options['toporadardatapath']=""

#======================================
# PLOTING PARAMETERS
#======================================

options['plot']=dict()
options['plot']['Enable']=True
options['plot']['Path']='./'
options['plot']['FigNameSufix']='.png'
options['plot']['VrMin']=-30
options['plot']['VrMax']=30
options['plot']['DbzMin']=0
options['plot']['DbzMax']=70
options['plot']['Elevs']=[0]             #List with elevations that will be ploted.
options['plot']['Show']=True
options['plot']['CmapWind']='pyart_NWSVel'
options['plot']['CmapDbz']='pyart_NWSRef'

#======================================
# QC PARAMETERS
#======================================

#General

options['name_ref'] ='ZH'    #'dBZ'              #Reflectivity
options['name_cref']='CZH'   #Corrected reflectivity (qc output)
options['name_v']   ='VRAD'  #'V'                #Dopper velocity
options['name_cv']  ='CVRAD' #Corrected wind (qc ouput)
options['name_rho'] ='RHOHV' #'RhoHV'            #Rho HV

options['name_model_ref_max']='dBZ_model_max'   #Maximum reflectivity from the model ensemble.
options['name_model_ref_min']='dBZ_model_min'   #Minimum reflectivity from the model ensemble.
options['name_model_ref_std']='dBZ_model_std'   #Standard deviation of reflectivity from the model ensemble.
options['name_model_v_max']  ='V_model_max'     #Maximum doppler velocity from the model ensemble.
options['name_model_v_min']  ='V_model_min'     #Minimum doppler velocity from the model ensemble.
options['name_model_v_std']  ='V_model_std'     #Standard deviation of reflectivity from the model ensemble.
options['radar_altitude_agl'] = 20.0            #Antena altitude to be used in case radar_altitude_agl is not defined.

options['is_rma'] = True                        #Wether this is an RMA file

options['norainrefval']=-0.1
options['undef']=-9.99e9


#Fuzzy logic parameter

options['w_tr']=0.5                  #When total normalized weight is greather than this value
                                     #the pixel is flagged as missing data.
#Each filter has associated an ify and ifx array that descrives the importance function.
#Each filter has also a w which determines its relative importance among the other filters.


#The weight is computed for each filter based on the importance function evaluated at the filter value.
#The relative importance of each filter is given by the importance_w.


#Dealiasing fliter parameters =============================================================
#Dealiasing parameters (pyart)
filter_name='Dealiasing'

options[filter_name]=dict()
options[filter_name]['flag']=True
options[filter_name]['interval_split']=3
options[filter_name]['skip_between_ray']=30
options[filter_name]['skip_along_ray']=200
options[filter_name]['nx']=3
options[filter_name]['ny']=3
options[filter_name]['nz']=0
options[filter_name]['code']=43
options[filter_name]['sequential']=True                    #Wheter this filter will affect the following filters.
options[filter_name]['order'] = [10]
options[filter_name]['var_update_list']=['v']              #Which variables will be filtered.

#Rho filter parameters  ===================================================================== 

filter_name='RhoFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                          #Enable / disable filter
options[filter_name]['nx']=2                                #NX
options[filter_name]['ny']=2                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['save']=False                          #Save filter aux fields to output?
options[filter_name]['ify']=np.array([1,1,0,0])             #Importance function y
options[filter_name]['ifx']=np.array([-2,0.5,0.8,2])        #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=10
options[filter_name]['force']=True                          #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['order'] = [20]
options[filter_name]['var_update_list']=['ref']             #Which variables will be filtered.
options[filter_name]['sequential']=True                     #Wheter this filter will affect the following filters.
options[filter_name]['fill_value']='min_ref'                  #Possible values, undef, min_ref or fill value

#Model filter parameters  ===================================================================== 

filter_name='ModelFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                          #Enable / disable filter
options[filter_name]['Enable_for_dBZ']=True                 #Enable / disable for reflectivity (filter has to be enabled)
options[filter_name]['Enable_for_Vr']=True                  #Enable / disable for doppler velocity (filter has to be enabled)
options[filter_name]['nx']=2                                #NX (smooth parameter)
options[filter_name]['ny']=2                                #NY (smooth parameter)
options[filter_name]['nz']=0                                #NZ (smooth parameter)
options[filter_name]['ref_error']=5.0                       #Estimated obs error for reflectivity
options[filter_name]['v_error']=2.0                         #Estimated obs error for doppler velocity
options[filter_name]['save']=False                           #Save filter aux fields to output?
options[filter_name]['uncertainty_metric']=0                # 0 - use range (max/min) , 1 - use spread
options[filter_name]['ify']=np.array([1,1,0,0,1,1])         #Importance function y
options[filter_name]['ifx']=np.array([-100,-2,-1,1,2,100])  #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=25
options[filter_name]['force']=False                         #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['order'] = [3]
options[filter_name]['var_update_list']=['ref']             #Which variables will be filtered.
options[filter_name]['sequential']=True                     #Wheter this filter will affect the following filters.
options[filter_name]['fill_value']='undef'                  #Possible values, undef, min_ref or fill value


#Echo top filter parameters ===================================================================

filter_name='EchoTopFilter'
options[filter_name]=dict()
options[filter_name]['flag']=True                          #Enable / disable filter
options[filter_name]['nx']=2                                #NX
options[filter_name]['ny']=2                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['fast_computation']=True               #Enable fast version of echo top computation.
options[filter_name]['save']=False                          #Save filter aux fields to output?
options[filter_name]['ify']=np.array([1,1,0,0])             #Importance function y
options[filter_name]['ifx']=np.array([0,2500,3000,20000])   #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=11
options[filter_name]['heigthtr']=3000                       #Do not use this filter if volume height 
                                                            #is lower than this threshold
options[filter_name]['force']=True                          #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['order'] = [15]
options[filter_name]['var_update_list']=['ref']             #Which variables will be filtered.
options[filter_name]['sequential']=True                     #Wheter this filter will affect the following filters.
options[filter_name]['fill_value']='min_ref'                #Possible values, undef, min_ref or fill value

#Echo depth filter parameters ===================================================================

filter_name='EchoDepthFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                          #Enable / disable filter
options[filter_name]['nx']=2                                #NX
options[filter_name]['ny']=2                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['save']=False                          #Save filter aux fields to output?
options[filter_name]['ify']=np.array([1,1,0,0])             #Importance function y
options[filter_name]['ifx']=np.array([0,2500,3000,20000])   #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=12
options[filter_name]['heigthtr']=3000                       #Do not use this filter if volume height 
options[filter_name]['force']=True                          #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['order'] = [5]
options[filter_name]['var_update_list']=['ref']             #Which variables will be filtered.
options[filter_name]['sequential']=True                     #Wheter this filter will affect the following filters.
options[filter_name]['fill_value']='undef'                  #Possible values, undef, min_ref or fill value

#Reflectivity speckle filter  parameters ==========================================================

filter_name='RefSpeckleFilter'
options[filter_name]=dict()
options[filter_name]['flag']=True                          #Enable / disable filter
options[filter_name]['nx']=2                                #NX
options[filter_name]['ny']=2                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['save']=False                          #Save filter aux fields to output?
options[filter_name]['ify']=np.array([1,1,0,0])             #Importance function y
options[filter_name]['ifx']=np.array([0,0.2,0.4,1])         #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=13
options[filter_name]['reftr']=5.0                           #Reflectivity threshold
options[filter_name]['force']=True                          #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['order'] = [30]
options[filter_name]['var_update_list']=['ref']             #Which variables will be filtered.
options[filter_name]['sequential']=True                     #Wheter this filter will affect the following filters.
options[filter_name]['fill_value']='min_ref'                  #Possible values, undef, min_ref or fill value

#Doppler speckle filter  parameters ==============================================================

filter_name='DopplerSpeckleFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                          #Enable / disable filter
options[filter_name]['nx']=2                                #NX
options[filter_name]['ny']=2                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['save']=False                          #Save filter aux fields to output?
options[filter_name]['ify']=np.array([1,1,0,0])             #Importance function y
options[filter_name]['ifx']=np.array([0,0.2,0.4,1])         #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=14
options[filter_name]['dvtr']=0.0                            #Wind threshold.
options[filter_name]['force']=True                          #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['order'] = [7]
options[filter_name]['var_update_list']=['v']               #Which variables will be filtered.
options[filter_name]['sequential']=True                     #Wheter this filter will affect the following filters.
options[filter_name]['fill_value']='undef'                  #Possible values, undef, min_ref or fill value

#Doppler texture filter  parameters ==============================================================

filter_name='DopplerTextureFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                          #Enable / disable filter
options[filter_name]['nx']=1                                #NX
options[filter_name]['ny']=2                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['save']=False                           #Save filter aux fields to output?
options[filter_name]['ify']=np.array([0,0,1,1])             #Importance function y
options[filter_name]['ifx']=np.array([0,2.5,10,200])        #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=15
options[filter_name]['force']=True                         #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['order'] = [8]
options[filter_name]['var_update_list']=['v']               #Which variables will be filtered.
options[filter_name]['sequential']=True                     #Wheter this filter will affect the following filters.
options[filter_name]['fill_value']='undef'                  #Possible values, undef, min_ref or fill value

#Reflectivity texture filter  parameters ==============================================================

filter_name='ReflectivityTextureFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                          #Enable / disable filter
options[filter_name]['nx']=3                                #NX for texture smoothing
options[filter_name]['ny']=5                                #NY for texture smoothing
options[filter_name]['nz']=0                                #NZ for texture smoothing
options[filter_name]['use_smooth_ref']=True                 #Wether high reflectivity cores will be protected.
options[filter_name]['smooth_ref_tr']=20.0                  #Reflectivity above this threshold wont be affected by this filter.
options[filter_name]['save']=False                          #Save filter aux fields to output?
options[filter_name]['ify']=np.array([0,0,1,1])             #Importance function y
options[filter_name]['ifx']=np.array([0,15,20,200])         #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=16
options[filter_name]['force']=True                          #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['order'] = [9]
options[filter_name]['var_update_list']=['ref']             #Which variables will be filtered.
options[filter_name]['sequential']=True                     #Wheter this filter will affect the following filters.
options[filter_name]['fill_value']='undef'                  #Possible values, undef, min_ref or fill value

#Attenuation parameters           ==============================================================
filter_name='AttenuationFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                           #Enable / Disable filter
options[filter_name]['nx']=0                                #NX
options[filter_name]['ny']=0                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['save']=False                          #Save filter aux fields to output?
options[filter_name]['ify']=np.array([0,0,1,1])             #Importance function y
options[filter_name]['ifx']=np.array([0,1,10,100])          #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=17
options[filter_name]['attcalerror']=1.0                     #Calibration error
options[filter_name]['is_power']=False                      #If input is in mm^6/m^3 set this to true.
options[filter_name]['att_coefs']=np.array([543,1.36,1.55e-3,1.30]) #Coefficients for the computation of attenuation (see below) 
options[filter_name]['force']=True                          #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['order'] = [2]
options[filter_name]['var_update_list']=['ref']             #Which variables will be filtered.
options[filter_name]['sequential']=True                     #Wheter this filter will affect the following filters.
options[filter_name]['fill_value']='undef'                  #Possible values, undef, min_ref or fill value

#Coefficients for C-Band radars based on
#Quantification of Path-Integrated Attenuation for X- and C-Band Weather
#Radar Systems Operating in Mediterranean Heavy Rainfall
#Delrieu, Andreiu, Creutin 1999 , Journal of Applied Meteorology
#a_coef=543
#b_coef=1.36
#c_coef=1.55e-3
#d_coef=1.30

#Blocking parameters               ==============================================================
#This filter is not part of the Fuzzy logic algorithm.
filter_name='BlockingFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                          #Enable / Disable filter
options[filter_name]['blocking_correction']=True            #Wether correction will be applied for partially blocked beams.
options[filter_name]['save']=False                          #Save blocking factor into qc_output dictionary.
options[filter_name]['code']=40                             #QC output code
options[filter_name]['ify']=np.array([0,0,1,1])             #Importance function y
options[filter_name]['ifx']=np.array([0,0.49,0.5,100])      #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=17                             #Filter QC code
options[filter_name]['force']=True                          #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['order'] = [40]
options[filter_name]['var_update_list']=['ref','v']         #Which variables will be filtered.
options[filter_name]['sequential']=True                     #Wheter this filter will affect the following filters.
options[filter_name]['fill_value']='undef'                  #Possible values, undef, min_ref or fill value

#Low elevation angles filter parameters ==============================================================
filter_name='LowElevFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                          #Enable / Disable filter
options[filter_name]['nx']=0                                #NX
options[filter_name]['ny']=0                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['save']=False                          #Save filter aux fields to output?
options[filter_name]['min_angle']=2.0                       #Minimun angle
options[filter_name]['height_thr']=20000.0                  #Heights above this wont be affected.
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=1
options[filter_name]['ify']=np.array([0,1])                 #Importance function y
options[filter_name]['ifx']=np.array([0,1])                 #Importance function x
options[filter_name]['force']=True                          #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['order'] = [12]
options[filter_name]['var_update_list']=['ref']             #Which variables will be filtered.
options[filter_name]['sequential']=True                     #Wheter this filter will affect the following filters.
options[filter_name]['fill_value']='undef'                  #Possible values, undef, min_ref or fill value

#Low doppler velocity filter            ==============================================================
filter_name='LowDopplerFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                          #Enable / Disable filter
options[filter_name]['nx']=0                                #NX
options[filter_name]['ny']=0                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['save']=False                          #Save filter aux fields to output?
options[filter_name]['ify']=np.array([1,1,0,0])             #Importance function y
options[filter_name]['ifx']=np.array([0,0.5,0.51,200])      #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']= 19
options[filter_name]['force']=True                          #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['height_thr']=1000                     #Height threshold.
options[filter_name]['order'] = [2]
options[filter_name]['var_update_list']=['v','ref']         #Which variables will be filtered.
options[filter_name]['sequential']=True                     #Wheter this filter will affect the following filters.
options[filter_name]['fill_value']='undef'                  #Possible values, undef, min_ref or fill value

#Interference filter            ==============================================================
#This filter is not included in the Fuzzy-logic approach.
filter_name='InterferenceFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                          #Enable / Disable filter
options[filter_name]['save']=False                          #Save filter aux fields to output?
options[filter_name]['nx']=0                                #NX
options[filter_name]['ny']=7                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['code']= 20
options[filter_name]['ify']=np.array([0,1.0])                 #Importance function y
options[filter_name]['ifx']=np.array([0,1.0])                 #Importance function x
options[filter_name]['force']=True                          #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['att']=0.01e-3                         #Atmospheric gases attenuation constant (dB / km )
options[filter_name]['Smooth_Ref']=True                     #Smooth reflectivity before applying robust regression
options[filter_name]['Power_Regression']=True               #Wether robust regression will be performed in Dbz scale or linear scale
options[filter_name]['offset']=100                          #Number of pixels from radar location that will be ignored in the regression.
options[filter_name]['AzimuthFilter']=True                  #Enable filter to remove isolated pixels in azimuth.
options[filter_name]['ElevationFilter']=False               #Enable filter to remove isolated pixels in elevation.
options[filter_name]['npass_filter']=3                      #Number of passes of the azimuthal continuity filter.
options[filter_name]['percent_valid_threshold']=0.2         #Rays with valid pixels over this percentaje will be examinated.
options[filter_name]['corr_threshold']=np.array([0.97,0.9,0.8,0.7])                 #Rays that correlates well with the interference pattern will be flagged as 
options[filter_name]['azimuth_ref_diff_threshold']=np.array([0.1,0.15,0.2,0.3])           #If more than this percent of the ray correlates well with the interference pattern, then
                                                            #the ray is flagged as contaminated by interference.#the ray is flagged as contaminated by interference.
options[filter_name]['ref_threshold']=7.0                   #Reflectivity threshold to count pixels which are close to the interference pattern.                 
options[filter_name]['percent_ref_threshold']=0.3           #Threshold to measure the azimuthal coherence of the reflectivity.
options[filter_name]['order'] = [14]
options[filter_name]['var_update_list']=['ref']             #Which variables will be filtered.
options[filter_name]['sequential']=True                     #Wheter this filter will affect the following filters.
options[filter_name]['fill_value']='undef'                  #Possible values, undef, min_ref or fill value

#Dealiasing border filter            ==============================================================
filter_name='DealiasingEdgeFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                          #Enable / Disable filter
options[filter_name]['nx']=3                                #NX
options[filter_name]['ny']=3                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['save']=False                          #Save filter aux fields to output?
options[filter_name]['ify']=np.array([0,1.0])               #Importance function y
options[filter_name]['ifx']=np.array([0,1.0])               #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=22                              
options[filter_name]['force']=True                          #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['order'] = [15]
options[filter_name]['var_update_list']=['v']               #Which variables will be filtered.
options[filter_name]['sequential']=True                     #Wheter this filter will affect the following filters.
options[filter_name]['fill_value']='undef'                  #Possible values, undef, min_ref or fill value

#Doppler Local Std Filter            ==============================================================
filter_name='DopplerLocalStdFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                          #Enable / Disable filter
options[filter_name]['nx']=1                                #NX
options[filter_name]['ny']=1                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['ify']=np.array([0,0,1,1])             #Importance function y
options[filter_name]['ifx']=np.array([0,3,4,10])            #Importance function x
options[filter_name]['save']=False                          #Save filter aux fields to output?
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=23
options[filter_name]['force']=True                         #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['order'] = [16]
options[filter_name]['var_update_list']=['v']               #Which variables will be filtered.
options[filter_name]['sequential']=True                     #Wheter this filter will affect the following filters.
options[filter_name]['fill_value']='undef'                  #Possible values, undef, min_ref or fill value

#Doppler Spatial Coherence Filter ==============================================================
filter_name='DopplerSpatialCoherenceFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                           #Enable / Disable filter
options[filter_name]['nx']=2                                #NX
options[filter_name]['ny']=2                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['threshold_undef']=0.1                 #Minimum percentage of valid points required
options[filter_name]['threshold_corr'] =0.6                 #Minimum correlation required to keep a ray.
options[filter_name]['threshold_coherence_index']=1         #Threshold to decied which pixels will be removed.
options[filter_name]['compute_horizontal_coherence']=True   #Flag to consider coherence in azimuth.
options[filter_name]['compute_vertical_coherence']=False    #Flag to consider coherence in elevation
options[filter_name]['npass_filter']=2                      #Number of applications of the issolated data filter.
options[filter_name]['azimuthfilter']=True                  #Apply issolated data filter in azimuth.
options[filter_name]['rangefilter']=True                    #Apply issolated data filter in range.
options[filter_name]['elevationfilter']=False               #Apply issolated data filter in elevation
options[filter_name]['enable_speckle']=True                 #Wether speckle filter will be applied to the remaining data.
options[filter_name]['speckle_threshold']=0.3               #Threshold to discard pixels based on speckle index.
options[filter_name]['consistency_metric']='constant'       #Possible values are 'Ransac' or 'Constant'
options[filter_name]['constant_consistency_threshold']=3.5  #Consecutive values which are farther than this threshold will be flagged.
options[filter_name]['ify']=np.array([0,0.9,20])            #Importance function y
options[filter_name]['ifx']=np.array([0,1,1])               #Importance function x
options[filter_name]['save']=False                          #Save filter aux fields to output?
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=34
options[filter_name]['force']=True                          #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['order'] = [10]
options[filter_name]['sequential']=True                        
options[filter_name]['var_update_list']=['v']               #Which variables will be filtered.
options[filter_name]['sequential']=True                     #Wheter this filter will affect the following filters.
options[filter_name]['fill_value']='undef'                  #Possible values, undef, min_ref or fill value

#Doppler Noise filter FIRST PASS, BEFORE DEALIASING    ==============================================================

filter_name='DopplerNoiseFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                          #Enable / Disable filter
options[filter_name]['nx']=[1,10]                           #NX
options[filter_name]['ny']=[1,10]                           #NY
options[filter_name]['nz']=[0,0]                            #NZ
options[filter_name]['threshold']=[2.0,15.0]
options[filter_name]['n_filter_pass']=[3,3]                 #Filter repetition
options[filter_name]['save']=False                          #Save filter aux fields to output?
options[filter_name]['ify']=np.array([0,1])                 #Importance function y
options[filter_name]['ifx']=np.array([0,1])                 #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']= 1
options[filter_name]['force']=True                          #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['order'] = [0]
options[filter_name]['var_update_list']=['v']               #Which variables will be filtered.
options[filter_name]['sequential']=True                     #Wheter this filter will affect the following filters.
options[filter_name]['fill_value']='undef'                  #Possible values, undef, min_ref or fill value


#Missing reflectivity filter ==================================================================
#Detects holes in high reflectivity regions. 
filter_name='MissingRefFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                           #Enable / Disable filter
options[filter_name]['threshold']=10                        #Threshold to detect sudden jumps in reflectivity between two consecutive pixels.
options[filter_name]['nmissing_max']=15                     #Maximum number of missing values in a radial beam.
options[filter_name]['save']=False                          #Save filter aux fields to output?
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['ify']=np.array([0,1.0])               #Importance function y
options[filter_name]['ifx']=np.array([0,1.0])               #Importance function x
options[filter_name]['code']= 32
options[filter_name]['force']=True                          #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['order'] = [60]
options[filter_name]['var_update_list']=['ref']             #Which variables will be filtered.
options[filter_name]['sequential']=True                     #Wheter this filter will affect the following  filters
options[filter_name]['fill_value']='undef'                  #Possible values, undef, min_ref or fill value

#Doppler Range Filter  ==================================================================
#Filter gates with out of range value for the Doppler field.
filter_name='DopplerRangeFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                           #Enable / Disable filter
options[filter_name]['min']=-100                            #Threshold to detect sudden jumps in reflectivity between two consecutive pixels.
options[filter_name]['max']=100                             #Maximum number of missing values in a radial beam.
options[filter_name]['save']=False                          #Save filter aux fields to output?
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['ify']=np.array([0,1.0])               #Importance function y
options[filter_name]['ifx']=np.array([0,1.0])               #Importance function x
options[filter_name]['code']= 33
options[filter_name]['force']=True                          #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['order'] = [61]
options[filter_name]['var_update_list']=['v']               #Which variables will be filtered.
options[filter_name]['sequential']=True                     #Wheter this filter will affect the following  filters
options[filter_name]['fill_value']='undef'                  #Possible values, undef, min_ref or fill value

#Doppler Range Filter  ==================================================================
#Filter gates with out of range value for the Doppler field.
filter_name='RefRangeFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                           #Enable / Disable filter
options[filter_name]['min']=-10.0                           #Threshold to detect sudden jumps in reflectivity between two consecutive pixels.
options[filter_name]['max']=80.0                            #Maximum number of missing values in a radial beam.
options[filter_name]['save']=False                          #Save filter aux fields to output?
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['ify']=np.array([0,1.0])               #Importance function y
options[filter_name]['ifx']=np.array([0,1.0])               #Importance function x
options[filter_name]['code']= 32
options[filter_name]['force']=True                          #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['order'] = [60]
options[filter_name]['var_update_list']=['ref']             #Which variables will be filtered.
options[filter_name]['sequential']=True                     #Wheter this filter will affect the following  filters
options[filter_name]['fill_value']='undef'                  #Possible values, undef, min_ref or fill value



