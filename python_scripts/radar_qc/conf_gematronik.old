#print __doc__

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

#======================================
# QC PARAMETERS
#======================================
options = {}

#General

options['name_ref'] ='dBZ'              #Reflectivity
options['name_cref']='CdBZ'             #Corrected reflectivity (qc output)
options['name_v']   ='V'                #Dopper velocity
options['name_cv']  ='CV'               #Corrected wind (qc ouput)
options['name_rho'] ='RhoHV'            #Rho HV

options['name_model_ref_max']='dBZ_model_max'   #Maximum reflectivity from the model ensemble.
options['name_model_ref_min']='dBZ_model_min'   #Minimum reflectivity from the model ensemble.
options['name_model_ref_std']='dBZ_model_std'   #Standard deviation of reflectivity from the model ensemble.
options['name_model_v_max']  ='V_model_max'     #Maximum doppler velocity from the model ensemble.
options['name_model_v_min']  ='V_model_min'     #Minimum doppler velocity from the model ensemble.
options['name_model_v_std']  ='V_model_std'     #Standard deviation of reflectivity from the model ensemble.
options['radar_altitude_agl'] = 20.0            #Antena altitude to be used in case radar_altitude_agl is not defined.

options['is_rma'] = False                        #Wether this is an RMA file

options['norainrefval']=-0.1
options['undef']=-9.99e9

#I/O options

options['qced_to_file']=False        #Write QCed data to the input file.


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
options[filter_name]['skip_between_ray']=10
options[filter_name]['skip_along_ray']=10
options[filter_name]['nx']=3
options[filter_name]['ny']=3
options[filter_name]['nz']=0
options[filter_name]['code']=43
options[filter_name]['sequential'] = True      #If sequential = true then following filters will use the dealiased data.
options[filter_name]['order'] = [1]

#Rho filter parameters  ===================================================================== 

filter_name='RhoFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                          #Enable / disable filter
options[filter_name]['nx']=2                                #NX
options[filter_name]['ny']=2                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['save']=True                           #Save filter aux fields to output?
options[filter_name]['ify']=np.array([1,1,1,0,0])           #Importance function y
options[filter_name]['ifx']=np.array([-2,0.5,0.8,0.85,2])   #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=10
options[filter_name]['force']=False                         #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['fill_value']='minref'                 #Possible values, undef, minref or fill value
options[filter_name]['order'] = [2]

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
options[filter_name]['save']=True                           #Save filter aux fields to output?
options[filter_name]['uncertainty_metric']=0                # 0 - use range (max/min) , 1 - use spread
options[filter_name]['ify']=np.array([1,1,0,0,1,1])         #Importance function y
options[filter_name]['ifx']=np.array([-100,-2,-1,1,2,100])  #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=25
options[filter_name]['force']=False                         #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['fill_value']='minref'                 #Possible values, undef, minref or fill value
options[filter_name]['order'] = [2]

#Echo top filter parameters ===================================================================

filter_name='EchoTopFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                          #Enable / disable filter
options[filter_name]['nx']=2                                #NX
options[filter_name]['ny']=2                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['save']=True                           #Save filter aux fields to output?
options[filter_name]['ify']=np.array([1,1,0,0])             #Importance function y
options[filter_name]['ifx']=np.array([0,2500,3000,20000])   #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=11
options[filter_name]['heigthtr']=3000                       #Do not use this filter if volume height 
                                                            #is lower than this threshold
options[filter_name]['force']=False                         #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['fill_value']='minref'                 #Possible values, undef, minref or fill value
options[filter_name]['order'] = [3]

#Echo depth filter parameters ===================================================================

filter_name='EchoDepthFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                          #Enable / disable filter
options[filter_name]['nx']=2                                #NX
options[filter_name]['ny']=2                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['save']=True                           #Save filter aux fields to output?
options[filter_name]['ify']=np.array([1,1,0,0])             #Importance function y
options[filter_name]['ifx']=np.array([0,2500,3000,20000])   #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=12
options[filter_name]['heigthtr']=3000                       #Do not use this filter if volume height 
options[filter_name]['force']=False                         #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['fill_value']='minref'                 #Possible values, undef, minref or fill value
options[filter_name]['order'] = [4]

#Reflectivity speckle filter  parameters ==========================================================

filter_name='RefSpeckleFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                          #Enable / disable filter
options[filter_name]['nx']=2                                #NX
options[filter_name]['ny']=2                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['save']=True                           #Save filter aux fields to output?
options[filter_name]['ify']=np.array([1,1,0,0])             #Importance function y
options[filter_name]['ifx']=np.array([0,0.2,0.4,1])         #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=13
options[filter_name]['reftr']=5                             #Reflectivity threshold
options[filter_name]['force']=False                         #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['fill_value']='minref'                 #Possible values, undef, minref or fill value
options[filter_name]['order'] = [5]

#Doppler speckle filter  parameters ==============================================================

filter_name='DopplerSpeckleFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                          #Enable / disable filter
options[filter_name]['nx']=2                                #NX
options[filter_name]['ny']=2                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['save']=True                           #Save filter aux fields to output?
options[filter_name]['ify']=np.array([1,1,0,0])             #Importance function y
options[filter_name]['ifx']=np.array([0,0.2,0.4,1])         #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=14
options[filter_name]['dvtr']=0.0                            #Wind threshold.
options[filter_name]['force']=False                         #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['fill_value']='undef'                  #Possible values, undef, minref or fill value
options[filter_name]['order'] = [6]

#Doppler texture filter  parameters ==============================================================

filter_name='DopplerTextureFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                          #Enable / disable filter
options[filter_name]['nx']=1                                #NX
options[filter_name]['ny']=2                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['save']=True                           #Save filter aux fields to output?
options[filter_name]['ify']=np.array([0,0,1,1])             #Importance function y
options[filter_name]['ifx']=np.array([0,2.5,10,200])        #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=15
options[filter_name]['force']=False                         #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['fill_value']='undef'                  #Possible values, undef, minref or fill value
options[filter_name]['order'] = [7]

#Reflectivity texture filter  parameters ==============================================================

filter_name='ReflectivityTextureFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                          #Enable / disable filter
options[filter_name]['nx']=3                                #NX for texture smoothing
options[filter_name]['ny']=5                                #NY for texture smoothing
options[filter_name]['nz']=0                                #NZ for texture smoothing
options[filter_name]['use_smooth_ref']=True                 #Wether high reflectivity cores will be protected.
options[filter_name]['smooth_ref_tr']=15.0                  #Reflectivity above this threshold wont be affected by this filter.
options[filter_name]['save']=True                           #Save filter aux fields to output?
options[filter_name]['ify']=np.array([0,0,1,1])             #Importance function y
options[filter_name]['ifx']=np.array([0,15,20,200])         #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=16
options[filter_name]['force']=False                         #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['fill_value']='undef'                  #Possible values, undef, minref or fill value
options[filter_name]['order'] = [8]

#Attenuation parameters           ==============================================================
filter_name='AttenuationFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                          #Enable / Disable filter
options[filter_name]['nx']=0                                #NX
options[filter_name]['ny']=0                                #NY
options[filter_name]['nz']=0                                #NZ
options[filter_name]['save']=True                           #Save filter aux fields to output?
options[filter_name]['ify']=np.array([0,0,1,1])             #Importance function y
options[filter_name]['ifx']=np.array([0,1,10,100])          #Importance function x
options[filter_name]['w']=1.0                               #Relative parameter weigth. 
options[filter_name]['code']=17
options[filter_name]['attcalerror']=1.0                     #Calibration error
options[filter_name]['is_power']=False                      #If input is in mm^6/m^3 set this to true.
options[filter_name]['att_coefs']=np.array([543,1.36,1.55e-3,1.30]) #Coefficients for the computation of attenuation (see below) 
options[filter_name]['force']=False                         #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                     #Threshold for force
options[filter_name]['fill_value']='undef'                  #Possible values, undef, minref or fill value
options[filter_name]['order'] = [9]

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
options[filter_name]['blocking_threshold']=0.5              #Beams with blocking above this threshold will be eliminated.
options[filter_name]['save']=True                           #Save blocking factor into qc_output dictionary.
options[filter_name]['code']=40                             #QC output code
options[filter_name]['fill_value']='undef'                  #Possible values, undef, minref or fill value
options[filter_name]['order'] = [10]

#Low elevation angles filter parameters ==============================================================
filter_name='LowElevFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                             #Enable / Disable filter
options[filter_name]['nx']=0                                   #NX
options[filter_name]['ny']=0                                   #NY
options[filter_name]['nz']=0                                   #NZ
options[filter_name]['save']=True                              #Save filter aux fields to output?
options[filter_name]['min_angle']=2.0                          #Minimun angle
options[filter_name]['height_thr']=20000.0                     #Heights above this wont be affected.
options[filter_name]['w']=1.0                                  #Relative parameter weigth. 
options[filter_name]['code']=18
options[filter_name]['force']=False                            #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                        #Threshold for force
options[filter_name]['fill_value']='undef'                     #Possible values, undef, minref or fill value
options[filter_name]['order'] = [11]

#Low doppler velocity filter            ==============================================================
filter_name='LowDopplerFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                             #Enable / Disable filter
options[filter_name]['nx']=0                                   #NX
options[filter_name]['ny']=0                                   #NY
options[filter_name]['nz']=0                                   #NZ
options[filter_name]['save']=True                              #Save filter aux fields to output?
options[filter_name]['ify']=np.array([0,0,1,1,0,0])            #Importance function y
options[filter_name]['ifx']=np.array([-200,-0.5,-0.25,0.25,0.5,200]) #Importance function x
options[filter_name]['w']=1.0                                  #Relative parameter weigth. 
options[filter_name]['code']= 19
options[filter_name]['force']=False                            #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                        #Threshold for force
options[filter_name]['height_thr']=1000                        #Height threshold.
options[filter_name]['use_terrain']=True                       #Wether AGL height is used.
options[filter_name]['fill_value']='undef'                     #Possible values, undef, minref or fill value
options[filter_name]['order'] = [12]

#Interference filter            ==============================================================
#This filter is not included in the Fuzzy-logic approach.
filter_name='InterferenceFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                             #Enable / Disable filter
options[filter_name]['save']=True                              #Save filter aux fields to output?
options[filter_name]['nx']=0                                   #NX
options[filter_name]['ny']=4                                   #NY
options[filter_name]['nz']=0                                   #NZ
options[filter_name]['code']= 20
options[filter_name]['att']=0.01e-3                            #Atmospheric gases attenuation constant (dB / km )
options[filter_name]['Smooth_Ref']=True                        #Smooth reflectivity before applying robust regression
options[filter_name]['Power_Regression']=True                  #Wether robust regression will be performed in Dbz scale or linear scale
options[filter_name]['offset']=100                             #Number of pixels from radar location that will be ignored in the regression.
options[filter_name]['AzimuthFilter']=True                     #Enable filter to remove isolated pixels in azimuth.
options[filter_name]['ElevationFilter']=True                   #Enable filter to remove isolated pixels in elevation.
options[filter_name]['npass_filter']=2                         #Number of passes of the azimuthal continuity filter.
options[filter_name]['percent_valid_threshold']=0.1            #Rays with valid pixels over this percentaje will be examinated.
options[filter_name]['corr_threshold']=0.5                     #Rays that correlates well with the interference pattern will be flagged as 
                                                               #contaminated.
options[filter_name]['ref_threshold']=5.0                      #Reflectivity threshold to count pixels which are close to the interference pattern.
options[filter_name]['percent_ref_threshold']=0.3              #If more than this percent of the ray correlates well with the interference pattern, then
                                                               #the ray is flagged as contaminated by interference.
options[filter_name]['fill_value']='minref'                    #Possible values, undef, minref or fill value
options[filter_name]['order'] = [13]

#Dealiasing border filter            ==============================================================
filter_name='DealiasingBorderFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                             #Enable / Disable filter
options[filter_name]['nx']=3                                   #NX
options[filter_name]['ny']=3                                   #NY
options[filter_name]['nz']=3                                   #NZ
options[filter_name]['save']=True                              #Save filter aux fields to output?
options[filter_name]['w']=1.0                                  #Relative parameter weigth. 
options[filter_name]['code']= 1
options[filter_name]['force']=False                            #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                        #Threshold for force
options[filter_name]['fill_value']='undef'                     #Possible values, undef, minref or fill value
options[filter_name]['order'] = [14]

#Doppler Local Std Filter            ==============================================================
filter_name='DopplerLocalStdFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                             #Enable / Disable filter
options[filter_name]['nx']=1                                   #NX
options[filter_name]['ny']=1                                   #NY
options[filter_name]['nz']=0                                   #NZ
options[filter_name]['ify']=np.array([0,0,1,1])                #Importance function y
options[filter_name]['ifx']=np.array([0,5,6,10])               #Importance function x
options[filter_name]['save']=True                              #Save filter aux fields to output?
options[filter_name]['w']=1.0                                  #Relative parameter weigth. 
options[filter_name]['code']= 1
options[filter_name]['force']=False                            #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                        #Threshold for force
options[filter_name]['fill_value']='undef'                     #Possible values, undef, minref or fill value
options[filter_name]['order'] = [15]

#Doppler Spatial Coherence Filter ==============================================================
filter_name='DopplerSpatialCoherenceFilter'
options[filter_name]=dict()
options[filter_name]['flag']=True                              #Enable / Disable filter
options[filter_name]['nx']=2                                   #NX
options[filter_name]['ny']=2                                   #NY
options[filter_name]['nz']=0                                   #NZ
options[filter_name]['threshold_undef']=0.1                    #Minimum percentage of valid points required
options[filter_name]['threshold_corr'] =0.7                    #Minimum correlation required to keep a ray.
options[filter_name]['threshold_coherence_index']=1            #Threshold to decied which pixels will be removed.
options[filter_name]['compute_horizontal_coherence']=True      #Flag to consider coherence in azimuth.
options[filter_name]['compute_vertical_coherence']=False       #Flag to consider coherence in elevation
options[filter_name]['npass_filter']=2                         #Number of applications of the issolated data filter.
options[filter_name]['azimuthfilter']=True                     #Apply issolated data filter in azimuth.
options[filter_name]['rangefilter']=True                       #Apply issolated data filter in range.
options[filter_name]['elevationfilter']=False                  #Apply issolated data filter in elevation
options[filter_name]['enable_speckle']=True                    #Wether speckle filter will be applied to the remaining data.
options[filter_name]['speckle_threshold']=0.3                  #Threshold to discard pixels based on speckle index.
options[filter_name]['consistency_metric']='constant'          #Possible values are 'Ransac' or 'Constant'
options[filter_name]['constant_consistency_threshold']=3.5     #Consecutive values which are farther than this threshold will be flagged.
options[filter_name]['ify']=np.array([0,0,1,1])                #Importance function y
options[filter_name]['ifx']=np.array([0,5,6,10])               #Importance function x
options[filter_name]['save']=True                              #Save filter aux fields to output?
options[filter_name]['w']=1.0                                  #Relative parameter weigth. 
options[filter_name]['code']= 1
options[filter_name]['force']=False                            #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                        #Threshold for force
options[filter_name]['fill_value']='undef'                     #Possible values, undef, minref or fill value
options[filter_name]['order'] = [16]

#Doppler Noise filter            ==============================================================
filter_name='DopplerNoiseFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                             #Enable / Disable filter
options[filter_name]['nx']=1                                   #NX
options[filter_name]['ny']=1                                   #NY
options[filter_name]['nz']=0                                   #NZ
options[filter_name]['nx2']=10                                 #NX
options[filter_name]['ny2']=10                                 #NY
options[filter_name]['nz2']=0                                  #NZ
options[filter_name]['threshold_1']=1.0
options[filter_name]['threshold_2']=15       
options[filter_name]['n_filter_pass']=3                        #Filter repetition
options[filter_name]['save']=True                              #Save filter aux fields to output?
options[filter_name]['ify_1']=np.array([0,0,1,1])              #Importance function y
options[filter_name]['ifx_1']=np.array([0,2,6,10])             #Importance function x
options[filter_name]['ify_2']=np.array([0,0,1,1])              #Importance function y
options[filter_name]['ifx_2']=np.array([0,10,15,20])           #Importance function x
options[filter_name]['w']=1.0                                  #Relative parameter weigth. 
options[filter_name]['code']= 1
options[filter_name]['force']=False                            #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                        #Threshold for force
options[filter_name]['fill_value']='undef'                     #Possible values, undef, minref or fill value
options[filter_name]['order'] = [17]

#Missing reflectivity filter ==================================================================
#Detects holes in high reflectivity regions. 
filter_name='MissingRefFilter'
options[filter_name]=dict()
options[filter_name]['flag']=False                             #Enable / Disable filter
options[filter_name]['threshold']=10                           #Threshold to detect sudden jumps in reflectivity between two consecutive pixels.
options[filter_name]['nmissing_max']=15                        #Maximum number of missing values in a radial beam.
options[filter_name]['save']=True                              #Save filter aux fields to output?
options[filter_name]['w']=1.0                                  #Relative parameter weigth. 
options[filter_name]['code']= 22
options[filter_name]['force']=False                            #Wether we will reject data based on this filter alone.
options[filter_name]['force_value']=0.5                        #Threshold for force
options[filter_name]['fill_value']='undef'                     #Possible values, undef, minref or fill value
options[filter_name]['order'] = [18]

#Topography parameters

#options['toporawdatapath']="/home/jruiz/share/radar_qc_da/data/terrain_data/raw/"
#options['toporadardatapath']="/home/jruiz/share/radar_qc_da/data/terrain_data/radar/"

options['toporawdatapath']="/media/jruiz/PAWR/Dropbox/DATA/radar_qc_da/data/terrain_data/raw/"
options['toporadardatapath']="/media/jruiz/PAWR/Dropbox/DATA/radar_qc_da/data/terrain_data/radar/"


