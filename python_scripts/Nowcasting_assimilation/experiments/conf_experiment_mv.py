import numpy as np 

#=================================================================

#=================================================================

GeneralConf=dict()

GeneralConf['DataBase'] = '/home/aarruti'
GeneralConf['radar'] = 'QC_ANGUIL'
GeneralConf['vol'] = '240'
GeneralConf['path_ref'] = 'QC_DATA'

GeneralConf['date_start_mv'] ='20100111093000'
GeneralConf['date_end_mv']   ='20100111113000'

GeneralConf['format_date'] = '%Y%m%d_%H%M'
GeneralConf['filedir_mv'] = 'asim_prueba'
GeneralConf['dt_entreimag'] = 600
#=================================================================
# MODEL SECTION : 
#=================================================================
#General model section

ModelConf=dict()

#General model section
ModelConf['nx'] = 240
ModelConf['ny'] = 240
ModelConf['box_size'] = 40 
ModelConf['dx'] = 2000
ModelConf['dt'] = 600
ModelConf['desp_max'] = 12
ModelConf['sigma'] = 20*ModelConf['dx']
ModelConf['sigmatrend'] = 10*ModelConf['dx']
ModelConf['Nx'] = ModelConf['nx']*ModelConf['ny'] 
ModelConf['nvar'] = 2


#Forcing section

ModelConf['Coef']=np.array([8])                     #Coefficient of parametrized forcing (polynom coefficients starting from coef[0]*x^0 + coef[1]*x ... ) 

ModelConf['NCoef']=np.size(ModelConf['Coef'])

#Space dependent parameter

ModelConf['FSpaceDependent']=False                     #If the forcing parameters will depend on the location.
ModelConf['FSpaceAmplitude']=np.array([1])             #Amplitude of space variantions (for each coefficient)
ModelConf['FSpaceFreq']     =np.array([1])             #Use integers >= 1

#Parameter random walk          

ModelConf['EnablePRF']=False                            #Activate Parameter random walk
ModelConf['CSigma']=np.array([0,0,0])                   #Parameter random walk sigma
ModelConf['CPhi'  ]=1.0                                 #Parameter random walk phi

#State random forcing parameters

ModelConf['EnableSRF']=False                            #Activate State random forcing.
ModelConf['XSigma']=0.0                                 #Amplitude of the random walk
ModelConf['XPhi'  ]=1.0                                 #Time autocorrelation parameter

ModelConf['XLoc'  ]=np.arange(1,ModelConf['nx']+1)      #Location of model grid points (1-nx)

#Two scale model parameters

ModelConf['TwoScaleParameters']=np.array([10,10,0])     #Small scale and coupling parameters C , B and Hint
                                                        #Set Hint /= 0 to enable two scale model

ModelConf['nxss']= ModelConf['nx'] * 8                  #Number of small scale variables
ModelConf['dtss']= ModelConf['dt'] / 5                  #Time step increment for the small scale variables

#=================================================================
#  DATA ASSIMILATION SECTION :
#=================================================================

DAConf=dict()

DAConf['nens'] = 20
DAConf['sigma_start'] = 0.5
DAConf['desv_gauss_start'] = 15 
DAConf['dalengh'] = 69            
DAConf['sigma_da'] = 0.2
DAConf['sigma_gauss_da'] = 15 
DAConf['LocScales']=np.array([15000.0,-1.0])             #Localization scale is space and time (negative means no localization)
DAConf['InfCoefs']=np.array([0.0,1.0,0.0,0.0,0.0])   #Mult inf, RTPS, RTPP, EPES, Additive inflation


 
#DAConf['Twin'] = True                                #When True, model configuration will be replaced by the model configuration in the nature run.

#DAConf['Freq'] = 4                                   #Assimilation frequency (in number of time steps)
#DAConf['TSFreq'] = 4                                 #Intra window ensemble output frequency (for 4D Data assimilation)

#DAConf['InfCoefs']=np.array([1.02,0.0,0.0,0.0,0.0])   #Mult inf, RTPS, RTPP, EPES, Additive inflation

#DAConf['LocScales']=np.array([4.0,-1.0])             #Localization scale is space and time (negative means no localization)

#Initial state ensemble.
#DAConf['InitialXSigma']=1                            #Initial ensemble spread for state variables.

#DAConf['UpdateSmoothCoef']=0.0                       #Data assimilation update smooth (for parameter estimation only)

#Parameter estimation/perturbation 

#DAConf['InitialPSigma']=np.array([0,0,0])            #Initial ensemble spread for the parameters. (0 means no parameter estimation)

#DAConf['InfCoefsP']=np.array([1.0,1.0,0.0,0.0,0.0])  #Mult inf, RTPS, RTPP, EPES, Additive inflation

#DAConf['UpdateSmoothCoefP']=0.0                      #Data assimilation update smooth (for parameter estimation only)

#DAConf['EstimateParameters']=False                   #Wether parameters will be estimated or not.

#DAConf['ParameterLocalizationType']=1                #1-Global parameter (no loc), 2-Averaged local estimation , 3-Full local estimation

#DAConf['LocScalesP']=np.array([4.0,-1.0])            #To be used with ParameterLocalizationTypes 2 or 3.


