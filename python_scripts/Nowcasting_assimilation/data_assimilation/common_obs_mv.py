import numpy as np 

#=======================================================================
# Model_to_Obs
#=======================================================================

def model_to_obs(nens, Ux, Vx, UObs, VObs, ObsLoc, ObsError, code=-999):

#    print( np.shape( UObs) , np.shape(VObs) , np.shape( Ux) , np.shape(Vx) )
    window_mask = np.logical_and(UObs[:]!=code,VObs[:]!=code)
    Ux = Ux[window_mask,:]
    Vx = Vx[window_mask,:]
    UObs = UObs[window_mask]
    VObs = VObs[window_mask]
    ObsLoc = ObsLoc[window_mask,:]
    ObsErrorW = ObsError[window_mask]

    #Add a quality control based on the distance between the first guess and the observations.
    #If the difference between the observation and the first guess is larger than 3 times the observation
    #error, then the observation is rejected

    Uqc_mask = np.abs( np.mean(Ux,1) - UObs ) < 3.0*ObsErrorW  
    Vqc_mask = np.abs( np.mean(Vx,1) - VObs ) < 3.0*ObsErrorW

    qc_mask = np.logical_and( Uqc_mask , Vqc_mask )

    Ux = Ux[qc_mask,:]
    Vx = Vx[qc_mask,:]
    UObs = UObs[qc_mask]
    VObs = VObs[qc_mask]
    ObsLoc = ObsLoc[qc_mask,:]
    ObsErrorW = ObsErrorW[qc_mask]

    return  Ux, Vx, UObs, VObs, ObsLoc, ObsErrorW
