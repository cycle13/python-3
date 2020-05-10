from common_functions import common_functions as cf
import matplotlib.pyplot as plt
import datetime as dt
import numpy as np
import pickle as pkl
import os
import time
import warnings
warnings.filterwarnings("ignore")

exp_path='/home/juan/tmp_scale_to_radar/output_data/LE_D1_1km_30sec_OFP_V2/'
rad_path='/home/juan/tmp_scale_to_radar/radar_data/'

init_freq=600   #Forecast initialization frequency in seconds.
out_freq =30    #Forecast output frequency in seconds.
for_len  =1800  #Forecast length in seconds.
time_tol =15    #Maximum time tolerance (if scale time and time in nearest radar data are larger than this the interpolation is not performed.)

minref=0.0      #Ref values below this threshold will be assumed equal to the threshold.

itime = dt.datetime(2013,7,13,5,5,0)  #Initial time.
etime = dt.datetime(2013,7,13,5,5,0)   #End time.

radar_files=None

nfor = int( np.floor( ( etime - itime ).seconds / init_freq ) )     #Number of forecasts.
nled = int( np.floor( for_len / out_freq ) )                        #Number of lead times.

first_read = True
scores = dict()

#Time loop to read forecast and interpolate them to nearest radar data.
ctime = itime
ifor = 0
while ( ctime <= etime ) :

    fi_time = ctime
    fe_time = ctime + dt.timedelta( seconds = for_len )
    cftime = fi_time

    forecast_time = 0
    while ( cftime <= fe_time ) :
        cpath = exp_path + dt.datetime.strftime( fi_time ,'%Y%m%d%H%M%S' ) + '/fcstrad/mean/'

        grid_file = cpath + '/fcstrad_grid' + dt.datetime.strftime( cftime ,'%Y%m%d%H%M%S' ) + '.pkl' 
        if os.path.exists( grid_file  ) :
            filehandler = open ( grid_file , "rb" )
            radar_grid = pickle.load(filehandler)

            if first_read  :
                [nz,ny,nx,n,nv] = np.shape( radar_grid['data_ave'] )
                #Allocate memory for the scores.
                scores['ref_ets']  = np.zeros( ( nt , nfor , nled ) )
                scores['ref_biasf']= np.zeros( ( nt , nfor , nled ) )
                scores['ref_csi']  = np.zeros( ( nt , nfor , nled ) )
                scores['ref_pod']= np.zeros( ( nt , nfor , nled ) )
                scores['ref_far']= np.zeros( ( nt , nfor , nled ) )
                scores['ref_bias']  = np.zeros( ( nfor , nled ) )
                scores['ref_rmse']= np.zeros( ( nfor , nled ) )
                scores['vr_bias']  = np.zeros( ( nfor , nled ) )
                scores['vr_rmse']= np.zeros( ( nfor , nled ) )
                first_read = False

            #Compute scores 
            undef = radar_grid['data_ave'].fill_value 
            
            [ets,csi,biasf,pod,far,ctable]=cf.cont_table( dfor= , dobs= , ndata= 
                    dfor=radar_grid['data_ave'][:,:,:,2].reshape( nz * ny * nx )   ,
                    dobs=radar_grid['data_ave'][:,:,:,0].reshape( nz * ny * nx )   ,
                    ndata=nz * ny * nx                                             ,
                    thresholds=thresholds_ref                                      ,
                    nthresholds=np.shape( thresholds_ref )                         ,
                    undef=undef) 

            scores['ref_ets'][:,ifor,forecast_time] = ets
            scores['ref_csi'][:,ifor,forecast_time] = csi
            scores['ref_biasf'][:,ifor,forecast_time] = biasf
            scores['ref_pod'][:,ifor,forecast_time] = pod
            scores['ref_far'][:,ifor,forecast_time] = far
            


        forecast_time = forecast_time + 1 
        cftime = cftime + dt.timedelta( seconds = out_freq )
    ifor  = ifor + 1
    ctime = ctime + dt.timedelta( seconds = init_freq )



