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
outfile=exp_path + '/time_mean/fcst_scores.pkl'

init_freq=600   #Forecast initialization frequency in seconds.
out_freq =30    #Forecast output frequency in seconds.
for_len  =1800  #Forecast length in seconds.
time_tol =15    #Maximum time tolerance (if scale time and time in nearest radar data are larger than this the interpolation is not performed.)

minref=0.0      #Ref values below this threshold will be assumed equal to the threshold.

thresholds_ref = np.array([10,20,30,40,50,60])
sizes_fss      = np.array([0,2,5,10]).astype(int)


nt=np.size( thresholds_ref )
ns=np.size( sizes_fss )

itime = dt.datetime(2013,7,13,5,5,0)  #Initial time.
etime = dt.datetime(2013,7,13,5,5,0)   #End time.

radar_files=None

nfor = int( np.floor( ( etime - itime ).seconds / init_freq ) ) + 1    #Number of forecasts.
nled = int( np.floor( for_len / out_freq ) )                    + 1    #Number of lead times.

first_read = True
scores = dict()
scores['thresholds_ref']=thresholds_ref
scores['sizes_fss'] = sizes_fss


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
            radar_grid = pkl.load(filehandler)
            #Filter invalid data 
            undef = radar_grid['data_ave'].fill_value
            radar_grid['data_ave'].mask = np.logical_or( radar_grid['data_ave'].mask , radar_grid['data_n'] == 0 )
            radar_grid['data_ave'].mask = np.logical_or( radar_grid['data_ave'].mask , np.isnan( radar_grid['data_ave'] ) )
            radar_grid['data_ave'].data[ radar_grid['data_ave'].mask  ] = undef 

            if first_read  :
                [nz,ny,nx,nv] = np.shape( radar_grid['data_ave'] )
                #Allocate memory for the scores.
                scores['ref_ets']  = np.zeros( ( nt , nfor , nled ) )
                scores['ref_biasf']= np.zeros( ( nt , nfor , nled ) )
                scores['ref_csi']  = np.zeros( ( nt , nfor , nled ) )
                scores['ref_pod']= np.zeros( ( nt , nfor , nled ) )
                scores['ref_far']= np.zeros( ( nt , nfor , nled ) )
                scores['mref_ets']  = np.zeros( ( nt , nfor , nled ) )
                scores['mref_biasf']= np.zeros( ( nt , nfor , nled ) )
                scores['mref_csi']  = np.zeros( ( nt , nfor , nled ) )
                scores['mref_pod']= np.zeros( ( nt , nfor , nled ) )
                scores['mref_far']= np.zeros( ( nt , nfor , nled ) )
                scores['ref_bias']  = np.zeros( ( nfor , nled ) )
                scores['ref_rmse']= np.zeros( ( nfor , nled ) )
                scores['ref_fss'] = np.zeros( ( nt , ns , nfor , nled ) )
                scores['rv_bias']  = np.zeros( ( nfor , nled ) )
                scores['rv_rmse']= np.zeros( ( nfor , nled ) )

                first_read = False

            #Compute scores based on contingency table
            undef = radar_grid['data_ave'].fill_value 

            #print( np.sum( (radar_grid['data_ave'].data)[:,:,:,2].reshape( nz * ny * nx ) == undef ) , nz * ny * nx , 'PEPE' )
            
            [ets,csi,biasf,pod,far,ctable]=cf.cont_table(  
                    dfor=(radar_grid['data_ave'].data)[:,:,:,2].reshape( nz * ny * nx )   ,
                    dobs=(radar_grid['data_ave'].data)[:,:,:,0].reshape( nz * ny * nx )   ,
                    ndata=nz * ny * nx                                             ,
                    thresholds=thresholds_ref                                      ,
                    nthresholds=np.size( thresholds_ref )                          ,
                    undef=undef                                                     ) 

            scores['ref_ets'][:,ifor,forecast_time] = ets
            scores['ref_csi'][:,ifor,forecast_time] = csi
            scores['ref_biasf'][:,ifor,forecast_time] = biasf
            scores['ref_pod'][:,ifor,forecast_time] = pod
            scores['ref_far'][:,ifor,forecast_time] = far


            [ets,csi,biasf,pod,far,ctable]=cf.cont_table(
                    dfor=np.max( radar_grid['data_ave'].data[:,:,:,2] , 0 ).reshape( ny * nx )   ,
                    dobs=np.max( radar_grid['data_ave'].data[:,:,:,0] , 0 ).reshape( ny * nx )   ,
                    ndata= ny * nx                                                 ,
                    thresholds=thresholds_ref                                      ,
                    nthresholds=np.size( thresholds_ref )                          ,
                    undef=undef                                                     )

            scores['mref_ets'][:,ifor,forecast_time] = ets
            scores['mref_csi'][:,ifor,forecast_time] = csi
            scores['mref_biasf'][:,ifor,forecast_time] = biasf
            scores['mref_pod'][:,ifor,forecast_time] = pod
            scores['mref_far'][:,ifor,forecast_time] = far


            #Compute Fraction Skill Score for max_dbz
            #fss   = cf.fss_det_2d( dfor=np.max( radar_grid['data_ave'][:,:,:,2] , 0 )          , 
            #                       dobs=np.max( radar_grid['data_ave'][:,:,:,0] , 0 )          ,
            #                       nx=nx                                                       , 
            #                       ny=ny                                                       , 
            #                       undef=undef                                                 , 
            #                       thresholds=thresholds_ref                                   ,
            #                       nthresholds=np.size( thresholds_ref )                       , 
            #                       sizes=sizes_fss                                             ,
            #                       nsizes=np.size( sizes_fss )                                  )
            #scores['ref_fss'][:,:,ifor,forecast_time] = fss

            #Compute RMSE for max_dbz.

            rmse = cf.com_rms(ndim = nz * ny * nx                                                                                                      ,
                              var  = (radar_grid['data_ave'].data)[:,:,:,2].reshape( nz * ny * nx ) - (radar_grid['data_ave'].data)[:,:,:,0].reshape( nz * ny * nx ) ,
                              undef= undef                                                                                                              )

            scores['ref_rmse'][ifor,forecast_time] = rmse
            scores['ref_bias'][ifor,forecast_time] = np.mean( radar_grid['data_ave'][:,:,:,2].reshape( nz * ny * nx ) - radar_grid['data_ave'][:,:,:,0].reshape( nz * ny * nx ) )

            rmse = cf.com_rms(ndim = nz * ny * nx                                                                                                      ,
                              var  = radar_grid['data_ave'][:,:,:,3].reshape( nz * ny * nx ) - radar_grid['data_ave'][:,:,:,1].reshape( nz * ny * nx ) ,
                              undef= undef                                                                                                              )

            scores['rv_rmse'][ifor,forecast_time] = rmse
            scores['rv_bias'][ifor,forecast_time] = np.mean( radar_grid['data_ave'][:,:,:,3].reshape( nz * ny * nx ) - radar_grid['data_ave'][:,:,:,1].reshape( nz * ny * nx ) )

        forecast_time = forecast_time + 1 
        cftime = cftime + dt.timedelta( seconds = out_freq )
    ifor  = ifor + 1
    ctime = ctime + dt.timedelta( seconds = init_freq )

    #Save the results
    os.makedirs( os.path.dirname( outfile ) , exist_ok=True)
    filehandler = open( outfile ,"wb")
    pkl.dump(scores,filehandler)

