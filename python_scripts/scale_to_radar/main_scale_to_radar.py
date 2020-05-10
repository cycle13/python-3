from src.python import io
from src.python import common_scale_to_radar as cs2r
from src.python import calc
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

proj = {
'type': 'LC',
'basepoint_lon': 135.523,
'basepoint_lat': 34.823,
'basepoint_x': 90000.0,
'basepoint_y': 90000.0,
'LC_lat1': 32.5,
'LC_lat2': 37.5
}

grid = {
        'xini':-60000.0,
        'nx'  :120,
        'dx'  :1000.0,
        'yini':-60000.0,
        'ny'  :120,
        'dy'  :1000.0,
        'zini':0.0,
        'nz'  :22,
        'dz'  :500.0,
        }



radar_files=None

#Read topo data
sioh = io.ScaleIO(exp_path + '/const/topo/topo' , verbose=1)
topo = sioh.readvar('TOPO',bufsize=2)
io.scale_close(sioh.rootgrps)

#Time loop to read forecast and interpolate them to nearest radar data.
ctime = itime
while ( ctime <= etime ) :

    fi_time = ctime
    fe_time = ctime + dt.timedelta( seconds = for_len )
    cftime = fi_time

    sio = io.ScaleIO( exp_path + dt.datetime.strftime( fi_time ,'%Y%m%d%H%M%S' ) + '/fcst/mean/history')

    forecast_time = 0
    while ( cftime <= fe_time ) :

        #Get radar data closest to the model time (within the tolerance)
        print('Looking for radar data close to ',dt.datetime.strftime( cftime , '%Y%m%d%H%M%S') )

        [radar , radar_files ] =  cs2r.get_radar_data( cftime , rad_path , radar_files = radar_files , time_tol = time_tol )
        
        #if False  :
        if not ( radar is None ) :  

            #Read model data and interpolate it to the radar grid.
            t0=time.time()
            radar = calc.radar_int( sio , proj , topo , radar , t=forecast_time )
            print("Model data was interpolated in {:.3f} seconds".format(time.time() - t0))
            #Mask model data according to radar data.
            radar['model_ref'][ radar['model_ref'] < radar['minref'] ] = radar['minref']
            radar['model_ref'][ radar['ref'].mask ] = radar['ref'].fill_value
            radar['model_rv'][ radar['rv'].mask ] = radar['rv'].fill_value
            radar['model_ref'] = np.ma.masked_values( radar['model_ref'] , radar['ref'].fill_value )
            radar['model_rv'] = np.ma.masked_values( radar['model_rv'] , radar['rv'].fill_value )
            #plt.pcolor( radar['lon_gate'][10,:,:] , radar['lat_gate'][10,:,:] , radar['model_ref'][10,:,:] )
            #plt.show()

            #Save the radar structure to disk
            outpath = exp_path + dt.datetime.strftime( fi_time ,'%Y%m%d%H%M%S' ) + '/fcstrad/mean/'
            os.makedirs( outpath , exist_ok=True)
            filehandler = open( outpath + '/fcstrad' + dt.datetime.strftime( cftime ,'%Y%m%d%H%M%S' ) + '.pkl',"wb")
            pkl.dump(radar,filehandler)

            #Perform data regriding.
            t0=time.time()
            radar_grid = calc.radar_regrid( radar , grid )
            print("Radar data was regrided in {:.3f} seconds".format(time.time() - t0))
            filehandler = open( outpath + '/fcstrad_grid' + dt.datetime.strftime( cftime ,'%Y%m%d%H%M%S' ) + '.pkl',"wb")
            pkl.dump(radar_grid,filehandler)
            #plt.pcolor( radar_grid['data_ave'][10,:,:,0] - radar_grid['data_ave'][10,:,:,2] );plt.colorbar();plt.show()
            #plt.pcolor( radar_grid['data_ave'][:,:,80,2]);plt.show()
            #plt.pcolor( radar_grid['data_ave'][:,:,80,0]);plt.show()
        #radar = None
        forecast_time = forecast_time + 1 
        cftime = cftime + dt.timedelta( seconds = out_freq )
    io.scale_close(sio.rootgrps)    
    ctime = ctime + dt.timedelta( seconds = init_freq )


#plt.pcolor( radar.fields['model_ref']['data'][10,:,:] )
#plt.show()

#plt.pcolor( radar.fields['model_rv']['data'][3000:3300,:] )
#plt.show()

