from src.python import io
from src.python import pawr_read
import calc
#import calc_tmp
import matplotlib.pyplot as plt
import pyart
import datetime as dt
import numpy as np


exp_path='/home/juan/tmp_scale_to_radar/output_data/LE_D1_1km_30sec_OFP_V2/'
rad_path='/home/juan/tmp_scale_to_radar/radar_data/'

init_freq=600   #Forecast initialization frequency in seconds.
out_freq =30    #Forecast output frequency in seconds.
for_len  =1800  #Forecast length in seconds.

minref=0.0      #Ref values below this threshold will be assumed equal to the threshold.

itime = dt.datetime(2013,7,13,5,10,0)  #Initial time.
etime = dt.datetime(2013,7,13,6,0,0)   #End time.



data_path='/home/juan/tmp_scale_to_radar/output_data/LE_D1_1km_30sec_OFP_V2/20130713050500/fcst/mean/history'
topo_path='/home/juan/tmp_scale_to_radar/output_data/LE_D1_1km_30sec_OFP_V2/const/topo/topo'
radar_file='/home/juan/tmp_scale_to_radar/radar_data/PAWR_LETKF_INPUTV4_20130713-145040.dat'
proj = {
'type': 'LC',
'basepoint_lon': 135.523,
'basepoint_lat': 34.823,
'basepoint_x': 90000.0,
'basepoint_y': 90000.0,
'LC_lat1': 32.5,
'LC_lat2': 37.5
}


#Time loop to read forecast and interpolate them to nearest radar data.
ctime = itime
while ( ctime <= etime ):

    fi_time = ctime
    fe_time = ctime + dt.timedelta( seconds = for_len )
    cftime = fi_time




    ctime = ctime + dt.timedelta( seconds = init_freq )

#Time loop to read forecast and interpolate them to nearest radar data.

#Generate file object.
sio = io.ScaleIO(data_path, verbose=1)
sioh = io.ScaleIO(topo_path , verbose=1)
#Read topo
topo = sioh.readvar('TOPO',bufsize=2) 
#Read radar
radar = pawr_read.pawr_read( radar_file )
#Read and interpolate model data
[model_ref , model_rv ]=calc.radar_int( sio , proj , topo , radar , t=10 ) 

model_ref[ model_ref < minref ] = minref

#Add interpolated data to radar object.
undef = radar.fields['ref']['data'].fill_value
model_ref[ radar.fields['ref']['data'].mask  ] = undef
model_rv [ radar.fields['rv']['data'].mask   ] = undef

ref_dict = pyart.config.get_metadata('reflectivity')
ref_dict['data'] = np.ma.masked_values(model_ref, undef)
rv_dict = pyart.config.get_metadata('velocity')
rv_dict['data'] = np.ma.masked_values(model_rv, undef)
radar.fields['model_ref']=ref_dict
radar.fields['model_rv'] =rv_dict


plt.pcolor( radar.fields['model_ref']['data'][3000:3300,:] )
plt.show()

plt.pcolor( radar.fields['model_rv']['data'][3000:3300,:] )
plt.show()

