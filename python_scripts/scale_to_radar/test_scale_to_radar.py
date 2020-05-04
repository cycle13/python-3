from src.python import io
import calc
import calc_tmp
import matplotlib.pyplot as plt
import pyart


data_path='/home/juan/tmp_scale_to_radar/output_data/LE_D1_1km_30sec_OFP_V2/20130713050500/fcst/mean/history'
topo_path='/home/juan/tmp_scale_to_radar/output_data/LE_D1_1km_30sec_OFP_V2/const/topo/topo'
proj = {
'type': 'LC',
'basepoint_lon': 135.523,
'basepoint_lat': 34.823,
'basepoint_x': 90000.0,
'basepoint_y': 90000.0,
'LC_lat1': 32.5,
'LC_lat2': 37.5
}


radar = pyart.testing.make_empty_ppi_radar( 600 , 300 , 90 )
radar.latitude['data'] = 34.41
radar.longitude['data'] = 135.5
radar.altitude['data'] = 0.0
radar.ngates=600
radar.nrays=300*90



sio = io.ScaleIO(data_path, verbose=1)
sioh = io.ScaleIO(topo_path , verbose=1)

bmap = calc.set_bmap(sio, proj, resolution=None, rtol=1.e-6)
topo = sioh.readvar('TOPO',bufsize=2) 
[dbz,rv,max_dbz]=calc.calc_ref_rv(sio, bmap , radar , topo ,  min_dbz=-20., rho=None, qr=None, qs=None, qg=None, t=10, dryrun=False)



plt.pcolor( rv[10,:,:] )
plt.colorbar()
plt.show()
