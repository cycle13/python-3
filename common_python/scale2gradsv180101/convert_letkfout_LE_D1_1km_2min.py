import numpy as np
import datetime as dt
import sys
sys.path.append('./')
from scale.letkf import letkfout_grads
import time


#--- use mpi4py
from mpi4py import MPI
comm = MPI.COMM_WORLD
#--- do not use mpi4py
#comm = None
#---


sim_read = 50
nprocs = comm.Get_size()
myrank = comm.Get_rank()

if nprocs > sim_read:
    raise ValueError('The maximum number of simultaneous I/O threads is set to ' + str(sim_read) + ', please use nprocs <= ' + str(sim_read))


vcoor = 'p'
hcoor = 'o'
plevels = [100000., 95000.,90000.,80000., 85000., 70000.,60000., 50000.,40000., 30000., 20000.]
varout_3d = ['u', 'v', 'w', 'tk', 'qv','dbz']  # variables for restart file
#varout_3d = ['u', 'v', 'w', 'p', 'tk', 'theta', 'z', 'qv', 'qc', 'qr', 'qi', 'qs', 'qg', 'qhydro', 'rh', 'dbz']                                   # variables for history file
varout_2d = ['slp']                                                                     # variables for restart file
#varout_2d = ['topo', 'slp', 'u10', 'v10', 't2', 'q2', 'rain', 'snow', 'max_dbz', 'olr', 'tsfc', 'tsfcocean', 'sst', 'glon', 'glat']               # variables for history file
proj = {
'type': 'LC',
'basepoint_lon': 135.523,
'basepoint_lat': 34.823,
'basepoint_x': 90000.0,
'basepoint_y': 90000.0,
'LC_lat1': 32.5,
'LC_lat2': 37.5
}
extrap = False
dlon = 0.1
dlat = 0.1

letkfoutdir = '/home/ra001011/a03471/data/output_data/LE_D1_1km_2min/'
topofile = letkfoutdir + 'const/topo/topo'

stime = dt.datetime(2013,  7, 13,  5,  50, 0)
etime = dt.datetime(2013,  7, 13,  6,   0, 0)
tint = dt.timedelta(seconds=120)

outtype = ['gues','anal']
member = 1000

start_pp = time.time()

letkfout_grads(letkfoutdir, topofile=topofile, proj=proj, stime=stime, etime=etime, tint=tint,
               outtype=outtype, member=member,
               vcoor=vcoor, hcoor=hcoor, plevels=plevels, dlon=dlon, dlat=dlat,
               varout_3d=varout_3d, varout_2d=varout_2d, extrap=extrap,
               comm=comm, sim_read=sim_read)

print('Post proceesing took', time.time()-start_pp,'seconds. Whit ',nprocs,' processors.' )
