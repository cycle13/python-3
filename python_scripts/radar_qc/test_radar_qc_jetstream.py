#print __doc__

# Author: Rapid Refresh Argentina Team
# License: BSD 3 clause

# History:  Created 10-2017

#Add path to additional modules.
import sys
sys.path.append('./src/python' )
sys.path.append('./src/fortran')

import numpy as np
import matplotlib.pyplot as plt
import pyart
import netCDF4
import radar_qc_module as rqc
import os

import conf

#=======================================

conf.options['toporawdatapath']="/home/jruiz/Dropbox/DATA/radar_qc/data/terrain_data/raw/"
conf.options['toporadardatapath']="/home/jruiz/Dropbox/DATA/radar_qc/data/terrain_data/radar/"


# read in the file, create a RadarMapDisplay object
filename = './cfrad.20091117_174348.000_to_20091117_174737.000_PAR_SUR.nc'

#Performs QC operations based on options
[radar , qc_output] = rqc.main_qc( filename , conf.options )

print('End of QC')

plt.figure()
plt.pcolor(qc_output['x'][:,:,0],qc_output['y'][:,:,0],qc_output['smooth_rho'][:,:,0])
plt.colorbar()
plt.show()

