# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:45:15 2016

@author:
"""

import sys
sys.path.append('../../common_python/common_modules/')

import numpy as np
import matplotlib.pyplot as plt
import os

import common_plot_functions as cpf
import common_gdal_functions as cgf

testfile='/volume64/data/ra001011/jruiz/input_data/raster_tiles/ASTGTM2_N34E133_dem.tif'


[lon,lat,data]=cgf.read_raster_data(testfile) 

plt.pcolor( lon[0:1000,0:1000] , lat[0:1000,0:1000] , np.squeeze( data[0:1000,0:1000] ) )

plt.show()



