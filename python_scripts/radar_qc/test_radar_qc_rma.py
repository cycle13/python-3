#print __doc__
# Author: Rapid Refresh Argentina Team
# License: BSD 3 clause
# History:  Created 10-2017

#Add path to additional modules.
import sys
import matplotlib.pyplot as plt
sys.path.append('./src/python' )
sys.path.append('./src/fortran')

import radar_qc_module as rqc
import conf_defaults as conf
import pyart

#=======================================

# read in the file, create a RadarMapDisplay object
#filename='./cfrad.20091117_174348.000_to_20091117_174737.000_PAR_SUR.nc3'    #Ejemplo Parana
conf.options['filename']='./cfrad.20170926_171335.0000_to_20170926_171446.0000_RMA1_0122_03.nc3'  #Ejemplo Cordoba

conf.options['filename_out']='./cfrad.20170926_171335.0000_to_20170926_171446.0000_RMA1_0122_03.corrected.nc3'

#Performs QC operations based on options

radar = pyart.io.read(conf.options['filename'])

[radar , qc_output] = rqc.main_qc( conf.options )


print('End of QC')



