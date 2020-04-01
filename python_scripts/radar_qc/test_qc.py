#print __doc__
# Author: Rapid Refresh Argentina Team
# License: BSD 3 clause
# History:  Created 10-2017

import numpy as np
import os 
import glob
import sys
sys.path.append('./src/python' )
sys.path.append('./src/fortran')

import radar_qc_module as rqc
import conf_defaults_test as conf

#La idea es hacer un loop sobre un grupo de carpetas / subcarpetas y archivos.
#En este caso asumo que hay 2 niveles de carpetas (por eso el /*/*/)
 
#datapath='/media/jruiz/PAWR/RADAR_DATA_RRA_RMA1/RMA1/'
 

##Generate a list with all the files contained in the folders and sub-folders.
#file_list=[]
     
#for (dirpath, dirnames, filenames) in os.walk(datapath):


#    for filename in filenames:
#        if ( os.path.basename(filename).find('.nc') > 0 ) and ( os.path.basename(filename).find('.corr.') < 0 ) :
#           f = '/'.join([dirpath,filename])
#           file_list.append(f)


#for filename in filenames:
#        if ( os.path.basename(filename).find('.nc') > 0 ) and ( os.path.basename(filename).find('.corr.') < 0 ) :
#           f = '/'.join([dirpath,filename])
#           file_list.append(f)
#file_list=['/media/jruiz/PAWR/RADAR_DATA_RRA_RMA1/RMA1/20170923/cfrad.20170923_114117.0000_to_20170923_114310.0000_RMA1_0122_02.nc']
#file_list=['/media/jruiz/PAWR/RADAR_DATA_RRA_RMA1/RMA1/20170926/cfrad.20170926_211512.0000_to_20170926_212043.0000_RMA1_0122_01.nc']
#file_list=['/media/jruiz/PAWR/RADAR_DATA_RRA_RMA1/ANG/20170930/cfrad.20170930_133002.000_to_20170930_133420.001_ANG_SUR.nc']
#file_list=['/media/jruiz/PAWR/RADAR_DATA_RRA_RMA1/PAR/20170930//cfrad.20170930_222434.000_to_20170930_222732.000_PAR_SUR.nc']

#file_list=['./cfrad.20170926_171335.0000_to_20170926_171446.0000_RMA1_0122_03.nc3']

#file_list=['./cfrad.20181008_181228.RMA2.nc']

file_list=['./RMA1_0200_02_TH_20181010T114222Z.H5']

toporawdatapath="/home/jruiz/share/radar_qc/data/terrain_data/raw/"
toporadardatapath="/home/jruiz/share/radar_qc/data/terrain_data/radar/"



options = conf.options #This is the default configuration.

#Proceed to perform qc for each file.
for ifile in file_list  :
   #Test if we have an RMA or GEMATRONIK RADAR
   if os.path.basename(ifile).find('RMA') > 0  :

      options['name_ref'] ='TH'              #Reflectivity
      options['name_v']   ='VRAD'            #Dopper velocity
      options['name_rho'] ='RHOHV'           #Rho HV
      options['is_rma'] = True               #Wether this is an RMA file

   else                   :

      options['name_ref'] ='dBZ'              #Reflectivity
      options['name_v']   ='V'                #Dopper velocity
      options['name_rho'] ='RhoHV'            #Rho HV
      options['is_rma'] = False                         #Wether this is an RMA file

   options['filename'] = ifile
   #Generate output file name
   options['filename_out'] = ifile[ 0:ifile.find('.nc') ] + '.corr' + '.nc' 
   options['plot']['Path'] = './'
   #os.makedirs( options['plot']['Path'],exist_ok=True)
   print('Performing qc for the following file')
   print(options['filename'])
   print('Output will be stored in the following file')
   print(options['filename_out'])
   if options['plot']['Enable']   :
      print('Figures will be created in the following file')
      print(options['plot']['Path'])
   print('Is this an RMA radar?')
   print(options['is_rma'])
   options['toporawdatapath']=toporawdatapath
   options['toporadardatapath']=toporadardatapath

   #Performs QC operations based on options
   [radar , qc_output] = rqc.main_qc( options )


