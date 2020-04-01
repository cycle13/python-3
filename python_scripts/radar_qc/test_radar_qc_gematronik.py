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

import qc_plot_tools as qcpt

import conf_gematronik as conf

#=======================================


# read in the file, create a RadarMapDisplay object
filename='./cfrad.20091117_174348.000_to_20091117_174737.000_PAR_SUR.nc3'    #Ejemplo Parana
#filename='./cfrad.20170926_171335.0000_to_20170926_171446.0000_RMA1_0122_03.nc3'  #Ejemplo Cordoba

#Performs QC operations based on options
[radar , qc_output] = rqc.main_qc( filename , conf.options )

print('End of QC')


elev=3

if conf.options['Dealiasing']['flag']  :
   qcpt.plot_dealiasing( qc_output , conf.options , figname='DealiasingTest.png' , elev=elev)
if conf.options['DopplerNoiseFilter']['flag']  :
   qcpt.plot_dopplernoisefilter( qc_output , conf.options , figname='DopplerNoiseFilterTest.png' , elev=elev )
if conf.options['DopplerLocalStdFilter']['flag']  :
   qcpt.plot_dopplerlocalstdfilter( qc_output , conf.options , figname='DopplerLocalStdFilterTest.png' , elev=elev )
if conf.options['DopplerSpatialCoherenceFilter']['flag']  :
   qcpt.plot_dopplerspatialcoherencefilter( qc_output , conf.options , figname='DopplerSpatialCoherenceFilterTest.png' , elev=elev , show = True)
if conf.options['RhoFilter']['flag']   :
   qcpt.plot_rhofilter( qc_output , conf.options , figname='RhoFilterTest.png' , elev=elev )
if conf.options['EchoTopFilter']['flag'] :
   qcpt.plot_echotopfilter( qc_output , conf.options , figname='EchoTopFilter.png' , elev=elev )
if conf.options['EchoDepthFilter']['flag'] :
   qcpt.plot_echodepthfilter( qc_output , conf.options , figname='EchoDepthFilter.png' , elev=elev )
if conf.options['RefSpeckleFilter']['flag'] :
   qcpt.plot_refspecklefilter( qc_output , conf.options , figname='RefSpeckleFilter.png' , elev=elev )
if conf.options['DopplerSpeckleFilter']['flag'] :
   qcpt.plot_dopplerspecklefilter( qc_output , conf.options , figname='DopplerSpeckleFilter.png' , elev=elev )
if conf.options['DopplerTextureFilter']['flag'] :
   qcpt.plot_dopplertexturefilter( qc_output , conf.options , figname='DopplerTextureFilter.png' , elev=elev )
if conf.options['ReflectivityTextureFilter']['flag'] :
   qcpt.plot_reflectivitytexturefilter( qc_output , conf.options , figname='ReflectivityTextureFilter.png' , elev=elev )
if conf.options['AttenuationFilter']['flag'] :
   qcpt.plot_attenuationfilter( qc_output , conf.options , figname='AttenuationFilter.png' , elev=elev )
if conf.options['BlockingFilter']['flag'] :
   qcpt.plot_blockingfilter( qc_output , conf.options , figname='BlockingFilter.png' , elev=elev )
if conf.options['LowElevFilter']['flag'] :
   qcpt.plot_lowelevfilter( qc_output , conf.options , figname='LowElevFilter.png' , elev=elev )
if conf.options['LowDopplerFilter']['flag'] :
   qcpt.plot_lowdopplerfilter( qc_output , conf.options , figname='LowDopplerFilter.png' , elev=elev )
if conf.options['InterferenceFilter']['flag'] :
   qcpt.plot_interferencefilter( qc_output , conf.options , figname='InterferenceFilter.png' , elev=elev )







