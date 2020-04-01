#!/home/qcradar/.conda/envs/da/bin/python

import os
import glob
from datetime import date, datetime, timedelta

import pyart

from busco_archivos_ref import busco_archivos_ref
from genero_obj_radar import armo_radar

qc_path = "/home/qcradar/scripts/radar_qc/"

import sys
sys.path.append( qc_path + '/src/python/' )
sys.path.append( qc_path + '/src/fortran/')

import radar_qc_module as rqc        #Radar qc python modules
import conf_defaults as conf         #Radar qc default configuration


"""
Esta parte del script es la que busca los archivos para una determinada hora

"""

PATH_DATOS = '/ms-36/mrugna/RMA/datos/'


# LA HORA DE SOL1 TIENE PROBLEMAS CON EL HUSO HORARIO
utc = timedelta(hours=3)
# 20 MINUTOS DE RETRASO ES MUCHO. 10 ESTA BIEN
ahora = (datetime.now() - timedelta(minutes=20) + utc).replace(second=0,
                                                               microsecond=0)

trunc = ahora - timedelta(minutes=ahora.minute%10)

print(ahora, trunc)


lista_radares = ['RMA1',# 'RMA2', 'RMA3', 'RMA4', 'RMA5', 'RMA6',
                 # 'RMA7', 'RMA8',
                 # 'PAR',
                 'PER', 'ANG']


"""
Ahora busco SOLAMENTE los TH/dBZ

"""

lista_ref_todos = []

for sitio_radar in lista_radares:
    lista_ref = busco_archivos_ref(sitio_radar, trunc, PATH_DATOS) # el path de datos es opcional

    lista_ref_todos.extend(lista_ref)
    lista_ref = []

print()
for j in lista_ref_todos:
    print(j)

"""
Para cada uno de los TH/dBZ busco el resto de las variables de interes
y genero un unico objeto radar

"""

lista_radares = armo_radar(lista_ref_todos)

for radar in lista_radares:

# llama a qc

    import conf_defaults as conf
    options = conf.options #This is the default configuration.

    options['name_ref'] ='ZH'         #Reflectivity
    options['name_v']   ='VRAD'       #Dopper velocity
    options['name_rho'] ='RHOHV'      #Rho HV
    options['toporawdatapath']= qc_path + "/data/terrain_data/raw/"
    options['toporadardatapath']= qc_path + "/data/terrain_data/radar/" 

    [ radar , qc_output ] = rqc.main_qc( options , radar )


# manda qc_ref a superobbing??
