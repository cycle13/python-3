#Este script plotea el empuje y sus componentes en un corte vertical.
#Las funciones estan pensadas para leer un diccionario (wrf_data) que contiene todos los datos necesarios.
#Cada funcion va agregando nuevas variables a ese diccionario.

#Las funciones buscan reducir el tiempo y no repetir calculos inecesariamente, por eso se fijan si el calculo ya esta hecho.
#Pero tambien esta la opcion de setear force=True para repetir el calulo a la fuerza.

import wrf_module as wrf
import wrf_plot   as wrfp
import numpy as np
import os


exp_path = '../celda_ordinaria_200m'    #Carpeta base del experimento.

data_path = exp_path + '/run/'                                                       #Carpeta donde estan los datos.

plot_path = exp_path + '/figuras/'                                                   #Carpeta donde se generan las figuras.
os.makedirs(plot_path,exist_ok=True)


wrf_data = wrf.get_profile( data_path , xp=0 , yp=0 , tp=0 ) 


wrfp.plot_sounding( wrf_data , plot_path , show=False )

