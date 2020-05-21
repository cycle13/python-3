#Este script plotea el empuje y sus componentes en un corte vertical.
#Las funciones estan pensadas para leer un diccionario (wrf_data) que contiene todos los datos necesarios.
#Cada funcion va agregando nuevas variables a ese diccionario.

#Las funciones buscan reducir el tiempo y no repetir calculos inecesariamente, por eso se fijan si el calculo ya esta hecho.
#Pero tambien esta la opcion de setear force=True para repetir el calulo a la fuerza.

import wrf_module as wrf
import wrf_plot   as wrfp
import numpy as np
import os


exp_path = '../cortante_cuartocirculo_fuerte_250m'    #Carpeta base del experimento.

data_path = exp_path + '/run/'                                                       #Carpeta donde estan los datos.

plot_path = exp_path + '/figuras/'                                                   #Carpeta donde se generan las figuras.
os.makedirs(plot_path,exist_ok=True)


#wrf_data = wrf.get_data_vslice( data_path , slice_type='vy' , slice_index = 80 , force=False ) 

#wrf_data = wrf.get_data_vslice( data_path , slice_type='vy' , slice_index = 80 , force=False ) 

wrf_data = wrf.get_data_vslice( data_path , slice_type='h' , slice_z = 500.0 , zres=200.0 , slice_index = 80 , force=False ) 

# wrf_data = wrf.get_moment_equation( wrf_data , force=False , save = True )

# wrf_data = wrf.get_termo_equation( wrf_data , force=False , save = True )

# wrf_data = wrf.get_water_equation( wrf_data , force=False , save = True )

# wrf_data = wrf.get_ppert_equation(  wrf_data , force=True , save=True )

wrf_data = wrf.get_vorticity_equation(  wrf_data , force=True , save=True )

#wrfp.plot_momentum_equation_v( wrf_data , plot_path  , ybound=[20000,40000] , xbound=[30000,50000] , scale_factor = 1.0 )

wrfp.plot_vortz_equation_v( wrf_data , plot_path  , ybound=[20000,40000] , xbound=[30000,50000] , scale_factor = 1.0 )

#wrfp.plot_momentum_equation_2_v( wrf_data , plot_path , ybound=[11900.0,19900.0] , xbound=[11900.0,19900.0], scale_factor = 1.0 )

#wrfp.plot_termo_equation_2_v( wrf_data , plot_path , ybound=[11900.0,19900.0] , xbound=[11900.0,19900.0], scale_factor = 1.0 )

#wrfp.plot_vapor_equation_2_v( wrf_data , plot_path , ybound=[11900.0,19900.0] , xbound=[11900.0,19900.0] , scale_factor = 1.5 )

#wrfp.plot_water_equation_2_v( wrf_data , plot_path , ybound=[11900.0,19900.0] , xbound=[11900.0,19900.0] , scale_factor = 1.5 )

#wrfp.plot_thetas_v( wrf_data , plot_path , ybound=[0,12500] , scale_factor = 1.5 )

#wrfp.plot_ppert_equation_v( wrf_data , plot_path , ybound=[11900.0,19900.0] , xbound=[11900.0,19900.0] , scale_factor = 1.0 )

#wrfp.plot_vortz_equation_v( wrf_data , plot_path , ybound=[11900.0,19900.0] , xbound=[11900.0,19900.0] , scale_factor = 1.0 )

