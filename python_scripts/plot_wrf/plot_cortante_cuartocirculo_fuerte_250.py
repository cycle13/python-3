#Este script plotea el empuje y sus componentes en un corte vertical.
#Las funciones estan pensadas para leer un diccionario (wrf_data) que contiene todos los datos necesarios.
#Cada funcion va agregando nuevas variables a ese diccionario.

#Las funciones buscan reducir el tiempo y no repetir calculos inecesariamente, por eso se fijan si el calculo ya esta hecho.
#Pero tambien esta la opcion de setear force=True para repetir el calulo a la fuerza.

import wrf_module as wrf
import wrf_plot   as wrfp
import numpy as np
import os


exp_path = '/home/jruiz/datosalertar1/SIMULACIONES_CYM_2017/cortante_cuartocirculo_fuerte_250m/'    #Carpeta base del experimento.

data_path = exp_path + '/run/'                                                       #Carpeta donde estan los datos.

plot_path = exp_path + '/figuras/'                                                   #Carpeta donde se generan las figuras.
os.makedirs(plot_path,exist_ok=True)

wrf_data = wrf.get_profile( data_path , xp=0 , yp=0 , tp=0 )

wrfp.plot_sounding( wrf_data , plot_path , show=False )

#Corte en Y
wrf_data = wrf.get_data_vslice( data_path , slice_type='vy' , slice_index = 160 , t_start=5 , t_end=25 ) 

wrf_data = wrf.get_moment_equation( wrf_data , save=True )

wrf_data = wrf.get_termo_equation( wrf_data , save=True )

wrf_data = wrf.get_water_equation( wrf_data , save=True )

wrf_data = wrf.get_ppert_equation(  wrf_data , save=True )

wrf_data = wrf.get_vorticity_equation(  wrf_data , save=True )

wrfp.plot_momentum_equation_v( wrf_data , plot_path , ybound=[0,15000] , xbound=[20000,50000]  )

wrfp.plot_termo_equation_v( wrf_data , plot_path , ybound=[0,15000] , xbound=[20000,50000] )

wrfp.plot_vapor_equation_v( wrf_data , plot_path , ybound=[0,15000] , xbound=[20000,50000] )

wrfp.plot_water_equation_v( wrf_data , plot_path , ybound=[0,15000] , xbound=[20000,50000] )

wrfp.plot_thetas_v( wrf_data , plot_path , ybound=[0,15000] , xbound=[20000,50000] )

wrfp.plot_ppert_equation_v( wrf_data , plot_path , ybound=[0,15000] , xbound=[20000,50000] )

wrfp.plot_vortz_equation_v( wrf_data , plot_path  , ybound=[0,15000] , xbound=[20000,50000] )


#Corte en X
wrf_data = wrf.get_data_vslice( data_path , slice_type='vx' , slice_index = 160 , t_start=5 , t_end=25 )              

wrf_data = wrf.get_moment_equation( wrf_data , save=True )

wrf_data = wrf.get_termo_equation( wrf_data , save=True )

wrf_data = wrf.get_water_equation( wrf_data , save=True )

wrf_data = wrf.get_ppert_equation(  wrf_data , save=True )

wrf_data = wrf.get_vorticity_equation(  wrf_data , save=True )

wrfp.plot_momentum_equation_v( wrf_data , plot_path , ybound=[0,15000] , xbound=[30000,50000]  )

wrfp.plot_termo_equation_v( wrf_data , plot_path , ybound=[0,15000] , xbound=[30000,50000] )

wrfp.plot_vapor_equation_v( wrf_data , plot_path , ybound=[0,15000] , xbound=[30000,50000] )

wrfp.plot_water_equation_v( wrf_data , plot_path , ybound=[0,15000] , xbound=[30000,50000] )

wrfp.plot_thetas_v( wrf_data , plot_path , ybound=[0,15000] , xbound=[30000,50000] )

wrfp.plot_ppert_equation_v( wrf_data , plot_path , ybound=[0,15000] , xbound=[30000,50000] )

wrfp.plot_vortz_equation_v( wrf_data , plot_path  , ybound=[0,15000] , xbound=[30000,50000] )

#Corte a 3000 metros
wrf_data = wrf.get_data_vslice( data_path , slice_type='h' , slice_z = 3000.0 , zres=200.0 , t_start=5 , t_end=25 )

wrf_data = wrf.get_moment_equation( wrf_data , save=True )

wrf_data = wrf.get_termo_equation( wrf_data , save=True )

wrf_data = wrf.get_water_equation( wrf_data , save=True )

wrf_data = wrf.get_ppert_equation(  wrf_data , save=True )

wrf_data = wrf.get_vorticity_equation(  wrf_data , save=True )

wrfp.plot_momentum_equation_v( wrf_data , plot_path , ybound=[25000,55000] , xbound=[20000,50000]  )

wrfp.plot_termo_equation_v( wrf_data , plot_path , ybound=[25000,55000] , xbound=[20000,50000] )

wrfp.plot_vapor_equation_v( wrf_data , plot_path , ybound=[25000,55000] , xbound=[20000,50000] )

wrfp.plot_water_equation_v( wrf_data , plot_path , ybound=[25000,55000] , xbound=[20000,50000] )

wrfp.plot_ppert_equation_v( wrf_data , plot_path , ybound=[25000,55000] , xbound=[20000,50000] )

wrfp.plot_vortz_equation_v( wrf_data , plot_path  , ybound=[25000,55000] , xbound=[20000,50000] )


#Corte a 5000 metros
wrf_data = wrf.get_data_vslice( data_path , slice_type='h' , slice_z = 5000.0 , zres=200.0 , t_start=5 , t_end=25 )

wrf_data = wrf.get_moment_equation( wrf_data , save=True )

wrf_data = wrf.get_termo_equation( wrf_data , save=True )

wrf_data = wrf.get_water_equation( wrf_data , save=True )

wrf_data = wrf.get_ppert_equation(  wrf_data , save=True )

wrf_data = wrf.get_vorticity_equation(  wrf_data , save=True )

wrfp.plot_momentum_equation_v( wrf_data , plot_path , ybound=[25000,55000] , xbound=[20000,50000]  )

wrfp.plot_termo_equation_v( wrf_data , plot_path , ybound=[25000,55000] , xbound=[20000,50000] )

wrfp.plot_vapor_equation_v( wrf_data , plot_path , ybound=[25000,55000] , xbound=[20000,50000] )

wrfp.plot_water_equation_v( wrf_data , plot_path , ybound=[25000,55000] , xbound=[20000,50000] )

wrfp.plot_ppert_equation_v( wrf_data , plot_path , ybound=[25000,55000] , xbound=[20000,50000] )

wrfp.plot_vortz_equation_v( wrf_data , plot_path  , ybound=[25000,55000] , xbound=[20000,50000] )









