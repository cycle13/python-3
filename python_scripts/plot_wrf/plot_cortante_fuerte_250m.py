#Este script plotea el empuje y sus componentes en un corte vertical.
#Las funciones estan pensadas para leer un diccionario (wrf_data) que contiene todos los datos necesarios.
#Cada funcion va agregando nuevas variables a ese diccionario.

#Las funciones buscan reducir el tiempo y no repetir calculos inecesariamente, por eso se fijan si el calculo ya esta hecho.
#Pero tambien esta la opcion de setear force=True para repetir el calulo a la fuerza.

import wrf_module as wrf
import wrf_plot   as wrfp
import numpy as np
import os

#Indicamos la carpeta donde estan los wrfout (el script asume que cada tiempo de salida del WRF esta en un archivo separado)
#va a buscar todos los archivos en wrfout en la carpeta y a ordenarlos cronologicamente.
exp_path = '/home/jruiz/datosalertar1/SIMULACIONES_CYM_2017/cortante_fuerte_250m/'    #Carpeta base del experimento.

data_path = exp_path + '/run/'                                                       #Carpeta donde estan los datos.

#A donde van a ir a parar las figuras.
plot_path = exp_path + '/figuras/'                                                   #Carpeta donde se generan las figuras.
os.makedirs(plot_path,exist_ok=True)

#Leo el primer archivo de la lista para obtener el perfil vertical de T en x=0, y=0 y con eso construir el grafico que representa el entorno.
wrf_data = wrf.get_profile( data_path , xp=0 , yp=0 , tp=0 )
#Grafica el perfil vertical de T, Td y viento en un emagrama.
wrfp.plot_sounding( wrf_data , plot_path , show=False )

#Corte en Y

#Esta funcion lee todos los wrfouts y saca todas las variables que se necesitan para hacer cuentas. 
#Estas variables se almacenan en una "slice" que es un corte vertical que puede ser en la direccion de y o de x (vx o vy).
#slice_index dice cual es el punto de reticula donde realiza el corte. t_start y t_end se especifican si se quiere graficar solo un rango
#de tiempos en lugar de toda las simulacion.
#Lo que se guarda para cada variable es un corte vertical centrado en slice_index pero que tiene un ancho de 3 puntos de reticula de manera
#tal de poder calcular sobre eso tanto derivadas en x como en y. Esta funcion tarda un poco la primera vez, pero guarda el corte que obtuvo en un
#archivo pickle con lo cual si le pedis nuevamente el mismo corte levanta el pickle en lugar de leer la salida del wrf.
#wrf_data es un diccionario que tiene como entradas las variables del wrf (y algunos parametros que guardan la metadata del corte y el tipo de corte).
#cada entrada correspondiente a una variable es un array de numpy y se preserva los niveles verticales originales del wrf asi como la resolucion horizontal.
wrf_data = wrf.get_data_vslice( data_path , slice_type='vy' , slice_index = 160 , t_start=5 , t_end=25 ) 

#Todo este grupete de funciones (get_xxxx_equation) trabajan sobre el wrf_data y le agregan entradas que tienen que ver con los terminos de las diferentes ecuaciones.
#ahi se calculan las derivadas verticales, horizontales, etc. El parametro save indica que guarde el wrf_data en un pickle (el mismo que genero get_data_vslice) pero 
#agregando las nuevas variables que se van calculando, de manera que si se quiere graficar un corte para el cual los calculos ya fueron hechos nos podemos ahorrar todo este paso.
wrf_data = wrf.get_moment_equation( wrf_data , save=True )

wrf_data = wrf.get_termo_equation( wrf_data , save=True )

wrf_data = wrf.get_water_equation( wrf_data , save=True )

wrf_data = wrf.get_ppert_equation(  wrf_data , save=True )

wrf_data = wrf.get_vorticity_equation(  wrf_data , save=True )


#Los plot_xxxx_equation son las funciones que generan las figuras. Todo se hace bastante automaticamente a partir del diccionario wrf_data. 
#ybound y xbound permiten controlar la extension de los ejes (en este caso como la slice es un corte vertical ybound controla la extension en z).
#Tambien hay un parametro optativo scale_factor que permite ajustar la escala, por ejemplo si pongo scale_factor = 2.0 usa una escala que tiene el doble de extensión que 
#la que viene por defecto y eso se aplica a todas las variables salvo la reflectividad.
wrfp.plot_momentum_equation_v( wrf_data , plot_path , ybound=[0,15000] , xbound=[20000,50000]  )

wrfp.plot_termo_equation_v( wrf_data , plot_path , ybound=[0,15000] , xbound=[20000,50000] )

wrfp.plot_vapor_equation_v( wrf_data , plot_path , ybound=[0,15000] , xbound=[20000,50000] )

wrfp.plot_water_equation_v( wrf_data , plot_path , ybound=[0,15000] , xbound=[20000,50000] )

wrfp.plot_thetas_v( wrf_data , plot_path , ybound=[0,15000] , xbound=[20000,50000] )

wrfp.plot_ppert_equation_v( wrf_data , plot_path , ybound=[0,15000] , xbound=[20000,50000] )

wrfp.plot_vortz_equation_v( wrf_data , plot_path  , ybound=[0,15000] , xbound=[20000,50000] )


#Corte en X. Este bloque es esencialmente lo mismo pero para un corte vertical a x constante en lugar de y constante.
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
#Esto es un ejemplo de corte a z constante. A diferencia de los otros cortes este usa el wrfpython para interpolar las variables del wrf
#a la altura solicitada (slice_z) y a una altura slice_z - zres y slice_z + zres que usa para calcular las derivadas verticales. Por ese motivo
#estos datos no respetan los niveles verticales originales del wrf, pero si la resolución en x e y. En este caso slice_type = 'h' indica que se trata
#de un corte horizontal a z constante.
wrf_data = wrf.get_data_vslice( data_path , slice_type='h' , slice_z = 3000.0 , zres=200.0 , t_start=5 , t_end=25 )

#Lo demas es todo bastante igual, las funciones que calculan las ecuaciones cambian poco dependiendo que sea un corte horizontal o vertical.
#el unico cambio tiene mas que ver con la eficiencia del calculo (las derivadas verticales son mucho mas faciles de calcular en este caso) que con otra cosa.
wrf_data = wrf.get_moment_equation( wrf_data , save=True )

wrf_data = wrf.get_termo_equation( wrf_data , save=True )

wrf_data = wrf.get_water_equation( wrf_data , save=True )

wrf_data = wrf.get_ppert_equation(  wrf_data , save=True )

wrf_data = wrf.get_vorticity_equation(  wrf_data , save=True )

#Para los graficos se usan las mismas funciones con los mismos parametros de entrada. Como wrf_data tiene la informacion de si es un corte
#horizontal o vertical la funcion se acomoda automaticamente a una cosa o la otra (en teoria...). 
wrfp.plot_momentum_equation_v( wrf_data , plot_path , ybound=[25000,55000] , xbound=[20000,50000]  )

wrfp.plot_termo_equation_v( wrf_data , plot_path , ybound=[25000,55000] , xbound=[20000,50000] )

wrfp.plot_vapor_equation_v( wrf_data , plot_path , ybound=[25000,55000] , xbound=[20000,50000] )

wrfp.plot_water_equation_v( wrf_data , plot_path , ybound=[25000,55000] , xbound=[20000,50000] )

wrfp.plot_ppert_equation_v( wrf_data , plot_path , ybound=[25000,55000] , xbound=[20000,50000] )

wrfp.plot_vortz_equation_v( wrf_data , plot_path  , ybound=[25000,55000] , xbound=[20000,50000] )


#Corte a 5000 metros
#Esto es un ejemplo de corte a z constante. A diferencia de los otros cortes este usa el wrfpython para interpolar las variables del wrf
#a la altura solicitada (slice_z) y a una altura slice_z - zres y slice_z + zres que usa para calcular las derivadas verticales. Por ese motivo
#estos datos no respetan los niveles verticales originales del wrf, pero si la resolución en x e y. En este caso slice_type = 'h' indica que se trata
#de un corte horizontal a z constante.
wrf_data = wrf.get_data_vslice( data_path , slice_type='h' , slice_z = 5000.0 , zres=200.0 , t_start=5 , t_end=25 )

#Lo demas es todo bastante igual, las funciones que calculan las ecuaciones cambian poco dependiendo que sea un corte horizontal o vertical.
#el unico cambio tiene mas que ver con la eficiencia del calculo (las derivadas verticales son mucho mas faciles de calcular en este caso) que con otra cosa.
wrf_data = wrf.get_moment_equation( wrf_data , save=True )

wrf_data = wrf.get_termo_equation( wrf_data , save=True )

wrf_data = wrf.get_water_equation( wrf_data , save=True )

wrf_data = wrf.get_ppert_equation(  wrf_data , save=True )

wrf_data = wrf.get_vorticity_equation(  wrf_data , save=True )

#Para los graficos se usan las mismas funciones con los mismos parametros de entrada. Como wrf_data tiene la informacion de si es un corte
#horizontal o vertical la funcion se acomoda automaticamente a una cosa o la otra (en teoria...). 
wrfp.plot_momentum_equation_v( wrf_data , plot_path , ybound=[25000,55000] , xbound=[20000,50000]  )

wrfp.plot_termo_equation_v( wrf_data , plot_path , ybound=[25000,55000] , xbound=[20000,50000] )

wrfp.plot_vapor_equation_v( wrf_data , plot_path , ybound=[25000,55000] , xbound=[20000,50000] )

wrfp.plot_water_equation_v( wrf_data , plot_path , ybound=[25000,55000] , xbound=[20000,50000] )

wrfp.plot_ppert_equation_v( wrf_data , plot_path , ybound=[25000,55000] , xbound=[20000,50000] )

wrfp.plot_vortz_equation_v( wrf_data , plot_path  , ybound=[25000,55000] , xbound=[20000,50000] )







