import glob

from datetime import datetime, timedelta

# archivos_ref = []
# filtro_tmp = []
# filtro_tmp_cond = []

def busco_archivos_ref(radar, hora_patron, path_datos='/ms-36/mrugna/RMA/datos/',
                       inicio_ventana=11, fin_ventana=0):
    """
    Entrada
    -------

    radar : str
        sitio de radar
    hora_patron : datetime obj
        hora final de la ventana de busqueda de archivos
    path_datos : str
        directorio con datos H5 o vol

    Devuelve
    --------

    lista_ref : lista
        conjunto de TH/dBZ para el sitio radar en la ventana

    """

    PATH_DATOS = path_datos # alias para compatibilidad
    filtro_tmp_cond = []

    if radar.startswith('RMA'):
        # aca hay que definir el tipo de archivo (estrategia, volumen) a levantar
        filtro_tmp = glob.glob(f'{PATH_DATOS}{radar}/*/*/*/*/*/*_TH_*.H5') # *_01_TH_*.H5

        for i in filtro_tmp:
            arch = i.split('/')[-1]
            hora = datetime.strptime(arch.split('_')[-1][:15], '%Y%m%dT%H%M%S')
            cond_0 = hora_patron-timedelta(minutes=inicio_ventana) < hora < hora_patron+timedelta(minutes=fin_ventana)
            if cond_0:
                filtro_tmp_cond.append(i)

        # En esta parte reviso la cantidad de archivos que hay
        # si no hay ninguno pasa, si hay uno toma ese y si hay mas de uno
        # elige el que este mas cerca de la hora de arriba (trunc)
        """
        if len(filtro_tmp_cond) == 1:  # ESTE ES IGUAL A 1
            archivos_ref.append(filtro_tmp_cond[0])
        elif len(filtro_tmp_cond) == 0:
            print(f'No habia datos de {radar}')
        elif len(filtro_tmp_cond) > 1: # ESTE ES MAYOR A 1
            print(f'Hay mas de 1 archivo de {radar}')
            minimo = []
            for i in filtro_tmp_cond:
                arch = i.split('/')[-1]
                hora = datetime.strptime(arch.split('_')[-1][:13],
                                         '%Y%m%dT%H%M')
                valor = abs(trunc.minute - hora.minute)
                if valor < 55:
                    minimo.append(valor)
                else:
                    minimo.append(valor%50)

            val, idx = min((val, idx) for (idx, val) in enumerate(minimo))
            print(minimo)
            print(val, idx)
            archivos_ref.append(filtro_tmp_cond[idx])
        """

    elif (radar == 'PAR') or (radar == 'PER') or (radar == 'ANG'):
        filtro_tmp = glob.glob(f'{PATH_DATOS}{radar}/*/*dBZ.vol')

        for i in filtro_tmp:
            if '/tmp/' in i:
                continue
            arch = i.split('/')[-1]
            hora = datetime.strptime(arch[:14], '%Y%m%d%H%M%S')
            cond_0 = hora_patron-timedelta(minutes=inicio_ventana) < hora < hora_patron+timedelta(minutes=fin_ventana)
            if cond_0:
                filtro_tmp_cond.append(i)

        """
        if len(filtro_tmp_cond) == 1:  # ESTE ES IGUAL A 1
            archivos_ref.append(filtro_tmp_cond[0])
        elif len(filtro_tmp_cond) == 0:
            print(f'No habia datos de {radar}')
        elif len(filtro_tmp_cond) > 1: # ESTE ES MAYOR A 1
            print(f'Hay mas de 1 archivo de {radar}')
            minimo = []
            for i in filtro_tmp_cond:
                arch = i.split('/')[-1]
                hora = datetime.strptime(arch[:12], '%Y%m%d%H%M')
                valor = abs(trunc.minute - hora.minute)
                if valor < 55:
                    minimo.append(valor)
                else:
                    minimo.append(valor%50)

            val, idx = min((val, idx) for (idx, val) in enumerate(minimo))
            print(minimo)
            print(val, idx)
            archivos_ref.append(filtro_tmp_cond[idx])
        """
    # filtro_tmp_cond = []
    # filtro_tmp = []

    lista_ref = filtro_tmp_cond

    return lista_ref
