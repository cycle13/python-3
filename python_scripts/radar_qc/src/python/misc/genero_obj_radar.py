import glob
import os
import sys
from datetime import datetime, timedelta

try:
    from netcdftime import utime
except ImportError:
    from cftime import utime
import numpy as np
import pyart

from pyart.aux_io.sinarame_h5 import read_sinarame_h5
from pyart.aux_io.rainbow_wrl import read_rainbow_wrl

from get_strat import get_strat


def armo_radar(lista_ref_todos, path_datos='/ms-36/mrugna/RMA/datos/'):
    """
    Parametros
    ----------
    lista_ref_todos : str
        Lista de archivos de reflectividad

    path_datos : str (opcional)
        Ruta donde estan los archivos

    Devuelve
    --------
    lista_radares : lista
        Lista con instancias de objeto Radar de PyART.

    """

    lista_radares = []

    for archivo_ref in lista_ref_todos:
        nombre_radar = archivo_ref.split('/')[5] # esto esta hardcodeado al path de ms-36

        if nombre_radar.startswith('RMA'):

            # try:
            hora = datetime.strptime(archivo_ref.split('_')[-1][:13], '%Y%m%dT%H%M')
            condicion = hora.strftime(f'{path_datos}{nombre_radar}/*/*/*/*/*/*_%Y%m%dT%H%M*Z.H5')
            files = glob.glob(condicion)
            # aca arriba hay que decidir que variables se agarran

            # Ordeno la lista para que siempre arranque con T
            files.sort(reverse=True)
            files.sort(key=lambda x: len(x.split('_')[-2]))

            # print('\tEl primer archivo es:', str(files[0]), end='\n\n')

            for j, file in enumerate(files):

                if j == 0:
                    try:
                        radar = read_sinarame_h5(file, file_field_names=True)
                        azi_todos, rango = radar.fields['TH']['data'].shape
                        # print('Creo objeto radar de', file, sep='\t')

                    # Si no puede crear radar
                    # Seria recomendable agarrar las excepciones de mejor manera
                    except ValueError:
                        print('x - No se pudo crear el objeto Radar de',
                              file, sep='\t')

                    # Puede venir mal convertido el archivo H5 desde el BUFR
                    except KeyError:
                        print('x - Se convirtio mal de BUFR', file, sep='\t')

                else:
                    try:
                        radar_prov = read_sinarame_h5(file, file_field_names=True)
                        # print('Creo objeto radar de', file, sep='\t')

                        campo = list(radar_prov.fields.keys())[0]

                        # Si el archivo tiene menos gates que TH hay que agregarlos
                        # esto es util en un volumen solo con staggered prf
                        if radar_prov.fields[campo]['data'].shape[1] != rango:
                            falta = (rango
                                     - radar_prov.fields[campo]['data'].shape[1])
                            resto = np.ma.masked_all((azi_todos, falta))
                            datos = np.ma.concatenate(
                                [radar_prov.fields[campo]['data'], resto], 1)
                            radar_prov.fields[campo]['data'] = datos

                        radar.fields.update(radar_prov.fields)

                    # Puede venir mal convertido el archivo H5 desde el BUFR
                    except KeyError:
                        print('x - Se convirtio mal de BUFR', file, sep='\t')

                    # Si no puede crear radar_prov
                    except ValueError:
                        print('x - No se puede crear el objeto radar_prov para',
                              file, sep='\t')

                    # Si no encuentra radar
                    except NameError:
                        print('No existia radar, lo creo')
                        radar = read_sinarame_h5(file, file_field_names=True)

            radar = get_strat(archivo_ref, radar)

            # except UnboundLocalError as e:
            #    print(e)
            #    continue

        elif (nombre_radar == 'PAR') or (nombre_radar == 'PER') or (nombre_radar == 'ANG'):

            # try:
            hora = datetime.strptime(archivo_ref.split('/')[-1][:12], '%Y%m%d%H%M')
            condicion = hora.strftime(f'{path_datos}{nombre_radar}/*/%Y%m%d%H%M*.vol') # /240/
            files = glob.glob(condicion)
            # files = [f in files if 'tmp' not in f]
            # aca arriba hay que decidir que variables se agarran

            # Ordeno la lista para que siempre arranque con T
            files.sort(reverse=True)

            # print('\tEl primer archivo es:', str(files[0]), end='\n\n')

            for j, file in enumerate(files):

                if j == 0:
                    try:
                        radar = read_rainbow_wrl(file, file_field_names=True)
                        azi_todos, rango = radar.fields['dBZ']['data'].shape

                        #radar.fields['ZH']=radar.fields.pop['dBZ']
                        # print('Creo objeto radar de', file, sep='\t')
                        #radar = get_strat(file, radar)

                    # Si no puede crear radar
                    # Seria recomendable agarrar las excepciones de mejor manera
                    except ValueError:
                        print('x - No se pudo crear el objeto Radar de',
                              file, sep='\t')

                    # Puede venir mal convertido el archivo H5 desde el BUFR
                    except KeyError:
                        print('x - Se convirtio mal de BUFR', file, sep='\t')

                else:
                    try:
                        radar_prov = read_rainbow_wrl(file, file_field_names=True)
                        # print('Creo objeto radar de', file, sep='\t')

                        campo = list(radar_prov.fields.keys())[0]

                        #if 'V' in radar_prov :
                        #   radar_prov.fields['VRAD']=radar_prov.fields.pop['V']
                        #if 'W' in radar_prov :
                        #   radar_prov.fields['WRAD']=radar_prov.fields.pop['W']

                        #radar_prov = get_strat(file , radar_prov) 

                        ##Check if we can merge this variable with the reflectivity.
                        #if np.shape( radar.azimuth['data'] )[0] == np.shape( radar_prov.azimuth['data'] )[0] )   :
                        #   if np.sum( radar.azimuth['data'] - radar_prov.azimuth['data'] ) == 0  :
                        #      radar.fields.update(radar_prov.fields)
                        #else 
                        #      lista_radares.append(radar)

                    # Puede venir mal convertido el archivo H5 desde el BUFR
                    except KeyError:
                        print('x - Se convirtio mal de BUFR', file, sep='\t')

                    # Si no puede crear radar_prov
                    except ValueError:
                        print('x - No se puede crear el objeto radar_prov para',
                              file, sep='\t')

                    # Si no encuentra radar
                    except NameError:
                        print('No existia radar, lo creo')
                        radar = read_rainbow_wrl(file, file_field_names=True)

            radar = get_strat(archivo_ref, radar)

            # except UnboundLocalError as e:
            #    print(e)
            #    continue

        """
        Acá se normalizan los nombres de las variables
        """

        # Reflectividad RMA
        try:
            radar.fields['ZH'] = radar.fields.pop('TH')
        except:
            pass

        # Viento RMA
        # try:
        #    radar.fields['VRAD'] = radar.fields.pop('VRAD')
        # except:
        #     pass

        # RHO RMA
        # try:
        #    radar.fields['VRAD'] = radar.fields.pop('VRAD')
        # except:
        #     pass

        # Reflectividad Gematronik
        try:
            radar.fields['ZH'] = radar.fields.pop('dBZ')
        except:
            pass

        # Viento Gematronik
        try:
            radar.fields['VRAD'] = radar.fields.pop('V')
        except:
            pass

        # Ancho espectral Gematronik
        try:
            radar.fields['WRAD'] = radar.fields.pop('W')
        except:
            pass


        lista_radares.append(radar)
	
        return lista_radares
