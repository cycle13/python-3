import numpy as np
import numpy.ma as ma
import os

#Este modulo lee un ctl de grads para obtener parametros imporantes para la lectura de los datos

def read_ctl( filename , coding='utf-8' )  :

	fp=open(filename,'rb')
	ctl_ori=fp.read().decode(coding)
	ctl=dict()

	#Read some parameters and dimensions.

	ctl['template']=False
	ctl['big_endian']=False
	ctl['yrev']=False
	ctl['sequential']=False
	use_xdef=True
	for line in ctl_ori.split('\n'):
		if ( 'options' in line  or 'OPTIONS' in line ) and  ( 'template' in line or 'TEMPLATE' in line )  :
			ctl['template']=True
		if ( 'options' in line or 'OPTIONS' in line ) and  ( 'big_endian' in line or 'BIG_ENDIAN' in line )  :
			ctl['big_endian']=True
		if ( 'options' in line or 'OPTIONS' in line ) and  ( 'byteswapped' in line or 'BYTESWAPPED' in line )  :
			ctl['big_endian']=True

		if ( 'options' in line or 'OPTIONS' in line ) and  ( 'sequential' in line or 'SEQUENTIAL' in line )  :
			ctl['sequential']=True
		if ( 'options' in line or 'OPTIONS' in line ) and  ( 'yrev' in line or 'YREV' in line )  :
			ctl['yrev']=True

		if 'pdef' in line or 'PDEF' in line :
			ctl['nx']=int(line.split()[1])
			ctl['ny']=int(line.split()[2])
			ctl['dx']=float(line.split()[11])
			ctl['dy']=float(line.split()[12])
			use_xdef=False
		if ( 'xdef' in line or 'XDEF' in line ) and use_xdef :
			ctl['nx']=int(line.split()[1])  
		if ( 'ydef' in line or 'YDEF' in line ) and use_xdef :
			ctl['ny']=int(line.split()[1])  
		if 'zdef' in line or 'ZDEF' in line  :
			ctl['nz']=int(line.split()[1])
			ctl['vcoordtype']=line.split()[2]

		if 'tdef' in line or 'TDEF' in line  :
			ctl['nt']=int(line.split()[1])
		if ( 'vars' in line or 'VARS' in line ) and not ( 'ENDVARS' in line or 'endvars' in line ) :
			ctl['nv']=int(line.split()[1])
		if 'undef' in line or 'UNDEF' in line :
			ctl['undef']=float(line.split()[1])

	#Read variable list
	ctl['var_list']=list()
	ctl['var_size']=list()
	ctl['var_desc']=list()

	var_section=False
	for line in ctl_ori.split('\n'):

		if ( 'ENDVARS' in line or 'endvars' in line ) :
			var_section=False
		if var_section                                :
			#Get the var name and store it in a var_list list	
			ctl['var_list'].append(line.split()[0])
			ctl['var_size'].append(line.split()[1])
			ctl['var_desc'].append(["".join(line.split()[2:])])
		if ( 'VARS' in line or 'vars' in line ) and not ( 'ENDVARS' in line or 'endvars' in line ) :
			var_section=True

	ctl['vlevels']=np.zeros(ctl['nz']) 
	vlev_section=False
	lev_counter=0
	for line in ctl_ori.split('\n'):
		if  vlev_section  :
			#Get  the level
			ctl['vlevels'][lev_counter]=float(line.split()[0])
			lev_counter=lev_counter + 1
			if lev_counter == ctl['nz']   :
				vlev_section=False


		if 'zdef' in line or 'ZDEF' in line  :
			if 'levels' in ctl['vcoordtype']    :  
				if np.size( line.split() ) == 3     :  #The levels are not in this line.
					#We need to run the following nz lines to input the levels.
					vlev_section=True
				else                            :  #The levels are in the same line.
					for ilev in range(0,ctl['nz'])   :
						ctl['vlevels'][ilev]=float( line.split()[ilev+3] )
 

	record=0
	ctl['ini_record']=list()
	ctl['end_record']=list()
	for var_size in ctl['var_size']   :
		ctl['ini_record'].append(record)
		if int(var_size) == 0    :
			var_size = 1
		record=record + int(var_size) 
		ctl['end_record'].append(record-1)

	ctl['ini_record']=np.array(ctl['ini_record'])
	ctl['end_record']=np.array(ctl['end_record'])

	#The following section will add two keys to the dictionary:
	#The added keys are:
	#full_var_list -> variable corresponding to each record
	#full_lev_list -> level corresponding to each record

	ctl['full_var_list']=list()
	ctl['full_lev_list']=list()
	ctl['full_rec_list']=np.arange(0,ctl['end_record'][-1])

	#record=0
	for ivar in range(0,int(ctl['nv']))  :
		ilev = 0 
		if int(ctl['var_size'][ivar]) >  0    :
			for ilev in range(0,int( ctl['var_size'][ivar] ))   :
				ctl['full_var_list'].append( ctl['var_list'][ivar] )
				ctl['full_lev_list'].append( ctl['vlevels'][ilev] )
  
				ilev = ilev + 1
				#record = record + 1
		elif int(ctl['var_size'][ivar]) == 0	 :	
			ctl['full_var_list'].append( ctl['var_list'][ivar] )
			ctl['full_lev_list'].append( ctl['vlevels'][0] )
         
	ctl['full_lev_list']=np.array(ctl['full_lev_list'])

	return ctl

def record_subset( ctl , subset_var , subset_lev )    :

	#Given a ctl file, a subset of variables and levels
	#this function returns the subset of records to be read and 
	#the variables and levels corresponding to each record.

	subset_var_list=list()
	subset_lev_list=list()
	subset_record_list=list()
                
	for irecord in ctl['full_rec_list']    :
		if ctl['full_var_list'][irecord] in subset_var  and np.in1d(ctl['full_lev_list'][irecord],subset_lev)  :

			subset_var_list.append( ctl['full_var_list'][irecord]  )
			subset_lev_list.append( ctl['full_lev_list'][irecord]  )
			subset_record_list.append( ctl['full_rec_list'][irecord] )

	return subset_var_list , np.asarray( subset_lev_list ) , np.asarray( subset_record_list )

def read_data_grads(filename,ctl,masked=False,undef2nan=False,singletime=True):
        #If single time == True then nt will be ignored and one time will be read.

	my_data=dict()

	nx=int(ctl['nx'])
	ny=int(ctl['ny'])
	if singletime :
		nt=1 
	else          :      
		nt=int(ctl['nt'])

	undef=ctl['undef']



	tmp_data=read_data( filename , ctl , undef2nan=undef2nan , singletime = singletime )  #Read all the data.

	#Loop over variables to create the dictionary. 
 
	for it in range(0,nt)      :
		ivar=0	
		for my_var in ctl['var_list']   :

			if it == 0        :
				nzvar=int(ctl['var_size'][ivar])
				if nzvar == 0    :
					nzvar=1
				my_data[my_var]=np.ones([ny,nx,nzvar,nt]).astype(np.float32)

			tmp_data_2=(tmp_data[:,:,ctl['ini_record'][ivar]:ctl['end_record'][ivar]+1,it])

			if masked                      :	
				rec=ctl['ini_record'][ivar]
				for iz in range(0,nzvar)  :
                        
					my_data[my_var][:,:,iz,it]=ma.masked_array( tmp_data[:,:,rec,it] , mask= tmp_data == undef )
					rec=rec+1
 
			else                           :
				rec=ctl['ini_record'][ivar]
				for iz in range(0,nzvar)  :

					my_data[my_var][:,:,iz,it]=tmp_data[:,:,rec,it]
					rec=rec+1

			ivar=ivar+1

	tmp_data_2=None

	tmp_data=None

	return my_data

def read_data(  inputfilename , ctl  , undef2nan=False , singletime = True ):
#dtypein is the native input data format.
#>f32 for big endian data format with single precission.
#f32 for little endian data with single precission ... and so on.

#If singletime = True only one time will be readed and nt from ctl will be ignored.

    nx=int(ctl['nx'])
    ny=int(ctl['ny'])
    nz=int(ctl['end_record'].max() + 1)
    if singletime :
       nt = 1 
    else          :
       nt=int(ctl['nt'])
    undef=np.float32( ctl['undef'] )
    seq_acces=ctl['sequential']

    field=np.float32( np.ones([ny,nx,nz,nt])*undef )

    if  ctl['big_endian']  :
        dtypein = '>f4'
    else                           :
        dtypein = 'f4'
    if  os.path.exists(inputfilename) :
        f=open(inputfilename,'r')
        for it in range(0,nt) :
            for ii in range(0,nz) :
                if seq_acces :
                    nada=np.fromfile(f,dtype='>i4',count=1)
                field[:,:,ii,it]= np.fromfile(f,dtype=dtypein,count=nx*ny).reshape(nx,ny) 
                if seq_acces :
                    nada=np.fromfile(f,dtype='>i4',count=1)
    else :
        print('Not found ',inputfilename)
        #If the file does not exist we will consider the entire data volume as missing data.

    if undef2nan   :
           field[ field == undef ] = np.nan

    return field

def read_data_records(  inputfilename , ctl , records ):
#dtypein is the native input data format.
#>f32 for big endian data format with single precission.
#f32 for little endian data with single precission ... and so on.

    nx=int(ctl['nx'])
    ny=int(ctl['ny'])
    undef=ctl['undef']
    seq_acces=ctl['sequential']
    max_record=int(np.max(records))
    n_record=int(np.size(records))
    records=records

    field=np.ones([ny,nx,n_record]) * undef 

    if  ctl['big_endian']  :
        dtypein = '>f4'
    else                           :
        dtypein = 'f4'
    if  os.path.exists(inputfilename) :
        f=open(inputfilename,'r')
        current_record = 0
        for ir in range(0,max_record + 1) :
                if seq_acces :
                    nada=np.fromfile(f,dtype='>i4',count=1)

                tmp = np.fromfile(f,dtype=dtypein,count=nx*ny).reshape(nx,ny)
                if seq_acces :
                    nada=np.fromfile(f,dtype='>i4',count=1)

                if np.any( records == ir )  :
                    field[:,:,current_record] = tmp
                    current_record = current_record + 1

    else :
        print('Not found ',inputfilename)

    return field


def write_data(inputfilename,field,dtypein):
#dtypein is the desired data format.
#>f32 for big endian data format with single precission.
#f32 for little endian data with single precission ... and so on.
    tmp_shape=field.shape
    nx=tmp_shape[0]
    ny=tmp_shape[1]
    nz=tmp_shape[2]

    #Reorder data to be consistent with input format.
    field=field.transpose(2,0,1)
    field=np.reshape(field,nx*ny*nz)
 
    #Write the data
    field.astype(dtypein).tofile(inputfilename)

def get_var_start_end( ctl , my_var )  :
#Given a var, if this var is in the ctl var list, return initial and end record 
#If the variable is not found a warning is printed and -1 is returned for both start and end variable record index.
    for ii , my_ctl_var in enumerate( ctl['var_list'] )  :
        if my_ctl_var == my_var :
            return ctl['ini_record'][ii] ,  ctl['end_record'][ii]
    print('Warning: variable '+my_var+' not in ctl file list')
    return -1.0 , -1.0


