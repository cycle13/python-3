MODULE common_scale_to_radar
!=======================================================================
!
! [PURPOSE:] Observational procedures
!
! [HISTORY:]
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_scale
  USE common_obs_scale
  USE common_namelist

  IMPLICIT NONE
  PUBLIC

  real(r_size),allocatable :: v3d(:,:,:,:)
  real(r_size),allocatable :: v2d(:,:,:)

  real(r_size)  :: min_ref_dbz
  
  CHARACTER*1000  :: model_file_prefix     !Model input file. T slots are for different times (or ensemble members)
  CHARACTER*1000  :: topo_file_prefix      !Model topography input file prefix.
  CHARACTER*1000  :: current_model_prefix 
      
CONTAINS


!-----------------------------------------------------------------------
! Interpolate model data to a radar (or part of a radar) domain

!-----------------------------------------------------------------------
SUBROUTINE model_to_radar( radar_na , radar_nr , radar_lat , radar_lon , radar_z , radar_az ,  & 
                           radar_el , model_time , model_ref , movel_rv )
  IMPLICIT NONE
  INTEGER      , INTENT(IN)   :: radar_na , radar_nr , model_time
  REAL(r_size) , INTENT(IN)   :: radar_lat(radar_na,radar_nr) , radar_lon(radar_na,radar_nr) 
  REAL(r_size) , INTENT(IN)   :: radar_z(radar_na,radar_nr) , radarar_az(radar_na) , radarar_el(radar_na)
  REAL(r_size) , INTENT(OUT)  :: model_ref(radar_na,radar_nr) , model_rv(radar_na,radar_nr)

  INTEGER              :: ia , ir , ie , i , j , k
  REAL(r_size)         :: ri , rj , rk , rlev
  REAL(r_size)         :: tmp_z ,  tmp
  REAL(r_size)         :: p , t
  REAL(r_size)         :: qv  , qc , qr
  REAL(r_size)         :: qi  , qs , qg 
  REAL(r_size)         :: u   , v  , w
  REAL(r_size)         :: ref , rv ,  cref , refdb 
  REAL(r_size)         :: pik !Path integrated attenuation coefficient.
  REAL(r_size)         :: prh !Pseudo relative humidity
  LOGICAL              :: ISALLOC
  REAL(r_size)         :: max_model_z
  REAL(r_size)         :: tmpref,tmperr(1)
  REAL(r_size)         :: test_lon , test_lat

  min_ref_dbz=10.0*log10(min_ref)
  ref_out = undef
  vr_out  = undef  

  !Begin with the interpolation. 

  max_model_z = 0.0d0
  do i = 1 , nlonh
    do j = 1 , nlath 
       if( v3d(nlev,i,j,iv3d_hgt) > max_model_z )then
         max_model_z = v3d(nlev,i,j,iv3d_hgt)
       endif
    enddo
  enddo


!$OMP PARALLEL DO DEFAULT(SHARED) FIRSTPRIVATE(ia,ir,pik,tmp_z,ri,rj,rk,qv,qc,qr,qci,qs,qg,t,p,u,v,w,ref,rv,cref,prh)

    DO ia=1,input_radar%na 
      pik=0.0d0 !Path integrated attenuation coefficient.
      DO ir=1,nr 
 
        !PAWR can have data up to 60.000 m we speed up the forward operator
        !by checking this simple condition.  
        IF( radar_z(ia,ir) .LE. max_model_z )THEN

        CALL  phys2ij(radar_lon(ia,ir),radar_lat(ia,ir),ri,rj)
        CALL  phys2ijkz(v3d(:,:,:,iv3d_hgt),ri,rj,radar_z(ia,ir),rk)
 

        if( rk /= undef )THEN
   
           !Interpolate qv,qc,qr,qci,qs,qg,t and p to the beam center.

           CALL itpl_3d(v3d(:,:,:,iv3d_u),rk,ri-0.5d0,rj,u)  
           CALL itpl_3d(v3d(:,:,:,iv3d_v),rk,ri,rj-0.5d0,v) 
           CALL itpl_3d(v3d(:,:,:,iv3d_w),rk-0.5d0,ri,rj,w) 
           CALL itpl_3d(v3d(:,:,:,iv3d_t),rk,ri,rj,t)
           CALL itpl_3d(v3d(:,:,:,iv3d_p),rk,ri,rj,p)
           CALL itpl_3d(v3d(:,:,:,iv3d_q),rk,ri,rj,qv)
           CALL itpl_3d(v3d(:,:,:,iv3d_qc),rk,ri,rj,qc)
           CALL itpl_3d(v3d(:,:,:,iv3d_qr),rk,ri,rj,qr)
           CALL itpl_3d(v3d(:,:,:,iv3d_qi),rk,ri,rj,qi)
           CALL itpl_3d(v3d(:,:,:,iv3d_qs),rk,ri,rj,qs)
           CALL itpl_3d(v3d(:,:,:,iv3d_qg),rk,ri,rj,qg)

           !TODO ROTATE WINDS
 
            !Compute reflectivity at the beam center.
           CALL calc_ref_rv(qv,qc,qr,qi,qs,qg,u,v,w,t,p,            &
               radar_az(ia),radar_el(ia)    &
                ,ref,rv)
   
           !ADD ERRORS TO THE OBSERVATIONS
           IF( ADD_OBS_ERROR )THEN
             !Add error to reflectivity.
             IF( ref > min_ref)THEN
              tmpref=10.0d0*log10(ref)
              CALL com_randn(1,tmperr)
              tmpref=tmpref+tmperr(1)*REFLECTIVITY_ERROR
              ref=10.0d0**(tmpref/10.0d0) 
              IF( ref < min_ref )ref=min_ref
             ENDIF
        
             !Add error to doppler velocity.
             CALL com_randn(1,tmperr)
             rv=rv+tmperr(1)*RADIALWIND_ERROR 

           ENDIF

          refdb=10.0d0*log10(ref)

          IF( ref .GT. min_ref )THEN
              ref_out(ir,ia)=refdb !refdb
          ELSE
              ref_out(ir,ia)=min_ref_dbz
          ENDIF

          ref_rv(ir,ia)=rv    !Radial wind


        ENDIF  !Domain check

       ENDIF   !Max height check

      ENDDO
    ENDDO
!$OMP END PARALLEL DO


  ENDIF
  


  RETURN
END SUBROUTINE model_to_radar


SUBROUTINE scale2rad( 
                      METHOD_REF_CALC_in ,
                      min_ref_in , compute_attenuation_in      ,
                      input_type_in , 
                      input_type_radar_in , 
                      MPRJ_type_in , 
                      MPRJ_basepoint_x_in , MPRJ_basepoint_y_in ,
                      MPRJ_rotation_in , 
                      MPRJ_LC_lat1_in , MPRJ_LC_lat2_in ,
                      MPRJ_PS_lat_in , MPRJ_M_lat_in , MPRJ_EC_lat_in ,
                      MPRJ_basepoint_lon_in , MPRJ_basepoint_lat_in   ,
                      model_time             !The corresponding time in the model input file.
                      radar_na , radar_nr ,  !Number of azimuths and ranges (cfradial format)
                      radar_lat , radar_lon , radar_z ,
                      radar_file_in , model_input_prefix_in , model_topo_prefix_in , outputh_path )


IMPLICIT NONE
!This is a python wrapper for the interp_to_real_radar subroutine. It can be used to call the scale to radar program from a python
!script directly.
!Namelist variables
INTEGER,INTENT(IN)            :: METHOD_REF_CALC_in      !=2        !Method used for reflectivity computaton.
REAL(r_size),INTENT(IN)       :: min_ref_in              ! = 1.0d0          !Minimum reflectivity value (in Power units)
INTEGER,INTENT(IN)            :: input_type_in          !=1  !1-restart files, 2-history files.
INTEGER,INTENT(IN)            :: input_type_radar_in    !=1  !1-cfradial, 2-binary files.

!Projection parameters
character(len=4),INTENT(IN) :: MPRJ_type_in       ! = 'NONE' !< map projection type
real(r_size),INTENT(IN)  :: MPRJ_basepoint_x_in        ! position of base point in the model [m]
real(r_size),INTENT(IN)  :: MPRJ_basepoint_y_in        ! position of base point in the model [m]
real(r_size),INTENT(IN)  :: MPRJ_rotation_in           ! =  0.0d0 ! rotation factor (only for 'NONE' type)
real(r_size),INTENT(IN)  :: MPRJ_LC_lat1_in            !  = 30.0d0 ! standard latitude1 for L.C. projection [deg]
real(r_size),INTENT(IN)  :: MPRJ_LC_lat2_in            !  = 60.0d0 ! standard latitude2 for L.C. projection [deg]
real(r_size),INTENT(IN)  :: MPRJ_PS_lat_in             !             ! standard latitude1 for P.S. projection [deg]
real(r_size),INTENT(IN)  :: MPRJ_M_lat_in              !    =  0.0d0 ! standard latitude1 for Mer. projection [deg] 
real(r_size),INTENT(IN)  :: MPRJ_EC_lat_in             !   =  0.0d0 ! standard latitude1 for E.C. projection [deg]
real(r_size),INTENT(IN)  :: MPRJ_basepoint_lon_in      ! = 135.221d0 ! position of base point
real(r_size),INTENT(IN)  :: MPRJ_basepoint_lat_in      ! =  34.653d0 ! position of base point

INTEGER, INTENT(IN)        :: model_time 

INTEGER, INTENT(IN)        :: radar_na , radar_nr
REAL(r_size) , INTENT(IN)  :: radar_lat(radar_na,radar_nr) , radar_lon(radar_na,radar_nr) 
REAL(r_size) , INTENT(IN)  :: radar_z(radar_na,radar_nr) , radarar_az(radar_na) , radarar_el(radar_na)
REAL(r_size) , INTENT(IN)  :: model_ref(radar_na,radar_nr) , model_rv(radar_na,radar_nr)

METHDO_REF_CALC=METHOD_REF_CALC_in 
min_ref=min_ref_in
compute_attenaution=compute_attenuation_in      
model_time = model_time_in 
  
input_type=input_type_in 
MPRJ_type=MPRJ_type_in 
MPRJ_basepoint_x=MPRJ_basepoint_x_in 
MPRJ_basepoint_y=MPRJ_basepoint_y_in 
MPRJ_rotation=MPRJ_rotation_in 
MPRJ_LC_lat1=MPRJ_LC_lat1_in 
MPRJ_LC_lat2=MPRJ_LC_lat2_in 
MPRJ_PS_lat=MPRJ_PS_lat_in 
MPRJ_M_lat=MPRJ_M_lat_in 
MPRJ_EC_lat=MPRJ_EC_lat_in 
MPRJ_basepoint=MPRJ_basepoint_lon_in 
MPRJ_basepoint_lat=MPRJ_basepoint_lat_in  
model_file_prefix = model_input_prefix_in              !Model input file. 
radar_file =  radar_file_in                            !Radar data input in cfradial format
topo_file_prefix = model_topo_prefix_in                !Model topography input file prefix.

!Get grid properties.
write(*,*)"Getting model info"

CALL set_common_scale( model_file_prefix , topo_file_prefix , input_type )

if ( model_time > ntime )then
   write(*,*)"Error: number of times in files is less than the requested time.
   stop
endif

ALLOCATE( v3d(nlev,nlonh,nlath,nv3d),v2d(nlonh,nlath,nv2d) )
ALLOCATE( model_ref(na,nr),model_rv(na,nr) )
write(*,*)"Reading model data"
CALL read_file(model_file_prefix,v3d,v2d,model_time,input_type)
write(*,*)"Interpolating model data to radar grid"
CALL model_to_radar( radar_na , radar_nr , radar_lat , radar_lon , radar_z , radar_az , &
                     radar_el , model_time , model_ref , movel_rv )


END SUBROUTINE scale2rad


END MODULE common_scale_to_radar
