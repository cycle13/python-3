MODULE QC_CONST
!MODULE QC
 IMPLICIT NONE
 PUBLIC
!=======================================================================
!
! [PURPOSE:]  Constants for QC
!
! [HISTORY:]
!   09/01/2014 Juan Ruiz created
!
!=======================================================================
!-----------------------------------------------------------------------
! Variable size definitions
!-----------------------------------------------------------------------
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_dble=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
!-----------------------------------------------------------------------
! Constants
!-----------------------------------------------------------------------
  REAL(r_size),PARAMETER :: pi=3.141592653589793d0
  REAL(r_size),PARAMETER :: gg=9.81d0
  REAL(r_size),PARAMETER :: rd=287.0d0
  REAL(r_size),PARAMETER :: cp=7.0d0 / 2.0d0 * rd
  REAL(r_size),PARAMETER :: re=6371.3d3
  REAL(r_size),PARAMETER :: r_omega=7.292d-5
  REAL(r_size),PARAMETER :: t0c=273.15d0
  REAL(r_size),PARAMETER :: deg2rad=3.1415926535d0/180d0
  REAL(r_size),PARAMETER :: rad2deg=180d0/3.1415926535d0
  REAL(r_size),PARAMETER :: clight=299792458.0d0 !Speed of light

!  INTEGER      , PARAMETER :: NPAR_ECHO_TOP_3D=6 , NPAR_ECHO_TOP_2D=7 !Number of parameters in output arrays.
!  REAL(r_size) , PARAMETER :: MAX_Z_ECHO_TOP=20.0d4 , MAX_R_ECHO_TOP=240.0d03
!  REAL(r_size) , PARAMETER :: DZ_ECHO_TOP = 500.0d0 , DX_ECHO_TOP = 500.0d0

!  INTEGER      , PARAMETER :: MAX_ECHO_TOP_LEVS=5 
!  REAL(r_size) , PARAMETER :: DBZ_THRESHOLD_ECHO_TOP=5.0d0  !Echo top detection value.

!  REAL(r_size)             :: undef 


!  REAL(r_size) , ALLOCATABLE  :: QCARRAY(:,:,:)  !Array to store qccodes

END MODULE QC_CONST

MODULE QC
!=======================================================================
!
! [PURPOSE:] Common quality control parameters computation
!
! [HISTORY:]
!   09/01/2014 Juan Ruiz created
!
!=======================================================================
!$USE OMP_LIB
  USE qc_const
  IMPLICIT NONE
  PUBLIC

 CONTAINS


SUBROUTINE SPECKLE_FILTER(var,na,nr,ne,undef,nx,ny,nz,threshold,speckle)
IMPLICIT NONE
INTEGER, INTENT(IN)        :: na,nr,ne,nx,ny,nz        !Var dims and box dims
REAL(r_size),INTENT(IN)    :: threshold                !Threshold 
REAL(r_size),INTENT(IN)    :: var(na,nr,ne)            !Input variable
REAL(r_size),INTENT(IN)    :: undef
REAL(r_size),INTENT(OUT)   :: speckle(na,nr,ne)        !Temporal array
INTEGER                    :: ia,ir,ie 

  CALL BOX_FUNCTIONS_2D(var,na,nr,ne,undef,nx,ny,nz,'COUN',threshold,speckle)

RETURN
END SUBROUTINE SPECKLE_FILTER

SUBROUTINE RHO_FILTER(var,na,nr,ne,undef,nx,ny,nz,rho_smooth)
IMPLICIT NONE
INTEGER, INTENT(IN)        :: na,nr,ne,nx,ny,nz        !Var dims and box dims
!REAL(r_size),INTENT(IN)    :: threshold                !Threshold 
REAL(r_size),INTENT(INOUT) :: var(na,nr,ne)            !Input variable
REAL(r_size),INTENT(IN)    :: undef                    !Undefined value.
REAL(r_size),INTENT(OUT)   :: rho_smooth(na,nr,ne)     !Temporal array
INTEGER                    :: ia,ir,ie

  CALL BOX_FUNCTIONS_2D(var,na,nr,ne,undef,nx,ny,nz,'MEAN',0.0d0,rho_smooth)


RETURN

END SUBROUTINE RHO_FILTER


SUBROUTINE GET_ATTENUATION(var,na,nr,ne,undef,beaml,cal_error,coefs,is_power,attenuation,corrected_var,mindbz)

!==========================================================================================
!  This function estimates the attenuation percentaje due to metereological
!  echoes and also computes a correction.
!  Input:
!  radar an structure containing radar information.
!  reflectivity: a 3D array (possibly 4D) containing reflectivity data with
!  most of the ground clutter already removed (we will assume that the
!  reflectivity is associatedi with weather echoes only.
!  Output:
!  attenuation which is the Path Integrated Attenuation (PIA) A. Berne and R.Uijlenhoet
!  2006.
!  correction is the correction (in dbz) that has to be added to each grid
!  point. 
!  To avoid the very well known instability of the forward attenuation
!  computation algorithm the algorithm is stopped when the attenuation
!  factor reaches a certain threshold. 
!===========================================================================================

IMPLICIT NONE
INTEGER     ,INTENT(IN)    :: na,nr,ne
REAL(r_size),INTENT(IN)    :: var(na,nr,ne) !Input reflectivity in dBZ
REAL(r_size),INTENT(OUT)   :: corrected_var(na,nr,ne) !Output corrected reflectivity in dBZ
REAL(r_size),INTENT(IN)    :: undef
LOGICAL     ,INTENT(IN)    :: is_power      !If input data is in mm^6/m^3 the set this to true.
REAL(r_size),INTENT(IN)    :: coefs(4)
REAL(r_size),INTENT(OUT)   :: attenuation(na,nr,ne) !Attenuation factor.
REAL(r_size),INTENT(IN)    :: beaml  !Beam length (m)
REAL(r_size),INTENT(IN)    :: mindbz 
REAL(r_size)               :: a_coef , b_coef , c_coef , d_coef , alfa , beta  !Attenuation parameters
REAL(r_size),INTENT(IN)    :: cal_error !Calibration erro (use 1.0 if we dont know it)
REAL(r_size)               :: power(nr) , attenuation_power(nr)
INTEGER                    :: ia,ir,ie,ir2
REAL(r_size)               :: mean_k
REAL(r_size),PARAMETER     :: max_dbz_correction = 10.0d0


!Coefficients for C-Band radars based on
!Quantification of Path-Integrated Attenuation for X- and C-Band Weather
!Radar Systems Operating in Mediterranean Heavy Rainfall
!Delrieu, Andreiu, Creutin 1999 , Journal of Applied Meteorology
!a_coef=543
!b_coef=1.36
!c_coef=1.55e-3
!d_coef=1.30

a_coef=coefs(1)
b_coef=coefs(2)
c_coef=coefs(3)
d_coef=coefs(4)

!Formulation of PIA is based on the equations presented in.
!Quantitative analysis of X-band weather radar attenuation
!correction accuracy
!A. Berne and R. Uijlenhoet , 2006
!Nat. Hazards Earth Syst. Sci.

alfa=(c_coef**d_coef)/(a_coef**(d_coef/b_coef))
beta=(d_coef/b_coef)

!tmp_data_3d=10.0d0**(var/10)
!where( var == undef )
!     tmp_data_3d=0.0d0
!endwhere

!==========================================================================
! Iterative algorithm to compute attenuation. Based on the forward
! attenuation estimation of HB.
!==========================================================================

!We iterate forward in the range direction.
attenuation(:,1,:)=1.0d0*cal_error

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(ia,ie,ir,ir2,mean_k,power,attenuation_power)
DO ia=1,na
 DO ie=1,ne
  power = 10.0d0**(var(ia,:,ie)/10.0d0)
  attenuation_power = 1.0d0
  where( var(ia,:,ie) == undef ) 
      power = 0.0d0
  endwhere

    DO ir=1,nr-1

       mean_k=1.0d-3*(alfa*( (0.5)*(power(ir)+power(ir+1)) )**beta ) 
       !Compute mean k between ir and ir+1 (k is dbz/m);
       attenuation_power(ir+1) = attenuation_power(ir) *exp(-0.46d0*mean_k*beaml)
       !attenuation(ia,ir+1,ie)=attenuation(ia,ir,ie)*exp(-0.46d0*mean_k*beaml)

       attenuation(ia,ir+1,ie) = 10.0d0*log10( attenuation_power(ir+1) )

       !If attenuation is les than max_dbz_correction, then correct attenuation.
       if ( attenuation(ia,ir+1,ie) > -1.0*max_dbz_correction )then
          do ir2=ir+1,nr
             if ( var(ia,ir2,ie) > mindbz ) then
                power(ir2) = power(ir2) / exp(-0.46d0*mean_k*beaml)
             endif
          enddo 
       endif 

    ENDDO
  
    DO ir = 1 , nr  
      IF ( power( ir ) > 0.0d0 )THEN
         corrected_var(ia,ir,ie) = 10.0d0*log10( power(ir) )
      ELSE
         corrected_var(ia,ir,ie) = undef 
      ENDIF
   ENDDO
 ENDDO
ENDDO
!$OMP END PARALLEL DO


END SUBROUTINE GET_ATTENUATION


SUBROUTINE COMPUTE_TEXTURE(var,na,nr,ne,undef,nx,ny,nz,texture)
!This routine performs the radar QC computing the requested fields.
IMPLICIT NONE
INTEGER     ,INTENT(IN) :: na , nr , ne    !Grid dimension
INTEGER     ,INTENT(IN) :: nx , ny , nz  !Box dimension
REAL(r_size),INTENT(IN) :: var(na,nr,ne) 
!REAL(r_size),INTENT(IN)     :: threshold
REAL(r_size),INTENT(IN)     :: undef
REAL(r_size),INTENT(OUT)    :: texture(na,nr,ne)
REAL(r_size)             :: tmp_data_3d(na,nr,ne) 
INTEGER                  :: ii , jj , kk

!Compute the difference along the radial direction.
tmp_data_3d=undef

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(ii,jj,kk)
 DO ii = 1,na
   DO jj = 1 ,nr-1
     DO kk = 1, ne
        IF( var(ii,jj,kk) /= undef .AND. var(ii,jj+1,kk) /= undef)THEN
          tmp_data_3d(ii,jj,kk) = ( var(ii,jj+1,kk)-var(ii,jj,kk) )**2
        ENDIF
     ENDDO
   ENDDO
 ENDDO 
!$OMP END PARALLEL DO

 !Average the squared radial differences.
 CALL BOX_FUNCTIONS_2D(tmp_data_3d,na,nr,ne,undef,nx,ny,nz,'MEAN',0.0d0,texture)

RETURN
END SUBROUTINE COMPUTE_TEXTURE

SUBROUTINE COMPUTE_SIGN(var,na,nr,ne,undef,nx,ny,nz,varsign)
!This routine computes the sign parameter
!Kessinger et al 2003
IMPLICIT NONE
INTEGER     ,INTENT(IN)     :: na , nr , ne    !Grid dimension
INTEGER     ,INTENT(IN)     :: nx , ny , nz  !Box dimension
REAL(r_size),INTENT(INOUT)  :: var(na,nr,ne)
REAL(r_size),INTENT(IN)     :: undef
!REAL(r_size),INTENT(IN)     :: threshold 
REAL(r_size)                :: tmp_data_3d(na,nr,ne) , diff
INTEGER                     :: ii , jj , kk
REAL(r_size),INTENT(OUT)    :: varsign(na,nr,ne)

!Compute the difference along the radial direction.
tmp_data_3d=undef

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(ii,jj,kk,diff)
 DO ii = 1,na
   DO jj = 1 ,nr-1
     DO kk = 1, ne
        IF( var(ii,jj,kk) /= undef .AND. var(ii,jj+1,kk) /= undef )THEN
           diff= var(ii,jj,kk) - var(ii,jj+1,kk) 
           IF( ABS(diff) > 0 )THEN
             tmp_data_3d(ii,jj,kk) = diff / ABS(diff)
           ELSE
             tmp_data_3d(ii,jj,kk) = 0.0d0
           ENDIF
        ENDIF
     ENDDO
   ENDDO
 ENDDO
!$OMP END PARALLEL DO

 !Average the squared radial differences.
 CALL BOX_FUNCTIONS_2D(tmp_data_3d,na,nr,ne,undef,nx,ny,nz,'MEAN',0.0d0,varsign)

! where( varsign > threshold .or. varsign == undef )
!    var=undef
!    qcarray=QCCODE_SIGN
! endwhere 


RETURN
END SUBROUTINE COMPUTE_SIGN

SUBROUTINE BOX_FUNCTIONS_2D(datain,na,nr,ne,undef,boxx,boxy,boxz,operation,threshold,dataout)

IMPLICIT NONE
INTEGER     ,INTENT(IN) :: na , nr , ne    !Grid dimension
INTEGER     ,INTENT(IN) :: boxx,boxy,boxz  !Box dimension
REAL(r_size),INTENT(IN) :: datain(na,nr,ne)
REAL(r_size),INTENT(IN) :: undef
CHARACTER(4),INTENT(IN) :: operation     
REAL(r_size),INTENT(IN) :: threshold
REAL(r_size),INTENT(OUT) :: dataout(na,nr,ne) !Result
REAL(r_size),ALLOCATABLE :: tmp_field(:) 
REAL(r_size)             :: tmp_mean , tmp_var
INTEGER                  :: NITEMS
INTEGER                  :: ii , jj , kk , bii , bjj , bkk , box_size , iin ,ii_index , data_count

box_size=(2*boxx+1)*(2*boxy+1)*(2*boxz+1)

!DO kk=1,ne
!WRITE(*,*)maxval( datain(:,:,kk) )
!ENDDO

ALLOCATE( tmp_field(box_size) )

dataout=undef

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(kk,ii,jj,bkk,bii,bjj,ii_index,NITEMS,tmp_field,tmp_mean,tmp_var,data_count,iin)
DO kk=1,ne
 DO ii=1,na
  DO jj=1,nr


     NITEMS=0
     tmp_field=0.0d0
      DO bkk=kk-boxz,kk+boxz
       DO bii=ii-boxx,ii+boxx
         DO bjj=jj-boxy,jj+boxy
            IF( bkk >= 1 .AND. bkk <= ne .AND. bjj >= 1 .AND. bjj <= nr )THEN
              !Boundary condition in X
              ii_index=bii
              IF( bii < 1 )ii_index=bii+na
              IF( bii > na)ii_index=bii-na 
              IF( OPERATION == 'MEA2' .AND. bii== ii .AND. bjj==jj .AND. bkk==kk)CYCLE !We will not count the center of the box.

              IF( OPERATION /= 'COU2' )THEN !For COU2 operation count all the data even those flagged as undef.
                IF( datain(ii_index,bjj,bkk) /= UNDEF )THEN
                  NITEMS=NITEMS + 1
                  tmp_field(NITEMS)=datain(ii_index,bjj,bkk)
                ENDIF
              ELSE
                NITEMS=NITEMS+1   
                tmp_field(NITEMS)=datain(ii_index,bjj,bkk)
              ENDIF
            ENDIF
         ENDDO
       ENDDO
      ENDDO 
      !Perform the operation and save the result in dataout
      IF( OPERATION .EQ. 'MEAN' .OR. OPERATION .EQ. 'MEA2' )THEN
        IF (NITEMS > 0 )THEN
          dataout(ii,jj,kk)=sum( tmp_field(1:NITEMS) )/ REAL( NITEMS , r_size )
          !WRITE(*,*)NITEMS,tmp_field
        ENDIF

      ELSEIF( OPERATION == 'SIGM')THEN
        !Undef values won't be considered.
        tmp_mean=0.0d0
        tmp_var=0.0d0
        data_count=0
        DO iin=1,NITEMS
          IF( tmp_field(iin) .ne. UNDEF )THEN
            data_count=data_count+1
            tmp_mean=tmp_mean+tmp_field(iin)
            tmp_var=tmp_var+tmp_field(iin) ** 2
          ENDIF
        ENDDO
        IF( data_count .GT. 0)THEN
         tmp_mean=tmp_mean/REAL(data_count,r_size)
         tmp_var=tmp_var/REAL(data_count,r_size)
         dataout(ii,jj,kk)=SQRT(tmp_var - tmp_mean**2 )
        ENDIF

      ELSEIF( OPERATION == 'COUN' .OR. OPERATION == 'COU2' )THEN
        !Count values over a certain threshold (note that undef values will be 
        !always below the threshold.
        data_count=0
        IF(  NITEMS > 0 )THEN
          DO iin=1,NITEMS
             IF( tmp_field( iin ) >= threshold .AND. tmp_field( iin ) /= UNDEF )THEN
               data_count = data_count + 1
             ENDIF
          ENDDO
          dataout(ii,jj,kk)= REAL( data_count , r_size ) / REAL( NITEMS , r_size )
        ELSE
          dataout(ii,jj,kk) = 0.0d0

        ENDIF

      ELSEIF( OPERATION == 'MAXN')THEN
        IF( NITEMS > 0 )THEN
           !Local maximum tacking care of undef values.
           dataout(ii,jj,kk)=maxval( tmp_field(1:NITEMS) )
        ELSE
           dataout(ii,jj,kk)=UNDEF
        ENDIF 

      ELSEIF( OPERATION == 'MINN')THEN
        IF( NITEMS > 0 )THEN
          !Local maximum tacking care of undef values.
          dataout(ii,jj,kk)=minval( tmp_field(1:NITEMS) )
        ELSE
          dataout(ii,jj,kk)=UNDEF
        ENDIF

      ENDIF

  ENDDO
 ENDDO
ENDDO 
!$OMP END PARALLEL DO

RETURN
END SUBROUTINE BOX_FUNCTIONS_2D

SUBROUTINE ECHO_TOP_FAST(reflectivity,heigth,rrange,na,nr,ne,undef,nx,ny,nz,echo_top_3d,echo_top_2d)
!Curretnly this routine:
!Compute 3D echo top, echo base , echo depth , max dbz and max dbz z
!Performs interpolation from original radar grid (r,elevation) to an uniform (r,z) grid, where
!the parameters are computed.
!Then the result is interpolated back to the original radar grid.
IMPLICIT NONE

!TODO SOME OF THESE PARAMETERS SHOULD BE INPUT ARGUMENTS.
REAL(r_size) , PARAMETER :: MAX_Z_ECHO_TOP=20.0d4 , MAX_R_ECHO_TOP=240.0d03
INTEGER      , PARAMETER :: MAX_ECHO_TOP_LEVS=5
REAL(r_size) , PARAMETER :: DBZ_THRESHOLD_ECHO_TOP=5.0d0  !Echo top detection value.
INTEGER     ,INTENT(IN)  :: na,nr,ne
INTEGER     ,INTENT(IN)  :: nx,ny,nz
REAL(r_size),INTENT(IN)  :: reflectivity(na,nr,ne) , heigth(nr,ne) , rrange(nr,ne)  
REAL(r_size),INTENT(IN)  :: undef
REAL(r_size),INTENT(OUT) :: echo_top_3d(na,nr,ne)  , echo_top_2d(na,nr)
!Reflectivity and echo top in the vertical Z grid.
REAL(r_size)             :: tmp_ref(nr,ne) , tmp_echo_top_3d(nr,ne) , tmp_z(nr,ne)
REAL(r_size)             :: tmp_echo_top_2d(nr)
REAL(r_size)             :: tmp_bufr_3d(na,nr,ne) , tmp_bufr_2d(na,nr)
INTEGER                  :: index_up(nr,ne),index_up_inv(nr,ne)
REAL(r_size)             :: w_up(nr,ne) , w_up_inv(nr,ne)

REAL(r_size)             :: range_et(nr,ne) , z_et(nr,ne)

INTEGER                  :: i, ii , iii , kk , jj , ia , ip 

  !Define new grid in the R direction 
  DO kk=1,ne
     !We will interpolate reflectivity in all ranges to the location of the first elevation pixels.
     range_et(:,kk) = rrange(:,1)  
  ENDDO

  !Compute the index and weights for fast interpolation between the original grid and the modified range grid.
  !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(kk,ii,iii)
  DO kk = 1,ne
     !Initialize variables with zeros.
     index_up(:,kk) = 0
     index_up_inv(:,kk) = 0
     w_up(:,kk) = 0.0d0
     w_up_inv(:,kk) = 0.0d0

     DO ii = 1 , nr
        DO iii = 2 , nr
           !Compute indexes and weight for the forward interpolation from original grid to the vertical z coordinate grid.
           IF ( (rrange(ii,kk) >= range_et(iii-1,kk)) .and. (rrange(ii,kk) <= range_et(iii,kk) ) )THEN
              index_up(ii,kk)=iii
              w_up(ii,kk) = ( rrange(ii,kk) - range_et(iii-1,kk) ) / ( range_et(iii,kk) - range_et(iii-1,kk) )

           ENDIF  
           !Compute indexes and weigths for the inverse interpolation from the vertical z coordinate grid to the original grid.
           IF ( (range_et(ii,kk) .ge. rrange(iii-1,kk)) .and. (range_et(ii,kk) .le. rrange(iii,kk) ) )THEN
              index_up_inv(ii,kk)=iii
              w_up_inv(ii,kk) = ( range_et(ii,kk) - rrange(iii-1,kk) ) / ( rrange(iii,kk) - rrange(iii-1,kk) )

           ENDIF

        ENDDO

        !This is for extrapolation of reflectivity and echo top        
        IF( rrange(ii,kk) .gt. range_et(nr,kk) )THEN !We are beyond the range of the first level.
          index_up(ii,kk) = nr
          w_up(ii,kk) = 1.0d0
        ELSEIF( rrange(ii,kk) .lt. range_et(1,kk) )THEN !We are closer to the radar than the first gate of the first level
          index_up(ii,kk) = 2
          w_up(ii,kk) = 0.0d0
        ENDIF        

        IF( range_et(ii,kk) .gt. rrange(nr,kk) )THEN !We are beyond the range of the first level.
          index_up_inv(ii,kk) = nr
          w_up_inv(ii,kk) = 1.0d0
        ELSEIF( range_et(ii,kk) .lt. rrange(1,kk) )THEN !We are closer to the radar than the first gate of the first level
          index_up_inv(ii,kk) = 2
          w_up_inv(ii,kk) = 0.0d0
        ENDIF        

        !Height interpolation
        IF( index_up(ii,kk) .gt. 1 ) THEN
            tmp_z(ii,kk)=heigth(index_up(ii,kk),kk)*w_up(ii,kk) + heigth(index_up(ii,kk)-1,kk)*(1.0d0-w_up(ii,kk))
        ENDIF

     ENDDO
  ENDDO
  !$OMP END PARALLEL DO

  !Forward interpolation from the original grid to a vertical z coordinate grid.
  tmp_ref=undef

  !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(ii,jj,kk,tmp_ref,tmp_echo_top_3d,tmp_echo_top_2d)
  DO jj=1,na
    !Forward reflectivity interpolation
    DO ii=1,nr
      DO kk=1,ne
        IF( index_up(ii,kk) .gt. 1 ) THEN
          tmp_ref(ii,kk)=reflectivity(jj,index_up(ii,kk),kk)*w_up(ii,kk) + reflectivity(jj,index_up(ii,kk)-1,kk)*(1.0d0-w_up(ii,kk))
        ENDIF
      ENDDO
    ENDDO
    !Echo top computation
    DO ii=1,nr
      CALL ECHO_TOP_FAST_SUB(tmp_ref(ii,:),tmp_z(ii,:),ne,undef,tmp_echo_top_3d(ii,:), & 
                             tmp_echo_top_2d(ii),MAX_ECHO_TOP_LEVS,DBZ_THRESHOLD_ECHO_TOP)
      !DO kk=1,ne
      !   IF( tmp_echo_top_3d(ii,kk) <= 0.0 .and. tmp_echo_top_3d(ii,kk) /= undef )THEN
      !     WRITE(*,*) tmp_echo_top_3d(ii,kk) 
      !   ENDIF
      !ENDDO
    ENDDO

    !Inverse interpolation of echo top (3d and 2d)
    DO ii=1,nr
      DO kk=1,ne
        IF( tmp_echo_top_3d(index_up(ii,kk),kk) /= undef .and.    &
            tmp_echo_top_3d(index_up(ii,kk)-1,kk) /= undef ) THEN
          echo_top_3d(jj,ii,kk)=tmp_echo_top_3d(index_up(ii,kk),kk)*w_up(ii,kk) +   &
                                tmp_echo_top_3d(index_up(ii,kk)-1,kk)*(1.0d0-w_up(ii,kk))
        ELSEIF( tmp_echo_top_3d( index_up(ii,kk) , kk ) /= undef )THEN
          echo_top_3d(jj,ii,kk)= tmp_echo_top_3d( index_up(ii,kk) , kk )
        ELSEIF( tmp_echo_top_3d( index_up(ii,kk) - 1 , kk ) /= undef )THEN  
          echo_top_3d(jj,ii,kk)= tmp_echo_top_3d( index_up(ii,kk) - 1 , kk )     
        ELSE
          echo_top_3d(jj,ii,kk) = undef    
        ENDIF
      ENDDO
        IF( tmp_echo_top_2d(index_up(ii,1))   /= undef .and.   &                             
            tmp_echo_top_2d(index_up(ii,1)-1) /= undef ) THEN
          echo_top_2d(jj,ii)=tmp_echo_top_2d(index_up(ii,1))*w_up(ii,1) + & 
                      tmp_echo_top_2d(index_up(ii,1)-1)*(1.0d0-w_up(ii,1))
        ELSEIF( tmp_echo_top_2d( index_up(ii,1) ) /= undef )THEN
          echo_top_2d(jj,ii)= tmp_echo_top_2d( index_up(ii,1)  )
        ELSEIF( tmp_echo_top_2d( index_up(ii,1) - 1 ) /= undef )THEN
          echo_top_2d(jj,ii)= tmp_echo_top_2d( index_up(ii,1) - 1 )
        ELSE
          echo_top_2d(jj,ii) = undef
        ENDIF

    ENDDO

  ENDDO
  !$OMP END PARALLEL DO

  CALL BOX_FUNCTIONS_2D(echo_top_3d,na,nr,ne,undef,nx,ny,nz,'MEAN',0.0d0,tmp_bufr_3d)
  CALL BOX_FUNCTIONS_2D(echo_top_2d,na,nr,1 ,undef,nx,ny,0 ,'MEAN',0.0d0,tmp_bufr_2d)
  
  echo_top_3d = tmp_bufr_3d
  echo_top_2d = tmp_bufr_2d

RETURN
END SUBROUTINE ECHO_TOP_FAST

SUBROUTINE ECHO_TOP_FAST_SUB(reflectivity,z,nz,undef,echo_top_3d,echo_top_2d,max_levs,threshold)
!Vertical columns calculations
!Compute the possition of multiple echo tops in a single reflectivity column.
!Compute echo depth of each echo layer
!compute echo base
!compute max dbz
!This routine returns a vertical profile of echo base, echo top , echo depth , max dbz and max dbz a fore 
!each cloud layer.
!It also returns the vertical profile of the vertical gradient of the reflectivity field.

IMPLICIT NONE
INTEGER, INTENT(IN)     :: nz , max_levs
REAL(r_size),INTENT(IN) :: reflectivity(nz) , z(nz)
REAL(r_size),INTENT(OUT):: echo_top_3d(nz) !Vertical profile of echo top.
REAL(r_size),INTENT(OUT):: echo_top_2d     !Max echo top.
REAL(r_size),INTENT(IN) :: undef
REAL(r_size),INTENT(IN) :: threshold    !Reflectivity threshold to detect echo top.
INTEGER, PARAMETER      :: Nlevelstop=2
REAL(r_size)            :: tmp(max_levs,2) !echo_top , echo_base for each cloud layer.
REAL(r_size)            :: ref(nz) , ave_ref , sum_z
INTEGER                 :: jj, iz , base_count , top_count , tmp_count , itop , imax
LOGICAL                 :: base_detected , top_detected
LOGICAL                 :: found_first_maximum
REAL(r_size), PARAMETER :: first_maximum_threshold = 10.0d0 
REAL(r_size), PARAMETER :: refmin =0.0d0  ! Reflectivity value that will be assumed for UNDEF values in gradient computation.

echo_top_3d=undef
echo_top_2d=undef
tmp=UNDEF

base_count=0
top_count=0

ref=reflectivity   !reflectivity is intent in.

base_detected=.false.
top_detected=.false.

!Before computation extend data one or two levels below the first echo. This is done to prevent the first level to fall outside the computation
!of these scores.
DO iz=1,nz
   IF( z(iz) > 3000 )EXIT

   IF( ref(iz) /= UNDEF )THEN
      IF( iz>= 2)THEN
       ref(iz-1)=ref(iz)
      ENDIF

     EXIT
   ENDIF
ENDDO

DO iz=1,nz
   !Look for an echo base
   IF( ref(iz) > threshold .AND.  .NOT. base_detected .AND. ref(iz) /= UNDEF )THEN
       !An echo base has been detected.
       IF( base_count < max_levs)THEN
       base_detected=.true.
       top_detected=.false.
       base_count=base_count+1
       tmp(base_count,2)=z(iz)   !Echo base
       ENDIF
   ENDIF
   !Look for an echo top.
   IF( iz > Nlevelstop )THEN
     tmp_count=0
     DO jj=iz-Nlevelstop+1,iz
        IF( ref(jj) < threshold .OR. ref(jj) == UNDEF )tmp_count=tmp_count+1
     ENDDO
     IF( tmp_count == Nlevelstop .AND. .NOT. top_detected .AND. base_detected )THEN
     !An echo top has been detected
        top_detected=.true.
        base_detected=.false.
        IF( base_count <= max_levs )THEN
           tmp(base_count,1)=z(iz-Nlevelstop)  !Echo top
        ENDIF
     ENDIF
   ENDIF
   !Echo top associated with top of the radar domain.
   IF( iz == nz .AND. base_detected .AND. .NOT. top_detected )THEN
   !Domain is over but echo top has not been found! :( 
   !Force echo top
       IF( base_count <= max_levs )THEN
           tmp(base_count,1)=z(iz)  !Echo top
       ENDIF
   ENDIF


ENDDO !End for loop over levels

DO itop=1,max_levs
   IF( tmp(itop,1) .NE. UNDEF  .AND.  tmp(itop,2) .NE. UNDEF )THEN  !Echo top and echo base
       DO iz=1,nz !Juan - 1
          IF( z(iz) >= tmp(itop,2) .AND. z(iz) <= tmp(itop,1))THEN
               echo_top_3d(iz)=tmp(itop,1)
          ENDIF
          IF( z(iz) > tmp(itop,1) )EXIT
       ENDDO
   ENDIF
   !Find maximum echo top
   IF ( tmp(itop,1) .NE. UNDEF )THEN
      IF( echo_top_2d .EQ. UNDEF )THEN
        echo_top_2d=tmp(itop,1)
        ELSE
         IF( tmp(itop,1) >= echo_top_2d )THEN
           echo_top_2d = tmp(itop,1)
         ENDIF
      ENDIF
   ENDIF


ENDDO


RETURN
END SUBROUTINE ECHO_TOP_FAST_SUB



SUBROUTINE  ECHO_TOP(reflectivity,heigth,rrange,na,nr,ne,undef,nx,ny,nz,output_data_3d,output_data_2d)
!Curretnly this routine:
!Compute 3D echo top, echo base , echo depth , max dbz and max dbz z
!Performs interpolation from original radar grid (r,elevation) to an uniform (r,z) grid, where
!the parameters are computed.
!Then the result is interpolated back to the original radar grid.



IMPLICIT NONE

!TODO SOME OF THESE PARAMETERS SHOULD BE INPUT ARGUMENTS.
INTEGER      , PARAMETER :: NPAR_ECHO_TOP_3D=6 , NPAR_ECHO_TOP_2D=7 !Number of parameters in output arrays.
REAL(r_size) , PARAMETER :: MAX_Z_ECHO_TOP=20.0d4 , MAX_R_ECHO_TOP=240.0d03
REAL(r_size) , PARAMETER :: DZ_ECHO_TOP = 500.0d0 , DX_ECHO_TOP = 500.0d0
INTEGER      , PARAMETER :: MAX_ECHO_TOP_LEVS=5
REAL(r_size) , PARAMETER :: DBZ_THRESHOLD_ECHO_TOP=5.0d0  !Echo top detection value.


INTEGER     ,INTENT(IN)  :: na,nr,ne
INTEGER     ,INTENT(IN)  :: nx,ny,nz
REAL(r_size),INTENT(IN)  :: reflectivity(na,nr,ne) , heigth(nr,ne) , rrange(nr,ne)  
REAL(r_size),INTENT(IN)  :: undef
REAL(r_size),INTENT(OUT) :: output_data_3d(na,nr,ne,NPAR_ECHO_TOP_3D)  !Echo top , echo base , echo depth , max_dbz , maz_dbz_z , vertical_z_gradient
REAL(r_size),INTENT(OUT) :: output_data_2d(na,nr,NPAR_ECHO_TOP_2D)  !Max echo top, max_echo_base, max_echo_depth, col_max, height weighted col_max
!REAL(r_size)             :: tmp_output_data_3d(na,nr,ne,NPAR_ECHO_TOP_3D)
REAL(r_size)              :: tmp_radgrid_3d(na,nr,ne,NPAR_ECHO_TOP_3D)
REAL(r_size)              :: tmp_radgrid_2d(na,nr,NPAR_ECHO_TOP_2D)
REAL(r_size), ALLOCATABLE :: Z(:,:) , R(:,:) 
INTEGER, ALLOCATABLE      :: REGJ(:,:,:) , REGI(:,:,:) , INVI(:,:,:) , INVJ(:,:,:) 
INTEGER, ALLOCATABLE      :: NEARESTN(:,:) , INVNEARESTN(:,:)
REAL(r_size),ALLOCATABLE  :: W(:,:,:),INVW(:,:,:)
INTEGER                   :: REGNZ , REGNR
LOGICAL                   :: INITIALIZED=.FALSE.
!----->
REAL(r_size), ALLOCATABLE      :: REGREF(:,:)
CHARACTER(4)                   :: OPERATION='MEAN'
INTEGER                        :: i, ii , jj , ia , ip 
REAL(r_size),ALLOCATABLE       :: tmp_data_3d(:,:,:) ,  tmp_data_2d(:,:)


!WRITE(6,*)'HELLO FROM COMPUTE_ECHO_TOP'
!IF( .NOT. INITIALIZED) THEN
!Perform this part only in the first call.

REGNZ=INT(MAX_Z_ECHO_TOP / DZ_ECHO_TOP)+1
REGNR=INT(MAX_R_ECHO_TOP / DX_ECHO_TOP)+1
!WRITE(6,*)'REGNZ = ',REGNZ,' REGNR = ',REGNR


  ALLOCATE( Z(REGNR,REGNZ) , R(REGNR,REGNZ) )
  ALLOCATE( REGI(REGNR,REGNZ,4) , REGJ(REGNR,REGNZ,4), NEARESTN(REGNR,REGNZ) )
  ALLOCATE( INVI(nr,ne,4) , INVJ(nr,ne,4), INVNEARESTN(nr,ne) )
  ALLOCATE( W(REGNR,REGNZ,4),INVW(nr,ne,4) )

  !Set interpolation from range-elevation to range-z grid

  !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(ii,jj)
  DO ii=1,REGNR
   DO jj=1,REGNZ
     Z(ii,jj)=REAL(jj-1,r_size)*DZ_ECHO_TOP
     R(ii,jj)=REAL(ii-1,r_size)*DX_ECHO_TOP
     CALL com_xy2ij(nr,ne,rrange,heigth,R(ii,jj),Z(ii,jj),REGI(ii,jj,:),REGJ(ii,jj,:),W(ii,jj,:),NEARESTN(ii,jj),undef)
   ENDDO
  ENDDO
  !$OMP END PARALLEL DO

  !Set interpolation from range-z to range-elevation grid.

  !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(ii,jj)
  DO ii=1,nr
   DO jj=1,ne
     CALL com_xy2ij(REGNR,REGNZ,R,Z,rrange(ii,jj),heigth(ii,jj),INVI(ii,jj,:),INVJ(ii,jj,:),INVW(ii,jj,:),INVNEARESTN(ii,jj),undef) 
   ENDDO
  ENDDO
  !$OMP END PARALLEL DO

!ENDIF !End of first call only section.

ALLOCATE( tmp_data_3d(REGNR,REGNZ,NPAR_ECHO_TOP_3D))
ALLOCATE( tmp_data_2d(REGNR,NPAR_ECHO_TOP_2D))

ALLOCATE( REGREF(REGNR,REGNZ) )

tmp_radgrid_3d=UNDEF
tmp_radgrid_2d=UNDEF

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(ia,ii,jj,REGREF,tmp_data_3d,tmp_data_2d)
DO ia=1,na
 REGREF=UNDEF
 !Interp reflectivity from elevation-range grid to z-range grid. (nearest neighbor)
 DO ii=1,REGNR
  DO jj=1,REGNZ 
     IF( NEARESTN(ii,jj) > 0)THEN
        REGREF(ii,jj)=reflectivity(ia,REGI(ii,jj,NEARESTN(ii,jj)),REGJ(ii,jj,NEARESTN(ii,jj)))
     ENDIF
  ENDDO

  CALL ECHO_TOP_SUB(REGREF(ii,:),Z(ii,:),REGNZ,undef,tmp_data_3d(ii,:,:),tmp_data_2d(ii,:),MAX_ECHO_TOP_LEVS,DBZ_THRESHOLD_ECHO_TOP)
 ENDDO

 DO ii=1,nr
  DO jj=1,ne
    IF( INVNEARESTN(ii,jj) .GT. 0 )THEN
      tmp_radgrid_3d(ia,ii,jj,:)=tmp_data_3d(INVI(ii,jj,INVNEARESTN(ii,jj)),INVJ(ii,jj,INVNEARESTN(ii,jj)),:)
    ENDIF
    !IF( tmp_radgrid_3d(ia,ii,jj,5) /= UNDEF)THEN
    !  tmp_radgrid_3d(ia,ii,jj,5)=tmp_radgrid_3d(ia,ii,jj,5) !-topography(ia,ii,jj)  (dejo para mas adelante)
    !ENDIF
  ENDDO

  IF( INVNEARESTN(ii,1) .GT. 0 )THEN !We interpolate the data to the lowest level.
     tmp_radgrid_2d(ia,ii,:)=tmp_data_2d(INVI(ii,1,INVNEARESTN(ii,1)),:)
  ENDIF
 ENDDO

ENDDO
!$OMP END PARALLEL DO

!WRITE(*,*)maxval( tmp_radgrid_3d(:,:,:,1) )

!DO ip=1,NPAR_ECHO_TOP_3D
CALL BOX_FUNCTIONS_2D(tmp_radgrid_3d(:,:,:,1),na,nr,ne,undef,nx,ny,nz,'MEAN',0.0d0,output_data_3d(:,:,:,1))
CALL BOX_FUNCTIONS_2D(tmp_radgrid_3d(:,:,:,2),na,nr,ne,undef,nx,ny,nz,'MEAN',0.0d0,output_data_3d(:,:,:,2))
CALL BOX_FUNCTIONS_2D(tmp_radgrid_3d(:,:,:,3),na,nr,ne,undef,nx,ny,nz,'MEAN',0.0d0,output_data_3d(:,:,:,3))
CALL BOX_FUNCTIONS_2D(tmp_radgrid_3d(:,:,:,4),na,nr,ne,undef,nx,ny,nz,'MEAN',0.0d0,output_data_3d(:,:,:,4))
CALL BOX_FUNCTIONS_2D(tmp_radgrid_3d(:,:,:,5),na,nr,ne,undef,0,1,0,'MINN',0.0d0,output_data_3d(:,:,:,5))
CALL BOX_FUNCTIONS_2D(tmp_radgrid_3d(:,:,:,6),na,nr,ne,undef,0,1,0,'MINN',0.0d0,output_data_3d(:,:,:,6))
!END DO

DO ip=1,NPAR_ECHO_TOP_2D
CALL BOX_FUNCTIONS_2D(tmp_radgrid_2d(:,:,ip),na,nr,1,undef,nx,ny,0,'MEAN',0.0d0,output_data_2d(:,:,ip))
END DO


!WRITE(*,*)maxval( output_data_3d(:,:,:,1) )

DEALLOCATE(  tmp_data_3d , tmp_data_2d )
DEALLOCATE( REGREF )

DEALLOCATE( Z, R )
DEALLOCATE( REGI , REGJ, NEARESTN )
DEALLOCATE( INVI , INVJ, INVNEARESTN )
DEALLOCATE( W , INVW )

!INITIALIZED=.TRUE.

RETURN
END SUBROUTINE ECHO_TOP

SUBROUTINE ECHO_TOP_SUB(reflectivity,z,nz,undef,output_3d,output_2d,max_levs,threshold)
!Vertical columns calculations
!Compute the possition of multiple echo tops in a single reflectivity column.
!Compute echo depth of each echo layer
!compute echo base
!compute max dbz
!This routine returns a vertical profile of echo base, echo top , echo depth , max dbz and max dbz a fore 
!each cloud layer.
!It also returns the vertical profile of the vertical gradient of the reflectivity field.


IMPLICIT NONE
INTEGER      , PARAMETER :: NPAR_ECHO_TOP_3D=6 , NPAR_ECHO_TOP_2D=7 !Number of parameters in output arrays.
INTEGER, INTENT(IN)  :: nz , max_levs
REAL(r_size),INTENT(IN) :: reflectivity(nz) , z(nz)
REAL(r_size),INTENT(OUT):: output_3d(nz,NPAR_ECHO_TOP_3D) !echo_top, echo_base, echo_depth , max_dbz , max_dbz_z , reflectivity gradient
REAL(r_size),INTENT(OUT):: output_2d(NPAR_ECHO_TOP_2D) !Max echo top, max_echo_base, max_echo_depth, col_max, height weighted col_max , first ref maximum height , intensity.
REAL(r_size),INTENT(IN) :: undef
REAL(r_size),INTENT(IN) :: threshold    !Reflectivity threshold to detect echo top.
INTEGER, PARAMETER      :: Nlevelstop=2
REAL(r_size)            :: tmp(max_levs,5) !echo_top, echo_base, echo_depth , max_dbz , max_dbz_z
REAL(r_size)            :: ref(nz) , ave_ref , sum_z
INTEGER                 :: jj, iz , base_count , top_count , tmp_count , itop , imax
LOGICAL                 :: base_detected , top_detected
LOGICAL                 :: found_first_maximum
REAL(r_size), PARAMETER :: first_maximum_threshold = 10.0d0 
INTEGER     , PARAMETER :: NDELTAZ=5      ! NDELTAZ * dz is the distance used to estimate vertical reflectivity gradient
REAL(r_size), PARAMETER :: refmin =0.0d0  ! Reflectivity value that will be assumed for UNDEF values in gradient computation.

output_3d=undef
output_2d=undef
tmp=UNDEF

base_count=0
top_count=0

ref=reflectivity   !reflectivity is intent in.

base_detected=.false.
top_detected=.false.

!Before computation extend data one or to levels below the first echo. This is done to prevent the first level to fall outside the computation
!of these scores.
DO iz=1,nz
   IF( z(iz) > 3000 )EXIT

   IF( ref(iz) /= UNDEF )THEN
      IF( iz>= 2)THEN
       ref(iz-1)=ref(iz)
      ENDIF

     EXIT
   ENDIF
ENDDO

DO iz=1,nz
   !Look for an echo base
   IF( ref(iz) > threshold .AND.  .NOT. base_detected .AND. ref(iz) /= UNDEF )THEN
       !An echo base has been detected.
       IF( base_count < max_levs)THEN
       base_detected=.true.
       top_detected=.false.
       base_count=base_count+1
       tmp(base_count,2)=z(iz)   !Echo base
       tmp(base_count,4)=ref(iz) !Max dbz
       tmp(base_count,5)=z(iz)   !Max dbz_z
       ENDIF
   ENDIF
   !Look for an echo top.
   IF( iz > Nlevelstop )THEN
     tmp_count=0
     DO jj=iz-Nlevelstop+1,iz
        IF( ref(jj) < threshold .OR. ref(jj) == UNDEF )tmp_count=tmp_count+1
     ENDDO
     IF( tmp_count == Nlevelstop .AND. .NOT. top_detected .AND. base_detected )THEN
     !An echo top has been detected
        top_detected=.true.
        base_detected=.false.
        IF( base_count <= max_levs )THEN
           tmp(base_count,1)=z(iz-Nlevelstop)  !Echo top
        ENDIF
     ENDIF
   ENDIF
   !Echo top associated with top of the radar domain.
   IF( iz == nz .AND. base_detected .AND. .NOT. top_detected )THEN
   !Domain is over but echo top has not been found! :( 
   !Force echo top
       IF( base_count <= max_levs )THEN
           tmp(base_count,1)=z(iz)  !Echo top
       ENDIF
   ENDIF
   !Compute max dbz
   IF( base_detected .AND. .NOT. top_detected .AND. ref(iz) /=UNDEF )THEN
       !We are within a cloud or an echo region. Compute max dbz.
       IF( ref(iz) > tmp(base_count,4) )THEN  !Max dbz
           tmp(base_count,4)=ref(iz)  !Max dbz
           tmp(base_count,5)=z(iz)    !Max dbz z
       ENDIF
   ENDIF
   !Compute vertical gradient of reflectivity.
   IF( iz <= nz-NDELTAZ)THEN
   IF( ref(iz) /= UNDEF )THEN
    IF(  ref( iz + NDELTAZ ) /= UNDEF )THEN
        output_3d(iz,6)= ( ref(iz+NDELTAZ) - ref(iz) ) /( z(iz+NDELTAZ) - z(iz) )
    ELSE
        output_3d(iz,6)= ( refmin          - ref(iz) ) /( z(iz+NDELTAZ) - z(iz) )
    ENDIF
   ENDIF
   ENDIF


ENDDO !End for loop over levels

DO itop=1,max_levs
   IF( tmp(itop,1) .NE. UNDEF  .AND.  tmp(itop,2) .NE. UNDEF )THEN  !Echo top and echo base
       DO iz=1,nz-1
          IF( z(iz) >= tmp(itop,2) .AND. z(iz) <= tmp(itop,1))THEN
               output_3d(iz,1:2)=tmp(itop,1:2)
               output_3d(iz,3)  =tmp(itop,1)-tmp(itop,2)
               output_3d(iz,4:5)=tmp(itop,4:5)
          ENDIF
          IF( z(iz) > tmp(itop,1) )EXIT
       ENDDO
   ENDIF
   !Find maximum echo top
   IF ( tmp(itop,1) .NE. UNDEF )THEN
      IF( output_2d(1) .EQ. UNDEF )THEN
        output_2d(1)=tmp(itop,1)
        ELSE
         IF( tmp(itop,1) >= output_2d(1) )THEN
           output_2d(1) = tmp(itop,1)
         ENDIF
      ENDIF
   ENDIF

   !Find maximum echo base
   IF ( tmp(itop,2) .NE. UNDEF )THEN
      IF( output_2d(2) .EQ. UNDEF )THEN
        output_2d(2)=tmp(itop,2)
        ELSE
         IF( tmp(itop,2) >= output_2d(2) )THEN
           output_2d(2) = tmp(itop,2)
         ENDIF
      ENDIF
   ENDIF

   !Find maximum echo depth
   IF ( tmp(itop,3) .NE. UNDEF )THEN
      IF( output_2d(3) .EQ. UNDEF )THEN
        output_2d(3)=tmp(itop,3)
        ELSE
         IF( tmp(itop,3) >= output_2d(3) )THEN
           output_2d(3) = tmp(itop,3)
         ENDIF
      ENDIF
   ENDIF

   !Find maximum reflectivity (colmax)
   IF ( tmp(itop,4) .NE. UNDEF )THEN
      IF( output_2d(4) .EQ. UNDEF )THEN
        output_2d(4)=tmp(itop,4)
        ELSE
         IF( tmp(itop,4) >= output_2d(4) )THEN
           output_2d(4) = tmp(itop,4)
         ENDIF
      ENDIF
   ENDIF

   IF ( tmp(itop,4) .NE. UNDEF )THEN
      IF( output_2d(4) .EQ. UNDEF )THEN
        output_2d(4)=tmp(itop,4)
        ELSE
         IF( tmp(itop,4) >= output_2d(4) )THEN
           output_2d(4) = tmp(itop,4)
         ENDIF
      ENDIF
   ENDIF

ENDDO


!Compute heigh weigthed averaged reflectivity, the height of the first reflectivity maximum and its intensity.


ave_ref=0
sum_z=0
found_first_maximum=.FALSE.
imax=0
DO iz=1,nz

 IF( reflectivity(iz) .NE. UNDEF )THEN
   ave_ref = z(iz) * reflectivity(iz)
   sum_z   = z(iz)
 ENDIF

 IF( reflectivity(iz) /= UNDEF .AND. output_2d(7) == UNDEF )THEN
   output_2d(7) = reflectivity(iz)    !Intensity of first reflectivity maximun
   output_2d(6) = z(iz)               !Height of first reflectivity maximum
 ENDIF
 IF( reflectivity(iz) /= UNDEF .AND. (reflectivity(iz) - output_2d(7)) <  & 
     first_maximum_threshold .AND. .NOT. found_first_maximum )THEN
   found_first_maximum=.TRUE. 
 ELSE
   IF( reflectivity(iz) > output_2d(7) )THEN
     output_2d(7) = reflectivity(iz) !Keep updating the maximum until we reach the first maximum.
     output_2d(6) = z(iz)            !Keep updating the height of the maximum 
   ENDIF
 ENDIF

ENDDO
IF( sum_z .GT. 0 )THEN
  output_2d(5) = ave_ref / sum_z
ELSE
  output_2d(5) = 0
ENDIF

!Max echo top, max_echo_base, max_echo_depth, col_max, height weighted col_max
RETURN
END SUBROUTINE ECHO_TOP_SUB

!-------------------------------------------------------------------------------------------
!
! Computation of blocking by the topography
!
!-------------------------------------------------------------------------------------------


SUBROUTINE COMPUTE_BLOCKING( radarz , topo , na , nr , ne , undef ,  & 
           &                 radar_beam_width_v , beam_length , radarrange , radarelev , blocking )

INTEGER     , INTENT(IN)  :: na,nr,ne
REAL(r_size), INTENT(IN)  :: radarz(na,nr,ne)
REAL(r_size), INTENT(IN)  :: topo(na,nr,ne)
REAL(r_size), INTENT(OUT) :: blocking(na,nr,ne)
REAL(r_size), INTENT(IN)  :: undef
REAL(r_size), INTENT(IN)  :: radar_beam_width_v , beam_length , radarrange(nr) , radarelev(ne)

REAL(r_size) :: alfa , beta , diag , effective_beam_width(nr,ne) , effective_height_factor(nr,ne)
REAL(r_size) :: max_blocking_factor
REAL(r_size) :: norm_h , min_norm_h 
INTEGER      :: ii , jj , kk
REAL(r_size) :: lthreshold

effective_height_factor = 0.0d0
effective_beam_width    = 0.0d0

DO kk=1,ne
 DO jj=1,nr
   !!!!alfa=atan(beam_length/vert_beam_width)
   !!!!diag =sqrt( beam_length**2 + vert_beam_width**2 )
   !!!!beta=alfa-radarelev(kk)*deg2rad
   !!!!max_vertical_extent(jj,kk)=diag*cos(beta)

   !We first compute the effective_height_factor which gives use the height of the lowest part of the beam.
   !taking into account the elevation angle (and assuming that the topography value is uniform over this volume.
   effective_height_factor(jj,kk) = -0.5 * beam_length * sin( radarelev(kk) * deg2rad  ) 

   !Then we compute the effective beam width. Tacking into account that the elevation angle
   effective_beam_width(jj,kk) =  ( radar_beam_width_v * radarrange(jj)*(deg2rad)/2.0d0 ) * cos( radarelev(kk) * deg2rad )

   !This procedure allow us to compute the height of the lowest volume vertex.
 ENDDO
ENDDO

!DO kk=1,ne
!  WRITE(*,*)'Effective height factor',kk,effective_height_factor(1,kk),effective_height_factor(nr,kk)
!  WRITE(*,*)'Effective beam width',kk,effective_beam_width(1,kk),effective_beam_width(nr,kk) 
!ENDDO

blocking=0.0d0
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(ii,jj,kk,norm_h,min_norm_h,max_blocking_factor)
DO ii=1,na
  DO kk=1,ne
    min_norm_h=1.0d0
    DO jj=1,nr
       !Compute heigth over the terrain normalized by the beam vertical extent.
       norm_h=( radarz(ii,jj,kk) + effective_height_factor(jj,kk) -topo(ii,jj,kk) ) / effective_beam_width(jj,kk)
       IF( norm_h < min_norm_h )THEN
          min_norm_h=norm_h
          IF( min_norm_h < 1.0d0  )THEN
            !We have some blocking, lets compute the blocking magnitude
            IF( min_norm_h > -1.0d0 )THEN
               blocking(ii,jj:nr,kk)=( min_norm_h * SQRT( 1.0d0 - min_norm_h ** 2 ) - ASIN( min_norm_h ) + pi/2.0d0 )/pi
            ELSE
               blocking(ii,jj:nr,kk)=1.0d0
               !This beam is totally blocked so we exit the loop.
               EXIT
            ENDIF
          ENDIF
       ENDIF
    ENDDO
  ENDDO
ENDDO
!$OMP END PARALLEL DO

!DO kk=1,ne
!   WRITE(*,*)'Blocking fortran',MINVAL(blocking(:,:,kk)),MAXVAL(blocking(:,:,kk))
!   WRITE(*,*)'Radar z fortran',MINVAL(radarz(:,:,kk)),MAXVAL(radarz(:,:,kk))
!   WRITE(*,*)'Topo fortran',MINVAL(topo(:,:,kk)),MAXVAL(topo(:,:,kk))
!ENDDO



RETURN
END SUBROUTINE COMPUTE_BLOCKING

SUBROUTINE SOBEL_FILTER( field , nx , ny , nz , undef , nboxx , nboxy , nboxz , edge_tr , edge_intensity , edge_mask )
!This subroutine detects edges in one field computing differences between nearby grid points.
!When these differences are over a certain threshold (nh_thresh) an edge is detected.
!The output is a logical mask (edge_mask) which is true where the routine detected an edge and false 
!otherwise.

IMPLICIT NONE
INTEGER       , INTENT(IN)       :: nx,ny,nz                  !Matrices size
INTEGER       , INTENT(IN)       :: nboxx,nboxy,nboxz         !expansion factors edges will be expanded nx,ny,nz grid points.
REAL(r_size)  , INTENT(IN)       :: field(nx,ny,nz)           !Original field
REAL(r_size)  , INTENT(IN)       :: edge_tr , undef 

REAL(r_size)  , INTENT(OUT)      :: edge_intensity(nx,ny,nz)
LOGICAL       , INTENT(OUT)      :: edge_mask(nx,ny,nz) 

REAL(r_size)                     :: sx , sy
INTEGER                          :: ix,iy,iz , maxx,maxy,maxz , minx,miny,minz , boxx,boxy,boxz


edge_mask=.False.

edge_intensity=0.0d0

boxx=max( nboxx , 0 )
boxy=max( nboxy , 0 )
boxz=max( nboxz , 0 )

!Loop over elevations

DO ix=2,nx-1
  !Loop over azimuths
  DO iy=2,ny-1
    !Loop over ranges
    DO iz=2,nz-1

        IF( field(ix,iy,iz) .ne. undef )THEN

        !TODO!!! CONSIDERAR MEJOR LOS UNDEF PORQUE NO ESTAN CONSIDERADOS EN ESTA CUENTA.

           !Compute the SOBEL operator.
           sx   =(field(ix-1,iy+1,iz) - field(ix+1,iy+1,iz) ) +  &
         & 2.0d0*(field(ix-1,iy  ,iz) - field(ix+1,iy  ,iz) ) +  &
         &       (field(ix-1,iy-1,iz) - field(ix+1,iy-1,iz) ) / 6.0

            sy  =(field(ix+1,iy-1,iz) - field(ix+1,iy+1,iz) ) +  &
         & 2.0d0*(field(ix  ,iy-1,iz) - field(ix  ,iy+1,iz) ) +  &
         &       (field(ix-1,iy-1,iz) - field(ix-1,iy+1,iz) ) / 6.0

           edge_intensity(ix,iy,iz)= sqrt(sx**2 + sy**2)

          
           IF(   edge_intensity(ix,iy,iz) > edge_tr  ) THEN  

             !Expand the grid points identified as edges using nx,ny,nz

             maxx=min( ix+boxx ,  nx )
             maxy=min( iy+boxy ,  ny )
             maxz=min( iz+boxz ,  nz )
             minx=max( ix-boxx ,  1  )
             miny=max( iy-boxy ,  1  )
             minz=max( iz-boxz ,  1  )

             edge_mask(minx:maxx,miny:maxy,minz:maxz)=.True.         
 

           ENDIF 
 
         ENDIF

     ENDDO  !End loop over ranges 
   ENDDO    !End loop over elevations
ENDDO       !End loop over azimuth


RETURN
END SUBROUTINE SOBEL_FILTER

SUBROUTINE DOPPLER_EDGE_FILTER( vdiff , v , nx , ny , nz , undef , nboxx , nboxy , nboxz , edge_tr , edge_mask )
!This subroutine detects edges in one field computing differences between nearby grid points.
!When these differences are over a certain threshold (nh_thresh) an edge is detected.
!The output is a logical mask (edge_mask) which is true where the routine detected an edge and false 
!otherwise.

IMPLICIT NONE
INTEGER       , INTENT(IN)       :: nx,ny,nz                       !Matrices size
INTEGER       , INTENT(IN)       :: nboxx,nboxy,nboxz              !expansion factors edges will be expanded nx,ny,nz grid points.
REAL(r_size)  , INTENT(IN)       :: vdiff(nx,ny,nz) , v(nx,ny,nz)  !Dealiased-aliased difference , corrected v
REAL(r_size)  , INTENT(IN)       :: edge_tr , undef

INTEGER       , INTENT(OUT)      :: edge_mask(nx,ny,nz)

INTEGER                          :: ix,iy,iz , maxx,maxy,maxz , minx,miny,minz , boxx,boxy,boxz , iix,iiy,iiz

edge_mask=0

boxx=max( nboxx , 0 )
boxy=max( nboxy , 0 )
boxz=max( nboxz , 0 )

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(maxx,maxy,maxz,minx,miny,minz,ix,iy,iz,iix,iiy,iiz)

DO ix=1,nx
  DO iy=1,ny
    DO iz=1,nz

        IF ( v(ix,iy,iz) .ne. undef .and. vdiff(ix,iy,iz) /= 0.0d0 ) THEN
               
           !Look for neighbors.
           maxx=min( ix+nboxx ,  nx )
           maxy=min( iy+nboxy ,  ny )
           maxz=min( iz+nboxz ,  nz )
           minx=max( ix-nboxx ,  1  )
           miny=max( iy-nboxy ,  1  )
           minz=max( iz-nboxz ,  1  )

           DO iix=minx,maxx
             DO iiy=miny,maxy
               DO iiz=minz,maxz
                 IF( vdiff(iix,iiy,iiz) == 0.0d0 .and. abs( v(ix,iy,iz) - v(iix,iiy,iiz) ) > edge_tr )THEN

                     edge_mask(iix,iiy,iiz) = 1
 
                 ENDIF 
               ENDDO
             ENDDO      
           ENDDO

        ENDIF

    ENDDO
  ENDDO
ENDDO

!$OMP END PARALLEL DO

RETURN

END SUBROUTINE DOPPLER_EDGE_FILTER

SUBROUTINE MULTIPLE_1D_INTERPOLATION( field , nx, ny , nz , undef , xx , yy , nxx , fieldo )
!Interpolate the values given in the 3-dimensional array to the 1d function defined by xx and yy.
!Linear interpolation is used and values ouside xx range are asigned to the max/min value of yy.
IMPLICIT NONE
INTEGER      , INTENT(IN)      :: nx , ny , nz , nxx 
REAL(r_size) , INTENT(IN)      :: field(nx,ny,nz)            !Input values to be interpolated.
REAL(r_size) , INTENT(IN)      :: undef            
REAL(r_size) , INTENT(IN)      :: xx(nxx) , yy(nxx)          !x and y that define the interpolation function.
REAL(r_size) , INTENT(OUT)     :: fieldo(nx,ny,nz)           !Field containing interpolated values.

INTEGER                        :: ix , iy , iz , iix

fieldo=0.0d0

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(ix,iy,iz,iix)
DO ix = 1 , nx
  DO iy = 1 , ny
    DO iz = 1 , nz
       IF( field(ix,iy,iz) == undef )THEN
         fieldo(ix,iy,iz) = undef 
         cycle
       ENDIF
       IF( field(ix,iy,iz) <= xx(1) )THEN
         fieldo(ix,iy,iz) = yy(1) 
       ELSEIF( field(ix,iy,iz) >= xx(nxx) )THEN
         fieldo(ix,iy,iz) = yy(nxx)
       ELSE
         DO iix = 1 , nxx-1
            IF( field(ix,iy,iz) >= xx(iix) .and. field(ix,iy,iz) < xx(iix+1) )THEN
                fieldo(ix,iy,iz)= yy(iix) + ( field(ix,iy,iz)-xx(iix) ) * ( ( yy(iix+1) - yy(iix) ) / ( xx(iix+1) - xx(iix) ) )
            ENDIF
         ENDDO
       ENDIF
    ENDDO
  ENDDO
ENDDO
!$OMP END PARALLEL DO

RETURN

END SUBROUTINE MULTIPLE_1D_INTERPOLATION

SUBROUTINE COMPUTE_DISTANCE(var,var2,nx,ny,nz,undef,nx_box,ny_box,nz_box,dist)

IMPLICIT NONE
INTEGER,INTENT(IN)       :: nx , ny , nz
INTEGER,INTENT(IN)       :: nx_box , ny_box , nz_box
REAL(r_size),INTENT(IN)  :: var(nx,ny,nz) 
REAL(r_size),INTENT(IN)  :: var2(nx,ny,nz)
REAL(r_size),INTENT(IN)  :: undef
REAL(r_size),INTENT(OUT) :: dist(nx,ny,nz)
INTEGER                  :: box_size , ii , jj , kk , ndata
INTEGER                  :: imin , imax , jmin , jmax , kmin , kmax
INTEGER                  :: i,j,k

dist=0.0d0

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(ii,jj,kk,imin,imax,jmin,jmax,kmin,kmax,ndata,i,j,k)
DO ii=1,nx
  DO jj=1,ny
    DO kk=1,nz

       imax=min(i+nx_box,nx)
       imin=max(i-nx_box,1 )

       jmax=min(j+ny_box,ny)
       jmin=max(j-ny_box,1 )

       kmax=min(k+nz_box,nz)
       kmin=max(k-nz_box,1 )

       ndata=0

       DO i=imin,imax
        DO j=jmin,jmax
         DO k=kmin,kmax
           IF( var2(i,j,k) /= undef )THEN
              ndata=ndata+1
              dist(ii,jj,kk) = dist(ii,jj,kk) + var2(i,j,k)
           ENDIF
         ENDDO
       ENDDO
     ENDDO

     !Compute the difference between the mean of var2 in a local box
     !and the value of var at the center point of the local box.
     IF( ndata >= 1 )THEN
        dist(ii,jj,kk)= abs( dist(ii,jj,kk)/ndata - var(i,j,k) )
     ELSE
        dist(ii,jj,kk)= undef
     ENDIF

    ENDDO
  ENDDO
ENDDO
!$OMP END PARALLEL DO
END SUBROUTINE COMPUTE_DISTANCE


SUBROUTINE DETECT_MISSING(ref , na , nr , ne , undef , min_ref , threshold , nmissing_max , missing_mask )
IMPLICIT NONE
INTEGER , INTENT(IN)       :: na , nr , ne , nmissing_max
REAL(r_size) , INTENT(IN)  :: ref(na,nr,ne) , undef , min_ref , threshold
LOGICAL , INTENT(OUT)      :: missing_mask(na,nr,ne)

INTEGER                    :: ia , ir , ie , nmissing
LOGICAL                    :: ismissing

missing_mask=.False.

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(ia,ie,ir,ismissing,nmissing)
DO ia = 1 , na
  DO ie = 1 , ne
     ismissing=.False.
     nmissing=1;
   
     DO ir = 2 , nr
         if( abs( ref(ia,ir-1,ie) - ref(ia,ir,ie) ) > threshold .and. ref(ia,ir,ie) <= min_ref )then
             ismissing=.True.
             nmissing=1
         endif
         if( ismissing .and. ref(ia,ir,ie) > min_ref )then
             ismissing=.False.
         endif
         if( nmissing > nmissing_max )then
             ismissing=.False.
         endif
         if( ismissing )then
            missing_mask(ia,ir,ie)=.True.
            nmissing=nmissing+1
         endif
    ENDDO
  ENDDO
ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE DETECT_MISSING

!-----------------------------------------------------------------------
! (X,Y) --> (i,j) conversion (General pourpuse interpolation)
!   [ORIGINAL AUTHOR:] Masaru Kunii
!-----------------------------------------------------------------------
SUBROUTINE com_xy2ij(nx,ny,fx,fy,datax,datay,dist_min_x,dist_min_y,ratio,nearestn,undef)


  IMPLICIT NONE
  ! --- inout variables
  INTEGER,INTENT(IN) :: nx,ny !number of grid points
  REAL(r_size),INTENT(IN) :: fx(nx,ny),fy(nx,ny) !(x,y) at (i,j)
  REAL(r_size),INTENT(IN) :: datax,datay !target (lon,lat)
  REAL(r_size),INTENT(IN) :: undef
  ! --- local work variables
  LOGICAL,PARAMETER :: detailout = .FALSE.
  INTEGER,PARAMETER :: num_grid_ave = 4  ! fix
  INTEGER :: ix,jy,ip,wk_maxp
  INTEGER :: iorder_we,iorder_sn
  INTEGER :: nxp,nyp
  INTEGER,PARAMETER :: order = 2
  REAL(r_size),PARAMETER :: max_dist = 2.0e+6
  REAL(r_size) :: rxmax, rxmin, rymax, rymin   
  REAL(r_size) :: dist(num_grid_ave)  , tmp_dist(num_grid_ave)
  
  INTEGER,INTENT(OUT) :: dist_min_x( num_grid_ave)
  INTEGER,INTENT(OUT) :: dist_min_y( num_grid_ave) 
  INTEGER,INTENT(OUT) :: nearestn(1)
  REAL(r_size) :: wk_dist, sum_dist
  REAL(r_size),INTENT(OUT) :: ratio(num_grid_ave)
  
  IF(detailout) THEN
    WRITE(6,'(A)') '====================================================='
    WRITE(6,'(A)') '      Detailed output of SUBROUTINE com_pos2ij       '
    WRITE(6,'(A)') '====================================================='    
  END IF
  ! ================================================================
  !   Check the Order of fx,fy 
  ! ================================================================   
  iorder_we = 1
  iorder_sn = 1
  IF(fx(1,1) > fx(2,1)) THEN
    iorder_we = -1
  END IF
  IF(fy(1,1) > fy(1,2)) THEN
    iorder_sn = -1
  END IF
  IF(detailout) THEN  
    WRITE(6,'(3X,A,I5)') 'X Order (WE) :',iorder_we 
    WRITE(6,'(3X,A,I5)') 'Y Order (SN) :',iorder_sn 

  END IF
   
  ratio=UNDEF
  dist_min_x=0
  dist_min_y=0
  nearestn=0
    ! ================================================================
    !   Nearest 4 Grid Points Interpolation
    ! ================================================================   
      ! ------------------------------------------------------------
      !    Search 4-Grid Points
      ! ------------------------------------------------------------      
      dist(1:num_grid_ave) = 1.D+10
      DO jy=1,ny-1
        DO ix=1,nx-1
          rxmax = MAXVAL(fx(ix:ix+1, jy:jy+1))
          rxmin = MINVAL(fx(ix:ix+1, jy:jy+1))
          rymax = MAXVAL(fy(ix:ix+1, jy:jy+1))
          rymin = MINVAL(fy(ix:ix+1, jy:jy+1))
         IF(rxmin <= datax .AND. rxmax >= datax .AND. &
           & rymin <= datay .AND. rymax >= datay ) THEN
          tmp_dist(1)=( fx(ix,jy) - datax )** order + ( fy(ix,jy) - datay )** order
          tmp_dist(2)=( fx(ix+1,jy) - datax )** order + ( fy(ix+1,jy) - datay )** order
          tmp_dist(3)=( fx(ix+1,jy+1) - datax )** order + ( fy(ix+1,jy+1) - datay )** order
          tmp_dist(4)=( fx(ix,jy+1) - datax )** order + ( fy(ix,jy+1) - datay )** order

   
          IF( maxval(tmp_dist) <= maxval(dist) )THEN
            nearestn=minloc(tmp_dist)
            dist=tmp_dist
            dist_min_x(1)=ix
            dist_min_x(2)=ix+1
            dist_min_x(3)=ix+1
            dist_min_x(4)=ix
            dist_min_y(1)=jy
            dist_min_y(2)=jy
            dist_min_y(3)=jy+1
            dist_min_y(4)=jy+1
          ENDIF 
         ENDIF

        END DO
      END DO

      IF( dist_min_x(1) > 0)THEN
      sum_dist = dist(1) + dist(2) + dist(3) + dist(4)
      ratio(1) = dist(1)/sum_dist
      ratio(2) = dist(2)/sum_dist
      ratio(3) = dist(3)/sum_dist
      ratio(4) = dist(4)/sum_dist
      ENDIF
      !IF(detailout) WRITE(6,'(2X,A,5F15.5)') 'ratio      :',ratio(1:4),SUM(ratio(1:4))

        

  RETURN
END SUBROUTINE com_xy2ij



END MODULE QC
