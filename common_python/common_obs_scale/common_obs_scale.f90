MODULE common_obs_scale
!=======================================================================
!
! [PURPOSE:] Observational procedures
!
! [HISTORY:]
!   01/23/2009 Takemasa MIYOSHI  created
!   10/04/2012 Guo-Yuan Lien     modified for GFS model
!   07/25/2014 Guo-Yuan Lien     modified for SCALE model
!   .......... See git history for the following revisions
!   05/02/2018 Modified to be used from python scripts.
!
!=======================================================================
!
! [LETKF observation format]
!   (In files, all variables are stored in single-precision float)
!
!  column  description
!     (1)  variable type (1..nid_obs; see 'id_*_obs' parameters)
!     (2)  longitude (degree)
!     (3)  latitude (degree)
!     (4)  level/height
!            u,v,t,tv,q,rh: level (hPa)
!            ps: station elevation (m)
!     (5)  observation value
!            wind (m/s)
!            temperature (K)
!            specific humidity (kg/kg)
!            relative humidity (%)
!            surface pressure (hPa)
!     (6)  observation error
!            unit same as observation value
!     (7)  observation platform type (1..nobtype+1; see 'obtypelist' array)
!     (8)  observation time relative to analysis time (sec)
!
!=======================================================================
!$USE OMP_LIB

  IMPLICIT NONE
  PUBLIC

  INTEGER,PARAMETER :: r_dble=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
  INTEGER,PARAMETER :: r_size=r_dble


  INTEGER,PARAMETER :: nid_obs_varlocal=9 !H08
!
! conventional observations
!
  INTEGER,PARAMETER :: id_u_obs=2819
  INTEGER,PARAMETER :: id_v_obs=2820
  INTEGER,PARAMETER :: id_t_obs=3073
  INTEGER,PARAMETER :: id_tv_obs=3074
  INTEGER,PARAMETER :: id_q_obs=3330
  INTEGER,PARAMETER :: id_rh_obs=3331
!
! surface observations codes > 9999
!
  INTEGER,PARAMETER :: id_ps_obs=14593
  INTEGER,PARAMETER :: id_rain_obs=19999
  INTEGER,PARAMETER :: id_tclon_obs=99991  ! TC vital
  INTEGER,PARAMETER :: id_tclat_obs=99992  ! TC vital
  INTEGER,PARAMETER :: id_tcmip_obs=99993  ! TC vital
!
! radar observations
!
  INTEGER,PARAMETER :: id_radar_ref_obs=4001
  INTEGER,PARAMETER :: id_radar_ref_zero_obs=4004
  INTEGER,PARAMETER :: id_radar_vr_obs=4002
  INTEGER,PARAMETER :: id_radar_prh_obs=4003
!
! Himawari-8 (H08) observations
!
  INTEGER,PARAMETER :: id_H08IR_obs=8800

CONTAINS

!-----------------------------------------------------------------------
! Basic modules for observation input
!-----------------------------------------------------------------------
SUBROUTINE get_nobs(cfile,nn,nrec,endian)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile , endian 
  INTEGER,INTENT(IN)  :: nrec
  INTEGER,INTENT(OUT) :: nn 
  INTEGER :: iunit
  LOGICAL :: ex
  INTEGER :: sz

  nn = 0
  iunit=91
  INQUIRE(FILE=cfile,EXIST=ex)
  IF(ex) THEN
   IF( endian == 'b' )THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential',convert='big_endian')
   ELSE
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential',convert='little_endian')
   ENDIF

!-----------------------------
    INQUIRE(UNIT=iunit, SIZE=sz)
    IF (MOD(sz, r_sngl * (nrec+2)) /= 0) THEN
      WRITE(6,'(3A)') '[Warning]',cfile,': Reading error -- skipped'
      RETURN
    END IF
    nn = sz / (r_sngl * (nrec+2))
!-----------------------------

    WRITE(6,'(I10,A)') nn,' OBSERVATIONS INPUT'
    CLOSE(iunit)
  ELSE
    WRITE(6,'(2A)') cfile,' does not exist -- skipped'
  END IF

  RETURN
END SUBROUTINE get_nobs

SUBROUTINE read_obs(cfile,nobs,obs,nrec,endian)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile , endian
  INTEGER , INTENT(IN)    :: nobs , nrec
  REAL(r_sngl),INTENT(OUT) :: obs(nobs,nrec)
  INTEGER :: n,iunit

  iunit=91
  IF( endian == 'b' )THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential',convert='big_endian')
  ELSE
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential',convert='little_endian')
  ENDIF
  DO n=1,nobs
    READ(iunit) obs(n,:)
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE read_obs

SUBROUTINE write_obs(cfile,nobs,obs,nrec,append,endian)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile , endian
  INTEGER , INTENT(IN)    :: nobs , nrec
  REAL(r_sngl),INTENT(IN) :: obs(nobs,nrec)
  LOGICAL,INTENT(IN)      :: append
  INTEGER :: n,iunit
  CHARACTER(20) :: converte

  iunit=92
  IF( endian == 'b' )THEN
    converte='big_endian'
  ELSE
    converte='little_endian'
  ENDIF

  IF(append) THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='append',convert=converte)
  ELSE
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential',convert=converte)
  END IF
  DO n=1,nobs
    WRITE(iunit) obs(n,:)
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE write_obs


SUBROUTINE get_nobs_radar(cfile,nn,radarlon,radarlat,radarz,nrec,endian)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile , endian
  INTEGER,INTENT(IN)  :: nrec
  INTEGER,INTENT(OUT) :: nn
  REAL(r_sngl) :: tmp
  INTEGER :: ios
  INTEGER :: iunit
  LOGICAL :: ex
  INTEGER :: sz
  REAL(r_sngl),INTENT(OUT) :: radarlon,radarlat,radarz
  CHARACTER(20) :: converte

  nn = 0
  iunit=91

  radarlon=0.0d0
  radarlat=0.0d0
  radarz  =0.0d0

  IF( endian == 'b' )THEN
    converte='big_endian'
  ELSE
    converte='little_endian'
  ENDIF


  INQUIRE(FILE=cfile,EXIST=ex)
  IF(ex) THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential',convert=converte)
    READ(iunit,IOSTAT=ios)tmp
    IF(ios /= 0) THEN
      WRITE(6,'(3A)') '[Warning]',cfile,': Reading error -- skipped'
      RETURN
    END IF
    radarlon=REAL(tmp,r_sngl)
    READ(iunit,IOSTAT=ios)tmp
    IF(ios /= 0) THEN
      WRITE(6,'(3A)') '[Warning]',cfile,': Reading error -- skipped'
      RETURN
    END IF
    radarlat=REAL(tmp,r_sngl)
    READ(iunit,IOSTAT=ios)tmp
    IF(ios /= 0) THEN
      WRITE(6,'(3A)') '[Warning]',cfile,': Reading error -- skipped'
      RETURN
    END IF
    radarz=REAL(tmp,r_sngl)

! get file size by INQUIRE statement... may not work for some older fortran compilers
!-----------------------------
    INQUIRE(UNIT=iunit, SIZE=sz)
    sz = sz - r_sngl * (1+2) * 3 ! substract the radar data header
    IF (MOD(sz, r_sngl * (nrec+2)) /= 0) THEN
      WRITE(6,'(3A)') '[Warning]',cfile,': Reading error -- skipped'
      RETURN
    END IF
    nn = sz / (r_sngl * (nrec+2))
!-----------------------------

    WRITE(6,*)' RADAR FILE ', cfile
    WRITE(6,*)' RADAR LON = ',radarlon
    WRITE(6,*)' RADAR LAT = ',radarlat
    WRITE(6,*)' RADAR Z   = ',radarz
    WRITE(6,'(I10,A)') nn,' OBSERVATIONS INPUT'
    CLOSE(iunit)
  ELSE
    WRITE(6,'(2A)') cfile,' does not exist -- skipped'
  END IF

  RETURN
END SUBROUTINE get_nobs_radar

SUBROUTINE read_obs_radar(cfile,nobs,obs,nrec,endian)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile , endian
  INTEGER , INTENT(IN)    ::  nobs , nrec
  REAL(r_sngl),INTENT(OUT) :: obs(nobs,nrec)
  REAL(r_sngl)      :: tmp
  INTEGER :: n,iunit,ios
  CHARACTER(20) :: converte

  iunit=91

  IF( endian == 'b' )THEN
    converte='big_endian'
  ELSE
    converte='little_endian'
  ENDIF


  OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential',convert=converte)
  READ(iunit, iostat=ios)tmp
  IF(ios /= 0) RETURN
  READ(iunit, iostat=ios)tmp
  IF(ios /= 0) RETURN
  READ(iunit, iostat=ios)tmp
  IF(ios /= 0) RETURN
  DO n=1,nobs
    READ(iunit) obs(n,:)
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE read_obs_radar

SUBROUTINE write_obs_radar(cfile,nobs,rlon,rlat,rz,obs,nrec,append,endian)
  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: cfile  , endian
  INTEGER , INTENT(IN)    ::  nobs , nrec
  REAL(r_sngl),INTENT(IN) :: obs(nobs,nrec) , rlon , rlat , rz
  LOGICAL,INTENT(IN)      :: append
  INTEGER :: n,iunit
  CHARACTER(20) :: converte

  IF( endian == 'b' )THEN
    converte='big_endian'
  ELSE
    converte='little_endian'
  ENDIF


  iunit=92
  IF(append) THEN
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='append',convert=converte)
  ELSE
    OPEN(iunit,FILE=cfile,FORM='unformatted',ACCESS='sequential',convert=converte)
  END IF
  WRITE(iunit) rlon
  WRITE(iunit) rlat
  WRITE(iunit) rz
  DO n=1,nobs
    WRITE(iunit) obs(n,:)
  END DO
  CLOSE(iunit)

  RETURN
END SUBROUTINE write_obs_radar

END MODULE common_obs_scale
