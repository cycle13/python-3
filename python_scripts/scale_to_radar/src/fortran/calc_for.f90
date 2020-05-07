MODULE calc_for
!=======================================================================
!$USE OMP_LIB
  IMPLICIT NONE
  PUBLIC
CONTAINS

!-----------------------------------------------------------------------
! Coordinate conversion (find rk in height levels)
!
! rk = 0.0d0  : surface observation
!-----------------------------------------------------------------------
SUBROUTINE get_k(z_full,nlev,nlon,nlat,ri,rj,rlev,nin,rk)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)

  INTEGER     ,INTENT(IN) :: nlev , nlon , nlat 
  INTEGER     ,INTENT(IN) :: nin
  REAL(r_size),INTENT(IN) :: z_full(nlev,nlon,nlat)
  REAL(r_size),INTENT(IN) :: ri(nin)
  REAL(r_size),INTENT(IN) :: rj(nin)
  REAL(r_size),INTENT(IN) :: rlev(nin) ! height levels
  REAL(r_size),INTENT(OUT) :: rk(nin)
  REAL(r_size) :: ak
  REAL(r_size) :: zlev(nlev)
  INTEGER :: i,j,k, ii, jj, ks  

!$OMP PARALLEL DO PRIVATE(zlev,k,ak,i,j)
DO ii=1,nin 

  if (ri(ii) < 1.0d0 .or. ri(ii) > nlon .or. rj(ii) < 1.0d0 .or. rj(ii) > nlat ) then
    rk(ii) = -1.0
    CYCLE
  end if

  call itpl_2d_column(z_full,nlev,nlon,nlat,ri(ii),rj(ii),zlev)
  !i = int( ri(ii) )
  !j = int( rj(ii) )
  !zlev = z_full( : , i , j )

  !
  ! determine if rlev is within bound.
  !
  IF(rlev(ii) > zlev(nlev)) THEN
    rk(ii) = -1.0
    CYCLE
  END IF
  IF(rlev(ii) < zlev(1)) THEN
    rk(ii) = -1.0
    CYCLE
  END IF
  !
  ! find rk
  !
  DO k=2,nlev
    IF(zlev(k) > rlev(ii) ) EXIT ! assuming ascending order of zlev
  END DO

  ak = (rlev(ii) - zlev(k-1)) / (zlev(k) - zlev(k-1))
  rk(ii) = REAL(k-1,r_size) + ak
  !write( *,* ) rk(ii) 

ENDDO
!$OMP END PARALLEL DO

RETURN
END SUBROUTINE get_k

!-----------------------------------------------------------------------
SUBROUTINE itpl_2d(var,nlon,nlat,ri,rj,var5)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)

  INTEGER     ,INTENT(IN) :: nlon , nlat
  REAL(r_size),INTENT(IN) :: var(nlon,nlat)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(OUT) :: var5
  REAL(r_size) :: ai,aj
  INTEGER :: i,j

  i = CEILING(ri)
  ai = ri - REAL(i-1,r_size)
  j = CEILING(rj)
  aj = rj - REAL(j-1,r_size)

  var5 = var(i-1,j-1) * (1-ai) * (1-aj) &
     & + var(i  ,j-1) *    ai  * (1-aj) &
     & + var(i-1,j  ) * (1-ai) *    aj  &
     & + var(i  ,j  ) *    ai  *    aj

  RETURN
END SUBROUTINE itpl_2d

SUBROUTINE itpl_2d_column(var,nlev,nlon,nlat,ri,rj,var5)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)

  INTEGER     ,INTENT(IN) :: nlev,nlon,nlat
  REAL(r_size),INTENT(IN) :: var(nlev,nlon,nlat)
  REAL(r_size),INTENT(IN) :: ri
  REAL(r_size),INTENT(IN) :: rj
  REAL(r_size),INTENT(OUT) :: var5(nlev)
  REAL(r_size) :: ai,aj
  INTEGER :: i,j

  i = CEILING(ri)
  ai = ri - REAL(i-1,r_size)
  j = CEILING(rj)
  aj = rj - REAL(j-1,r_size)

  var5(:) = var(:,i-1,j-1) * (1-ai) * (1-aj) &
        & + var(:,i  ,j-1) *    ai  * (1-aj) &
        & + var(:,i-1,j  ) * (1-ai) *    aj  &
        & + var(:,i  ,j  ) *    ai  *    aj

  RETURN
END SUBROUTINE itpl_2d_column


END MODULE calc_for
