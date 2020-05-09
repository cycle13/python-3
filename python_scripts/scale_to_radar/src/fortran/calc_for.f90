MODULE calc_for
!=======================================================================
!$USE OMP_LIB
  IMPLICIT NONE

  PUBLIC
CONTAINS


!-----------------------------------------------------------------------
! Coordinate conversion (find rk in height levels)
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

  rk=0.0d0
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

SUBROUTINE itpl_3d_vec(var,nlev,nlon,nlat,ri,rj,rk,fill_value,nin,var5)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)

  INTEGER     ,INTENT(IN) :: nlev,nlon,nlat,nin
  REAL(r_size),INTENT(IN) :: var(nlev,nlon,nlat)
  REAL(r_size),INTENT(IN) :: ri(nin)
  REAL(r_size),INTENT(IN) :: rj(nin)
  REAL(r_size),INTENT(IN) :: rk(nin)
  REAL(r_size),INTENT(IN) :: fill_value 
  REAL(r_size),INTENT(OUT) :: var5(nin)
  REAL(r_size) :: ai,aj,ak
  INTEGER :: i,j,k,ii

  !$OMP PARALLEL DO PRIVATE(i,j,k,ai,aj,ak)
  DO ii=1,nin
     
     i = CEILING(ri(ii))
     if ( i <= 1 .or. i > nlon )THEN
        var5(ii) = fill_value
        cycle
     endif
     j = CEILING(rj(ii))
     if ( j <= 1 .or. j > nlat )THEN
        var5(ii) = fill_value
        cycle
     endif
     k = CEILING(rk(ii))
     if ( k <= 1 .or. k > nlev )THEN
        var5(ii) = fill_value
        cycle
     endif
         
     ai = ri(ii) - REAL(i-1,r_size)
     aj = rj(ii) - REAL(j-1,r_size)
     ak = rk(ii) - REAL(k-1,r_size)

     var5(ii) = var(k-1,i-1,j-1) * (1-ai) * (1-aj) * (1-ak) &
     & + var(k-1,i  ,j-1) *    ai  * (1-aj) * (1-ak) &
     & + var(k-1,i-1,j  ) * (1-ai) *    aj  * (1-ak) &
     & + var(k-1,i  ,j  ) *    ai  *    aj  * (1-ak) &
     & + var(k  ,i-1,j-1) * (1-ai) * (1-aj) *    ak  &
     & + var(k  ,i  ,j-1) *    ai  * (1-aj) *    ak  &
     & + var(k  ,i-1,j  ) * (1-ai) *    aj  *    ak  &
     & + var(k  ,i  ,j  ) *    ai  *    aj  *    ak

  END DO
  !$OMP END PARALLEL DO


  RETURN
END SUBROUTINE itpl_3d_vec



!2D interpolation using box average. Destination grid is assumed to be regular.
SUBROUTINE com_interp_boxavereg(xini,dx,nx,yini,dy,ny,zini,dz,nz,nvar,xin,yin,zin,datain,undef,nin    &
               &                ,data_ave,data_max,data_min,data_std,data_n)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)

  INTEGER , INTENT(IN)           :: nx , ny , nz , nvar , nin
  REAL(r_size),INTENT(IN)        :: dx , dy , dz , xini , yini , zini
  REAL(r_size),INTENT(IN)        :: xin(nin),yin(nin),zin(nin),datain(nin,nvar)
  REAL(r_size),INTENT(IN)        :: undef
  REAL(r_size),INTENT(OUT)       :: data_ave(nz,ny,nx,nvar)
  REAL(r_size),INTENT(OUT)       :: data_max(nz,ny,nx,nvar)
  REAL(r_size),INTENT(OUT)       :: data_min(nz,ny,nx,nvar)
  REAL(r_size),INTENT(OUT)       :: data_std(nz,ny,nx,nvar)
  INTEGER,INTENT(OUT)            :: data_n(nz,ny,nx,nvar)
  REAL(r_size)                   :: tmp_data 

  INTEGER                        :: ii , ix , iy , iz , iv , tmp_n

data_max=-1.0d9
data_min=1.0d9
data_ave=0.0d0
data_std=0.0d0
data_n  =0

!$OMP PARALLEL DO PRIVATE( tmp_data,ii,ix,iy,iz )
DO iv = 1,nvar !Loop over the variables (we can perform OMP over this loop)

  DO ii = 1,nin  !Loop over the input data 

    !Compute the location of the current point with in grid coordinates (rx,ry)
    ix = int( ( xin(ii) - xini ) / dx ) + 1
    iy = int( ( yin(ii) - yini ) / dy ) + 1
    iz = int( ( zin(ii) - zini ) / dz ) + 1

    !Check is the data is within the grid.
    IF( ix <= nx .and. ix >= 1 .and. iy <= ny .and. iy >= 1 .and.   &
        iz <= nz .and. iz >= 1 .and. datain(ii,iv) /= undef )THEN

        data_n(iz,iy,ix,iv) = data_n(iz,iy,ix,iv) + 1
        data_ave(iz,iy,ix,iv) = data_ave(iz,iy,ix,iv) + datain(ii,iv) 
        data_std(iz,iy,ix,iv) = data_std(iz,iy,ix,iv) + ( datain(ii,iv) ** 2 )

        IF( datain(ii,iv) > data_max(iz,iy,ix,iv) )THEN
          data_max(iz,iy,ix,iv) = datain(ii,iv)
        ENDIF
        IF( datain(ii,iv) < data_min(iz,iy,ix,iv) )THEN
          data_min(iz,iy,ix,iv) = datain(ii,iv)
        ENDIF

    ENDIF

  ENDDO

  WHERE( data_n(:,:,:,iv) > 0)
       data_ave(:,:,:,iv) = data_ave(:,:,:,iv) / REAL( data_n(:,:,:,iv) , r_sngl )
       data_std(:,:,:,iv) = SQRT( data_std(:,:,:,iv)/REAL( data_n(:,:,:,iv) , r_sngl ) - data_ave(:,:,:,iv) ** 2 )
  ENDWHERE
  WHERE( data_n(:,:,:,iv) == 0)
       data_ave(:,:,:,iv) = undef
       data_std(:,:,:,iv) = undef
       data_min(:,:,:,iv) = undef
       data_max(:,:,:,iv) = undef
  ENDWHERE

ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE com_interp_boxavereg



END MODULE calc_for
