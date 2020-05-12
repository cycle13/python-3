!MODULE common
!=======================================================================
!
! [PURPOSE:] General constants and procedures
!
! [HISTORY:]
!
!s=======================================================================
!  IMPLICIT NONE
!  PUBLIC
!-----------------------------------------------------------------------
! Variable size definitions
!-----------------------------------------------------------------------

!CONTAINS
!-----------------------------------------------------------------------
! Mean
!-----------------------------------------------------------------------
MODULE common_functions

!  INTEGER,PARAMETER :: r_size=kind(0.0d0)
!  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)

!-----------------------------------------------------------------------
! Constants
!-----------------------------------------------------------------------
!  REAL(r_size),PARAMETER :: pi=3.141592653589793d0
!  REAL(r_size),PARAMETER :: gg=9.81d0
!  REAL(r_size),PARAMETER :: rd=287.0d0
!  REAL(r_size),PARAMETER :: cp=7.0d0 / 2.0d0 * rd
!  REAL(r_size),PARAMETER :: re=6371.3d3
!  REAL(r_size),PARAMETER :: r_omega=7.292d-5
!  REAL(r_size),PARAMETER :: t0c=273.15d0
!  REAL(r_size),PARAMETER :: undef=9.99d33
!  REAL(r_sngl),PARAMETER :: undefs=9.99e33
!  REAL(r_size),PARAMETER :: deg2rad=3.1415926535d0/180d0
!  REAL(r_size),PARAMETER :: rad2deg=180d0/3.1415926535d0
!  REAL(r_size),PARAMETER :: clight=299792458.0d0 !Speed of light

PUBLIC

contains

SUBROUTINE com_mean(ndim,var,amean)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)

  INTEGER,INTENT(IN) :: ndim
  REAL(r_size),INTENT(IN) :: var(ndim)
  REAL(r_size),INTENT(OUT) :: amean

  INTEGER :: i

  amean = 0.0d0
  DO i=1,ndim
    amean = amean + var(i)
  END DO
  amean = amean / REAL(ndim,r_size)

  RETURN
END SUBROUTINE com_mean

SUBROUTINE com_mean_sngl(ndim,var,amean,w)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)

  INTEGER,INTENT(IN) :: ndim
  REAL(r_sngl),INTENT(IN) :: var(ndim)
  REAL(r_sngl),INTENT(IN) :: w(ndim)
  REAL(r_sngl),INTENT(OUT) :: amean

  INTEGER :: i

  !IF( .not. present( w ) )THEN
  !   w=1.0e0
  !ENDIF

  amean = 0.0d0
  DO i=1,ndim
    amean = amean + var(i) * w(i)
  END DO
  amean = amean / SUM( w )

  RETURN
END SUBROUTINE com_mean_sngl


!-----------------------------------------------------------------------
! Standard deviation
!-----------------------------------------------------------------------
SUBROUTINE com_stdev(ndim,var,aout)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)

  INTEGER,INTENT(IN) :: ndim
  REAL(r_size),INTENT(IN) :: var(ndim)
  REAL(r_size),INTENT(OUT) :: aout

  REAL(r_size) :: amean
  REAL(r_size) :: dev(ndim)

  CALL com_mean(ndim,var,amean)

  dev(:) = var(:) - amean

  aout = SQRT( SUM(dev*dev) / REAL(ndim-1,r_size) )

  RETURN
END SUBROUTINE com_stdev

SUBROUTINE com_stdev_sngl(ndim,var,aout,w)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)

  INTEGER,INTENT(IN) :: ndim
  REAL(r_sngl),INTENT(IN) :: var(ndim)
  REAL(r_sngl),INTENT(IN) :: w(ndim)
  REAL(r_sngl),INTENT(OUT) :: aout

  REAL(r_sngl) :: amean
  REAL(r_sngl) :: dev(ndim)

  !IF( .not. present( w ) )THEN
  !  w=1.0e0
  !ENDIF

  CALL com_mean_sngl(ndim,var,amean,w)

  dev(:) = var(:) - amean 

  aout = SQRT( SUM(w*dev*dev) / ( sum(w)*(ndim-1)/real(ndim,r_sngl) ) )

  RETURN
END SUBROUTINE com_stdev_sngl

!-----------------------------------------------------------------------
! Covariance
!-----------------------------------------------------------------------
SUBROUTINE com_covar(ndim,var1,var2,cov)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)

  INTEGER,INTENT(IN) :: ndim
  REAL(r_size),INTENT(IN) :: var1(ndim)
  REAL(r_size),INTENT(IN) :: var2(ndim)
  REAL(r_size),INTENT(OUT) :: cov

  REAL(r_size) :: amean1,amean2
  REAL(r_size) :: dev1(ndim),dev2(ndim)

  CALL com_mean(ndim,var1,amean1)
  CALL com_mean(ndim,var2,amean2)

  dev1(:) = var1(:) - amean1
  dev2(:) = var2(:) - amean2

  cov = SUM( dev1*dev2 ) / REAL(ndim-1,r_size)

  RETURN
END SUBROUTINE com_covar

SUBROUTINE com_covar_sngl(ndim,var1,var2,cov,w)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)

  INTEGER,INTENT(IN) :: ndim
  REAL(r_sngl),INTENT(IN) :: var1(ndim)
  REAL(r_sngl),INTENT(IN) :: var2(ndim)
  REAL(r_sngl),INTENT(IN) :: w(ndim)
  REAL(r_sngl),INTENT(OUT) :: cov

  REAL(r_sngl) :: amean1,amean2
  REAL(r_sngl) :: dev1(ndim),dev2(ndim)

  CALL com_mean_sngl(ndim,var1,amean1,w)
  CALL com_mean_sngl(ndim,var2,amean2,w)

  dev1(:) = var1(:) - amean1
  dev2(:) = var2(:) - amean2

  cov = SUM( dev1*dev2*w ) / ( SUM(w) )

  RETURN
END SUBROUTINE com_covar_sngl

SUBROUTINE com_covar_sngl_sample(ndim,var1,var2,sampleindex,cov,w)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)

  INTEGER,INTENT(IN) :: ndim
  REAL(r_sngl),INTENT(IN) :: var1(ndim)
  REAL(r_sngl),INTENT(IN) :: var2(ndim)
  INTEGER     ,INTENT(IN) :: sampleindex(ndim)
  REAL(r_sngl),INTENT(IN) :: w(ndim)
  REAL(r_sngl),INTENT(OUT) :: cov
  INTEGER  :: ii
  REAL(r_sngl) :: amean1,amean2
  REAL(r_sngl) :: dev1(ndim),dev2(ndim),devw(ndim)

  !Reorder the sample
  DO ii=1,ndim
    dev1(ii) = var1(sampleindex(ii))
    dev2(ii) = var2(sampleindex(ii))
    devw(ii) = w(sampleindex(ii))

  ENDDO

  CALL com_mean_sngl(ndim,dev1,amean1,devw)
  CALL com_mean_sngl(ndim,dev2,amean2,devw)

  dev1=dev1-amean1
  dev2=dev2-amean2

  cov = SUM( dev1*dev2*devw ) / ( SUM(w) * REAL(ndim-1,r_sngl) )

  RETURN
END SUBROUTINE com_covar_sngl_sample


!-----------------------------------------------------------------------
! Correlation
!-----------------------------------------------------------------------
SUBROUTINE com_correl(ndim,var1,var2,cor)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)

  INTEGER,INTENT(IN) :: ndim
  REAL(r_size),INTENT(IN) :: var1(ndim)
  REAL(r_size),INTENT(IN) :: var2(ndim)
  REAL(r_size),INTENT(OUT) :: cor

  REAL(r_size) :: cov,stdev1,stdev2

  CALL com_stdev(ndim,var1,stdev1)
  CALL com_stdev(ndim,var2,stdev2)
  CALL com_covar(ndim,var1,var2,cov)

  cor = cov/stdev1/stdev2

  RETURN
END SUBROUTINE com_correl

SUBROUTINE com_correl_sngl(ndim,var1,var2,cor,w)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)

  INTEGER,INTENT(IN) :: ndim
  REAL(r_sngl),INTENT(IN) :: var1(ndim)
  REAL(r_sngl),INTENT(IN) :: var2(ndim)
  REAL(r_sngl),INTENT(OUT) :: cor
  REAL(r_sngl),INTENT(IN)  :: w(ndim)

  REAL(r_sngl) :: cov,stdev1,stdev2

  CALL com_stdev_sngl(ndim,var1,stdev1,w)
  CALL com_stdev_sngl(ndim,var2,stdev2,w)
  CALL com_covar_sngl(ndim,var1,var2,cov,w)

  cor = cov/stdev1/stdev2

  RETURN
END SUBROUTINE com_correl_sngl

!-----------------------------------------------------------------------
! Mode computation
!-----------------------------------------------------------------------

RECURSIVE SUBROUTINE com_mode_sngl(ndim,var,mode,w,is_sorted)
IMPLICIT NONE 
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
INTEGER,  INTENT(IN)     :: ndim
REAL(r_sngl), INTENT(IN) :: var(ndim)
REAL(r_sngl), INTENT(OUT):: mode
REAL(r_sngl), INTENT(IN) :: w(ndim,1)
LOGICAL     , INTENT(IN) :: is_sorted
REAL(r_sngl)             :: vars(ndim) , ws(ndim,1)
REAL(r_sngl)             :: max_density , sum_weigth , delta_x , density , undef
INTEGER                  :: mdi , imax , ii

!This algoritm provides a recursive estimation of the mode.
!Bickel and Fruhwirth 2005 https://arxiv.org/pdf/math/0505419.pdf
!This implementation considers a weigthed sample through the optional w
!array that contains the weigth corresponding to each sample element.
undef=-9999e9
!first sort the array in ascending order.
vars=var
ws=w
mode=undef

IF ( .not. is_sorted )THEN
   !sort the variable and the weights
   CALL com_cosort_sngl(ndim,var,w,1,vars,ws)
   !write(*,*) minval( ws ) , maxval( ws ) , minval( w ) , maxval( w )
   
ENDIF


!Special cases ndim from 1 to 3.
IF( ndim == 1 )THEN
  mode = vars(1)
  RETURN
ENDIF
IF( ndim == 2 )THEN
  mode = (vars(1) + vars(2) )/2.0e0
  RETURN
ENDIF
IF( ndim == 3 )THEN
  IF( abs( vars(1) - vars(2) ) < abs( vars(2) - vars(3) ) )THEN
     mode = ( vars(1) + vars(2) ) / 2.0
     RETURN
  ELSE
     mode = ( vars(2) + vars(3) ) / 2.0
  ENDIF
ENDIF


IF ( ndim >= 4 )THEN

   max_density = 0

   imax = ceiling( ndim / real(2.0,r_sngl) ) 

   !Search for the interval with the maximum density.
   DO ii = 1,ndim-imax+1

      density = sum( ws(ii:ii+imax-1,1) )/(vars(ii+imax-1)-vars(ii))

      IF( density > max_density )THEN
         max_density = density
         mdi = ii 
      ENDIF

   ENDDO

   !Recursively call this subroutine but only on the segment that shows
   !the largest pdf density.
   CALL com_mode_sngl(imax,vars(mdi:mdi+imax-1),mode,ws(mdi:mdi+imax-1,1),.True.)

ENDIF

END SUBROUTINE com_mode_sngl

!-----------------------------------------------------------------------
! Array sort
!-----------------------------------------------------------------------

! quicksort.f -*-f90-*-
! Author: t-nissie, some tweaks by 1AdAstra1
! License: GPLv3
! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!!
RECURSIVE SUBROUTINE com_sort_sngl( ndim , var , vars )
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
  INTEGER     , INTENT(IN) :: ndim
  REAL(r_sngl), INTENT(IN) :: var(ndim)
  REAL(r_sngl), INTENT(OUT):: vars(ndim)
  REAL                     :: x , t
  INTEGER                  :: first = 1, last
  INTEGER                  :: i, j

  vars=var
  last = ndim
  x = vars( (first+last) / 2 )
  i = first
  j = last
  
  do
     do while (vars(i) < x)
        i=i+1
     end do
     do while (x < vars(j))
        j=j-1
     end do
     if (i >= j) exit
     t = vars(i);  vars(i) = vars(j);  vars(j) = t
     i=i+1
     j=j-1
  end do
  
  if (first < i - 1) CALL com_sort_sngl(i-first,vars(first : i - 1),vars(first : i - 1))
  if (j + 1 < last)  CALL com_sort_sngl(last-j,vars(j + 1 : last),vars(j + 1 : last))

END SUBROUTINE com_sort_sngl

!-----------------------------------------------------------------------
! Array co - sort
!-----------------------------------------------------------------------
! co sort based on 
! quicksort.f -*-f90-*-
! Author: t-nissie, some tweaks by 1AdAstra1
! License: GPLv3
! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!!
RECURSIVE SUBROUTINE com_cosort_sngl( ndim , var , covar , ncovar , vars , covars )
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
  INTEGER     , INTENT(IN) :: ndim
  REAL(r_sngl), INTENT(IN) :: var(ndim)
  REAL(r_sngl), INTENT(OUT):: vars(ndim)
  INTEGER                  :: ncovar             !Number of co-sorted variables
  REAL(r_sngl), INTENT(IN) :: covar(ndim,ncovar) !Variables to be co-sorted
  REAL(r_sngl), INTENT(OUT):: covars(ndim,ncovar)!Co sorted variables
  REAL                     :: x , t , cot(ncovar)
  INTEGER                  :: first = 1, last
  INTEGER                  :: i, j

  vars=var
  covars=covar
  last = ndim
  x = vars( (first+last) / 2 )
  i = first
  j = last

  do
     do while (vars(i) < x)
        i=i+1
     end do
     do while (x < vars(j))
        j=j-1
     end do
     if (i >= j) exit
     t = vars(i);  vars(i) = vars(j);  vars(j) = t
     cot = covars(i,:); covars(i,:) = covars(j,:) ; covars(j,:)=cot 
     i=i+1
     j=j-1
  end do

  if (first < i - 1) CALL com_cosort_sngl(i-first,vars(first : i - 1),covars(first:i-1,:), & 
                                        ncovar,vars(first : i - 1),covars(first:i-1,:))
  if (j + 1 < last)  CALL com_cosort_sngl(last-j,vars(j + 1 : last),covars(j+1: last,:),   & 
                                        ncovar,vars(j + 1 : last),covars(j+1:last,:))

END SUBROUTINE com_cosort_sngl



!-----------------------------------------------------------------------
! Anomaly Correlation
!-----------------------------------------------------------------------
SUBROUTINE com_anomcorrel(ndim,var1,var2,varmean,cor)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
  INTEGER,INTENT(IN) :: ndim
  REAL(r_size),INTENT(IN) :: var1(ndim)
  REAL(r_size),INTENT(IN) :: var2(ndim)
  REAL(r_size),INTENT(IN) :: varmean(ndim)
  REAL(r_size),INTENT(OUT) :: cor

  REAL(r_size) :: dev1(ndim),dev2(ndim)

  dev1 = var1 - varmean
  dev2 = var2 - varmean

  cor = SUM( dev1*dev2 ) / SQRT( SUM(dev1*dev1) * SUM(dev2*dev2) )

  RETURN
END SUBROUTINE com_anomcorrel
!-----------------------------------------------------------------------
! L2 Norm
!-----------------------------------------------------------------------
SUBROUTINE com_l2norm(ndim,var,anorm)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
  INTEGER,INTENT(IN) :: ndim
  REAL(r_size),INTENT(IN) :: var(ndim)
  REAL(r_size),INTENT(OUT) :: anorm

  anorm = SQRT( SUM(var*var) )

  RETURN
END SUBROUTINE com_l2norm
!-----------------------------------------------------------------------
! RMS (root mean square)
!-----------------------------------------------------------------------
SUBROUTINE com_rms(ndim,var,undef,rmsv)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
  INTEGER,INTENT(IN) :: ndim
  REAL(r_sngl),INTENT(IN) :: var(ndim)
  REAL(r_sngl),INTENT(IN) :: undef
  REAL(r_sngl),INTENT(OUT) :: rmsv
  REAL(r_sngl)             :: rmsv_tmp
  INTEGER                  :: nt
  INTEGER                  :: ii

  rmsv=undef
  rmsv_tmp=0.0d0
  nt=0
  
  !$OMP PARALLEL DO REDUCTION(+:rmsv_tmp,nt)
  DO ii = 1 , ndim
     IF( var(ii) /= undef )THEN
        rmsv_tmp = rmsv_tmp + var(ii)**2
        nt = nt + 1
     ENDIF
  ENDDO
  !$OMP END PARALLEL DO
  IF( nt > 0 )THEN
    rmsv = SQRT( rmsv_tmp / REAL( nt , r_sngl ) )
  ENDIF

  RETURN
END SUBROUTINE com_rms
!-----------------------------------------------------------------------
! Lanczos Filter (Low-pass) with cyclic boundary
!-----------------------------------------------------------------------
SUBROUTINE com_filter_lanczos(ndim,fc,var)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
  REAL(r_size),PARAMETER :: pi=3.141592653589793d0
  REAL(r_size),PARAMETER :: deg2rad=3.1415926535d0/180d0
  REAL(r_size),PARAMETER :: rad2deg=180d0/3.1415926535d0
  REAL(r_size),PARAMETER :: re=6371.3d3
  INTEGER,INTENT(IN) :: ndim
  REAL(r_size),INTENT(IN) :: fc    ! critical frequency in [0,pi]
  REAL(r_size),INTENT(INOUT) :: var(ndim)


  INTEGER,PARAMETER :: lresol=10

  REAL(r_size) :: weight(-lresol:lresol)
  REAL(r_size) :: varwk(1-lresol:ndim+lresol)
  REAL(r_size) :: rl,rlresol
  INTEGER :: i,l
!
! Weight
!
  rlresol = REAL(lresol,r_size)
  DO l=-lresol,-1
    rl = REAL(l,r_size)
    weight(l) = SIN(fc*rl) * SIN(pi*rl/rlresol) &
      & * rlresol / pi / rl / pi / rl
  END DO
  DO l=1,lresol
    rl = REAL(l,r_size)
    weight(l) = SIN(fc*rl) * SIN(pi*rl/rlresol) &
      & * rlresol / pi / rl / pi / rl
  END DO
  weight(0) = fc / pi
!
! Cyclic boundary
!
  DO i=0,1-lresol,-1
    varwk(i) = var(ndim+i)
  END DO
  DO i=ndim+1,ndim+lresol
    varwk(i) = var(i-ndim)
  END DO
  varwk(1:ndim) = var(1:ndim)
!
! Filter
!
  var = 0.0d0
  DO i=1,ndim
    DO l=-lresol,lresol
      var(i) = var(i) + weight(l) * varwk(i+l)
    END DO
  END DO

  RETURN
END SUBROUTINE com_filter_lanczos

!-----------------------------------------------------------------------
! DISTANCE BETWEEN TWO POINTS (LONa,LATa)-(LONb,LATb)
!-----------------------------------------------------------------------
SUBROUTINE com_distll(ndim,alon,alat,blon,blat,dist)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
  REAL(r_size),PARAMETER :: pi=3.141592653589793d0
  REAL(r_size),PARAMETER :: deg2rad=3.1415926535d0/180d0
  REAL(r_size),PARAMETER :: rad2deg=180d0/3.1415926535d0
  REAL(r_size),PARAMETER :: re=6371.3d3
  INTEGER,INTENT(IN) :: ndim
  REAL(r_size),INTENT(IN) :: alon(ndim)
  REAL(r_size),INTENT(IN) :: alat(ndim)
  REAL(r_size),INTENT(IN) :: blon(ndim)
  REAL(r_size),INTENT(IN) :: blat(ndim)
  REAL(r_size),INTENT(OUT) :: dist(ndim)
  REAL(r_size),PARAMETER :: r180=1.0d0/180.0d0
  REAL(r_size) :: lon1,lon2,lat1,lat2
  REAL(r_size) :: cosd(ndim)
  INTEGER :: i


  DO i=1,ndim
    lon1 = alon(i) * pi * r180
    lon2 = blon(i) * pi * r180
    lat1 = alat(i) * pi * r180
    lat2 = blat(i) * pi * r180

    cosd(i) = SIN(lat1)*SIN(lat2) + COS(lat1)*COS(lat2)*COS(lon2-lon1)
    cosd(i) = MIN( 1.d0,cosd(i))
    cosd(i) = MAX(-1.d0,cosd(i))

    dist(i) = ACOS( cosd(i) ) * re
  END DO

  RETURN
END SUBROUTINE com_distll
!-----------------------------------------------------------------------
! DISTANCE BETWEEN TWO POINTS (LONa,LATa)-(LONb,LATb)
!-----------------------------------------------------------------------
SUBROUTINE com_distll_1(alon,alat,blon,blat,dist)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
  REAL(r_size),PARAMETER :: pi=3.141592653589793d0
  REAL(r_size),PARAMETER :: deg2rad=3.1415926535d0/180d0
  REAL(r_size),PARAMETER :: rad2deg=180d0/3.1415926535d0
  REAL(r_size),PARAMETER :: re=6371.3d3
  REAL(r_size),INTENT(IN) :: alon
  REAL(r_size),INTENT(IN) :: alat
  REAL(r_size),INTENT(IN) :: blon
  REAL(r_size),INTENT(IN) :: blat
  REAL(r_size),INTENT(OUT) :: dist
  REAL(r_size),PARAMETER :: r180=1.0d0/180.0d0
  REAL(r_size) :: lon1,lon2,lat1,lat2
  REAL(r_size) :: cosd
 
  lon1 = alon * pi * r180
  lon2 = blon * pi * r180
  lat1 = alat * pi * r180
  lat2 = blat * pi * r180

  cosd = SIN(lat1)*SIN(lat2) + COS(lat1)*COS(lat2)*COS(lon2-lon1)
  cosd = MIN( 1.d0,cosd)
  cosd = MAX(-1.d0,cosd)

  dist = ACOS( cosd ) * re

  RETURN
END SUBROUTINE com_distll_1

SUBROUTINE com_distll_sngl(alon,alat,blon,blat,dist)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
  REAL(r_size),PARAMETER :: pi=3.141592653589793d0
  REAL(r_size),PARAMETER :: deg2rad=3.1415926535d0/180d0
  REAL(r_size),PARAMETER :: rad2deg=180d0/3.1415926535d0
  REAL(r_size),PARAMETER :: re=6371.3d3
  REAL(r_sngl),INTENT(IN) :: alon
  REAL(r_sngl),INTENT(IN) :: alat
  REAL(r_sngl),INTENT(IN) :: blon
  REAL(r_sngl),INTENT(IN) :: blat
  REAL(r_sngl),INTENT(OUT) :: dist
  REAL(r_sngl),PARAMETER :: r180=1.0d0/180.0d0
  REAL(r_size) :: lon1,lon2,lat1,lat2,cosd
  !Double precission is esential in internal computations.
  !So interface is single precission, but internal computations
  !are performed with double precission.

  lon1 = REAL(alon,r_size) * deg2rad
  lon2 = REAL(blon,r_size) * deg2rad
  lat1 = REAL(alat,r_size) * deg2rad
  lat2 = REAL(blat,r_size) * deg2rad

  cosd = SIN(lat1)*SIN(lat2) + COS(lat1)*COS(lat2)*COS(lon2-lon1)
  cosd = MIN( 1.d0,cosd)
  cosd = MAX(-1.d0,cosd)

  dist = REAL( ACOS( cosd ) * re , r_sngl)

  RETURN
END SUBROUTINE com_distll_sngl

SUBROUTINE com_compassbearing_sngl(alon,alat,blon,blat,bearing)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
  REAL(r_size),PARAMETER :: pi=3.141592653589793d0
  REAL(r_size),PARAMETER :: deg2rad=3.1415926535d0/180d0
  REAL(r_size),PARAMETER :: rad2deg=180d0/3.1415926535d0
  REAL(r_size),PARAMETER :: re=6371.3d3
  REAL(r_sngl),INTENT(IN) :: alon
  REAL(r_sngl),INTENT(IN) :: alat
  REAL(r_sngl),INTENT(IN) :: blon
  REAL(r_sngl),INTENT(IN) :: blat
  REAL(r_sngl),INTENT(OUT) :: bearing
  REAL(r_size) :: lon1,lon2,lat1,lat2,dlon,x,y

    lat1 = REAL(alat,r_size) * deg2rad
    lat2 = REAL(blat,r_size) * deg2rad

    dlon = REAL(blon - alon , r_size) * deg2rad

    x = sin(dlon) * cos(lat2)
    y = cos(lat1) * sin(lat2) - (sin(lat1) * cos(lat2) * cos(dlon))

    bearing = atan2(x, y)

    bearing = REAL( bearing * rad2deg , r_sngl )

    if ( bearing < 0.0e0 ) then
       bearing = bearing + 360.0e0
    endif

END SUBROUTINE com_compassbearing_sngl

!-----------------------------------------------------------------------
! Lon , lat moving in a certain direction.
!-----------------------------------------------------------------------

!-------------------------------------------------------------------------------------
! Compute the lat and lon reached by moving a certain distance in a certain
! direction
!        Formula from Map Projections - a working manual.  USGS paper
!        1395.  Equations (5-5) and (5-6).
!-------------------------------------------------------------------------------------
! Input azimuth in degrees with respect to North
! Input distance in meters
! Input and output lat lon in degrees

SUBROUTINE com_ll_arc_distance_sngl(ini_lon,ini_lat,distance,azimuth,final_lon,final_lat)
IMPLICIT NONE
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
REAL(r_size),PARAMETER :: pi=3.141592653589793d0
REAL(r_size),PARAMETER :: deg2rad=3.1415926535d0/180d0
REAL(r_size),PARAMETER :: rad2deg=180d0/3.1415926535d0
REAL(r_size),PARAMETER :: re=6371.3d3
REAL(r_sngl), INTENT(IN)  :: ini_lon,ini_lat,distance,azimuth
REAL(r_sngl), INTENT(OUT) :: final_lon,final_lat
REAL(r_size)  :: cdist,sdist,sinll1,cosll1,cosaz,sinaz
!So interface is single precission, but internal computations
!are performed with double precission.


 IF ( distance .EQ. 0e0 )THEN
    final_lon=ini_lon
    final_lat=ini_lat
 ELSE

 cdist = cos(REAL(distance,r_size)/re)
 sdist = sin(REAL(distance,r_size)/re)

 sinll1 = sin(REAL(ini_lat,r_size)*deg2rad)
 cosll1 = cos(REAL(ini_lat,r_size)*deg2rad)

 cosaz  = cos(REAL(azimuth,r_size)*deg2rad)
 sinaz  = sin(REAL(azimuth,r_size)*deg2rad)

 final_lat=asin(sinll1*cdist + cosll1*sdist*cosaz )*rad2deg
 final_lon=ini_lon + rad2deg *atan2(sdist*sinaz,cosll1*cdist - &
                  &  sinll1*sdist*cosaz)

 ENDIF

END SUBROUTINE com_ll_arc_distance_sngl

!Convert lat,lon to range,azimuth (with respect to a certain location)
SUBROUTINE com_ll_to_ra(cen_lon,cen_lat,lon,lat,nx,ny,r,az)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
  INTEGER     ,INTENT(IN)    :: nx,ny
  REAL(r_sngl),INTENT(IN)    :: cen_lon , cen_lat
  REAL(r_sngl),INTENT(IN)    :: lon(nx,ny) , lat(nx,ny)
  REAL(r_sngl),INTENT(OUT)   :: r(nx,ny)   , az(nx,ny)
  INTEGER                    :: ii , jj
  !Compute the azimuth and range from for each element of lon and lat arrays
  !as seen from the location defined by cen_lon , cen_lat
 
  
  !$OMP PARALLEL DO PRIVATE(ii,jj)
 
  DO ii=1,nx
    DO jj=1,ny

     CALL com_distll_sngl(cen_lon,cen_lat,lon(ii,jj),lat(ii,jj),r(ii,jj))
     CALL com_compassbearing_sngl(cen_lon,cen_lat,lon(ii,jj),lat(ii,jj),az(ii,jj))

    ENDDO
  ENDDO

  !$OMP END PARALLEL DO

END SUBROUTINE com_ll_to_ra

!Convert range,azimuth with respect to a certain location to lat,lon
SUBROUTINE com_ra_to_ll(cen_lon,cen_lat,r,a,nr,na,lon,lat)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
  INTEGER     ,INTENT(IN)    :: nr,na
  REAL(r_sngl),INTENT(IN)    :: cen_lon , cen_lat
  REAL(r_sngl),INTENT(OUT)   :: lon(nr,na) , lat(nr,na)
  REAL(r_sngl),INTENT(IN)    :: r(nr,na)   , a(nr,na)
  INTEGER                    :: ii , jj
  !Compute the azimuth and range from for each element of lon and lat arrays
  !as seen from the location defined by cen_lon , cen_lat

  !$OMP PARALLEL DO PRIVATE(ii,jj)

  DO ii=1,nr
    DO jj=1,na
      
     CALL com_ll_arc_distance_sngl(cen_lon,cen_lat,r(ii,jj),a(ii,jj),lon(ii,jj),lat(ii,jj) )

    ENDDO
  ENDDO

  !$OMP END PARALLEL DO

END SUBROUTINE com_ra_to_ll

!Complete missing values in a 2D field by successive interpolation from the
!neigbhor values.
SUBROUTINE com_complete_missing_2d(field,missing_mask,nx,ny,npass,cfield) 
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
  INTEGER  , INTENT(IN)     :: nx , ny
  REAL(r_sngl),INTENT(IN)   :: field(nx,ny)         !Data field
  REAL(r_sngl),INTENT(OUT)  :: cfield(nx,ny)        !Completed data field.
  LOGICAL     ,INTENT(IN)   :: missing_mask(nx,ny)  !Mask indicating missing values
                                                    !True means missing , false valid
  LOGICAL                   :: missing_maskc(nx,ny)
  INTEGER     ,INTENT(IN)   :: npass                !Number of filter passes.

  INTEGER                   :: ii , jj , iii , jjj  , ipass
  INTEGER                   :: imin , imax , jmin , jmax , n_num
  REAL(r_sngl)              :: n_mean
 
 cfield = field 
 
 missing_maskc = missing_mask

 DO ipass = 1 , npass

 !$OMP PARALLEL DO PRIVATE( ii , jj , iii , jjj , n_num , n_mean , imin , imax , jmin , jmax )

  DO ii = 1 , nx
    DO jj = 1 , ny

       IF( missing_maskc(ii,jj) ) THEN
         !We have a missing value
         imin=max(1,ii-1)
         imax=min(nx,ii+1)
         jmin=max(1,jj-1)
         jmax=min(ny,jj+1)
       
         n_num=0
         n_mean=0.0e0 
         DO iii=imin,imax
           DO jjj=jmin,jmax
              IF( .not. missing_maskc(iii,jjj) )THEN
                n_num = n_num + 1
                n_mean = n_mean + cfield(iii,jjj)
              ENDIF
           ENDDO
         ENDDO 
         IF( n_num > 0 )THEN
           cfield(ii,jj)=n_mean/n_num
           missing_maskc(ii,jj)=.false.
         ENDIF

       ENDIF

    ENDDO
  ENDDO

  !$OMP END PARALLEL DO

 ENDDO
  


END SUBROUTINE com_complete_missing_2d

!2D interpolation using box average. Destination grid is assumed to be regular.
SUBROUTINE com_interp_boxavereg(xini,dx,nx,yini,dy,ny,xin,yin,datain,nin    &
               &                ,data_sum,data_max,data_min,data_n,undef)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
  INTEGER , INTENT(IN)        :: nx , ny , nin
  REAL(r_sngl),INTENT(IN)     :: dx , dy , xini , yini
  REAL(r_sngl),INTENT(IN)     :: undef
  REAL(r_sngl),INTENT(IN)     :: xin(nin),yin(nin),datain(nin)
  REAL(r_sngl),INTENT(INOUT)  :: data_sum(nx,ny)
  REAL(r_sngl),INTENT(INOUT)  :: data_max(nx,ny)
  REAL(r_sngl),INTENT(INOUT)  :: data_min(nx,ny)
  INTEGER     ,INTENT(INOUT)  :: data_n(nx,ny)
  INTEGER                     :: ii , ix , iy

 !I can not use OPEN MP on this loop since several threads may
 !need to access to the same memory location at the same time.
 !An OMP CRITICAL can be included but that probably will take away
 !all the advantage of paralelization in the first place.

  DO ii = 1,nin  !Loop over the input data 

    !Compute the location of the current point with in grid coordinates (rx,ry)
    ix = int( ( xin(ii) - xini ) / dx ) + 1
    iy = int( ( yin(ii) - yini ) / dy ) + 1

    !Check is the data is within the grid.
    IF( ix <= nx .and. ix >= 1 .and. iy <= ny .and. iy >= 1 .and. datain(ii) /= undef )THEN

      IF(  data_n(ix,iy) == 0 )THEN
        data_max(ix,iy) = datain(ii)
        data_min(ix,iy) = datain(ii)
        data_sum(ix,iy) = datain(ii)
        
      ELSE
        data_sum(ix,iy) = data_sum(ix,iy) + datain(ii)

        IF( datain(ii) > data_max(ix,iy) )THEN
          data_max(ix,iy) = datain(ii)
        ENDIF
        IF( datain(ii) < data_min(ix,iy) )THEN
          data_min(ix,iy) = datain(ii)
        ENDIF  
 
      ENDIF

       data_n(ix,iy) = data_n(ix,iy) + 1

    ENDIF

  ENDDO

END SUBROUTINE com_interp_boxavereg

!-----------------------------------------------------------------------
! Cubic spline interpolation
!   [Reference:] Akima, H., 1970: J. ACM, 17, 589-602.
!-----------------------------------------------------------------------
SUBROUTINE com_interp_spline(ndim,x,y,n,x5,y5)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
  INTEGER,INTENT(IN) :: ndim         ! number of grid points
  REAL(r_size),INTENT(IN) :: x(ndim) ! coordinate
  REAL(r_size),INTENT(IN) :: y(ndim) ! variable
  INTEGER,INTENT(IN) :: n            ! number of targets
  REAL(r_size),INTENT(IN) :: x5(n)   ! target coordinates
  REAL(r_size),INTENT(OUT) :: y5(n)  ! target values
  INTEGER :: i,j,m
  REAL(r_size) :: dydx(5),ddydx(4),t(2),dx21,dx
  REAL(r_size) :: wk

  TGT: DO j=1,n
    DO i=1,ndim
      IF(x5(j) == x(i)) THEN
        y5(j) = y(i)
        CYCLE TGT
      END IF
      IF(x5(j) < x(i)) EXIT
    END DO
!       i-3   i-2   i-1    i    i+1   i+2
!     ---+-----+-----+---*-+-----+-----+---
!dydx       1     2     3     4     5
!ddydx         1     2     3     4
!t                   1     2
    IF(i==2) THEN
      DO m=3,5
        dydx(m) = (y(i-3+m)-y(i-4+m)) / (x(i-3+m)-x(i-4+m))
      END DO
      dydx(2) = 2.0d0*dydx(3) - dydx(4)
      dydx(1) = 2.0d0*dydx(2) - dydx(3)
    ELSE IF(i==3) THEN
      DO m=2,5
        dydx(m) = (y(i-3+m)-y(i-4+m)) / (x(i-3+m)-x(i-4+m))
      END DO
      dydx(1) = 2.0d0*dydx(2) - dydx(3)
    ELSE IF(i==ndim) THEN
      DO m=1,3
        dydx(m) = (y(i-3+m)-y(i-4+m)) / (x(i-3+m)-x(i-4+m))
      END DO
      dydx(4) = 2.0d0*dydx(3) - dydx(2)
      dydx(5) = 2.0d0*dydx(4) - dydx(3)
    ELSE IF(i==ndim-1) THEN
      DO m=1,4
        dydx(m) = (y(i-3+m)-y(i-4+m)) / (x(i-3+m)-x(i-4+m))
      END DO
      dydx(5) = 2.0d0*dydx(4) - dydx(3)
    ELSE
      DO m=1,5
        dydx(m) = (y(i-3+m)-y(i-4+m)) / (x(i-3+m)-x(i-4+m))
      END DO
    END IF
    DO m=1,4
      ddydx(m) = ABS(dydx(m+1) - dydx(m))
    END DO
    DO m=1,2
      wk = ddydx(m+2) + ddydx(m)
      IF(wk == 0) THEN
        t(m) = 0.0d0
      ELSE
        t(m) = (ddydx(m+2)*dydx(m+1)+ddydx(m)*dydx(m+2))/wk
      END IF
    END DO
    dx21 = x(i)-x(i-1)
    dx = x5(j) - x(i-1)
    y5(j) = y(i-1) &
        & + dx*t(1) &
        & + dx*dx*(3.0d0*dydx(3)-2.0d0*t(1)-t(2))/dx21 &
        & + dx*dx*dx*(t(1)+t(2)-2.0d0*dydx(3))/dx21/dx21
  END DO TGT

  RETURN
END SUBROUTINE com_interp_spline
!-----------------------------------------------------------------------
! (LON,LAT) --> (i,j) conversion
!   [ORIGINAL AUTHOR:] Masaru Kunii
!-----------------------------------------------------------------------
SUBROUTINE com_pos2ij(msw,nx,ny,flon,flat,num_obs,olon,olat,oi,oj)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
  ! --- inout variables
  INTEGER,INTENT(IN) :: msw   !MODE SWITCH: 1: fast, 2: accurate
  INTEGER,INTENT(IN) :: nx,ny !number of grid points
  REAL(r_size),INTENT(IN) :: flon(nx,ny),flat(nx,ny) !(lon,lat) at (i,j)
  INTEGER,INTENT(IN) :: num_obs !repetition number of conversion
  REAL(r_size),INTENT(IN) :: olon(num_obs),olat(num_obs) !target (lon,lat)
  REAL(r_size),INTENT(OUT) :: oi(num_obs),oj(num_obs) !target (i,j)
  ! --- local work variables
  LOGICAL,PARAMETER :: detailout = .FALSE.
  INTEGER,PARAMETER :: num_grid_ave = 4  ! fix
  INTEGER :: inum,ix,jy,ip,wk_maxp
  INTEGER :: iorder_we,iorder_sn
  INTEGER :: nxp,nyp
  REAL(r_size),PARAMETER :: miss = -32768 
  REAL(r_size),PARAMETER :: max_dist = 2.0e+6
  REAL(r_size) :: rlat_max, rlat_min, rlon_max, rlon_min   
  REAL(r_size) :: dist(num_grid_ave) 
  REAL(r_size) :: dist_min_x(num_obs, num_grid_ave)
  REAL(r_size) :: dist_min_y(num_obs, num_grid_ave) 
  REAL(r_size) :: wk_dist, sum_dist
  REAL(r_size) :: ratio(num_grid_ave)
  IF(detailout) THEN
    WRITE(6,'(A)') '====================================================='
    WRITE(6,'(A)') '      Detailed output of SUBROUTINE com_pos2ij       '
    WRITE(6,'(A)') '====================================================='    
  END IF
  ! ================================================================
  !   Check the Order of flon, flat
  ! ================================================================   
  iorder_we = 1
  iorder_sn = 1
  IF(flon(1,1) > flon(2,1)) THEN
    iorder_we = -1
  END IF
  IF(flat(1,1) > flat(1,2)) THEN
    iorder_sn = -1
  END IF
  IF(detailout) THEN  
    WRITE(6,'(3X,A,I5)') 'Obs Order (WE) :',iorder_we 
    WRITE(6,'(3X,A,I5)') 'Obs Order (SN) :',iorder_sn 
  END IF
  ! ================================================================
  !  FAST MODE
  ! ================================================================   
  IF(msw == 1) THEN
    ! ==============================================================
    !   Surrounding 4 Grid Points Interpolation
    ! ==============================================================   
    Obs_Loop_1 : DO inum=1,num_obs 
      IF(detailout) WRITE(6,'(A,I5,2F15.5)') '*** START OBS ',inum,olon(inum),olat(inum) 
      ! ------------------------------------------------------------
      !    Search Basic Point
      ! ------------------------------------------------------------ 
      nxp = miss
      nyp = miss
      DO jy=1,ny-1
        DO ix=1,nx-1
          rlon_max = MAXVAL(flon(ix:ix+1, jy:jy+1))
          rlon_min = MINVAL(flon(ix:ix+1, jy:jy+1))
          rlat_max = MAXVAL(flat(ix:ix+1, jy:jy+1))
          rlat_min = MINVAL(flat(ix:ix+1, jy:jy+1))
          IF(rlon_min <= olon(inum) .AND. rlon_max >= olon(inum) .AND. &
           & rlat_min <= olat(inum) .AND. rlat_max >= olat(inum)) THEN
            nxp = ix
            nyp = jy
            EXIT
          END IF
        END DO
      END DO
      IF(detailout) WRITE(6,'(3X,A,2I7)') 'nxp, nyp =',nxp,nyp
      IF(nxp == miss .OR. nyp == miss) THEN
        WRITE(6,'(A)') '!!WARNING(com_pos2ij): obs position cannot be detected'
        oi(inum) = miss
        oj(inum) = miss
        CYCLE Obs_Loop_1
      END IF
      ! ------------------------------------------------------------
      !    Interpolation
      ! ------------------------------------------------------------    
      CALL com_distll_1(flon(nxp  ,nyp  ),flat(nxp  ,nyp  ),&
                      & olon(inum),olat(inum),dist(1))
      CALL com_distll_1(flon(nxp+1,nyp  ),flat(nxp+1,nyp  ),&
                      & olon(inum),olat(inum),dist(2))
      CALL com_distll_1(flon(nxp  ,nyp+1),flat(nxp  ,nyp+1),&
                      & olon(inum),olat(inum),dist(3))      
      CALL com_distll_1(flon(nxp+1,nyp+1),flat(nxp+1,nyp+1),&
                      & olon(inum),olat(inum),dist(4))      
      dist(1:4) = dist(1:4) * 1.D-3  
      IF(detailout) WRITE(6,'(3X,A,4F15.5)') 'distance :',dist(1:4) 
      sum_dist = dist(1) * dist(1) * dist(2) * dist(2) * dist(3) * dist(3) &
             & + dist(2) * dist(2) * dist(3) * dist(3) * dist(4) * dist(4) &
             & + dist(3) * dist(3) * dist(4) * dist(4) * dist(1) * dist(1) &
             & + dist(4) * dist(4) * dist(1) * dist(1) * dist(2) * dist(2)
      ratio(1) = (dist(2)*dist(2)*dist(3)*dist(3)*dist(4)*dist(4))/sum_dist
      ratio(2) = (dist(3)*dist(3)*dist(4)*dist(4)*dist(1)*dist(1))/sum_dist
      ratio(3) = (dist(4)*dist(4)*dist(1)*dist(1)*dist(2)*dist(2))/sum_dist
      ratio(4) = (dist(1)*dist(1)*dist(2)*dist(2)*dist(3)*dist(3))/sum_dist
      IF(detailout) WRITE(6,'(3X,A,5F15.5)') 'ratio    :',ratio(1:4),SUM(ratio(1:4))
      oi(inum) = ratio(1) *  nxp    + ratio(2) * (nxp+1) &
             & + ratio(3) *  nxp    + ratio(4) * (nxp+1)
      oj(inum) = ratio(1) *  nyp    + ratio(2) *  nyp    &
             & + ratio(3) * (nyp+1) + ratio(4) * (nyp+1)    
      IF(detailout) WRITE(6,'(3X,A,2F15.5)') 'position :',oi(inum), oj(inum)
 
    END DO Obs_Loop_1
  ! ================================================================
  !  ACCURATE MODE
  ! ================================================================   
  ELSE IF(msw == 2) THEN
    ! ================================================================
    !   Nearest 4 Grid Points Interpolation
    ! ================================================================   
    Obs_Loop_2 : DO inum=1,num_obs
      IF(detailout) WRITE(6,'(A,I5,2F15.5)') '*** START OBS ',inum,olon(inum),olat(inum) 
      ! ------------------------------------------------------------
      !    Search 4-Grid Points
      ! ------------------------------------------------------------      
      dist(1:num_grid_ave) = 1.D+10
      wk_maxp = num_grid_ave    
      DO jy=1,ny
        DO ix=1,nx
          CALL com_distll_1(flon(ix,jy),flat(ix,jy),&
                          & olon(inum) ,olat(inum) ,wk_dist)
          IF(wk_dist > max_dist) CYCLE
          IF(wk_dist < dist(wk_maxp)) THEN
            dist(wk_maxp) = wk_dist
            dist_min_x(inum, wk_maxp) = ix
            dist_min_y(inum, wk_maxp) = jy
            DO ip = 1, num_grid_ave
              IF(dist(ip) == maxval(dist(1:num_grid_ave))) THEN
                wk_maxp = ip
                EXIT
              END IF
            END DO
          END IF
        END DO
      END DO
      IF(detailout) WRITE(6,'(A,4(A,I4,A,I4,A))')  '  Intp Grids : ', &
        & '(', INT(dist_min_x(inum, 1)), ',', INT(dist_min_y(inum, 1)), ') ', &
        & '(', INT(dist_min_x(inum, 2)), ',', INT(dist_min_y(inum, 2)), ') ', &
        & '(', INT(dist_min_x(inum, 3)), ',', INT(dist_min_y(inum, 3)), ') ', &
        & '(', INT(dist_min_x(inum, 4)), ',', INT(dist_min_y(inum, 4)), ') '
      ! ------------------------------------------------------------
      !    Interpolation
      ! ------------------------------------------------------------ 
      dist(1:num_grid_ave) =  dist(1:num_grid_ave) * 1.0D-3
      sum_dist = dist(1) * dist(1) * dist(2) * dist(2) * dist(3) * dist(3)  &
             & + dist(2) * dist(2) * dist(3) * dist(3) * dist(4) * dist(4)  &
             & + dist(3) * dist(3) * dist(4) * dist(4) * dist(1) * dist(1)  &
             & + dist(4) * dist(4) * dist(1) * dist(1) * dist(2) * dist(2)
      ratio(1) = (dist(2)*dist(2)*dist(3)*dist(3)*dist(4)*dist(4))/sum_dist
      ratio(2) = (dist(3)*dist(3)*dist(4)*dist(4)*dist(1)*dist(1))/sum_dist
      ratio(3) = (dist(4)*dist(4)*dist(1)*dist(1)*dist(2)*dist(2))/sum_dist
      ratio(4) = (dist(1)*dist(1)*dist(2)*dist(2)*dist(3)*dist(3))/sum_dist
      IF(detailout) WRITE(6,'(2X,A,5F15.5)') 'ratio      :',ratio(1:4),SUM(ratio(1:4))
      oi(inum) = SUM(ratio(1:num_grid_ave) * dist_min_x(inum, 1:num_grid_ave))
      oj(inum) = SUM(ratio(1:num_grid_ave) * dist_min_y(inum, 1:num_grid_ave))
      IF(detailout) WRITE(6,'(2X,A,2F15.5)') 'position   :',oi(inum),oj(inum)
    END DO Obs_Loop_2
  END IF

  RETURN
END SUBROUTINE com_pos2ij
!-----------------------------------------------------------------------
! UTC to TAI93
!-----------------------------------------------------------------------
SUBROUTINE com_utc2tai(iy,im,id,ih,imin,sec,tai93)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
  INTEGER,INTENT(IN) :: iy,im,id,ih,imin
  REAL(r_size),INTENT(IN) :: sec
  REAL(r_size),INTENT(OUT) :: tai93
  REAL(r_size),PARAMETER :: mins = 60.0d0
  REAL(r_size),PARAMETER :: hour = 60.0d0*mins
  REAL(r_size),PARAMETER :: day = 24.0d0*hour
  REAL(r_size),PARAMETER :: year = 365.0d0*day
  INTEGER,PARAMETER :: mdays(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  INTEGER :: days,i

  tai93 = REAL(iy-1993,r_size)*year + FLOOR(REAL(iy-1993)/4.0,r_size)*day
  days = id -1
  DO i=1,12
    IF(im > i) days = days + mdays(i)
  END DO
  IF(MOD(iy,4) == 0 .AND. im > 2) days = days + 1 !leap year
  tai93 = tai93 + REAL(days,r_size)*day + REAL(ih,r_size)*hour &
              & + REAL(imin,r_size)*mins + sec
  IF(iy > 1993 .OR. (iy==1993 .AND. im > 6)) tai93 = tai93 + 1.0d0 !leap second
  IF(iy > 1994 .OR. (iy==1994 .AND. im > 6)) tai93 = tai93 + 1.0d0 !leap second
  IF(iy > 1995) tai93 = tai93 + 1.0d0 !leap second
  IF(iy > 1997 .OR. (iy==1997 .AND. im > 6)) tai93 = tai93 + 1.0d0 !leap second
  IF(iy > 1998) tai93 = tai93 + 1.0d0 !leap second
  IF(iy > 2005) tai93 = tai93 + 1.0d0 !leap second
  IF(iy > 2008) tai93 = tai93 + 1.0d0 !leap second
  IF(iy > 2012 .OR. (iy==2012 .AND. im > 6)) tai93 = tai93 + 1.0d0 !leap second

  RETURN
END SUBROUTINE com_utc2tai
!-----------------------------------------------------------------------
! TAI93 to UTC
!-----------------------------------------------------------------------
SUBROUTINE com_tai2utc(tai93,iy,im,id,ih,imin,sec)
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)

  INTEGER,PARAMETER :: n=8 ! number of leap seconds after Jan. 1, 1993
  INTEGER,PARAMETER :: leapsec(n) = (/  15638399,  47174400,  94608001,&
                                  &    141868802, 189302403, 410227204,&
                                  &    504921605, 615254406/)
  REAL(r_size),INTENT(IN) :: tai93
  INTEGER,INTENT(OUT) :: iy,im,id,ih,imin
  REAL(r_size),INTENT(OUT) :: sec
  REAL(r_size),PARAMETER :: mins = 60.0d0
  REAL(r_size),PARAMETER :: hour = 60.0d0*mins
  REAL(r_size),PARAMETER :: day = 24.0d0*hour
  REAL(r_size),PARAMETER :: year = 365.0d0*day
  INTEGER,PARAMETER :: mdays(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  REAL(r_size) :: wk,tai
  INTEGER :: days,i,leap

  tai = tai93
  sec = 0.0d0
  DO i=1,n
    IF(FLOOR(tai93) == leapsec(i)+1) sec = 60.0d0 + tai93-FLOOR(tai93,r_size)
    IF(FLOOR(tai93) > leapsec(i)) tai = tai -1.0d0
  END DO
  iy = 1993 + FLOOR(tai /year)
  wk = tai - REAL(iy-1993,r_size)*year - FLOOR(REAL(iy-1993)/4.0,r_size)*day
  IF(wk < 0.0d0) THEN
    iy = iy -1
    wk = tai - REAL(iy-1993,r_size)*year - FLOOR(REAL(iy-1993)/4.0,r_size)*day
  END IF
  days = FLOOR(wk/day)
  wk = wk - REAL(days,r_size)*day
  im = 1
  DO i=1,12
    leap = 0
    IF(im == 2 .AND. MOD(iy,4)==0) leap=1
    IF(im == i .AND. days >= mdays(i)+leap) THEN
      im = im + 1
      days = days - mdays(i)-leap
    END IF
  END DO
  id = days +1

  ih = FLOOR(wk/hour)
  wk = wk - REAL(ih,r_size)*hour
  imin = FLOOR(wk/mins)
  IF(sec < 60.0d0) sec = wk - REAL(imin,r_size)*mins

  RETURN
END SUBROUTINE com_tai2utc

!Get the missing value mask.
SUBROUTINE getmask(my_ensemble,my_undefmask,nx,ny,nz,nbv,undef)
IMPLICIT NONE
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
REAL(r_sngl) , INTENT(IN)  :: undef
REAL(r_sngl) , INTENT(IN)  :: my_ensemble(nx,ny,nz,nbv) 
LOGICAL      , INTENT(OUT) :: my_undefmask(nx,ny,nz)
INTEGER     ,INTENT(IN)    :: nx,ny,nz,nbv
INTEGER                    :: ii,jj,kk,im

my_undefmask=.false.

DO  ii = 1,nx
  DO jj = 1,ny
    DO kk = 1,nz
      DO im = 1,nbv
         if( abs(my_ensemble(ii,jj,kk,im)) == undef )then
           my_undefmask(ii,jj,kk)=.true.
           exit
         endif
      ENDDO 
    ENDDO
  ENDDO
ENDDO


END SUBROUTINE getmask

!compute different pdf metrics.
SUBROUTINE compute_pdf_metric(my_ensemble,my_undefmask,nx,ny,nz,nbv,nbins,pdf_metrics,undef)
IMPLICIT NONE  
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
REAL(r_size),PARAMETER :: pi=3.141592653589793d0
REAL(r_size),PARAMETER :: deg2rad=3.1415926535d0/180d0
REAL(r_size),PARAMETER :: rad2deg=180d0/3.1415926535d0
INTEGER, INTENT(IN)       :: nbins
INTEGER, INTENT(IN)       :: nx,ny,nz,nbv !Ensemble dimensions.
REAL(r_sngl) , INTENT(IN) :: my_ensemble(nx,ny,nz,nbv)
LOGICAL      , INTENT(IN) :: my_undefmask(nx,ny,nz)
!REAL(r_sngl), INTENT(IN)  :: ensemble(nx,ny,nz,nbv) !Ensemble data.
REAL(r_sngl), INTENT(IN)  :: undef        !Missing data code.
!LOGICAL, INTENT(IN)       :: undefmask(nx,ny,nz) !Valid grids.
REAL(r_sngl), INTENT(OUT) :: pdf_metrics(nx,ny,nz,7) !1-KL divergence
                                                     !2-Hellinger distance
                                                     !3-Total variation distance
                                                     !4-Renyi Divergence
                                                     !5-Jensen Shannon Divergence
                                                     !6-Bhattacharyy distance
                                                     !7-Energy distance.
REAL(r_sngl)              :: tmp(nbv)
INTEGER                   :: ii , jj , kk , ib 
INTEGER                   :: histogram(nx,ny,nz,nbins)
REAL(r_sngl)              :: varmin(nx,ny,nz) , varmax(nx,ny,nz)
REAL(r_sngl)              :: moments(nx,ny,nz,2)

!Double precission will be used for computation of the metrics.
REAL(r_size)              :: p , q , dvar , gconst , gconst2 , x
REAL(r_size)              :: cp , cq , m  , alpha , logp , logq , logm

pdf_metrics=0.0e0
alpha=2.0d0  !Parameter in the computation of the Renyi divergence

 !Compute the histogram.
 !write(*,*)"Computing histogram"
 CALL compute_histogram(my_ensemble,my_undefmask,nx,ny,nz,nbv,nbins,undef,varmin,varmax,histogram) 

 !write(*,*)"Computing moments"
 !Compute the mean and standard deviation.
 CALL compute_moments(my_ensemble,my_undefmask,nx,ny,nz,nbv,2,moments,undef)

 !write(*,*)"Computing kl distance"
 !$OMP PARALLEL DO PRIVATE(ii,jj,kk,ib,x,dvar,p,q,cp,cq,m,gconst,gconst2,logp,logq,logm)
  DO ii=1,nx
    DO jj=1,ny
     DO kk=1,nz
      if( .not. my_undefmask(ii,jj,kk) )then
       dvar=real(varmax(ii,jj,kk)-varmin(ii,jj,kk),r_size)/real(nbins,r_size) 
       gconst=sqrt(1.d0/real(pi*2.0d0*moments(ii,jj,kk,2),r_size))
       gconst2=(1.d0/real(2.0d0*moments(ii,jj,kk,2),r_size))

       cp=0.0d0
       cq=0.0d0
       DO ib=1,nbins
         x=varmin(ii,jj,kk)+(ib-0.5)*dvar
         p=real(histogram(ii,jj,kk,ib),r_size)/real(nbv,r_size)                         !Ensemble pdf
         q=dvar*gconst*exp(-((x-real(moments(ii,jj,kk,1),r_size))**2)*gconst2)          !Gaussian pdf.
         m=0.5*(p + q )                                                         !Mean distribution
         !Compute cdf
         cp = cp + p
         cq = cq + q
         logp=log(p)
         logq=log(q)
         logm=log(m)

         !Accumulation for Kullback Liebler Divergence.
         if( p > 0 )then
           !TODO when q(ib) is 0 replace by the machine epsilon
           pdf_metrics(ii,jj,kk,1)=pdf_metrics(ii,jj,kk,1)+p*( logp - logq )
         endif
         !Accumulation for Hellinger distance
         pdf_metrics(ii,jj,kk,2)=pdf_metrics(ii,jj,kk,2)+ ( sqrt(p)-sqrt(q) )**2 
         !Accumulation for Total Variation Distance
         if( abs( p - q ) > pdf_metrics(ii,jj,kk,3) ) then
           pdf_metrics(ii,jj,kk,3) = abs( p - q )
         endif
         !Accumulation for Renyi divergence
         if ( p > 0 )then
           pdf_metrics(ii,jj,kk,4) = pdf_metrics(ii,jj,kk,4) + (p**alpha)/(q ** (alpha-1))
         endif
         !Accumulation for Jensen-Shannon divergence
         if( p > 0 )then
           pdf_metrics(ii,jj,kk,5) = pdf_metrics(ii,jj,kk,5) + 0.5*( p*( logp - logm ) )
         endif
         if( q > 0 )then
           pdf_metrics(ii,jj,kk,5) = pdf_metrics(ii,jj,kk,5) + 0.5*( q*( logq - logm ) )
         endif 
         !Accumulation for Bhattacharyya distance
         pdf_metrics(ii,jj,kk,6) = pdf_metrics(ii,jj,kk,6) + sqrt( p * q )
         !Accumulation of Energy Distance
         pdf_metrics(ii,jj,kk,7) = pdf_metrics(ii,jj,kk,7) + ( cp - cq )**2
       ENDDO

       !Normalization of Hellinger distance
       pdf_metrics(ii,jj,kk,2) = sqrt( pdf_metrics(ii,jj,kk,2) )/sqrt( 2.0d0 )
       !Normalization of Renyi divergence
       pdf_metrics(ii,jj,kk,4) = (1.0d0/(alpha-1.0d0))*log( pdf_metrics(ii,jj,kk,4) )
       !Normalization of Jensen-Shannon divergence
       pdf_metrics(ii,jj,kk,5) = pdf_metrics(ii,jj,kk,5) / ( log( 10.0d0 ) / log( 2.0d0 ) ) 
       !Normalization of Bhattacharyya distance
       pdf_metrics(ii,jj,kk,6) = -log( pdf_metrics(ii,jj,kk,6) )
       !Normalization of Energy Distance
       pdf_metrics(ii,jj,kk,7) = sqrt( pdf_metrics(ii,jj,kk,7) )

      else
       pdf_metrics(ii,jj,kk,:)=undef
      endif
     ENDDO
    ENDDO
  ENDDO
 !$OMP END PARALLEL DO

END SUBROUTINE compute_pdf_metric

!compute the kld distance
SUBROUTINE compute_kld(my_ensemble,my_undefmask,nx,ny,nz,nbv,nbins,kld,undef)
IMPLICIT NONE
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
REAL(r_size),PARAMETER :: pi=3.141592653589793d0
INTEGER, INTENT(IN)       :: nbins
INTEGER, INTENT(IN)       :: nx,ny,nz,nbv !Ensemble dimensions.
REAL(r_sngl) , INTENT(IN) :: my_ensemble(nx,ny,nz,nbv)
LOGICAL      , INTENT(IN) :: my_undefmask(nx,ny,nz)
!REAL(r_sngl), INTENT(IN)  :: ensemble(nx,ny,nz,nbv) !Ensemble data.
REAL(r_sngl), INTENT(IN)  :: undef        !Missing data code.
!LOGICAL, INTENT(IN)       :: undefmask(nx,ny,nz) !Valid grids.
REAL(r_sngl), INTENT(OUT) :: kld(nx,ny,nz)
REAL(r_sngl)              :: tmp(nbv)
INTEGER                   :: ii , jj , kk , ib
INTEGER                   :: histogram(nx,ny,nz,nbins)
REAL(r_sngl)              :: varmin(nx,ny,nz) , varmax(nx,ny,nz)
REAL(r_sngl)              :: moments(nx,ny,nz,2)

REAL(r_sngl)              :: p(nbins) , q(nbins) , dvar , gconst , gconst2 , x

kld=undef

 !Compute the histogram.
 !write(*,*)"Computing histogram"
 CALL compute_histogram(my_ensemble,my_undefmask,nx,ny,nz,nbv,nbins,undef,varmin,varmax,histogram)

 !write(*,*)"Computing moments"
 !Compute the mean and standard deviation.
 CALL compute_moments(my_ensemble,my_undefmask,nx,ny,nz,nbv,2,moments,undef)

 !write(*,*)"Computing kl distance"
!$OMP PARALLEL DO PRIVATE(ii,jj,kk,ib,x,dvar,p,q,gconst,gconst2)
  DO ii=1,nx
    DO jj=1,ny
     DO kk=1,nz
      if( .not. my_undefmask(ii,jj,kk) )then
       dvar=( varmax(ii,jj,kk)-varmin(ii,jj,kk) )/real(nbins,r_sngl)
       gconst=sqrt(1./(pi*2*moments(ii,jj,kk,2)))
       gconst2=(1./(2*moments(ii,jj,kk,2)))
       kld(ii,jj,kk)=0.0e0
       DO ib=1,nbins
         x=varmin(ii,jj,kk)+(ib-0.5)*dvar
         p(ib)=real(histogram(ii,jj,kk,ib),r_sngl)/real(nbv,r_sngl)   !Ensemble pdf
         q(ib)=dvar*gconst*exp(-((x-moments(ii,jj,kk,1))**2)*gconst2) !Gaussian pdf.

         if( p(ib) > 0 .and. q(ib) > 0 )then
         !if( p(ib) > 0 ) then
           !write(*,*)p(ib),q(ib),log(p(ib)),log(q(ib))
           !TODO when q(ib) is 0 replace by the machine epsilon
           kld(ii,jj,kk)=kld(ii,jj,kk)+p(ib)*( log(p(ib)) - log(q(ib) ) )
         endif
       ENDDO
       !stop
       !if( kld(ii,jj,kk) < 0 )kld(ii,jj,kk)=0.0e0
      !else
      ! kld(ii,jj,kk)=undef
      endif
     ENDDO
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE compute_kld

SUBROUTINE compute_kld_sub( pdf_1 , pdf_2 , nbin , kld ) 

IMPLICIT NONE
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
INTEGER      , INTENT(IN) :: nbin                       !Number of bins for pdf discretization.
REAL(r_sngl) , INTENT(IN) :: pdf_1(nbin) , pdf_2(nbin)  !Input pdf to be compared.
REAL(r_sngl) , INTENT(OUT):: kld
INTEGER                   :: ib

 kld=0.0e0
 DO ib=1,nbin
    if( pdf_1(ib) > 0 .and. pdf_2(ib) > 0 )then
      !TODO when q(ib) is 0 replace by the machine epsilon
      kld=kld+pdf_1(ib)*( log(pdf_1(ib)) - log(pdf_2(ib) ) )
    endif
 ENDDO

END SUBROUTINE compute_kld_sub 


!compute the firtst N moments of the PDF centered around the mean.
SUBROUTINE compute_moments(my_ensemble,my_undefmask,nx,ny,nz,nbv,nmoments,moments,undef)
IMPLICIT NONE  
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
INTEGER, INTENT(IN)       :: nx,ny,nz,nbv !Ensemble dimensions.
REAL(r_sngl) , INTENT(IN) :: my_ensemble(nx,ny,nz,nbv)
LOGICAL      , INTENT(IN) :: my_undefmask(nx,ny,nz)
!REAL(r_sngl), INTENT(IN)  :: ensemble(nx,ny,nz,nbv) !Ensemble data.
REAL(r_sngl), INTENT(IN)  :: undef        !Missing data code.
INTEGER, INTENT(IN)       :: nmoments     !Number of moments to be computed.
!LOGICAL, INTENT(IN)       :: undefmask(nx,ny,nz) !Valid grids.
REAL(r_sngl), INTENT(OUT) :: moments(nx,ny,nz,nmoments)
REAL(r_sngl)              :: tmp(nbv)
INTEGER                   :: ii , jj , kk , im
REAL(r_sngl)              :: w(nbv)

w=1.0e0


moments=0.0e0

!$OMP PARALLEL DO PRIVATE(ii,jj,kk,im,tmp)
  DO ii=1,nx
   DO im=1,nmoments
    DO jj=1,ny
     DO kk=1,nz
      if( .not. my_undefmask(ii,jj,kk) )then
       tmp=( ( my_ensemble(ii,jj,kk,:)-moments(ii,jj,kk,1) ) )**im
       CALL com_mean_sngl(nbv,tmp,moments(ii,jj,kk,im),w)
      else
       moments(ii,jj,kk,im)=undef
      endif
     ENDDO
    ENDDO
   ENDDO
  ENDDO
!$OMP END PARALLEL DO



END SUBROUTINE compute_moments

!compute the sample linear correlation coefficient between variable 1 and variable 2 at each grid point.
SUBROUTINE compute_correlation(my_ensemble1,my_ensemble2,my_undefmask1,my_undefmask2  &
     &                          ,nx,ny,nz,nbv,correlation,undef)
IMPLICIT NONE
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
INTEGER, INTENT(IN)       :: nx,ny,nz,nbv !Ensemble dimensions.
REAL(r_sngl) , INTENT(IN) :: my_ensemble1(nx,ny,nz,nbv) , my_ensemble2(nx,ny,nz,nbv)
LOGICAL      , INTENT(IN) :: my_undefmask1(nx,ny,nz) , my_undefmask2(nx,ny,nz)
REAL(r_sngl), INTENT(IN)  :: undef        !Missing data code.
REAL(r_sngl), INTENT(OUT) :: correlation(nx,ny,nz)
INTEGER                   :: ii , jj , kk 
REAL(r_sngl)              :: w(nbv)

w=1.0e0

correlation=0.0e0

!$OMP PARALLEL DO PRIVATE(ii,jj,kk)
  DO ii=1,nx
   DO jj=1,ny
    DO kk=1,nz
     if( ( .not. my_undefmask1(ii,jj,kk) ) .and. ( .not. my_undefmask2(ii,jj,kk) )  )then
       CALL com_correl_sngl(nbv,my_ensemble1(ii,jj,kk,:),my_ensemble2(ii,jj,kk,:),correlation(ii,jj,kk),w)
      else
       correlation(ii,jj,kk)=undef
     endif
    ENDDO
   ENDDO
  ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE compute_correlation

!Compute the histogram.
SUBROUTINE compute_histogram(my_ensemble,my_undefmask,nx,ny,nz,nbv,nbins,undef,varmin,varmax,histogram)
IMPLICIT NONE
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
INTEGER, INTENT(IN)       :: nx,ny,nz,nbv           !Ensemble dimensions.
REAL(r_sngl) , INTENT(IN) :: my_ensemble(nx,ny,nz,nbv)
LOGICAL      , INTENT(IN) :: my_undefmask(nx,ny,nz)
INTEGER, INTENT(IN)       :: nbins                  !Number of moments to be computed.
REAL(r_sngl), INTENT(IN)  :: undef
INTEGER  , INTENT(OUT)    :: histogram(nx,ny,nz,nbins) 
INTEGER                   :: tmpindex , reclength , iunit
REAL(r_sngl),INTENT(OUT)  :: varmin(nx,ny,nz),varmax(nx,ny,nz)
REAL(r_sngl)              :: tmp(nbv) , dvar
INTEGER                   :: ii , jj , kk , im , irec

histogram=0

!$OMP PARALLEL DO PRIVATE(ii,jj,kk)
  DO ii=1,nx
   DO jj=1,ny
    DO kk=1,nz
       varmin(ii,jj,kk)=MINVAL( my_ensemble(ii,jj,kk,:) )
       varmax(ii,jj,kk)=MAXVAL( my_ensemble(ii,jj,kk,:) )
    ENDDO
   ENDDO
  ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(ii,jj,kk,im,tmp,dvar,tmpindex)
  DO ii=1,nx
   DO jj=1,ny
    DO kk=1,nz
      if( .not. my_undefmask(ii,jj,kk) ) then

       dvar=( varmax(ii,jj,kk)-varmin(ii,jj,kk) )/real(nbins,r_sngl)
  
       IF( dvar > 0 )THEN
         DO im=1,nbv
           tmpindex=floor( ( my_ensemble(ii,jj,kk,im) - varmin(ii,jj,kk) )/dvar ) + 1
           IF( tmpindex > nbins )tmpindex=nbins
           IF( tmpindex < 1     )tmpindex=1
           !WRITE(*,*)histogram(ii,jj,kk,tmpindex),tmpindex
           histogram(ii,jj,kk,tmpindex)=histogram(ii,jj,kk,tmpindex) + 1
         ENDDO
       ENDIF
      endif
    ENDDO
   ENDDO
  ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE compute_histogram

!Compute the histogram but using weights and assuming that the output is real. 
SUBROUTINE compute_histogram_sub(ens , nbv , nbins , varmin , varmax , histogram , w )
IMPLICIT NONE
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
INTEGER, INTENT(IN)       :: nbv           !Ensemble dimensions.
REAL(r_sngl) , INTENT(IN) :: ens(nbv)      !Ensemble size.
INTEGER, INTENT(IN)       :: nbins         !Number of bins that will be used.
REAL(r_sngl), INTENT(OUT) :: histogram(nbins)
INTEGER                   :: tmpindex 
REAL(r_sngl),INTENT(IN)   :: varmin,varmax
REAL(r_sngl),INTENT(IN)   :: w(nbv)  
REAL(r_sngl)              :: dvar
INTEGER                   :: im 

histogram=0.0

dvar=( varmax-varmin )/real(nbins,r_sngl)

IF( dvar > 0 )THEN
 DO im=1,nbv
   tmpindex=floor( ( ens(im) - varmin )/dvar ) + 1
   IF( tmpindex > nbins )tmpindex=nbins
   IF( tmpindex < 1     )tmpindex=1
   histogram(tmpindex)=histogram(tmpindex) + w(im)
 ENDDO
ENDIF


END SUBROUTINE compute_histogram_sub

!Compute the histogram bimodality index for a histograms computed from a 3D
!array
SUBROUTINE compute_histogram_bimodality(histogram,my_undefmask,nx,ny,nz,nbins,smoothw,bmindex)
IMPLICIT NONE
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
INTEGER, INTENT(IN)       :: nx,ny,nz,nbins           !Dimensions.
INTEGER, INTENT(IN)       :: smoothw                  !Smoothing filter length for the histogram.
LOGICAL      , INTENT(IN) :: my_undefmask(nx,ny,nz)
INTEGER      , INTENT(IN) :: histogram(nx,ny,nz,nbins)
INTEGER                   :: ii , jj , kk 
REAL(r_sngl) , INTENT(OUT):: bmindex(nx,ny,nz)        !Bimodality index for each grid point.

bmindex=0.0e0

!$OMP PARALLEL DO PRIVATE(ii,jj,kk)
DO ii = 1 , nx
  DO jj = 1 , ny
    DO kk = 1 , nz
      IF( .NOT. my_undefmask(ii,jj,kk) )THEN
         CALL  histogram_bimodality_index(histogram(ii,jj,kk,:),nbins,smoothw,bmindex(ii,jj,kk))
      ENDIF
    ENDDO
  ENDDO
ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE compute_histogram_bimodality

!Compute the histogram bimodality index for a histograms computed from a 3D
!array
SUBROUTINE compute_histogram_bimodality_improved(histogram,my_undefmask,nx,ny,nz,nbins,smoothw,bmindex)
IMPLICIT NONE
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
INTEGER, INTENT(IN)       :: nx,ny,nz,nbins           !Dimensions.
INTEGER, INTENT(IN)       :: smoothw                  !Smoothing filter length for the histogram.
LOGICAL      , INTENT(IN) :: my_undefmask(nx,ny,nz)
INTEGER      , INTENT(IN) :: histogram(nx,ny,nz,nbins)
INTEGER                   :: ii , jj , kk
REAL(r_sngl) , INTENT(OUT):: bmindex(nx,ny,nz)        !Bimodality index for each grid point.

bmindex=0.0e0

!$OMP PARALLEL DO PRIVATE(ii,jj,kk)
DO ii = 1 , nx
  DO jj = 1 , ny
    DO kk = 1 , nz
      IF( .NOT. my_undefmask(ii,jj,kk) )THEN
         CALL histogram_bimodality_index_improved(histogram(ii,jj,kk,:),nbins,smoothw,bmindex(ii,jj,kk))
      ENDIF
    ENDDO
  ENDDO
ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE compute_histogram_bimodality_improved

SUBROUTINE histogram_bimodality_index(histogramin,nb,smoothw,bmindex)
IMPLICIT NONE
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
INTEGER    ,   INTENT(IN) :: histogramin(nb) !histogram values and center location of each bin.
INTEGER      ,   INTENT(IN) :: nb ! number of bins
INTEGER      ,   INTENT(IN) :: smoothw  !Length of the smoothing filter for the histogram.
REAL(r_sngl) ,   INTENT(OUT):: bmindex  !Bimodality index.
REAL(r_sngl)                :: mean , std !Mean and std of the sample that produced the histogram.
REAL(r_sngl)                ::  histogram(nb)
INTEGER                     ::  ii , imin , imax , im 
REAL(r_sngl)                ::  hist_max(nb)                !P value at the local maximum location.
INTEGER                     ::  hist_max_loc(nb) , n_max    !Max location and number.
REAL(r_sngl)                ::  in_between_min(nb) , global_max_dist(nb) , max_index(nb)              
REAL(r_sngl)                ::  global_max                  !Global max P value
INTEGER                     ::  global_max_loc              !Global max location
REAL(r_sngl)                ::  hmean , hstd                !Mean and std estimated from the histogram

bmindex=0.0e0
!histogram=REAL(histogramin,r_sngl)

IF ( SUM(histogramin) == 0 )THEN
   !If the histogram in 0 (which happens for reflectivity when all the members
   !have the same value)
   RETURN
ENDIF


!Step 1 : Smooth the histogram (optional)

IF( smoothw > 0 )THEN !We will smooth the histogram.

  DO ii = 1 , nb
     !Box window smoothing.
     imin=max(1,ii-smoothw)
     imax=min(nb,ii+smoothw)

     histogram(ii) = SUM(histogramin( imin:imax ))/REAL( imax-imin+1,r_sngl)

  ENDDO

ELSE

  histogram = REAL( histogramin , r_sngl )

ENDIF

!Step 2 : Normalize the frequency and compute the approximated mean and spread.

histogram = histogram / sum(histogram) 

hmean=0.0e0
hstd =0.0e0

DO ii = 1 , nb
   hmean = hmean + ii * histogram(ii)
   hstd  = hstd  + (ii**2) * histogram(ii)
ENDDO

hstd = sqrt( hstd - mean**2 )

!Step 3 : Find the total number of local maxima.

n_max=0

DO ii = 2 , nb - 1

   IF( ii == 2 )THEN
     IF( histogram(1) > histogram(2) )THEN
       n_max=n_max+1
       hist_max(n_max) = histogram(1)
       hist_max_loc(n_max) = 1
     ENDIF
   ENDIF
   IF( ii == nb - 1 )THEN
     IF( histogram(nb) > histogram(nb-1) )THEN
       n_max=n_max+1
       hist_max(n_max) = histogram(nb)
       hist_max_loc(n_max) = nb
     ENDIF
   ENDIF

   IF( histogram(ii) >= histogram(ii+1) .and. histogram(ii) >= histogram(ii-1) )THEN
     n_max=n_max+1
     hist_max(n_max) = histogram(ii) 
     hist_max_loc(n_max) = ii
   ENDIF
ENDDO

!If there is only one minimum the return a value of 0.
IF( n_max  == 1 )THEN
  !There is only one maximum which is the global one. Then the index is 0.
  RETURN
ENDIF

!Step 4 : Identify the global maximum location and value.
global_max=0.0
DO im = 1 , n_max
   IF( hist_max( im ) > global_max )THEN
     global_max = hist_max(im)
     global_max_loc = hist_max_loc(im)
   ENDIF
ENDDO

!Step 5 : For each local maximum measure the "in between minimum strength" and
!the distance with respect to the local maximum.

in_between_min=0.0e0
global_max_dist=0.0e0
max_index=0.0e0

DO im = 1 , n_max

   imin=min(global_max_loc,hist_max_loc(im) )
   imax=max(global_max_loc,hist_max_loc(im) )

   !Get the minimum P between the current local minima and 
   !the global minima.  
   in_between_min(im)=minval( histogram(imin:imax) )

   !Compute the distance between the current local minima and the
   !global minima normalized by the histogram standard deviation.
   global_max_dist(im) = abs( imax - imin ) / hstd

   !Now compute the index for each local minima.
   !The index is the normalized distance times the minimum probability
   !difference between any of the maximums (the global and the local) and the
   !minimum probability in between.
   max_index(im) = min( abs( global_max - in_between_min(im) ) , abs( hist_max(im) - in_between_min(im) ) )

   max_index(im) = max_index(im) * global_max_dist(im)

ENDDO

!Step 6 : Given the value of the bimodality index computed for all the local
!minima present in the histogram, return the maximum of all these index values.

bmindex = maxval( max_index(1:n_max) ) 

!if( bmindex < 0 )then
!      write(*,*)bmindex
!      write(*,*)histogram
!      write(*,*)hmean,hstd,n_max
!      STOP
!endif

RETURN 

END SUBROUTINE histogram_bimodality_index



SUBROUTINE histogram_bimodality_index_improved(histogramin,nb,smoothw,bmindex)
IMPLICIT NONE
!This improved version computes the total probability mass that is in the basin
!within two local maxima.
!This formulation is more robust to sampling noise than the previous one.
!It is also easier to justify since the index is based on the "fluid" that can
!be contained by the basin or valley
!created in between two relative probability modes. If one of the mode is too
!small, the the resulting basin will be shallow 
!and will contain less fluid than two strong maxima which are well separated by
!a deep valley in the histogram. 
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
INTEGER    ,   INTENT(IN) :: histogramin(nb) !histogram values and center location of each bin.
INTEGER      ,   INTENT(IN) :: nb ! number of bins
INTEGER      ,   INTENT(IN) :: smoothw  !Length of the smoothing filter for the histogram.
REAL(r_sngl) ,   INTENT(OUT):: bmindex  !Bimodality index.
REAL(r_sngl)                :: mean , std !Mean and std of the sample that produced the histogram.
REAL(r_sngl)                ::  histogram(nb)
INTEGER                     ::  ii , imin , imax , im , ic 
REAL(r_sngl)                ::  hist_max(nb)                !P value at the local maximum location.
INTEGER                     ::  hist_max_loc(nb) , n_max    !Max location and number.
REAL(r_sngl)                ::  valley_volume(nb) , valley_volume_left , valley_volume_right 
LOGICAL                     ::  cont              

bmindex=0.0e0
!histogram=REAL(histogramin,r_sngl)

IF ( SUM(histogramin) == 0 )THEN
   !If the histogram in 0 (which happens for reflectivity when all the members
   !have the same value)
   RETURN
ENDIF


!Step 1 : Smooth the histogram (optional)

IF( smoothw > 0 )THEN !We will smooth the histogram.

  DO ii = 1 , nb
     !Box window smoothing.
     imin=max(1,ii-smoothw)
     imax=min(nb,ii+smoothw)

     histogram(ii) = SUM(histogramin( imin:imax ))/REAL( imax-imin+1,r_sngl)

  ENDDO

ELSE

  histogram = REAL( histogramin , r_sngl )

ENDIF

!Step 2 : Normalize the frequency.

histogram = histogram / sum(histogram) 

!Step 3 : Find the total number of local maxima.

n_max=0

DO ii = 2 , nb - 1

   IF( ii == 2 )THEN
     IF( histogram(1) > histogram(2) )THEN
       n_max=n_max+1
       hist_max(n_max) = histogram(1)
       hist_max_loc(n_max) = 1
     ENDIF
   ENDIF
   IF( ii == nb - 1 )THEN
     IF( histogram(nb) > histogram(nb-1) )THEN
       n_max=n_max+1
       hist_max(n_max) = histogram(nb)
       hist_max_loc(n_max) = nb
     ENDIF
   ENDIF

   IF( histogram(ii) >= histogram(ii+1) .and. histogram(ii) >= histogram(ii-1) )THEN
     n_max=n_max+1
     hist_max(n_max) = histogram(ii) 
     hist_max_loc(n_max) = ii
   ENDIF
ENDDO

!If there is only one minimum the return a value of 0.
IF( n_max  == 1 )THEN
  !There is only one maximum which is the global one. Then the index is 0.
  RETURN
ENDIF

!Step 4 : For each local maximum measure the "in between minimum strength" and
!the distance with respect to the local maximum.

valley_volume = 0.0e0

DO im = 1 , n_max

   valley_volume_right=0.0e0
   valley_volume_left =0.0e0

   !Compute the volume of the basin to the right of the maximum
   ic = hist_max_loc(im) + 1
   cont=.true.
   do while ( cont )
      if ( histogram(ic) < hist_max(im) ) then
         valley_volume_right = valley_volume_right + ( hist_max(im) - histogram(ic) )
      else
         cont=.false.
      endif
       
      if ( ic >= nb .and. cont )  then
         !We reached the end of the histogram and we could not find the other
         !side of the basin.
         valley_volume_right = 0.0e0
         cont = .false. 
      endif
      ic = ic + 1
  
   enddo 

   !Compute the volume of the basin to the left of the maximum
   ic = hist_max_loc(im) -1
   cont=.true.
   do while ( cont )
      if ( histogram(ic) < hist_max(im) ) then
         valley_volume_left = valley_volume_left + ( hist_max(im) - histogram(ic) )
      else
         cont=.false.
      endif
       
      if ( ic <= 1 .and. cont )  then
         !We reached the end of the histogram and we could not find the other
         !side of the basin.
         valley_volume_left = 0.0e0
         cont = .false. 
      endif
      ic = ic - 1
  
   enddo 
      
   valley_volume(im) = max( valley_volume_right , valley_volume_left )

ENDDO

!Step 5 : Given the value of the bimodality index computed for all the local
!minima present in the histogram, return the maximum of all these index values.

!The bmindex is defined as the volume corresponding to the largest basin.
!Another possible formulation is to use the sum of the volumes of all possible
!basins 
!this is ok because the different basin detected on each maximum are not
!overlaped, however
!taking the sum instead of the maximum might increase the impact of sampling
!noise (many small basin
!generated by sampling noise can add up to a volume similar to the one produced
!by one large basin
!associated to a truly bimodal distribution).
bmindex = maxval( valley_volume(1:n_max) ) 

RETURN 

END SUBROUTINE histogram_bimodality_index_improved



SUBROUTINE compare_update_sub( hxf , ens  , nens , yo , obs_error , mean_diff , mode_diff , std_diff  &
                             , kld , updated_mean_kf , updated_std_kf , updated_mode_kf , weights )  
IMPLICIT NONE
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
INTEGER     , INTENT(IN)  :: nens
REAL(r_sngl), INTENT(IN)  :: hxf(nens) , ens(nens)  !Ensemble in obs space (1D) and ensemble in model space (1D)
REAL(r_sngl), INTENT(IN)  :: yo                     !Observation value.
REAL(r_sngl), INTENT(IN)  :: obs_error              !Observation error.
REAL(r_sngl), INTENT(IN)  :: weights(nens)          !Input weights for pf update.
REAL(r_sngl), INTENT(OUT) :: mean_diff , mode_diff , std_diff , kld !Difference in mean , mode and var of the updated ensembles and kld between the two posterior distributions.
REAL(r_sngl) :: weights_kf(nens)
REAL(r_sngl) :: updated_ens_kf(nens) , updated_weights_pf(nens) 
REAL(r_sngl), INTENT(OUT) :: updated_mean_kf , updated_std_kf , updated_mode_kf
REAL(r_sngl) :: updated_mean_pf , updated_std_pf , updated_mode_pf

REAL(r_sngl) , ALLOCATABLE :: updated_hist_kf(:) , updated_hist_pf(:)
INTEGER      :: nbins 
REAL(r_sngl) :: var_min_hist , var_max_hist , varmin , varmax

   updated_weights_pf = weights/sum(weights)

   !Assume uniform weigths for the KF prior.
   weights_kf = 1.0
   

!Initialization
mean_diff=0.0e0
mode_diff=0.0e0
std_diff =0.0e0

CALL compute_enkf_update(hxf,ens,nens,yo,obs_error,updated_ens_kf)
CALL compute_pf_update(hxf,ens,nens,yo,obs_error,weights,updated_weights_pf)

CALL com_mean_sngl(nens,updated_ens_kf,updated_mean_kf,weights_kf)
CALL com_mean_sngl(nens,ens           ,updated_mean_pf,updated_weights_pf)
mean_diff = updated_mean_kf - updated_mean_pf 

CALL com_stdev_sngl(nens,updated_ens_kf,updated_std_kf,weights_kf)
CALL com_stdev_sngl(nens,ens           ,updated_std_pf,updated_weights_pf)
std_diff = updated_std_kf - updated_std_pf

CALL com_mode_sngl(nens,updated_ens_kf,updated_mode_kf,weights_kf,.false.)
CALL com_mode_sngl(nens,ens           ,updated_mode_pf,updated_weights_pf,.false.)
mode_diff = updated_mode_kf - updated_mode_pf 

nbins = INT( sqrt( REAL( nens , r_sngl ) ) )

ALLOCATE( updated_hist_kf( nbins ) , updated_hist_pf( nbins ) )

varmax = MAX( MAXVAL( ens ) , MAXVAL( updated_ens_kf ) )
varmin = MIN( MINVAL( ens ) , MINVAL( updated_ens_kf ) ) 


CALL compute_histogram_sub( updated_ens_kf , nens , nbins , varmin , varmax , updated_hist_kf , weights_kf ) 
CALL compute_histogram_sub( ens           , nens , nbins , varmin , varmax , updated_hist_pf ,  updated_weights_pf ) 

!write(*,*) updated_hist_kf , updated_hist_pf

updated_hist_kf = updated_hist_kf / sum( updated_hist_kf ) 
updated_hist_pf = updated_hist_pf / sum( updated_hist_pf )

CALL compute_kld_sub( updated_hist_kf , updated_hist_pf , nbins , kld )

!write(*,*) kld 

DEALLOCATE( updated_hist_kf , updated_hist_pf )

END SUBROUTINE compare_update_sub


SUBROUTINE compare_update( obs_var , state_var , undefmask , nens , nx , ny , nz , obs_increment &
                         , obs_error , mean_diff , mode_diff ,std_diff , kld  &
                         , updated_mean_kf , updated_std_kf , updated_mode_kf )

IMPLICIT NONE
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
INTEGER, INTENT(IN)       :: nx,ny,nz,nens               !Dimensions.
REAL(r_sngl), INTENT(IN) :: obs_var(nx,ny,nz,nens)       !Observed variable.
REAL(r_sngl), INTENT(IN) :: state_var(nx,ny,nz,nens)     !State variable to be updated.
LOGICAL     , INTENT(IN) :: undefmask(nx,ny,nz)
REAL(r_sngl), INTENT(IN) :: obs_increment                !Diff between obs_var ensemble mean and the observation.
REAL(r_sngl), INTENT(IN) :: obs_error                    !Observational error standard deviation.
REAL(r_sngl), INTENT(OUT) :: mean_diff(nx,ny,nz)         !Difference in posterior mean between pf and kf.
REAL(r_sngl), INTENT(OUT) :: mode_diff(nx,ny,nz)         !Difference in posterior mode between pf and kf.
REAL(r_sngl), INTENT(OUT) :: std_diff(nx,ny,nz)          !Difference in posterior std in pf and kf.
REAL(r_sngl), INTENT(OUT) :: kld(nx,ny,nz)               !Difference in kld in pf and kf.
REAL(r_sngl), INTENT(OUT) :: updated_mean_kf(nx,ny,nz)   !Updated ensemble mean by the EnKF
REAL(r_sngl), INTENT(OUT) :: updated_std_kf(nx,ny,nz)    !Updated ensemble spread by the EnKF
REAL(r_sngl), INTENT(OUT) :: updated_mode_kf(nx,ny,nz)   !Updated ensemble mode by the EnKF
INTEGER      :: ii , jj , kk , im
REAL(r_sngl) :: yo , w(nens)

w = 1.0
!$OMP PARALLEL DO PRIVATE(ii,jj,kk,yo)
DO ii = 1 , nx
  DO jj = 1 , ny
    DO kk = 1 , nz
      IF ( .NOT. undefmask(ii,jj,kk) )THEN
         CALL com_mean_sngl( nens , obs_var(ii,jj,kk,:) , yo , w )
         yo = yo + obs_increment
         CALL compare_update_sub( obs_var(ii,jj,kk,:) , state_var(ii,jj,kk,:) ,                       &
                               nens , yo , obs_error , mean_diff(ii,jj,kk) , mode_diff(ii,jj,kk) ,    &
                               std_diff(ii,jj,kk) , kld(ii,jj,kk) , updated_mean_kf(ii,jj,kk) ,       &
                               updated_std_kf(ii,jj,kk) , updated_mode_kf(ii,jj,kk) , w )
      ENDIF
    ENDDO
  ENDDO
ENDDO
!OMP END PARALLEL DO

END SUBROUTINE compare_update

SUBROUTINE compute_single_obs_update(hxf,ens,mask,nx,ny,nz,nens,yo,obs_error,undef,update_type, &
                    &     updated_mean , updated_std , updated_mode )
IMPLICIT NONE
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
REAL(r_sngl), INTENT(IN)  :: hxf(nens)   !Ensemble in observation space (1D)
REAL(r_sngl), INTENT(IN)  :: ens(nx,ny,nz,nens)   !Ensemble in model space 
REAL(r_sngl), INTENT(IN)  :: undef
LOGICAL     , INTENT(IN)  :: mask(nx,ny,nz)
INTEGER     , INTENT(IN)  :: nx,ny,nz,nens        !Dimensions
REAL(r_sngl), INTENT(IN)  :: yo          !Observed value.
REAL(r_sngl), INTENT(IN)  :: obs_error   !Observation error standard deviation
INTEGER     , INTENT(IN)  :: update_type !1-EnKF , 2-PF
REAL(r_sngl), INTENT(OUT) :: updated_mean(nx,ny,nz)  !Updated ensemble mean
REAL(r_sngl), INTENT(OUT) :: updated_std(nx,ny,nz)   !Updated ensemble var
REAL(r_sngl), INTENT(OUT) :: updated_mode(nx,ny,nz)  !Updated ensemble mode
REAL(r_sngl)              :: alpha
REAL(r_sngl)              :: w(nens)
INTEGER                   :: iens , ii , jj , kk 
REAL(r_sngl)              :: updated_ens(nens)

updated_mean = 0.0e0
updated_std  = 0.0e0
updated_mode = 0.0e0

w=1.0e0

!$OMP PARALLEL DO PRIVATE(ii,jj,kk,w,updated_ens)

DO ii = 1 , nx
  DO jj = 1 , ny
    DO kk = 1 , nz
       IF( .not. mask( ii , jj , kk ) )THEN
         IF ( update_type == 1 ) THEN
            CALL compute_enkf_update(hxf,ens(ii,jj,kk,:),nens,yo,obs_error,updated_ens)
         ENDIF 
         IF ( update_type == 2 ) THEN
            w=1.0d0
            updated_ens=ens(ii,jj,kk,:)
            CALL compute_pf_update(hxf,ens(ii,jj,kk,:),nens,yo,obs_error,w,w)
         ENDIF
         CALL com_mean_sngl(nens,updated_ens,updated_mean(ii,jj,kk),w)
         CALL com_stdev_sngl(nens,updated_ens,updated_std(ii,jj,kk),w)
         CALL com_mode_sngl(nens,updated_ens,updated_mode(ii,jj,kk),w,.false.)
       ELSE
          updated_mean( ii , jj , kk ) = undef
          updated_std( ii , jj , kk ) = undef
          updated_mode( ii , jj , kk ) = undef
       ENDIF
    ENDDO
  ENDDO
ENDDO

!$OMP END PARALLEL DO


END SUBROUTINE compute_single_obs_update


SUBROUTINE compute_enkf_update(hxf,ens,nens,yo,obs_error,updated_ens)
IMPLICIT NONE
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
REAL(r_sngl), INTENT(IN)  :: hxf(nens)   !Ensemble in observation space (1D)
REAL(r_sngl), INTENT(IN)  :: ens(nens)   !Ensemble in model space (1D)
INTEGER     , INTENT(IN)  :: nens        !Ensemble size.
REAL(r_sngl), INTENT(IN)  :: yo          !Observed value.
REAL(r_sngl), INTENT(IN)  :: obs_error   !Observation error standard deviation
REAL(r_sngl), INTENT(OUT) :: updated_ens(nens)  !Posterior ensemble in model space (1D)
REAL(r_sngl)              :: alpha
REAL(r_sngl)              :: hxf_mean , hxf_std , hxf_var 
REAL(r_sngl)              :: ens_mean , ens_std , ens_var , cov_ens
REAL(r_sngl)              :: updated_mean , w(nens) 
INTEGER                   :: iens

!Computhe the enkf update of ens (Following the ensemble square root filter of W&H 2002 MWR)
w=1.0e0

CALL com_mean_sngl(nens,hxf,hxf_mean,w)
CALL com_mean_sngl(nens,ens,ens_mean,w)
CALL com_stdev_sngl(nens,hxf,hxf_std,w)
CALL com_stdev_sngl(nens,ens,ens_std,w)
ens_var = ens_std ** 2
hxf_var = hxf_std ** 2
CALL com_covar_sngl(nens,hxf,ens,cov_ens,w)

alpha = 1.0 / ( 1.0 + sqrt( ( obs_error ** 2 ) / ( hxf_var + ( obs_error ** 2 ) ) ) )

updated_mean = ens_mean + ( cov_ens / ( hxf_var + ( obs_error ** 2 ) ) ) * ( yo - hxf_mean )



DO iens = 1 , nens
   updated_ens( iens ) = ( ens( iens ) - ens_mean ) -  &
   alpha * ( cov_ens / ( hxf_var + ( obs_error ** 2 ) ) ) * ( hxf( iens ) - hxf_mean )  + updated_mean
ENDDO

END SUBROUTINE compute_enkf_update

SUBROUTINE compute_pf_update(hxf,ens,nens,yo,obs_error,weights_in,weights_out)

IMPLICIT NONE
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
REAL(r_sngl), INTENT(IN)  :: hxf(nens)   !Ensemble in observation space (1D)
REAL(r_sngl), INTENT(IN)  :: ens(nens)   !Ensemble in model space (1D)
INTEGER     , INTENT(IN)  :: nens        !Ensemble size.
REAL(r_sngl), INTENT(IN)  :: yo          !Observed value.
REAL(r_sngl), INTENT(IN)  :: obs_error   !Observation error standard deviation
REAL(r_sngl), INTENT(IN)  :: weights_in(nens)   !Prior weights (1D)
REAL(r_sngl), INTENT(OUT) :: weights_out(nens)  !Posterior weights (1D)
REAL(r_sngl)              :: hxf_mean
INTEGER                   :: ii

do ii = 1 , nens 
   weights_out(ii) = weights_in(ii) * exp( -0.5 * ( ( hxf(ii) - yo )**2 ) / ( obs_error ** 2 )  )
enddo

!Normalize weights
weights_out = weights_out / sum( weights_out )

END SUBROUTINE compute_pf_update

SUBROUTINE read_ensemble(path,ensemble,undef_mask,nx,ny,nbv,selected_fields,n_selected_fields,undef,ie,acc)
IMPLICIT NONE
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
INTEGER , INTENT(IN)       :: selected_fields(n_selected_fields)
INTEGER , INTENT(IN)       :: nx , ny , nbv , n_selected_fields 
REAL(r_sngl) , INTENT(IN)  :: undef
REAL(r_sngl) , INTENT(OUT) :: ensemble(nx,ny,n_selected_fields,nbv)
LOGICAL      , INTENT(OUT) :: undef_mask(nx,ny,n_selected_fields)
REAL(r_sngl)               :: bufr(nx,ny)
CHARACTER(*) , INTENT(IN)  :: path , ie , acc                           !Data pathand input endian
INTEGER                    :: ifield , maxfield , i, ibv , iunit , reclength , ii , jj , isf , base_unit
INTEGER                    :: field_counter
CHARACTER(40)              :: filename 
CHARACTER(200)             :: filenamewithpath

INQUIRE(IOLENGTH=reclength)reclength
reclength=nx*ny*reclength

undef_mask=.false.
ensemble=0.0e0
base_unit=200

!filename='____.grd'

!$OMP PARALLEL DO PRIVATE(ibv,iunit,filename,filenamewithpath,ifield,bufr,ii,jj,maxfield,field_counter,isf)
DO ibv = 1,nbv

   iunit = base_unit + ibv

   filename='____.grd'
   write(filename(1:4),'(I4.4)')ibv
   WRITE(*,*)"Reading file ", filename
   filenamewithpath= path // filename

   IF( acc == 'direct' )THEN

     OPEN(iunit,FILE=filenamewithpath,FORM='unformatted',ACCESS='direct',RECL=reclength,CONVERT=ie)

     !Read only the selected fields.
     DO ifield = 1,n_selected_fields
        READ(iunit,rec=selected_fields(ifield)) bufr
        ensemble(:,:,ifield,ibv)=bufr
        DO ii=1,nx
         DO jj=1,ny
           if( ensemble(ii,jj,ifield,ibv)==undef)then
             undef_mask(ii,jj,ifield)=.true.
           endif
         ENDDO
        ENDDO
     ENDDO

   ENDIF
   IF( acc == 'sequential' )THEN

     OPEN(iunit,FILE=filenamewithpath,FORM='unformatted',ACCESS='sequential',RECL=reclength,CONVERT=ie)

     !Read only the selected fields.
     maxfield=maxval( selected_fields )
     field_counter = 1 
     DO ifield = 1,maxfield
        READ(iunit) bufr
        DO isf = 1 , n_selected_fields
           IF( ifield == selected_fields(isf) )THEN
               ensemble(:,:,field_counter,ibv)=bufr
               field_counter = field_counter + 1
               DO ii=1,nx
                 DO jj=1,ny
                   IF( ensemble(ii,jj,ifield,ibv)==undef)THEN
                     undef_mask(ii,jj,ifield)=.false.
                   ENDIF
                 ENDDO
               ENDDO
               EXIT
           ENDIF
        ENDDO        
     ENDDO
   ENDIF

ENDDO
!OMP END PARALLEL DO

END SUBROUTINE read_ensemble


SUBROUTINE write_data(outfile,mydata,nx,ny,nz,ie,acc)
IMPLICIT NONE
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
INTEGER , INTENT(IN)       :: nx , ny , nz
REAL(r_sngl) , INTENT(IN)  :: mydata(nx,ny,nz)
CHARACTER(*) , INTENT(IN)  :: outfile , ie , acc                           !Data
INTEGER                    :: iz , iunit , reclength 

INQUIRE(IOLENGTH=reclength)reclength
reclength=nx*ny*reclength

WRITE(*,*)"Writing a file ", outfile

IF( acc == 'direct' )THEN

  OPEN(iunit,FILE=outfile,FORM='unformatted',ACCESS='direct',RECL=reclength,CONVERT=ie)
  !Read only the selected fields.
  DO iz = 1,nz
     WRITE( iunit,rec=iz ) mydata(:,:,iz) 
  ENDDO

ENDIF

IF( acc == 'sequential' )THEN

  OPEN(iunit,FILE=outfile,FORM='unformatted',ACCESS='sequential',RECL=reclength,CONVERT=ie)
  !Read only the selected fields.
  DO iz = 1,nz
     WRITE(iunit) mydata(:,:,iz) 
  ENDDO
ENDIF

END SUBROUTINE write_data


!Perform a 2D smoothing of the ensemble using a Lanczos filter.
SUBROUTINE smooth_2d(mydata,myfdata,my_undefmask,nx,ny,nz,dx,xfs)
IMPLICIT NONE
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
INTEGER      , INTENT(IN)       :: nx,ny,nz               !Grid dimensions
LOGICAL      , INTENT(IN)       :: my_undefmask(nx,ny,nz)
REAL(r_sngl) , INTENT(IN)       :: mydata(nx,ny,nz)       !Model Data
REAL(r_sngl) , INTENT(OUT)      :: myfdata(nx,ny,nz)      !Filtered Model Data
REAL(r_sngl) , INTENT(IN)       :: dx      !Grid resolution.
REAL(r_sngl) , INTENT(IN)       :: xfs     !Filter scale in x and y.
!LOGICAL      , INTENT(IN)       :: undefmask(nx,ny,nz)
INTEGER                         :: kk
INTEGER                         :: integermask(nx,ny,nz)
INTEGER                         :: filter_size_x
REAL(r_sngl)                    :: tmpdata(nx,ny,nz)

integermask=0
WHERE( my_undefmask )
  integermask=1
END WHERE

filter_size_x = NINT( xfs / dx )     

tmpdata=REAL(mydata,r_size)

!$OMP PARALLEL DO PRIVATE(kk)
DO kk=1,nz  
    CALL lanczos_2d(tmpdata(:,:,kk),tmpdata(:,:,kk),integermask(:,:,kk),filter_size_x,nx,ny)
ENDDO
!$OMP END PARALLEL DO

myfdata=REAL(tmpdata,r_sngl)

END SUBROUTINE smooth_2d

SUBROUTINE lanczos_2d(inputvar2d,outputvar2d,mask,lambda,nx,ny)
!This is a single-routine lanczos filter. This routine has been prepared in
!order to be used
  IMPLICIT NONE
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
  INTEGER,INTENT(IN)        :: nx,ny 
  REAL(r_sngl) , INTENT(IN) :: inputvar2d(nx,ny) 
  REAL(r_sngl) , INTENT(OUT):: outputvar2d(nx,ny)
  INTEGER      , INTENT(IN) :: mask(nx,ny)
  INTEGER      , INTENT(IN) :: lambda
  REAL(r_sngl)                         :: fn, rspval        ! cutoff freq/wavelength, spval
  REAL(r_sngl), DIMENSION(0:2*lambda)  :: dec, de           ! weight in r8, starting index 0:nband
  INTEGER                              :: ji, jj, jmx, jkx, jk, jt, jvar ! dummy loop index
  INTEGER                              :: ik1x, ik2x, ikkx
  INTEGER                              :: ifrst=0
  INTEGER                              :: inxmin, inxmaxi
  INTEGER                              :: inymin, inymaxi
  REAL(r_sngl), DIMENSION(nx,ny)       :: dl_tmpx, dl_tmpy
  REAL(r_sngl)                         :: dl_yy, dl_den, dl_pi, dl_ey, dl_coef

  !   PRINT *,'        ncut     : number of grid step to be filtered'
  !  remark: for a spatial filter, fn=dx/lambda where dx is spatial step, lamda
  !  is cutting wavelength

  fn    = 1./lambda

  !Filter init -------------------------------------  
 
    dl_pi   = ACOS(-1.d0)
    dl_coef = 2.0d0*dl_pi*fn

    de(0) = 2.d0*fn
    DO  ji=1,2*lambda
       de(ji) = SIN(dl_coef*ji)/(dl_pi*ji)
    END DO
    !
    dec(0) = 2.d0*fn
    DO ji=1,2*lambda
       dl_ey   = dl_pi*ji/(2*lambda)
       dec(ji) = de(ji)*SIN(dl_ey)/dl_ey
    END DO

  !End of filter init --------------------------------
  
  !FILTER-----------------------------------------------

   IF ( lambda /= 0 )THEN

    inxmin   =  2*lambda
    inxmaxi  =  nx-2*lambda+1
    inymin   =  2*lambda
    inymaxi  =  ny-2*lambda+1


    DO jj=1,ny
       DO  jmx=1,nx
          ik1x = -2*lambda
          ik2x =  2*lambda
          !
          IF (jmx <= inxmin ) ik1x = 1-jmx
          IF (jmx >= inxmaxi) ik2x = nx-jmx
          !
          dl_yy  = 0.d0
          dl_den = 0.d0
          !
          DO jkx=ik1x,ik2x
             ikkx=ABS(jkx)
             IF (mask(jkx+jmx,jj)  ==  1) THEN
                dl_den = dl_den + dec(ikkx)
                dl_yy  = dl_yy  + dec(ikkx)*inputvar2d(jkx+jmx,jj)
             END IF
          END DO
          !
          dl_tmpx(jmx,jj)=dl_yy/dl_den
       END DO
    END DO

    DO ji=1,nx
       DO  jmx=1,ny
          ik1x = -2*lambda
          ik2x =  2*lambda
          !
          IF (jmx <= inymin ) ik1x = 1-jmx
          IF (jmx >= inymaxi) ik2x = ny-jmx
          !
          dl_yy  = 0.d0
          dl_den = 0.d0
          !
          DO jkx=ik1x,ik2x
             ikkx=ABS(jkx)
             IF (mask(ji,jkx+jmx)  ==  1) THEN
                dl_den = dl_den + dec(ikkx)
                dl_yy  = dl_yy  + dec(ikkx)*dl_tmpx(ji,jkx+jmx)
             END IF
          END DO
          outputvar2d(ji,jmx)=0.
          IF (dl_den /=  0.) outputvar2d(ji,jmx) = dl_yy/dl_den
       END DO
    END DO
    !


   ENDIF 

  !END OF FILTER-----------------------------------------


   IF ( lambda == 0 ) outputvar2d=inputvar2d

   WHERE( mask == 0 )outputvar2d = inputvar2d


  END SUBROUTINE lanczos_2d

SUBROUTINE GAUSSIAN_FILTER(field0,field1,dx,sigma,nx,ny,undef)

IMPLICIT NONE
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
INTEGER            :: ii , jj , iii , jjj
REAL(r_sngl), INTENT(IN) :: undef 
REAL(r_sngl), INTENT(IN) :: field0(nx,ny),dx,sigma
REAL(r_sngl), INTENT(OUT):: field1(nx,ny)
INTEGER , INTENT(IN) :: nx,ny
REAL(r_sngl)             :: tmp , suma , sumaw , dist

!WRITE(*,*) " INSIDE GAUSSIAN FILTER "
!WRITE(*,*) " DELTA X                  ",dx
!WRITE(*,*) " SIGMA                    ",sigma
!WRITE(*,*) " NX ,  NY                 ",nx,ny

!$OMP PARALLEL DEFAULT(SHARED)  PRIVATE(ii,jj,iii,jjj,dist,tmp,suma,sumaw) 
!$OMP DO 
DO ii = 1,nx
  DO jj = 1,ny
        !COMPUTE AVERAGED VALUE.
        suma=0.0d0
        sumaw=0.0d0
        DO iii=1,nx
          DO jjj=1,ny
               dist=( ((ii-iii)*dx) ** 2 + ((jj-jjj)*dx) ** 2 )/( sigma**2)
               IF( dist <= 10)THEN
                 tmp=exp(-dist)
                 suma=suma+field0(iii,jjj)*tmp
                 sumaw=sumaw+tmp
               ENDIF
          ENDDO
        ENDDO

  IF( sumaw > 0.0d0)THEN
   field1(ii,jj)=suma/sumaw
  ELSE
   field1(ii,jj)=UNDEF
  ENDIF

  ENDDO
ENDDO

!$END OMP DO
!$OMP END PARALLEL

!WRITE(*,*)"FINISH GAUSSIAN FILTER"
END SUBROUTINE GAUSSIAN_FILTER

SUBROUTINE cont_table( dfor , dobs , ndata , thresholds , nthresholds , undef ,  &
               ets , csi , bias , pod , far , ctable )
IMPLICIT NONE
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
INTEGER      , INTENT(IN)  :: ndata , nthresholds
REAL(r_size) , INTENT(IN)  :: dfor(ndata) , dobs(ndata)
REAL(r_size) , INTENT(IN)  :: thresholds(nthresholds)
REAL(r_size) , INTENT(IN)  :: undef 
INTEGER      , INTENT(OUT) :: ctable(2,2,nthresholds)
REAL(r_size) , INTENT(OUT) :: ets(nthresholds) , csi(nthresholds) , bias(nthresholds)
REAL(r_size) , INTENT(OUT) :: pod(nthresholds) , far(nthresholds)
INTEGER                    :: hi , mi , fa , cn , nt 
INTEGER                    :: it , ii
REAL(r_size)               :: hi_r

ctable=0
ets=undef
csi=undef
bias=undef
pod=undef
far=undef

DO it = 1 , nthresholds 
   hi=0
   mi=0
   fa=0
   cn=0
   nt=0

   !$OMP PARALLEL DO REDUCTION(+:hi,mi,fa,cn,nt)
   DO ii = 1,ndata
      IF( dfor(ii) /= undef .and. dobs(ii) /= undef )THEN
        nt = nt + 1
        IF( dfor(ii) > thresholds(it) .and. dobs(ii) > thresholds(it) )THEN
          hi = hi+1
        ENDIF  
        IF( dfor(ii) <= thresholds(it) .and. dobs(ii) <= thresholds(it) )THEN
          cn = cn+1
        ENDIF
        IF( dfor(ii) > thresholds(it) .and. dobs(ii) <= thresholds(it) )THEN
          fa = fa+1
        ENDIF
        IF( dfor(ii) <= thresholds(it) .and. dobs(ii) > thresholds(it) )THEN
          mi = mi+1
        ENDIF
      ENDIF
   ENDDO
   !$OMP END PARALLEL DO
   ctable(1,1,it)=hi
   ctable(1,2,it)=fa
   ctable(2,1,it)=mi
   ctable(2,2,it)=cn

   !Compute scores
   IF ( nt > 0 ) THEN
      bias(it) = REAL( hi + fa , r_sngl )/REAL( hi + mi , r_sngl )
      csi(it)  = REAL( hi , r_sngl )/REAL( hi + mi + fa , r_sngl )
      pod(it)  = REAL( hi , r_sngl )/REAL( hi + mi , r_sngl )
      far(it)  = REAL( fa , r_sngl )/REAL( hi + fa , r_sngl )
      hi_r = REAL( ( hi + mi ) * ( hi + fa ) , r_sngl )/REAL( nt , r_sngl )
      ets(it) = ( REAL( hi , r_sngl ) - hi_r )/( REAL( hi + mi + fa , r_sngl ) - hi_r )
   ENDIF

ENDDO

WRITE(*,*)'Contingency table computed, ',REAL(nt,r_sngl)/REAL(ndata),' % valida data'

END SUBROUTINE cont_table

SUBROUTINE fss_det_2d( dfor , dobs , nx , ny , undef , thresholds , nthresholds , sizes , nsizes , fss )
   !Computes the Fraction Skill Score ussuming 2D fields and deterministic forecasts.
   !Computation is done following Skok and Roberts (2016) QJRMS
   !Undef are allowed and a check is performed to detect them on each box.
   !If the proportion of undef within the box is larger than (1-per_threshold) then the box is rejected
   !and is not considered for the computation of the FSS.
   !This criteria is also apply to handle boxes which are close to the domain boundaries.  
   !
   IMPLICIT NONE
   INTEGER,PARAMETER   :: r_size=kind(0.0d0)
   INTEGER,PARAMETER   :: r_sngl=kind(0.0e0)
   INTEGER,INTENT(IN)  :: nx , ny , nthresholds , nsizes 
   REAL(r_size), INTENT(IN) :: dfor(nx,ny) , dobs(nx,ny) 
   REAL(r_size), INTENT(IN) :: thresholds(nthresholds) 
   INTEGER     , INTENT(IN) :: sizes(nsizes) 
   REAL(r_size), INTENT(IN) :: undef
   REAL(r_size), INTENT(OUT):: fss(nthresholds,nsizes)
   REAL(r_size)             :: tmp_dfor(nx,ny,1) , tmp_dobs(nx,ny,1) 
   REAL(r_size)             :: tmp_sdfor(nx,ny,1) , tmp_sdobs(nx,ny,1)
   REAL(r_size)             :: tmp_n(nx,ny,1)
   REAL(r_size)             :: fbs , fbs_ref
   INTEGER                  :: validn , is , it
   LOGICAL                  :: mask( nx , ny , 1)
   INTEGER                  :: maski(nx , ny , 1)
   REAL(r_size),PARAMETER   :: per_thresh = 0.8


   fss = undef 
   WHERE( tmp_dfor == undef )
           tmp_dobs = undef
   ENDWHERE
   WHERE( tmp_dobs == undef )
           tmp_dfor = undef 
   ENDWHERE

   tmp_dfor(:,:,1) = dfor   !box_functions assumes 3D inputs.
   tmp_dobs(:,:,1) = dobs   !box_functions assumes 3D inputs.

   DO it = 1 , nthresholds
     DO is = 1 , nsizes 
        CALL box_functions( tmp_dfor,nx,ny,1,undef,sizes(is),sizes(is),0,'COUN',thresholds(it),tmp_sdfor)
        CALL box_functions( tmp_dobs,nx,ny,1,undef,sizes(is),sizes(is),0,'COUN',thresholds(it),tmp_sdobs)
        CALL box_functions( tmp_dfor,nx,ny,1,undef,sizes(is),sizes(is),0,'VALI',0.0d0,tmp_n)
        tmp_sdfor = tmp_sdfor / REAL( tmp_n , r_size )
        tmp_sdobs = tmp_sdobs / REAL( tmp_n , r_size ) 
        tmp_n = tmp_n / REAL( (sizes(is)*2+1) ** 2 ,r_size)  !Compute percentage of valid data within each box.
        mask  = tmp_n >= per_thresh         !Define the mask of valid elements.
        maski = mask
        validn = SUM( maski )
        fbs = SUM( (tmp_sdfor - tmp_sdobs)**2 , mask )
        fbs_ref = SUM( (tmp_sdfor )**2 , mask ) + SUM( (tmp_sdobs)**2 , mask )
        fss(it,is) = 1.0d0 - ( fbs / fbs_ref )
     ENDDO
   ENDDO


END SUBROUTINE fss_det_2d

SUBROUTINE BOX_FUNCTIONS(datain,nx,ny,nz,undef,boxx,boxy,boxz,operation,threshold,dataout)

IMPLICIT NONE
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
INTEGER     ,INTENT(IN) :: nx,ny,nz    !Grid dimension
INTEGER     ,INTENT(IN) :: boxx,boxy,boxz  !Box dimension
REAL(r_size),INTENT(IN) :: datain(nx,ny,nz)
REAL(r_size),INTENT(IN) :: undef
CHARACTER(4),INTENT(IN) :: operation     
REAL(r_size),INTENT(IN) :: threshold
REAL(r_size),INTENT(OUT) :: dataout(nx,ny,nz) !Result
REAL(r_size),ALLOCATABLE :: tmp_field(:) 
REAL(r_size)             :: tmp_mean , tmp_var
INTEGER                  :: NITEMS
INTEGER                  :: ii , jj , kk , bii , bjj , bkk , box_size , iin ,ii_index , data_count

box_size=(2*boxx+1)*(2*boxy+1)*(2*boxz+1)

!DO kk=1,nz
!WRITE(*,*)maxval( datain(:,:,kk) )
!ENDDO

ALLOCATE( tmp_field(box_size) )

dataout=undef

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(kk,ii,jj,bkk,bii,bjj,ii_index,NITEMS,tmp_field,tmp_mean,tmp_var,data_count,iin)
DO kk=1,nz
 DO ii=1,nx
  DO jj=1,ny


     NITEMS=0
     tmp_field=0.0d0
      DO bkk=kk-boxz,kk+boxz
       DO bii=ii-boxx,ii+boxx
         DO bjj=jj-boxy,jj+boxy
            IF( bkk >= 1 .AND. bkk <= nz .AND. bjj >= 1 .AND. bjj <= ny )THEN
              !Boundary condition in X
              ii_index=bii
              IF( bii < 1 )ii_index=bii+nx
              IF( bii > nx)ii_index=bii-nx
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
      ELSEIF( OPERATION == 'VALI' )THEN
        !Count all valid data (i.e. those not equal to undef)
        dataout(ii,jj,kk) = REAL( NITEMS , r_size )

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
END SUBROUTINE BOX_FUNCTIONS





END MODULE common_functions

