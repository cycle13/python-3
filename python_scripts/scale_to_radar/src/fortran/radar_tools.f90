MODULE radar_tools
!=======================================================================
!
! [PURPOSE:] Common radar procedures and variables.
!
! [HISTORY:]
!   02/19/2014 Juan Ruiz created
!
!=======================================================================
!$USE OMP_LIB
  IMPLICIT NONE
  PUBLIC
  CONTAINS
!-----------------------------------------------------------------------
! Compute the radar grid in LAT , LON , Z
!-----------------------------------------------------------------------
SUBROUTINE radar_georeference( rrange,relev,razim,rlon0,rlat0,rz0,ne,nr,na,   &
                radar_z , radar_lon , radar_lat , radar_x , radar_y ,         & 
                radar_local_elevation , distance_to_radar )
IMPLICIT NONE
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)

INTEGER      , INTENT(IN) :: ne,nr,na
REAL(r_size) , INTENT(IN) :: rlon0 , rlat0 , rz0 
REAL(r_size) , INTENT(IN) :: rrange(nr) ,relev(ne), razim(na)
REAL(r_size) , INTENT(OUT)::  radar_z(ne,nr,na)
REAL(r_size) , INTENT(OUT)::  radar_lat(ne,nr,na)
REAL(r_size) , INTENT(OUT)::  radar_lon(ne,nr,na)
REAL(r_size) , INTENT(OUT)::  radar_x(ne,nr,na)
REAL(r_size) , INTENT(OUT)::  radar_y(ne,nr,na)
REAL(r_size) , INTENT(OUT)::  radar_local_elevation(ne,nr)
REAL(r_size) , INTENT(OUT)::  distance_to_radar(ne,nr)
REAL(r_size)  :: tmp1 , tmp2 
INTEGER       :: ia , ir , ie
REAL(r_size),PARAMETER :: pi=3.141592653589793d0
REAL(r_size),PARAMETER :: gg=9.81d0
REAL(r_size),PARAMETER :: rd=287.0d0
REAL(r_size),PARAMETER :: cp=7.0d0 / 2.0d0 * rd
REAL(r_size),PARAMETER :: re=6371.3d3
REAL(r_size),PARAMETER :: r_omega=7.292d-5
REAL(r_size),PARAMETER :: t0c=273.15d0
REAL(r_size),PARAMETER :: undef=9.99d33
REAL(r_sngl),PARAMETER :: undefs=9.99e33
REAL(r_size),PARAMETER :: deg2rad=3.1415926535d0/180d0
REAL(r_size),PARAMETER :: rad2deg=180d0/3.1415926535d0
REAL(r_size),PARAMETER :: ke=(4d0/3d0)

  !Perform standard height beam heigth computation.
  !$OMP PARALLEL DO PRIVATE(ie,ia,tmp1,tmp2)
  DO ir=1,nr
     DO ie=1,ne
      tmp1= rrange(ir)**2 + (ke*Re)**2 + 2*rrange(ir)*ke*Re*sin( relev(ie)*deg2rad )
      radar_z(ie,ir,:)=rz0 + SQRT( tmp1 )-ke*Re;
      !Compute the local elevation angle tacking into account the efective earth radius.
      !tmp1= rrange(ir) * cos( relev(ie)*deg2rad )
      !tmp2= rrange(ir) * sin( relev(ie)*deg2rad ) + ke*Re
      !local_elevation(ie,ir) = relev(ie) + atan(tmp1/tmp2)
      
     ENDDO

     !Compute the latitude and longitude corresponding to each radar point.
     DO ie=1,ne
       distance_to_radar(ie,ir)=ke*Re*asin( rrange(ir) * cos(relev(ie)*deg2rad) / (ke*Re) )
       
       DO ia=1,na
        CALL com_ll_arc_distance(rlon0,rlat0,               &
             distance_to_radar(ie,ir),razim(ia),radar_lon(ie,ir,ia), &
             radar_lat(ie,ir,ia))
        radar_x(ie,ir,ia) = distance_to_radar(ie,ir) * sin(razim(ia)*deg2rad)
        radar_y(ie,ir,ia) = distance_to_radar(ie,ir) * cos(razim(ia)*deg2rad)

       ENDDO
     ENDDO
  ENDDO
  !$OMP END PARALLEL DO

END SUBROUTINE radar_georeference

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

SUBROUTINE com_ll_arc_distance(ini_lon,ini_lat,distance,azimuth,final_lon,final_lat)
IMPLICIT NONE
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
REAL(r_size), INTENT(IN)  :: ini_lon,ini_lat,distance,azimuth
REAL(r_size), INTENT(OUT) :: final_lon,final_lat
REAL(r_size)  :: cdist,sdist,sinll1,cosll1,cosaz,sinaz
REAL(r_size),PARAMETER :: deg2rad=3.1415926535d0/180d0
REAL(r_size),PARAMETER :: rad2deg=180d0/3.1415926535d0
REAL(r_size),PARAMETER :: re=6371.3d3

!So interface is single precission, but internal computations
!are performed with double precission.

 IF ( distance .EQ. 0e0 )THEN
    final_lon=ini_lon
    final_lat=ini_lat
 ELSE

 cdist = cos(distance/re)
 sdist = sin(distance/re)

 sinll1 = sin(ini_lat*deg2rad)
 cosll1 = cos(ini_lat*deg2rad)

 cosaz  = cos(azimuth*deg2rad)
 sinaz  = sin(azimuth*deg2rad)

 final_lat=asin(sinll1*cdist + cosll1*sdist*cosaz )*rad2deg
 final_lon=ini_lon + rad2deg *atan2(sdist*sinaz,cosll1*cdist - &
                  &  sinll1*sdist*cosaz)

 ENDIF

END SUBROUTINE com_ll_arc_distance


END MODULE radar_tools
