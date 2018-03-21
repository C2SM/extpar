!+  Fortran module for bilinear interpolation routines
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_8         2013/03/12 Frank Brenner
!  fixed bug in get_4_surrounding_raw_data_indices regarding negative longitudes
! V2_0_3       2015-01-12 Juergen Helmert
!  Bugfix correction covers CSCS SVN r5907-r6359
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module for bilinear interpolation routines
!> to calculate a bilinear interpolation from a regular raw data grid (geographical projection, 
!> rectangular grid) to a given point with geographical coordinates
!> \author Hermann Asensio
MODULE mo_bilinterpol


!> kind parameters are defined in MODULE data_parameters
USE mo_kind, ONLY: wp, &
                   i4, &
                   i8


!> abort_extpar defined in MODULE utilities_extpar
USE mo_utilities_extpar, ONLY: abort_extpar


IMPLICIT NONE

PRIVATE

PUBLIC :: get_4_surrounding_raw_data_indices, &
          calc_weight_bilinear_interpol, &
          calc_value_bilinear_interpol
    CONTAINS

       !> calculate the 4 surrounding raw data indeces for a given target point
       !!
       !! the definition of the regular lon-lat grid requires 
       !! - the coordinates of the north-western point of the domain ("upper left") startlon_reg_lonlat and startlat_reg_lonlat
       !! - the increment dlon_reg_lonlat and dlat_reg_lonlat(implict assuming that the 
       !! grid definiton goes from the west to the east and from the north to the south)
       !! - the number of grid elements nlon_reg_lonlat and nlat_reg_lonlat for both directions
       SUBROUTINE get_4_surrounding_raw_data_indices(   reg_data_grid, &
                                                        reg_lon,           &
                                                        reg_lat,           &
                                                        gldata,            &
                                                        point_lon_geo,      &
                                                        point_lat_geo,      &
                                                        western_column,     &
                                                        eastern_column,     &
                                                        northern_row,       &
                                                        southern_row)
       USE mo_grid_structures, ONLY : reg_lonlat_grid

       USE mo_search_ll_grid, ONLY : find_reg_lonlat_grid_element_index

       TYPE(reg_lonlat_grid), INTENT(IN)  :: reg_data_grid               !< structure with the definition of the raw data grid
       REAL (KIND=wp), INTENT(IN)         :: reg_lon(1:reg_data_grid%nlon_reg)  !< longitude of raw data in geographical system
       REAL (KIND=wp), INTENT(IN)          :: reg_lat(1:reg_data_grid%nlat_reg)  !< latitude of raw date in geographical system
       LOGICAL, INTENT(IN) :: gldata !< logical switch to indicate whether the raw data has global coverage or not

       REAL (KIND=wp), INTENT(IN) :: point_lon_geo       !< longitude coordinate in geographical system of input point 
       REAL (KIND=wp), INTENT(IN) :: point_lat_geo       !< latitude coordinate in geographical system of input point

       
      INTEGER (KIND=i8), INTENT(OUT) :: western_column     !< the index of the western_column of raw data 
      INTEGER (KIND=i8), INTENT(OUT) :: eastern_column     !< the index of the eastern_column of raw data 
      INTEGER (KIND=i8), INTENT(OUT) :: northern_row       !< the index of the northern_row of raw data 
      INTEGER (KIND=i8), INTENT(OUT) :: southern_row       !< the index of the southern_row of raw data 



       ! local variables
       INTEGER  (KIND=i8) :: point_lon_index !< longitude index of point for regular lon-lat grid
       INTEGER  (KIND=i8) :: point_lat_index !< latitude index of point for regular lon-lat grid

       INTEGER (KIND=i8) :: undefined_integer   !< value for undefined integer

       REAL (KIND=wp) :: dist_lon_m1, dist_lon_p1
       INTEGER  (KIND=i8) :: point_lon_index_m1,  point_lon_index_p1

       REAL (KIND=wp) :: dist_lat_m1, dist_lat_p1
       INTEGER  (KIND=i8) :: point_lat_index_m1,  point_lat_index_p1

       
       INTEGER  (KIND=i8) :: c_m1, c_p1
       INTEGER  (KIND=i8) :: r_p1, r_m1
       
       undefined_integer = 0 ! set undefined to zero
        
       point_lon_index = 0 ! initialise with zero
       point_lat_index = 0 ! initialise with zero
       
       ! find  grid element indices for the point
       CALL find_reg_lonlat_grid_element_index(point_lon_geo, &
         &                                      point_lat_geo, &
         &                                      reg_data_grid, &
         &                                      point_lon_index,&
         &                                      point_lat_index)


       IF ((point_lon_index == 0).OR.(point_lat_index == 0) )   THEN  ! point is out of data grid range

       ! !HA debug
       ! PRINT *,'get_4_surrounding_raw_data_indices, point_lon_geo: ', point_lon_geo
       ! PRINT *,'get_4_surrounding_raw_data_indices, eg_data_grid%start_lon_re: ', reg_data_grid%start_lon_reg
       ! PRINT *,'gldata: ',gldata
        
         IF (gldata) THEN
           point_lat_index = NINT(( point_lat_geo - reg_data_grid%start_lat_reg)/reg_data_grid%dlat_reg) + 1 
           point_lon_index = NINT( (point_lon_geo - reg_data_grid%start_lon_reg)/reg_data_grid%dlon_reg) + 1

           IF (point_lat_index < 1) THEN ! point out of range of regular lon-lat grid
             point_lat_index = undefined_integer
           ELSE IF (point_lat_index > reg_data_grid%nlat_reg) THEN ! point out of range of regular lon-lat grid
             point_lat_index = undefined_integer
           ENDIF

           IF (point_lon_index < 1 ) THEN ! point out of range of regular lon-lat grid
             point_lon_index = undefined_integer
           ELSE IF (point_lon_index > reg_data_grid%nlon_reg) THEN ! point out of range of regular lon-lat grid
             point_lon_index = undefined_integer
           ENDIF

           IF (point_lon_index == undefined_integer) THEN 
             IF ((point_lon_geo >=-360.0).AND.(point_lon_geo<=reg_data_grid%start_lon_reg)) THEN
                 point_lon_index = 1
             ELSEIF (( point_lon_geo <=360.0).AND.(point_lon_geo>=reg_data_grid%end_lon_reg)) THEN
                  point_lon_index = reg_data_grid%nlon_reg
             ELSE
                 western_column  = undefined_integer
                 eastern_column  = undefined_integer
                 northern_row    = undefined_integer
                 southern_row    = undefined_integer
                 RETURN   ! don't proceed in this subroutine
             ENDIF
           ENDIF

           IF (point_lat_index == undefined_integer) THEN
             IF ((point_lat_geo <=90.0).AND.(point_lat_geo>=reg_data_grid%start_lat_reg)) THEN
               point_lat_index = 1
             ELSEIF (( point_lat_geo >=-90.0).AND.(point_lat_geo<=reg_data_grid%end_lat_reg)) THEN
                point_lat_index = reg_data_grid%nlat_reg
             ELSE
               western_column  = undefined_integer
               eastern_column  = undefined_integer
               northern_row    = undefined_integer
               southern_row    = undefined_integer
               RETURN   ! don't proceed in this subroutine
             ENDIF  
           ENDIF

         ELSE ! gldata - no global data
           western_column  = undefined_integer
           eastern_column  = undefined_integer
           northern_row    = undefined_integer
           southern_row    = undefined_integer
           RETURN   ! don't proceed in this subroutine
         ENDIF !gldata
       ENDIF ! out of range

       ! get startindex of surrounding boxes to interpolate
       ! first longitude index
       point_lon_index_m1 = point_lon_index - 1
       point_lon_index_p1 = point_lon_index + 1

       IF ( point_lon_index_m1 <= 0) THEN ! point is at (western) boundary
         c_m1 = point_lon_index
         c_p1 = point_lon_index_p1
       ELSE IF ( point_lon_index_p1 > reg_data_grid%nlon_reg) THEN ! point is at (eastern) boundary 
!DWD MERGE  ELSE IF ( point_lon_index_p1 >= reg_data_grid%nlon_reg) THEN 
         c_m1 = point_lon_index_m1
         c_p1 = point_lon_index
       ELSE  
         dist_lon_m1 = abs(point_lon_geo - reg_lon(point_lon_index_m1))
         dist_lon_p1 = abs(point_lon_geo - reg_lon(point_lon_index_p1))
           IF (dist_lon_m1 < dist_lon_p1 ) THEN
              c_m1 = point_lon_index_m1
              c_p1 = point_lon_index
           ELSE
              c_m1 = point_lon_index
              c_p1 = point_lon_index_p1
           ENDIF
       ENDIF  

       ! check for grid direction (west-east) or (east-west)
       IF ((reg_lon(c_p1) - reg_lon(c_m1)) > 0. ) THEN ! grid in west-east direction
           western_column = c_m1
           eastern_column = c_p1
       ELSE ! this is very unlikly to happen, but you never know what kind of grids you get...
          western_column = c_p1
          eastern_column = c_m1
       ENDIF 

       !-----


       ! latitude index
       point_lat_index_m1 = point_lat_index - 1
       point_lat_index_p1 = point_lat_index + 1

       IF ( point_lat_index_m1 <= 0) THEN ! point is at (northern/southern) boundary
         r_m1 = point_lat_index
         r_p1 = point_lat_index_p1
       ELSE IF ( point_lat_index_p1 >= reg_data_grid%nlat_reg) THEN ! point is at (southern/northern) boundary 
         r_m1 = point_lat_index_m1
         r_p1 = point_lat_index
       ELSE  
         dist_lat_m1 = abs(point_lat_geo - reg_lat(point_lat_index_m1))
         dist_lat_p1 = abs(point_lat_geo - reg_lat(point_lat_index_p1))
           IF (dist_lat_m1 < dist_lat_p1 ) THEN
              r_m1 = point_lat_index_m1
              r_p1 = point_lat_index
           ELSE
              r_m1 = point_lat_index
              r_p1 = point_lat_index_p1
           ENDIF
       ENDIF  

       ! check for grid direction (north-south) or (south-north)
       IF ((reg_lat(r_p1) - reg_lat(r_m1)) > 0. ) THEN ! grid in south-north direction
           northern_row = r_p1
           southern_row = r_m1
       ELSE  ! grid in north-south direction
          northern_row = r_m1
          southern_row = r_p1
       ENDIF   

        

       END SUBROUTINE get_4_surrounding_raw_data_indices

    !> ELEMENTAL SUBROUTINE to calculate weights bwlon and bwlat for bilinear interploation with a regular lonlat grid as input
       ELEMENTAL SUBROUTINE calc_weight_bilinear_interpol(pixel_lon, pixel_lat,&
                                                        & west_lon, east_lon, north_lat, south_lat, &
                                                        & bwlon, bwlat)
       
       REAL (KIND=wp), INTENT(in) :: pixel_lon !< longitude coordinate in geographical system of input point 
       REAL (KIND=wp), INTENT(in) :: pixel_lat !< latitude coordinate in geographical system of input point 
       REAL (KIND=wp), INTENT(in) :: west_lon  !< longitude of western pixel
       REAL (KIND=wp), INTENT(in) :: east_lon  !< longitude of eastern pixel
       REAL (KIND=wp), INTENT(in) :: north_lat !< latitude of northern pixel
       REAL (KIND=wp), INTENT(in) :: south_lat !< latitude of southern pixel
       REAL (KIND=wp), INTENT(out):: bwlon     !< weight bwlon
       REAL (KIND=wp), INTENT(out):: bwlat     !< weight bwlat

       REAL (KIND=wp) :: pixel_lon2, east_lon2

       east_lon2 = east_lon
       pixel_lon2 = pixel_lon

       DO WHILE (ABS(east_lon2+360._wp-west_lon) < ABS(east_lon2-west_lon))
         east_lon2 = east_lon2 + 360._wp
       ENDDO
       
       IF (east_lon2 < west_lon) THEN   ! other than expected east_lon is smaller than west_lon
         DO WHILE (ABS(pixel_lon2+360._wp-east_lon2) < ABS(pixel_lon2-east_lon2))
           pixel_lon2 = pixel_lon2 + 360._wp
         ENDDO
         IF (pixel_lon2 < east_lon2) THEN
           bwlon = 0._wp
         ELSE IF (pixel_lon2 > west_lon) THEN
           bwlon = 1._wp
         ELSE                            ! pixel_lon2 lies between east_lon2 and west_lon
           bwlon = ABS((pixel_lon2 - west_lon)/(east_lon2 - west_lon))
       ENDIF
       ELSE IF (east_lon2 > west_lon) THEN
         DO WHILE (ABS(pixel_lon2+360._wp-west_lon) < ABS(pixel_lon2-west_lon))
           pixel_lon2 = pixel_lon2 + 360._wp
         ENDDO
         IF (pixel_lon2 > east_lon2) THEN
           bwlon = 0._wp
         ELSE IF (pixel_lon2 < west_lon) THEN
           bwlon = 1._wp
         ELSE                            ! pixel_lon2 lies between east_lon2 and west_lon
           bwlon = ABS((pixel_lon2 - west_lon)/(east_lon2 - west_lon))
         ENDIF
       ELSE
         bwlon = 0.5
       ENDIF

       bwlat = (pixel_lat - south_lat)/(north_lat - south_lat)

       END  SUBROUTINE calc_weight_bilinear_interpol

       !> ELEMENTAL SUBROUTINE to calculate the bilinerar interpolated value
       ELEMENTAL FUNCTION calc_value_bilinear_interpol(bwlon, bwlat, &
                                                      & data_value_sw, data_value_se, data_value_ne, data_value_nw) &
                                                      &  RESULT (target_value)
       REAL (KIND=wp), INTENT(IN) :: bwlon     !< weight bwlon
       REAL (KIND=wp), INTENT(IN) :: bwlat     !< weight bwlat
       REAL (KIND=wp), INTENT(IN) :: data_value_sw !< data value at south-western point
       REAL (KIND=wp), INTENT(IN) :: data_value_se !< data value at south-eastern point
       REAL (KIND=wp), INTENT(IN) :: data_value_ne !< data value at north-eastern point
       REAL (KIND=wp), INTENT(IN) :: data_value_nw !< data value at north-western point
       REAL (KIND=wp) :: target_value !< interpolated value, return value

       !calculate bilinear interpolation
        target_value = (1-bwlon) * (1-bwlat) * data_value_sw +  &
                      & bwlon    * (1-bwlat) * data_value_se +  &
                      & bwlon    *  bwlat    * data_value_ne +  &
                      &(1-bwlon) *  bwlat    * data_value_nw 


       END FUNCTION calc_value_bilinear_interpol


END MODULE mo_bilinterpol
