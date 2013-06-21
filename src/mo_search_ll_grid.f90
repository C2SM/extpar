!+ Fortran module to find grid elements in a regular or rotated lon-lat grid
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V2_0         2013/06/07 Martina Messmer
!  The routine 'find_reg_lonlat_grid_element_index' must be adapted, as the Globcover 2009 data set consists of 6 tiles.
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module to find grid elements in a regular or rotated lon-lat grid
!> \author Hermann Asensio
MODULE mo_search_ll_grid

!> kind parameters are defined in MODULE data_parameters
USE mo_kind, ONLY: wp
USE mo_kind, ONLY: i8
USE mo_kind, ONLY: i4

!> abort_extpar defined in MODULE utilities_extpar
USE mo_utilities_extpar, ONLY: abort_extpar


!> data type structures form module GRID_structures
USE mo_GRID_structures, ONLY: reg_lonlat_grid
USE mo_GRID_structures, ONLY: rotated_lonlat_grid




IMPLICIT NONE

PRIVATE

PUBLIC :: find_reg_lonlat_grid_element_index
PUBLIC :: find_rotated_lonlat_grid_element_index


    CONTAINS


     !> subroutine to find in a prescribed lon-lat grid the grid element indices which inclose given point in geographical system (nearest neighbour)
       !!
       !! for a regular lon-lat grid the index of the data grid element which incloses the given point is directly computed with the help of the grid definition
       !! if the point lies outside the boundary of the data grid, the index is set to an undefinied value ('0')
       !! the definition of the regular lon-lat grid requires 
       !! - the coordinates of the north-western point of the domain ("upper left") regular_grid_info%start_lon_reg and regular_grid_info%start_lat_reg
       !! - the increment regular_grid_info%dlon_reg and regular_grid_info%dlat_reg
       !!   (implict assuming that the grid axis goes from the west to the east and from the south to the north, 
       !!   for the oreder north to south a negative regular_grid_info%dlat_reg value )
       !! - the number of grid elements regular_grid_info%nlon_reg and regular_grid_info%nlat_reg for both directions
       !! the result of the routine could also be interpreted as a result of a nearest neighbour search
       SUBROUTINE find_reg_lonlat_grid_element_index(point_lon_geo,&
                                               point_lat_geo,      &
                                               regular_grid_info,  &
                                               point_lon_index,    &
                                               point_lat_index,    &
                                               tile,               &
                                               regular_tiles_grid_info)

       USE mo_globcover_data,   ONLY: ntiles_globcover

       REAL (KIND=wp), INTENT(in) :: point_lon_geo       !< longitude coordinate in geographical system of input point 
       REAL (KIND=wp), INTENT(in) :: point_lat_geo       !< latitude coordinate in geographical system of input point

       TYPE(reg_lonlat_grid), INTENT(in) :: regular_grid_info !< structure with the definition of regular grid (startlon, startlat etc)

        
       INTEGER  (KIND=i8), INTENT(out):: point_lon_index !< longitude index of point for regular lon-lat grid
       INTEGER  (KIND=i8), INTENT(out):: point_lat_index !< latitude index of point for regular lon-lat grid
       TYPE(reg_lonlat_grid), INTENT(in), OPTIONAL:: regular_tiles_grid_info(:)!< structure with the definition of regular tile grid (startlon, startlat etc), only used for GLOBCOVER, as 6 tiles are present
       INTEGER (KIND=i4), INTENT(out), OPTIONAL:: tile

       ! local variables
       INTEGER (KIND=i4) :: undefined_integer   !< value for undefined integer
       INTEGER (KIND=i4) :: k                   !< counter
       REAL (KIND=wp)    :: point_lon_geo_var   !< longitude coordinate in geographical system of input point, variable for eventual shift


       undefined_integer = 0 ! set undefined to zero

       ! the regular grid domain could cross the 180/-180 meridian or the 0/360 meridian (depending on the chosen value domain for the geographical coordinates)
       ! therefore the longitude coordinate of the point is shiftet to a value which is in any case greater than the startlon_reg_lonlat value, 
       ! implicit assuming that dlon_reg_lonlat is always positiv 
       point_lon_geo_var = point_lon_geo             ! set point_lon_geo_var to value point_lon_geo
       if (point_lon_geo_var < regular_grid_info%start_lon_reg ) then 
           point_lon_geo_var = point_lon_geo + 360.  ! shift coordinate value of point_lon_geo (see above)
       endif

! >mes

       if(present(regular_tiles_grid_info))  then
         tiles: DO k = 1, ntiles_globcover

                   point_lon_index = NINT( (point_lon_geo_var - regular_tiles_grid_info(k)%start_lon_reg)/&
                        &                   regular_tiles_grid_info(k)%dlon_reg) + 1
           
                   point_lat_index = NINT(( point_lat_geo - regular_tiles_grid_info(k)%start_lat_reg)/&
                        &                   regular_tiles_grid_info(k)%dlat_reg) + 1



                   tile = k
                   if(point_lat_index.lt.1 .or. point_lon_index.lt.1 .or. &
                    & point_lon_index.gt.regular_tiles_grid_info(k)%nlon_reg .or. point_lat_index.gt.regular_tiles_grid_info(k)%nlat_reg) then
                     CYCLE tiles
                   else
                     EXIT tiles
                   endif
                   
                 ENDDO tiles

        else
! <mes

        point_lon_index = NINT( (point_lon_geo_var - regular_grid_info%start_lon_reg)/regular_grid_info%dlon_reg) + 1  ! calculate index for regular lon-lat grid
        ! point_lon_geo = regular_grid_info%start_lon_reg + regular_grid_info%dlon_reg * (point_lon_index - 1)

        point_lat_index = NINT(( point_lat_geo - regular_grid_info%start_lat_reg)/regular_grid_info%dlat_reg) + 1 ! calculate index for regular lon-lat grid
        ! point_lat_geo = regular_grid_info%start_lat_reg + regular_grid_info%dlat_reg * (point_lat_index - 1)
      endif


        if (point_lat_index < 1) then ! point out of range of regular lon-lat grid
            point_lat_index = undefined_integer
            point_lat_index = undefined_integer
 
        else if (point_lon_index < 1 ) then ! point out of range of regular lon-lat grid
           point_lon_index = undefined_integer
           point_lat_index = undefined_integer
       
        else if (point_lon_index > regular_grid_info%nlon_reg) then ! point out of range of regular lon-lat grid
           point_lon_index = undefined_integer
           point_lat_index = undefined_integer
       
        else if (point_lat_index > regular_grid_info%nlat_reg) then ! point out of range of regular lon-lat grid
           point_lon_index = undefined_integer
           point_lat_index = undefined_integer
        endif

        
       END SUBROUTINE find_reg_lonlat_grid_element_index

       
       !> subroutine to find in a prescribed rotated lon-lat grid the grid element indices which inclose given point in geographical system (nearest neighbour)
       !!
       !! for a rotated lon-lat grid the index of the data grid element which incloses the given point is directly computed with the help of the grid definition
       !! if the point lies outside the boundary of the data grid, the index is set to an undefinied value ('0')
       !! the result of the routine could also be interpreted as a result of a nearest neighbour search
       SUBROUTINE find_rotated_lonlat_grid_element_index(point_lon_geo,      &
                                               point_lat_geo,      &
                                               rot_grid_info, &
                                               point_rot_lon_index,    &
                                               point_rot_lat_index)
       
       USE mo_utilities_extpar, ONLY: rla2rlarot, &
                                   phi2phirot



      REAL (KIND=wp), INTENT(in) :: point_lon_geo       !< longitude coordinate in geographical system of input point 
      REAL (KIND=wp), INTENT(in) :: point_lat_geo       !< latitude coordinate in geographical system of input point

      TYPE(rotated_lonlat_grid), INTENT(IN) :: rot_grid_info !< !< structure which contains the definition of the rotated grid


      INTEGER  (KIND=i8), INTENT(out):: point_rot_lon_index !< longitude index of point for rotated lon-lat grid
      INTEGER  (KIND=i8), INTENT(out):: point_rot_lat_index !< latitude index of point for rotated lon-lat grid

      ! local variables

      REAL (KIND=wp) :: point_lon_rot  !< longitude coordinate in rotated system of input point
      REAL (KIND=wp) :: point_lat_rot  !< latitude coordinate in rotated system of input point


      INTEGER (KIND=i4) :: undefined_integer   !< value for undefined integer

      undefined_integer = 0 ! set undefined to zero


      ! convert coordinates of given point from geographical system to rotated system
      point_lon_rot = rla2rlarot(point_lat_geo, point_lon_geo, rot_grid_info%pollat, rot_grid_info%pollon, rot_grid_info%polgam)
      point_lat_rot = phi2phirot(point_lat_geo, point_lon_geo, rot_grid_info%pollat, rot_grid_info%pollon)

      
      ! the rotated grid domain could cross the 180/-180 meridian or the 0/360 meridian (depending on the chosen value domain for the rotated coordinates)
      ! therefore the longitude coordinate of the point is shiftet to a value which is in any case greater than the startlon_reg_lonlat value, 
      ! implicit assuming that dlon_reg_lonlat is always positiv 
      if (point_lon_rot < rot_grid_info%startlon_rot ) then 
          point_lon_rot = point_lon_rot + 360.  ! shift coordinate value of point_lon_rot (see above)
      endif

      ! calculate the index of the grid element which incloses the given point
       point_rot_lon_index = NINT( (point_lon_rot - rot_grid_info%startlon_rot)/rot_grid_info%dlon_rot) + 1  ! calculate index for rotated lon-lat grid
       ! point_lon_rot = startlon_rot + dlon_rot * (point_rot_lon_index - 1)

       point_rot_lat_index = NINT(( point_lat_rot - rot_grid_info%startlat_rot)/ rot_grid_info%dlat_rot) + 1 ! calculate index for rotated lon-lat grid
       ! point_lat_rot = startlat_rot + dlat_rot * (point_rot_lat_index - 1)


       if (point_rot_lon_index < 1) then ! point out of range of rotated lon-lat grid
           point_rot_lon_index = undefined_integer
           point_rot_lat_index = undefined_integer
      

        else if (point_rot_lat_index < 1) then ! point out of range of rotated lon-lat grid
           point_rot_lon_index = undefined_integer
           point_rot_lat_index = undefined_integer
       
        else if (point_rot_lon_index > rot_grid_info%nlon_rot) then ! point out of range of rotated lon-lat grid
           point_rot_lon_index = undefined_integer
           point_rot_lat_index = undefined_integer
       
        else if (point_rot_lat_index > rot_grid_info%nlat_rot) then ! point out of range of rotated lon-lat grid
          point_rot_lon_index = undefined_integer
          point_rot_lat_index = undefined_integer
       endif

        
       END SUBROUTINE find_rotated_lonlat_grid_element_index
!----------------------------------------------------------------------------------------------------------------






END MODULE mo_search_ll_grid
