!+  Fortran module to aggregate the CRU near surface climatology to target grid
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_7         2013/01/25 Guenther Zaengl 
!   Parallel threads for ICON and COSMO using Open-MP, 
!   Several bug fixes and optimizations for ICON search algorithm, 
!   particularly for the special case of non-contiguous domains; 
!   simplified namelist control for ICON  
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module to aggregate the CRU near surface climatology to target grid
!!
!> \author Hermann Asensio
MODULE mo_agg_cru

  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i8
  USE mo_kind, ONLY: i4

  !> abort_extpar defined in MODULE utilities_extpar
  USE mo_utilities_extpar, ONLY: abort_extpar

  USE mo_bilinterpol, ONLY:  get_4_surrounding_raw_data_indices, &
    &                        calc_weight_bilinear_interpol, &
    &                        calc_value_bilinear_interpol



IMPLICIT NONE

PRIVATE

PUBLIC :: agg_cru_data_to_target_grid


CONTAINS

  !> Subroutine to aggregate CRU temperature data to the target grid
  SUBROUTINE agg_cru_data_to_target_grid(nrows,ncolumns,ntime)
  !-------------------------------------------------------------------------------------
  ! list of modules which are used as "input"
    USE mo_grid_structures, ONLY: target_grid_def   !< type definition of structure for tg
    !> data type structures form module GRID_structures
    USE mo_GRID_structures, ONLY: reg_lonlat_grid, &
      &                            rotated_lonlat_grid


    ! USE structure which contains the definition of the ICON grid
    USE  mo_icon_grid_data, ONLY: ICON_grid !< structure which contains the definition of the ICON grid

    ! USE structure which contains the definition of the COSMO grid
    USE  mo_cosmo_grid, ONLY: COSMO_grid !< structure which contains the definition of the COSMO grid


    USE mo_base_geometry,   ONLY: geographical_coordinates
    USE mo_base_geometry,   ONLY: cartesian_coordinates
    USE mo_math_constants, ONLY: pi, rad2deg, eps
    USE mo_physical_constants, ONLY: re
    USE mo_additional_geometry,   ONLY: cc2gc,                  &
      &                                 gc2cc
       
     USE mo_search_ll_grid, ONLY: find_reg_lonlat_grid_element_index

     ! USE global data fields (coordinates)
     USE mo_target_grid_data, ONLY: lon_geo, & !< longitude coordinates of the COSMO grid in the geographical system 
                                     lat_geo !< latitude coordinates of the COSMO grid in the geographical system
     USE mo_target_grid_data, ONLY: tg  !< structure with target grid description


     ! "input" fields
     USE mo_cru_data, ONLY : lon_cru, &
       &                   lat_cru, &
       &                   cru_raw_data, &
       &                   cru_grid

     !-------------------------------------------------------------------------------------
     ! list of modules and fields for "output"

     USE mo_cru_target_fields, ONLY: crutemp

                                  
     INTEGER (KIND=i8), INTENT(IN) :: nrows !< number of rows
     INTEGER (KIND=i8), INTENT(IN) :: ncolumns !< number of columns
     INTEGER (KIND=i8), INTENT(IN) :: ntime !< number of times



     REAL (KIND=wp) :: undefined            !< undef value

     INTEGER (KIND=i8) :: undefined_integer ! undef value
     REAL (KIND=wp)    :: default_real

     REAL (KIND=wp) :: deg2rad ! degree to radian


     INTEGER :: i,j,k,l, t ! counters
     INTEGER (KIND=i8) :: ie, je, ke  ! indices for target grid elements

     INTEGER :: idom  ! counter

     INTEGER :: nlon

     REAL(KIND=wp)   :: point_lon, point_lat
     TYPE(geographical_coordinates) :: target_geo_co  !< structure for geographical coordinates of raw data pixel
     TYPE(cartesian_coordinates)  :: target_cc_co     !< coordinates in cartesian system of point 
                                                      !< for which the nearest ICON grid cell is to be determined

     INTEGER        :: k_error     ! error return code

      REAL (KIND=wp) :: bound_north_cosmo !< northern boundary for COSMO target domain
      REAL (KIND=wp) :: bound_south_cosmo !< southern boundary for COSMO target domain

      REAL (KIND=wp) :: point_lon_geo       !< longitude coordinate in geographical system of input point 
      REAL (KIND=wp) :: point_lat_geo       !< latitude coordinate in geographical system of input point
      
       
      INTEGER (KIND=i8) :: western_column     !< the index of the western_column of raw data 
      INTEGER (KIND=i8) :: eastern_column     !< the index of the eastern_column of raw data 
      INTEGER (KIND=i8) :: northern_row       !< the index of the northern_row of raw data 
      INTEGER (KIND=i8) :: southern_row       !< the index of the southern_row of raw data 

      REAL (KIND=wp) :: data_sw
      REAL (KIND=wp) :: data_se
      REAL (KIND=wp) :: data_ne
      REAL (KIND=wp) :: data_nw

      REAL (KIND=wp) :: target_value

      ! global data flag
      LOGICAL :: gldata=.TRUE. ! CRU data are global




      INTEGER (KIND=i4) :: igrid_type  !< target grid type, 1 for ICON, 2 for COSMO, 3 for GME grid

      REAL (KIND=wp) :: bwlon !< weight for bilinear interpolation
      REAL (KIND=wp) :: bwlat !< weight for bilinear interpolation

      REAL (KIND=wp) :: tem_clim_raw(ncolumns,nrows)
      INTEGER :: dpm(12)
      DATA dpm / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /

       tem_clim_raw = 0.0
       DO j=1, nrows
       DO i=1,ncolumns

       DO t=1, ntime
       tem_clim_raw(i,j) = tem_clim_raw(i,j) +  dpm(t) * cru_raw_data(i,j,t)
       ENDDO

       tem_clim_raw(i,j) = 273.15 +  tem_clim_raw(i,j) / 365    ! unit in K instad degC, and yearly mean instead of mohntly means

       ENDDO
       ENDDO
      



      ! loop through all target grid elements
     DO i=1,tg%ie
     DO j=1,tg%je
     DO k=1,tg%ke

!     DO i=200,202
!     DO j=300,301
!     DO k=1,tg%ke


     point_lon_geo = lon_geo(i,j,k) 
     point_lat_geo = lat_geo(i,j,k)

     ! get four surrounding raw data indices

       CALL  get_4_surrounding_raw_data_indices(   cru_grid, &
                                                   lon_cru,           &
                                                   lat_cru,           &
                                                   gldata,            &
                                                   point_lon_geo,      &
                                                   point_lat_geo,      &
                                                   western_column,     &
                                                   eastern_column,     &
                                                   northern_row,       &
                                                   southern_row)
     ! for the bilinear interpolation, use the four data pixels
     ! data(western_column,northern_row)
     ! data(western_column,southern_row)
     ! data(eastern_column,southern_row)
     ! data(eastern_column,northern_row)

   ! PRINT *,'western_column: ', western_column
   ! PRINT *,'eastern_column: ', eastern_column
   ! PRINT *,'northern_row: ', northern_row
   ! PRINT *,'southern_row: ', southern_row



     ! calculate weight for bilinear interpolation
    target_value = -999.
    if ((western_column /= 0) ) then  ! point is not out of data grid range

     CALL calc_weight_bilinear_interpol(point_lon_geo, &
                                        point_lat_geo, &
                                        lon_cru(western_column),      &
                                        lon_cru(eastern_column),      &
                                        lat_cru(northern_row),     &
                                        lat_cru(southern_row),     &
                                        bwlon,         &
                                        bwlat)
     ! the weights are bwlon and bwlat
   !  PRINT *,'bwlon, bwlat', bwlon,bwlat
     ! put all relevent data to an array



     data_sw = tem_clim_raw(western_column,southern_row) 
     data_se = tem_clim_raw(eastern_column,southern_row)
     data_ne = tem_clim_raw(eastern_column,northern_row)
     data_nw = tem_clim_raw(western_column,northern_row)



     ! perform the interpolation

     target_value = calc_value_bilinear_interpol(bwlon, bwlat, &
                                       data_sw, data_se, data_ne, data_nw)
     ELSE
    !  PRINT *,'point_lon_geo: ', point_lon_geo
    !  PRINT *,'point_lat_geo: ', point_lat_geo
    !  PRINT *,'HA debug: western_column: ', western_column
    !  PRINT *,'HA debug: eastern_column: ', eastern_column
    !  PRINT *,'northern_row: ', northern_row
    !  PRINT *,'southern_row: ', southern_row
!
!
!      PRINT *,'Grid element i,j,k: ', i,j,k
!      PRINT *,' lon_geo(i,j,k) : ', lon_geo(i,j,k)
!      PRINT *,' lat_geo(i,j,k) : ', lat_geo(i,j,k)

     ENDIF
     crutemp(i,j,k) = target_value

     ENDDO
     ENDDO
     ENDDO ! loop through all target grid elements





    END SUBROUTINE agg_cru_data_to_target_grid

END MODULE mo_agg_cru
