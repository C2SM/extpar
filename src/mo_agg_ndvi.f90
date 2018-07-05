!+ Fortran module to aggregate NDVI data to the target grid
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
! V2.0_3         2015-01-12 Juergen Helmert
!  Bugfix correction covers CSCS SVN r5907-r6359
!  Bugfix in loop index boundaries
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module to aggregate NDVI data to the target grid
!> \author Hermann Asensio
MODULE mo_agg_ndvi

  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i4
  USE mo_kind, ONLY: i8


  !> abort_extpar defined in MODULE utilities_extpar
  USE mo_utilities_extpar, ONLY: abort_extpar



  USE mo_grid_structures, ONLY: reg_lonlat_grid, &
       &                           rotated_lonlat_grid, &
       &                           target_grid_def, &
       &                           icosahedral_triangular_grid

  USE mo_grid_structures, ONLY: igrid_icon
  USE mo_grid_structures, ONLY: igrid_cosmo


  USE mo_search_ll_grid, ONLY: find_reg_lonlat_grid_element_index, &
       &                          find_rotated_lonlat_grid_element_index

  USE mo_search_target_grid, ONLY: find_nearest_target_grid_element

  IMPLICIT NONE

  PRIVATE
  !
  PUBLIC :: agg_ndvi_data_to_target_grid

CONTAINS

  !> Subroutine to aggregate NDVI data to target grid
  SUBROUTINE agg_ndvi_data_to_target_grid(tg,undefined, path_ndvi_file)

    USE mo_ndvi_data, ONLY: ndvi_raw_data_grid, &
         ndvi_field_row_mom, &
         ndvi_field_row, &
         lon_ndvi, &
         lat_ndvi, &
         ntime_ndvi

    USE mo_ndvi_tg_fields, ONLY: ndvi_field, &
         ndvi_max, &
         ndvi_field_mom, &
         ndvi_ratio_mom

    USE mo_ndvi_routines, ONLY: open_netcdf_NDVI_data, &
         close_netcdf_NDVI_data, &
         get_one_row_NDVI_data, &
         get_pixel_NDVI_data

    USE mo_target_grid_data, ONLY: lon_geo, &
         lat_geo, &
         no_raw_data_pixel

    USE mo_target_grid_data, ONLY: search_res !< resolution of ICON grid search index list


    ! USE structure which contains the definition of the COSMO grid
    USE  mo_cosmo_grid, ONLY: COSMO_grid !< structure which contains the definition of the COSMO


    ! USE structure which contains the definition of the ICON grid
    USE  mo_icon_grid_data, ONLY: ICON_grid !< structure which contains the definition of the ICON grid

    USE mo_bilinterpol, ONLY: get_4_surrounding_raw_data_indices, &
         calc_weight_bilinear_interpol, &
         calc_value_bilinear_interpol





    IMPLICIT NONE

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description

    REAL (KIND=wp), INTENT(IN) :: undefined  !< undefined value

    CHARACTER (len=*), INTENT(in) :: path_ndvi_file         !< filename with path for NDVI raw data


    ! local variables

    INTEGER (KIND=i4)    :: time_index            !< the index of the time (month) to read in

    INTEGER :: ncid_ndvi

    REAL (KIND=wp) :: default_value

    INTEGER (KIND=i8) :: ie   !< index value for longitude
    INTEGER (KIND=i8)  :: je   !< index value for latitude
    INTEGER (KIND=i8)  :: ke   !< counter
    INTEGER (KIND=i8) :: start_cell_id !< ID of starting cell for ICON search
    INTEGER (KIND=i8) :: i,j,k !< counter
    INTEGER (KIND=i8) :: i1, i2

    INTEGER :: row_index !< counter for NDVI data row
    INTEGER :: column_index !< counter for NDVI data column
    REAL (KIND=wp) :: northern_bound !< northern boundary for NDVI data to read for COSMO grid domain
    REAL (KIND=wp) :: southern_bound !< southern boundary for NDVI data to read for COSMO grid domain

    INTEGER :: northern_bound_index !< northern boundary for NDVI data to read for COSMO grid domain
    INTEGER :: southern_bound_index !< southern boundary for NDVI data to read for COSMO grid domain


    REAL (KIND=wp) ::  ndvi_sum(1:tg%ie,1:tg%je,1:tg%ke) !< field of target grid with sum of NDVI values

    INTEGER (KIND=i4) :: point_rot_lon_index          !< longitude index of point for rotated lon-lat grid
    INTEGER (KIND=i4) :: point_rot_lat_index          !< latitude index of point for rotated lon-lat grid


    INTEGER (KIND=i4) :: point_reg_lon_index          !< longitude index of point for regular lon-lat grid
    INTEGER (KIND=i4) :: point_reg_lat_index          !< latitude index of point for regular lon-lat grid

    INTEGER (KIND=i4) :: nlon_reg !< number of columns
    INTEGER (KIND=i4) :: nlat_reg !< number of rows

    REAL (KIND=wp) :: point_lon_geo       !< longitude coordinate in geographical system of input point
    REAL (KIND=wp) :: point_lat_geo       !< latitude coordinate in geographical system of input point

    REAL(KIND=wp)   :: point_lon, point_lat

    INTEGER (KIND=i8) :: western_column     !< the index of the western_column of raw data
    INTEGER (KIND=i8) :: eastern_column     !< the index of the eastern_column of raw data
    INTEGER (KIND=i8) :: northern_row       !< the index of the northern_row of raw data
    INTEGER (KIND=i8) :: southern_row       !< the index of the southern_row of raw data

    REAL (KIND=wp) :: bwlon !< weight for bilinear interpolation
    REAL (KIND=wp) :: bwlat !< weight for bilinear interpolation

    REAL (KIND=wp)   :: ndvi_point_sw       !< value of the NDVI raw data pixel south west
    REAL (KIND=wp)   :: ndvi_point_se       !< value of the NDVI raw data pixel south east
    REAL (KIND=wp)   :: ndvi_point_ne       !< value of the NDVI raw data pixel north east
    REAL (KIND=wp)   :: ndvi_point_nw       !< value of the NDVI raw data pixel north west

    REAL (KIND=wp)   :: target_value


    ! matrix to save search results

    INTEGER (KIND=i8) :: map_ie(ndvi_raw_data_grid%nlon_reg, ndvi_raw_data_grid%nlat_reg)
    INTEGER (KIND=i8) :: map_je(ndvi_raw_data_grid%nlon_reg, ndvi_raw_data_grid%nlat_reg)
    INTEGER (KIND=i8) :: map_ke(ndvi_raw_data_grid%nlon_reg, ndvi_raw_data_grid%nlat_reg)

    ! buffer for ndvi data for one month
    REAL (KIND=wp)   :: ndvi_raw_data(ndvi_raw_data_grid%nlon_reg, ndvi_raw_data_grid%nlat_reg)

    ! global data flag
    LOGICAL :: gldata=.TRUE. ! NDVI data are global





    default_value = -1.
    ndvi_field = default_value
    ndvi_field_mom = default_value

    !HA debug:
    print *,' ndvi_raw_data_grid: ', ndvi_raw_data_grid
    print *,' ndvi_raw_data_grid%start_lat_reg: ', ndvi_raw_data_grid%start_lat_reg
    print *,' ndvi_raw_data_grid%dlat_reg: ', ndvi_raw_data_grid%dlat_reg
    print *,' ndvi_raw_data_grid%nlat_reg: ', ndvi_raw_data_grid%nlat_reg

    IF (tg%igrid_type == igrid_cosmo) THEN


      ! determine northern and southern boundary for NDVI data to read for COSMO grid domain
      northern_bound =  MAXVAL(lat_geo) + 1.* ABS(ndvi_raw_data_grid%dlat_reg)
      ! One row of NDVI data north of northern boundary of COSMO grid domain
      northern_bound = MIN(90._wp, northern_bound)                                 ! Check for the poles

      northern_bound_index = NINT(ABS((northern_bound - ndvi_raw_data_grid%start_lat_reg)/&
           ABS(ndvi_raw_data_grid%dlat_reg))) +1  ! calculate index for regular lon-lat NDVI grid

      if (northern_bound_index < 1) then
        northern_bound_index = 1 !< check for bounds
      else if (northern_bound_index > ndvi_raw_data_grid%nlat_reg) then
        northern_bound_index=ndvi_raw_data_grid%nlat_reg
      endif
      southern_bound = MINVAL(lat_geo) - 1.* ABS(ndvi_raw_data_grid%dlat_reg)
      ! One row of NDVI data south of southern boundary of COSMO grid domain
      southern_bound = MAX(-90._wp, southern_bound)                               ! Check for the poles

      southern_bound_index = NINT(ABS((southern_bound - ndvi_raw_data_grid%start_lat_reg)/&
           ABS(ndvi_raw_data_grid%dlat_reg))) + 2  ! calculate index for regular lon-lat NDVI grid
      ! DWD  ABS(ndvi_raw_data_grid%dlat_reg))) +1
      if (southern_bound_index < 1) then
        southern_bound_index = 1 !< check for bounds
      else if (southern_bound_index > ndvi_raw_data_grid%nlat_reg) then
        southern_bound_index=ndvi_raw_data_grid%nlat_reg
      endif

    END IF
    IF (tg%igrid_type == igrid_icon) THEN
      northern_bound_index=1
      southern_bound_index=ndvi_raw_data_grid%nlat_reg
    END IF

    !HA debug:
    print *,' MAXVAL(lat_geo): ', MAXVAL(lat_geo)
    print *,' MINVAL(lat_geo): ', MINVAL(lat_geo)
    print *, 'northern_bound: ', northern_bound
    print *, 'southern_bound: ', southern_bound

    print *, 'northern_bound_index: ', northern_bound_index
    print *, 'southern_bound_index: ', southern_bound_index
    print *,'lat_ndvi(northern_bound_index): ', lat_ndvi(northern_bound_index)
    print *,'lat_ndvi(southern_bound_index): ', lat_ndvi(southern_bound_index)


    nlon_reg = ndvi_raw_data_grid%nlon_reg
    nlat_reg = ndvi_raw_data_grid%nlat_reg
    start_cell_id = 1

    ! open netcdf file with NDVI data
    CALL open_netcdf_NDVI_data(path_ndvi_file, &
         ncid_ndvi)

    ! read in NDVI data row by row and assign NDVI raw data pixel to COSMO grid
    ! start loop over NDVI raw data
    time_loop: DO time_index=1,12
      no_raw_data_pixel = 0 ! set count to 0
      ndvi_sum = 0.  ! set sum to 0

      data_rows: DO row_index=northern_bound_index,southern_bound_index
        ! get NDVI raw data row

        CALL get_one_row_NDVI_data(ncid_ndvi,      &
             nlon_reg,      &
             nlat_reg,     &
             ntime_ndvi,     &
             row_index,           &
             time_index,         &
             ndvi_field_row)

        ! store ndvi data for subsequent filling algorithm
        ndvi_raw_data(1:ndvi_raw_data_grid%nlon_reg, row_index) = ndvi_field_row(1:ndvi_raw_data_grid%nlon_reg)


        column: DO column_index=1, ndvi_raw_data_grid%nlon_reg

          IF (time_index == 1) THEN
            point_lon = lon_ndvi(column_index)
            point_lat = lat_ndvi(row_index)

            ! Reset start cell when entering a new row or when the previous data point was outside
            ! the model domain
            IF (tg%igrid_type == igrid_icon .AND. (column_index == 1 .OR. start_cell_id == 0)) THEN
              i1 = NINT(point_lon*search_res)
              i2 = NINT(point_lat*search_res)
              start_cell_id = tg%search_index(i1,i2)
              IF (start_cell_id == 0) EXIT column ! in this case, the whole row is empty
            ENDIF

            CALL find_nearest_target_grid_element( point_lon, &
                 &      point_lat, &
                 &      tg,        &
                 &      start_cell_id, &
                 &      ie,      &
                 &      je,      &
                 &      ke)

            map_ie(column_index, row_index) = ie
            map_je(column_index, row_index) = je
            map_ke(column_index, row_index) = ke
          ELSE
            ie = map_ie(column_index, row_index)
            je = map_je(column_index, row_index)
            ke = map_ke(column_index, row_index)
          ENDIF

          IF ((ie /= 0).AND.(je/=0).AND.(ke/=0))THEN

            no_raw_data_pixel(ie,je,ke) = no_raw_data_pixel(ie,je,ke) + 1
            ! count raw data pixel within COSMO/ICON grid element
            ndvi_sum(ie,je,ke)  = ndvi_sum(ie,je,ke) + ndvi_field_row(column_index) ! sum data values
          ENDIF


        ENDDO column

      END DO data_rows

      print *,'MAXVAL of NDVI_sum is ', MAXVAL(ndvi_sum)






      !HA debug:
      print *,'MAXVAL(no_raw_data_pixel): ',MAXVAL(no_raw_data_pixel)
      print *,'MINVAL(no_raw_data_pixel): ',MINVAL(no_raw_data_pixel)

      !WHERE (no_raw_data_pixel > 0)
      !    ndvi_field = ndvi_sum / no_raw_data_pixel   ! calculate arithmetic mean
      !END WHERE

      PRINT *,'target grid tg: ',tg%ie, tg%je, tg%ke, tg%minlon, tg%maxlon, tg%minlat, tg%maxlat

      DO k=1, tg%ke
        DO j=1, tg%je
          DO i=1, tg%ie

            !print *,' i,j,k: ',i,j,k
            !print *,'no_raw_data_pixel(i,j,k): ', no_raw_data_pixel(i,j,k)

            IF (no_raw_data_pixel(i,j,k) /= 0) THEN
              !---------------------------------------------------------------------------------------------------------------------
              !   PRINT *,'no_raw_data_pixel(i,j,k) /= 0'
              !   print *,'no_raw_data_pixel(i,j,k): ', no_raw_data_pixel(i,j,k)
              !   print *,'ndvi_sum(i,j,k): ', ndvi_sum(i,j,k)
              !   PRINT *,'SIZE(ndvi_field): ', SIZE(ndvi_field)
              !   PRINT *,'SHAPE(ndvi_field): ', SHAPE(ndvi_field)
              !
              !   print *,'ndvi_field(i,j,k): ', ndvi_field(i,j,k)
              ndvi_field(i,j,k) = ndvi_sum(i,j,k) / REAL(no_raw_data_pixel(i,j,k),wp)   ! calculate arithmetic mean
              !    PRINT *,' ndvi_field(i,j,k) = ndvi_sum(i,j,k) / REAL(no_raw_data_pixel(i,j,k))'

            ENDIF

          ENDDO
        ENDDO
      ENDDO

      print *,'MAXVAL/MINVAL of NDVI_field for no_raw_data_pixel >0 is ', MAXVAL(ndvi_field), MINVAL(ndvi_field)

      DO k=1, tg%ke
        DO j=1, tg%je
          DO i=1, tg%ie

            !print *,' i,j,k: ',i,j,k
            !print *,'no_raw_data_pixel(i,j,k): ', no_raw_data_pixel(i,j,k)

            IF (no_raw_data_pixel(i,j,k) == 0) THEN
              !     ELSE ! bilinear interpolation
              !--------------------------------------------------------------------------------------------------------------------

              !  print *,'start bilinear interpolation'

              point_lon_geo = lon_geo(i,j,k)
              point_lat_geo = lat_geo(i,j,k)

              ! get four surrounding raw data indices

!DIR$ NOINLINE
              CALL  get_4_surrounding_raw_data_indices(   ndvi_raw_data_grid, &
                   lon_ndvi,           &
                   lat_ndvi,           &
                   gldata,             &
                   point_lon_geo,      &
                   point_lat_geo,      &
                   western_column,     &
                   eastern_column,     &
                   northern_row,       &
                   southern_row)
              !        print *,'western_column, eastern_column, northern_row, southern_row'
              !        print *, western_column, eastern_column, northern_row, southern_row

              target_value = -999.

              IF ( (western_column /= 0) .AND. &
                   (eastern_column /= 0) .AND. &
                   (northern_row /= 0)   .AND. &
                   (southern_row /= 0)        ) THEN


                ! get NDVI data for the pixel south west
                point_reg_lon_index = western_column
                point_reg_lat_index = southern_row
                CALL  get_pixel_NDVI_data( ncid_ndvi,           &
                     nlon_reg,           &
                     nlat_reg,           &
                     ntime_ndvi,          &
                     point_reg_lon_index,     &
                     point_reg_lat_index,     &
                     time_index,          &
                     ndvi_point_sw)


                ! get NDVI data for the pixel south east
                point_reg_lon_index = eastern_column
                point_reg_lat_index = southern_row
                CALL  get_pixel_NDVI_data( ncid_ndvi,           &
                     nlon_reg,           &
                     nlat_reg,           &
                     ntime_ndvi,          &
                     point_reg_lon_index,     &
                     point_reg_lat_index,     &
                     time_index,          &
                     ndvi_point_se)

                ! get NDVI data for the pixel north east
                point_reg_lon_index = eastern_column
                point_reg_lat_index = northern_row
                CALL  get_pixel_NDVI_data( ncid_ndvi,           &
                     nlon_reg,           &
                     nlat_reg,           &
                     ntime_ndvi,          &
                     point_reg_lon_index,     &
                     point_reg_lat_index,     &
                     time_index,          &
                     ndvi_point_ne)


                ! get NDVI data for the pixel north west
                point_reg_lon_index = western_column
                point_reg_lat_index = northern_row
                CALL  get_pixel_NDVI_data( ncid_ndvi,           &
                     nlon_reg,           &
                     nlat_reg,           &
                     ntime_ndvi,          &
                     point_reg_lon_index,     &
                     point_reg_lat_index,     &
                     time_index,          &
                     ndvi_point_nw)



                ! calculate weight for bilinear interpolation

                CALL calc_weight_bilinear_interpol(point_lon_geo, &
                     point_lat_geo, &
                     lon_ndvi(western_column),      &
                     lon_ndvi(eastern_column),      &
                     lat_ndvi(northern_row),     &
                     lat_ndvi(southern_row),     &
                     bwlon,         &
                     bwlat)

                ndvi_point_sw = ndvi_raw_data(western_column, southern_row)
                ndvi_point_se = ndvi_raw_data(eastern_column, southern_row)
                ndvi_point_ne = ndvi_raw_data(eastern_column, northern_row)
                ndvi_point_nw = ndvi_raw_data(western_column, northern_row)



                ! the weights are bwlon and bwlat
                !  PRINT *,'bwlon, bwlat', bwlon,bwlat

                ! perform the interpolation

                target_value = calc_value_bilinear_interpol(bwlon, bwlat, &
                     ndvi_point_sw, ndvi_point_se, ndvi_point_ne, ndvi_point_nw)

                if (target_value < 0.) print*,'Caution target_value < 0: ', bwlon, bwlat, &
                     ndvi_point_sw, ndvi_point_se, ndvi_point_ne, ndvi_point_nw


                !        ndvi_field(i,j,k) = MAX(0._wp,target_value)
                ndvi_field(i,j,k) =target_value
              ELSE ! grid element outside target grid
                ndvi_field(i,j,k) = default_value
              ENDIF

              !---------------------------------------------------------------------------------------------------------------------
            ENDIF

          ENDDO
        ENDDO
      ENDDO

      print *,'ndvi_field determined'
      print *,'time_index:', time_index
      print *,'MAXVAL of NDVI is ', MAXVAL(ndvi_field)

      ndvi_field_mom(:,:,:,time_index) = ndvi_field(:,:,:)

    END DO time_loop

    CALL  close_netcdf_NDVI_data(ncid_ndvi)


    ! calculate NDVI Max and NDVI ratio

    ndvi_max = MAXVAL(ndvi_field_mom,4) ! maximum for the "time" dimension

    ndvi_ratio_mom = default_value

    DO k=1, tg%ke
      DO j=1, tg%je
        DO i=1, tg%ie
          IF(ndvi_max(i,j,k) /= 0.) THEN
            DO time_index=1,12
              ndvi_ratio_mom(i,j,k,time_index) = ndvi_field_mom(i,j,k,time_index) / ndvi_max(i,j,k)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    ! print *,' subroutine aggregate_ndvi_raw_data_to_cosmo_grid: MAXVAL(cosmo_ndvi_field): ',MAXVAL(cosmo_ndvi_field)


  END SUBROUTINE agg_ndvi_data_to_target_grid


END MODULE mo_agg_ndvi
