!+ Fortran module to aggregate albedo data to the target grid
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_8         2013-03-12 Frank Brenner
!  introduced MODIS albedo dataset(s) as new external parameter(s)
! V1_9         2013/03/15 Frank Brenner
!  minor bug fix
! V1_13        2013-05-29 Frank Brenner
!  missing values fixed
! V1_14        2014-07-18 Juergen Helmert
!  Combined COSMO Release
! V2_0_3       2015-01-12 Juergen Helmert
!  Bugfix correction covers CSCS SVN r5907-r6359
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module to aggregate albedo data to the target grid
!> \author Frank Brenner, Hermann Asensio
MODULE mo_agg_albedo

  USE mo_kind, ONLY: wp, i4, i8

  USE mo_grid_structures, ONLY: target_grid_def

  USE mo_grid_structures, ONLY: igrid_icon

  USE mo_search_target_grid, ONLY: find_nearest_target_grid_element

  USE mo_albedo_data, ONLY: alb_raw_data_grid, &
       &                    alb_field_row, &
       &                    lon_alb, &
       &                    lat_alb, &
       &                    ntime_alb

  USE mo_albedo_data, ONLY: ialb_type

  USE mo_albedo_routines, ONLY: open_netcdf_ALB_data, &
       &                        close_netcdf_ALB_data, &
       &                        get_one_row_ALB_data, &
       &                        get_pixel_ALB_data

  USE mo_target_grid_data, ONLY: lon_geo, &
       &                         lat_geo, &
       &                         no_raw_data_pixel

  USE mo_target_grid_data, ONLY: search_res !< resolution of ICON grid search index list


  USE mo_bilinterpol, ONLY: get_4_surrounding_raw_data_indices, &
       &                    calc_weight_bilinear_interpol, &
       &                    calc_value_bilinear_interpol

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: agg_alb_data_to_target_grid

CONTAINS

  !> Subroutine to aggregate albedo data to target grid
  SUBROUTINE agg_alb_data_to_target_grid(tg,undefined, path_alb_file, &
       &                     alb_source, alb_field_mom_d)

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    REAL (wp),             INTENT(IN) :: undefined  !< undefined value

    CHARACTER (len=*),     INTENT(in)  :: path_alb_file         !< filename with path for albedo raw data
    CHARACTER (len=*),     INTENT(in)  :: alb_source    !< albedo variable name inside input file

    REAL (wp),             INTENT(OUT) :: alb_field_mom_d(:,:,:,:)   !< monthly mean albedo, output variable

    ! local variables

    INTEGER (i4) :: time_index            !< the index of the time (month) to read in

    INTEGER :: ncid_alb

    INTEGER (i8) :: ie   !< index value for longitude
    INTEGER (i8) :: je   !< index value for latitude
    INTEGER (i8) :: ke   !< counter
    INTEGER (i8) :: start_cell_id !< ID of starting cell for ICON search
    INTEGER (i8) :: i, j, k !< counter
    INTEGER (i8) :: i1, i2

    INTEGER :: row_index !< counter for data row
    INTEGER :: column_index !< counter for data column

    REAL (wp) :: northern_bound !< northern boundary for input data to read for COSMO grid domain
    REAL (wp) :: southern_bound !< southern boundary for input data to read for COSMO grid domain

    INTEGER :: northern_bound_index !< northern boundary for input data to read for COSMO grid domain
    INTEGER :: southern_bound_index !< southern boundary for input data to read for COSMO grid domain

    REAL (wp) ::  alb_sum(1:tg%ie,1:tg%je,1:tg%ke) !< field of target grid with sum of albedo values

    INTEGER (i4) :: point_reg_lon_index          !< longitude index of point for regular lon-lat grid
    INTEGER (i4) :: point_reg_lat_index          !< latitude index of point for regular lon-lat grid

    INTEGER (i4) :: nlon_reg !< number of columns
    INTEGER (i4) :: nlat_reg !< number of rows

    REAL (wp) :: point_lon_geo       !< longitude coordinate in geographical system of input point
    REAL (wp) :: point_lat_geo       !< latitude coordinate in geographical system of input point

    REAL(wp)   :: point_lon, point_lat

    INTEGER (i8) :: western_column     !< the index of the western_column of raw data
    INTEGER (i8) :: eastern_column     !< the index of the eastern_column of raw data
    INTEGER (i8) :: northern_row       !< the index of the northern_row of raw data
    INTEGER (i8) :: southern_row       !< the index of the southern_row of raw data

    REAL (wp) :: bwlon !< weight for bilinear interpolation
    REAL (wp) :: bwlat !< weight for bilinear interpolation

    REAL (wp)   :: alb_point_sw       !< albedo value of the raw data pixel south west
    REAL (wp)   :: alb_point_se       !< albedo value of the raw data pixel south east
    REAL (wp)   :: alb_point_ne       !< albedo value of the raw data pixel north east
    REAL (wp)   :: alb_point_nw       !< albedo value of the raw data pixel north west

    REAL (wp)   :: target_value

    REAL (wp), PARAMETER :: undef_raw     = -0.01_wp
    REAL (wp), PARAMETER :: default_value =  0.07_wp

    ! matrix to save search results

    INTEGER (i8) :: map_ie(alb_raw_data_grid%nlon_reg, alb_raw_data_grid%nlat_reg)
    INTEGER (i8) :: map_je(alb_raw_data_grid%nlon_reg, alb_raw_data_grid%nlat_reg)
    INTEGER (i8) :: map_ke(alb_raw_data_grid%nlon_reg, alb_raw_data_grid%nlat_reg)

    ! global data flag

    LOGICAL :: gldata=.TRUE. ! input data is global

    INTEGER (i4), ALLOCATABLE :: no_valid_raw_data_pixel(:,:,:)

    alb_field_mom_d = undefined

#ifdef DEBUG
    PRINT '(a,f18.10)','  alb_raw_data_grid%start_lat_reg: ', alb_raw_data_grid%start_lat_reg
    PRINT '(a,f18.10)','  alb_raw_data_grid%dlat_reg: ',      alb_raw_data_grid%dlat_reg
    PRINT '(a,f18.10)','  alb_raw_data_grid%start_lon_reg: ', alb_raw_data_grid%start_lon_reg
    PRINT '(a,f18.10)','  alb_raw_data_grid%dlon_reg: ',      alb_raw_data_grid%dlon_reg
    PRINT '(a,2i6)',   '  alb_raw_data_grid nlon, nlat: ',    alb_raw_data_grid%nlon_reg, alb_raw_data_grid%nlat_reg
#endif
    
    ! determine northern and southern boundary for input data to read for COSMO grid domain
    northern_bound =  MAXVAL(lat_geo) + 1.0_wp * alb_raw_data_grid%dlat_reg
    ! One row of albedo data north of northern boundary of COSMO grid domain
    northern_bound = MIN(90._wp, northern_bound)                                 ! Check for the poles
    ! calculate index for regular lon-lat grid
    northern_bound_index = NINT((northern_bound - alb_raw_data_grid%start_lat_reg) / alb_raw_data_grid%dlat_reg) + 1
    !< check for bounds
    IF (northern_bound_index < 1) THEN
      northern_bound_index = 1
    ELSE IF (northern_bound_index > alb_raw_data_grid%nlat_reg) THEN
      northern_bound_index=alb_raw_data_grid%nlat_reg
    ENDIF

    southern_bound = MINVAL(lat_geo) - 1.0_wp * alb_raw_data_grid%dlat_reg
    ! One row of albedo data south of southern boundary of COSMO grid domain
    southern_bound = MAX(-90._wp, southern_bound)                               ! Check for the poles
    ! calculate index for regular lon-lat grid
    southern_bound_index = NINT((southern_bound - alb_raw_data_grid%start_lat_reg) / alb_raw_data_grid%dlat_reg) + 1
    !< check for bounds
    IF (southern_bound_index < 1) THEN
      southern_bound_index = 1
    ELSE IF (southern_bound_index > alb_raw_data_grid%nlat_reg) THEN
      southern_bound_index=alb_raw_data_grid%nlat_reg
    ENDIF

#ifdef DEBUG
    PRINT *,' MAXVAL(lat_geo): ', MAXVAL(lat_geo)
    PRINT *,' MINVAL(lat_geo): ', MINVAL(lat_geo)
    PRINT *, 'northern_bound: ', northern_bound
    PRINT *, 'southern_bound: ', southern_bound

    PRINT *, 'northern_bound_index: ', northern_bound_index
    PRINT *, 'southern_bound_index: ', southern_bound_index
    PRINT *,'lat_alb(northern_bound_index): ', lat_alb(northern_bound_index)
    PRINT *,'lat_alb(southern_bound_index): ', lat_alb(southern_bound_index)
#endif
    
    nlon_reg = alb_raw_data_grid%nlon_reg
    nlat_reg = alb_raw_data_grid%nlat_reg
    start_cell_id = 1

    ! open netcdf file with albedo data
#ifdef DEBUG
    PRINT *, 'mo_agg_albedo: read from ', TRIM(path_alb_file)
#endif
    CALL open_netcdf_ALB_data(path_alb_file, ncid_alb)
#ifdef DEBUG
    PRINT *, 'mo_agg_albedo: albedo nc input file opened'
#endif
    ALLOCATE(no_valid_raw_data_pixel(tg%ie, tg%je, tg%ke))

    ! read in albedo data row by row and assign albedo raw data pixel to target grid
    ! start loop over albedo raw data
    time_loop: DO time_index = 1, ntime_alb

      no_raw_data_pixel       = 0
      no_valid_raw_data_pixel = 0
      alb_sum(:,:,:)          = 0.0_wp

      data_rows: DO row_index = southern_bound_index, northern_bound_index
#ifdef DEBUG
        PRINT '(a,i6)','  row index ', row_index
#endif        
        ! get input raw data row
        CALL get_one_row_ALB_data(ncid_alb,      &
             &                    nlon_reg,      &
             &                    nlat_reg,      &
             &                    ntime_alb,     &
             &                    row_index,     &
             &                    time_index,    &
             &                    alb_field_row, &
             &                    alb_source)

        column: DO column_index = 1, alb_raw_data_grid%nlon_reg

          IF (time_index == 1) THEN

            point_lon = lon_alb(column_index)
            point_lat = lat_alb(row_index)

            ! Reset start cell when entering a new row or when the previous data point was outside
            ! the model domain
            IF (tg%igrid_type == igrid_icon .AND. (column_index == 1 .OR. start_cell_id == 0)) THEN
              i1 = NINT(point_lon*search_res)
              i2 = NINT(point_lat*search_res)
              start_cell_id = tg%search_index(i1,i2)
              IF (start_cell_id == 0) EXIT column ! in this case, the whole row is empty
            ENDIF

            CALL find_nearest_target_grid_element( point_lon,     &
                 &                                 point_lat,     &
                 &                                 tg,            &
                 &                                 start_cell_id, &
                 &                                 ie,            &
                 &                                 je,            &
                 &                                 ke)

            map_ie(column_index, row_index) = ie
            map_je(column_index, row_index) = je
            map_ke(column_index, row_index) = ke

          ELSE

            ie = map_ie(column_index, row_index)
            je = map_je(column_index, row_index)
            ke = map_ke(column_index, row_index)

          ENDIF

          IF ((ie /= 0).AND.(je /= 0).AND.(ke /= 0))THEN
            no_raw_data_pixel(ie,je,ke) = no_raw_data_pixel(ie,je,ke) + 1
            IF (alb_field_row(column_index) > 0.02_wp) THEN
              no_valid_raw_data_pixel(ie,je,ke) = no_valid_raw_data_pixel(ie,je,ke) + 1
              ! count raw data pixel within COSMO grid element
              alb_sum(ie,je,ke)  = alb_sum(ie,je,ke) + alb_field_row(column_index) ! sum data values
            ENDIF
          ENDIF
        ENDDO column
      END DO data_rows

#ifdef DEBUG
      PRINT *,'MAXVAL(no_raw_data_pixel): ',MAXVAL(no_raw_data_pixel)
      PRINT *,'MINVAL(no_raw_data_pixel): ',MINVAL(no_raw_data_pixel)
      PRINT *,'tg: ',tg%ie, tg%je, tg%ke, tg%minlon, tg%maxlon, tg%minlat, tg%maxlat
#endif
      DO k = 1, tg%ke
        DO j = 1, tg%je
          DO i = 1, tg%ie

            IF (no_valid_raw_data_pixel(i,j,k) /= 0) THEN

              alb_field_mom_d(i,j,k,time_index) = alb_sum(i,j,k) / REAL(no_valid_raw_data_pixel(i,j,k))   ! calculate arithmetic mean

            ELSEIF (no_valid_raw_data_pixel(i,j,k) == 0) THEN

#ifdef DEBUG
              PRINT *,' problem getting a value ...'
#endif
              point_lon_geo = lon_geo(i,j,k)
              point_lat_geo = lat_geo(i,j,k)
#ifdef DEBUG
              PRINT *,' ... at lon, lat ', point_lon_geo, point_lat_geo
#endif
              ! get four surrounding raw data indices
              CALL  get_4_surrounding_raw_data_indices(alb_raw_data_grid, &
                   &                                   lon_alb,           &
                   &                                   lat_alb,           &
                   &                                   gldata,             &
                   &                                   point_lon_geo,      &
                   &                                   point_lat_geo,      &
                   &                                   western_column,     &
                   &                                   eastern_column,     &
                   &                                   northern_row,       &
                   &                                   southern_row)
#ifdef DEBUG
              PRINT *,' ... indices (wesn) ', western_column, eastern_column, northern_row, southern_row
#endif              
              IF ( (western_column /= 0) .AND. (eastern_column /= 0) .AND. &
                   (northern_row /= 0)   .AND. (southern_row /= 0) ) THEN

                ! get albedo data for the pixel south west
                point_reg_lon_index = western_column
                point_reg_lat_index = southern_row
                CALL  get_pixel_ALB_data(ncid_alb,            &
                     &                   nlon_reg,            &
                     &                   nlat_reg,            &
                     &                   ntime_alb,           &
                     &                   point_reg_lon_index, &
                     &                   point_reg_lat_index, &
                     &                   time_index,          &
                     &                   alb_point_sw,        &
                     &                   alb_source)

                ! get albedo data for the pixel south east
                point_reg_lon_index = eastern_column
                point_reg_lat_index = southern_row
                CALL  get_pixel_ALB_data( ncid_alb,           &
                     &                   nlon_reg,            &
                     &                   nlat_reg,            &
                     &                   ntime_alb,           &
                     &                   point_reg_lon_index, &
                     &                   point_reg_lat_index, &
                     &                   time_index,          &
                     &                   alb_point_se,        &
                     &                   alb_source)

                ! get albedo data for the pixel north east
                point_reg_lon_index = eastern_column
                point_reg_lat_index = northern_row
                CALL  get_pixel_ALB_data( ncid_alb,           &
                     &                   nlon_reg,            &
                     &                   nlat_reg,            &
                     &                   ntime_alb,           &
                     &                   point_reg_lon_index, &
                     &                   point_reg_lat_index, &
                     &                   time_index,          &
                     &                   alb_point_ne,        &
                     &                   alb_source)

                ! get albedo data for the pixel north west
                point_reg_lon_index = western_column
                point_reg_lat_index = northern_row
                CALL  get_pixel_ALB_data( ncid_alb,           &
                     &                   nlon_reg,            &
                     &                   nlat_reg,            &
                     &                   ntime_alb,           &
                     &                   point_reg_lon_index, &
                     &                   point_reg_lat_index, &
                     &                   time_index,          &
                     &                   alb_point_nw,        &
                     &                   alb_source)

                ! calculate weight for bilinear interpolation
                CALL calc_weight_bilinear_interpol(point_lon_geo, &
                     &                   point_lat_geo,           &
                     &                   lon_alb(western_column), &
                     &                   lon_alb(eastern_column), &
                     &                   lat_alb(northern_row),   &
                     &                   lat_alb(southern_row),   &
                     &                   bwlon, bwlat)

                ! the weights are bwlon and bwlat
#ifdef DEBUG
                PRINT *, bwlon
                PRINT *, bwlat
#endif
                IF (ialb_type == 2) THEN
                  IF ((alb_point_ne > undef_raw) .OR. (alb_point_nw > undef_raw) .OR. &
                      (alb_point_se > undef_raw) .OR. (alb_point_sw > undef_raw)) THEN
                    IF (alb_point_ne < undef_raw) THEN
                      IF (alb_point_nw > undef_raw) THEN
                        alb_point_ne = alb_point_nw
                      ELSEIF (alb_point_se > undef_raw) THEN
                        alb_point_ne = alb_point_se
                      ELSE
                        alb_point_ne = alb_point_sw
                      ENDIF
                    ENDIF
                    IF (alb_point_nw < undef_raw) THEN
                      IF (alb_point_ne > undef_raw) THEN
                        alb_point_nw = alb_point_ne
                      ELSEIF (alb_point_sw > undef_raw) THEN
                        alb_point_nw = alb_point_sw
                      ELSE
                        alb_point_nw = alb_point_se
                      ENDIF
                    ENDIF
                    IF (alb_point_se < undef_raw) THEN
                      IF (alb_point_sw > undef_raw) THEN
                        alb_point_se = alb_point_sw
                      ELSEIF (alb_point_ne > undef_raw) THEN
                        alb_point_se = alb_point_ne
                      ELSE
                        alb_point_se = alb_point_nw
                      ENDIF
                    ENDIF
                    IF (alb_point_sw < undef_raw) THEN
                      IF (alb_point_se > undef_raw) THEN
                        alb_point_sw = alb_point_se
                      ELSEIF (alb_point_nw > undef_raw) THEN
                        alb_point_sw = alb_point_nw
                      ELSE
                        alb_point_sw = alb_point_ne
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF

                ! perform the interpolation

                IF (alb_point_sw > 0.02_wp .AND. alb_point_se > 0.02_wp .AND. &
                    alb_point_ne > 0.02_wp .AND. alb_point_nw > 0.02_wp) THEN
                  target_value = calc_value_bilinear_interpol(bwlon, bwlat, alb_point_sw, alb_point_se, alb_point_ne, alb_point_nw)
                  IF (target_value < 0.02_wp) THEN
#ifdef DEBUG
                    PRINT *,'Interpolation gone wrong! ', target_value,alb_point_sw, alb_point_se, alb_point_ne, &
                         &                                alb_point_nw,bwlon, bwlat
#endif
                    target_value = -999.0_wp
                  ENDIF
                ELSE
                  ! assume missing value - will be fixed later in the cross check
                  target_value = -999.0_wp
                ENDIF
                alb_field_mom_d(i,j,k,time_index) = target_value

              ELSE ! grid element outside target grid
                alb_field_mom_d(i,j,k,time_index) = default_value
              ENDIF

              IF (bwlon > 1) THEN !calculation of bwlon gone wrong (eastern border)
                alb_field_mom_d(i,j,k,time_index) = -999.0_wp
              ENDIF

              IF (bwlon < 0) THEN !calculation of bwlon gone wrong (western border)
                alb_field_mom_d(i,j,k,time_index) = -999.0_wp
              ENDIF

            ELSE
              ! assume missing value if no valid data point was present, will be fixed later in the cross check
            ENDIF

          ENDDO ! i
        ENDDO ! j
      ENDDO ! k

    ENDDO time_loop

    CALL  close_netcdf_ALB_data(ncid_alb)

    DEALLOCATE(no_valid_raw_data_pixel)

  END SUBROUTINE agg_alb_data_to_target_grid

END MODULE mo_agg_albedo
