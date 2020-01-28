!+ Fortran module to aggregate EMISS data to the target grid
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
!> Fortran module to aggregate EMISS data to the target grid
!> \author Hermann Asensio
MODULE mo_agg_emiss

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4

  USE mo_grid_structures,       ONLY: target_grid_def, & 
       &                              igrid_icon, & 
       &                              igrid_cosmo

  USE mo_search_target_grid,    ONLY: find_nearest_target_grid_element

  USE mo_emiss_data,            ONLY: emiss_raw_data_grid, &
       &                              emiss_field_row, &
       &                              lon_emiss, &
       &                              lat_emiss, &
       &                              ntime_emiss

  USE mo_emiss_tg_fields,       ONLY: emiss_field, &
       &                              emiss_max, &
       &                              emiss_field_mom, &
       &                              emiss_ratio_mom
 
  USE mo_emiss_routines,        ONLY: open_netcdf_EMISS_data, &
       &                              close_netcdf_EMISS_data, &
       &                              get_one_row_EMISS_data, &
       &                              get_pixel_EMISS_data
 
  USE mo_target_grid_data,      ONLY: lon_geo, &
                                      lat_geo, &
                                      no_raw_data_pixel
        
  USE mo_target_grid_data,      ONLY: search_res !< resolution of ICON grid search index list

  USE mo_bilinterpol,           ONLY: get_4_surrounding_raw_data_indices, &
      &                               calc_weight_bilinear_interpol, &
      &                               calc_value_bilinear_interpol

  IMPLICIT NONE

  PRIVATE
 
  PUBLIC :: agg_emiss_data_to_target_grid

  CONTAINS

  !> Subroutine to aggregate EMISS data to target grid
  SUBROUTINE agg_emiss_data_to_target_grid(tg, path_emiss_file)

    IMPLICIT NONE

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description

    CHARACTER (len=*), INTENT(in)     :: path_emiss_file         !< filename with path for EMISS raw data

    ! local variables
    INTEGER (KIND=i4)                 :: time_index, &             !< the index of the time (month) to read in
         &                               ie, &    !< index value for longitude
         &                               je, &    !< index value for latitude
         &                               ke, &    !< counter
         &                               start_cell_id, &  !< ID of starting cell for ICON search
         &                               i,j,k, &  !< counter
         &                               i1, i2, & 
         &                               point_reg_lon_index, &           !< longitude index of point for regular lon-lat grid
         &                               point_reg_lat_index, &           !< latitude index of point for regular lon-lat grid
         &                               nlon_reg, &  !< number of columns
         &                               nlat_reg, &  !< number of rows
         &                               western_column, &      !< the index of the western_column of raw data 
         &                               eastern_column, &      !< the index of the eastern_column of raw data 
         &                               northern_row, &        !< the index of the northern_row of raw data 
         &                               southern_row, &        !< the index of the southern_row of raw data 
         &                               map_ie(emiss_raw_data_grid%nlon_reg, emiss_raw_data_grid%nlat_reg), & 
         &                               map_je(emiss_raw_data_grid%nlon_reg, emiss_raw_data_grid%nlat_reg), & 
         &                               map_ke(emiss_raw_data_grid%nlon_reg, emiss_raw_data_grid%nlat_reg), & 
         &                               ncid_emiss, & 
         &                               row_index, &  !< counter for EMISS data row
         &                               column_index, &  !< counter for EMISS data column
         &                               northern_bound_index, &  !< northern boundary for EMISS data to read for COSMO grid domain
         &                               southern_bound_index !< southern boundary for EMISS data to read for COSMO grid domain
                                      
    REAL (KIND=wp)                    :: default_value, & 
         &                               northern_bound, &  !< northern boundary for EMISS data to read for COSMO grid domain
         &                               southern_bound, &  !< southern boundary for EMISS data to read for COSMO grid domain
         &                               emiss_sum(1:tg%ie,1:tg%je,1:tg%ke), &  !< field of target grid with sum of EMISS values
         &                               point_lon_geo, &        !< longitude coordinate in geographical system of input point 
         &                               point_lat_geo, &        !< latitude coordinate in geographical system of input point
         &                               point_lon, point_lat, & 
         &                               bwlon, &  !< weight for bilinear interpolation
         &                               bwlat, &  !< weight for bilinear interpolation
         &                               emiss_point_sw, &        !< value of the EMISS raw data pixel south west
         &                               emiss_point_se, &        !< value of the EMISS raw data pixel south east
         &                               emiss_point_ne, &        !< value of the EMISS raw data pixel north east
         &                               emiss_point_nw, &        !< value of the EMISS raw data pixel north west
         &                               target_value, & 
         &                               emiss_raw_data(emiss_raw_data_grid%nlon_reg, emiss_raw_data_grid%nlat_reg)

    ! global data flag
    LOGICAL                           :: gldata=.TRUE. ! EMISS data are global

    CALL logging%info('Enter routine: agg_emiss_data_to_target_grid')

    default_value = -1.
    emiss_field = default_value
    emiss_field_mom = default_value

    IF (tg%igrid_type == igrid_cosmo) THEN

      ! determine northern and southern boundary for EMISS data to read for COSMO grid domain
      northern_bound =  MAXVAL(lat_geo) + 1.* ABS(emiss_raw_data_grid%dlat_reg) 
      ! One row of EMISS data north of northern boundary of COSMO grid domain
      northern_bound = MIN(90._wp, northern_bound)                                 ! Check for the poles
        
      northern_bound_index = NINT(ABS((northern_bound - emiss_raw_data_grid%start_lat_reg)/&
                                   ABS(emiss_raw_data_grid%dlat_reg))) +1  ! calculate index for regular lon-lat EMISS grid
      
      IF (northern_bound_index < 1) THEN 
        northern_bound_index = 1 !< check for bounds
      ELSE IF (northern_bound_index > emiss_raw_data_grid%nlat_reg) THEN
        northern_bound_index=emiss_raw_data_grid%nlat_reg
      ENDIF
      southern_bound = MINVAL(lat_geo) - 1.* ABS(emiss_raw_data_grid%dlat_reg)
      ! One row of EMISS data south of southern boundary of COSMO grid domain
      southern_bound = MAX(-90._wp, southern_bound)                               ! Check for the poles

      southern_bound_index = NINT(ABS((southern_bound - emiss_raw_data_grid%start_lat_reg)/&
                              ABS(emiss_raw_data_grid%dlat_reg))) + 2  ! calculate index for regular lon-lat EMISS grid
      IF (southern_bound_index < 1) THEN 
        southern_bound_index = 1 !< check for bounds
      ELSE IF (southern_bound_index > emiss_raw_data_grid%nlat_reg) THEN
        southern_bound_index=emiss_raw_data_grid%nlat_reg
      ENDIf

    END IF

    IF (tg%igrid_type == igrid_icon) THEN
      northern_bound_index=1
      southern_bound_index=emiss_raw_data_grid%nlat_reg
    END IF

    nlon_reg = emiss_raw_data_grid%nlon_reg
    nlat_reg = emiss_raw_data_grid%nlat_reg
    start_cell_id = 1

    ! open netcdf file with EMISS data
    CALL open_netcdf_EMISS_data(path_emiss_file, &
                               ncid_emiss)

    ! read in EMISS data row by row and assign EMISS raw data pixel to COSMO grid
    ! start loop over EMISS raw data
    time_loop: DO time_index=1,12
      no_raw_data_pixel = 0 ! set count to 0
      emiss_sum = 0.  ! set sum to 0

      data_rows: DO row_index=northern_bound_index,southern_bound_index
                    ! get EMISS raw data row
                    
        CALL get_one_row_EMISS_data(ncid_emiss,      &
                  &                 nlon_reg,      &
                  &                 nlat_reg,     &  
                  &                 ntime_emiss,     &
                  &                 row_index,           &
                  &                 time_index,         &
                  &                 emiss_field_row)

        ! store emiss data for subsequent filling algorithm
        emiss_raw_data(1:emiss_raw_data_grid%nlon_reg, row_index) = emiss_field_row(1:emiss_raw_data_grid%nlon_reg)/1.E3 

        column: DO column_index=1, emiss_raw_data_grid%nlon_reg

          IF (time_index == 1) THEN
            point_lon = lon_emiss(column_index)
            point_lat = lat_emiss(row_index)

            ! Reset start cell when entering a new row or when the previous data point was outside
            ! the model domain
            IF (tg%igrid_type == igrid_icon .AND. (column_index == 1 .OR. start_cell_id == 0)) THEN
              i1 = NINT(point_lon*search_res)
              i2 = NINT(point_lat*search_res)
              start_cell_id = tg%search_index(i1,i2)
              IF (start_cell_id == 0) EXIT column ! in this case, the whole row is empty
            ENDIF

            CALL find_nearest_target_grid_element( point_lon, &
                 &                                 point_lat, &
                 &                                 tg,        &
                 &                                 start_cell_id, &
                 &                                 ie,      &
                 &                                 je,      &
                 &                                 ke)

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
            emiss_sum(ie,je,ke)  = emiss_sum(ie,je,ke) + emiss_field_row(column_index)/1.E3 ! sum data values
          ENDIF
        ENDDO column
      END DO data_rows

      DO k=1, tg%ke
        DO j=1, tg%je
          DO i=1, tg%ie
            IF (no_raw_data_pixel(i,j,k) /= 0) THEN 
              emiss_field(i,j,k) = emiss_sum(i,j,k) / REAL(no_raw_data_pixel(i,j,k),wp)   ! calculate arithmetic mean
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      DO k=1, tg%ke
        DO j=1, tg%je
          DO i=1, tg%ie
            IF (no_raw_data_pixel(i,j,k) == 0) THEN 
              point_lon_geo = lon_geo(i,j,k) 
              point_lat_geo = lat_geo(i,j,k)

              ! get four surrounding raw data indices
!DIR$ NOINLINE
              CALL  get_4_surrounding_raw_data_indices( emiss_raw_data_grid, &
                   &                                    lon_emiss,           &
                   &                                    lat_emiss,           &
                   &                                    gldata,             &
                   &                                    point_lon_geo,      &
                   &                                    point_lat_geo,      &
                   &                                    western_column,     &
                   &                                    eastern_column,     &
                   &                                    northern_row,       &
                   &                                    southern_row)
              
              target_value = -999.
              
              IF ( (western_column /= 0) .AND. &
                   (eastern_column /= 0) .AND. &
                   (northern_row /= 0)   .AND. &
                   (southern_row /= 0)        ) THEN
              
                ! get EMISS data for the pixel south west
                point_reg_lon_index = western_column
                point_reg_lat_index = southern_row
                CALL  get_pixel_EMISS_data( ncid_emiss,           &
                                           nlon_reg,           & 
                                           nlat_reg,           &
                                           ntime_emiss,          &
                                           point_reg_lon_index,     &
                                           point_reg_lat_index,     &
                                           time_index,          &
                                           emiss_point_sw)
              
                ! get EMISS data for the pixel south east
                point_reg_lon_index = eastern_column
                point_reg_lat_index = southern_row
                CALL  get_pixel_EMISS_data( ncid_emiss,           &
                                           nlon_reg,           & 
                                           nlat_reg,           &
                                           ntime_emiss,          &
                                           point_reg_lon_index,     &
                                           point_reg_lat_index,     &
                                           time_index,          &
                                           emiss_point_se)
              
                ! get EMISS data for the pixel north east
                point_reg_lon_index = eastern_column
                point_reg_lat_index = northern_row
                CALL  get_pixel_EMISS_data( ncid_emiss,           &
                                           nlon_reg,           & 
                                           nlat_reg,           &
                                           ntime_emiss,          &
                                           point_reg_lon_index,     &
                                           point_reg_lat_index,     &
                                           time_index,          &
                                           emiss_point_ne)
              
                ! get EMISS data for the pixel north west
                point_reg_lon_index = western_column
                point_reg_lat_index = northern_row
                CALL  get_pixel_EMISS_data( ncid_emiss,           &
                                           nlon_reg,           & 
                                           nlat_reg,           &
                                           ntime_emiss,          &
                                           point_reg_lon_index,     &
                                           point_reg_lat_index,     &
                                           time_index,          &
                                           emiss_point_nw)
              
                ! calculate weight for bilinear interpolation
                CALL calc_weight_bilinear_interpol(point_lon_geo, &
                                                   point_lat_geo, &
                                                   lon_emiss(western_column),      &
                                                   lon_emiss(eastern_column),      &
                                                   lat_emiss(northern_row),     &
                                                   lat_emiss(southern_row),     &
                                                   bwlon,         &
                                                   bwlat)
              
                emiss_point_sw = emiss_raw_data(western_column, southern_row) 
                emiss_point_se = emiss_raw_data(eastern_column, southern_row) 
                emiss_point_ne = emiss_raw_data(eastern_column, northern_row) 
                emiss_point_nw = emiss_raw_data(western_column, northern_row) 
              
              
                ! perform the interpolation
                target_value = calc_value_bilinear_interpol(bwlon, bwlat, &
                     &                                      emiss_point_sw, emiss_point_se, emiss_point_ne, emiss_point_nw)
              
                IF (target_value < 0.)THEN 
                  WRITE(message_text,*)'bwlon: ', bwlon, 'bwlat: ', bwlat, &
                       &               'emiss_point_sw: ', emiss_point_sw, 'emiss_point_se: ',&
                       &                emiss_point_se, 'emiss_point_ne: ', emiss_point_ne,&
                       &               'emiss_point_nw: ', emiss_point_nw
              
                  CALL logging%warning(message_text)
                  WRITE(message_text,*)'Target value < 0 '
                  CALL logging%warning(message_text)
                ENDIF
              
                emiss_field(i,j,k) =target_value
              ELSE ! grid element outside target grid
                emiss_field(i,j,k) = default_value
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      WRITE(message_text,*)'emiss_field determined for time_index: ', time_index,&
           &               'with MAXVAL of EMISS:', MAXVAL(emiss_field)
      CALL logging%info(message_text)
      
      emiss_field_mom(:,:,:,time_index) = emiss_field(:,:,:)

    END DO time_loop
   
    CALL  close_netcdf_EMISS_data(ncid_emiss)
   
    ! calculate EMISS Max and EMISS ratio
    emiss_max = MAXVAL(emiss_field_mom,4) ! maximum for the "time" dimension
    
    emiss_ratio_mom = default_value
    
    DO k=1, tg%ke
      DO j=1, tg%je
        DO i=1, tg%ie
          IF(emiss_max(i,j,k) /= 0.) THEN
            DO time_index=1,12
              emiss_ratio_mom(i,j,k,time_index) = emiss_field_mom(i,j,k,time_index) / emiss_max(i,j,k)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    
    CALL logging%info('Exit routine: agg_emiss_data_to_target_grid')

  END SUBROUTINE agg_emiss_data_to_target_grid

END MODULE mo_agg_emiss
