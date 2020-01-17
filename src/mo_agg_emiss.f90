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

  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i4
  USE mo_kind, ONLY: i8

  USE mo_grid_structures, ONLY: target_grid_def
  USE mo_grid_structures, ONLY: igrid_icon
  USE mo_grid_structures, ONLY: igrid_cosmo

  USE mo_search_target_grid, ONLY: find_nearest_target_grid_element

IMPLICIT NONE

PRIVATE
!
PUBLIC :: agg_emiss_data_to_target_grid

    CONTAINS

    !> Subroutine to aggregate EMISS data to target grid
    SUBROUTINE agg_emiss_data_to_target_grid(tg, path_emiss_file)

       USE mo_emiss_data, ONLY: emiss_raw_data_grid, &
                               emiss_field_row, &
                               lon_emiss, &
                               lat_emiss, &
                               ntime_emiss
                               
       USE mo_emiss_tg_fields, ONLY: emiss_field, &
                                    emiss_max, &
                                    emiss_field_mom, &
                                    emiss_ratio_mom

       USE mo_emiss_routines, ONLY: open_netcdf_EMISS_data, &
                                   close_netcdf_EMISS_data, &
                                   get_one_row_EMISS_data, &
                                   get_pixel_EMISS_data

       USE mo_target_grid_data, ONLY: lon_geo, &
                                      lat_geo, &
                                      no_raw_data_pixel
        
    USE mo_target_grid_data, ONLY: search_res !< resolution of ICON grid search index list

    USE mo_bilinterpol, ONLY: get_4_surrounding_raw_data_indices, &
                              calc_weight_bilinear_interpol, &
                              calc_value_bilinear_interpol





       IMPLICIT NONE

       TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description

       CHARACTER (len=*), INTENT(in) :: path_emiss_file         !< filename with path for EMISS raw data


      ! local variables

       INTEGER (KIND=i4)    :: time_index            !< the index of the time (month) to read in

       INTEGER :: ncid_emiss

       REAL (KIND=wp) :: default_value

    INTEGER (KIND=i8) :: ie   !< index value for longitude
    INTEGER (KIND=i8)  :: je   !< index value for latitude
    INTEGER (KIND=i8)  :: ke   !< counter
    INTEGER (KIND=i8) :: start_cell_id !< ID of starting cell for ICON search
    INTEGER (KIND=i8) :: i,j,k !< counter
    INTEGER (KIND=i8) :: i1, i2

    INTEGER :: row_index !< counter for EMISS data row
    INTEGER :: column_index !< counter for EMISS data column
    REAL (KIND=wp) :: northern_bound !< northern boundary for EMISS data to read for COSMO grid domain
    REAL (KIND=wp) :: southern_bound !< southern boundary for EMISS data to read for COSMO grid domain

    INTEGER :: northern_bound_index !< northern boundary for EMISS data to read for COSMO grid domain
    INTEGER :: southern_bound_index !< southern boundary for EMISS data to read for COSMO grid domain


     REAL (KIND=wp) ::  emiss_sum(1:tg%ie,1:tg%je,1:tg%ke) !< field of target grid with sum of EMISS values

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

    REAL (KIND=wp)   :: emiss_point_sw       !< value of the EMISS raw data pixel south west
    REAL (KIND=wp)   :: emiss_point_se       !< value of the EMISS raw data pixel south east
    REAL (KIND=wp)   :: emiss_point_ne       !< value of the EMISS raw data pixel north east
    REAL (KIND=wp)   :: emiss_point_nw       !< value of the EMISS raw data pixel north west

    REAL (KIND=wp)   :: target_value


    ! matrix to save search results

    INTEGER (KIND=i8) :: map_ie(emiss_raw_data_grid%nlon_reg, emiss_raw_data_grid%nlat_reg)
    INTEGER (KIND=i8) :: map_je(emiss_raw_data_grid%nlon_reg, emiss_raw_data_grid%nlat_reg)
    INTEGER (KIND=i8) :: map_ke(emiss_raw_data_grid%nlon_reg, emiss_raw_data_grid%nlat_reg)

    ! buffer for emiss data for one month
    REAL (KIND=wp)   :: emiss_raw_data(emiss_raw_data_grid%nlon_reg, emiss_raw_data_grid%nlat_reg)

    ! global data flag
    LOGICAL :: gldata=.TRUE. ! EMISS data are global





    default_value = -1.
    emiss_field = default_value
    emiss_field_mom = default_value
    
    !HA debug:
    print *,' emiss_raw_data_grid: ', emiss_raw_data_grid
    print *,' emiss_raw_data_grid%start_lat_reg: ', emiss_raw_data_grid%start_lat_reg
    print *,' emiss_raw_data_grid%dlat_reg: ', emiss_raw_data_grid%dlat_reg
    print *,' emiss_raw_data_grid%nlat_reg: ', emiss_raw_data_grid%nlat_reg

     IF (tg%igrid_type == igrid_cosmo) THEN

    
    ! determine northern and southern boundary for EMISS data to read for COSMO grid domain
    northern_bound =  MAXVAL(lat_geo) + 1.* ABS(emiss_raw_data_grid%dlat_reg) 
! One row of EMISS data north of northern boundary of COSMO grid domain
    northern_bound = MIN(90._wp, northern_bound)                                 ! Check for the poles
      
    northern_bound_index = NINT(ABS((northern_bound - emiss_raw_data_grid%start_lat_reg)/&
                                 ABS(emiss_raw_data_grid%dlat_reg))) +1  ! calculate index for regular lon-lat EMISS grid
    
  if (northern_bound_index < 1) then 
        northern_bound_index = 1 !< check for bounds
    else if (northern_bound_index > emiss_raw_data_grid%nlat_reg) then
        northern_bound_index=emiss_raw_data_grid%nlat_reg
    endif
    southern_bound = MINVAL(lat_geo) - 1.* ABS(emiss_raw_data_grid%dlat_reg)
! One row of EMISS data south of southern boundary of COSMO grid domain
    southern_bound = MAX(-90._wp, southern_bound)                               ! Check for the poles

    southern_bound_index = NINT(ABS((southern_bound - emiss_raw_data_grid%start_lat_reg)/&
                            ABS(emiss_raw_data_grid%dlat_reg))) + 2  ! calculate index for regular lon-lat EMISS grid
! DWD  ABS(emiss_raw_data_grid%dlat_reg))) +1 
    if (southern_bound_index < 1) then 
        southern_bound_index = 1 !< check for bounds
    else if (southern_bound_index > emiss_raw_data_grid%nlat_reg) then
        southern_bound_index=emiss_raw_data_grid%nlat_reg
    endif

END IF
IF (tg%igrid_type == igrid_icon) THEN
northern_bound_index=1
southern_bound_index=emiss_raw_data_grid%nlat_reg
END IF

    !HA debug:
    print *,' MAXVAL(lat_geo): ', MAXVAL(lat_geo)
    print *,' MINVAL(lat_geo): ', MINVAL(lat_geo)
    print *, 'northern_bound: ', northern_bound
    print *, 'southern_bound: ', southern_bound

    print *, 'northern_bound_index: ', northern_bound_index
    print *, 'southern_bound_index: ', southern_bound_index
    print *,'lat_emiss(northern_bound_index): ', lat_emiss(northern_bound_index)
    print *,'lat_emiss(southern_bound_index): ', lat_emiss(southern_bound_index)


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
                                nlon_reg,      &
                                nlat_reg,     &  
                                ntime_emiss,     &
                                row_index,           &
                                time_index,         &
                                emiss_field_row)

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
                        emiss_sum(ie,je,ke)  = emiss_sum(ie,je,ke) + emiss_field_row(column_index)/1.E3 ! sum data values
                    ENDIF


                  ENDDO column

    END DO data_rows

    print *,'MAXVAL of EMISS_sum is ', MAXVAL(emiss_sum)

    !HA debug:
    print *,'MAXVAL(no_raw_data_pixel): ',MAXVAL(no_raw_data_pixel)
    print *,'MINVAL(no_raw_data_pixel): ',MINVAL(no_raw_data_pixel)

    !WHERE (no_raw_data_pixel > 0)
    !    emiss_field = emiss_sum / no_raw_data_pixel   ! calculate arithmetic mean
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
      !   print *,'emiss_sum(i,j,k): ', emiss_sum(i,j,k)
      !   PRINT *,'SIZE(emiss_field): ', SIZE(emiss_field)
      !   PRINT *,'SHAPE(emiss_field): ', SHAPE(emiss_field)
      !   
      !   print *,'emiss_field(i,j,k): ', emiss_field(i,j,k)
          emiss_field(i,j,k) = emiss_sum(i,j,k) / REAL(no_raw_data_pixel(i,j,k),wp)   ! calculate arithmetic mean
      !    PRINT *,' emiss_field(i,j,k) = emiss_sum(i,j,k) / REAL(no_raw_data_pixel(i,j,k))'

       ENDIF

     ENDDO
     ENDDO
     ENDDO

 print *,'MAXVAL/MINVAL of EMISS_field for no_raw_data_pixel >0 is ', MAXVAL(emiss_field), MINVAL(emiss_field)

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
         CALL  get_4_surrounding_raw_data_indices(   emiss_raw_data_grid, &
                                                     lon_emiss,           &
                                                     lat_emiss,           &
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



         ! the weights are bwlon and bwlat
       !  PRINT *,'bwlon, bwlat', bwlon,bwlat

        ! perform the interpolation

         target_value = calc_value_bilinear_interpol(bwlon, bwlat, &
                                           emiss_point_sw, emiss_point_se, emiss_point_ne, emiss_point_nw)

         if (target_value < 0.) print*,'Caution target_value < 0: ', bwlon, bwlat, &
                                           emiss_point_sw, emiss_point_se, emiss_point_ne, emiss_point_nw


 !        emiss_field(i,j,k) = MAX(0._wp,target_value)
         emiss_field(i,j,k) =target_value
       ELSE ! grid element outside target grid
         emiss_field(i,j,k) = default_value
      ENDIF

       !---------------------------------------------------------------------------------------------------------------------
   ENDIF

     ENDDO
     ENDDO
     ENDDO

     print *,'emiss_field determined'
     print *,'time_index:', time_index
     print *,'MAXVAL of EMISS is ', MAXVAL(emiss_field)

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

   ! print *,' subroutine aggregate_emiss_raw_data_to_cosmo_grid: MAXVAL(cosmo_emiss_field): ',MAXVAL(cosmo_emiss_field)


       END SUBROUTINE agg_emiss_data_to_target_grid


END MODULE mo_agg_emiss
