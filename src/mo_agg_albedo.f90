!# Comment from Merge: Different approaches used! Check the code and results!
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
  USE mo_grid_structures, ONLY: igrid_gme


  USE mo_search_ll_grid, ONLY: find_reg_lonlat_grid_element_index, &
    &                          find_rotated_lonlat_grid_element_index

  USE mo_search_target_grid, ONLY: find_nearest_target_grid_element

IMPLICIT NONE

PRIVATE
!
PUBLIC :: agg_alb_data_to_target_grid

    CONTAINS

    !> Subroutine to aggregate albedo data to target grid
    SUBROUTINE agg_alb_data_to_target_grid(tg,undefined, path_alb_file, &
         &                     alb_source, alb_field_mom_d)

       USE mo_albedo_data, ONLY: alb_raw_data_grid, &
                               alb_field_row_mom, &
                               alnid_field_row_mom, &
                               aluvd_field_row_mom, &
                               alb_field_row, &
                               alnid_field_row, &
                               aluvd_field_row, &
                               lon_alb, &
                               lat_alb, &
                               ntime_alb
                               
       USE mo_albedo_tg_fields, ONLY: alb_field, ialb_type, &
                                      alnid_field, &
                                      aluvd_field

       USE mo_albedo_routines, ONLY: open_netcdf_ALB_data, &
                                   close_netcdf_ALB_data, &
                                   get_one_row_ALB_data, &
                                   get_pixel_ALB_data

       USE mo_target_grid_data, ONLY: lon_geo, &
                                      lat_geo, &
                                      no_raw_data_pixel

       USE mo_target_grid_data, ONLY: search_res !< resolution of ICON grid search index list
        
    USE mo_gme_grid, ONLY: gme_grid
    USE mo_gme_grid, ONLY: sync_diamond_edge
    USE mo_gme_grid, ONLY: gme_real_field, gme_int_field
    USE mo_gme_grid, ONLY: cp_buf2gme, cp_gme2buf


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

       CHARACTER (len=*), INTENT(in)  :: path_alb_file         !< filename with path for albedo raw data  
       CHARACTER (len=*), INTENT(in)  :: alb_source    !< albedo variable name inside input file
       REAL (KIND=wp),    INTENT(OUT) :: alb_field_mom_d(:,:,:,:)   !< monthly mean albedo, output variable


      ! local variables

       INTEGER (KIND=i4)    :: time_index            !< the index of the time (month) to read in

       INTEGER :: ncid_alb

       REAL (KIND=wp) :: default_value

    INTEGER (KIND=i8)  :: ie   !< index value for longitude
    INTEGER (KIND=i8)  :: je   !< index value for latitude
    INTEGER (KIND=i8)  :: ke   !< counter
    INTEGER (KIND=i8) :: start_cell_id !< ID of starting cell for ICON search
    INTEGER (KIND=i8)  :: i,j,k !< counter
    INTEGER (KIND=i8) :: i1, i2

    INTEGER :: row_index !< counter for data row
    INTEGER :: column_index !< counter for data column
    REAL (KIND=wp) :: northern_bound !< northern boundary for input data to read for COSMO grid domain
    REAL (KIND=wp) :: southern_bound !< southern boundary for input data to read for COSMO grid domain

    INTEGER :: northern_bound_index !< northern boundary for input data to read for COSMO grid domain
    INTEGER :: southern_bound_index !< southern boundary for input data to read for COSMO grid domain

    REAL (KIND=wp) ::  alb_sum(1:tg%ie,1:tg%je,1:tg%ke) !< field of target grid with sum of albedo values

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

    REAL (KIND=wp)   :: alb_point_sw       !< albedo value of the raw data pixel south west
    REAL (KIND=wp)   :: alb_point_se       !< albedo value of the raw data pixel south east
    REAL (KIND=wp)   :: alb_point_ne       !< albedo value of the raw data pixel north east
    REAL (KIND=wp)   :: alb_point_nw       !< albedo value of the raw data pixel north west

    REAL (KIND=wp)   :: target_value
    REAL (KIND=wp)   :: lon_alt
    REAL (KIND=wp), PARAMETER   :: undef_raw=-0.01_wp


    ! matrix to save search results

    INTEGER (KIND=i8) :: map_ie(alb_raw_data_grid%nlon_reg, alb_raw_data_grid%nlat_reg)
    INTEGER (KIND=i8) :: map_je(alb_raw_data_grid%nlon_reg, alb_raw_data_grid%nlat_reg)
    INTEGER (KIND=i8) :: map_ke(alb_raw_data_grid%nlon_reg, alb_raw_data_grid%nlat_reg)

    ! buffer for albedo raw data for one month
    REAL (KIND=wp)   :: albedo_raw_data(alb_raw_data_grid%nlon_reg, alb_raw_data_grid%nlat_reg)

    ! global data flag
    LOGICAL :: gldata=.TRUE. ! input data is global

    INTEGER (KIND=i4), ALLOCATABLE :: no_valid_raw_data_pixel(:,:,:)

!$   INTEGER :: omp_get_max_threads, omp_get_thread_num, thread_id, num_threads
!$   INTEGER (KIND=i8), ALLOCATABLE :: start_cell_arr(:)
     LOGICAL :: l_new_row

!      num_threads = 1
!      call OMP_SET_NUM_THREADS(num_threads)
!$     print*, ' Number of threads used : ', omp_get_max_threads()

    alb_field = undefined ! set cosmo_alb_field to undefined value at the start
    alb_field_mom_d = undefined
    default_value = 0.07
    
    !HA debug:
    print *,' alb_raw_data_grid: ', alb_raw_data_grid
    print *,' alb_raw_data_grid%start_lat_reg: ', alb_raw_data_grid%start_lat_reg
    print *,' alb_raw_data_grid%dlat_reg: ', alb_raw_data_grid%dlat_reg
    print *,' alb_raw_data_grid%nlat_reg: ', alb_raw_data_grid%nlat_reg
    
    ! determine northern and southern boundary for input data to read for COSMO grid domain
    northern_bound =  MAXVAL(lat_geo) + 1.* alb_raw_data_grid%dlat_reg
 ! One row of albedo data north of northern boundary of COSMO grid domain
    northern_bound = MIN(90._wp, northern_bound)                                 ! Check for the poles
   
   print *, 'northern_bound: ', northern_bound
 
    northern_bound_index = NINT((northern_bound - alb_raw_data_grid%start_lat_reg)/&
                                 alb_raw_data_grid%dlat_reg) +1  ! calculate index for regular lon-lat grid
  
  print *, 'northern_bound_index: ', northern_bound_index

  if (northern_bound_index < 1) then 
        northern_bound_index = 1 !< check for bounds
    else if (northern_bound_index > alb_raw_data_grid%nlat_reg) then
        northern_bound_index=alb_raw_data_grid%nlat_reg
    endif
    southern_bound = MINVAL(lat_geo) - 1.* alb_raw_data_grid%dlat_reg
 ! One row of albedo data south of southern boundary of COSMO grid domain
    southern_bound = MAX(-90._wp, southern_bound)                               ! Check for the poles

    print *, 'southern_bound: ', southern_bound


    southern_bound_index = NINT((southern_bound - alb_raw_data_grid%start_lat_reg)/&
                                 alb_raw_data_grid%dlat_reg) +1  ! calculate index for regular lon-lat grid
    print *, 'southern_bound_index: ', southern_bound_index
    if (southern_bound_index < 1) then 
        southern_bound_index = 1 !< check for bounds
    else if (southern_bound_index > alb_raw_data_grid%nlat_reg) then
        southern_bound_index=alb_raw_data_grid%nlat_reg
    endif

    !HA debug:
    print *,' MAXVAL(lat_geo): ', MAXVAL(lat_geo)
    print *,' MINVAL(lat_geo): ', MINVAL(lat_geo)
    print *, 'northern_bound: ', northern_bound
    print *, 'southern_bound: ', southern_bound

    print *, 'northern_bound_index: ', northern_bound_index
    print *, 'southern_bound_index: ', southern_bound_index
    print *,'lat_alb(northern_bound_index): ', lat_alb(northern_bound_index)
    print *,'lat_alb(southern_bound_index): ', lat_alb(southern_bound_index)


    nlon_reg = alb_raw_data_grid%nlon_reg
    nlat_reg = alb_raw_data_grid%nlat_reg
    start_cell_id = 1

!$   num_threads = omp_get_max_threads()
!$   ALLOCATE(start_cell_arr(num_threads))
!$   start_cell_arr(:) = 1

   ! open netcdf file with albedo data
    print *, 'In mo_agg_albedo.f90'
    print *, TRIM(path_alb_file)
    CALL open_netcdf_ALB_data(path_alb_file, &
                               ncid_alb)

    ALLOCATE(no_valid_raw_data_pixel(tg%ie, tg%je, tg%ke))

    PRINT *, 'mo_agg_albedo, albedo nc input file opened'

    ! read in albedo data row by row and assign albedo raw data pixel to target grid
    ! start loop over albedo raw data
    time_loop: DO time_index=1,ntime_alb
    no_raw_data_pixel       = 0 ! set count to 0
    no_valid_raw_data_pixel = 0 ! set count to 0
    alb_sum = 0.  ! set sum to 0

!    data_rows: DO row_index=northern_bound_index,southern_bound_index
    data_rows: DO row_index=southern_bound_index,northern_bound_index
                 ! get input raw data row
                    
                 CALL get_one_row_ALB_data(ncid_alb,      &
                                           nlon_reg,      &
                                           nlat_reg,      &  
                                           ntime_alb,     &
                                           row_index,     &
                                           time_index,    &
                                           alb_field_row, &
                                           alb_source)

                 IF (time_index == 1) THEN

!$OMP PARALLEL PRIVATE(thread_id,start_cell_id,l_new_row)
!$    thread_id = omp_get_thread_num()+1
!$    start_cell_id = start_cell_arr(thread_id)
      l_new_row = .TRUE.

!$OMP DO PRIVATE(column_index,point_lon,point_lat,i1,i2,ie,je,ke)
                   column: DO column_index=1, alb_raw_data_grid%nlon_reg

                       IF (time_index == 1) THEN
                     point_lon = lon_alb(column_index)
                     point_lat = lat_alb(row_index)

                     ! Reset start cell when entering a new row or when the previous data point was outside
                     ! the model domain
                     IF (tg%igrid_type == igrid_icon .AND. (l_new_row .OR. start_cell_id == 0)) THEN
                       IF (point_lon > 180.) THEN
                         i1 = NINT((point_lon-360.)*search_res)
                       ELSE
                         i1 = NINT(point_lon*search_res)
                       ENDIF
                       i2 = NINT(point_lat*search_res)
                       start_cell_id = tg%search_index(i1,i2)
                       l_new_row = .FALSE.
                       IF (start_cell_id == 0) CYCLE 
                     ENDIF

                     CALL find_nearest_target_grid_element( point_lon, &
                                                            point_lat, &
                                                            tg,        &
                                                            start_cell_id, &
                                                            ie,      &
                                                            je,      &
                                                            ke)
!                         print *, 'FB debug: find_nearest_target_grid_element done '

                     map_ie(column_index, row_index) = ie
                     map_je(column_index, row_index) = je
                     map_ke(column_index, row_index) = ke

                   ENDDO column

!$OMP END DO
!$    start_cell_arr(thread_id) = start_cell_id
!$OMP END PARALLEL
                 ENDIF

                   column2: DO column_index=1, alb_raw_data_grid%nlon_reg

                     ! store albedo data for subsequent filling algorithm
                     albedo_raw_data(column_index, row_index) = alb_field_row(column_index)

                     ! aggregate albedo data on model grid
                     ie = map_ie(column_index, row_index)
                     je = map_je(column_index, row_index)
                     ke = map_ke(column_index, row_index)
                        ENDIF

                     IF ((ie /= 0).AND.(je/=0).AND.(ke/=0))THEN
                       no_raw_data_pixel(ie,je,ke) = no_raw_data_pixel(ie,je,ke) + 1
                       IF (alb_field_row(column_index)>0.02) THEN
                         no_valid_raw_data_pixel(ie,je,ke) = no_valid_raw_data_pixel(ie,je,ke) + 1
            ! count raw data pixel within COSMO grid element
                         alb_sum(ie,je,ke)  = alb_sum(ie,je,ke) + alb_field_row(column_index) ! sum data values
                     ENDIF
                    ENDIF
                  ENDDO column
    END DO data_rows

     SELECT CASE(tg%igrid_type)
     CASE(igrid_gme)  ! in GME grid the diamond edges need to be synrchonized

       ! alb_sum
       CALL cp_buf2gme(tg,gme_grid,alb_sum,gme_real_field)
       CALL sync_diamond_edge(gme_grid, gme_real_field)
       CALL cp_gme2buf(tg,gme_grid,gme_real_field,alb_sum)

      ! no_raw_data_pixel
      CALL cp_buf2gme(tg,gme_grid,no_raw_data_pixel,gme_int_field)
      CALL sync_diamond_edge(gme_grid, gme_int_field)
      CALL cp_gme2buf(tg,gme_grid,gme_int_field,no_raw_data_pixel)
     END SELECT

    !HA debug:
    PRINT *,'MAXVAL(no_raw_data_pixel): ',MAXVAL(no_raw_data_pixel)
    PRINT *,'MINVAL(no_raw_data_pixel): ',MINVAL(no_raw_data_pixel)

    PRINT *,'tg: ',tg%ie, tg%je, tg%ke, tg%minlon, tg%maxlon, tg%minlat, tg%maxlat

! OMP was deactivated here due to race conditions !!!
!OMP PARALLEL   
    DO k=1, tg%ke
    DO j=1, tg%je
!OMP DO PRIVATE(i,point_lon_geo,point_lat_geo,western_column,eastern_column,northern_row,southern_row, &
!OMP            alb_point_sw,alb_point_se,alb_point_nw,alb_point_ne,lon_alt,bwlon,target_value )
    DO i=1, tg%ie

     IF (no_valid_raw_data_pixel(i,j,k) /= 0) THEN 
          alb_field_mom_d(i,j,k,time_index) = alb_sum(i,j,k) / REAL(no_valid_raw_data_pixel(i,j,k))   ! calculate arithmetic mean
!roa sorry but the next line is not correct! Please ask DWD to repair! I just do a quick fix... >
  !   ELSE IF (no_raw_data_pixel(ie,je,ke) == 0) THEN! bilinear interpolation
     ELSEIF (no_valid_raw_data_pixel(i,j,k) == 0) THEN
!roa <
       !--------------------------------------------------------------------------------------------------
       point_lon_geo = lon_geo(i,j,k) 
       point_lat_geo = lat_geo(i,j,k)

       ! get four surrounding raw data indices
       CALL  get_4_surrounding_raw_data_indices(   alb_raw_data_grid, &
                                                   lon_alb,           &
                                                   lat_alb,           &
                                                   gldata,             &
                                                   point_lon_geo,      &
                                                   point_lat_geo,      &
                                                   western_column,     &
                                                   eastern_column,     &
                                                   northern_row,       &
                                                   southern_row)

         IF ( (western_column /= 0) .AND. &
              (eastern_column /= 0) .AND. &
              (northern_row /= 0)   .AND. &
              (southern_row /= 0)        ) THEN

           ! get albedo data for the pixel south west
           point_reg_lon_index = western_column
           point_reg_lat_index = southern_row
           CALL  get_pixel_ALB_data( ncid_alb,           &
                                      nlon_reg,           & 
                                      nlat_reg,           &
                                      ntime_alb,          &
                                      point_reg_lon_index,     &
                                      point_reg_lat_index,     &
                                      time_index,          &
                                      alb_point_sw,    &
                                      alb_source)

           ! get albedo data for the pixel south east
           point_reg_lon_index = eastern_column
           point_reg_lat_index = southern_row
           CALL  get_pixel_ALB_data( ncid_alb,           &
                                      nlon_reg,           & 
                                      nlat_reg,           &
                                      ntime_alb,          &
                                      point_reg_lon_index,     &
                                      point_reg_lat_index,     &
                                      time_index,          &
                                      alb_point_se,      &
                                      alb_source)

           ! get albedo data for the pixel north east
           point_reg_lon_index = eastern_column
           point_reg_lat_index = northern_row
           CALL  get_pixel_ALB_data( ncid_alb,           &
                                      nlon_reg,           & 
                                      nlat_reg,           &
                                      ntime_alb,          &
                                      point_reg_lon_index,     &
                                      point_reg_lat_index,     &
                                      time_index,          &
                                      alb_point_ne,       &
                                      alb_source)

           ! get albedo data for the pixel north west
           point_reg_lon_index = western_column
           point_reg_lat_index = northern_row
           CALL  get_pixel_ALB_data( ncid_alb,           &
                                      nlon_reg,           & 
                                      nlat_reg,           &
                                      ntime_alb,          &
                                      point_reg_lon_index,     &
                                      point_reg_lat_index,     &
                                      time_index,          &
                                      alb_point_nw,       &
                                      alb_source)

           ! calculate weight for bilinear interpolation
           CALL calc_weight_bilinear_interpol(point_lon_geo, &
                                              point_lat_geo, &
                                              lon_alb(western_column),      &
                                              lon_alb(eastern_column),      &
                                              lat_alb(northern_row),     &
                                              lat_alb(southern_row),     &
                                              bwlon,         &
                                              bwlat)
           ! the weights are bwlon and bwlat

        IF (ialb_type == 2) THEN
         IF ((alb_point_ne>undef_raw).OR.(alb_point_nw>undef_raw).OR. &
            (alb_point_se>undef_raw).OR.(alb_point_sw>undef_raw)) THEN
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
           IF (alb_point_sw > 0.02 .AND. alb_point_se > 0.02 .AND. alb_point_ne > 0.02 .AND. alb_point_nw > 0.02) THEN
             target_value = calc_value_bilinear_interpol(bwlon, bwlat, &
                                           alb_point_sw, alb_point_se, alb_point_ne, alb_point_nw)
          IF (target_value < 0.02) THEN
!               print *,'Interpolation gone wrong! ',target_value,alb_point_sw, alb_point_se, alb_point_ne, &
!                        alb_point_nw,bwlon, bwlat 
               target_value = -999.
          ENDIF
           ELSE
             target_value = -999. ! assume missing value - will be fixed later in the cross check
           ENDIF

!         IF (target_value.gt.1) THEN
!           PRINT *, 'alb_field interpolation error: ',i,j,k,target_value
!           PRINT *, alb_point_sw,alb_point_se,alb_point_ne,alb_point_nw
!           target_value = 0.25 * (alb_point_sw + alb_point_se + alb_point_ne + &
!  &                       alb_point_nw)
!           PRINT *, point_lon_geo, point_lat_geo, bwlon, bwlat
!           PRINT *, lon_alb(western_column),lon_alb(eastern_column)
!           PRINT *, 'new value: ',target_value
!         ENDIF
         alb_field_mom_d(i,j,k,time_index) = target_value

         ELSE ! grid element outside target grid
         alb_field_mom_d(i,j,k,time_index) = default_value
         ENDIF
       IF (bwlon>1) THEN !calculation of bwlon gone wrong (eastern border)
         alb_field_mom_d(i,j,k,time_index) = -999.
       ENDIF
       IF (bwlon<0) THEN !calculation of bwlon gone wrong (western border)
         alb_field_mom_d(i,j,k,time_index) = -999.
       ENDIF

       ELSE
       alb_field_mom_d(i,j,k,time_index) = -999. !assume missing value if no valid data point was present, 
                                                 !will be fixed later in the cross check
       ENDIF

     ENDDO !i
!OMP END DO
     ENDDO !j
     ENDDO !k
!OMP END PARALLEL
     print *,'alb_field determined, time_index:', time_index

    alb_field_mom_d(:,:,:,time_index) = alb_field(:,:,:)

    END DO time_loop
    CALL  close_netcdf_ALB_data(ncid_alb)
    DEALLOCATE(no_valid_raw_data_pixel)
  
  END SUBROUTINE agg_alb_data_to_target_grid

END MODULE mo_agg_albedo
