!+ Fortran module to aggregate AHF data to the target grid
!
!
! Description:
! Fortran module to aggregate AHF data to the target grid
!
! Current Code Owner: DWD, <Name of person responsible for this code>
!    <smail, phone, fax and email>
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V2_4         2015-05-21 Hendrik Wouters
!  Initial release based on mo_ndvi_ahf.f90 (V1_14)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================
!> Fortran module to aggregate AHF data to the target grid
!> \author Hermann Asensio
MODULE mo_agg_ahf

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
PUBLIC :: agg_ahf_data_to_target_grid

    CONTAINS

    !> Subroutine to aggregate AHF data to target grid
    SUBROUTINE agg_ahf_data_to_target_grid(tg,undefined, path_ahf_file)

       USE mo_ahf_data, ONLY: ahf_raw_data_grid, &
                               ahf_field_row, &
                               lon_ahf, &
                               lat_ahf
                               
       USE mo_ahf_tg_fields, ONLY: ahf_field




       USE mo_ahf_routines, ONLY: open_netcdf_AHF_data, &
                                   close_netcdf_AHF_data, &
                                   get_one_row_AHF_data, &
                                   get_pixel_AHF_data

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

       CHARACTER (len=*), INTENT(in) :: path_ahf_file         !< filename with path for AHF raw data


      ! local variables


       INTEGER :: ncid_ahf

       REAL (KIND=wp) :: default_value

    INTEGER (KIND=i8) :: ie   !< index value for longitude
    INTEGER (KIND=i8)  :: je   !< index value for latitude
    INTEGER (KIND=i8)  :: ke   !< counter
    INTEGER (KIND=i8) :: start_cell_id !< ID of starting cell for ICON search
    INTEGER (KIND=i8) :: i,j,k !< counter
    INTEGER (KIND=i8) :: i1, i2

    INTEGER :: row_index !< counter for AHF data row
    INTEGER :: column_index !< counter for AHF data column
    REAL (KIND=wp) :: northern_bound !< northern boundary for AHF data to read for COSMO grid domain
    REAL (KIND=wp) :: southern_bound !< southern boundary for AHF data to read for COSMO grid domain

    INTEGER :: northern_bound_index !< northern boundary for AHF data to read for COSMO grid domain
    INTEGER :: southern_bound_index !< southern boundary for AHF data to read for COSMO grid domain


     REAL (KIND=wp) ::  ahf_sum(1:tg%ie,1:tg%je,1:tg%ke) !< field of target grid with sum of AHF values

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

    REAL (KIND=wp)   :: ahf_point_sw       !< value of the AHF raw data pixel south west
    REAL (KIND=wp)   :: ahf_point_se       !< value of the AHF raw data pixel south east
    REAL (KIND=wp)   :: ahf_point_ne       !< value of the AHF raw data pixel north east
    REAL (KIND=wp)   :: ahf_point_nw       !< value of the AHF raw data pixel north west

    REAL (KIND=wp)   :: target_value


    ! matrix to save search results

    INTEGER (KIND=i8) :: map_ie(ahf_raw_data_grid%nlon_reg, ahf_raw_data_grid%nlat_reg)
    INTEGER (KIND=i8) :: map_je(ahf_raw_data_grid%nlon_reg, ahf_raw_data_grid%nlat_reg)
    INTEGER (KIND=i8) :: map_ke(ahf_raw_data_grid%nlon_reg, ahf_raw_data_grid%nlat_reg)

    ! buffer for ahf data 
    REAL (KIND=wp)   :: ahf_raw_data(ahf_raw_data_grid%nlon_reg, ahf_raw_data_grid%nlat_reg)

    ! global data flag
    LOGICAL :: gldata=.TRUE. ! AHF data are global


    ahf_field = undefined ! set cosmo_ahf_field to undefined value at the start

    default_value = 0.
   ! ahf_field = default_value
   ! ahf_field_mom = default_value
    
    !HA debug:
    print *,' ahf_raw_data_grid: ', ahf_raw_data_grid
    print *,' ahf_raw_data_grid%start_lat_reg: ', ahf_raw_data_grid%start_lat_reg
    print *,' ahf_raw_data_grid%dlat_reg: ', ahf_raw_data_grid%dlat_reg
    print *,' ahf_raw_data_grid%nlat_reg: ', ahf_raw_data_grid%nlat_reg
    
    ! determine northern and southern boundary for AHF data to read for COSMO grid domain
    northern_bound =  MAXVAL(lat_geo) + 1.* ABS(ahf_raw_data_grid%dlat_reg) 

! One row of AHF data north of northern boundary of COSMO grid domain
    northern_bound = MIN(90._wp, northern_bound)                                 ! Check for the poles
      

    northern_bound_index = NINT(ABS((northern_bound - ahf_raw_data_grid%start_lat_reg)/&
                                 ABS(ahf_raw_data_grid%dlat_reg))) +1  ! calculate index for regular lon-lat AHF grid

  if (northern_bound_index < 1) then 
        northern_bound_index = 1 !< check for bounds
    else if (northern_bound_index > ahf_raw_data_grid%nlat_reg) then
        northern_bound_index=ahf_raw_data_grid%nlat_reg
    endif
    southern_bound = MINVAL(lat_geo) - 1.* ABS(ahf_raw_data_grid%dlat_reg) 
! One row of AHF data south of southern boundary of COSMO grid domain
    southern_bound = MAX(-90._wp, southern_bound)                               ! Check for the poles


    southern_bound_index = NINT(ABS((southern_bound - ahf_raw_data_grid%start_lat_reg)/&
                                 ABS(ahf_raw_data_grid%dlat_reg))) +1  ! calculate index for regular lon-lat AHF grid

    if (southern_bound_index < 1) then 
        southern_bound_index = 1 !< check for bounds
    else if (southern_bound_index > ahf_raw_data_grid%nlat_reg) then
        southern_bound_index=ahf_raw_data_grid%nlat_reg
    endif

IF (tg%igrid_type == igrid_icon) THEN

northern_bound_index=1
southern_bound_index=ahf_raw_data_grid%nlat_reg

END IF


    !HA debug:
    print *,' MAXVAL(lat_geo): ', MAXVAL(lat_geo)
    print *,' MINVAL(lat_geo): ', MINVAL(lat_geo)
    print *, 'northern_bound: ', northern_bound
    print *, 'southern_bound: ', southern_bound

    print *, 'northern_bound_index: ', northern_bound_index
    print *, 'southern_bound_index: ', southern_bound_index
    print *,'lat_ahf(northern_bound_index): ', lat_ahf(northern_bound_index)
    print *,'lat_ahf(southern_bound_index): ', lat_ahf(southern_bound_index)


    nlon_reg = ahf_raw_data_grid%nlon_reg
    nlat_reg = ahf_raw_data_grid%nlat_reg
    start_cell_id = 1

   ! open netcdf file with AHF data
    CALL open_netcdf_AHF_data(path_ahf_file, &
                               ncid_ahf)

    ! read in AHF data row by row and assign AHF raw data pixel to COSMO grid
    ! start loop over AHF raw data
    no_raw_data_pixel = 0 ! set count to 0
    ahf_sum = 0.  ! set sum to 0

    
    !print *, 'lat_ahf',lat_ahf

    data_rows: DO row_index=northern_bound_index,southern_bound_index
                    ! get AHF raw data row
                     
                    CALL get_one_row_AHF_data(ncid_ahf,      &
                                nlon_reg,      &
                                nlat_reg,     &  
                                row_index,           &
                                ahf_field_row)

                    ! store ahf data for subsequent filling algorithm
                    ahf_raw_data(1:ahf_raw_data_grid%nlon_reg, row_index) = ahf_field_row(1:ahf_raw_data_grid%nlon_reg)


                    column: DO column_index=1, ahf_raw_data_grid%nlon_reg


   
                           point_lon = lon_ahf(column_index)
                           point_lat = lat_ahf(row_index)











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

                    IF ((ie /= 0).AND.(je/=0).AND.(ke/=0))THEN

                        no_raw_data_pixel(ie,je,ke) = no_raw_data_pixel(ie,je,ke) + 1  
                        ! count raw data pixel within COSMO/ICON grid element
                        ahf_sum(ie,je,ke)  = ahf_sum(ie,je,ke) + ahf_field_row(column_index) ! sum data values
                    ENDIF


                  ENDDO column

    END DO data_rows






    !HA debug:
    print *,'MAXVAL(no_raw_data_pixel): ',MAXVAL(no_raw_data_pixel)
    print *,'MINVAL(no_raw_data_pixel): ',MINVAL(no_raw_data_pixel)

    !WHERE (no_raw_data_pixel > 0)
    !    ahf_field = ahf_sum / no_raw_data_pixel   ! calculate arithmetic mean
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
      !   print *,'ahf_sum(i,j,k): ', ahf_sum(i,j,k)
      !   PRINT *,'SIZE(ahf_field): ', SIZE(ahf_field)
      !   PRINT *,'SHAPE(ahf_field): ', SHAPE(ahf_field)
      !   
      !   print *,'ahf_field(i,j,k): ', ahf_field(i,j,k)
          ahf_field(i,j,k) = ahf_sum(i,j,k) / REAL(no_raw_data_pixel(i,j,k),wp)   ! calculate arithmetic mean
      !    PRINT *,' ahf_field(i,j,k) = ahf_sum(i,j,k) / REAL(no_raw_data_pixel(i,j,k))'
     ELSE ! bilinear interpolation
       !--------------------------------------------------------------------------------------------------------------------

       !  print *,'start bilinear interpolation'
         
         point_lon_geo = lon_geo(i,j,k) 
         point_lat_geo = lat_geo(i,j,k)

        ! get four surrounding raw data indices

!DIR$ NOINLINE
         CALL  get_4_surrounding_raw_data_indices(   ahf_raw_data_grid, &
                                                     lon_ahf,           &
                                                     lat_ahf,           &
                                                     gldata,             &
                                                     point_lon_geo,      &
                                                     point_lat_geo,      &
                                                     western_column,     &
                                                     eastern_column,     &
                                                     northern_row,       &
                                                     southern_row)
         !print *,'western_column, eastern_column, northern_row, southern_row'  
         !print *, western_column, eastern_column, northern_row, southern_row  

         IF ( (western_column /= 0) .AND. &
              (eastern_column /= 0) .AND. &
              (northern_row /= 0)   .AND. &
              (southern_row /= 0)        ) THEN

           ahf_point_sw = ahf_raw_data(western_column, southern_row) 
           ahf_point_se = ahf_raw_data(eastern_column, southern_row) 
           ahf_point_ne = ahf_raw_data(eastern_column, northern_row) 
           ahf_point_nw = ahf_raw_data(western_column, northern_row) 


         ! calculate weight for bilinear interpolation

         CALL calc_weight_bilinear_interpol(point_lon_geo, &
                                            point_lat_geo, &
                                            lon_ahf(western_column),      &
                                            lon_ahf(eastern_column),      &
                                            lat_ahf(northern_row),     &
                                            lat_ahf(southern_row),     &
                                            bwlon,         &
                                            bwlat)
         ! the weights are bwlon and bwlat
       !  PRINT *,'bwlon, bwlat', bwlon,bwlat

        ! perform the interpolation

         target_value = calc_value_bilinear_interpol(bwlon, bwlat, &
                                           ahf_point_sw, ahf_point_se, ahf_point_ne, ahf_point_nw)


         ahf_field(i,j,k) = MAX(0._wp,target_value)

       ELSE ! grid element outside target grid
         ahf_field(i,j,k) = default_value

       ENDIF

       !---------------------------------------------------------------------------------------------------------------------
       ENDIF

     ENDDO
     ENDDO
     ENDDO

     print *,'ahf_field determined'

    !ahf_field_mom(:,:,:,time_index) = ahf_field(:,:,:)


    CALL  close_netcdf_AHF_data(ncid_ahf)

      ! calculate AHF Max and AHF ratio

       END SUBROUTINE agg_ahf_data_to_target_grid


END MODULE mo_agg_ahf
