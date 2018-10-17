!+ Fortran main program to aggregate isa data to a target grid 
!
!
! Description:
! Fortran main program to aggregate isa data to a target grid
!
! Current Code Owner: DWD, <Name of person responsible for this code>
!    <smail, phone, fax and email>
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V2_4         2015-05-21 Hendrik Wouters
!  Initial release based on extpar_landuse_to_buffer.f90 V1_14
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================
!> Fortran main program to aggregate isa data to a target grid
!!
!! @par extpar_isa_to_buffer 
!!
!! 
!! This program aggregates the ISA (Impervious Surface Area) to a given target grid (COSMO/ICON).
!!
!! @author
!!     Hermann Asensio
!!     (DWD)
!!
!!
!!
PROGRAM extpar_isa_to_buffer
  
  USE info_extpar, ONLY: info_print
  USE mo_logging

  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i8

  USE mo_target_grid_data, ONLY: lon_geo, &
    &                            lat_geo
  
  USE mo_target_grid_data, ONLY: tg
  
  USE mo_target_grid_routines, ONLY: init_target_grid

  USE mo_io_units,          ONLY: filename_max

  USE mo_utilities_extpar,  ONLY: abort_extpar

  USE mo_isa_routines, ONLY: read_namelists_extpar_isa

  USE mo_isa_routines, ONLY:  get_isa_tiles_grid

  !USE mo_isa_tg_fields, ONLY :  i_isa_isa, i_isa_glc2000, i_isa_glcc
  USE mo_isa_tg_fields, ONLY: allocate_isa_target_fields, allocate_add_isa_fields
  USE mo_isa_tg_fields, ONLY: isa_field,       &
  &                          isa_tot_npixel



  USE mo_isa_output_nc, ONLY: write_netcdf_buffer_isa
           

  USE mo_isa_routines, ONLY: get_dimension_isa_data, &
    &                            get_lonlat_isa_data

  USE mo_isa_data, ONLY: isa_grid,                &
    &                          lon_isa,                 &
    &                          lat_isa,                 &
    &                          isa_tiles_grid,          &
    &                          ntiles_isa,              &
    &                          max_tiles_isa,                  &
    &                          isa_tiles_lon_min,              &
    &                          isa_tiles_lon_max,              &
    &                          isa_tiles_lat_min,              &
    &                          isa_tiles_lat_max,              &
    &                          nc_tiles_isa,                   &
    &                          allocate_raw_isa_fields, &
    &                          allocate_isa_data,       &
    &                          fill_isa_data,           &
    &                          deallocate_isa_data

  USE mo_isa_data, ONLY : isa_type !_br 14.04.16

 USE mo_agg_isa, ONLY : agg_isa_data_to_target_grid

  IMPLICIT NONE

  CHARACTER(len=filename_max) :: netcdf_filename

  CHARACTER(len=filename_max) :: namelist_grid_def


  CHARACTER(len=filename_max) :: input_isa_namelist_file
  CHARACTER(len=filename_max), ALLOCATABLE:: isa_file(:)

  CHARACTER (len=filename_max) :: raw_data_isa_path        !< path to raw data
  CHARACTER (len=filename_max) :: raw_data_isa_filename(1:max_tiles_isa) !< filename glc2000 raw data


  CHARACTER (len=filename_max) :: isa_buffer_file !< name for glc2000 buffer file

  INTEGER :: i,k !< counter
  INTEGER :: errorcode

  REAL (KIND=wp) :: undefined

  INTEGER :: undef_int


  INTEGER (KIND=i8) :: nlon_isa !< number of grid elements in zonal direction for isa data
  INTEGER (KIND=i8) :: nlat_isa !< number of grid elements in meridional direction for isa data



  !--------------------------------------------------------------------------------------
  CALL initialize_logging("expar_isa_to_buffer.log", stdout_level=debug)  
  CALL info_print ()
  !--------------------------------------------------------------------------------------------------------

  namelist_grid_def = 'INPUT_grid_org'
  CALL init_target_grid(namelist_grid_def)

  CALL allocate_isa_target_fields(tg)
  print *,'Grid defined, isa target fields allocated'

  !------------------------------------------------------------------------------------

  ! get information about isa data

  ! get info on raw data file
  input_isa_namelist_file = 'INPUT_ISA'

  !---------------------------------------------------------------------------
  CALL read_namelists_extpar_isa(input_isa_namelist_file, &
    &                                 isa_type,    & !_br 14.04.16
    &                                 raw_data_isa_path,       &
    &                                 raw_data_isa_filename,   &
    &                                 ntiles_isa,              &
    &                                 isa_buffer_file          )

! >mes
!     print*,input_isa_namelist_file,   raw_data_isa_path,       &
!     &                                 raw_data_isa_filename,   &
!     &                                 ntiles_isa,       &
!     &                                 isa_buffer_file


       
      CALL allocate_isa_data(ntiles_isa)                  ! allocates the data using ntiles
      CALL fill_isa_data(raw_data_isa_path,     &
                          raw_data_isa_filename, &  ! the allocated vectors need to be filled with the respective value.
                          isa_tiles_lon_min, &
                          isa_tiles_lon_max, &    
                          isa_tiles_lat_min, &
                          isa_tiles_lat_max, &
                          nc_tiles_isa)

   print*, 'ISA TILES, LON, LAT (MIN,MAX): ' 
     DO i = 1,ntiles_isa
       WRITE(*,998)  i, isa_tiles_lon_min(i), isa_tiles_lon_max(i), &
                     isa_tiles_lat_min(i), isa_tiles_lat_max(i) 
998    FORMAT(I1,1X,4(F9.4,1X))      
     END DO

   print*, 'MODEL DOMAIN, LON, LAT (MIN,MAX): ' 

       WRITE(*,999)  lon_geo(1,1,1), lon_geo(tg%ie,tg%je,tg%ke), &
                     lat_geo(1,1,1), lat_geo(tg%ie,tg%je,tg%ke)
999    FORMAT(4(F9.4,1X)) 

!      print*, 'lon_min: ', isa_tiles_lon_min
!      print*, 'lon_max: ', isa_tiles_lon_max
!      print*, 'lat_min: ', isa_tiles_lat_min
!      print*, 'lat_max: ', isa_tiles_lat_max

       DO i = 1,ntiles_isa
   IF(isa_tiles_lon_min(i) < lon_geo(1,1,1).AND.isa_tiles_lon_max(i) > lon_geo(tg%ie,tg%je,tg%ke).AND. &
        isa_tiles_lat_min(i) < lat_geo(1,1,1).AND.isa_tiles_lat_max(i) > lat_geo(tg%ie,tg%je,tg%ke)) THEN

      print*,'MODEL DOMAIN COVERED BY ISA TILE ',i

   END IF
       END DO

    ALLOCATE(isa_file(1:ntiles_isa), STAT= errorcode)
    IF(errorcode /= 0) CALL abort_extpar('Cant allocate isa_file')
! <mes
!    print*, 'ntiles_isa: ', ntiles_isa

! >mes
  DO k = 1,ntiles_isa
    isa_file(k) = TRIM(raw_data_isa_path) // TRIM(raw_data_isa_filename(k))
    PRINT *,'isa_file: ', TRIM(isa_file(k))

  END DO
! <mes



      CALL get_dimension_isa_data(nlon_isa, &
        &                                  nlat_isa)
      CALL allocate_raw_isa_fields(nlat_isa,nlon_isa)
      CALL allocate_add_isa_fields(tg)
      CALL get_lonlat_isa_data( &
        &                              nlon_isa, &
        &                              nlat_isa, &
        &                              lon_isa,  &
        &                              lat_isa,  &
        &                              isa_grid)
        !HA debug
        PRINT *,'isa_grid: ',isa_grid

! >mes
      CALL get_isa_tiles_grid(isa_tiles_grid)
      print*,'isa_tiles_grid(1): ', isa_tiles_grid(1)
! <mes



  undefined = 0.0_wp
  PRINT *,'aggregate isa data to target grid'

    CALL agg_isa_data_to_target_grid(isa_file,                &  
    &                                        undefined,            &
    &                                        isa_tiles_grid, &
    &                                        tg,                   &
    &                                        isa_tot_npixel,        &
    &                                        isa_field    )


    PRINT *,'aggregation of isa data done'

  !--------------------------------------------------------------------------------
  ! output
   undefined = -999.0_wp
   undef_int = -999


   netcdf_filename = TRIM(isa_buffer_file)
   print *, 'Land-use buffer filename: ',TRIM(netcdf_filename)

   !print*,  isa_field
   CALL write_netcdf_buffer_isa(TRIM(netcdf_filename),  &
    &                                     tg,         &
    !&                                     i_isa_data, &
    &                                     undefined, &
    &                                     undef_int,   &
    &                                     lon_geo,     &
    &                                     lat_geo, &
    &                                     isa_tot_npixel, &
    &                                     isa_field)


   CALL deallocate_isa_data()

  PRINT *,'============= isa_to_buffer done ==============='

END PROGRAM extpar_isa_to_buffer

