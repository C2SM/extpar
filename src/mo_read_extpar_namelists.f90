!+ Fortran modules with namelist definitions for the external parameters software extpar
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_1         2011/01/20 Hermann Asensio
!  small bug fixes accroding to Fortran compiler warnings
! V1_3         2011/04/19 Hermann Asensio
! introduce Globcover 2009 land use data set for external parameters
! add support for GRIB1 and GRIB2
! V1_8         2013-03-12 Frank Brenner
!  introduced MODIS albedo dataset(s) as new external parameter(s)         
! V1_11        2013/04/16 Juergen Helmert
!  Adaptions for using special points and external land-sea-mask
! V1_14        2014-07-18 Juergen Helmert
!  Combined COSMO Release
! V2_1         2015-01-12 Juergen Helmert
!  Bugfix correction covers CSCS SVN r5907-r6359
! V2_3         2015-05-18 Juergen Helmert
!  Change tile_mode switch to integer         
! V2_10        2018-02-19 Juergen Helmert 
!  lsubtract_mean_slope, ERA-I surface temp for land points         
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran modules with namelist definitions for the external parameters software extpar
!> and input routines
MODULE mo_read_extpar_namelists

  USE mo_kind, ONLY: wp, i4, i8
  USE mo_logging
  USE mo_io_units, ONLY: filename_max

  PUBLIC :: read_namelists_extpar_grid_def
  PUBLIC :: read_namelists_extpar_check_icon
  PUBLIC :: read_namelists_extpar_check_cosmo
  PUBLIC :: read_namelists_extpar_special_points

CONTAINS

  !---------------------------------------------------------------------------
  !> subroutine to read namelist for grid settings for EXTPAR
  SUBROUTINE read_namelists_extpar_grid_def(namelist_grid_def, &
       igrid_type, &
       domain_def_namelist, &
       domain_refinement_opt)

    USE mo_utilities_extpar, ONLY: free_un ! function to get free unit number

    CHARACTER (len=*), INTENT(IN) :: namelist_grid_def !< filename with namelists for grid settings for EXTPAR

    INTEGER (KIND=i4), INTENT(OUT)            :: igrid_type       !< target grid type, 1 for ICON, 2 for COSMO, 3 for GME grid
    CHARACTER (len=filename_max), INTENT(OUT) :: domain_def_namelist !< namelist file with domain definition
    CHARACTER (len=filename_max),OPTIONAL, INTENT(OUT) :: domain_refinement_opt   
    !< namelist file with domain refinement defintion (e.g. for the ICON grid)

    ! local variables
    CHARACTER (len=filename_max) :: domain_refinement

    !> namelist with grid defintion
    NAMELIST /grid_def/ igrid_type, domain_def_namelist, domain_refinement


    INTEGER           :: nuin !< unit number
    INTEGER (KIND=i4) :: ierr !< error flag


    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_grid_def), IOSTAT=ierr)
    READ(nuin, NML=grid_def, IOSTAT=ierr)

    CLOSE(nuin)

    ! If optional argument is present for output, copy the value from the local variable to the output argument variable
    IF (PRESENT(domain_refinement_opt)) domain_refinement_opt = TRIM(domain_refinement)

  END SUBROUTINE read_namelists_extpar_grid_def

  !---------------------------------------------------------------------------
  !> subroutine to read namelist for consitency check settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_check_icon(namelist_file,         &
       grib_output_filename,  &
       grib_sample,           &
       netcdf_output_filename,&
       orography_buffer_file, &
       soil_buffer_file,      &
       lu_buffer_file,        &
       glcc_buffer_file,      &
       flake_buffer_file,     &
       ndvi_buffer_file,      &
       sst_icon_file,         &
       t2m_icon_file,         &
       t_clim_buffer_file,    &
       aot_buffer_file,       &
       alb_buffer_file,       &
       i_lsm_data,            &
       land_sea_mask_file,    &
       lwrite_netcdf,         &
       lwrite_grib,           &
       number_special_points, tile_mode )

    USE mo_utilities_extpar, ONLY: free_un ! function to get free unit number


    CHARACTER (len=*), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings


    CHARACTER (len=filename_max) :: grib_output_filename  !< name for grib output filename
    CHARACTER (len=filename_max) :: grib_sample  !< name for grib sample  (sample to be found in $GRIB_SAMPLES_PATH)
    CHARACTER (len=filename_max) :: netcdf_output_filename!< name for netcdf output filename
    CHARACTER (len=filename_max) :: orography_buffer_file  !< name for orography buffer file
    CHARACTER (len=filename_max) :: soil_buffer_file !< name for soil buffer file
    CHARACTER (len=filename_max) :: lu_buffer_file  !< name for glc2000 buffer file

    CHARACTER (len=filename_max) :: glcc_buffer_file  !< name for glcc buffer file

    CHARACTER (len=filename_max) :: flake_buffer_file  !< name for flake buffer file

    CHARACTER (len=filename_max) :: ndvi_buffer_file  !< name for ndvi buffer file
    CHARACTER (len=filename_max) :: sst_icon_file  !< name for sst file
    CHARACTER (len=filename_max) :: t2m_icon_file  !< name for sst file
    CHARACTER (len=filename_max) :: t_clim_buffer_file  !< name for t_clim buffer file
    CHARACTER (len=filename_max) :: aot_buffer_file  !< name for aot buffer file
    CHARACTER (len=filename_max) :: alb_buffer_file  !< name for albedo buffer file
    CHARACTER (len=filename_max) :: land_sea_mask_file  !< name for land-sea mask file
    INTEGER                      :: number_special_points, i_lsm_data
    INTEGER                      :: tile_mode
    LOGICAL                      :: lwrite_netcdf, lwrite_grib

    !> namelist with filenames for output of soil data
    NAMELIST /extpar_consistency_check_io/ grib_output_filename, &
         grib_sample, &
         netcdf_output_filename, &
         orography_buffer_file, &
         soil_buffer_file, &
         lu_buffer_file, &
         glcc_buffer_file, &
         flake_buffer_file, &
         ndvi_buffer_file, &
         sst_icon_file, &
         t2m_icon_file, &
         t_clim_buffer_file, &
         aot_buffer_file, &
         alb_buffer_file, &
         i_lsm_data, &
         land_sea_mask_file,&
         lwrite_netcdf, &
         lwrite_grib, &
         number_special_points, tile_mode


    INTEGER           :: nuin !< unit number
    INTEGER (KIND=i4) :: ierr !< error flag


    nuin = free_un()  ! functioin free_un returns free Fortran unit number

    tile_mode = 0
    lwrite_netcdf = .TRUE.
    lwrite_grib   = .FALSE.

    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)

    READ(nuin, NML=extpar_consistency_check_io, IOSTAT=ierr)

    CLOSE(nuin)

    IF (lwrite_grib) THEN
      CALL logging%info('Direct Grib output is not supported anymore, but has been moved to an post-processing step!', __FILE__, __LINE__)
      lwrite_grib=.FALSE.
    END IF

    print*, "soil_buffer_file = ", TRIM(soil_buffer_file)
    print*, "ndvi_buffer_file = ", TRIM(ndvi_buffer_file)
    print*, "sst_icon_file    = ", TRIM(sst_icon_file)
    print*, "number_special_points, tile_mode ", number_special_points, tile_mode

  END SUBROUTINE read_namelists_extpar_check_icon

  !---------------------------------------------------------------------------
  !> subroutine to read namelist for consitency check settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_check_cosmo(namelist_file,         &
       grib_output_filename,  &
       grib_sample,           &
       netcdf_output_filename,&
       orography_buffer_file, &
       soil_buffer_file,      &
       lu_buffer_file,        &
       glcc_buffer_file,      &
       flake_buffer_file,     &
       ndvi_buffer_file,      &
       t_clim_buffer_file,    &
       aot_buffer_file,       &
       alb_buffer_file,       &
       i_lsm_data,            &
       land_sea_mask_file,    &
       lwrite_netcdf,         &
       lwrite_grib,           &
       number_special_points, &
       tile_mode,             &
       lflake_correction)

    USE mo_utilities_extpar, ONLY: free_un ! function to get free unit number


    CHARACTER (len=*), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings


    CHARACTER (len=filename_max) :: grib_output_filename  !< name for grib output filename
    CHARACTER (len=filename_max) :: grib_sample  !< name for grib sample  (sample to be found in $GRIB_SAMPLES_PATH)
    CHARACTER (len=filename_max) :: netcdf_output_filename!< name for netcdf output filename
    CHARACTER (len=filename_max) :: orography_buffer_file  !< name for orography buffer file
    CHARACTER (len=filename_max) :: soil_buffer_file !< name for soil buffer file
    CHARACTER (len=filename_max) :: lu_buffer_file  !< name for glc2000 buffer file

    CHARACTER (len=filename_max) :: glcc_buffer_file  !< name for glcc buffer file

    CHARACTER (len=filename_max) :: flake_buffer_file  !< name for flake buffer file

    CHARACTER (len=filename_max) :: ndvi_buffer_file  !< name for ndvi buffer file
    CHARACTER (len=filename_max) :: t_clim_buffer_file  !< name for t_clim buffer file
    CHARACTER (len=filename_max) :: aot_buffer_file  !< name for aot buffer file
    CHARACTER (len=filename_max) :: alb_buffer_file  !< name for albedo buffer file
    CHARACTER (len=filename_max) :: land_sea_mask_file  !< name for land-sea mask file
    INTEGER                      :: number_special_points, i_lsm_data
    INTEGER                      :: tile_mode
    LOGICAL                      :: lwrite_netcdf, lwrite_grib, lflake_correction

    !> namelist with filenames for output of soil data
    NAMELIST /extpar_consistency_check_io/ grib_output_filename, &
         grib_sample, &
         netcdf_output_filename, &
         orography_buffer_file, &
         soil_buffer_file, &
         lu_buffer_file, &
         glcc_buffer_file, &
         flake_buffer_file, &
         ndvi_buffer_file, &
         t_clim_buffer_file, &
         aot_buffer_file, &
         alb_buffer_file, &
         i_lsm_data, &
         land_sea_mask_file,&
         lwrite_netcdf, &
         lwrite_grib, &
         tile_mode, &
         number_special_points, &
         lflake_correction


    INTEGER           :: nuin !< unit number
    INTEGER (KIND=i4) :: ierr !< error flag


    nuin = free_un()  ! functioin free_un returns free Fortran unit number

    lwrite_netcdf = .TRUE.
    lwrite_grib   = .FALSE.
    lflake_correction = .TRUE.
    tile_mode = 0

    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)

    READ(nuin, NML=extpar_consistency_check_io, IOSTAT=ierr)

    CLOSE(nuin)

    IF (lwrite_grib) THEN
      CALL logging%info('Direct Grib output is not supported anymore, but has been moved to an post-processing step!', __FILE__, __LINE__)
      lwrite_grib=.FALSE.
    END IF

    print*, "soil_buffer_file = ", soil_buffer_file
    print*, "number_special_points, tile_mode ", number_special_points, tile_mode

  END SUBROUTINE read_namelists_extpar_check_cosmo
  !---------------------------------------------------------------------------
  SUBROUTINE read_namelists_extpar_special_points(namelist_file,        &
       lon_geo_sp,           &
       lat_geo_sp,           &
       soiltype_sp,          &
       z0_sp,                &
       rootdp_sp,            &
       plcovmn_sp,           &
       plcovmx_sp,           &
       laimn_sp,             &
       laimx_sp,             &
       for_d_sp,             &
       for_e_sp,             &
       fr_land_sp            )


    USE mo_utilities_extpar, ONLY: free_un, & ! function to get free unit number
         abort_extpar

    CHARACTER (len=*), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings
    ! orography smoothing



    REAL(KIND=wp),    INTENT(OUT) ::              lon_geo_sp,           &
         lat_geo_sp,           &
         soiltype_sp,          &
         z0_sp,                &
         rootdp_sp,            &
         plcovmn_sp,           &
         plcovmx_sp,           &
         laimn_sp,             &
         laimx_sp,             &
         for_d_sp,             &
         for_e_sp,             &
         fr_land_sp

    !> local variables
    INTEGER           :: nuin     !< unit number
    INTEGER (KIND=i4) :: ierr     !< error flag


    !> define the namelist group
    NAMELIST /special_points/ &
         lon_geo_sp, lat_geo_sp, soiltype_sp, z0_sp, rootdp_sp, plcovmn_sp, plcovmx_sp, &
         laimn_sp, laimx_sp,for_d_sp,for_e_sp,fr_land_sp

    !> initialization
    ierr     = 0



    !> read namelist  
    print *,"special points namelist:", TRIM(namelist_file)
    nuin = free_un()  ! function free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL abort_extpar('read_namelists_extpar_special_points: cannot open file')
    ENDIF
    READ(nuin, NML=special_points, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL abort_extpar(&
           'read_namelists_extpar_special_points: cannot read file: please compare number of special points and defined files!')
    ENDIF
    CLOSE(nuin)

    !> check values for consistency

!!$    IF ((ifill_valley < 1) .OR. (ifill_valley > 2)) THEN
!!$      PRINT *,' Warning  *** ifill valley has to be 1 or 2 *** '
!!$      PRINT *,'          *** set ifill valley = 1 (default value)! *** '
!!$      ifill_valley = 1
!!$    ENDIF    




  END SUBROUTINE read_namelists_extpar_special_points

END MODULE mo_read_extpar_namelists


