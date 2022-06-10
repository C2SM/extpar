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

  USE mo_logging                
  USE mo_kind,                  ONLY: wp, i4
  USE mo_utilities_extpar,      ONLY: check_input_file 
  USE mo_io_units,              ONLY: filename_max
  USE info_extpar,              ONLY: INFO_CompilerVersion

  IMPLICIT NONE

  PRIVATE
  
  PUBLIC :: read_namelists_extpar_grid_def, &
       &    read_namelists_extpar_check_icon, &
       &    read_namelists_extpar_check_cosmo, &
       &    read_namelists_extpar_special_points

  CONTAINS

  !> subroutine to read namelist for grid settings for EXTPAR
  SUBROUTINE read_namelists_extpar_grid_def(namelist_grid_def, &
       &                                    igrid_type, &
       &                                    domain_def_namelist, &
       &                                    domain_refinement_opt)

    CHARACTER (len=*), INTENT(IN)       :: namelist_grid_def !< filename with namelists for grid settings for EXTPAR
    INTEGER(KIND=i4), INTENT(OUT)                :: igrid_type       !< target grid type, 1 for ICON, 2 for COSMO, 3 GME
    CHARACTER (len=filename_max), INTENT(OUT) &
      &                                 :: domain_def_namelist !< namelist file with domain definition
    CHARACTER (len=*),OPTIONAL, INTENT(OUT)   &
      &                                 :: domain_refinement_opt   

    ! local variables
    CHARACTER (len=filename_max)        :: domain_refinement
    INTEGER                             :: nuin, ierr

    !> namelist with grid defintion
    NAMELIST /grid_def/ igrid_type, domain_def_namelist, domain_refinement

    CALL logging%info('Enter routine: read_namelists_extpar_grid_def')

    OPEN(NEWUNIT=nuin, FILE=TRIM(namelist_grid_def), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_grid_def)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=grid_def, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist grid_def',__FILE__, __LINE__) 
    ENDIF
    CLOSE(nuin)

    ! If optional argument is present for output, copy the value from the local variable to the output argument variable
    IF (PRESENT(domain_refinement_opt)) domain_refinement_opt = TRIM(domain_refinement)

  END SUBROUTINE read_namelists_extpar_grid_def

  !> subroutine to read namelist for consitency check settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_check_icon(namelist_file,         &
       &                                      grib_output_filename,  &
       &                                      grib_sample,           &
       &                                      netcdf_output_filename,&
       &                                      sst_icon_file,         &
       &                                      t2m_icon_file,         &
       &                                      i_lsm_data,            &
       &                                      land_sea_mask_file,    &
       &                                      lwrite_netcdf,         &
       &                                      lwrite_grib,           &
       &                                      number_special_points, &
       &                                      tile_mode,             &
       &                                      l_use_glcc,            &
       &                                      l_use_array_cache      )

    CHARACTER (len=*), INTENT(IN)             :: namelist_file !< filename with namelists for for EXTPAR settings

    CHARACTER (len=filename_max), INTENT(OUT) :: grib_output_filename, &  !< name for grib output filename
         &                                       grib_sample, &  !< name for grib sample  (sample to be found in $GRIB_SAMPLES_PATH)
         &                                       netcdf_output_filename, &!< name for netcdf output filename
         &                                       sst_icon_file,         &
         &                                       t2m_icon_file,         &
         &                                       land_sea_mask_file  !< name for land-sea mask file

    INTEGER(KIND=i4),INTENT(OUT)              :: number_special_points, i_lsm_data, &
         &                                       tile_mode

    LOGICAL,INTENT(OUT)                       :: lwrite_netcdf, lwrite_grib, &
         &                                       l_use_glcc, l_use_array_cache

    !local variables
    INTEGER(KIND=i4)                          :: nuin, ierr

    !> namelist with filenames for output of soil data
    NAMELIST /extpar_consistency_check_io/ grib_output_filename, &
         &                                 grib_sample, &
         &                                 netcdf_output_filename, &
         &                                 sst_icon_file,         &
         &                                 t2m_icon_file,         &
         &                                 i_lsm_data, &
         &                                 land_sea_mask_file,&
         &                                 lwrite_netcdf, &
         &                                 lwrite_grib, &
         &                                 number_special_points, &
         &                                 tile_mode, &
         &                                 l_use_glcc, &
         &                                 l_use_array_cache 


    CALL logging%info('Enter routine: read_namelists_extpar_check_icon')

    sst_icon_file = '_'
    t2m_icon_file = '_'
    tile_mode = 0
    l_use_array_cache  = .FALSE.   ! Might be slower, but required for really high resolution
    lwrite_netcdf = .TRUE.
    lwrite_grib   = .FALSE.
    l_use_glcc    = .TRUE. ! Assume that GLCC land-use data file exists!

    OPEN(NEWUNIT=nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=extpar_consistency_check_io, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist extpar_consistency_check_io',__FILE__, __LINE__) 
    ENDIF

    CLOSE(nuin)

    IF (lwrite_grib) THEN
      CALL logging%warning('Direct Grib output is not supported anymore, but has been moved to a post-processing step!')
      lwrite_grib=.FALSE.
    END IF

    ! only check era-files if specified in namelist
    IF (TRIM(sst_icon_file) /= '_' .AND. TRIM(t2m_icon_file) /= '_') THEN
      CALL check_input_file(TRIM(sst_icon_file), __FILE__, __LINE__)
      CALL check_input_file(TRIM(t2m_icon_file), __FILE__, __LINE__)
    ENDIF

    WRITE(message_text,'(a,i0)') 'No of special points: ', number_special_points
    CALL logging%info(message_text)
    WRITE(message_text,'(a,i0)') 'Tile mode: ',  tile_mode
    CALL logging%info(message_text)

    IF (l_use_array_cache) THEN
      CALL logging%info('array caching (less memory consumption) is activated!')
      
      ! list of compilers not working with array-caching
      IF ( (index(INFO_CompilerVersion,'intel') /= 0) .OR. &
           (index(INFO_CompilerVersion,'INTEL') /= 0) .OR. &
           (index(INFO_CompilerVersion,'Intel') /= 0) ) THEN

        WRITE(message_text,*) 'Array-caching not supported for', &
             &                TRIM(INFO_CompilerVersion),' compiler!'
        CALL logging%error(message_text,__FILE__,__LINE__)
      ENDIF
    ENDIF

    ! Currently we only support GCC for Extpar
    IF ( (index(INFO_CompilerVersion,'gcc') == 0) .AND. &
         (index(INFO_CompilerVersion,'GCC') == 0) .AND. &
         (index(INFO_CompilerVersion,'Gcc') == 0) ) THEN

      WRITE(message_text,*) 'EXTPAR IS NOT TESTED WITH ', &
           &                TRIM(INFO_CompilerVersion),' COMPILER!'
      CALL logging%error(message_text,__FILE__,__LINE__)
    ENDIF
    
    CALL logging%info('Exit routine: read_namelists_extpar_check_icon')

  END SUBROUTINE read_namelists_extpar_check_icon

  !> subroutine to read namelist for consitency check settings for EXTPAR 
  SUBROUTINE read_namelists_extpar_check_cosmo(namelist_file,         &
       &                                       grib_output_filename,  &
       &                                       grib_sample,           &
       &                                       netcdf_output_filename,&
       &                                       i_lsm_data,            &
       &                                       land_sea_mask_file,    &
       &                                       lwrite_netcdf,         &
       &                                       lwrite_grib,           &
       &                                       number_special_points, &
       &                                       tile_mode,             &
       &                                       lflake_correction,     &
       &                                       l_use_array_cache)

    CHARACTER (len=*), INTENT(IN)             :: namelist_file !< filename with namelists for for EXTPAR settings

    CHARACTER (len=filename_max), INTENT(OUT) :: grib_output_filename, &  !< name for grib output filename
         &                                       grib_sample, &  !< name for grib sample  (sample to be found in $GRIB_SAMPLES_PATH)
         &                                       netcdf_output_filename, &!< name for netcdf output filename
         &                                       land_sea_mask_file  !< name for land-sea mask file

    INTEGER(KIND=i4), INTENT(OUT)             :: number_special_points, i_lsm_data, &
         &                                       tile_mode

    LOGICAL, INTENT(OUT)                      :: lwrite_netcdf, lwrite_grib, lflake_correction, &
                                                 l_use_array_cache

    INTEGER(KIND=i4)                          :: nuin, ierr

    !> namelist with filenames for output of soil data
    NAMELIST /extpar_consistency_check_io/ grib_output_filename, &
         grib_sample, &
         netcdf_output_filename, &
         i_lsm_data, &
         land_sea_mask_file,&
         lwrite_netcdf, &
         lwrite_grib, &
         tile_mode, &
         number_special_points, &
         lflake_correction, &
         l_use_array_cache 

    lwrite_netcdf = .TRUE.
    lwrite_grib   = .FALSE.
    lflake_correction = .TRUE.
    tile_mode = 0
    l_use_array_cache  = .FALSE.   ! Might be slower, but required for really high resolution

    CALL logging%info('Enter routine: read_namelists_extpar_check_cosmo')

    OPEN(NEWUNIT=nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=extpar_consistency_check_io, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist extpar_consistency_check_io',__FILE__, __LINE__) 
    ENDIF

    CLOSE(nuin)
 
    IF (lwrite_grib) THEN
      CALL logging%warning('Direct Grib output is not supported anymore, but has been moved to an post-processing step!')
      lwrite_grib=.FALSE.
    END IF

    WRITE(message_text,'(a,i0)') 'Number of special points: ', number_special_points
    CALL logging%info(message_text)
    WRITE(message_text,'(a,i0)') 'Tile mode: ', tile_mode
    CALL logging%info(message_text)
    WRITE(message_text,*) 'Flake corrrection: ', lflake_correction
    CALL logging%info(message_text)

    IF (l_use_array_cache) THEN
      CALL logging%info('array caching (less memory consumption) is activated!')
      
      ! list of compilers not working with array-caching
      IF ( (index(INFO_CompilerVersion,'intel') /= 0) .OR. &
           (index(INFO_CompilerVersion,'INTEL') /= 0) .OR. &
           (index(INFO_CompilerVersion,'Intel') /= 0) ) THEN

        WRITE(message_text,*) 'Array-caching not supported for', &
             &                TRIM(INFO_CompilerVersion),' compiler!'
        CALL logging%error(message_text,__FILE__,__LINE__)
      ENDIF
    ENDIF

    ! Currently we only support GCC for Extpar
    IF ( (index(INFO_CompilerVersion,'gcc') == 0) .AND. &
         (index(INFO_CompilerVersion,'GCC') == 0) .AND. &
         (index(INFO_CompilerVersion,'Gcc') == 0) ) THEN

      WRITE(message_text,*) 'EXTPAR IS NOT TESTED WITH ', &
           &                TRIM(INFO_CompilerVersion),' COMPILER!'
      CALL logging%error(message_text,__FILE__,__LINE__)
    ENDIF

    CALL logging%info('Exit routine: read_namelists_extpar_check_cosmo')

  END SUBROUTINE read_namelists_extpar_check_cosmo
                                                                              
  SUBROUTINE read_namelists_extpar_special_points(namelist_file,        &
       &                                          lon_geo_sp,           &
       &                                          lat_geo_sp,           &
       &                                          soiltype_sp,          &
       &                                          z0_sp,                &
       &                                          rootdp_sp,            &
       &                                          plcovmn_sp,           &
       &                                          plcovmx_sp,           &
       &                                          laimn_sp,             &
       &                                          laimx_sp,             &
       &                                          for_d_sp,             &
       &                                          for_e_sp,             &
       &                                          fr_land_sp            )

    CHARACTER (len=*), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

    REAL(KIND=wp),    INTENT(OUT) :: lon_geo_sp,           &
         &                           lat_geo_sp,           &
         &                           soiltype_sp,          &
         &                           z0_sp,                &
         &                           rootdp_sp,            &
         &                           plcovmn_sp,           &
         &                           plcovmx_sp,           &
         &                           laimn_sp,             &
         &                           laimx_sp,             &
         &                           for_d_sp,             &
         &                           for_e_sp,             &
         &                           fr_land_sp

    !> local variables
    INTEGER(KIND=i4)              :: nuin, ierr     !< unit number

    !> define the namelist group
    NAMELIST /special_points/ &
         lon_geo_sp, lat_geo_sp, soiltype_sp, z0_sp, rootdp_sp, plcovmn_sp, plcovmx_sp, &
         laimn_sp, laimx_sp,for_d_sp,for_e_sp,fr_land_sp

    CALL logging%info('Enter routine: read_namelists_extpar_special_points') 

    !> initialization
    ierr     = 0
    lon_geo_sp  = -999.0_wp
    lat_geo_sp  = -999.0_wp
    soiltype_sp = -999.0_wp
    z0_sp       = -999.0_wp
    rootdp_sp   = -999.0_wp
    plcovmn_sp  = -999.0_wp
    plcovmx_sp  = -999.0_wp
    laimn_sp    = -999.0_wp
    for_d_sp    = -999.0_wp
    for_e_sp    = -999.0_wp
    fr_land_sp  = -999.0_wp
    

    OPEN(NEWUNIT=nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text)
    ENDIF

    READ(nuin, NML=special_points, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%warning('Please compare number of special points and defined files!')
      CALL logging%error('Cannot read in namelist extpar_special_point',__FILE__, __LINE__) 
    ENDIF

    CLOSE(nuin)

    CALL logging%info('Exit routine: read_namelists_extpar_special_points') 

  END SUBROUTINE read_namelists_extpar_special_points

END MODULE mo_read_extpar_namelists
