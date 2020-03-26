!+ Fortran main program to read in GLOBE orography data and aggregate to target grid
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_1         2011/01/20 Hermann Asensio
!  small bug fixes accroding to Fortran compiler warnings
! V1_2         2011/03/25 Hermann Asensio
!  update to support ICON refinement grids
! V1_4         2011/04/21 Anne Roches
!  implementation of orography smoothing
! V1_7         2013/01/25 Guenther Zaengl
!  Parallel threads for ICON and COSMO using Open-MP,
!  Several bug fixes and optimizations for ICON search algorithm,
!  particularly for the special case of non-contiguous domains;
!  simplified namelist control for ICON
! V2_0         2013/06/04 Martina Messmer, Anne Roches
!  introduction of the ASTER topography raw data set for external parameters
!  switch to choose if SSO parameters are desired or not
! V2_0         2013/06/04 Anne Roches
!  Implementation of the topographical corrected radiation parameters
! V1_14        2014-07-18 Juergen Helmert
!  Combined COSMO Release
! V2_1         2015-01-12 Juergen Helmert
!  Bugfix correction covers CSCS SVN r5907-r6359
! V2_6         2016-10-07 Juergen Helmert
!  Add namelist switch lfilter_topo
! V2_10        2018-02-19 Juergen Helmert
!  lsubtract_mean_slope, ERA-I surface temp for land points
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran main program to read in GLOBE orography data and aggregate to target grid
!>
!! @par extpar_topo_to_buffer
!!
!! This program reads in the GLOBE/ASTER orography data set and aggregates the orographic height to the target grid
!! and computes the subgrid-scale orography parameters (SSO) required by the SSO-parameterization.
!!
!> Purpose: read in GLOBE/ASTER orography data and aggregate to COSMO grid
!> \author Hermann Asensio
PROGRAM extpar_topo_to_buffer

  USE mo_logging
  USE info_extpar,              ONLY: info_print
  USE mo_kind,                  ONLY: wp, i4
                                
  USE mo_target_grid_data,      ONLY: lon_geo,           &
       &                              lat_geo,           &
       &                              no_raw_data_pixel, &
       &                              tg  !< structure with target grid description
                                
  USE mo_target_grid_routines,  ONLY: init_target_grid
                                
  USE mo_grid_structures,       ONLY: igrid_icon, igrid_cosmo
                                
  USE mo_cosmo_grid,            ONLY: COSMO_grid
                                
  USE mo_icon_grid_data,        ONLY: ICON_grid, icon_grid_region
                                
  USE mo_io_units,              ONLY: filename_max
                                
  USE mo_topo_routines,         ONLY: read_namelists_extpar_orography, &
       &                              det_topo_tiles_grid, &
       &                              det_topo_grid, &
       &                              read_namelists_extpar_scale_sep
                                
  USE mo_topo_tg_fields,        ONLY: fr_land_topo,                &
       &                              hh_topo,                     &
       &                              hh_topo_max,                 &
       &                              hh_topo_min,                 &
       &                              stdh_topo,                   &
       &                              theta_topo,                  &
       &                              aniso_topo,                  &
       &                              slope_topo,                  &
       &                              z0_topo,                     &
       &                              sgsl,                        &
       &                              allocate_topo_target_fields, &
       &                              slope_asp_topo,              &
       &                              slope_ang_topo,              &
       &                              horizon_topo,                &
       &                              skyview_topo,                &
       &                              vertex_param,                &
       &                              allocate_additional_param
                                
  USE mo_topo_data,             ONLY:  topo_aster,        &
       &                               itopo_type,        &
       &                               topo_tiles_grid,   &
       &                               topo_grid,         &
       &                               ntiles,            &
       &                               max_tiles,         &
       &                               nc_tot,            &
       &                               nr_tot,            &
       &                               nc_tile,           &
       &                               tiles_lon_min,     &
       &                               tiles_lon_max,     &
       &                               tiles_lat_min,     &
       &                               tiles_lat_max,     &
       &                               aster_lat_min,     &
       &                               aster_lat_max,     &
       &                               aster_lon_min,     &
       &                               aster_lon_max,     &
       &                               num_tiles,         &
       &                               allocate_topo_data,&
       &                               fill_topo_data,    &
       &                               lradtopo,          &
       &                               nhori,             &
       &                               deallocate_topo_fields
                                
                                
  USE mo_agg_topo_icon,         ONLY: agg_topo_data_to_target_grid_icon
  USE mo_agg_topo_cosmo,        ONLY: agg_topo_data_to_target_grid_cosmo
                                
  USE mo_topo_output_nc,        ONLY: write_netcdf_buffer_topo,    &
       &                              write_netcdf_icon_grid_topo, &
       &                              write_netcdf_cosmo_grid_topo
                                
  USE mo_oro_filter,            ONLY: read_namelists_extpar_orosmooth
  USE mo_lradtopo,              ONLY: read_namelists_extpar_lradtopo, &
       &                              compute_lradtopo

  USE mo_preproc_for_sgsl,      ONLY: preproc_orography

  USE mo_agg_sgsl,              ONLY: agg_sgsl_data_to_target_grid
  USE mo_preproc_topo,          ONLY: reduce_grid

  IMPLICIT NONE

  CHARACTER (len=filename_max)   :: netcdf_filename, &
       &                            namelist_grid_def, &
       &                            namelist_topo_data_input, &     !< file with input namelist with GLOBE data information
       &                            namelist_scale_sep_data_input, &!< file with input namelist with scale separated data information
       &                            namelist_oro_smooth, &          !< file with orography smoothing information (switches)
       &                            namelist_lrad, &                !< file with opo information (switches)
       &                            topo_files(1:max_tiles), &      !< filenames globe raw data
       &                            sgsl_files(1:max_tiles), &      !< filenames subgrid-slope
       &                            orography_buffer_file, &        !< name for orography buffer file
       &                            orography_output_file, &        !< name for orography output file
       &                            sgsl_output_file,      &        !< name for sgsl output file
       &                            raw_data_orography_path, &      !< path to raw data
       &                            raw_data_scale_sep_orography_path, & !< path to raw data
       &                            scale_sep_files(1:max_tiles) !< filenames globe raw data
                               
  REAL(KIND=wp)                  :: undefined                     !< value to indicate undefined grid elements
                               
  INTEGER (KIND=i4)              :: &
       &                            k,ie,je,ke, &
       &                            igrid_type, &           !< target grid type, 1 for ICON, 2 for COSMO, 3 for GME grid
       &                            nvertex, &  !< total number of vertices
       &                            ntiles_column, &        !< number of tile columns in total domain
       &                            ntiles_row, &           !< number of tile rows in total domain
       &                            ilow_pass_oro,   &
       &                            numfilt_oro,     &
       &                            ifill_valley,    &
       &                            ilow_pass_xso,   &
       &                            numfilt_xso

  INTEGER (KIND=i4), ALLOCATABLE :: topo_startrow(:), &     !< startrow indeces for each GLOBE tile
       &                            topo_endrow(:), &       !< endrow indeces for each GLOBE tile
       &                            topo_startcolumn(:), &  !< starcolumn indeces for each GLOBE tile
       &                            topo_endcolumn(:)    !< endcolumn indeces for each GLOBE tile

  REAL (KIND=wp)                  :: eps_filter, &
       &                             rfill_valley,    &
       &                             rxso_mask

  LOGICAL                         :: lsso_param, &
       &                             lcompute_sgsl=.FALSE., & !compute subgrid slope
       &                             lpreproc_oro = .FALSE., & !TRUE: preproc raw oro data FALSE: read directly from NetCDF
       &                             lreduce_grid = .FALSE., & 
       &                             lscale_separation=.FALSE., &
       &                             lscale_file= .FALSE., &
       &                             lsubtract_mean_slope, &
       &                             lfilter_oro,     &
       &                             lxso_first

  namelist_grid_def                = 'INPUT_grid_org'
  namelist_scale_sep_data_input    = 'INPUT_SCALE_SEP'
  namelist_lrad                    = 'INPUT_RADTOPO'
  namelist_topo_data_input         = 'INPUT_ORO'
  namelist_oro_smooth              = 'INPUT_OROSMOOTH'
  
  CALL initialize_logging("extpar_topo_to_buffer.log")
  CALL info_print ()


  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= start topo_to_buffer =============')
  CALL logging%info( '')

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  CALL logging%info( '')
  CALL logging%info( '============= read namelist and init grid ======')
  CALL logging%info( '')

  CALL read_namelists_extpar_orography(namelist_topo_data_input,  &
       &                               raw_data_orography_path,   &
       &                               topo_files,                &
       &                               sgsl_files,                &
       &                               ntiles_column,             &
       &                               ntiles_row,                &
       &                               itopo_type,                &
       &                               lcompute_sgsl,             &
       &                               lreduce_grid,              &
       &                               lpreproc_oro,              &
       &                               lsso_param,                &
       &                               lsubtract_mean_slope,      &
       &                               orography_buffer_file,     &
       &                               orography_output_file,     &
       &                               sgsl_output_file)

  IF (lcompute_sgsl) THEN 
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    CALL logging%info( '')
    CALL logging%warning( 'Subgrid-slope (SGSL) active')
    CALL logging%info( '')

    IF (lpreproc_oro) THEN !preprocess raw oro data to get s_oro field
      !--------------------------------------------------------------------------
      !--------------------------------------------------------------------------
      CALL logging%info( '')
      CALL logging%info( '======= SGSL: preprocess raw oro data ==========')
      CALL logging%info( '')

      CALL preproc_orography(raw_data_orography_path, &
           &                 topo_files, &
           &                 sgsl_files, &
           &                 itopo_type, &
           &                 ntiles_row, &
           &                 ntiles_column)

      !--------------------------------------------------------------------------
      !--------------------------------------------------------------------------
      CALL logging%info( '')
      CALL logging%info( '======= SGSL: end preprocess raw oro data =======')
      CALL logging%info( '')

    ELSE ! read s_oro field from files defined in namelist sgsl_io_extpar
      !--------------------------------------------------------------------------
      !--------------------------------------------------------------------------
      CALL logging%info( '')
      CALL logging%info( '======= SGSL: read S_ORO from netcdf-files ======')
      CALL logging%info( '')
    ENDIF

  ENDIF

  INQUIRE(file=TRIM(namelist_scale_sep_data_input),exist=lscale_file)
  IF (lscale_file) THEN
    CALL read_namelists_extpar_scale_sep(namelist_scale_sep_data_input,        &
         &                                  raw_data_scale_sep_orography_path, &
         &                                  scale_sep_files,                   &
         &                                  lscale_separation)
  ENDIF

  CALL read_namelists_extpar_lradtopo(namelist_lrad,lradtopo,nhori)

  ! get information on target grid
  CALL init_target_grid(namelist_grid_def,lrad=lradtopo)

  igrid_type = tg%igrid_type

  ! Checks
  IF (igrid_type == igrid_cosmo) THEN
    WRITE(message_text,*) 'lradtopo = ', lradtopo,' lfilter_oro = ', lfilter_oro
    CALL logging%info(message_text)
  ELSE
    lradtopo    = .FALSE.
    lfilter_oro = .FALSE.
  END IF

  IF (lscale_separation .AND. itopo_type == 2) THEN
    lscale_separation = .FALSE.
    CALL logging%warning('Scale separation can only be used with GLOBE as raw topography')
  ENDIF

  IF (igrid_type == igrid_cosmo) THEN ! cosmo grid
    IF (itopo_type == 1 .AND. &
    &  (cosmo_grid%dlon_rot <= 0.01 .OR. cosmo_grid%dlat_rot <= 0.01 )) &
    &  THEN
      CALL logging%warning('GLOBE raw topography data is used for horizontal grid &
           & resolution smaller than 1km')
    ENDIF
  ELSE !icon grid
    IF (itopo_type == 1 ) THEN
      CALL logging%warning('GLOBE raw topography should not be used for &
           & horizontal grid resolution smaller than 1km -> please check &
           & your icon grid to fulfill this condition.')
    ENDIF
  ENDIF

  CALL read_namelists_extpar_orosmooth(namelist_oro_smooth,  &
       &                               lfilter_oro,          &
       &                               ilow_pass_oro,        &
       &                               numfilt_oro,          &
       &                               eps_filter,           &
       &                               ifill_valley,         &
       &                               rfill_valley,         &
       &                               ilow_pass_xso,        &
       &                               numfilt_xso,          &
       &                               lxso_first,           &
       &                               rxso_mask)

  IF (lradtopo .AND. (.NOT. lfilter_oro)) THEN
    CALL logging%warning('lradtopo should not be used without orography filtering!')
  ENDIF

  ! gives back the number of tiles that are available 16 for GLOBE or 36 for ASTER
  CALL num_tiles(ntiles_column, ntiles_row, ntiles)

  ! need to be allocated after ntiles is known!
  ALLOCATE (topo_startrow(1:ntiles), topo_endrow(1:ntiles),topo_startcolumn(1:ntiles),topo_endcolumn(1:ntiles))

  CALL allocate_topo_data(ntiles)      ! allocates the data using ntiles

  CALL fill_topo_data(raw_data_orography_path,topo_files, &! the allocated vectors need to be filled with the respective value.
       &              tiles_lon_min, &
       &              tiles_lon_max, &
       &              tiles_lat_min, &
       &              tiles_lat_max, &
       &              nc_tot,        &
       &              nr_tot,        &
       &              nc_tile)

  SELECT CASE(itopo_type)
    CASE(topo_aster)
      WRITE(message_text,*)'Edges of domain: ', aster_lon_min,' ', aster_lon_max,' ', aster_lat_min,' ',aster_lat_max
      CALL logging%info(message_text)
      IF (lon_geo (tg%ie,tg%je,tg%ke) > aster_lon_max .OR. lon_geo(1,1,1) < aster_lon_min) THEN
        WRITE(message_text,*) 'ASTER min lon is: ', aster_lon_min, ' and ASTER max lon is: ', aster_lon_max
        CALL logging%warning(message_text)
        CALL logging%error('The chosen longitude edges are not within the ASTER domain.',__FILE__,__LINE__)
      END IF
      IF (lat_geo(tg%ie,tg%je,tg%ke) > aster_lat_max .OR. lat_geo(1,1,1) < aster_lat_min) THEN
        WRITE(message_text,*) 'ASTER min lat is: ', aster_lat_min, ' and ASTER max lat is: ', aster_lat_max
        CALL logging%warning(message_text)
        CALL logging%error('The chosen latitude edges are not within the ASTER domain.',__FILE__,__LINE__)
      END IF
  END SELECT

  CALL det_topo_tiles_grid(topo_tiles_grid)

  CALL logging%info('Topo input files:')
  DO k = 1, ntiles
    WRITE(message_text,'(3x,a,a,4f7.1,2i6)')  &
         &      TRIM(topo_files(k)),              &
         &      ' Tile'//char(64+k),              &
         &      topo_tiles_grid(k)%start_lat_reg, &
         &      topo_tiles_grid(k)%end_lat_reg,   &
         &      topo_tiles_grid(k)%start_lon_reg, &
         &      topo_tiles_grid(k)%end_lon_reg,   &
         &      topo_tiles_grid(k)%nlon_reg,      &
         &      topo_tiles_grid(k)%nlat_reg
    CALL logging%info(message_text)
  ENDDO

  CALL det_topo_grid(topo_grid)

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info( '============= allocate fields ==================')
  CALL logging%info( '')

  CALL allocate_topo_target_fields(tg,nhori,lcompute_sgsl)

  ! allocate additional fields for icon grid
  SELECT CASE(igrid_type)
  CASE(igrid_icon) ! ICON GRID
    ! allocate addtional target fields
    nvertex = icon_grid_region%nverts
    CALL  allocate_additional_param(nvertex, lcompute_sgsl)
  END SELECT

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  IF (lreduce_grid) THEN
    CALL logging%info( '')
    CALL logging%info('===== Grid reduction: start reducing grid ======')
    CALL logging%info( '')

    CALL reduce_grid(topo_tile_grid, &
         &           topo_files, &
         &           raw_dat_orography_path)
  ENDIF


  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= start aggregation ================')
  CALL logging%info( '')

  IF (lcompute_sgsl) THEN
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    CALL logging%info( '')
    CALL logging%info( '======= SGSL: start aggregation ================')
    CALL logging%info( '')

    CALL agg_sgsl_data_to_target_grid(topo_tiles_grid, &
         &                           topo_grid,        &
         &                           tg,               &
         &                           sgsl_files,       &
         &                           sgsl,         &
         &                           no_raw_data_pixel,    &
         &                           raw_data_sgsl_path=raw_data_orography_path)

    !jj_tmp: reset field for topo aggregation
    no_raw_data_pixel= 0

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    CALL logging%info( '')
    CALL logging%info( '======= SGSL: all calculations done ============')
    CALL logging%info( '')

  ENDIF

  IF (igrid_type == igrid_icon) THEN ! ICON GRID
      CALL agg_topo_data_to_target_grid_icon(topo_tiles_grid,  &
           &                                 topo_grid,        &
           &                                 tg,               &
           &                                 topo_files,       &
           &                                 lsso_param,       &
           &                                 lsubtract_mean_slope, &
           &                                 lfilter_oro,      &
           &                                 ilow_pass_oro,    &
           &                                 numfilt_oro,      &
           &                                 eps_filter,       &
           &                                 ifill_valley,     &
           &                                 rfill_valley,     &
           &                                 ilow_pass_xso,    &
           &                                 numfilt_xso,      &
           &                                 lxso_first,       &
           &                                 rxso_mask,        &
           &                                 hh_topo,          &
           &                                 hh_topo_max,      &
           &                                 hh_topo_min,      &
           &                                 stdh_topo,        &
           &                                 fr_land_topo,     &
           &                                 z0_topo,          &
           &                                 no_raw_data_pixel,&
           &                                 theta_topo,       &
           &                                 aniso_topo,       &
           &                                 slope_topo,       &
           &                                 raw_data_orography_path)

  ELSE  ! COSMO/GME GRID
    CALL agg_topo_data_to_target_grid_cosmo(topo_tiles_grid,   &
         &                                  topo_grid,         &
         &                                  tg,                &
         &                                  topo_files,        &
         &                                  lsso_param,        &
         &                                  lscale_separation, &
         &                                  lfilter_oro,       &
         &                                  ilow_pass_oro,     &
         &                                  numfilt_oro,       &
         &                                  eps_filter,        &
         &                                  ifill_valley,      &
         &                                  rfill_valley,      &
         &                                  ilow_pass_xso,     &
         &                                  numfilt_xso,       &
         &                                  lxso_first,        &
         &                                  rxso_mask,         &
         &                                  hh_topo,           &
         &                                  stdh_topo,         &
         &                                  fr_land_topo,      &
         &                                  z0_topo,           &
         &                                  no_raw_data_pixel, &
         &                                  theta_topo,        &
         &                                  aniso_topo,        &
         &                                  slope_topo,        &
         &                                  raw_data_orography_path, &
         &                                  raw_data_scale_sep_orography_path, &
         &                                  scale_sep_files)

  END IF !igrid_type

  ! if the target domain has a higher resolution of than the GLOBE data set (30'') some grid elements might not
  ! be set by the routine agg_topo_data_to_target_grid, (no_raw_data_pixel(ie,je,ke) == 0 in this case
  ! loop overa all grid elements to check and perform a bilinear interplation if necessary
  k = 0
  undefined = -999.9_wp

  ! consistency for small grid sizes, do not use estimates of variance for small sample size
  IF (MAXVAL(no_raw_data_pixel) < 10) THEN 
    IF (lsso_param) THEN
      stdh_topo  = 0.0_wp
      theta_topo = 0.0_wp
      aniso_topo = 0.0_wp
      slope_topo = 0.0_wp
    ENDIF
    z0_topo     = 0.0_wp
  ENDIF

  DO ke=1,tg%ke
    DO je=1,tg%je
      DO ie=1,tg%ie
        IF (no_raw_data_pixel(ie,je,ke) < 1 ) THEN
          IF (lsso_param) THEN
            stdh_topo(ie,je,ke)  = 0.0_wp
            theta_topo(ie,je,ke) = 0.0_wp
            aniso_topo(ie,je,ke) = 0.0_wp
            slope_topo(ie,je,ke) = 0.0_wp
          ENDIF
          z0_topo(ie,je,ke)     = 0.0_wp
        ENDIF
      ENDDO
    ENDDO
  ENDDO

  ! compute the lradtopo parameters if needed
  IF ( lradtopo ) THEN
    CALL compute_lradtopo(nhori,tg,hh_topo,slope_asp_topo,slope_ang_topo,horizon_topo,skyview_topo)
  ENDIF

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= write data to netcdf==============')
  CALL logging%info( '')

  ! output to netcdf file
  undefined = -999.9_wp



  netcdf_filename = TRIM(orography_buffer_file)

  CALL write_netcdf_buffer_topo(netcdf_filename, &
       &                        tg,              &
       &                        undefined,       &
       &                        igrid_type,      &
       &                        lon_geo,         &
       &                        lat_geo,         &
       &                        fr_land_topo,    &
       &                        hh_topo,         &
       &                        stdh_topo,       &
       &                        z0_topo,         &
       &                        lradtopo,        &
       &                        lsso_param,      &
       &                        lcompute_sgsl,   &
       &                        nhori,           &
       &                        hh_topo_max,     &
       &                        hh_topo_min,     &
       &                        theta_topo,      &
       &                        aniso_topo,      &
       &                        slope_topo,      &
       &                        slope_asp_topo,  &
       &                        slope_ang_topo,  &
       &                        horizon_topo,    &
       &                        skyview_topo,    &
       &                        vertex_param,    &
       &                        sgsl)


  netcdf_filename = TRIM(orography_output_file)

  SELECT CASE(igrid_type)

  CASE(igrid_icon)
    CALL write_netcdf_icon_grid_topo(netcdf_filename,         &
         &                           icon_grid,               &
         &                           tg,                      &
         &                           undefined,               &
         &                           lon_geo,                 &
         &                           lat_geo,                 &
         &                           fr_land_topo,            &
         &                           hh_topo,                 &
         &                           stdh_topo,               &
         &                           z0_topo,                 &
         &                           lsso_param,              &
         &                           vertex_param,            &
         &                           hh_topo_max,             &
         &                           hh_topo_min,             &
         &                           theta_topo,              &
         &                           aniso_topo,              &
         &                           slope_topo)          

  CASE(igrid_cosmo) ! COSMO grid
    CALL write_netcdf_cosmo_grid_topo(netcdf_filename, &
         &                        cosmo_grid,      &
         &                        tg,              &
         &                        undefined,       &
         &                        lon_geo,         &
         &                        lat_geo,         &
         &                        fr_land_topo,    &
         &                        hh_topo,         &
         &                        stdh_topo,       &
         &                        z0_topo,         &
         &                        lradtopo,        &
         &                        lsso_param,      &
         &                        lcompute_sgsl,   &
         &                        nhori,           &
         &                        theta_topo,      &
         &                        aniso_topo,      &
         &                        slope_topo,      &
         &                        slope_asp_topo,  &
         &                        slope_ang_topo,  &
         &                        horizon_topo,    &
         &                        skyview_topo,    &
         &                        sgsl)


  END SELECT

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info('============= deallocate fields =================')
  CALL logging%info( '')

  CALL deallocate_topo_fields(lcompute_sgsl)

  DEALLOCATE (topo_startrow, topo_endrow, topo_startcolumn, topo_endcolumn)

  CALL logging%info( '')
  CALL logging%info('============= topo_to_buffer done ===============')

END PROGRAM extpar_topo_to_buffer
