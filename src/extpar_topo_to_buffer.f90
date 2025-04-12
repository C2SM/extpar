!+ Fortran main program to read in GLOBE orography data and aggregate to target grid
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_1         2011/01/20 Hermann Asensio
!  small bug fixes according to Fortran compiler warnings
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
! V2_11        2024-12-10 Christian R. Steger
!  introduction of the COPERNICUS topography raw data set for external parameters
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran main program to read in GLOBE orography data and aggregate to target grid
!>
!! @par extpar_topo_to_buffer
!!
!! This program reads in the GLOBE/ASTER/MERIT/COPERNICUS orography data set
!! and aggregates the orographic height to the target grid
!! and computes the subgrid-scale orography parameters (SSO) required by the SSO-parameterization.
!!
!> Purpose: read in GLOBE/ASTER orography data and aggregate to COSMO grid
!> \author Hermann Asensio
PROGRAM extpar_topo_to_buffer

  USE, INTRINSIC :: iso_c_binding !, ONLY: c_loc, c_f_pointer ########## adjust later (todo)

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
                                
  USE mo_icon_grid_data,        ONLY: ICON_grid, &
       &                              icon_grid_region
                                
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
       &                              skyview_topo
                                
  USE mo_topo_data,             ONLY:  topo_aster,        &
       &                               topo_merit,        &
       &                               topo_copernicus,   &
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
       &                               merit_lat_min,     &
       &                               merit_lat_max,     &
       &                               merit_lon_min,     &
       &                               merit_lon_max,     &
       &                               copernicus_lat_min,  &
       &                               copernicus_lat_max,  &
       &                               copernicus_lon_min,  &
       &                               copernicus_lon_max,  &
       &                               num_tiles,         &
       &                               allocate_topo_data,&
       &                               fill_topo_data,    &
       &                               lradtopo,          &
       &                               nhori,             &
       &                               radius,            &
       &                               min_circ_cov,      &
       &                               max_missing,       &
       &                               itype_scaling,     &
       &                               deallocate_topo_fields
                                
                                
  USE mo_agg_topo_icon,         ONLY: agg_topo_data_to_target_grid_icon
  USE mo_agg_topo_cosmo,        ONLY: agg_topo_data_to_target_grid_cosmo
                                
  USE mo_topo_output_nc,        ONLY: write_netcdf_buffer_topo,    &
       &                              write_netcdf_icon_grid_topo, &
       &                              write_netcdf_cosmo_grid_topo
                                
  USE mo_oro_filter,            ONLY: read_namelists_extpar_orosmooth
  USE mo_lradtopo,              ONLY: read_namelists_extpar_lradtopo, &
       &                              compute_lradtopo, &
       &                              lradtopo_icon

  USE mo_preproc_for_sgsl,      ONLY: preproc_orography

  USE mo_agg_sgsl,              ONLY: agg_sgsl_data_to_target_grid

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

  INTEGER(c_int)                 :: num_cell_c, num_vertex_c, num_hori_c
  REAL(c_double), ALLOCATABLE    :: clon_c(:), &
       &                            clat_c(:), &
       &                            hsurf_c(:), &
       &                            vlon_c(:), &
       &                            vlat_c(:), &
       &                            horizon_topo_c(:, :), &
       &                            skyview_topo_c(:)
  INTEGER(c_int), ALLOCATABLE    :: cells_of_vertex_c(:, :)

  REAL (KIND=wp)                  :: eps_filter, &
       &                             rfill_valley,    &
       &                             rxso_mask

  LOGICAL                         :: lsso_param, &
       &                             lcompute_sgsl=.FALSE., & !compute subgrid slope
       &                             lpreproc_oro = .FALSE., & !TRUE: preproc raw oro data FALSE: read directly from NetCDF
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
       &                               lpreproc_oro,              &
       &                               lsso_param,                &
       &                               lsubtract_mean_slope,      &
       &                               orography_buffer_file,     &
       &                               orography_output_file,     &
       &                               sgsl_output_file)

  IF (lcompute_sgsl) THEN

    IF (itopo_type == 2) THEN
      CALL logging%error('ASTER topography currently not supported for SGSL', &
             __FILE__, __LINE__)
    ENDIF

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

  CALL read_namelists_extpar_lradtopo(namelist_lrad,lradtopo,nhori, radius,min_circ_cov, max_missing, itype_scaling)

  ! get information on target grid
  CALL init_target_grid(namelist_grid_def,lrad=lradtopo)

  igrid_type = tg%igrid_type

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
    IF (lcompute_sgsl) THEN
      CALL logging%error('SGSL not supported for ICON',__FILE__,__LINE__)
    ENDIF
  ENDIF

  CALL read_namelists_extpar_orosmooth(namelist_oro_smooth,  &
       &                               igrid_type,           &
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

  IF (lradtopo .AND. (.NOT. lfilter_oro) .AND. igrid_type == igrid_cosmo) THEN
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
    CASE(topo_merit)
      WRITE(message_text,*)'Edges of domain: ', merit_lon_min,' ', merit_lon_max,' ', merit_lat_min,' ',merit_lat_max
      CALL logging%info(message_text)
      IF (lon_geo (tg%ie,tg%je,tg%ke) > merit_lon_max .OR. lon_geo(1,1,1) < merit_lon_min) THEN
        WRITE(message_text,*) 'MERIT min lon is: ', merit_lon_min, ' and MERIT max lon is: ', merit_lon_max
        CALL logging%warning(message_text)
        CALL logging%error('The chosen longitude edges are not within the MERIT domain.',__FILE__,__LINE__)
      END IF
      IF (lat_geo(tg%ie,tg%je,tg%ke) > merit_lat_max .OR. lat_geo(1,1,1) < merit_lat_min) THEN
        WRITE(message_text,*) 'MERIT min lat is: ', merit_lat_min, ' and MERIT max lat is: ', merit_lat_max
        CALL logging%warning(message_text)
        CALL logging%error('The chosen latitude edges are not within the MERIT domain.',__FILE__,__LINE__)
      END IF
    CASE(topo_copernicus)
      WRITE(message_text,*)'Edges of domain: ', copernicus_lon_min,' ', copernicus_lon_max,' ', copernicus_lat_min,' ' &
           & ,copernicus_lat_max
      CALL logging%info(message_text)
      IF (lon_geo (tg%ie,tg%je,tg%ke) > copernicus_lon_max .OR. lon_geo(1,1,1) < copernicus_lon_min) THEN
        WRITE(message_text,*) 'COPERNICUS min lon is: ', copernicus_lon_min, ' and COPERNICUS max lon is: ', copernicus_lon_max
        CALL logging%warning(message_text)
        CALL logging%error('The chosen longitude edges are not within the COPERNICUS domain.',__FILE__,__LINE__)
      END IF
      IF (lat_geo(tg%ie,tg%je,tg%ke) > copernicus_lat_max .OR. lat_geo(1,1,1) < copernicus_lat_min) THEN
        WRITE(message_text,*) 'COPERNICUS min lat is: ', copernicus_lat_min, ' and COPERNICUS max lat is: ', copernicus_lat_max
        CALL logging%warning(message_text)
        CALL logging%error('The chosen latitude edges are not within the COPERNICUS domain.',__FILE__,__LINE__)
      END IF
  END SELECT

  CALL det_topo_tiles_grid(topo_tiles_grid)

  CALL logging%info('Topo input files:')
  DO k = 1, ntiles
    WRITE(message_text,'(3x,a,a,4f7.1,2i6,2x,2f8.6)')  &
         &      TRIM(topo_files(k)),              &
         &      ' Tile'//char(64+k),              &
         &      topo_tiles_grid(k)%start_lat_reg, &
         &      topo_tiles_grid(k)%end_lat_reg,   &
         &      topo_tiles_grid(k)%start_lon_reg, &
         &      topo_tiles_grid(k)%end_lon_reg,   &
         &      topo_tiles_grid(k)%nlon_reg,      &
         &      topo_tiles_grid(k)%nlat_reg,      &
         &      topo_tiles_grid(k)%dlon_reg,      &
         &      topo_tiles_grid(k)%dlat_reg       
    CALL logging%info(message_text)
  ENDDO

  CALL det_topo_grid(topo_grid)

  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------

  CALL logging%info( '')
  CALL logging%info( '============= allocate fields ==================')
  CALL logging%info( '')

  CALL logging%info('l_use_array_cache=.FALSE. -> can only be used in consistency_check')

  CALL allocate_topo_target_fields(tg,nhori,lcompute_sgsl, l_use_array_cache=.FALSE.)

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
    IF ( igrid_type == igrid_cosmo ) THEN
      CALL compute_lradtopo(nhori,tg,hh_topo,slope_asp_topo,slope_ang_topo, &
           &                horizon_topo,skyview_topo)
    ELSEIF ( igrid_type == igrid_icon ) THEN
      ! CALL lradtopo_icon(nhori, radius, min_circ_cov,tg, hh_topo, horizon_topo, &
      !      &             skyview_topo, max_missing, itype_scaling)

      ! temporary -------------------------------------------------------------
      CALL logging%info("icon_grid_region%cells%center(:)%lon")
      WRITE(message_text,*) 'Is contiguous: ', IS_CONTIGUOUS(icon_grid_region%cells%center(:)%lon)
      CALL logging%info(message_text)
      CALL logging%info("icon_grid_region%verts%vertex(:)%lat")
      WRITE(message_text,*) 'Is contiguous: ', IS_CONTIGUOUS(icon_grid_region%verts%vertex(:)%lat)
      CALL logging%info(message_text)
      CALL logging%info("icon_grid_region%verts%cell_index")
      WRITE(message_text,*) 'Is contiguous: ', IS_CONTIGUOUS(icon_grid_region%verts%cell_index)
      CALL logging%info(message_text)
      ! temporary -------------------------------------------------------------

      ! Cast arrays that are non-contiguous in memory to C types
      ALLOCATE(clon_c(icon_grid_region%ncells))
      ALLOCATE(clat_c(icon_grid_region%ncells))
      DO k = 1, icon_grid_region%ncells
        clon_c(k) = REAL(icon_grid_region%cells%center(k)%lon, KIND=c_double)
        clat_c(k) = REAL(icon_grid_region%cells%center(k)%lat, KIND=c_double)
      END DO
      ! temporary -------------------------------------------------------------
      WRITE(message_text,*) 'clat_c(7): ', clat_c(7)
      CALL logging%info(message_text)
      ! temporary -------------------------------------------------------------
      ALLOCATE(vlon_c(icon_grid_region%nverts))
      ALLOCATE(vlat_c(icon_grid_region%nverts))
      DO k = 1, icon_grid_region%nverts
        vlon_c(k) = REAL(icon_grid_region%verts%vertex(k)%lon, KIND=c_double)
        vlat_c(k) = REAL(icon_grid_region%verts%vertex(k)%lat, KIND=c_double)
      END DO
      ! temporary -------------------------------------------------------------
      WRITE(message_text,*) 'vlon_c(33): ', vlon_c(33)
      CALL logging%info(message_text)
      ! temporary -------------------------------------------------------------

      ! Cast remaining data to C types
      ALLOCATE(hsurf_c(icon_grid_region%ncells))
      hsurf_c = REAL(hh_topo(:,1,1), KIND=c_double) ! cast correct (:,1,1) -> (:)?
      ALLOCATE(cells_of_vertex_c(icon_grid_region%nverts, 6))
      cells_of_vertex_c = INT(icon_grid_region%verts%cell_index, KIND=c_int) ! shape: (575572, 6), size: 3453432
      num_cell_c = INT(icon_grid_region%ncells, KIND=c_int)
      num_vertex_c = INT(icon_grid_region%nverts, KIND=c_int)
      num_hori_c = INT(nhori, KIND=c_int)
      ! temporary -------------------------------------------------------------
      WRITE(message_text,*) 'num_cell: ', num_cell_c
      CALL logging%info(message_text)
      WRITE(message_text,*) 'num_vertex: ', num_vertex_c
      CALL logging%info(message_text)
      ! temporary -------------------------------------------------------------

      ! Allocate output arrays
      ALLOCATE(horizon_topo_c(icon_grid_region%ncells, nhori)) ! order of dim correct?
      ALLOCATE(skyview_topo_c(icon_grid_region%ncells)) ! order of dim correct?
      horizon_topo_c = 2.3 ! temporary
      skyview_topo_c = 4.7 ! temporary

      ! #######################################################################
      ! -> call C++/Embree implementation of horizon computation

      ! -----------------------------------------------------------------------
      ! The following arrays needs to be passsed:
      ! input:
      ! - lon_cell_centre (==clon) -> g%cells%center(:)%lon [radian]
      ! - lat_cell_centre (==clat) -> g%cells%center(:)%lat [radian]
      ! - longitude_vertices (== vlon) -> g%verts%vertex(:)%lon [radian]
      ! - latitude_vertices (== vlat) -> g%verts%vertex(:)%lat [radian]
      ! - cells_of_vertex(6, num_vertex) -> g%verts%cell_index (cell idx, 1 to noOfNeigbors)
           ! start with 1!, -1: "empty", no connected cell
      ! output:
      ! - horizon_topo (input/output) (num_cell,1,1,num_azim) [degree]
      ! - skyview_topo (input/output) (num_cell,1,1) [-]
      ! -----------------------------------------------------------------------
      ! Definition of interface (at top or here?):

      ! INTERFACE
      !   SUBROUTINE horayzon(clon_c, clat_c, hsurf_c, vlon_c, vlat_c, &
      !     & cells_of_vertex_c, horizon_topo_c, skyview_topo_c, &
      !     & num_cell_c, num_vertex_c, num_hori_c) bind(C, name="horayzon")
      !     USE iso_c_binding
      !     IMPLICIT NONE
      !     REAL(c_double), DIMENSION(*), INTENT(IN) :: clon_c, clat_c
      !     REAL(c_double), DIMENSION(*), INTENT(IN) :: hsurf_c
      !     REAL(c_double), DIMENSION(*), INTENT(IN) :: vlon_c, vlat_c
      !     REAL(c_int), DIMENSION(*), INTENT(IN) :: cells_of_vertex_c
      !     REAL(c_double), DIMENSION(*), INTENT(INOUT) :: horizon_topo_c
      !     REAL(c_double), DIMENSION(*), INTENT(INOUT) :: skyview_topo_c
      !     INTEGER(c_int), value :: num_cell_c, num_vertex_c, num_hori_c
      !   END SUBROUTINE horayzon
      ! END INTERFACE

      ! -----------------------------------------------------------------------

      ! CALL horayzon(clon_c, clat_c, hsurf_c, vlon_c, vlat_c, &
      !     & cells_of_vertex_c, horizon_topo_c, skyview_topo_c, &
      !     & num_cell_c, num_vertex_c, num_hori_c)

      ! #######################################################################

      ! Cast output to Fortran types
      horizon_topo(:,1,1,:) = REAL(horizon_topo_c) ! cast correct (:,:) -> (:,1,1,:)?
      skyview_topo(:,1,1) = REAL(skyview_topo_c) ! cast correct (:) -> (:,1,1)?

      ! temporary -------------------------------------------------------------
      WRITE(message_text,*) 'horizon_topo(32,1,1,5): ', horizon_topo(32,1,1,5)
      CALL logging%info(message_text)
      WRITE(message_text,*) 'skyview_topo(11,1,1): ', skyview_topo(11,1,1)
      CALL logging%info(message_text)
      ! temporary -------------------------------------------------------------

      DEALLOCATE(clon_c)
      DEALLOCATE(clat_c)
      DEALLOCATE(vlon_c)
      DEALLOCATE(vlat_c)
      DEALLOCATE(hsurf_c)
      DEALLOCATE(cells_of_vertex_c)
      DEALLOCATE(horizon_topo_c)
      DEALLOCATE(skyview_topo_c)

      ! -----------------------------------------------------------------------
    ENDIF
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
         &                           lradtopo,                &
         &                           nhori,                   &
         &                           hh_topo_max,             &
         &                           hh_topo_min,             &
         &                           horizon_topo,            &
         &                           skyview_topo,            &
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
