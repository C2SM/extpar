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

  ! Load the library information data:
  USE info_extpar,             ONLY: info_define, info_print

  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind,                 ONLY: wp, i4, i8

  USE mo_target_grid_data,     ONLY: lon_geo,           &
       &                             lat_geo,           &
       &                             no_raw_data_pixel, &
       &                             tg  !< structure with target grid description

  USE mo_target_grid_routines, ONLY: init_target_grid

  USE mo_grid_structures,      ONLY: igrid_icon, igrid_cosmo

  USE mo_cosmo_grid,           ONLY: COSMO_grid

  !< structure which contains the definition of the ICON grid
  USE mo_icon_grid_data,       ONLY: ICON_grid, icon_grid_region

  USE mo_io_units,             ONLY: filename_max

  USE mo_utilities_extpar,     ONLY: abort_extpar

  USE mo_topo_routines,        ONLY: read_namelists_extpar_orography, &
       &                             read_namelists_extpar_scale_sep

  USE mo_topo_tg_fields,       ONLY: fr_land_topo,                &
       &                             hh_topo,                     &
       &                             hh_topo_max,                 &
       &                             hh_topo_min,                 &
       &                             stdh_topo,                   &
       &                             theta_topo,                  &
       &                             aniso_topo,                  &
       &                             slope_topo,                  &
       &                             z0_topo,                     &
       &                             allocate_topo_target_fields, &
       &                             slope_asp_topo,              &
       &                             slope_ang_topo,              &
       &                             horizon_topo,                &
       &                             skyview_topo,                &
       &                             vertex_param,                &
       &                             allocate_additional_hh_param

  USE mo_topo_data,            ONLY:  topo_aster,        &
       &                              itopo_type,        &
       &                              topo_tiles_grid,   &
       &                              topo_grid,         &
       &                              ntiles,            &
       &                              max_tiles,         &
       &                              nc_tot,            &
       &                              nr_tot,            &
       &                              nc_tile,           &
       &                              tiles_lon_min,     &
       &                              tiles_lon_max,     &
       &                              tiles_lat_min,     &
       &                              tiles_lat_max,     &
       &                              aster_lat_min,     &
       &                              aster_lat_max,     &
       &                              aster_lon_min,     &
       &                              aster_lon_max,     &
       &                              num_tiles,         &
       &                              allocate_topo_data,&
       &                              fill_topo_data,    &
       &                              lradtopo,          &
       &                              nhori,             &
       &                              deallocate_topo_fields

  USE mo_topo_routines,        ONLY:  det_topo_tiles_grid, &
       &                              det_topo_grid

  USE mo_agg_topo_icon,        ONLY: agg_topo_data_to_target_grid_icon
  USE mo_agg_topo_cosmo,       ONLY: agg_topo_data_to_target_grid_cosmo

  USE mo_topo_output_nc,       ONLY: write_netcdf_buffer_topo,    &
       &                             write_netcdf_icon_grid_topo, &
       &                             write_netcdf_cosmo_grid_topo

  USE mo_oro_filter,           ONLY: read_namelists_extpar_orosmooth
  USE mo_lradtopo,             ONLY: read_namelists_extpar_lradtopo, &
       &                             compute_lradtopo

  IMPLICIT NONE

  CHARACTER (len=filename_max) :: netcdf_filename
  CHARACTER (len=filename_max) :: namelist_grid_def

  CHARACTER (len=filename_max) :: namelist_topo_data_input     !< file with input namelist with GLOBE data information
  CHARACTER (len=filename_max) :: namelist_scale_sep_data_input!< file with input namelist with scale separated data information
  !roa >
  CHARACTER (len=filename_max) :: namelist_oro_smooth          !< file with orography smoothing information (switches)
  CHARACTER (len=filename_max) :: namelist_lrad                !< file with opo information (switches)
  !roa <

  CHARACTER (LEN=filename_max) :: topo_files(1:max_tiles)      !< filenames globe raw data

  CHARACTER (len=filename_max) :: orography_buffer_file        !< name for orography buffer file
  CHARACTER (len=filename_max) :: orography_output_file        !< name for orography output file

  CHARACTER (len=filename_max) :: raw_data_orography_path      !< path to raw data
  CHARACTER (len=filename_max) :: raw_data_scale_sep_orography_path !< path to raw data
  CHARACTER (LEN=filename_max) :: scale_sep_files(1:max_tiles) !< filenames globe raw data

  REAL(wp)                :: undefined                    !< value to indicate undefined grid elements
  INTEGER (i4)            :: undef_int                    !< value for undefined integer

  INTEGER (i4), ALLOCATABLE :: topo_startrow(:)     !< startrow indeces for each GLOBE tile
  INTEGER (i4), ALLOCATABLE :: topo_endrow(:)       !< endrow indeces for each GLOBE tile
  INTEGER (i4), ALLOCATABLE :: topo_startcolumn(:)  !< starcolumn indeces for each GLOBE tile
  INTEGER (i4), ALLOCATABLE :: topo_endcolumn(:)    !< endcolumn indeces for each GLOBE tile

  INTEGER :: k !< counter
  INTEGER(i8) :: ie !< counter
  INTEGER(i8) :: je !< counter
  INTEGER(i8) :: ke !< counter

  INTEGER (i4) :: igrid_type           !< target grid type, 1 for ICON, 2 for COSMO, 3 for GME grid

  ! variables for the ICON grid
  INTEGER :: nvertex  !< total number of vertices

  REAL :: timestart
  REAL :: timeend
  REAL :: timediff

  INTEGER (i4) :: ntiles_column        !< number of tile columns in total domain
  INTEGER (i4) :: ntiles_row           !< number of tile rows in total domain
  LOGICAL           :: lsso_param
  LOGICAL           :: lscale_separation
  LOGICAL           :: lsubtract_mean_slope

  LOGICAL           ::  &
       lfilter_oro,     &
       lxso_first

  INTEGER(i4)  ::       &
       ilow_pass_oro,   &
       numfilt_oro,     &
       ifill_valley,    &
       ilow_pass_xso,   &
       numfilt_xso

  REAL(wp)     ::       &
       eps_filter,      &
       rfill_valley,    &
       rxso_mask
  
  ! Print the default information to stdout:
  CALL info_define ()
  CALL info_print ()
  !--------------------------------------------------------------------------------------------------------

  namelist_grid_def = 'INPUT_grid_org'

  namelist_lrad     = 'INPUT_RADTOPO'

  CALL read_namelists_extpar_lradtopo(namelist_lrad,lradtopo,nhori)

  ! get information on target grid
  CALL init_target_grid(namelist_grid_def,lrad=lradtopo)

  print*,'Target grid tg: ',tg%ie, tg%je, tg%ke, tg%minlon, tg%maxlon, tg%minlat, tg%maxlat
  igrid_type = tg%igrid_type

  ! Checks
  IF (igrid_type /= igrid_cosmo) THEN
    lradtopo    = .FALSE.
    lfilter_oro = .FALSE.
    print*,"The model grid is not COSMO: Assume "
    print*,"   lradtopo = .FALSE. "
    print*,"   lfilter_oro = .FALSE. "
  END IF

  !--------------------------------------------------------------------------------------------------------
  ! get GLOBE raw data information
  !
  ! read namelist with globe data information

  namelist_topo_data_input = 'INPUT_ORO'
  CALL read_namelists_extpar_orography(namelist_topo_data_input,  &
       &                               raw_data_orography_path,   &
       &                               topo_files,                &
       &                               ntiles_column,             &
       &                               ntiles_row,                &
       &                               itopo_type,                &
       &                               lsso_param,                &
       &                               lsubtract_mean_slope,      &
       &                               orography_buffer_file,     &
       &                               orography_output_file)

  namelist_scale_sep_data_input = 'INPUT_SCALE_SEP'
  CALL read_namelists_extpar_scale_sep(namelist_scale_sep_data_input,        &
       &                                  raw_data_scale_sep_orography_path, &
       &                                  scale_sep_files,                   &
       &                                  lscale_separation)

  IF (lscale_separation .AND. itopo_type.eq.2) THEN   !_br 21.02.14 replaced eq by eqv
    lscale_separation = .FALSE.
    print*, '*** Scale separation can only be used with GLOBE as raw topography ***'
  ENDIF

  CALL num_tiles(itopo_type, ntiles_column, ntiles_row, ntiles)
  ! gives back the number of tiles that are available 16 for GLOBE or 36 for ASTER

  ! need to be allocated after ntiles is known!
  ALLOCATE (topo_startrow(1:ntiles), topo_endrow(1:ntiles),topo_startcolumn(1:ntiles),topo_endcolumn(1:ntiles))  !_br 21.02.14

  CALL allocate_topo_data(ntiles)                  ! allocates the data using ntiles

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
    print*,'   edges of domain: ', aster_lon_min,' ', aster_lon_max,' ', aster_lat_min,' ',aster_lat_max
    IF (lon_geo (tg%ie,tg%je,tg%ke) > aster_lon_max .OR. lon_geo(1,1,1) < aster_lon_min) THEN
      print*, '   ASTER min lon is: ', aster_lon_min, ' and ASTER max lon is: ', aster_lon_max
      CALL abort_extpar('The chosen longitude edges are not within the ASTER domain.')
    END IF
    IF (lat_geo(tg%ie,tg%je,tg%ke) > aster_lat_max .OR. lat_geo(1,1,1) < aster_lat_min) THEN
      print*, '   ASTER min lat is: ', aster_lat_min, ' and ASTER max lat is: ', aster_lat_max
      CALL abort_extpar('The chosen latitude edges are not within the ASTER domain.')
    END IF
  END SELECT

  namelist_oro_smooth = 'INPUT_OROSMOOTH'

  CALL read_namelists_extpar_orosmooth(namelist_oro_smooth,  &
&                                               lfilter_oro,          &
&                                               ilow_pass_oro,        &
&                                               numfilt_oro,          &
&                                               eps_filter,           &
&                                               ifill_valley,         &
&                                               rfill_valley,         &
&                                               ilow_pass_xso,        &
&                                               numfilt_xso,          &
&                                               lxso_first,           &
&                                               rxso_mask)

  IF (lradtopo .AND. (.NOT. lfilter_oro)) THEN
    print *,' Warning *** lradtopo should not be used without orography filtering *** '
    print *,'                            (consistency problem)                        '
  ENDIF

  CALL det_topo_tiles_grid(topo_tiles_grid)

  print*,'Topo input files:'
  DO k = 1, ntiles
    print '(3x,a,a,4f7.1,2i6)',                   &
         &      TRIM(topo_files(k)),              &
         &      ' Tile'//char(64+k),              &
         &      topo_tiles_grid(k)%start_lat_reg, &
         &      topo_tiles_grid(k)%end_lat_reg,   &
         &      topo_tiles_grid(k)%start_lon_reg, &
         &      topo_tiles_grid(k)%end_lon_reg,   &
         &      topo_tiles_grid(k)%nlon_reg,      &
         &      topo_tiles_grid(k)%nlat_reg
  ENDDO

  CALL det_topo_grid(topo_grid)
  print '(3x,a,4f7.1,2i6)','Full grid size: ', &
       &        topo_grid%start_lat_reg,       &
       &        topo_grid%end_lat_reg,         &
       &        topo_grid%start_lon_reg,       &
       &        topo_grid%end_lon_reg,         &
       &        topo_grid%nlon_reg,            &
       &        topo_grid%nlat_reg

  ! allocate globe fields for target grid
  !--------------------------------------------------------------------------------------------------------
  CALL allocate_topo_target_fields(tg,nhori)

  ! allocate additional fields for icon grid
  !--------------------------------------------------------------------------------------------------------

  SELECT CASE(igrid_type)
  CASE(igrid_icon) ! ICON GRID
    ! allocate addtional target fields
    nvertex = icon_grid_region%nverts
    CALL  allocate_additional_hh_param(nvertex)
  END SELECT

  ! call the aggregation routine
  !--------------------------------------------------------------------------------------------------------
  IF (igrid_type == igrid_icon) THEN ! ICON GRID

    IF (lsso_param) THEN
      print *,'CALL agg_topo_data_to_target_grid with SSO'
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
           &                                 raw_data_orography_path=raw_data_orography_path)
    ELSE
      PRINT *,'CALL agg_topo_data_to_target_grid without SSO'
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
           &                                 no_raw_data_pixel, &
           &                                 raw_data_orography_path=raw_data_orography_path)
    ENDIF

  ELSE  ! COSMO/GME GRID

    IF (lsso_param) THEN

      IF (lscale_separation) THEN
        PRINT *,'CALL agg_topo_data_to_target_grid_cosmo with SSO'
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
             &                                  raw_data_orography_path=raw_data_orography_path,                     & !_br 17.09.14
             &                                  raw_data_scale_sep_orography_path=raw_data_scale_sep_orography_path, & !_br 17.09.14
             &                                  scale_sep_files = scale_sep_files)
      ELSE
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
             &                                  raw_data_orography_path=raw_data_orography_path)
      ENDIF

    ELSE

      IF (lscale_separation) THEN
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
             &                                  raw_data_orography_path=raw_data_orography_path,                     & 
             &                                  raw_data_scale_sep_orography_path=raw_data_scale_sep_orography_path, & 
             &                                  scale_sep_files = scale_sep_files)
        !
      ELSE
        PRINT *,'CALL agg_topo_data_to_target_grid_cosmo without SSO'
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
             &                                  raw_data_orography_path=raw_data_orography_path) 
      ENDIF

    ENDIF

  END IF !igrid_type

  !HA debug
  !PRINT *,'agg_topo_data_to_target_grid finished'

  ! if the target domain has a higher resolution of than the GLOBE data set (30'') some grid elements might not
  ! be set by the routine agg_topo_data_to_target_grid, (no_raw_data_pixel(ie,je,ke) == 0 in this case
  ! loop overa all grid elements to check and perform a bilinear interplation if necessary
  k = 0
  undefined = -999.9_wp

  PRINT*,'lsubtract_mean_slope is set to ',lsubtract_mean_slope

  PRINT *,'Maximum number of TOPO raw data pixel in a target grid element: '
  PRINT *,'MAXVAL(no_raw_data_pixel): ', MAXVAL(no_raw_data_pixel)

  PRINT *,'Minimal number of TOPO raw data pixel in a target grid element: '
  PRINT *,'MINVAL(no_raw_data_pixel): ', MINVAL(no_raw_data_pixel)

  PRINT *,'Maximum stdh_topo: '
  PRINT *,'MAXVAL(stdh_topo): ', MAXVAL(stdh_topo)

  PRINT *,'Minimal stdh_topo: '
  PRINT *,'MINVAL(stdh_topo): ', MINVAL(stdh_topo)

  ! consistency for small grid sizes, do not use estimates of variance for small sample size
  !   IF ( (MAXVAL(no_raw_data_pixel)< 10).OR. (MINVAL(no_raw_data_pixel)==0)) THEN
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


  !   SELECT CASE(igrid_type)
  !     CASE(igrid_icon) ! ICON GRID
  !     nvertex = icon_domain_grid%nverts
  !     je=1
  !     ke=1
  !     DO nv=1, nvertex
  !       k=0
  !       IF (vertex_param%npixel_vert(nv,je,ke) == 0) THEN ! interpolate from raw data in this case
  !         k = k + 1 ! count number of grid element for which a bilinear interpolation is done
  !          point_lon_geo =  rad2deg * icon_domain_grid%verts%vertex(nv)%lon
  !          point_lat_geo =  rad2deg * icon_domain_grid%verts%vertex(nv)%lat
  !          CALL bilinear_interpol_topo_to_target_point(topo_files, &
  !            &                                 topo_grid,          &
  !            &                                 topo_files,         &
  !            &                                 undefined,          &
  !            &                                 point_lon_geo,      &
  !            &                                 point_lat_geo,      &
  !            &                                 nrows,              &
  !            &                                 h_band,             &
  !            &                                 startrow_index,     &
  !            &                                 endrow_index,       &
  !            &                                 topo_target_value)
  !
  !          vertex_param%hh_vert(nv,je,ke) = topo_target_value
  !       ENDIF
  !     ENDDO
  !     print *,'bilinear interpolation of GLOBE data used for no of grid vertices: ',k
  !   END SELECT

  ! compute the lradtopo parameters if needed
  IF ( lradtopo ) THEN
    CALL CPU_TIME(timestart)
    CALL compute_lradtopo(nhori,tg,hh_topo,slope_asp_topo,slope_ang_topo,horizon_topo,skyview_topo)
    CALL CPU_TIME(timeend)
    timediff = timeend - timestart
    PRINT *,'lradtopo calculations were done in: ', timediff
  ENDIF

  ! output to netcdf file
  undefined = -999.9_wp
  undef_int = -999

  netcdf_filename = TRIM(orography_buffer_file)
  print *, 'Output filename: ',TRIM(netcdf_filename)

  SELECT CASE(igrid_type)
  CASE(igrid_icon)
    IF (lsso_param) THEN
      CALL write_netcdf_buffer_topo(netcdf_filename,         &
           &                        tg,                      &
           &                        undefined,               &
           &                        undef_int,               &
           &                        igrid_type,              &
           &                        lon_geo,                 &
           &                        lat_geo,                 &
           &                        fr_land_topo,            &
           &                        hh_topo,                 &
           &                        stdh_topo,               &
           &                        z0_topo,                 &
           &                        lradtopo,                &
           &                        nhori,                   &
           &                        hh_topo_max=hh_topo_max, &
           &                        hh_topo_min=hh_topo_min, &
           &                        theta_topo=theta_topo,   &
           &                        aniso_topo=aniso_topo,   &
           &                        slope_topo=slope_topo,   &
           &                        vertex_param=vertex_param)
    ELSE
      CALL write_netcdf_buffer_topo(netcdf_filename,         &
           &                        tg,                      &
           &                        undefined,               &
           &                        undef_int,               &
           &                        igrid_type,              &
           &                        lon_geo,                 &
           &                        lat_geo,                 &
           &                        fr_land_topo,            &
           &                        hh_topo,                 &
           &                        stdh_topo,               &
           &                        z0_topo,                 &
           &                        lradtopo,                &
           &                        nhori,                   &
           &                        hh_topo_max=hh_topo_max, &
           &                        hh_topo_min=hh_topo_min, &
           &                        vertex_param=vertex_param)
    ENDIF

  CASE DEFAULT

    IF (lradtopo) THEN

      IF (lsso_param) THEN
        CALL write_netcdf_buffer_topo(netcdf_filename,                &
             &                        tg,                             &
             &                        undefined,                      &
             &                        undef_int,                      &
             &                        igrid_type,                     &
             &                        lon_geo,                        &
             &                        lat_geo,                        &
             &                        fr_land_topo,                   &
             &                        hh_topo,                        &
             &                        stdh_topo,                      &
             &                        z0_topo,                        &
             &                        lradtopo,                       &
             &                        nhori,                          &
             &                        theta_topo=theta_topo,          &
             &                        aniso_topo=aniso_topo,          &
             &                        slope_topo=slope_topo,          &
             &                        slope_asp_topo=slope_asp_topo,  &
             &                        slope_ang_topo= slope_ang_topo, &
             &                        horizon_topo=horizon_topo,      &
             &                        skyview_topo=skyview_topo)
      ELSE
        CALL write_netcdf_buffer_topo(netcdf_filename,                &
             &                        tg,                             &
             &                        undefined,                      &
             &                        undef_int,                      &
             &                        igrid_type,                     &
             &                        lon_geo,                        &
             &                        lat_geo,                        &
             &                        fr_land_topo,                   &
             &                        hh_topo,                        &
             &                        stdh_topo,                      &
             &                        z0_topo,                        &
             &                        lradtopo,                       &
             &                        nhori,                          &
             &                        slope_asp_topo=slope_asp_topo,  &
             &                        slope_ang_topo= slope_ang_topo, &
             &                        horizon_topo=horizon_topo,      &
             &                        skyview_topo=skyview_topo)
      ENDIF

    ELSE

      IF (lsso_param) THEN
        CALL write_netcdf_buffer_topo(netcdf_filename,       &
             &                        tg,                    &
             &                        undefined,             &
             &                        undef_int,             &
             &                        igrid_type,            &
             &                        lon_geo,               &
             &                        lat_geo,               &
             &                        fr_land_topo,          &
             &                        hh_topo,               &
             &                        stdh_topo,             &
             &                        z0_topo,               &
             &                        lradtopo,              &
             &                        nhori,                 &
             &                        theta_topo=theta_topo, &
             &                        aniso_topo=aniso_topo, &
             &                        slope_topo=slope_topo)
      ELSE
        CALL write_netcdf_buffer_topo(netcdf_filename, &
             &                        tg,              &
             &                        undefined,       &
             &                        undef_int,       &
             &                        igrid_type,      &
             &                        lon_geo,         &
             &                        lat_geo,         &
             &                        fr_land_topo,    &
             &                        hh_topo,         &
             &                        stdh_topo,       &
             &                        z0_topo,         &
             &                        lradtopo,        &
             &                        nhori)
      ENDIF

    ENDIF

  END SELECT


  SELECT CASE(igrid_type)
  CASE(igrid_icon)

    netcdf_filename = TRIM(orography_output_file)
    PRINT *,'write out ', TRIM(netcdf_filename)

    IF (lsso_param) THEN
      CALL write_netcdf_icon_grid_topo(netcdf_filename,         &
           &                           icon_grid,               &
           &                           tg,                      &
           &                           undefined,               &
           &                           undef_int,               &
           &                           lon_geo,                 &
           &                           lat_geo,                 &
           &                           fr_land_topo,            &
           &                           hh_topo,                 &
           &                           stdh_topo,               &
           &                           z0_topo,                 &
           &                           vertex_param,            &
           &                           hh_topo_max=hh_topo_max, &
           &                           hh_topo_min=hh_topo_min, &
           &                           theta_topo=theta_topo,   &
           &                           aniso_topo=aniso_topo,   &
           &                           slope_topo=slope_topo)
    ELSE
      CALL write_netcdf_icon_grid_topo(netcdf_filename, &
           &                           icon_grid,       &
           &                           tg,              &
           &                           undefined,       &
           &                           undef_int,       &
           &                           lon_geo,         &
           &                           lat_geo,         &
           &                           fr_land_topo,    &
           &                           hh_topo,         &
           &                           stdh_topo,       &
           &                           z0_topo,         &
           &                           vertex_param)
    ENDIF

  CASE(igrid_cosmo) ! COSMO grid

    netcdf_filename = TRIM(orography_output_file)
    PRINT *,'write out ', TRIM(netcdf_filename)

    IF(lradtopo) THEN

      IF (lsso_param) THEN
        CALL write_netcdf_cosmo_grid_topo(netcdf_filename,                 &
             &                            cosmo_grid,                      &
             &                            tg,                              &
             &                            undefined,                       &
             &                            undef_int,                       &
             &                            lon_geo,                         &
             &                            lat_geo,                         &
             &                            fr_land_topo,                    &
             &                            hh_topo,                         &
             &                            stdh_topo,                       &
             &                            z0_topo,                         &
             &                            lradtopo,                        &
             &                            nhori,                           &
             &                            theta_topo=theta_topo,           &
             &                            aniso_topo=aniso_topo,           &
             &                            slope_topo=slope_topo,           &
             &                            slope_asp_topo=slope_asp_topo,   &
             &                            slope_ang_topo=slope_ang_topo,   &
             &                            horizon_topo=horizon_topo,       &
             &                            skyview_topo=skyview_topo)
      ELSE
        CALL write_netcdf_cosmo_grid_topo(netcdf_filename,                 &
             &                            cosmo_grid,                      &
             &                            tg,                              &
             &                            undefined,                       &
             &                            undef_int,                       &
             &                            lon_geo,                         &
             &                            lat_geo,                         &
             &                            fr_land_topo,                    &
             &                            hh_topo,                         &
             &                            stdh_topo,                       &
             &                            z0_topo,                         &
             &                            lradtopo,                        &
             &                            nhori,                           &
             &                            slope_ang_topo=slope_ang_topo,   &
             &                            horizon_topo=horizon_topo,       &
             &                            skyview_topo=skyview_topo)
      ENDIF

    ELSE

      IF (lsso_param) THEN
        CALL write_netcdf_cosmo_grid_topo(netcdf_filename,       &
             &                            cosmo_grid,            &
             &                            tg,                    &
             &                            undefined,             &
             &                            undef_int,             &
             &                            lon_geo,               &
             &                            lat_geo,               &
             &                            fr_land_topo,          &
             &                            hh_topo,               &
             &                            stdh_topo,             &
             &                            z0_topo,               &
             &                            lradtopo,              &
             &                            nhori,                 &
             &                            theta_topo=theta_topo, &
             &                            aniso_topo=aniso_topo, &
             &                            slope_topo=slope_topo)
      ELSE
        CALL write_netcdf_cosmo_grid_topo(netcdf_filename, &
             &                            cosmo_grid,      &
             &                            tg,              &
             &                            undefined,       &
             &                            undef_int,       &
             &                            lon_geo,         &
             &                            lat_geo,         &
             &                            fr_land_topo,    &
             &                            hh_topo,         &
             &                            stdh_topo,       &
             &                            z0_topo,         &
             &                            lradtopo,        &
             &                            nhori)
      ENDIF

    ENDIF

  END SELECT

  CALL deallocate_topo_fields()

  DEALLOCATE (topo_startrow, topo_endrow, topo_startcolumn, topo_endcolumn)

  PRINT *,'============= topo_to_buffer done ==============='

END PROGRAM extpar_topo_to_buffer
