!+ Fortran module to aggregate GLOBE orogrphy data to the target grid 
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_1         2011/01/20 Hermann Asensio
!  correct calculation of standard deviation of height
! V1_2         2011/03/25 Hermann Asensio
!  update to support ICON refinement grids
!  bug fix in interpolation routine and handling of undefined GLOBE values
! V1_4         2011/04/21 Anne Roches
!  implementation of orography smoothing
! V1_7         2013/01/25 Guenther Zaengl 
!   Parallel threads for ICON and COSMO using Open-MP, 
!   Several bug fixes and optimizations for ICON search algorithm, 
!   particularly for the special case of non-contiguous domains; 
!   simplified namelist control for ICON
!   Potential bugfix in treatment of GLOBE data for COSMO -
!   no impact on results detected so far
! V2_0        2013/04/09 Martina Messmer
!  introduce ASTER topography for external parameters
!  Change all 'globe' to topo in globe_files, remove all 'globe' in 
!  change mo_GLOBE_data to mo_topo_data, globe_tiles_grid to 
!  topo_tiles_grid, globe_files to topo_files, globe_grid to
!  topo_grid and change ntiles_gl to ntiles to obtain a more 
!  dynamical code.  
! V2_0_3       2014/09/17 Burkhardt Rockel
!  Added use of directory information to access raw data files
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module to aggregate GLOBE orogrphy data to the target grid
!> \author Hermann Asensio
MODULE mo_agg_topo_icon

  USE mo_kind,                ONLY: wp, i4, i8
  USE mo_utilities_extpar,    ONLY: abort_extpar
  USE mo_io_units,            ONLY: filename_max
  USE mo_grid_structures,     ONLY: reg_lonlat_grid,             &
       &                            target_grid_def,             &
       &                            igrid_icon
  USE mo_topo_data,           ONLY: ntiles,         & !< there are 16/240 GLOBE/ASTER tiles 
       &                            max_tiles,      &
       &                            nc_tot,         & !< number of total GLOBE/ASTER columns un a latitude circle
       &                            nr_tot,         & !< total number of rows in GLOBE/ASTER data
       &                            get_fill_value, &
       &                            itopo_type,     &
       &                            topo_gl,        & 
       &                            topo_aster
  USE mo_topo_sso,            ONLY: calculate_sso, auxiliary_sso_parameter_icon
  USE mo_topo_routines,       ONLY: open_netcdf_TOPO_tile,       &
       &                            close_netcdf_TOPO_tile,      &
       &                            get_topo_data_block,         &
       &                            det_band_gd
  USE mo_target_grid_data,    ONLY: lon_geo, & !< longitude coordinates of the grid in the geographical system 
       &                            lat_geo, & !< latitude coordinates of the grid in the geographical system
       &                            search_res !< resolution of ICON grid search index list
  USE mo_icon_grid_data,      ONLY: icon_grid, icon_grid_region
  USE mo_icon_domain,         ONLY: icon_domain
  USE mo_topo_tg_fields,      ONLY: vertex_param
  USE mo_search_icongrid,     ONLY: walk_to_nc, find_nearest_vert
  USE mo_base_geometry,       ONLY: geographical_coordinates, cartesian_coordinates
  USE mo_additional_geometry, ONLY: gc2cc

  USE mo_math_constants,      ONLY: rad2deg, deg2rad
  USE mo_physical_constants,  ONLY: re !< av. radius of the earth [m]
  USE mo_bilinterpol,         ONLY: get_4_surrounding_raw_data_indices, &
       &                            calc_weight_bilinear_interpol, &
       &                            calc_value_bilinear_interpol
  USE mo_oro_filter,          ONLY: do_orosmooth

  USE mo_bilinterpol,         ONLY: get_4_surrounding_raw_data_indices, &
       &                            calc_weight_bilinear_interpol, &
       &                            calc_value_bilinear_interpol
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: agg_topo_data_to_target_grid_icon

CONTAINS

  SUBROUTINE agg_topo_data_to_target_grid_icon(topo_tiles_grid,                   &
       &                                       topo_grid,                         &
       &                                       tg,                                &
       &                                       topo_files,                        &
       &                                       lsso_param,                        &
       &                                       lscale_separation,                 &
       &                                       lfilter_oro,                       &
       &                                       ilow_pass_oro,                     &
       &                                       numfilt_oro,                       &
       &                                       eps_filter,                        &
       &                                       ifill_valley,                      &
       &                                       rfill_valley,                      &
       &                                       ilow_pass_xso,                     &
       &                                       numfilt_xso,                       &
       &                                       lxso_first,                        &
       &                                       rxso_mask,                         & 
       &                                       hh_target,                         &
       &                                       stdh_target,                       &
       &                                       fr_land_topo,                      &
       &                                       z0_topo,                           &
       &                                       no_raw_data_pixel,                 &
       &                                       theta_target,                      &
       &                                       aniso_target,                      &
       &                                       slope_target,                      &
       &                                       raw_data_orography_path,           &
       &                                       raw_data_scale_sep_orography_path, &
       &                                       scale_sep_files)



    TYPE(reg_lonlat_grid) ::                    topo_tiles_grid(1:ntiles)  !< raw data grid for the 16/36 GLOBE/ASTER tiles
    TYPE(target_grid_def),        INTENT(in) :: tg                         !< target grid description
    TYPE(reg_lonlat_grid) ::                    topo_grid                  !< raw data grid for the whole GLOBE/ASTER
    CHARACTER (len=filename_max), INTENT(IN) :: topo_files(1:max_tiles)    !< filenames globe/aster raw data
    LOGICAL,                      INTENT(in) :: lsso_param
    LOGICAL,                      INTENT(in) :: lscale_separation
    LOGICAL,                      INTENT(in) :: lfilter_oro                !< oro smoothing to be performed? (TRUE/FALSE) 
    INTEGER(i4),                  INTENT(in) :: ilow_pass_oro              !< type of oro smoothing and stencil width (1,4,5,6,8)
    INTEGER(i4),                  INTENT(in) :: numfilt_oro                !< number of applications of the filter
    REAL(wp),                     INTENT(in) :: eps_filter                 !< smoothing param ("strength" of the filtering)
    INTEGER(i4),                  INTENT(in) :: ifill_valley               !< fill valleys before or after oro smoothing 
    !  (1: before, 2: after)
    REAL(wp),                     INTENT(in) :: rfill_valley               !< mask for valley filling (threshold value)
    INTEGER(i4),                  INTENT(in) :: ilow_pass_xso              !< type of oro eXtra SmOothing for steep
    !  orography and stencil width (1,4,5,6,8)
    INTEGER(i4),                  INTENT(in) :: numfilt_xso                !< number of applications of the eXtra filter
    LOGICAL,                      INTENT(in) :: lxso_first                 !< eXtra SmOothing before or after oro
    !  smoothing? (TRUE/FALSE)
    REAL(wp),                     INTENT(in) :: rxso_mask                  !< mask for eXtra SmOothing (threshold value)

    REAL(wp),                     INTENT(out) :: hh_target(1:tg%ie,1:tg%je,1:tg%ke)
    REAL(wp),                     INTENT(out) :: stdh_target(1:tg%ie,1:tg%je,1:tg%ke)
    REAL(wp),                     INTENT(out) :: z0_topo(1:tg%ie,1:tg%je,1:tg%ke) 
    REAL(wp),                     INTENT(out) :: fr_land_topo(1:tg%ie,1:tg%je,1:tg%ke) 
    INTEGER (i8),                 INTENT(out) :: no_raw_data_pixel(1:tg%ie,1:tg%je,1:tg%ke)  


    REAL(wp), INTENT(out), OPTIONAL:: theta_target(1:tg%ie,1:tg%je,1:tg%ke) !< sso parameter, angle of principal axis
    REAL(wp), INTENT(out), OPTIONAL:: aniso_target(1:tg%ie,1:tg%je,1:tg%ke) !< sso parameter, anisotropie factor
    REAL(wp), INTENT(out), OPTIONAL:: slope_target(1:tg%ie,1:tg%je,1:tg%ke) !< sso parameter, mean slope

    CHARACTER(len=filename_max), INTENT(in), OPTIONAL :: scale_sep_files(1:max_tiles)  !< files globe/aster raw scale separated data
    CHARACTER(len=filename_max), INTENT(in), OPTIONAL :: raw_data_orography_path
    CHARACTER(len=filename_max), INTENT(in), OPTIONAL :: raw_data_scale_sep_orography_path

    ! local variables

    REAL (wp)    :: lon_topo(1:nc_tot)   !< longitude coordinates of the GLOBE grid
    REAL (wp)    :: lat_topo(1:nr_tot)   !< latititude coordinates of the GLOBE grid

    INTEGER (i4) :: nc_tot_p1
    !< ncid for the GLOBE/ASTER tiles, the netcdf files have to be opened by a previous call of open_netcdf_topo_tile
    INTEGER      :: ncids_topo(1:ntiles)  
    !< ncid for the GLOBE/ASTER scale separated tiles, the netcdf files have to be opened by a previous call of open_netcdf_topo_tile
    INTEGER      :: ncids_scale(1:ntiles)  

    INTEGER (i4) :: h_parallel(1:nc_tot)        !< one line with GLOBE/ASTER data
    INTEGER (i4) :: h_parallel_scale(1:nc_tot)  !< one line with GLOBE/ASTER scale separated data
    INTEGER (i4) :: h_3rows(1:nc_tot,1:3)       !< three rows with GLOBE/ASTER data
    INTEGER (i4) :: h_3rows_scale(1:nc_tot,1:3) !< three rows with GLOBE/ASTER scale separated data
    INTEGER (i4) :: hh(0:nc_tot+1,1:3)          !< topographic height for gradient calculations
    INTEGER (i4) :: hh_scale(0:nc_tot+1,1:3)    !< scale separated topographic height for gradient calculations

    REAL (wp)   :: hh_target_scale(1:tg%ie,1:tg%je,1:tg%ke)

    REAL (wp)   :: dhdxdx(1:nc_tot)  !< x-gradient square for one latitude row
    REAL (wp)   :: dhdydy(1:nc_tot)  !< y-gradient square for one latitude row
    REAL (wp)   :: dhdxdy(1:nc_tot)  !< dxdy for one latitude row

    REAL (wp)   :: hh1_target(1:tg%ie,1:tg%je,1:tg%ke)        !< mean height of grid element
    REAL (wp)   :: hh2_target(1:tg%ie,1:tg%je,1:tg%ke)        !< square mean height of grid element
    REAL (wp)   :: hh2_target_scale(1:tg%ie,1:tg%je,1:tg%ke)  !< square mean scale separated height of grid element
    !< squared difference between the filtered (scale separated) and original topography
    REAL (wp)   :: hh_sqr_diff(1:tg%ie,1:tg%je,1:tg%ke) 
    REAL (wp)   :: hsmooth(1:tg%ie,1:tg%je,1:tg%ke)  !< mean smoothed height of grid element

    REAL (wp)   :: h11(1:tg%ie,1:tg%je,1:tg%ke) !< help variables
    REAL (wp)   :: h12(1:tg%ie,1:tg%je,1:tg%ke) !< help variables
    REAL (wp)   :: h22(1:tg%ie,1:tg%je,1:tg%ke) !< help variables

    INTEGER (i8) :: ndata(1:tg%ie,1:tg%je,1:tg%ke)  !< number of raw data pixel with land point

    TYPE(geographical_coordinates) :: target_geo_co  !< structure for geographical coordinates of raw data pixel

    INTEGER (i4) :: undef_topo
    INTEGER (i4) :: default_topo

    INTEGER :: i,j ! counters

    INTEGER (i8) :: ie, je, ke  ! indices for grid elements
    INTEGER (i8), ALLOCATABLE :: ie_vec(:), iev_vec(:)  ! indices for target grid elements
    INTEGER (i8) :: i_vert, j_vert, k_vert ! indeces for ICON grid vertices
    INTEGER (i8) :: i1, i2

    INTEGER :: nv            ! counter
    INTEGER :: nt            ! counter
    INTEGER :: j_n, j_c, j_s ! counter for northern, central and southern row
    INTEGER :: j_new         ! counter for swapping indices j_n, j_c, j_s
    INTEGER :: mlat          ! row number for GLOBE data

    REAL (wp) :: dx, dy, dx0    !  grid distance for gradient calculation (in [m])
    REAL (wp) :: d2x, d2y       ! 2 times grid distance for gradient calculation (in [m])
    REAL (wp) :: row_lat(1:3)    ! latitude of the row for the topographic height array hh
    REAL (wp) :: znorm, znfi2sum, zarg ! help variables for the estiamtion of the variance
    REAL (wp) :: znorm_z0, zarg_z0 ! help variables for the estiamtion of the variance   
    REAL (wp) :: stdh_z0(1:tg%ie,1:tg%je,1:tg%ke)

    ! Some stuff for OpenMP parallelization

    INTEGER :: num_blocks, ib, il, blk_len, istartlon, iendlon, nlon_sub, ishift
    !$ INTEGER :: omp_get_max_threads, omp_get_thread_num, thread_id
    !$ INTEGER (i8), ALLOCATABLE :: start_cell_arr(:)

    TYPE(reg_lonlat_grid) :: ta_grid
    TYPE(icon_domain)     :: domain

    !< structure with definition of the target area grid (dlon must be the same as for the whole GLOBE/ASTER dataset)
    INTEGER (i4), ALLOCATABLE :: h_block(:,:) !< a block of GLOBE/ASTER altitude data
    INTEGER (i4), ALLOCATABLE :: h_block_scale(:,:) !< a block of GLOBE/ASTER altitude scale separated data
    INTEGER :: block_row_start
    INTEGER :: block_row
    INTEGER :: errorcode !< error status variable
    ! test with walk_to_nc at start
    INTEGER (i8) :: start_cell_id
    !< coordinates in cartesian system of point for which the nearest ICON grid cell is to be determined
    TYPE(cartesian_coordinates) :: target_cc_co     

    ! global data flag
    LOGICAL :: lskip
    REAL (wp) :: point_lon_geo         !< longitude coordinate in geographical system of input point 
    REAL (wp) :: point_lat_geo         !< latitude coordinate in geographical system of input point
    REAL (wp) :: point_lon, point_lat
    REAL (wp) :: topo_target_value     !< interpolated altitude from GLOBE data
    REAL (wp) :: fr_land_pixel         !< interpolated fr_land from GLOBE data
    ! variables for the "Erdmann Heise Formel"
    REAL (wp) :: dnorm                 !< scale factor 
    REAL (wp) :: zlnorm = 2250.        !< scale factor [m]
    REAL (wp) :: alpha  = 1.E-05       !< scale factor [1/m] 
    REAL (wp) :: factor                !< factor
    REAL (wp) :: zhp = 10.0            !< height of Prandtl-layer [m]
    REAL (wp) :: z0_topography         !< rougness length according to Erdmann Heise Formula

    CHARACTER(len=filename_max) :: topo_file_1
    CHARACTER(len=filename_max) :: scale_sep_file_1
    CHARACTER(len=filename_max) :: data_path, data_path_scale_sep
    
    nc_tot_p1 = nc_tot + 1

    IF (PRESENT(raw_data_orography_path)) THEN
      data_path = TRIM(raw_data_orography_path)
    ELSE
      data_path = "./"
    ENDIF

    IF (PRESENT(raw_data_scale_sep_orography_path)) THEN
      data_path_scale_sep = TRIM(raw_data_scale_sep_orography_path)
    ELSE
      data_path_scale_sep = "./"
    ENDIF
      
    topo_file_1 = TRIM(data_path)//TRIM(topo_files(1))
    IF (lscale_separation) THEN
      scale_sep_file_1 = TRIM(data_path_scale_sep)//TRIM(scale_sep_files(1))
    ENDIF

    ke = 1

    j_n = 1 ! index for northern row
    j_c = 2 ! index for central row
    j_s = 3 ! index for southern row

    CALL get_fill_value(topo_file_1,undef_topo)

    default_topo = 0

    SELECT CASE(itopo_type)
    CASE(topo_aster)
      hh = default_topo
      h_3rows = default_topo
    CASE(topo_gl)
      hh = undef_topo
      h_3rows = undef_topo
    END SELECT

    ! initialize some variables
    no_raw_data_pixel = 0
    ndata      = 0
    z0_topo     = 0.0
    hh_target   = 0.0
    hh1_target  = 0.0
    hh2_target  = 0.0
    stdh_target = 0.0

    IF (lscale_separation) THEN
      hh_target_scale  = 0.0
      hh2_target_scale = 0.0
      stdh_z0          = 0.0
    ENDIF

    IF (lsso_param) THEN
      theta_target = 0.0
      aniso_target = 0.0
      slope_target = 0.0
    ENDIF
    h11         = 0.0
    h12         = 0.0
    h22         = 0.0

    hsmooth     = 0.0

    IF (tg%igrid_type == igrid_icon) THEN ! Icon grid
      vertex_param%hh_vert = 0.0
      vertex_param%npixel_vert = 0
    ENDIF

    ! calculate the longitude coordinate of the GLOBE columns
    DO i = 1, nc_tot
      lon_topo(i) = topo_grid%start_lon_reg + (i-1) * topo_grid%dlon_reg
    ENDDO

    ! calculate the latitiude coordinate of the GLOBE columns
    DO j=1,nr_tot
      lat_topo(j) = topo_grid%start_lat_reg + (j-1) * topo_grid%dlat_reg
    ENDDO
    !HA debug:
    print *,'lat_topo(1): ', lat_topo(1)
    print *,'lat_topo(nr_tot) ', lat_topo(nr_tot)

    ALLOCATE(ie_vec(nc_tot),iev_vec(nc_tot))
    ie_vec(:)     = 0
    iev_vec(:)    = 0
    start_cell_id = 1

    nt = 1
    ! longitudinal distance between to topo grid elemtens at equator
    dx0 =  topo_tiles_grid(nt)%dlon_reg * deg2rad * re 
    print *, 'dx0: ',dx0
    ! latitudinal distance  between to topo grid elemtens ! note the negative increment, as direction of data from north to south
    dy = -1.0 * topo_tiles_grid(nt)%dlat_reg * deg2rad * re
    print *,'dy: ',dy
    d2y = 2. * dy

    print *,'open TOPO netcdf files'
    DO nt=1,ntiles
      CALL open_netcdf_TOPO_tile(TRIM(data_path)//TRIM(topo_files(nt)), ncids_topo(nt))
    ENDDO
    IF (lscale_separation) THEN
      DO nt=1,ntiles
        print*, 'scale_sep_files(nt): ', TRIM(scale_sep_files(nt))
        CALL open_netcdf_TOPO_tile(TRIM(data_path_scale_sep)//TRIM(scale_sep_files(nt)), ncids_scale(nt))
      ENDDO
    ENDIF

    mlat = 1
    block_row_start = mlat

    CALL det_band_gd(topo_grid,block_row_start, ta_grid)
    PRINT *,'first call of det_band_gd'
    PRINT *,'ta_grid: ',ta_grid

    IF (ALLOCATED(h_block)) THEN
      DEALLOCATE(h_block, STAT=errorcode)
      IF (errorcode /= 0) CALL abort_extpar('Cant deallocate the h_block')
    ENDIF
    ALLOCATE (h_block(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg), STAT=errorcode)
    IF (errorcode /= 0) CALL abort_extpar('Cant allocate h_block')

    CALL get_topo_data_block(topo_file_1,       &
         &                   ta_grid,           &
         &                   topo_tiles_grid,   &
         &                   ncids_topo,        &
         &                   h_block)

    IF (lscale_separation) THEN
      IF (ALLOCATED(h_block_scale)) THEN
        DEALLOCATE(h_block_scale, STAT=errorcode)
        IF (errorcode/=0) CALL abort_extpar('Cant deallocate the h_block_scale')
      ENDIF
      ALLOCATE (h_block_scale(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg), STAT=errorcode)
      IF (errorcode/=0) CALL abort_extpar('Cant allocate h_block_scale')

      CALL get_topo_data_block(scale_sep_file_1, &
           &                   ta_grid,          &
           &                   topo_tiles_grid,  &
           &                   ncids_scale,      &
           &                   h_block_scale)
    ENDIF

    block_row = 0 

    ! Determine start and end longitude of search

    istartlon = 1
    iendlon = nc_tot

    DO i = 1, nc_tot
      point_lon = lon_topo(i)
      IF (point_lon < tg%minlon) istartlon = i + 1
      IF (point_lon > tg%maxlon) THEN
        iendlon = i - 1
        EXIT
      ENDIF
    ENDDO
    nlon_sub = iendlon - istartlon + 1

    num_blocks = 1
    !$ num_blocks = omp_get_max_threads()
    IF (MOD(nlon_sub,num_blocks)== 0) THEN
      blk_len = nlon_sub/num_blocks
    ELSE
      blk_len = nlon_sub/num_blocks + 1
    ENDIF
    !$ ALLOCATE(start_cell_arr(num_blocks))
    !$ start_cell_arr(:) = 1
    PRINT*, 'nlon_sub, num_blocks, blk_len: ',nlon_sub, num_blocks, blk_len


    print *,'Start loop over TOPO rows'
    !-----------------------------------------------------------------------------
    topo_rows: DO mlat = 1, nr_tot
      !-----------------------------------------------------------------------------
      !-----------------------------------------------------------------------------
      IF (MOD(mlat,100)==0) PRINT *, 'TOPO row:', mlat
      block_row= block_row + 1
      IF ((block_row > ta_grid%nlat_reg).AND.(mlat<nr_tot)) THEN ! read in new block
        block_row_start = mlat + 1
        block_row = 1
        CALL det_band_gd(topo_grid,block_row_start, ta_grid)
        PRINT *,'next call of det_band_gd'
        PRINT *,'ta_grid: ',ta_grid
        IF (ALLOCATED(h_block)) THEN
          DEALLOCATE(h_block, STAT=errorcode)
          IF (errorcode/=0) CALL abort_extpar('Cant deallocate the h_block')
        ENDIF
        ALLOCATE (h_block(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg), STAT=errorcode)
        IF (errorcode/=0) CALL abort_extpar('Cant allocate h_block')
        CALL get_topo_data_block(topo_file_1,     &
             &                   ta_grid,         &
             &                   topo_tiles_grid, &
             &                   ncids_topo,      &
             &                   h_block)
        IF (lscale_separation) THEN
          IF (ALLOCATED(h_block_scale)) THEN
            DEALLOCATE(h_block_scale, STAT=errorcode)
            IF (errorcode/=0) CALL abort_extpar('Cant deallocate the h_block_scale')
          ENDIF
          ALLOCATE (h_block_scale(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg), STAT=errorcode)
          IF (errorcode/=0) CALL abort_extpar('Cant allocate h_block_scale')
          CALL get_topo_data_block(scale_sep_file_1,    &
               &                       ta_grid,         &
               &                       topo_tiles_grid, &
               &                       ncids_scale,     &
               &                       h_block_scale)
        ENDIF
      ENDIF

      IF (mlat==1) THEN
        h_parallel(1:nc_tot) = h_block(1:nc_tot,block_row)
        row_lat(j_c) = topo_grid%start_lat_reg + (mlat-1) * topo_grid%dlat_reg

        h_3rows(1:nc_tot,j_c) = h_parallel(1:nc_tot)  ! put data to "central row"
        hh(1:nc_tot,j_c)      = h_parallel(1:nc_tot)  ! put data to "central row"
        hh(0,j_c)             = h_parallel(nc_tot)    ! western wrap at -180/180 degree longitude
        hh(nc_tot_p1,j_c)     = h_parallel(1)         ! eastern wrap at -180/180 degree longitude

        IF (lscale_separation) THEN
          h_parallel_scale(1:nc_tot)  = h_block_scale(1:nc_tot,block_row)
          h_3rows_scale(1:nc_tot,j_c) = h_parallel_scale(1:nc_tot)  ! put data to "central row"
          hh_scale(1:nc_tot,j_c)      = h_parallel_scale(1:nc_tot)  ! put data to "central row"
          hh_scale(0,j_c)             = h_parallel_scale(nc_tot)    ! western wrap at -180/180 degree longitude
          hh_scale(nc_tot_p1,j_c)     = h_parallel_scale(1)         ! eastern wrap at -180/180 degree longitude
        ENDIF

        block_row = block_row + 1

      ENDIF

      row_lat(j_s) = topo_grid%start_lat_reg + mlat * topo_grid%dlat_reg  !  ((mlat+1)-1)

      lskip = .FALSE.
      IF (row_lat(j_s) > tg%maxlat .OR. row_lat(j_s) < tg%minlat) lskip = .TRUE.
      ! read raw data south of "central" row except when you are at the most southern raw data line
      IF (mlat /= nr_tot) THEN 
        h_parallel(1:nc_tot)  = h_block(1:nc_tot,block_row)
        h_3rows(1:nc_tot,j_s) = h_parallel(1:nc_tot)
        hh(1:nc_tot,j_s)      = h_parallel(1:nc_tot) ! put data to "southern row"
        hh(0,j_s)             = h_parallel(nc_tot)   ! western wrap at -180/180 degree longitude
        hh(nc_tot_p1,j_s)     = h_parallel(1)        ! eastern wrap at -180/180 degree longitude

        IF (lscale_separation) THEN
          h_parallel_scale(1:nc_tot) = h_block_scale(1:nc_tot,block_row)
          h_3rows_scale(1:nc_tot,j_s) = h_parallel_scale(1:nc_tot)
          hh_scale(1:nc_tot,j_s)      = h_parallel_scale(1:nc_tot)  ! put data to "central row"
          hh_scale(0,j_s)             = h_parallel_scale(nc_tot)    ! western wrap at -180/180 degree longitude
          hh_scale(nc_tot_p1,j_s)     = h_parallel_scale(1)         ! eastern wrap at -180/180 degree longitude
        ENDIF
      ENDIF

      IF (lskip) THEN
        ! swap indices of the hh array for next data row before skipping the loop
        j_new = j_n   ! the new data will be written in the former "northern" array
        j_n   = j_c   ! the "center" row will become "northern" row
        j_c   = j_s   ! the "southern" row will become "center" row
        j_s   = j_new ! the new data will be written in the "southern" row
        CYCLE topo_rows
      ENDIF

      dx  = dx0 * COS(row_lat(j_c) * deg2rad)  ! longitudinal distance between to GLOBE grid elemtens
      d2x = 2. * dx
      d2y = 2. * dy
      IF (mlat==1) THEN ! most northern row of raw data
        ! j_n = j_c  ! put the index of "northern row" to the same index as "central row"
        hh(:,j_n) = hh(:,j_c)
        IF (lscale_separation) THEN
          hh_scale(:,j_n) = hh_scale(:,j_c)
        ENDIF
        d2y = dy   ! adjust d2y in this case too
      ELSEIF (mlat==nr_tot) THEN ! most southern row of raw data
        j_s = j_c  ! put the index of "southern row" to the same index as "central row"
        d2y = dy   ! adjust d2y in this case too
      ENDIF

      ! set undefined values to 0 altitude (default)
      WHERE (hh == undef_topo)  
        hh = default_topo
      END WHERE
      WHERE (h_parallel == undef_topo)
        h_parallel = default_topo
      END WHERE

      IF (lscale_separation) THEN
        WHERE (hh_scale == undef_topo)  
          hh_scale = default_topo
        END WHERE
        WHERE (h_parallel_scale == undef_topo)
          h_parallel_scale = default_topo
        END WHERE
      ENDIF

      IF (lsso_param) THEN
        IF (lscale_separation) THEN
          CALL auxiliary_sso_parameter_icon(d2x,d2y,j_n,j_c,j_s,hh_scale,dhdxdx,dhdydy,dhdxdy)
        ELSE
          CALL auxiliary_sso_parameter_icon(d2x,d2y,j_n,j_c,j_s,hh,dhdxdx,dhdydy,dhdxdy)
        ENDIF
      ENDIF

      ie_vec(istartlon:iendlon)  = 0
      iev_vec(istartlon:iendlon) = 0

      point_lat = row_lat(j_c)

!$OMP PARALLEL DO PRIVATE(ib,il,i,i1,i2,ishift,point_lon,thread_id,start_cell_id,target_geo_co,target_cc_co)
      DO ib = 1, num_blocks

        !$   thread_id = omp_get_thread_num()+1
        !$   start_cell_id = start_cell_arr(thread_id)
        ishift = istartlon-1+(ib-1)*blk_len

        ! loop over one latitude circle of the raw data
        columns1: DO il = 1, blk_len
          i = ishift+il
          IF (i > iendlon) CYCLE columns1

          ! find the corresponding target grid indices
          point_lon = lon_topo(i) 

          ! Reset start cell when entering a new row or when the previous data point was outside
          ! the model domain
          IF (il == 1 .OR. start_cell_id == 0) THEN
            i1 = NINT(point_lon*search_res)
            i2 = NINT(point_lat*search_res)
            start_cell_id = tg%search_index(i1,i2)
            IF (start_cell_id == 0) EXIT ! in this case, the whole row is empty; may happen with merged (non-contiguous) domains
          ENDIF

          target_geo_co%lon = point_lon * deg2rad ! note that the ICON coordinates do not have the unit degree but radians
          target_geo_co%lat = point_lat * deg2rad
          target_cc_co = gc2cc(target_geo_co)
          CALL walk_to_nc(icon_grid_region,            &
               &          target_cc_co,                &
               &          start_cell_id,               &
               &          icon_grid%nvertex_per_cell,  &
               &          icon_grid%nedges_per_vertex, &
               &          ie_vec(i))

          ! additional get the nearest vertex index for accumulating height values there
          IF (ie_vec(i) /= 0_i8) THEN
            CALL  find_nearest_vert(icon_grid_region,           &
                 &                  target_cc_co,               &
                 &                  ie_vec(i),                  &
                 &                  icon_grid%nvertex_per_cell, &
                 &                  iev_vec(i))
          ENDIF

        ENDDO columns1
        !$   start_cell_arr(thread_id) = start_cell_id
      ENDDO
!$OMP END PARALLEL DO


!$OMP PARALLEL DO PRIVATE(i,ie,i_vert,point_lon)
      DO i = istartlon, iendlon

        ! call here the attribution of raw data pixel to target grid for different grid types
        ie = ie_vec(i)
        je = 1
        ke = 1

        ! get the nearest vertex index for accumulating height values there

        ! aggregate the vertex parameter here
        i_vert = iev_vec(i)
        j_vert = 1
        k_vert = 1
        IF ((i_vert /=0)) THEN ! raw data pixel within target grid
          vertex_param%npixel_vert(i_vert,j_vert,k_vert) = vertex_param%npixel_vert(i_vert,j_vert,k_vert) + 1
          vertex_param%hh_vert(i_vert,j_vert,k_vert)     = vertex_param%hh_vert(i_vert,j_vert,k_vert) +  h_parallel(i)
        ENDIF

        IF ((ie /= 0).AND.(je/=0).AND.(ke/=0))THEN 
          ! raw data pixel within target grid, see output of routine find_rotated_lonlat_grid_element_index
          no_raw_data_pixel(ie,je,ke) = no_raw_data_pixel(ie,je,ke) + 1
          !- summation of variables

          SELECT CASE(itopo_type)
          CASE(topo_aster)
            IF (h_3rows(i,j_c) /= default_topo) THEN       
              ndata(ie,je,ke)      = ndata(ie,je,ke) + 1
              hh_target(ie,je,ke)  = hh_target(ie,je,ke)  +  h_3rows(i,j_c)
              hh2_target(ie,je,ke) = hh2_target(ie,je,ke) + (h_3rows(i,j_c) * h_3rows(i,j_c))

              IF (lscale_separation) THEN
                hh_target_scale(ie,je,ke)  = hh_target_scale(ie,je,ke)  +  h_3rows_scale(i,j_c)
                hh2_target_scale(ie,je,ke) = hh2_target_scale(ie,je,ke) + (h_3rows_scale(i,j_c) * h_3rows_scale(i,j_c)) 
                hh_sqr_diff(ie,je,ke)      = hh_sqr_diff(ie,je,ke)      + (h_3rows(i,j_c)       - h_3rows_scale(i,j_c))**2
              ENDIF

              IF (lsso_param) THEN
                h11(ie,je,ke)        = h11(ie,je,ke) + dhdxdx(i)
                h12(ie,je,ke)        = h12(ie,je,ke) + dhdxdy(i)
                h22(ie,je,ke)        = h22(ie,je,ke) + dhdydy(i)
              ENDIF
            ENDIF
          CASE(topo_gl)
            IF (h_3rows(i,j_c) /= undef_topo) THEN            
              ndata(ie,je,ke)      = ndata(ie,je,ke) + 1
              hh_target(ie,je,ke)  = hh_target(ie,je,ke) + h_3rows(i,j_c)
              hh2_target(ie,je,ke) = hh2_target(ie,je,ke) + (h_3rows(i,j_c) * h_3rows(i,j_c))

              IF (lscale_separation) THEN
                hh_target_scale(ie,je,ke)  = hh_target_scale(ie,je,ke)  +  h_3rows_scale(i,j_c)
                hh2_target_scale(ie,je,ke) = hh2_target_scale(ie,je,ke) + (h_3rows_scale(i,j_c) * h_3rows_scale(i,j_c))
                hh_sqr_diff(ie,je,ke)      = hh_sqr_diff(ie,je,ke)      + (h_3rows(i,j_c)       - h_3rows_scale(i,j_c))**2
              ENDIF

              IF (lsso_param) THEN
                h11(ie,je,ke)        = h11(ie,je,ke) + dhdxdx(i)
                h12(ie,je,ke)        = h12(ie,je,ke) + dhdxdy(i)
                h22(ie,je,ke)        = h22(ie,je,ke) + dhdydy(i)
              ENDIF
            ENDIF
          END SELECT
        ENDIF

      ENDDO ! loop over one latitude circle of the raw data
!$OMP END PARALLEL DO

      ! swap indices of the hh array for next data row
      j_new = j_n ! the new data will be written in the former "northern" array
      j_n = j_c   ! the "center" row will become "northern" row
      j_c = j_s   ! the "southern" row will become "center" row
      j_s = j_new ! the new data will be written in the "southern" row

      !-----------------------------------------------------------------------------
    ENDDO topo_rows
    !-----------------------------------------------------------------------------

    DEALLOCATE(ie_vec,iev_vec)
    !$     DEALLOCATE(start_cell_arr)

    print *,'loop over topo_rows done'

    PRINT *,'Maximum number of TOPO raw data pixel in a target grid element: '
    PRINT *,'MAXVAL(no_raw_data_pixel): ', MAXVAL(no_raw_data_pixel)
    PRINT *,'Index of target grid element: ', MAXLOC(no_raw_data_pixel)

    PRINT *,'Maximum number of TOPO land pixel in a target grid element: '
    PRINT *,'MAXVAL(ndata): ', MAXVAL(ndata)
    PRINT *,'Index of target grid element: ', MAXLOC(ndata)

    PRINT *,'Minimal number of TOPO raw data pixel in a target grid element: '
    PRINT *,'MINVAL(no_raw_data_pixel): ', MINVAL(no_raw_data_pixel)
    PRINT *,'Index of target grid element: ', MINLOC(no_raw_data_pixel)

    PRINT *,'Minimal number of TOPO land pixel in a target grid element: '
    PRINT *,'MINVAL(ndata): ', MINVAL(ndata)
    PRINT *,'Index of target grid element: ', MINLOC(ndata)

    hh1_target = hh_target ! save values of hh_target for computations of standard deviation

    print *,'Average height'
    ! Average height
    DO ke=1, tg%ke
      DO je=1, tg%je
        DO ie=1, tg%ie
          IF (no_raw_data_pixel(ie,je,ke) /= 0) THEN ! avoid division by zero for small target grids
            hh_target(ie,je,ke) = hh_target(ie,je,ke)/no_raw_data_pixel(ie,je,ke) 
            ! average height, oceans point counted as 0 height
            fr_land_topo(ie,je,ke) =  REAL(ndata(ie,je,ke),wp) / REAL(no_raw_data_pixel(ie,je,ke),wp) ! fraction land
          ELSE
            hh_target(ie,je,ke) = REAL(default_topo)
            fr_land_topo(ie,je,ke) = 0.0
          ENDIF
#ifdef DEBUG
          IF (  fr_land_topo(ie,je,ke) > 0 .AND. domain%cells%sea_land_mask(ie) < 0 ) THEN
             WRITE(*,*) "Warning: GLOBE land on ICON water at ICON cell index ", ie 
          ENDIF
          IF (  fr_land_topo(ie,je,ke) == 0 .AND. domain%cells%sea_land_mask(ie) > 0 ) THEN
             WRITE(*,*) "Warning: GLOBE water on ICON land at ICON cell index ", ie 
          ENDIF
#endif
        ENDDO
      ENDDO
    ENDDO

    hsmooth = hh_target

    IF (lfilter_oro) THEN
      CALL do_orosmooth(tg,               &
           &            hh_target,        &
           &            fr_land_topo,     &
           &            lfilter_oro,      &
           &            ilow_pass_oro,    &
           &            numfilt_oro,      &
           &            eps_filter,       &
           &            ifill_valley,     &
           &            rfill_valley,     &
           &            ilow_pass_xso,    &
           &            numfilt_xso,      &
           &            lxso_first,       &
           &            rxso_mask,        &
           &            hsmooth)
    ENDIF

    print *,'Average height for vertices'
    ! Average height for vertices
    DO ke=1, 1
      DO je=1, 1
        DO ie=1, icon_grid_region%nverts          
          IF (vertex_param%npixel_vert(ie,je,ke) /= 0) THEN ! avoid division by zero for small target grids
            vertex_param%hh_vert(ie,je,ke) = vertex_param%hh_vert(ie,je,ke)/vertex_param%npixel_vert(ie,je,ke) ! average height
          ELSE
            vertex_param%hh_vert(ie,je,ke) = REAL(default_topo)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    
    print *,'Standard deviation of height'
    !     Standard deviation of height.
    DO ke=1, tg%ke
      DO je=1, tg%je
        DO ie=1, tg%ie
          ! estimation of variance
          IF (no_raw_data_pixel(ie,je,ke) > 1) THEN
            znorm_z0 = 1.0/(no_raw_data_pixel(ie,je,ke)-1)
            znorm    = 1.0/(no_raw_data_pixel(ie,je,ke)*(no_raw_data_pixel(ie,je,ke)-1))
          ELSE
            znorm_z0 = 0.0
            znorm    = 0.0
          ENDIF

          IF (lscale_separation) THEN
            ! Standard deviation between filtred and un-filtred raw data
            ! (used to compute z0 later on)
            zarg_z0 = znorm_z0 * hh_sqr_diff(ie,je,ke)
            zarg_z0 = MAX(zarg_z0,0.0_wp) ! truncation errors may cause zarg_sso < 0.0
            stdh_z0(ie,je,ke) = SQRT(zarg_z0)

            ! Standard deviation between target grid and filtered raw data
            ! (used to compute SSO parameters later on)
            IF (lfilter_oro) THEN
              zarg = znorm_z0 * (hh2_target_scale(ie,je,ke) -               &   
                   & 2.0 * hsmooth(ie,je,ke) * hh_target_scale(ie,je,ke) +  &
                   & no_raw_data_pixel(ie,je,ke) * hsmooth(ie,je,ke)**2)
            ELSE
              zarg = znorm_z0 * (hh2_target_scale(ie,je,ke) -                &   
                   & 2.0 * hh_target(ie,je,ke) * hh_target_scale(ie,je,ke) + &
                   & no_raw_data_pixel(ie,je,ke) * hh_target(ie,je,ke)**2)
            ENDIF

          ELSE
            ! Standard deviation between target grid and raw data
            ! (used to compute both z0 and SSO parameters later on)
            IF (lfilter_oro) THEN
!!!!! standard deviation of height using oro filt !!!!!
              zarg = znorm_z0 * (hh2_target(ie,je,ke) -                &
                   & 2.0 * hsmooth(ie,je,ke) * hh1_target(ie,je,ke) +  &
                   & no_raw_data_pixel(ie,je,ke) * hsmooth(ie,je,ke)**2)
            ELSE
              znfi2sum = no_raw_data_pixel(ie,je,ke) * hh2_target(ie,je,ke) 
              zarg     = ( znfi2sum - (hh1_target(ie,je,ke)*hh1_target(ie,je,ke))) * znorm
            ENDIF
          ENDIF
          zarg = MAX(zarg,0.0_wp) ! truncation errors may cause zarg < 0.0
          stdh_target(ie,je,ke) = SQRT(zarg)
        ENDDO
      ENDDO
    ENDDO
    
    IF (lsso_param) THEN
      
      CALL calculate_sso(tg,no_raw_data_pixel,    &
           &             h11,h12,h22,stdh_target, &
           &             theta_target,            &
           &             aniso_target,            &
           &             slope_target)

    ENDIF
    !----------------------------------------------------------------------------------
    ! calculate roughness length
    ! first zo_topo with "Erdmann Heise formula"
    !----------------------------------------------------------------------------------

    dnorm = 60000.         ! dummy value for normation of Erdmann Heise formula
    !---------------------------------------------------------------------------------
    ! Erdman Heise Formel
    !---------------------------------------------------------------------------------
    factor= alpha*ATAN(dnorm/zlnorm) !  alpha  = 1.E-05 [1/m] ,  zlnorm = 2250 [m]  
    DO ke=1, tg%ke
      DO je=1, tg%je
        DO ie=1, tg%ie
          IF (lscale_separation) THEN
            z0_topography = factor*stdh_z0(ie,je,ke)**2
          ELSE
            z0_topography = factor*stdh_target(ie,je,ke)**2
          ENDIF
          z0_topography = MIN(z0_topography,zhp-1.0_wp)
          z0_topo(ie,je,ke) = z0_topography
        ENDDO
      ENDDO
    ENDDO

    ! set the orography variable hh_target to the smoothed orography variable
    ! hsmooth in case of orogrpahy smoothing in extpar

    IF (lfilter_oro) THEN
      hh_target (:,:,:) = hsmooth (:,:,:)
    ENDIF

    ! bilinear interpolation of the orography in case of target grid points having
    ! no corresponding points in GLOBE

    DO ke=1, tg%ke
      DO je=1, tg%je
        DO ie=1, tg%ie
          IF (no_raw_data_pixel(ie,je,ke) == 0) THEN  ! bilinear interpolation to target grid
            
            point_lon_geo = lon_geo(ie,je,ke)
            point_lat_geo = lat_geo(ie,je,ke)
            
            CALL bilinear_interpol_topo_to_target_point(data_path,               &
                 &                                      topo_files,              &
                 &                                      topo_grid,               &     
                 &                                      topo_tiles_grid,         &
                 &                                      ncids_topo,              &
                 &                                      lon_topo,                &
                 &                                      lat_topo,                &
                 &                                      point_lon_geo,           &
                 &                                      point_lat_geo,           &
                 &                                      fr_land_pixel,           &
                 &                                      topo_target_value)

            fr_land_topo(ie,je,ke) = fr_land_pixel
            hh_target(ie,je,ke) = topo_target_value
            
            IF (lsso_param) THEN            
              theta_target(ie,je,ke) = 0.0
              aniso_target(ie,je,ke) = 0.0
              slope_target(ie,je,ke) = 0.0
            ENDIF
            stdh_target(ie,je,ke)  = 0.0             
            z0_topo(ie,je,ke)      = 0.0
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    je=1
    ke=1
    DO nv=1, icon_grid_region%nverts
      IF (vertex_param%npixel_vert(nv,je,ke) == 0) THEN ! interpolate from raw data in this case
        point_lon_geo =  rad2deg * icon_grid_region%verts%vertex(nv)%lon
        point_lat_geo =  rad2deg * icon_grid_region%verts%vertex(nv)%lat
        
        CALL bilinear_interpol_topo_to_target_point(data_path,               &
             &                                      topo_files,              &
             &                                      topo_grid,               &
             &                                      topo_tiles_grid,         &
             &                                      ncids_topo,              &
             &                                      lon_topo,                &
             &                                      lat_topo,                &
             &                                      point_lon_geo,           &
             &                                      point_lat_geo,           &
             &                                      fr_land_pixel,           &
             &                                      topo_target_value)
        
        vertex_param%hh_vert(nv,je,ke) = topo_target_value
      ENDIF
    ENDDO
    
    DO nt = 1, ntiles
      CALL close_netcdf_TOPO_tile(ncids_topo(nt))
    ENDDO
    IF (lscale_separation) THEN
      DO nt = 1, ntiles
        CALL close_netcdf_TOPO_tile(ncids_scale(nt))
      ENDDO
    ENDIF
    PRINT *,'TOPO netcdf files closed'
    PRINT *,'Subroutine agg_topo_data_to_target_grid done'
  END SUBROUTINE agg_topo_data_to_target_grid_icon

  !----------------------------------------------------------------------------------------------------------------

  !> subroutine for bilenar interpolation from GLOBE data (regular lonlat grid) to a single target point
  !!
  !! the GLOBE data are passed to the subroutine in the topo_data_block 2D-Array, which is re-read 
  !! from the raw data file if the target point is out of the range of the data block. 
  !! (If the data block is not too small, repeated I/O to the hard disk is avoided, reading from memory is much faster.)
  !! 
  !! the definition of the regular lon-lat grid requires 
  !! - the coordinates of the north-western point of the domain ("upper left") startlon_reg_lonlat and startlat_reg_lonlat
  !! - the increment dlon_reg_lonlat and dlat_reg_lonlat(implict assuming that the grid definiton goes 
  !!   from the west to the east and from the north to the south)
  !! - the number of grid elements nlon_reg_lonlat and nlat_reg_lonlat for both directions
  SUBROUTINE bilinear_interpol_topo_to_target_point(raw_data_orography_path, &
       &                                            topo_files,              &
       &                                            topo_grid,               &
       &                                            topo_tiles_grid,         &
       &                                            ncids_topo,              &
       &                                            lon_topo,                &
       &                                            lat_topo,                &
       &                                            point_lon_geo,           &
       &                                            point_lat_geo,           &
       &                                            fr_land_pixel,           &
       &                                            topo_target_value)

    CHARACTER(len=filename_max), INTENT(in)     :: topo_files(1:max_tiles)
    TYPE(reg_lonlat_grid), INTENT(in) :: topo_grid                 !< raw data grid for the whole GLOBE/ASTER dataset
    TYPE(reg_lonlat_grid), INTENT(in) :: topo_tiles_grid(1:ntiles) !< raw data grid for the 16/36 GLOBE/ASTER tiles
    !< ncid for the topo tiles, the netcdf files have to be opened by a previous call of open_netcdf_GLOBE_tile
    INTEGER (i4), INTENT(in)     :: ncids_topo(1:ntiles)  
    
    REAL (wp), INTENT(in) :: lon_topo(1:nc_tot)   !< longitude coordinates of the GLOBE grid
    REAL (wp), INTENT(in) :: lat_topo(1:nr_tot)   !< latititude coordinates of the GLOBE grid
    REAL (wp), INTENT(in) :: point_lon_geo       !< longitude coordinate in geographical system of input point 
    REAL (wp), INTENT(in) :: point_lat_geo       !< latitude coordinate in geographical system of input point
    REAL (wp), INTENT(out) :: fr_land_pixel  !< interpolated fr_land from GLOBE data
    REAL (wp), INTENT(out) :: topo_target_value  !< interpolated altitude from GLOBE data

    ! local variables
    INTEGER (i4), ALLOCATABLE :: h_block(:,:) !< a block of GLOBE altitude data
    TYPE(reg_lonlat_grid) :: ta_grid 
    !< structure with definition of the target area grid (dlon must be the same as for the whole GLOBE dataset)
    INTEGER (i8) :: western_column     !< the index of the western_column of data to read in
    INTEGER (i8) :: eastern_column     !< the index of the eastern_column of data to read in
    INTEGER (i8) :: northern_row       !< the index of the northern_row of data to read in
    INTEGER (i8) :: southern_row       !< the index of the southern_row of data to read in
    REAL (wp)   :: bwlon  !< weight for bilinear interpolation
    REAL (wp)   :: bwlat  !< weight for bilinear interpolation
    REAL (wp)   :: topo_point_sw       !< value of the GLOBE raw data pixel south west
    REAL (wp)   :: topo_point_se       !< value of the GLOBE raw data pixel south east
    REAL (wp)   :: topo_point_ne       !< value of the GLOBE raw data pixel north east
    REAL (wp)   :: topo_point_nw       !< value of the GLOBE raw data pixel north west
    INTEGER :: errorcode
    LOGICAL :: gldata=.TRUE. ! GLOBE data are global
    INTEGER (i4) :: undef_topo
    INTEGER (i4) :: default_topo

    CHARACTER(len=filename_max) :: topo_file_1   
    CHARACTER(len=filename_max) :: raw_data_orography_path !_br 26.09.14

    topo_file_1 = TRIM(raw_data_orography_path)//TRIM(topo_files(1))

    CALL get_fill_value(topo_file_1,undef_topo)

    default_topo = 0

    ! get four surrounding raw data indices
    CALL  get_4_surrounding_raw_data_indices(topo_grid,      &
         &                                   lon_topo,       &
         &                                   lat_topo,       &
         &                                   gldata,         &
         &                                   point_lon_geo,  &
         &                                   point_lat_geo,  &
         &                                   western_column, &
         &                                   eastern_column, &
         &                                   northern_row,   &
         &                                   southern_row)
    !print *,'western_column, eastern_column, northern_row, southern_row'  
    !print *, western_column, eastern_column, northern_row, southern_row  


    ta_grid%dlon_reg = topo_grid%dlon_reg
    ta_grid%dlat_reg = topo_grid%dlat_reg
    ta_grid%nlon_reg = eastern_column - western_column + 1
    ta_grid%nlat_reg = southern_row - northern_row + 1
    ta_grid%start_lon_reg = lon_topo(western_column)
    ta_grid%end_lon_reg = lon_topo(eastern_column)
    ta_grid%start_lat_reg = lat_topo(northern_row)
    ta_grid%end_lat_reg  = lat_topo(southern_row) 

    ! calculate weight for bilinear interpolation
    CALL calc_weight_bilinear_interpol(point_lon_geo,            &
         &                             point_lat_geo,            &
         &                             lon_topo(western_column), &
         &                             lon_topo(eastern_column), &
         &                             lat_topo(northern_row),   &
         &                             lat_topo(southern_row),   &
         &                             bwlon,                    &
         &                             bwlat)

    ALLOCATE (h_block(western_column:eastern_column,northern_row:southern_row), STAT=errorcode)
    IF (errorcode/=0) CALL abort_extpar('Cant allocate h_block')
    CALL get_topo_data_block(topo_file_1,     &
         &                   ta_grid,         &
         &                   topo_tiles_grid, &
         &                   ncids_topo,      & 
         &                   h_block)
    ! check for undefined GLOBE data, which indicate ocean grid element

    IF ( h_block(western_column,southern_row) == undef_topo) THEN
      topo_point_sw = 0.0
      h_block(western_column,southern_row) = default_topo
    ELSE
      topo_point_sw = 1.0
    ENDIF

    IF ( h_block(eastern_column,southern_row) == undef_topo) THEN
      topo_point_se = 0.0
      h_block(eastern_column,southern_row) = default_topo
    ELSE
      topo_point_se = 1.0
    ENDIF

    IF ( h_block(eastern_column,northern_row) == undef_topo) THEN
      topo_point_ne = 0.0
      h_block(eastern_column,northern_row) = default_topo
    ELSE
      topo_point_ne = 1.0
    ENDIF

    IF ( h_block(western_column,northern_row) == undef_topo) THEN
      topo_point_nw = 0.0
      h_block(western_column,northern_row) = default_topo 
    ELSE
      topo_point_nw = 1.0
    ENDIF


    ! perform the interpolation, first for fraction land
    fr_land_pixel = calc_value_bilinear_interpol(bwlon,         &
         &                                       bwlat,         &
         &                                       topo_point_sw, &
         &                                       topo_point_se, &
         &                                       topo_point_ne, &
         &                                       topo_point_nw)

    topo_point_sw = h_block(western_column,southern_row)
    topo_point_se = h_block(eastern_column,southern_row)
    topo_point_ne = h_block(eastern_column,northern_row)
    topo_point_nw = h_block(western_column,northern_row)

    ! perform the interpolation for height
    topo_target_value = calc_value_bilinear_interpol(bwlon,         &
         &                                           bwlat,         &
         &                                           topo_point_sw, &
         &                                           topo_point_se, &
         &                                           topo_point_ne, &
         &                                           topo_point_nw)

  END SUBROUTINE bilinear_interpol_topo_to_target_point

END MODULE mo_agg_topo_icon


