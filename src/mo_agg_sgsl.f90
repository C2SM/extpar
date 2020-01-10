!+ Fortran module to aggregate subgrid-scale slope data to the target grid
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V4_0         2016/07/28 Daniel Luethi
!  Added use of directory information to access raw data files
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module to aggregate subgrid-scale slope data to the target grid
!> \author Daniel Luethi
MODULE mo_agg_sgsl
  !> kind parameters are defined in MODULE data_parameters
  USE mo_kind, ONLY: wp
  USE mo_kind, ONLY: i4
  USE mo_kind, ONLY: i8
  !> abort_extpar defined in MODULE utilities_extpar
  USE mo_utilities_extpar, ONLY: abort_extpar
  USE mo_io_units,          ONLY: filename_max
  !> data type structures form module GRID_structures
  USE mo_grid_structures, ONLY: reg_lonlat_grid, &
       &                           rotated_lonlat_grid
  USE mo_grid_structures, ONLY: igrid_icon
  USE mo_grid_structures, ONLY: igrid_cosmo
  USE mo_search_ll_grid, ONLY: find_reg_lonlat_grid_element_index
  USE mo_search_ll_grid, ONLY: find_rotated_lonlat_grid_element_index
  USE mo_logging

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: agg_sgsl_data_to_target_grid

CONTAINS
  !> aggregate DEM slope to target grid
  SUBROUTINE agg_sgsl_data_to_target_grid(sgsl_tiles_grid,       &
       &                                      sgsl_grid,            &
       &                                      tg,                   &
       &                                      sgsl_files,           &
       &                                      sgsl,                 &
       &                                      no_raw_data_pixel,    &
       &                                      raw_data_sgsl_path)

    USE mo_sgsl_data, ONLY : ntiles,   & !< there are 16/240 GLOBE/ASTER tiles
         max_tiles

    USE mo_sgsl_data, ONLY: nc_tot !< number of total GLOBE/ASTER columns un a latitude circle
    USE mo_sgsl_data, ONLY: nr_tot !< total number of rows in GLOBE/ASTER data
    USE mo_sgsl_data, ONLY: get_fill_value   !< determines the _FillValue of either GLOBE or ASTER
    USE mo_sgsl_data, ONLY: idem_type
    USE mo_sgsl_data, ONLY: dem_gl
    USE mo_sgsl_data, ONLY: dem_aster

    USE mo_grid_structures, ONLY: reg_lonlat_grid  !< Definition of Data Type to describe a regular (lonlat) grid
    USE mo_grid_structures, ONLY: rotated_lonlat_grid !< Definition of Data Type to describe a rotated lonlat grid
    USE mo_grid_structures, ONLY: target_grid_def  !< Definition of data type with target grid definition

    USE mo_sgsl_routines, ONLY: open_netcdf_sgsl_tile
    USE mo_sgsl_routines, ONLY: close_netcdf_sgsl_tile
    USE mo_sgsl_routines, ONLY: get_sgsl_data_block
    USE mo_sgsl_routines, ONLY: det_band_gd

    ! USE global data fields (coordinates)
    USE mo_target_grid_data, ONLY: lon_geo, & !< longitude coordinates of the grid in the geographical system
         &                            lat_geo !< latitude coordinates of the grid in the geographical system
    USE mo_target_grid_data, ONLY: search_res ! resolution of ICON grid search index list

    ! USE structure which contains the definition of the COSMO grid
    USE  mo_cosmo_grid, ONLY: COSMO_grid !< structure which contains the definition of the COSMO grid
    ! USE structure which contains the definition of the ICON grid
    USE  mo_icon_grid_data, ONLY: icon_grid !< structure which contains the definition of the ICON grid
    ! USE icon domain structure wich contains the ICON coordinates (and parent-child indices etc)
    USE mo_icon_grid_data, ONLY: icon_grid_region
    ! use additional parameters for height on vertices
    ! as a test the fields are loaded from a module instead of passing in the subroutine call
    USE mo_sgsl_tg_fields, ONLY: add_parameters_domain !< data structure
    USE mo_sgsl_tg_fields, ONLY: vertex_param          !< this structure contains the fields
    !! vertex_param%npixel_vert
    ! USE modules to search in ICON grid
    USE mo_search_icongrid, ONLY:   walk_to_nc, find_nearest_vert

    USE mo_icon_domain,     ONLY: icon_domain

    ! structure for geographica coordintaes
    USE mo_base_geometry,   ONLY: geographical_coordinates
    USE mo_base_geometry,   ONLY: cartesian_coordinates
    USE mo_additional_geometry,   ONLY: cc2gc,                  &
         &                                  gc2cc

    USE mo_math_constants, ONLY: pi, rad2deg, deg2rad
    USE mo_physical_constants, ONLY: re !< av. radius of the earth [m]

    USE mo_bilinterpol, ONLY: get_4_surrounding_raw_data_indices, &
         &                        calc_weight_bilinear_interpol, &
         &                        calc_value_bilinear_interpol

    TYPE(reg_lonlat_grid) :: sgsl_tiles_grid(1:ntiles)        !< raw data grid for the 16/36 GLOBE/ASTER tiles
    TYPE(target_grid_def), INTENT(IN)      :: tg              !< !< structure with target grid description

    TYPE(reg_lonlat_grid) :: sgsl_grid                        !< raw data grid for the whole GLOBE/ASTER dataset
    CHARACTER (LEN=filename_max), INTENT(IN) :: sgsl_files(1:max_tiles)  !< filenames globe/aster raw data

    REAL(KIND=wp), INTENT(OUT)          :: sgsl(1:tg%ie,1:tg%je,1:tg%ke)
    !< mean subgrid-scale slope of target grid element
    INTEGER (KIND=i8), INTENT(OUT)      :: no_raw_data_pixel(1:tg%ie,1:tg%je,1:tg%ke)
    !< number of raw data pixel for a target grid element
    CHARACTER(LEN=filename_max), INTENT(IN), OPTIONAL :: raw_data_sgsl_path !< path to raw data !_br 17.09.14

    ! local variables
    REAL (KIND=wp)    :: lon_sgsl(1:nc_tot)   !< longitude coordinates of the GLOBE grid
    REAL (KIND=wp)    :: lat_sgsl(1:nr_tot)   !< latititude coordinates of the GLOBE grid
    INTEGER (KIND=i4) :: nc_tot_p1
    INTEGER  :: ncids_sgsl(1:ntiles)
    !< ncid for the GLOBE/ASTER tiles, the netcdf files have to be opened by a previous call of open_netcdf_sgsl_tile
    INTEGER (KIND=i8) :: ndata(1:tg%ie,1:tg%je,1:tg%ke)  !< number of raw data pixel with land point

    TYPE(geographical_coordinates) :: target_geo_co  !< structure for geographical coordinates of raw data pixel
    INTEGER :: i,j,k,l ! counters
    INTEGER (KIND=i8) :: ie, je, ke  ! indices for grid elements
    INTEGER (KIND=i8), ALLOCATABLE :: ie_vec(:), iev_vec(:)  ! indices for target grid elements
    INTEGER (KIND=i8) :: i_vert, j_vert, k_vert ! indeces for ICON grid vertices
    INTEGER (KIND=i8) :: i1, i2
    INTEGER :: nv ! counter
    INTEGER :: nt      ! counter
    INTEGER :: mlat, mlat_start ! row number for GLOBE data

    REAL (KIND=wp) :: undef_sgsl
    REAL (KIND=wp) :: default_sgsl
    REAL(KIND=wp)  :: sl_parallel(1:nc_tot)!< one line with GLOBE/ASTER data
    REAL(KIND=wp)  :: sl(0:nc_tot+1) !< slope
    REAL(KIND=wp), ALLOCATABLE :: sl_block(:,:) !< a block of GLOBE/ASTER slope data
    REAL(KIND=wp)  :: row_lat     ! latitude of the row for the slope row
    REAL (KIND=wp) :: bound_north_cosmo !< northern boundary for COSMO target domain
    REAL (KIND=wp) :: bound_south_cosmo !< southern boundary for COSMO target domain
    REAL (KIND=wp) :: bound_west_cosmo  !< western  boundary for COSMO target domain
    REAL (KIND=wp) :: bound_east_cosmo  !< eastern  boundary for COSMO target domain

    ! Some stuff for OpenMP parallelization
    INTEGER :: num_blocks, ib, il, blk_len, istartlon, iendlon, nlon_sub, ishift
    !$ INTEGER :: omp_get_max_threads, omp_get_thread_num, thread_id
    !$ INTEGER (KIND=i8), ALLOCATABLE :: start_cell_arr(:)

    TYPE(reg_lonlat_grid) :: ta_grid
    !< structure with definition of the target area grid (dlon must be the same as for the whole GLOBE/ASTER dataset)
    INTEGER :: block_row_start
    INTEGER :: block_row
    INTEGER :: errorcode !< error status variable
    ! test with walk_to_nc at start
    INTEGER (KIND=i8) :: start_cell_id  !< start cell id
    TYPE(cartesian_coordinates)  :: target_cc_co
    !< coordinates in cartesian system of point for which the nearest ICON grid cell is to be determined
    !variables for GME search
    LOGICAL :: ldebug=.FALSE.
    ! global data flag
    LOGICAL :: gldata=.TRUE. ! GLOBE data are global
    LOGICAL :: lskip
    REAL (KIND=wp) :: point_lon_geo       !< longitude coordinate in geographical system of input point
    REAL (KIND=wp) :: point_lat_geo       !< latitude coordinate in geographical system of input point
    REAL(KIND=wp)   :: point_lon, point_lat
    INTEGER (KIND=i8) :: western_column     !< the index of the western_column of raw data
    INTEGER (KIND=i8) :: eastern_column     !< the index of the eastern_column of raw data
    INTEGER (KIND=i8) :: northern_row       !< the index of the northern_row of raw data
    INTEGER (KIND=i8) :: southern_row       !< the index of the southern_row of raw data
    REAL (KIND=wp) :: sgsl_point_sw
    REAL (KIND=wp) :: sgsl_point_se
    REAL (KIND=wp) :: sgsl_point_ne
    REAL (KIND=wp) :: sgsl_point_nw
    REAL (KIND=wp) :: bwlon !< weight for bilinear interpolation
    REAL (KIND=wp) :: bwlat !< weight for bilinear interpolation
    REAL (KIND=wp) :: sgsl_target_value  !< interpolated altitude from GLOBE data

    CHARACTER(LEN=filename_max) :: sgsl_file_1
    nc_tot_p1 = nc_tot + 1
    sgsl_file_1 = TRIM(raw_data_sgsl_path)//TRIM(sgsl_files(1)) !_br 17.09.14

    !mes <

    SELECT CASE(tg%igrid_type)
    CASE(igrid_icon)  ! ICON GRID
      ke = 1
    CASE(igrid_cosmo)  ! COSMO GRID
      ke = 1
      bound_north_cosmo = MAXVAL(lat_geo) + 0.05_wp  ! add some "buffer"
      bound_north_cosmo = MIN(bound_north_cosmo,90.0_wp)
      bound_south_cosmo = MINVAL(lat_geo) - 0.05_wp  ! add some "buffer"
      bound_south_cosmo = MAX(bound_south_cosmo,-90.0_wp)
      bound_east_cosmo  = MAXVAL(lon_geo) + 0.25_wp  ! add some "buffer"
      bound_east_cosmo  = MIN(bound_east_cosmo,180.0_wp)
      bound_west_cosmo  = MINVAL(lon_geo) - 0.25_wp  ! add some "buffer"
      bound_west_cosmo  = MAX(bound_west_cosmo,-180.0_wp)
    END SELECT

    CALL get_fill_value(sgsl_file_1,undef_sgsl)
    default_sgsl = 0.0

    ! initialize some variables
    no_raw_data_pixel = 0
    ndata      = 0
    sgsl       = 0.0

    IF (tg%igrid_type == igrid_icon) THEN ! Icon grid
      vertex_param%sgsl_vert = 0.0
      vertex_param%npixel_vert = 0
    ENDIF

    ! calculate the longitude coordinate of the GLOBE columns
    DO i=1,nc_tot
      lon_sgsl(i) = sgsl_grid%start_lon_reg + (i-1) * sgsl_grid%dlon_reg
    ENDDO

    ! calculate the latitiude coordinate of the GLOBE columns
    DO j=1,nr_tot
      lat_sgsl(j) = sgsl_grid%start_lat_reg + (j-1) * sgsl_grid%dlat_reg
    ENDDO
    IF (verbose >= idbg_low ) THEN
      WRITE(logging%fileunit,*)'lat_sgsl(1): ', lat_sgsl(1)
      WRITE(logging%fileunit,*)'lat_sgsl(nr_tot) ', lat_sgsl(nr_tot)
    ENDIF

    ALLOCATE(ie_vec(nc_tot),iev_vec(nc_tot))
    ie_vec(:) = 0
    iev_vec(:) = 0
    start_cell_id = 1

    nt = 1

    ! first open the slope netcdf files
    DO nt=1,ntiles
      CALL open_netcdf_sgsl_tile(TRIM(raw_data_sgsl_path)//TRIM(sgsl_files(nt)), ncids_sgsl(nt))
    ENDDO

    mlat_start = 1
    block_row_start = mlat

    CALL det_band_gd(sgsl_grid,block_row_start,ta_grid)
    IF (verbose >= idbg_high ) THEN
      WRITE(logging%fileunit,*)'first call of det_band_gd'
      WRITE(logging%fileunit,*)'ta_grid: ',ta_grid
    ENDIF

    ! Determine start and end longitude of search
    istartlon = 1
    iendlon = nc_tot
    IF (tg%igrid_type == igrid_icon) THEN
      DO i = 1, nc_tot
        point_lon = lon_sgsl(i)
        IF (point_lon < tg%minlon) istartlon = i + 1
        IF (point_lon > tg%maxlon) THEN
          iendlon = i - 1
          EXIT
        ENDIF
      ENDDO
    ELSE IF (tg%igrid_type == igrid_cosmo) THEN
      DO i = 1, nc_tot
        point_lon = lon_sgsl(i)
        IF (point_lon < bound_west_cosmo) istartlon = i + 1
        IF (point_lon > bound_east_cosmo) THEN
          iendlon = i - 1
          EXIT
        ENDIF
      ENDDO
    ENDIF
    nlon_sub = iendlon - istartlon + 1

    num_blocks = 1
    !$ num_blocks = omp_get_max_threads()
    IF (MOD(nlon_sub,num_blocks)== 0) THEN
      blk_len = nlon_sub/num_blocks
    ELSE
      blk_len = nlon_sub/num_blocks + 1
    ENDIF
    IF (verbose >= idbg_low ) WRITE(logging%fileunit,*) 'nlon_sub, num_blocks, blk_len: ',nlon_sub, num_blocks, blk_len

    IF (tg%igrid_type == igrid_cosmo) THEN ! CASE COSMO grid
      DO WHILE (ta_grid%end_lat_reg > bound_north_cosmo)
        block_row_start = block_row_start + ta_grid%nlat_reg
        mlat_start = mlat_start + ta_grid%nlat_reg
        CALL det_band_gd(sgsl_grid,block_row_start, ta_grid)
        PRINT *,block_row_start,ta_grid%end_lat_reg,bound_north_cosmo
      ENDDO
    ELSE IF (tg%igrid_type == igrid_icon) THEN
      DO WHILE (ta_grid%end_lat_reg > tg%maxlat)
        block_row_start = block_row_start + ta_grid%nlat_reg
        mlat_start = mlat_start + ta_grid%nlat_reg
        CALL det_band_gd(sgsl_grid,block_row_start, ta_grid)
      ENDDO
    ENDIF ! grid type

    IF(ALLOCATED(sl_block)) THEN
      DEALLOCATE(sl_block, STAT=errorcode)
      IF(errorcode/=0) CALL abort_extpar('Cant deallocate the sl_block')
    ENDIF
    ALLOCATE (sl_block(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg), STAT=errorcode)
    IF(errorcode/=0) CALL abort_extpar('Cant allocate sl_block')

    CALL get_sgsl_data_block(sgsl_file_1,       &
         &                       ta_grid,         &
         &                       sgsl_tiles_grid, &
         &                       ncids_sgsl,      &
         &                       sl_block)

    block_row = 0

    WRITE(logging%fileunit,*)'INFO: In routine agg_sgsl_data_to_target_grid:'
    WRITE(logging%fileunit,*)'                    Start loop over slope rows'
    !-----------------------------------------------------------------------------
    slope_rows: DO mlat=mlat_start,nr_tot
      !-----------------------------------------------------------------------------
      !-----------------------------------------------------------------------------
      IF (verbose >= idbg_high ) THEN
        IF (MOD(mlat,100)==0) WRITE(logging%fileunit,*) 'slope row:', mlat
      ENDIF
      block_row= block_row + 1
      IF((block_row > ta_grid%nlat_reg).AND.(mlat<nr_tot)) THEN ! read in new block
        block_row_start = mlat + 1
        block_row = 1
        CALL det_band_gd(sgsl_grid,block_row_start, ta_grid)
        IF (verbose >= idbg_high ) THEN
          WRITE(logging%fileunit,*)'next call of det_band_gd'
          WRITE(logging%fileunit,*)'ta_grid: ',ta_grid
        ENDIF

        IF (tg%igrid_type == igrid_cosmo) THEN ! CASE COSMO grid
          IF (ta_grid%start_lat_reg < bound_south_cosmo) THEN
            EXIT slope_rows
          ENDIF
        ELSE IF (tg%igrid_type == igrid_icon) THEN
          IF (ta_grid%start_lat_reg < tg%minlat) THEN
            EXIT slope_rows
          ENDIF
        ENDIF ! grid type

        IF(ALLOCATED(sl_block)) THEN
          DEALLOCATE(sl_block, STAT=errorcode)
          IF(errorcode/=0) CALL abort_extpar('Cant deallocate the sl_block')
        ENDIF
        ALLOCATE (sl_block(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg), STAT=errorcode)
        IF(errorcode/=0) CALL abort_extpar('Cant allocate sl_block')
        CALL get_sgsl_data_block(sgsl_file_1,     &
             &                       ta_grid,         &
             &                       sgsl_tiles_grid, &
             &                       ncids_sgsl,     &
             &                       sl_block)
      ENDIF

      sl_parallel(1:nc_tot) = sl_block(1:nc_tot,block_row)
      sl(1:nc_tot) = sl_parallel(1:nc_tot)  ! put data to "central row"
      sl(0)        = sl_parallel(nc_tot) ! western wrap at -180/180 degree longitude
      sl(nc_tot_p1) = sl_parallel(1)      ! eastern wrap at -180/180 degree longi
      row_lat = sgsl_grid%start_lat_reg + (mlat-1) * sgsl_grid%dlat_reg
      lskip = .FALSE.
      IF (tg%igrid_type == igrid_cosmo) THEN ! CASE COSMO grid
        IF ((row_lat > bound_north_cosmo).OR.(row_lat < bound_south_cosmo) ) THEN ! raw data out of target grid
          lskip = .TRUE.
        ENDIF
      ELSE IF (tg%igrid_type == igrid_icon) THEN
        IF (row_lat > tg%maxlat .OR. row_lat < tg%minlat) lskip = .TRUE.
      ENDIF ! grid type

      IF (lskip) THEN
        CYCLE slope_rows
      ENDIF

      ! set undefined values to 0 altitude (default)
      WHERE (sl == undef_sgsl)
        sl = default_sgsl
      END WHERE
      WHERE (sl_parallel == undef_sgsl)
        sl_parallel = default_sgsl
      END WHERE

      IF (tg%igrid_type == igrid_icon) THEN
        ie_vec(istartlon:iendlon) = 0
        iev_vec(istartlon:iendlon) = 0
      ENDIF

      point_lat = row_lat

      IF (tg%igrid_type == igrid_icon) THEN
!$OMP PARALLEL DO PRIVATE(ib,il,i,i1,i2,ishift,point_lon,thread_id,start_cell_id,target_geo_co,target_cc_co)
        DO ib = 1, num_blocks

          !$   thread_id = omp_get_thread_num()+1
          !$   start_cell_id = start_cell_arr(thread_id)
          ishift = istartlon-1+(ib-1)*blk_len

          ! loop over one latitude circle of the raw data
          columns1: DO il = 1,blk_len
            i = ishift+il
            IF (i > iendlon) CYCLE columns1

            ! find the corresponding target grid indices
            point_lon = lon_sgsl(i)

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
            CALL walk_to_nc(icon_grid_region,   &
                 target_cc_co,     &
                 start_cell_id,    &
                 icon_grid%nvertex_per_cell, &
                 icon_grid%nedges_per_vertex, &
                 ie_vec(i))

            ! additional get the nearest vertex index for accumulating height values there
            IF (ie_vec(i) /= 0_i8) THEN
              CALL  find_nearest_vert(icon_grid_region, &
                   target_cc_co,                  &
                   ie_vec(i),        &
                   icon_grid%nvertex_per_cell,    &
                   iev_vec(i))
            ENDIF

          ENDDO columns1
          !$   start_cell_arr(thread_id) = start_cell_id
        ENDDO ! ib
!$OMP END PARALLEL DO
      ENDIF ! ICON only

!$OMP PARALLEL DO PRIVATE(i,ie,je,ke,i_vert,point_lon)
      DO i=istartlon,iendlon

        ! call here the attribution of raw data pixel to target grid for different grid types
        SELECT CASE(tg%igrid_type)
        CASE(igrid_icon)  ! ICON GRID

          ie = ie_vec(i)
          je = 1
          ke = 1

          ! get the nearest vertex index for accumulating height values there

          ! aggregate the vertex parameter here
          i_vert = iev_vec(i)
          j_vert = 1
          k_vert = 1
          IF ((i_vert /=0)) THEN ! raw data pixel within target grid
            vertex_param%npixel_vert(i_vert,j_vert,k_vert) =  &
                 vertex_param%npixel_vert(i_vert,j_vert,k_vert) + 1

            vertex_param%sgsl_vert(i_vert,j_vert,k_vert) =  &
                 vertex_param%sgsl_vert(i_vert,j_vert,k_vert) +  sl_parallel(i)
          ENDIF

        CASE(igrid_cosmo)  ! COSMO GRID

          point_lon = lon_sgsl(i)

          CALL find_rotated_lonlat_grid_element_index(point_lon,  &
               point_lat,  &
               COSMO_grid, &
               ie,         &
               je)
          ke = 1

          IF ((ie /= 0).AND.(je/=0).AND.(ke/=0))THEN
!$OMP CRITICAL
            ! raw data pixel within target grid, see output of routine find_rotated_lonlat_grid_element_index
            no_raw_data_pixel(ie,je,ke) = no_raw_data_pixel(ie,je,ke) + 1
            !- summation of variables
            SELECT CASE(idem_type)
            CASE(dem_aster)
              IF (sl(i) /= default_sgsl) THEN
                ndata(ie,je,ke)      = ndata(ie,je,ke) + 1
                sgsl(ie,je,ke)  = sgsl(ie,je,ke) + sl(i)
              ENDIF
            CASE(dem_gl)
              IF (sl(i) /= undef_sgsl) THEN
                ndata(ie,je,ke)      = ndata(ie,je,ke) + 1
                sgsl(ie,je,ke)  = sgsl(ie,je,ke) + sl(i)
              ENDIF
            END SELECT
!$OMP END CRITICAL
          ENDIF

        END SELECT
      ENDDO ! loop over one latitude circle of the raw data
!$OMP END PARALLEL DO

      !-----------------------------------------------------------------------------
      !-----------------------------------------------------------------------------
    ENDDO slope_rows
    !-----------------------------------------------------------------------------

    DEALLOCATE(ie_vec,iev_vec)

    WRITE(logging%fileunit,*)'INFO: In routine agg_sgsl_data_to_target_grid:'
    WRITE(logging%fileunit,*)'                      loop over topo rows done'

    IF (verbose >= idbg_low ) THEN
      WRITE(logging%fileunit,*)'Maximum number of slope raw data pixel in a target grid element: '
      WRITE(logging%fileunit,*)'MAXVAL(no_raw_data_pixel): ', MAXVAL(no_raw_data_pixel)
      WRITE(logging%fileunit,*)'Index of target grid element: ', MAXLOC(no_raw_data_pixel)
      WRITE(logging%fileunit,*)'Maximum number of slope land pixel in a target grid element: '
      WRITE(logging%fileunit,*)'MAXVAL(ndata): ', MAXVAL(ndata)
      WRITE(logging%fileunit,*)'Index of target grid element: ', MAXLOC(ndata)
      WRITE(logging%fileunit,*)'Minimal number of slope raw data pixel in a target grid element: '
      WRITE(logging%fileunit,*)'MINVAL(no_raw_data_pixel): ', MINVAL(no_raw_data_pixel)
      WRITE(logging%fileunit,*)'Index of target grid element: ', MINLOC(no_raw_data_pixel)
      WRITE(logging%fileunit,*)'Minimal number of slope land pixel in a target grid element: '
      WRITE(logging%fileunit,*)'MINVAL(ndata): ', MINVAL(ndata)
      WRITE(logging%fileunit,*)'Index of target grid element: ', MINLOC(ndata)
    ENDIF




    IF (verbose >= idbg_low ) WRITE(logging%fileunit,*)'Average slope'
    ! Average height
    DO ke=1, tg%ke
      DO je=1, tg%je
        DO ie=1, tg%ie
          IF (no_raw_data_pixel(ie,je,ke) /= 0) THEN
            ! avoid division by zero for small target grids
            sgsl(ie,je,ke) = sgsl(ie,je,ke)/no_raw_data_pixel(ie,je,ke)
            ! average slope, oceans point counted as 0
          ELSE
            sgsl(ie,je,ke) = REAL(default_sgsl)
          ENDIF
        ENDDO
      ENDDO
    ENDDO


    IF (tg%igrid_type == igrid_icon) THEN ! CASE ICON grid
      IF (verbose >= idbg_low ) WRITE(logging%fileunit,*)'Average slope for vertices'
      ! Average slope for vertices
      DO ke=1, 1
        DO je=1, 1
          DO ie=1, icon_grid_region%nverts
            IF (vertex_param%npixel_vert(ie,je,ke) /= 0) THEN
              ! avoid division by zero for small target grids
              vertex_param%sgsl_vert(ie,je,ke) =  &
                   vertex_param%sgsl_vert(ie,je,ke)/vertex_param%npixel_vert(ie,je,ke)              ! average slope
            ELSE
              vertex_param%sgsl_vert(ie,je,ke) = REAL(default_sgsl)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    DO ke=1, tg%ke
      DO je=1, tg%je
        DO ie=1, tg%ie
          IF (no_raw_data_pixel(ie,je,ke) == 0) THEN
            ! bilinear interpolation to target grid

            point_lon_geo = lon_geo(ie,je,ke)
            point_lat_geo = lat_geo(ie,je,ke)

            CALL bilinear_interpol_sgsl_to_target_point(raw_data_sgsl_path, &
                 &                                      sgsl_files, &
                 &                                      sgsl_grid,       &
                 &                                      sgsl_tiles_grid, &
                 &                                      ncids_sgsl,     &
                 &                                      lon_sgsl,       &
                 &                                      lat_sgsl,       &
                 &                                      point_lon_geo,   &
                 &                                      point_lat_geo,   &
                 &                                      sgsl_target_value)

            sgsl(ie,je,ke) = sgsl_target_value

          ENDIF
        ENDDO
      ENDDO
    ENDDO

    SELECT CASE(tg%igrid_type)
    CASE(igrid_icon) ! ICON GRID
      je=1
      ke=1
      DO nv=1, icon_grid_region%nverts
        IF (vertex_param%npixel_vert(nv,je,ke) == 0) THEN
          ! interpolate from raw data in this case
          point_lon_geo =  rad2deg * icon_grid_region%verts%vertex(nv)%lon
          point_lat_geo =  rad2deg * icon_grid_region%verts%vertex(nv)%lat

          CALL bilinear_interpol_sgsl_to_target_point(raw_data_sgsl_path, &
               &                                      sgsl_files, &
               &                                      sgsl_grid,       &
               &                                      sgsl_tiles_grid, &
               &                                      ncids_sgsl,     &
               &                                      lon_sgsl,       &
               &                                      lat_sgsl,       &
               &                                      point_lon_geo,   &
               &                                      point_lat_geo,   &
               &                                      sgsl_target_value)

          vertex_param%sgsl_vert(nv,je,ke) = sgsl_target_value
        ENDIF
      ENDDO
    END SELECT
    ! close the raw data netcdf files
    DO nt=1,ntiles
      CALL close_netcdf_sgsl_tile(ncids_sgsl(nt))
    ENDDO
    WRITE(logging%fileunit,*)'INFO: routine agg_sgsl_data_to_target_grid done'
  END SUBROUTINE agg_sgsl_data_to_target_grid

  !----------------------------------------------------------------------------------------------------------------

  !> subroutine for bilenar interpolation from GLOBE data (regular lonlat grid) to a single target point
  !!
  !! the slope data are passed to the subroutine in the sgsl_data_block 2D-Array, which is re-read
  !! from the raw data file if the target point is out of the range of the data block.
  !! (If the data block is not too small, repeated I/O to the hard disk is avoided, reading from memory is much faster.)
  !!
  !! the definition of the regular lon-lat grid requires
  !! - the coordinates of the north-western point of the domain ("upper left") startlon_reg_lonlat and startlat_reg_lonlat
  !! - the increment dlon_reg_lonlat and dlat_reg_lonlat(implict assuming that the grid definiton goes
  !!   from the west to the east and from the north to the south)
  !! - the number of grid elements nlon_reg_lonlat and nlat_reg_lonlat for both directions
  SUBROUTINE bilinear_interpol_sgsl_to_target_point(raw_data_sgsl_path, &
       sgsl_files,  &
       sgsl_grid,      &
       sgsl_tiles_grid,&
       ncids_sgsl,    &
       lon_sgsl,      &
       lat_sgsl,      &
       point_lon_geo,  &
       point_lat_geo,  &
       sgsl_target_value)

    USE mo_sgsl_data, ONLY: ntiles !< there are 16/36 GLOBE/ASTER tiles
    USE mo_sgsl_data, ONLY: nc_tot !< number of total DEM columns at a latitude circle
    USE mo_sgsl_data, ONLY: nr_tot !< number of total DEM rows at a latitude circle

    USE mo_sgsl_data, ONLY: get_fill_value  ! determines the corresponding _FillValue of slope data
    USE mo_sgsl_data, ONLY: max_tiles

    USE mo_grid_structures, ONLY: reg_lonlat_grid  !< Definition of Data Type to describe a regular (lonlat) grid
    USE mo_grid_structures, ONLY: rotated_lonlat_grid !< Definition of Data Type to describe a rotated lonlat grid
    USE mo_sgsl_routines, ONLY: open_netcdf_sgsl_tile
    USE mo_sgsl_routines, ONLY: close_netcdf_sgsl_tile
    USE mo_sgsl_routines, ONLY: get_sgsl_data_block
    USE mo_bilinterpol, ONLY: get_4_surrounding_raw_data_indices, &
         &                        calc_weight_bilinear_interpol, &
         &                        calc_value_bilinear_interpol
    CHARACTER(len=filename_max), INTENT(IN)     :: sgsl_files(1:max_tiles)
    TYPE(reg_lonlat_grid), INTENT(IN) :: sgsl_grid
    !< structure with definition of the raw data grid for the whole slope dataset
    TYPE(reg_lonlat_grid), INTENT(IN) :: sgsl_tiles_grid(1:ntiles) !< raw data grid for the 16/36 GLOBE/ASTER tiles
    !< ncid for the sgsl tiles, the netcdf files have to be opened by a previous call of open_netcdf_GLOBE_tile
    INTEGER (KIND=i4), INTENT(IN)     :: ncids_sgsl(1:ntiles)

    REAL (KIND=wp), INTENT(IN) :: lon_sgsl(1:nc_tot)   !< longitude coordinates of the GLOBE grid
    REAL (KIND=wp), INTENT(IN) :: lat_sgsl(1:nr_tot)   !< latititude coordinates of the GLOBE grid
    REAL (KIND=wp), INTENT(IN) :: point_lon_geo       !< longitude coordinate in geographical system of input point
    REAL (KIND=wp), INTENT(IN) :: point_lat_geo       !< latitude coordinate in geographical system of input point
    REAL (KIND=wp), INTENT(OUT) :: sgsl_target_value  !< interpolated slope value from DEM

    ! local variables
    REAL (KIND=wp), ALLOCATABLE :: sl_block(:,:) !< a block of slope data
    TYPE(reg_lonlat_grid) :: ta_grid
    !< structure with definition of the target area grid (dlon must be the same as for the whole GLOBE dataset)
    INTEGER :: nt      ! counter
    INTEGER  (KIND=i8) :: point_lon_index !< longitude index of point for regular lon-lat grid
    INTEGER  (KIND=i8) :: point_lat_index !< latitude index of point for regular lon-lat grid
    INTEGER (KIND=i8) :: western_column     !< the index of the western_column of data to read in
    INTEGER (KIND=i8) :: eastern_column     !< the index of the eastern_column of data to read in
    INTEGER (KIND=i8) :: northern_row       !< the index of the northern_row of data to read in
    INTEGER (KIND=i8) :: southern_row       !< the index of the southern_row of data to read in
    REAL (KIND=wp)   :: bwlon  !< weight for bilinear interpolation
    REAL (KIND=wp)   :: bwlat  !< weight for bilinear interpolation
    REAL (KIND=wp) :: south_lat !< southern latitude of DEM data pixel for bilinear interpolation
    REAL (KIND=wp) :: west_lon  !< western longitude of DEM data pixel for bilinear interpolation
    REAL (KIND=wp) :: pixel_lon !< longitude coordinate in geographical system of input point
    REAL (KIND=wp) :: pixel_lat !< latitude coordinate in geographical system of input point
    REAL (KIND=wp) :: sgsl_pixel_lon !< longitude coordinate in geographical system of globe raw data point
    REAL (KIND=wp) :: sgsl_pixel_lat !< latitude coordinate in geographical system of globe raw data point
    REAL (KIND=wp)   :: sgsl_point_sw       !< value of the DEM raw data pixel south west
    REAL (KIND=wp)   :: sgsl_point_se       !< value of the DEM raw data pixel south east
    REAL (KIND=wp)   :: sgsl_point_ne       !< value of the DEM raw data pixel north east
    REAL (KIND=wp)   :: sgsl_point_nw       !< value of the DEM raw data pixel north west
    INTEGER :: errorcode
    LOGICAL :: gldata=.TRUE. ! DEM data are global
    INTEGER :: ndata
    INTEGER :: nland
    REAL (KIND=wp) :: undef_sgsl
    REAL (KIND=wp) :: default_sgsl
    CHARACTER(len=filename_max) :: sgsl_file_1
    CHARACTER(len=filename_max) :: raw_data_sgsl_path

    sgsl_file_1 = TRIM(raw_data_sgsl_path)//TRIM(sgsl_files(1))

    CALL get_fill_value(sgsl_file_1,undef_sgsl)

    default_sgsl = 0

    ! get four surrounding raw data indices
    CALL  get_4_surrounding_raw_data_indices(sgsl_grid,     &
         &                                      lon_sgsl,     &
         &                                      lat_sgsl,     &
         &                                      gldata,        &
         &                                      point_lon_geo, &
         &                                      point_lat_geo, &
         &                                      western_column,&
         &                                      eastern_column,&
         &                                      northern_row,  &
         &                                      southern_row)

    ta_grid%dlon_reg = sgsl_grid%dlon_reg
    ta_grid%dlat_reg = sgsl_grid%dlat_reg
    ta_grid%nlon_reg = eastern_column - western_column + 1
    ta_grid%nlat_reg = southern_row - northern_row + 1
    ta_grid%start_lon_reg = lon_sgsl(western_column)
    ta_grid%end_lon_reg = lon_sgsl(eastern_column)
    ta_grid%start_lat_reg = lat_sgsl(northern_row)
    ta_grid%end_lat_reg  = lat_sgsl(southern_row)

    ! calculate weight for bilinear interpolation
    CALL calc_weight_bilinear_interpol(point_lon_geo,            &
         &                                point_lat_geo,            &
         &                                lon_sgsl(western_column), &
         &                                lon_sgsl(eastern_column), &
         &                                lat_sgsl(northern_row),   &
         &                                lat_sgsl(southern_row),   &
         &                                bwlon,                    &
         &                                bwlat)

    ALLOCATE (sl_block(western_column:eastern_column,northern_row:southern_row), STAT=errorcode)
    IF(errorcode/=0) CALL abort_extpar('Cant allocate h_block')
    CALL get_sgsl_data_block(sgsl_file_1,     &
         &                       ta_grid,         &
         &                       sgsl_tiles_grid, &
         &                       ncids_sgsl,     &
         &                       sl_block)
    ! check for undefined GLOBE data, which indicate ocean grid element

    sgsl_point_sw = sl_block(western_column,southern_row)
    sgsl_point_se = sl_block(eastern_column,southern_row)
    sgsl_point_ne = sl_block(eastern_column,northern_row)
    sgsl_point_nw = sl_block(western_column,northern_row)

    ! perform the interpolation for height
    sgsl_target_value = calc_value_bilinear_interpol(bwlon,         &
         &       bwlat,         &
         &       sgsl_point_sw, &
         &       sgsl_point_se, &
         &       sgsl_point_ne, &
         &       sgsl_point_nw)

  END SUBROUTINE bilinear_interpol_sgsl_to_target_point

  !----------------------------------------------------------------------------------------------------------------

END MODULE mo_agg_sgsl
