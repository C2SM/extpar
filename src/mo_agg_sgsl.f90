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

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4
  USE mo_io_units,              ONLY: filename_max
                                
  USE mo_grid_structures,       ONLY: igrid_icon, &
       &                              igrid_cosmo, &
       &                              reg_lonlat_grid, &  !< Definition of Data Type to describe a regular (lonlat) grid
       &                              target_grid_def  !< Definition of data type with target grid definition
                                
  USE mo_search_ll_grid,        ONLY: find_rotated_lonlat_grid_element_index 
                                
  USE mo_sgsl_routines,          ONLY: open_netcdf_sgsl_tile, &
       &                               close_netcdf_sgsl_tile, &
       &                               get_sgsl_data_block, &
       &                               det_band_gd
                                
  USE mo_target_grid_data,       ONLY: lon_geo, & !< longitude coordinates of the grid in the geographical system
       &                               lat_geo, & !< latitude coordinates of the grid in the geographical system
       &                               search_res ! resolution of ICON grid search index list
                                
  USE mo_cosmo_grid,             ONLY: COSMO_grid !< structure which contains the definition of the COSMO grid
                                
  USE mo_icon_grid_data,         ONLY: icon_grid, & !< structure which contains the definition of the ICON grid
       &                               icon_grid_region 
                                
  USE mo_search_icongrid,        ONLY: walk_to_nc
                                
  USE mo_base_geometry,          ONLY: geographical_coordinates
  USE mo_base_geometry,          ONLY: cartesian_coordinates
  USE mo_additional_geometry,    ONLY: gc2cc

  USE mo_math_constants,         ONLY: deg2rad

  USE mo_bilinterpol,            ONLY: get_4_surrounding_raw_data_indices, &
       &                               calc_weight_bilinear_interpol, &
       &                               calc_value_bilinear_interpol

  USE mo_topo_data,              ONLY: ntiles, &
       &                               nc_tot, & !< number of total GLOBE/ASTER columns un a latitude circle
       &                               max_tiles, &
       &                               itopo_type, &
       &                               topo_gl, &
       &                               topo_aster, &
       &                               topo_merit, &
       &                               topo_copernicus, &
       &                               get_fill_value_sgsl, &
       &                               nr_tot !< total number of rows in GLOBE/ASTER data

   USE mo_io_utilities,         ONLY: join_path

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


    TYPE(target_grid_def), INTENT(IN)                 :: tg              !< !< structure with target grid description
                                                      
    CHARACTER (LEN=*), INTENT(IN)                     :: sgsl_files(1:max_tiles), &  !< filenames globe/aster raw data
         &                                               raw_data_sgsl_path !< path to raw data !_br 17.09.14

    REAL(KIND=wp), INTENT(OUT)                        :: sgsl(1:tg%ie,1:tg%je,1:tg%ke)

    INTEGER (KIND=i4), INTENT(OUT)                    :: no_raw_data_pixel(1:tg%ie,1:tg%je,1:tg%ke)

    ! local variables
    TYPE(reg_lonlat_grid)                             :: sgsl_tiles_grid(1:ntiles), &        !< raw data grid for the 16/36 GLOBE/ASTER tiles
         &                                               sgsl_grid, &                        !< raw data grid for the whole GLOBE/ASTER dataset
         &                                               ta_grid

    TYPE(geographical_coordinates)                    :: target_geo_co  !< structure for geographical coordinates of raw data pixel
    TYPE(cartesian_coordinates)                       :: target_cc_co

    INTEGER (KIND=i4)                                 :: ie, je, ke, &  ! indices for grid elements
         &                                               nc_tot_p1, &
         &                                               ndata(1:tg%ie,1:tg%je,1:tg%ke), &  !< number of raw data pixel with land point
         &                                               i1, i2, &
         &                                               start_cell_id, &  !< start cell id
         &                                               ncids_sgsl(1:ntiles), &
         &                                               i,j, & ! counters
         &                                               nt, &      ! counter
         &                                               mlat, mlat_start, & ! row number for GLOBE data
         &                                               block_row_start, &
         &                                               block_row, &
         &                                               errorcode !< error status variable

    INTEGER (KIND=i4), ALLOCATABLE                    :: ie_vec(:), iev_vec(:)  ! indices for target grid elements

    REAL (KIND=wp)                                    :: lon_sgsl(1:nc_tot), &   !< longitude coordinates of the GLOBE grid
         &                                               lat_sgsl(1:nr_tot), &   !< latititude coordinates of the GLOBE grid
         &                                               undef_sgsl, &
         &                                               default_sgsl, &
         &                                               sl_parallel(1:nc_tot), &!< one line with GLOBE/ASTER data
         &                                               sl(0:nc_tot+1), & !< slope
         &                                               row_lat, &     ! latitude of the row for the slope row
         &                                               bound_north_cosmo, & !< northern boundary for COSMO target domain
         &                                               bound_south_cosmo, & !< southern boundary for COSMO target domain
         &                                               bound_west_cosmo, &  !< western  boundary for COSMO target domain
         &                                               bound_east_cosmo, &  !< eastern  boundary for COSMO target domain
         &                                               point_lon_geo, &       !< longitude coordinate in geographical system of input point
         &                                               point_lat_geo, &       !< latitude coordinate in geographical system of input point
         &                                               point_lon, point_lat, &
         &                                               sgsl_target_value  !< interpolated altitude from GLOBE data

    REAL(KIND=wp), ALLOCATABLE                        :: sl_block(:,:) !< a block of GLOBE/ASTER slope data

    ! Some stuff for OpenMP parallelization
    INTEGER(KIND=i4)                                  :: num_blocks, ib, il, blk_len, istartlon, iendlon, nlon_sub, ishift
    !$ INTEGER :: omp_get_max_threads, omp_get_thread_num, thread_id
    !$ INTEGER (KIND=i4), ALLOCATABLE :: start_cell_arr(:)

    ! global data flag
    LOGICAL                                           :: lskip
    CHARACTER(LEN=filename_max)                       :: sgsl_file_1

    CALL logging%info('Enter routine: agg_sgsl_data_to_target_grid')

    nc_tot_p1 = nc_tot + 1
    sgsl_file_1 = join_path(raw_data_sgsl_path,sgsl_files(1)) !_br 17.09.14

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

    CALL get_fill_value_sgsl(sgsl_file_1,undef_sgsl)
    default_sgsl = 0.0

    ! initialize some variables
    no_raw_data_pixel = 0
    ndata      = 0
    sgsl       = 0.0

    ! calculate the longitude coordinate of the GLOBE columns
    DO i=1,nc_tot
      lon_sgsl(i) = sgsl_grid%start_lon_reg + (i-1) * sgsl_grid%dlon_reg
    ENDDO

    ! calculate the latitiude coordinate of the GLOBE columns
    DO j=1,nr_tot
      lat_sgsl(j) = sgsl_grid%start_lat_reg + (j-1) * sgsl_grid%dlat_reg
    ENDDO

    ALLOCATE(ie_vec(nc_tot),iev_vec(nc_tot))
    ie_vec(:) = 0
    iev_vec(:) = 0
    start_cell_id = 1

    nt = 1
    ! first open the slope netcdf files
    DO nt=1,ntiles
      CALL open_netcdf_sgsl_tile(join_path(raw_data_sgsl_path,sgsl_files(nt)), ncids_sgsl(nt))

    ENDDO

    mlat_start = 1
    block_row_start = mlat

    CALL det_band_gd(sgsl_grid,block_row_start,ta_grid)

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

    IF (tg%igrid_type == igrid_cosmo) THEN ! CASE COSMO grid
      DO WHILE (ta_grid%end_lat_reg > bound_north_cosmo)
        block_row_start = block_row_start + ta_grid%nlat_reg
        mlat_start = mlat_start + ta_grid%nlat_reg
        CALL det_band_gd(sgsl_grid,block_row_start, ta_grid)
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
      IF(errorcode/=0) CALL logging%error('Cant deallocate the sl_block',__FILE__,__LINE__)
    ENDIF
    ALLOCATE (sl_block(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg), STAT=errorcode)
    IF(errorcode/=0) CALL logging%error('Cant allocate sl_block',__FILE__,__LINE__)

    CALL get_sgsl_data_block(sgsl_file_1,       &
         &                       ta_grid,         &
         &                       sgsl_tiles_grid, &
         &                       ncids_sgsl,      &
         &                       sl_block)

    block_row = 0

    WRITE(message_text,*) 'nlon_sub: ',nlon_sub,' num_blocks: ',num_blocks, ' blk_len: ',blk_len
    CALL logging%info(message_text)
    
    CALL logging%info('Start loop over slope rows...')
    !-----------------------------------------------------------------------------
    slope_rows: DO mlat=mlat_start,nr_tot
      !-----------------------------------------------------------------------------
      !-----------------------------------------------------------------------------
      block_row= block_row + 1
      IF((block_row > ta_grid%nlat_reg).AND.(mlat<nr_tot)) THEN ! read in new block
        block_row_start = mlat + 1
        block_row = 1
        CALL det_band_gd(sgsl_grid,block_row_start, ta_grid)

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
          IF(errorcode/=0) CALL logging%error('Cant deallocate the sl_block',__FILE__,__LINE__)
        ENDIF
        ALLOCATE (sl_block(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg), STAT=errorcode)
        IF(errorcode/=0) CALL logging%error('Cant allocate sl_block',__FILE__,__LINE__)
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

          ENDDO columns1
          !$   start_cell_arr(thread_id) = start_cell_id
        ENDDO ! ib
!$OMP END PARALLEL DO
      ENDIF ! ICON only

!$OMP PARALLEL DO PRIVATE(i,ie,je,ke,point_lon)
      DO i=istartlon,iendlon

        ! call here the attribution of raw data pixel to target grid for different grid types
        SELECT CASE(tg%igrid_type)
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
              SELECT CASE(itopo_type)
                CASE(topo_aster)
                  IF (sl(i) /= default_sgsl) THEN
                    ndata(ie,je,ke)      = ndata(ie,je,ke) + 1
                    sgsl(ie,je,ke)  = sgsl(ie,je,ke) + sl(i)
                  ENDIF
                CASE(topo_gl)
                  IF (sl(i) /= undef_sgsl) THEN
                    ndata(ie,je,ke)      = ndata(ie,je,ke) + 1
                    sgsl(ie,je,ke)  = sgsl(ie,je,ke) + sl(i)
                  ENDIF
                CASE(topo_merit)
                  IF (sl(i) /= undef_sgsl) THEN
                    ndata(ie,je,ke)      = ndata(ie,je,ke) + 1
                    sgsl(ie,je,ke)  = sgsl(ie,je,ke) + sl(i)
                  ENDIF
                CASE(topo_copernicus)
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

    CALL logging%info('...done')

    DEALLOCATE(ie_vec,iev_vec)

    CALL logging%info('Compute Average slope')
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

    ! close the raw data netcdf files
    DO nt=1,ntiles
      CALL close_netcdf_sgsl_tile(ncids_sgsl(nt))
    ENDDO
    
  CALL logging%info('Exit routine: agg_sgsl_data_to_target_grid')

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

    CHARACTER(len=*), INTENT(IN)            :: sgsl_files(1:max_tiles)

    TYPE(reg_lonlat_grid), INTENT(IN)       :: sgsl_grid, &
         &                                     sgsl_tiles_grid(1:ntiles) !< raw data grid for the 16/36 GLOBE/ASTER tiles

    INTEGER (KIND=i4), INTENT(IN)           :: ncids_sgsl(1:ntiles)

    REAL (KIND=wp), INTENT(IN)              :: lon_sgsl(1:nc_tot), &   !< longitude coordinates of the GLOBE grid
         &                                     lat_sgsl(1:nr_tot), &   !< latititude coordinates of the GLOBE grid
         &                                     point_lon_geo, &       !< longitude coordinate in geographical system of input point
         &                                     point_lat_geo       !< latitude coordinate in geographical system of input point

    REAL (KIND=wp), INTENT(OUT)             :: sgsl_target_value  !< interpolated slope value from DEM
                                            
    ! local variables                       
    TYPE(reg_lonlat_grid)                   :: ta_grid
    REAL (KIND=wp), ALLOCATABLE             :: sl_block(:,:) !< a block of slope data
                                            
    INTEGER (KIND=i4)                       :: western_column, &     !< the index of the western_column of data to read in
         &                                     eastern_column, &     !< the index of the eastern_column of data to read in
         &                                     northern_row, &       !< the index of the northern_row of data to read in
         &                                     errorcode, &
         &                                     southern_row       !< the index of the southern_row of data to read in
                                            
    REAL (KIND=wp)                          :: bwlon, &  !< weight for bilinear interpolation
         &                                     bwlat, &  !< weight for bilinear interpolation
         &                                     sgsl_point_sw, &       !< value of the DEM raw data pixel south west
         &                                     sgsl_point_se, &       !< value of the DEM raw data pixel south east
         &                                     sgsl_point_ne, &       !< value of the DEM raw data pixel north east
         &                                     sgsl_point_nw, &       !< value of the DEM raw data pixel north west
         &                                     undef_sgsl

    LOGICAL                                 :: gldata=.TRUE. ! DEM data are global
                                            
    CHARACTER(len=filename_max)             :: sgsl_file_1, &
         &                                     raw_data_sgsl_path

    sgsl_file_1 = join_path(raw_data_sgsl_path,sgsl_files(1))

    CALL get_fill_value_sgsl(sgsl_file_1,undef_sgsl)

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
    IF(errorcode/=0) CALL logging%error('Cant allocate h_block',__FILE__,__LINE__)
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

END MODULE mo_agg_sgsl
