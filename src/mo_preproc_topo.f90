MODULE mo_preproc_topo

  USE mo_logging
  USE mo_kind,                  ONLY: wp,i4

  USE mo_io_units,              ONLY: filename_max
  USE mo_io_utilities,          ONLY: check_netcdf

  USE mo_grid_structures,       ONLY: reg_lonlat_grid, &
       &                              rotated_lonlat_grid, &
       &                              icosahedral_triangular_grid, &
       &                              target_grid_def, &
       &                              igrid_icon

  USE mo_icon_domain,           ONLY: icon_domain, grid_cells

  USE mo_topo_data,             ONLY: ntiles,         & !< there are 16/240 GLOBE/ASTER tiles
       &                              max_tiles,      &
       &                              nr_tot,         & !< total number of rows in GLOBE/ASTER data
       &                              get_fill_value, &
       &                              itopo_type,     &
       &                              topo_gl,        &
       &                              topo_aster,     &
       &                              get_varname

  USE mo_topo_routines,         ONLY: open_netcdf_TOPO_tile, &
       &                              close_netcdf_TOPO_tile, &
       &                              get_topo_data_block, &
       &                              det_band_gd

  USE mo_target_grid_data,      ONLY: lon_geo, & !< longitude coordinates of the grid in the geographical system
       &                              lat_geo, &    !< latitude coordinates of the grid in the geographical system
       &                              search_res ! resolution of ICON grid search index list

  USE mo_icon_grid_data,        ONLY: icon_grid, & !< structure which contains the definition of the ICON grid
       &                              icon_grid_region

  USE mo_topo_tg_fields,        ONLY: vertex_param          !< this structure contains the fields
  USE mo_search_icongrid,       ONLY: walk_to_nc, find_nearest_vert

  USE mo_base_geometry,         ONLY: geographical_coordinates, &
       &                              cartesian_coordinates

  USE mo_additional_geometry,   ONLY: gc2cc

  USE mo_math_constants,        ONLY: rad2deg, deg2rad
  USE mo_physical_constants,    ONLY: re !< av. radius of the earth [m]
  USE netcdf

  IMPLICIT NONE

  PUBLIC :: reduce_grid

CONTAINS

  SUBROUTINE reduce_grid(topo_tiles_grid,  &
       &                 topo_files,       &
       &                 raw_data_orography_path)

    CHARACTER (LEN=*), INTENT(IN)            :: topo_files(1:max_tiles), &  !< filenames globe/aster raw data
         &                                      raw_data_orography_path

    TYPE(reg_lonlat_grid), INTENT(IN)        :: topo_tiles_grid(1:ntiles) !< structure w/ def of the raw data grid

    ! local variables
    REAL (KIND=wp), ALLOCATABLE              :: lon_tile(:), &   !< longitude coordinates of the GLOBE grid
         &                                      lon_red(:),& 
         &                                      lat_tile(:), &   !< latititude coordinates of the GLOBE grid
         &                                      hh_red(:), &   !< topographic height on reduced grid
         &                                      hh_red_to_write(:,:)
    
    REAL (KIND=wp)                           :: row_lat, &    ! latitude of the row for the topographic height array hh
         &                                      lon_diff, lontopo, &
         &                                      dxrat, &                                       ! ratio of dy to dx when given in m
         &                                      dlon0, &
         &                                      d2y,dx0,dy,&
         &                                      wgt, wgtsum   ! filter weights
       
    INTEGER(KIND=i4), ALLOCATABLE            :: ijlist(:), &
         &                                      hh(:), & !< 
         &                                      h_block(:,:)

    INTEGER(KIND=i4)                         :: undef_topo, &
         &                                      default_topo, &
         &                                      ncids_topo(ntiles), &
         &                                      nc_tile, &
         &                                      nr_tile, &
         &                                      nc_red, &
         &                                      i, j, & ! counters
         &                                      errorcode, &
         &                                      nt, &
         &                                      ij, np, istart, iend, &
         &                                      varid, &
         &                                      block_fix, &
         &                                      block_dyn, &
         &                                      idx, &
         &                                      mlat ! row number for GLOBE data

    CHARACTER (LEN=80)                       :: varname_topo  !< name of variable for topo data
    CHARACTER(len=filename_max)              :: topo_file

    CALL logging%info('Grid reduction: Enter routine: reduce grid')

    ! first open the raw topography  netcdf files
    DO nt=1,ntiles
      CALL open_netcdf_TOPO_tile(TRIM(raw_data_orography_path)//''//TRIM(topo_files(nt)), ncids_topo(nt))
    ENDDO

    CALL logging%info('Start processing topo tiles...')

    ! loop over all tiles
    DO nt=1,ntiles

      topo_file = TRIM(raw_data_orography_path)//TRIM(topo_files(nt))

      WRITE(message_text,*)TRIM(topo_file), ' --> to outfile'
      CALL logging%info(message_text)

      ! number of columns and rows of tile
      nc_tile= topo_tiles_grid(nt)%nlon_reg
      nr_tile= topo_tiles_grid(nt)%nlat_reg
      nc_red = 0

      ! block_fix -> nr of data rows loaded at once
      block_fix= 3000
      IF (block_fix > nr_tile) THEN
        CALL logging%warning('block size is bigger than nlat of tile -> set block_size to nlat of tile')
        block_fix = nr_tile
      ENDIF
      ! block_dyn used for last block
      block_dyn=block_fix

      ! allocate data
      ALLOCATE ( lon_tile(1:nc_tile), lon_red(1:nc_tile), ijlist(1:nc_tile) )
      ALLOCATE ( hh(0:nc_tile+1), hh_red(0:nc_tile+1) )
      ALLOCATE ( lat_tile(1:nr_tile) )
      ALLOCATE ( hh_red_to_write(0:nc_tile+1, nr_tile) )


      CALL get_fill_value(topo_file,undef_topo)
      default_topo = 0

      SELECT CASE(itopo_type)
        CASE(topo_aster)
          hh = default_topo
          hh_red = default_topo
        CASE(topo_gl)
          hh = undef_topo
          hh_red = undef_topo
      END SELECT

      ! calculate the longitude coordinate of the tile
      DO i =1,nc_tile
        lon_tile(i) = topo_tiles_grid(nt)%start_lon_reg + (i-1) * topo_tiles_grid(nt)%dlon_reg
      ENDDO

      ! calculate the latitiude coordinate of the tile
      DO j = 1, nr_tile
        lat_tile(j) = topo_tiles_grid(nt)%start_lat_reg + (j-1) * topo_tiles_grid(nt)%dlat_reg
      ENDDO

      dx0 =  topo_tiles_grid(nt)%dlon_reg * deg2rad * re ! longitudinal distance between to topo grid elemtens at equator
      dy = topo_tiles_grid(nt)%dlat_reg * deg2rad * re
      d2y = 2.0_wp * dy

      CALL get_varname(topo_file,varname_topo)

      !-----------------------------------------------------------------------------
      topo_rows: DO mlat=1,nr_tile

        ! load data in chunks of block_fix rows
        IF ( mlat == 1 .OR. MOD(mlat,block_fix) == 1 ) THEN

          ! check if last block -> adjust block site in that case
          IF (mlat + block_fix > nr_tile) THEN
            block_dyn= (nr_tile -mlat) +1
          ENDIF

          ! allocate array for block data
          IF(ALLOCATED(h_block)) THEN 
            DEALLOCATE(h_block, stat=errorcode)
            IF(errorcode/=0) CALL logging%error('cant deallocate the h_block',__FILE__,__LINE__)
          ENDIF
          ALLOCATE (h_block(1:nc_tile, 1:block_dyn), stat=errorcode)
          IF(errorcode/=0) CALL logging%error('cant allocate h_block',__FILE__,__LINE__)

          ! load data block
          CALL check_netcdf(nf90_inq_varid(ncids_topo(nt),TRIM(varname_topo),varid), __FILE__, __LINE__)
          CALL check_netcdf(nf90_get_var(ncids_topo(nt), varid, h_block, &
            &     start=(/1,mlat/),count=(/nc_tile,block_dyn/)), __FILE__, __LINE__)
        ENDIF

        ! determine row index from block
        IF(MOD(mlat,block_fix) > 0)THEN !normal case
          idx=MOD(mlat,block_fix)
        ELSE                            !last index before new block is created
          idx=block_dyn
        ENDIF

        hh(1:nc_tile) = h_block(1:nc_tile,idx) 
        hh(0)         = h_block(nc_tile,MOD(mlat,block_fix)) ! western wrap at -180/180 degree longitude
        hh(nc_tile+1) = h_block(1,MOD(mlat,block_fix))      ! eastern wrap at -180/180 degree longitude

        !dr test
        ! set undefined values to 0 altitude (default)
        WHERE (hh(:) == undef_topo)
          hh(:) = default_topo
        END WHERE

        ! absolute values of latitude at index mlat
        row_lat = topo_tiles_grid(nt)%start_lat_reg + (mlat-1) * topo_tiles_grid(nt)%dlat_reg

        ! compute hh_red
        dxrat = 1.0_wp/(COS(row_lat*deg2rad))
        nc_red = NINT(REAL(nc_tile,wp)/dxrat)
        dxrat = REAL(nc_tile,wp)/REAL(nc_red,wp)
        np = INT((dxrat-1)/2.0_wp) +1
        dlon0 = ABS(topo_tiles_grid(nt)%dlat_reg)*dxrat
        ijlist(:) = 0

        ! loop over columns of reduced grid
        DO i = 1, nc_red
          lon_red(i) = lon_tile(1)+(lon_tile(i)-lon_tile(1))*dxrat
          wgtsum = 0.0_wp
          hh_red(i) = 0.0_wp
          ij = INT(1+(i-1)*dxrat)
          ijlist(ij) = i
          istart = ij-np
          iend   = ij+np
          DO j = istart-1,iend+1
            ij = j
            IF (ij > nc_tile) THEN
              ij = ij - nc_tile
              lontopo = lon_tile(ij) + 360.0_wp
            ELSE IF (ij < 1) THEN
              ij = ij + nc_tile
              lontopo = lon_tile(ij) - 360.0_wp
            ELSE
              lontopo = lon_tile(ij)
            ENDIF
            lon_diff = ABS(lon_red(i)-lontopo)
            IF (lon_diff < 0.5_wp*dlon0) THEN
              wgt = MIN(5.0_wp,2.0_wp*dxrat,0.5_wp*dlon0/MAX(1.e-6_wp,lon_diff))
            ELSE IF (lon_diff-ABS(topo_tiles_grid(nt)%dlon_reg) < 0.5_wp*dlon0) THEN
              wgt = ABS(lon_diff-0.5_wp*dlon0)/ABS(topo_tiles_grid(nt)%dlon_reg)
            ELSE
              wgt = 0.0_wp
            ENDIF
            hh_red(i) = hh_red(i) + wgt*hh(ij)
            wgtsum = wgtsum + wgt
          ENDDO
          hh_red(i) = hh_red(i)/wgtsum
        ENDDO

        hh_red(0)        = hh_red(nc_red) ! western wrap at -180/180 degree longitude
        hh_red(nc_red+1) = hh_red(1)      ! eastern wrap at -180/180 degree longitude

        hh_red_to_write(:,mlat)=hh_red(:)
        

      ENDDO topo_rows

      WRITE(message_text,*)'MINAL hh_red_to_write: ' , MINVAL(hh_red_to_write)
      CALL logging%info(message_text)

      WRITE(message_text,*)'MAXVAL hh_red_to_write: ', MAXVAL(hh_red_to_write)
      CALL logging%info(message_text)

      !******************UWE TO DO******************
      ! make an interface to write hh_red_to_write
      ! and all the related meta-data for each tile
      !********************************************

      DEALLOCATE ( lon_tile ,lon_red, ijlist, &
           &       hh, hh_red, lat_tile,hh_red_to_write )

    ENDDO ! loop over topo tiles

    CALL logging%info('                       ...done')

    ! abort extpar during developing
    CALL logging%error('debug exit')

    CALL logging%info('Grid reduction: Exit routine: reduce grid')

    END SUBROUTINE reduce_grid

END MODULE mo_preproc_topo
