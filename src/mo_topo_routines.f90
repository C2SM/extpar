!+ Fortran module with routines and settings for GLOBE orography data
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V2_0         2013/06/04 Martina Messmer
!  Change all 'globe' to topo in globe_files, remove all 'globe' in
!  globe_tiles_lat/lon_min/max, change mo_GLOBE_data to mo_topo_data,
!  globe_tiles_grid to topo_tiles_grid, globe_tiles_ncolumns/nrow to
!  tiles_ncolumns/nrows, globe_files to topo_files, globe_grid to
!  topo_grid and change ntiles_gl to ntiles to obtain a more
!  dynamical code.
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module with routines and settings for GLOBE orography data
!> \author Hermann Asensio
!>
MODULE mo_topo_routines

  USE mo_logging
  USE mo_kind,                  ONLY: wp, i4
  USE mo_io_utilities,          ONLY: check_netcdf
                                
  USE netcdf,                   ONLY: nf90_open,       &
       &                              nf90_close,      &
       &                              nf90_inq_varid,  &
       &                              nf90_nowrite,    &
       &                              nf90_get_var
                                
  USE mo_utilities_extpar,      ONLY: free_un

  USE mo_grid_structures,       ONLY: reg_lonlat_grid
  USE mo_base_geometry,         ONLY: geographical_coordinates

  USE mo_topo_data,             ONLY: max_tiles, &
       &                              ntiles ,        &    !< GLOBE raw data has 16 tiles and ASTER has 13
       &                              tiles_lon_min,  &
       &                              tiles_lon_max,  &
       &                              tiles_lat_min,  &
       &                              tiles_lat_max,  &
       &                              tiles_ncolumns, &
       &                              tiles_nrows,    &
       &                              nc_tot,        &
       &                              nr_tot,        &
       &                              itopo_type,    &
       &                              topo_aster,    &
       &                              topo_gl,       &
       &                              aster_lat_min, &
       &                              aster_lat_max, &
       &                              aster_lon_min, &
       &                              get_varname, &
       &                              h_tile_row, &
       &                              aster_lon_max

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_topo_data_input_namelist,   &
       &    read_namelists_extpar_orography, &
       &    read_namelists_extpar_scale_sep, &
       &    det_topo_tiles_grid,             &
       &    det_topo_grid,                   &
       &    get_topo_tile_nr,                &
       &    get_topo_tile_block_indices,     &
       &    open_netcdf_topo_tile,           &
       &    close_netcdf_topo_tile,          &
       &    get_topo_data_band,              &
       &    get_topo_data_parallel

  PUBLIC :: det_band_gd, get_topo_data_block,get_topo_data_block_cosmo

  CONTAINS

  !---------------------------------------------------------------------------

  !> subroutine to read namelist for orography data settings for EXTPAR
  SUBROUTINE read_namelists_extpar_orography(namelist_file,          &
       &                                     raw_data_orography_path,&
       &                                     topo_files,             &  !mes>
       &                                     sgsl_files,             &
       &                                     ntiles_column,          &
       &                                     ntiles_row,             &
       &                                     itopo_type,             &
       &                                     lcompute_sgsl,          &
       &                                     lpreproc_oro,           &
       &                                     lsso_param,             &
       &                                     lsubtract_mean_slope,   &
       &                                     orography_buffer_file,  &
       &                                     orography_output_file,  &
       &                                     sgsl_buffer_file)


    CHARACTER (LEN=*), INTENT(IN)     :: namelist_file !< filename with namelists for for EXTPAR settings
    ! orography
    CHARACTER (LEN=1024), INTENT(OUT) :: raw_data_orography_path, &        !< path to raw data
         &                               topo_files(1:max_tiles), &        !< filenames globe raw data
         &                               sgsl_files(1:max_tiles)

    INTEGER (KIND=i4), INTENT(OUT)    :: ntiles_column, &      !< number of tile columns
         &                               ntiles_row, &         !< number of tile rows
         &                               itopo_type

    LOGICAL, INTENT(OUT)              :: lsso_param, &
         &                               lcompute_sgsl, &
         &                               lpreproc_oro, &
         &                               lsubtract_mean_slope

    CHARACTER (len=1024), INTENT(OUT) :: orography_buffer_file, &!< name for orography buffer file
         &                               orography_output_file, &!< name for orography output file
         &                               sgsl_buffer_file        !< name for sgsl output file

    INTEGER(KIND=i4)                  :: nuin, ierr, nzylen
    
    !> namelist for enable/disable subgrid-slope (SGSL)calculation
    NAMELIST /oro_runcontrol/      lcompute_sgsl
    
    !> namelist with information on orography data input
    NAMELIST /orography_raw_data/  itopo_type, lsso_param, lsubtract_mean_slope, &
         &                         raw_data_orography_path, ntiles_column, ntiles_row, & 
         &                         topo_files
       
    !> namelist with filenames for orography data output
    NAMELIST /orography_io_extpar/ orography_buffer_file, orography_output_file

    !> namelist with filenames for subgrid-slope (SGSL) data output
    NAMELIST /sgsl_io_extpar/      sgsl_files, sgsl_buffer_file, lpreproc_oro

    nuin = free_un()  ! function free_un returns free Fortran unit number

    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=oro_runcontrol, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist oro_runcontrol',__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=orography_io_extpar, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist orography_io_extpar',__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=orography_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist orography_raw_data',__FILE__, __LINE__) 
    ENDIF

    IF (lcompute_sgsl) THEN
      READ(nuin, NML=sgsl_io_extpar, IOSTAT=ierr)
      IF (ierr /= 0) THEN
        CALL logging%error('Cannot read in namelist sgsl_io_extpar',__FILE__, __LINE__) 
      ENDIF
    ENDIF
    
    CLOSE(nuin, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot close ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    nzylen=LEN_TRIM(raw_data_orography_path)
    IF( nzylen > 0 ) THEN
      IF( raw_data_orography_path(nzylen:nzylen) /= '/') THEN
        IF( nzylen < LEN(raw_data_orography_path) ) THEN
          raw_data_orography_path = raw_data_orography_path (1:nzylen)//'/'
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE read_namelists_extpar_orography

  !---------------------------------------------------------------------------

  !> subroutine to read namelist for scale separated data settings for EXTPAR
  SUBROUTINE read_namelists_extpar_scale_sep(namelist_file,           &
       &                                     raw_data_scale_sep_path, &
       &                                     scale_sep_files,         &
       &                                     lscale_separation)


    CHARACTER (len=1024), INTENT(IN)  :: namelist_file                !< filename with namelists for for EXTPAR settings

    CHARACTER (len=1024), INTENT(OUT) :: raw_data_scale_sep_path, &      !< path to raw data
         &                               scale_sep_files(1:max_tiles) !< filenames globe raw data

    LOGICAL, INTENT(OUT)              :: lscale_separation

    INTEGER(KIND=i4)                  :: nuin, ierr, nzylen

    !> namelist with information on scale separated data input
    NAMELIST /scale_separated_raw_data/ lscale_separation, raw_data_scale_sep_path, scale_sep_files

    nuin = free_un()  ! function free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=scale_separated_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist scale_separated_raw_data',__FILE__, __LINE__) 
    ENDIF

    CLOSE(nuin, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot close ', TRIM(namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    nzylen=LEN_TRIM(raw_data_scale_sep_path)
    IF( nzylen > 0 ) THEN
      IF( raw_data_scale_sep_path(nzylen:nzylen) /= '/') THEN
        IF( nzylen < LEN(raw_data_scale_sep_path) ) THEN
          raw_data_scale_sep_path = raw_data_scale_sep_path (1:nzylen)//'/'
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE read_namelists_extpar_scale_sep

  !-----------------------------------------------------------------------------

  !> read namelist with settings for GLOBE raw data grid
  !> or ASTER raw data grid                               !mes
  !> \author Hermann Asensio
  SUBROUTINE read_topo_data_input_namelist(input_namelist_file, topo_files)


    CHARACTER (LEN=1024), INTENT(IN)  :: input_namelist_file     !< file with input namelist
    CHARACTER (LEN=1024), INTENT(OUT) :: topo_files(1:max_tiles) !< filenames globe raw data

    INTEGER (KIND=i4)                 :: ierr, &   !< error flag
         &                               nuin, &   !< unit number
         &                               nfiles ! number of files

    !>Define the namelist group
    NAMELIST /GLOBE_files_info/ nfiles, topo_files

    nuin = free_un()  ! functioin free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(input_namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(input_namelist_file)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=GLOBE_files_info, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist GLOBE_files_info',__FILE__, __LINE__) 
    ENDIF
    CLOSE(nuin)

  END SUBROUTINE read_topo_data_input_namelist

  !-----------------------------------------------------------------------------

  !> determine GLOBE raw data grid
  !> \author Hermann Asensio
  SUBROUTINE det_topo_tiles_grid(topo_tiles_grid)


    TYPE(reg_lonlat_grid), INTENT(OUT) :: topo_tiles_grid(1:ntiles)

    INTEGER(KIND=i4)                   :: k ! counter

    REAL (KIND=wp)                     :: dlon, dlat

    ! determine the globe_tile_grid information from the namelist information
    DO k = 1, ntiles
      dlon =           (tiles_lon_max(k) - tiles_lon_min(k)) / REAL(tiles_ncolumns(k),wp)
      dlat = -1.0_wp * (tiles_lat_max(k) - tiles_lat_min(k)) / REAL(tiles_nrows(k),wp)

      ! latitude from north to south, negative increment
      topo_tiles_grid(k)%start_lon_reg  = tiles_lon_min(k) + 0.5 * dlon
      topo_tiles_grid(k)%end_lon_reg    = tiles_lon_max(k) - 0.5 * dlon

      ! latitude from north to south, note the negative increment!
      topo_tiles_grid(k)%start_lat_reg  = tiles_lat_max(k) + 0.5 * dlat
      ! latitude from north to south, note the negative increment!
      topo_tiles_grid(k)%end_lat_reg    = tiles_lat_min(k) - 0.5 * dlat

      topo_tiles_grid(k)%dlon_reg = dlon
      topo_tiles_grid(k)%dlat_reg = dlat

      topo_tiles_grid(k)%nlon_reg = tiles_ncolumns(k)
      topo_tiles_grid(k)%nlat_reg = tiles_nrows(k)
    ENDDO

  END SUBROUTINE det_topo_tiles_grid

  !-----------------------------------------------------------------------------

  !> determine complete(global) GLOBE raw data grid
  !> \author Hermann Asensio
  SUBROUTINE det_topo_grid(topo_grid)

    TYPE(reg_lonlat_grid), INTENT(OUT) :: topo_grid !< structure with definition of the global data grid of the GLOBE data

    REAL (KIND=wp)                     :: dlon, dlat

    !mes > as ASTER does not cover the whole globe until now different procedures must be chosen for ASTER and GLOBE
    SELECT CASE(itopo_type)
      CASE(topo_aster)

        dlon = (aster_lon_max - aster_lon_min) / REAL(nc_tot,wp)

        dlat = -1. * (aster_lat_max - aster_lat_min) / REAL(nr_tot,wp)
        ! latitude from north to south, negative increment

        topo_grid%start_lon_reg  =  aster_lon_min + 0.5_wp * dlon
        topo_grid%end_lon_reg    =  aster_lon_max - 0.5_wp * dlon

        topo_grid%start_lat_reg = aster_lat_max + 0.5_wp * dlat ! latitude from north to south, note the negative increment!
        topo_grid%end_lat_reg  =  aster_lat_min - 0.5_wp * dlat ! latitude from north to south, note the negative increment!

      CASE(topo_gl)

        dlon = 360._wp / REAL(nc_tot,wp)
        dlat = -1. * 180._wp / REAL(nr_tot,wp)
        ! latitude from north to south, negative increment

        topo_grid%start_lon_reg  = -180._wp + 0.5_wp * dlon
        topo_grid%end_lon_reg    =  180._wp - 0.5_wp * dlon

        topo_grid%start_lat_reg = 90._wp + 0.5_wp * dlat ! latitude from north to south, note the negative increment!
        topo_grid%end_lat_reg  = -90._wp - 0.5_wp * dlat ! latitude from north to south, note the negative increment!

    END SELECT

    topo_grid%dlon_reg = dlon
    topo_grid%dlat_reg = dlat

    topo_grid%nlon_reg = nc_tot
    topo_grid%nlat_reg = nr_tot

  END SUBROUTINE det_topo_grid

  !-----------------------------------------------------------------------------

  !> determine grid description of band for GLOBE I/O
  !> \author Hermann Asensio
  SUBROUTINE det_band_gd(topo_grid,start_topo_row, ta_grid)

    TYPE(reg_lonlat_grid), INTENT(IN)  :: topo_grid      !< structure with definition of the global data grid of the GLOBE data
    INTEGER(KIND=i4), INTENT(IN)       :: start_topo_row !< number of the start row of band of topo_grid (global domain)
    TYPE(reg_lonlat_grid), INTENT(OUT) :: ta_grid        !< structure with defenition of the target area grid

    INTEGER(KIND=i4)                   :: nrows = 500 !< number of rows, set to 1000 as default

    ! band from east to west for the whole globe, like the complete topo_grid
    ta_grid%dlon_reg = topo_grid%dlon_reg
    ta_grid%dlat_reg = topo_grid%dlat_reg

    ta_grid%start_lon_reg = topo_grid%start_lon_reg
    ta_grid%end_lon_reg =  topo_grid%end_lon_reg
    ta_grid%nlon_reg = topo_grid%nlon_reg

    ! latitude from north to south, negative increment
    ta_grid%nlat_reg = nrows
    ta_grid%start_lat_reg = topo_grid%start_lat_reg + ta_grid%dlat_reg * (start_topo_row - 1)
    ! latitude from north to south, note the negative increment!
    ta_grid%end_lat_reg  =  ta_grid%start_lat_reg + ta_grid%dlat_reg * (nrows - 1)

    ! latitude from north to south, note the negative increment!
    ! check for south pole
    IF (ta_grid%end_lat_reg < topo_grid%end_lat_reg) THEN ! band is at south pole
      ta_grid%end_lat_reg =  topo_grid%end_lat_reg
      ta_grid%nlat_reg =  NINT(((ta_grid%end_lat_reg - ta_grid%start_lat_reg) / ta_grid%dlat_reg)) + 1
    ENDIF

  END SUBROUTINE det_band_gd

  !-----------------------------------------------------------------------------

  !> find GLOBE tile for given geographical coordinates
  ELEMENTAL FUNCTION get_topo_tile_nr(point_geo) RESULT (index_k)

    TYPE(geographical_coordinates), INTENT(IN) :: point_geo !< geographical coordinats of a point [degrees]

    INTEGER (KIND=i4) :: index_k !< index of GLOBE tile which contains point_geo

    ! local variables

    INTEGER (KIND=i4) :: t_i
    INTEGER (KIND=i4) :: t_j

    REAL (KIND=wp)    :: lon0_t
    REAL (KIND=wp)    :: lat0_t
    REAL (KIND=wp)    :: dlon_t
    REAL (KIND=wp)    :: dlat_t

    REAL (KIND=wp)    :: point_lon_coor

    ! the GLOBE data are diveded in 16 tiles,
    ! this defines a "dummy grid" to determine the index with a function
    lon0_t = -180._wp
    lat0_t =  100._wp
    dlon_t =   90._wp
    dlat_t =  -50._wp

    point_lon_coor = point_geo%lon
    IF (point_lon_coor > 180._wp) THEN  ! shift longitude range
      point_lon_coor = point_lon_coor -360._wp
    ENDIF

    t_i = INT((point_lon_coor - lon0_t)/dlon_t) + 1 ! get the tile index for the column

    t_j = INT((lat0_t - point_geo%lat)/dlat_t) + 1  ! get the tile index for the row,
    !note the negative increment (rows from north to south

    index_k = (t_j - 1) * 4 + t_i ! the way the 16 element array is sorted (columns first)

  END FUNCTION get_topo_tile_nr
  !----------------------------------------------------------------------------------------------------------------

  !> get startrow, endrow, startcolumn and endcolumn of each GLOBE tile (raw data) for a
  !! given target area (ta_grid) and
  !! get start_indices (lon, lat) and end_indices of the target area for each GLOBE tile
  !! The GLOBE raw data are split in 16 tiles (ASTER in 36), so the target area may overlap several tiles.
  !! This subroutine determines the necesarry indices to read in the GLOBE/ASTER data into the
  !! target area.
  !! GLOBE/ASTER tiles which are outside the target block will get indices with the value '0'

  SUBROUTINE get_topo_tile_block_indices(ta_grid,          &
       &                                 topo_tiles_grid,  &
       &                                 topo_startrow,    &
       &                                 topo_endrow,      &
       &                                 topo_startcolumn, &
       &                                 topo_endcolumn,   &
       &                                 ta_start_ie,      &
       &                                 ta_end_ie,        &
       &                                 ta_start_je,      &
       &                                 ta_end_je)


    !< structure with definition of the target area grid (dlon must be the same as for the whole GLOBE dataset)
    TYPE(reg_lonlat_grid), INTENT(IN) :: ta_grid, &
        &                                topo_tiles_grid(1:ntiles)

    INTEGER (KIND=i4), INTENT(OUT)    :: topo_startrow(1:ntiles), &    !< startrow indices for each GLOBE tile
        &                                topo_endrow(1:ntiles), &      !< endrow indices for each GLOBE tile
        &                                topo_startcolumn(1:ntiles), &  !< starcolumn indices for each GLOBE tile
        &                                topo_endcolumn(1:ntiles), &   !< endcolumn indices for each GLOBE tile
        &                                ta_start_ie(1:ntiles), &
        &                                ta_end_ie(1:ntiles), &
        &                                ta_start_je(1:ntiles), &
        &                                ta_end_je(1:ntiles)

    ! local variables
    INTEGER(KIND=i4)                  :: i,j,k,m,n,o,undefined, &
        &                                startrow, & ! startrow for tile
        &                                endrow, &
        &                                startcolumn, &
        &                                endcolumn

    REAL (KIND=wp)                    :: dlon, dlat

    undefined         = 0
    topo_startrow     = undefined
    topo_endrow       = undefined
    topo_startcolumn  = undefined
    topo_endcolumn    = undefined
    ta_start_ie       = undefined
    ta_end_ie         = undefined
    ta_start_je       = undefined
    ta_end_je         = undefined

    k    = 1 ! determin dlon and dlat (are the same for all tiles)
    dlon = ta_grid%dlon_reg
    dlat = ta_grid%dlat_reg

    !mes > SELECT CASE as the two DEMs do not have the same amount of tiles.
    SELECT CASE(itopo_type)
      CASE(topo_aster)
        m = 1
        n = 1
        o = ntiles
      CASE(topo_gl)
        m = 1
        n = ntiles/4
        o = ntiles/4
    END SELECT

    DO j = m, n
      DO i = 1, o
        k = (j - 1) * 4 + i ! the way the 16 element array (13 element array) is sorted (columns first/ only rows)

        ! get startcolumn for tile k
        startcolumn = NINT((ta_grid%start_lon_reg - topo_tiles_grid(k)%start_lon_reg)/dlon) + 1
        !< here I want nearest index (NINT)
        IF (startcolumn < 1) THEN
          topo_startcolumn(k) = 1
          ! get the start index of the subtile for the target area block
          ta_start_ie(k) = NINT ((topo_tiles_grid(k)%start_lon_reg - ta_grid%start_lon_reg)/dlon) + 1
          !< index of target area block
        ELSE IF (startcolumn > tiles_ncolumns(k)) THEN
          topo_startcolumn(k) = 0
          ta_start_ie(k) = 0
        ELSE
          topo_startcolumn(k) = startcolumn
          ta_start_ie(k) = 1
        ENDIF

        ! get endcolumn for tile k
        endcolumn = NINT((ta_grid%end_lon_reg - topo_tiles_grid(k)%start_lon_reg)/dlon) +1
        IF (endcolumn > tiles_ncolumns(k)) THEN
          topo_endcolumn(k) = tiles_ncolumns(k)
          ! get the end index of the subtile for the target area block
          ta_end_ie(k) = NINT ((topo_tiles_grid(k)%end_lon_reg - ta_grid%start_lon_reg)/dlon) + 1
          !< index of target area block
        ELSE IF (endcolumn < 1) THEN
          topo_endcolumn(k) = 0
          ta_end_ie(k) = 0
        ELSE
          topo_endcolumn(k) = endcolumn
          ta_end_ie(k) = ta_grid%nlon_reg
        ENDIF

        ! get startrow for tile k
        startrow = NINT((ta_grid%start_lat_reg - topo_tiles_grid(k)%start_lat_reg)/dlat) + 1

        IF (startrow < 1) THEN
          topo_startrow(k) = 1
          ! get the start index of the subtile for the target area block
          ta_start_je(k) = NINT ((topo_tiles_grid(k)%start_lat_reg  - ta_grid%start_lat_reg)/dlat) + 1
          !< index of target area block

        ELSE IF (startrow > tiles_nrows(k)) THEN
          topo_startrow(k) = 0
          ta_start_je(k) = 0
        ELSE
          topo_startrow(k) = startrow
          ta_start_je(k) = 1
        ENDIF

        ! get endrow for tile k
        endrow   = NINT(( ta_grid%end_lat_reg - topo_tiles_grid(k)%start_lat_reg )/dlat)  + 1

        IF (endrow > tiles_nrows(k)) THEN
          topo_endrow(k) = tiles_nrows(k)
          ! get the start index of the subtile for the target area block
          ta_end_je(k) = NINT ((topo_tiles_grid(k)%end_lat_reg -  ta_grid%start_lat_reg )/dlat) + 1
          !< index of target area block

        ELSE IF (endrow < 1) THEN
          topo_endrow(k) = 0
          ta_end_je(k) = 0
        ELSE
          topo_endrow(k) = endrow
          ta_end_je(k) =  ta_grid%nlat_reg
        ENDIF

      ENDDO
    ENDDO  ! loop over the tiles

  END SUBROUTINE get_topo_tile_block_indices

  !----------------------------------------------------------------------------------------------------------------

  !> open netcdf-file and get netcdf unit file number
  SUBROUTINE open_netcdf_topo_tile(path_topo_tile, ncid)

    CHARACTER (len=*), INTENT(in) :: path_topo_tile         !< filename with path to GLOBE tile
    INTEGER, INTENT(out)          :: ncid                   !< netcdf unit file number

    !! open netcdf file
    call check_netcdf( nf90_open(TRIM(path_topo_tile),NF90_NOWRITE, ncid), __FILE__, __LINE__)

  END SUBROUTINE open_netcdf_TOPO_tile

  !> close netcdf-file
  SUBROUTINE close_netcdf_TOPO_tile(ncid)
    INTEGER, INTENT(in) :: ncid                             !< netcdf unit file number

    !! close netcdf file
    call check_netcdf( nf90_close( ncid))

  END SUBROUTINE close_netcdf_TOPO_tile

  !----------------------------------------------------------------------------------------------------------------

  !> get globe data on a single circle of latitude
  SUBROUTINE get_topo_data_parallel(mlat,       &
       &                            ncids_topo, &
       &                            h_parallel)


    INTEGER(KIND=i4) , INTENT(IN)  :: mlat, &  !< global index of raw data line
         &                            ncids_topo(1:ntiles)
    !< ncid for the GLOBE tiles, the netcdf files have to be opened by a previous call of open_netcdf_GLOBE_tile
    INTEGER (KIND=i4), INTENT(OUT) :: h_parallel(1:nc_tot)     !< GLOBE altitude data along a parallel

    ! local variables
    INTEGER(KIND=i4)               :: tile_start, &
         &                            tile_end, &
         &                            tile_row, &
         &                            varid, &    !< id of variable
         &                            k, &  !< counter
         &                            os, & !< counter
         &                            nt ! counter

    CHARACTER (LEN=80)            :: varname  !< name of variable

    varname = 'altitude'  ! I know that in the GLOBE netcdf files the height data are stored in a variable "altitude"

    SELECT CASE(mlat)
      CASE (1:4800)
        tile_start = 1    ! GLOBE TILE A, or 1
        tile_row   = mlat ! row in the Tiles A, B, C, D
      CASE (4801:10800)
        tile_start = 5    ! GLOBE TILE E, or 5
        tile_row   = mlat - 4800  ! row in the tiles E, F, G, H
      CASE (10801:16800)
        tile_start = 9    ! GLOBE TILE I, or 9
        tile_row   = mlat - 10800 ! row in the tiles I, J, K, L
      CASE (16801:21600)
        tile_start = 13
        tile_row   = mlat - 16800 ! row in the tiles M, N, O, P
      CASE DEFAULT
        CALL logging%error('get_topo_data_parallel: mlat not in data range of TOPO tiles',__FILE__,__LINE__)
    END SELECT

    tile_end = tile_start + 3 ! numbering of GLOBE tiles
    nt = 0
    DO k=tile_start,tile_end
      nt = nt + 1 ! count the number of tiles
      CALL check_netcdf(nf90_inq_varid(ncids_topo(k),TRIM(varname),varid)) ! get the varid of the altitude variable
      CALL check_netcdf(nf90_get_var(ncids_topo(k), varid,  h_tile_row, & ! get the data of one tile row
           &    start=(/1,tile_row/),count=(/nc_tile,1/)))
      os = (nt-1) * nc_tile ! offset for array with data along latitude circle
      h_parallel(os+1:os+nc_tile) = h_tile_row(1:nc_tile)
    ENDDO

  END SUBROUTINE get_topo_data_parallel

  !> get GLOBE data block for a given target area from the tile block indices
  SUBROUTINE get_topo_data_block(str_topo,     &   !mes ><
       &                         ta_grid,         &
       &                         topo_tiles_grid, &
       &                         ncids_topo,      &
       &                         h_block,         &
       &                         lis_icon)

    CHARACTER (len=*), INTENT(IN)     :: str_topo
    LOGICAL, INTENT(IN)               :: lis_icon
    TYPE(reg_lonlat_grid), INTENT(IN) :: ta_grid, &
         &                               topo_tiles_grid(1:ntiles)

    INTEGER(KIND=i4) , INTENT(IN)     :: ncids_topo(1:ntiles)

    INTEGER(KIND=i4), INTENT(OUT)     :: h_block(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg)

    !local variables
    INTEGER (KIND=i4)                 :: topo_startrow(1:ntiles), &    !< startrow indices for each GLOBE tile
         &                               topo_endrow(1:ntiles), &      !< endrow indices for each GLOBE tile
         &                               topo_startcolumn(1:ntiles), & !< starcolumn indices for each GLOBE tile
         &                               topo_endcolumn(1:ntiles), &   !< endcolumn indices for each GLOBE tile
         &                               ta_start_ie(1:ntiles), &      !< indices of target area block for first column of each GLOBE tile
         &                               ta_end_ie(1:ntiles), &        !< indices of target area block for last column of each GLOBE tile
         &                               ta_start_je(1:ntiles), &      !< indices of target area block for first row of each GLOBE tile
         &                               ta_end_je(1:ntiles), &        !< indices of target area block for last row of each GLOBE tile
         &                               varid, nrows, ncolumns, k, errorcode

    INTEGER (KIND=i4), ALLOCATABLE    :: raw_topo_block(:,:) !< a block with GLOBE data

    CHARACTER (LEN=80):: varname                    !< name of variable

    IF (lis_icon) THEN
      varname = str_topo
    ELSE
      CALL get_varname(str_topo,varname)
    ENDIF

    CALL get_topo_tile_block_indices( ta_grid,          &
         &                            topo_tiles_grid,  &
         &                            topo_startrow,    &
         &                            topo_endrow,      &
         &                            topo_startcolumn, &
         &                            topo_endcolumn,   &
         &                            ta_start_ie,      &
         &                            ta_end_ie,        &
         &                            ta_start_je,      &
         &                            ta_end_je)

    DO k = 1, ntiles
      IF ((topo_startrow(k)/=0).AND.(topo_startcolumn(k)/=0)) THEN
        nrows = topo_endrow(k) - topo_startrow(k) + 1
        ncolumns = topo_endcolumn(k) - topo_startcolumn(k) + 1
        ALLOCATE (raw_topo_block(1:ncolumns,1:nrows), STAT=errorcode)
        IF(errorcode/=0) CALL logging%error('Cant allocate the array raw_topo_block',__FILE__,__LINE__)

        CALL check_netcdf(nf90_inq_varid(ncids_topo(k),TRIM(varname),varid), __FILE__, __LINE__)
        ! get the data into the raw_topo_block
        CALL check_netcdf(nf90_get_var(ncids_topo(k), varid,  raw_topo_block,     &
             &     start=(/topo_startcolumn(k),topo_startrow(k)/),count=(/ncolumns,nrows/)), __FILE__, __LINE__)
        h_block(ta_start_ie(k):ta_end_ie(k),ta_start_je(k):ta_end_je(k)) = raw_topo_block(1:ncolumns,1:nrows)

        DEALLOCATE (raw_topo_block, STAT=errorcode)
        IF(errorcode/=0) CALL logging%error('Cant deallocate the array raw_topo_block',__FILE__,__LINE__)

      ENDIF
    ENDDO

  END SUBROUTINE get_topo_data_block

  !----------------------------------------------------------------------------------------------------------------

  !> get GLOBE data block for a given target area from the tile block indices
  SUBROUTINE get_topo_data_block_cosmo(topo_file_1,     &   !mes ><
       &                         ta_grid,         &
       &                         topo_tiles_grid, &
       &                         ncids_topo,      &
       &                         h_block)


    CHARACTER (len=*), INTENT(IN)     :: topo_file_1

    TYPE(reg_lonlat_grid), INTENT(IN) :: ta_grid
    !< structure with definition of the target area grid (dlon must be the same as for the whole GLOBE dataset)
 
    TYPE(reg_lonlat_grid), INTENT(IN) :: topo_tiles_grid(1:ntiles)
    !< structure with defenition of the raw data grid for the 16 GLOBE tiles

    INTEGER , INTENT(IN)              :: ncids_topo(1:ntiles)
    !< ncid for the GLOBE tiles, the netcdf files have to be opened by a previous call of open_netcdf_GLOBE_tile

    INTEGER (KIND=i4), INTENT(OUT)    :: h_block(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg)
    !< a block of GLOBE altitude data

    !local variables
    INTEGER (KIND=i4)                 :: topo_startrow(1:ntiles), &    !< startrow indices for each GLOBE tile
         &                               topo_endrow(1:ntiles), &      !< endrow indices for each GLOBE tile
         &                               topo_startcolumn(1:ntiles), & !< starcolumn indices for each GLOBE tile
         &                               topo_endcolumn(1:ntiles), &   !< endcolumn indices for each GLOBE tile
         &                               ta_start_ie(1:ntiles), &      !< indices of target area block for first column of each GLOBE tile
         &                               ta_end_ie(1:ntiles), &        !< indices of target area block for last column of each GLOBE tile
         &                               ta_start_je(1:ntiles), &      !< indices of target area block for first row of each GLOBE tile
         &                               ta_end_je(1:ntiles), &        !< indices of target area block for last row of each GLOBE tile
         &                               nrows, &                      !< number of rows ! dimensions for raw_topo_block
         &                               ncolumns, &                   !< number of columns ! dimensions for raw_topo_block
         &                               k, &                          !< counter
         &                               errorcode, &                  !< error status variable
         &                               varid                      !< id of variable

    INTEGER (KIND=i4), ALLOCATABLE    :: raw_topo_block(:,:) !< a block with GLOBE data

    CHARACTER (LEN=80)                :: varname                    !< name of variable


    CALL get_varname(topo_file_1,varname)
    !       varname = 'altitude'  ! I know that in the GLOBE netcdf files the height data are stored in a variable "altitude"

    CALL get_topo_tile_block_indices( ta_grid,          &
         &                            topo_tiles_grid,  &
         &                            topo_startrow,    &
         &                            topo_endrow,      &
         &                            topo_startcolumn, &
         &                            topo_endcolumn,   &
         &                            ta_start_ie,      &
         &                            ta_end_ie,        &
         &                            ta_start_je,      &
         &                            ta_end_je)
    !  allocate_raw_topo_fields(nrows,ncolumns)
    !  raw_topo_block

    DO k = 1, ntiles
      IF ((topo_startrow(k)/=0).AND.(topo_startcolumn(k)/=0)) THEN
        nrows = topo_endrow(k) - topo_startrow(k) + 1
        ncolumns = topo_endcolumn(k) - topo_startcolumn(k) + 1
        ALLOCATE (raw_topo_block(1:ncolumns,1:nrows), STAT=errorcode)
        IF(errorcode/=0) CALL logging%error('Cant allocate the array raw_topo_block',__FILE__,__LINE__)

        CALL check_netcdf(nf90_inq_varid(ncids_topo(k),TRIM(varname),varid), __FILE__, __LINE__)
        ! get the data into the raw_topo_block
        CALL check_netcdf(nf90_get_var(ncids_topo(k), varid,  raw_topo_block,     &
             &     start=(/topo_startcolumn(k),topo_startrow(k)/),count=(/ncolumns,nrows/)), __FILE__, __LINE__)
        h_block(ta_start_ie(k):ta_end_ie(k),ta_start_je(k):ta_end_je(k)) = raw_topo_block(1:ncolumns,1:nrows)

        DEALLOCATE (raw_topo_block, STAT=errorcode)
        IF(errorcode/=0) CALL logging%error('Cant deallocate the array raw_topo_block',__FILE__,__LINE__)

      ENDIF
    ENDDO

  END SUBROUTINE get_topo_data_block_cosmo

  !> get globe data band on a circle of latitude
  SUBROUTINE get_topo_data_band(mstart,     &
       &                        nrows,      &
       &                        ncids_topo, &
       &                        h_band)


    INTEGER(KIND=i4) , INTENT(IN) :: mstart, &      !< global index of first raw data line
         &                           nrows, &       !< total number or row data rows to read in
         &                           ncids_topo(1:ntiles)

    !< ncid for the GLOBE tiles, the netcdf files have to be opened by a previous call of open_netcdf_TOPO_tile
    INTEGER (KIND=i4), INTENT(OUT):: h_band(1:nc_tot,1:nrows)
    !< GLOBE altitude data along a parallel

    ! local variables
    INTEGER(KIND=i4)              :: tile_start, &
         &                           tile_end, &
         &                           tile_row, &
         &                           varid, &               !< id of variable
         &                           k, &      !< counter
         &                           os, &     !< counter
         &                           nt, &     !< counter
         &                           n_row, &  !< counter
         &                           mlat, &   !< global index of GLOBE raw data row to read in
         &                           m_end  !< global index of last raw data line

    CHARACTER (LEN=80)            :: varname  !< name of variable

    m_end = mstart+nrows

    varname = 'altitude'  ! I know that in the GLOBE netcdf files the height data are stored in a variable "altitude"
    SELECT CASE(mstart)
      CASE (1:4800)
        tile_start = 1    ! GLOBE TILE A, or 1
        !tile_row   = mlat ! row in the Tiles A, B, C, D
      CASE (4801:10800)
        tile_start = 5    ! GLOBE TILE E, or 5
        !tile_row   = mlat - 4800  ! row in the tiles E, F, G, H
      CASE (10801:16800)
        tile_start = 9    ! GLOBE TILE I, or 9
        !tile_row   = mlat - 10800 ! row in the tiles I, J, K, L
      CASE (16801:21600)
        tile_start = 13
        !tile_row   = mlat - 16800 ! row in the tiles M, N, O, P
      CASE DEFAULT
        CALL logging%error('get_topo_data_band: mlat not in data range of TOPO tiles',__FILE__,__LINE__)
    END SELECT

    SELECT CASE(m_end)
      CASE (1:4800)
        tile_end = 4    ! row in the Tiles A, B, C, D
      CASE (4801:10800)
        tile_end = 8    ! row in the tiles E, F, G, H
      CASE (10801:16800)
        tile_end = 12   ! row in the tiles I, J, K, L
      CASE (16801:21600)
        tile_end = 16   ! row in the tiles M, N, O, P
      CASE DEFAULT
        CALL logging%error('get_topo_data_band: mlat not in data range of TOPO tiles',__FILE__,__LINE__)
    END SELECT

    DO n_row=1,nrows
      mlat= mstart + n_row -1 ! global index of GLOBE row

      SELECT CASE(mlat)
        CASE (1:4800)
          tile_start = 1    ! GLOBE TILE A, or 1
          tile_row   = mlat ! row in the Tiles A, B, C, D
        CASE (4801:10800)
          tile_start = 5    ! GLOBE TILE E, or 5
          tile_row   = mlat - 4800  ! row in the tiles E, F, G, H
        CASE (10801:16800)
          tile_start = 9    ! GLOBE TILE I, or 9
          tile_row   = mlat - 10800 ! row in the tiles I, J, K, L
        CASE (16801:21600)
          tile_start = 13
          tile_row   = mlat - 16800 ! row in the tiles M, N, O, P
        CASE DEFAULT
          CALL logging%error('get_topo_data_band: mlat not in data range of TOPO tiles',__FILE__,__LINE__)
      END SELECT

      tile_end = tile_start + 3 ! numbering of GLOBE tiles
      nt = 0
      DO k=tile_start,tile_end
        nt = nt + 1 ! count the number of tiles
        CALL check_netcdf(nf90_inq_varid(ncids_topo(k),TRIM(varname),varid)) ! get the varid of the altitude variable
        CALL check_netcdf(nf90_get_var(ncids_topo(k), varid,  h_tile_row,     & ! get the data of one tile row
             start=(/1,tile_row/),count=(/nc_tile,1/)))
        os = (nt-1) * nc_tile ! offset for array with data along latitude circle
        h_band(os+1:os+nc_tile,n_row) = h_tile_row(1:nc_tile)
      ENDDO
    ENDDO

  END SUBROUTINE get_topo_data_band

END MODULE mo_topo_routines
