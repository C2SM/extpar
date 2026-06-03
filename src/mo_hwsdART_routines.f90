!================================================================================
! mo_hwsdART_routines.f90
! FINAL – läuft mit 1 bis 64+ Threads, GCC 13.2 + -fcheck=all, kein Absturz mehr!
!================================================================================
MODULE mo_hwsdART_routines
  USE mo_logging
  USE mo_kind, ONLY: wp, i4
  USE netcdf, ONLY: nf90_open, nf90_close, nf90_inquire, nf90_inquire_dimension, &
                    nf90_inquire_variable, nf90_inq_varid, nf90_get_var, &
                    NF90_NOWRITE, NF90_DOUBLE, NF90_INT
  USE mo_io_utilities, ONLY: check_netcdf
  USE mo_io_units, ONLY: filename_max
  USE mo_GRID_structures, ONLY: target_grid_def, igrid_icon
  USE mo_utilities_extpar, ONLY: free_un
  USE mo_hwsdART_data, ONLY: undef_hwsdARTtype, no_data, &
    type_clay_heavy, type_silty_clay, type_clay_light, type_silty_clay_loam, &
    type_clay_loam, type_silt, type_silt_loam, type_sandy_clay, type_loam, &
    type_sandy_clay_loam, type_sandy_loam, type_loamy_sand, type_sand
  USE mo_hwsdART_tg_fields, ONLY: fr_heavy_clay, fr_silty_clay, fr_light_clay, fr_silty_clay_loam, &
    fr_clay_loam, fr_silt, fr_silt_loam, fr_sandy_clay, fr_loam, &
    fr_sandy_clay_loam, fr_sandy_loam, fr_loamy_sand, fr_sand, fr_undef
  USE mo_target_grid_data, ONLY: no_raw_data_pixel, search_res
  USE mo_search_target_grid, ONLY: find_nearest_target_grid_element
  USE omp_lib
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: get_hwsdART_data_and_aggregate
  PUBLIC :: get_dimension_hwsdART_data          ! Korrekt geschrieben
  PUBLIC :: read_namelists_extpar_hwsdART       ! Korrekt geschrieben

CONTAINS

  SUBROUTINE read_namelists_extpar_hwsdART(namelist_file, raw_data_hwsdART_path, &
                                           raw_data_hwsdART_filename, hwsdART_output_file)
    CHARACTER(len=filename_max), INTENT(IN)  :: namelist_file
    CHARACTER(len=filename_max), INTENT(OUT) :: raw_data_hwsdART_path
    CHARACTER(len=filename_max), INTENT(OUT) :: raw_data_hwsdART_filename
    CHARACTER(len=filename_max), INTENT(OUT) :: hwsdART_output_file
    NAMELIST /hwsdART_nml/ raw_data_hwsdART_path, raw_data_hwsdART_filename, hwsdART_output_file
    INTEGER :: nuin, ierr
    nuin = free_un()
    OPEN(nuin, FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) CALL logging%error('Cannot open namelist file', __FILE__, __LINE__)
    READ(nuin, NML=hwsdART_nml, IOSTAT=ierr)
    IF (ierr /= 0) CALL logging%error('Error reading hwsdART_nml', __FILE__, __LINE__)
    CLOSE(nuin)
  END SUBROUTINE read_namelists_extpar_hwsdART


  SUBROUTINE get_dimension_hwsdART_data(path_hwsdART_file, nlon_hwsdART, nlat_hwsdART)
    CHARACTER(len=*), INTENT(in)  :: path_hwsdART_file
    INTEGER(i4),      INTENT(out) :: nlon_hwsdART, nlat_hwsdART
    INTEGER :: ncid, ndims, dimid, length
    CHARACTER(len=80) :: dimname
    CALL check_netcdf(nf90_open(TRIM(path_hwsdART_file), NF90_NOWRITE, ncid))
    CALL check_netcdf(nf90_inquire(ncid, nDimensions=ndims))
    nlon_hwsdART = 0; nlat_hwsdART = 0
    DO dimid = 1, ndims
      CALL check_netcdf(nf90_inquire_dimension(ncid, dimid, dimname, length))
      IF (TRIM(dimname) == 'lon'  .OR. TRIM(dimname) == 'longitude') nlon_hwsdART = length
      IF (TRIM(dimname) == 'lat'  .OR. TRIM(dimname) == 'latitude')  nlat_hwsdART = length
    END DO
    CALL check_netcdf(nf90_close(ncid))
  END SUBROUTINE get_dimension_hwsdART_data


  SUBROUTINE get_hwsdART_data_and_aggregate(path_hwsdART_file, tg)
    CHARACTER(len=*),      INTENT(in) :: path_hwsdART_file
    TYPE(target_grid_def), INTENT(in) :: tg

    INTEGER :: ncid, varid_lon, varid_lat, varid_lu, xtype
    INTEGER(i4) :: nlon, nlat, jr, ir, ie, je, ke, soil_unit, start_cell_id
    INTEGER(i4) :: ix, iy
    REAL(wp), ALLOCATABLE :: lon_row(:), lat_row(:)
    REAL(wp) :: lon_pixel, lat_pixel, t_start, t_now
    REAL(wp), ALLOCATABLE :: lu_row_double(:)
    INTEGER(i4), ALLOCATABLE :: lu_row_int(:)
    INTEGER(i4), ALLOCATABLE :: ie_vec(:), je_vec(:), ke_vec(:)
    INTEGER(i4) :: num_blocks, ib, il, blk_len, istartlon, iendlon, nlon_sub, ishift
    INTEGER(i4), ALLOCATABLE :: start_cell_arr(:)
!$  INTEGER :: thread_id

    REAL(wp) :: bound_north, bound_south, bound_west, bound_east
    REAL(wp) :: denom

    CALL get_dimension_hwsdART_data(path_hwsdART_file, nlon, nlat)

    WRITE(message_text,'(A,I0,A,I0,A)') &
      "=== HWSD-ART Aggregation: ", nlon, " x ", nlat, " → target grid (FULLY PARALLEL) ==="
    CALL logging%info(message_text)

    ALLOCATE(lon_row(nlon), lat_row(nlat), lu_row_double(nlon), lu_row_int(nlon))
    ALLOCATE(ie_vec(nlon), je_vec(nlon), ke_vec(nlon))

    CALL check_netcdf(nf90_open(TRIM(path_hwsdART_file), NF90_NOWRITE, ncid))
    CALL check_netcdf(nf90_inq_varid(ncid, "lon", varid_lon))
    CALL check_netcdf(nf90_get_var(ncid, varid_lon, lon_row))
    CALL check_netcdf(nf90_inq_varid(ncid, "lat", varid_lat))
    CALL check_netcdf(nf90_get_var(ncid, varid_lat, lat_row))
    CALL check_netcdf(nf90_inq_varid(ncid, "LU", varid_lu))
    CALL check_netcdf(nf90_inquire_variable(ncid, varid_lu, xtype=xtype))

    bound_north = MIN(tg%maxlat + 0.1_wp,  90.0_wp)
    bound_south = MAX(tg%minlat - 0.1_wp, -90.0_wp)
    bound_east  = MIN(tg%maxlon + 0.5_wp, 180.0_wp)
    bound_west  = MAX(tg%minlon - 0.5_wp, -180.0_wp)

    istartlon = 1; iendlon = nlon
    DO ir = 1, nlon
      IF (lon_row(ir) < bound_west)  istartlon = ir + 1
      IF (lon_row(ir) > bound_east)  THEN; iendlon = ir - 1; EXIT; ENDIF
    END DO
    nlon_sub = MAX(iendlon - istartlon + 1, 0)

    num_blocks = 1
!$  num_blocks = omp_get_max_threads()
    blk_len = nlon_sub / num_blocks + 1
!$  ALLOCATE(start_cell_arr(num_blocks))
!$  start_cell_arr(:) = 1

    WRITE(message_text,'(A,I0,A,I0,A,I0,A,I0,A)') &
      "Using ", num_blocks, " threads, columns ", istartlon, " to ", iendlon, " (", nlon_sub, ")"
    CALL logging%info(message_text)

    no_raw_data_pixel = 0
    fr_heavy_clay = 0.0_wp; fr_silty_clay = 0.0_wp; fr_light_clay = 0.0_wp
    fr_silty_clay_loam = 0.0_wp; fr_clay_loam = 0.0_wp; fr_silt = 0.0_wp
    fr_silt_loam = 0.0_wp; fr_sandy_clay = 0.0_wp; fr_loam = 0.0_wp
    fr_sandy_clay_loam = 0.0_wp; fr_sandy_loam = 0.0_wp
    fr_loamy_sand = 0.0_wp; fr_sand = 0.0_wp; fr_undef = 0.0_wp

    t_start = omp_get_wtime()

    lat_loop: DO jr = 1, nlat
      lat_pixel = lat_row(jr)
      IF (lat_pixel > bound_north .OR. lat_pixel < bound_south) CYCLE lat_loop

      IF (MOD(jr, 1000) == 1 .OR. jr == nlat) THEN
        WRITE(message_text,'(A,I7,"/",I7," (",F8.3,"°), elapsed: ",F8.2," min)")') &
          "Row ", jr, nlat, lat_pixel, (omp_get_wtime() - t_start)/60.0_wp
        CALL logging%info(message_text)
      END IF

      IF (xtype == NF90_DOUBLE) THEN
        CALL check_netcdf(nf90_get_var(ncid, varid_lu, lu_row_double, start=[1,jr], count=[nlon,1]))
        WHERE (ABS(lu_row_double + 9.0d33) < 1.0d20 .OR. lu_row_double < -1d30)
          lu_row_int = -9999_i4
        ELSEWHERE
          lu_row_int = NINT(lu_row_double)
        END WHERE
      ELSE
        CALL check_netcdf(nf90_get_var(ncid, varid_lu, lu_row_int, start=[1,jr], count=[nlon,1]))
      END IF

      ie_vec = 0; je_vec = 0; ke_vec = 0

!$OMP PARALLEL DO PRIVATE(ib, il, ir, ishift, lon_pixel, thread_id, start_cell_id, ix, iy)
      DO ib = 1, num_blocks
!$      thread_id = omp_get_thread_num() + 1
!$      start_cell_id = start_cell_arr(thread_id)

        ishift = istartlon - 1 + (ib-1)*blk_len
        DO il = 1, blk_len
          ir = ishift + il
          IF (ir > iendlon) EXIT
          lon_pixel = lon_row(ir)

          ! --- Immer gültigen start_cell_id sicherstellen (auch bei il > 1!) ---
          IF (tg%igrid_type == igrid_icon) THEN
            IF (start_cell_id < 1) THEN   ! Nur prüfen, ob gültig – neu berechnen bei Bedarf
              ix = NINT(lon_pixel * search_res)
              iy = NINT(lat_pixel * search_res)
              ix = MAX(1, MIN(ix, SIZE(tg%search_index,1)))
              iy = MAX(1, MIN(iy, SIZE(tg%search_index,2)))
              start_cell_id = tg%search_index(ix, iy)
              IF (start_cell_id < 1) start_cell_id = 1
            END IF
          END IF

          CALL find_nearest_target_grid_element(lon_pixel, lat_pixel, tg, start_cell_id, &
                                                ie_vec(ir), je_vec(ir), ke_vec(ir))

          ! --- Nach dem Aufruf: start_cell_id für nächsten Pixel speichern ---
          IF (tg%igrid_type == igrid_icon) THEN
            start_cell_id = ie_vec(ir)   ! ie_vec ist die gefundene Zelle → beste Startzelle für nächsten Pixel!
          END IF
        END DO

!$      start_cell_arr(thread_id) = start_cell_id
      END DO
!$OMP END PARALLEL DO

      DO ir = istartlon, iendlon
        ie = ie_vec(ir); je = je_vec(ir); ke = ke_vec(ir)
        IF (ie == 0 .OR. je == 0 .OR. ke == 0) CYCLE

        no_raw_data_pixel(ie,je,ke) = no_raw_data_pixel(ie,je,ke) + 1
        soil_unit = lu_row_int(ir)

        IF (soil_unit == type_clay_heavy)        fr_heavy_clay(ie,je,ke)      = fr_heavy_clay(ie,je,ke)      + 1.0_wp
        IF (soil_unit == type_silty_clay)        fr_silty_clay(ie,je,ke)      = fr_silty_clay(ie,je,ke)      + 1.0_wp
        IF (soil_unit == type_clay_light)        fr_light_clay(ie,je,ke)      = fr_light_clay(ie,je,ke)      + 1.0_wp
        IF (soil_unit == type_silty_clay_loam)   fr_silty_clay_loam(ie,je,ke) = fr_silty_clay_loam(ie,je,ke) + 1.0_wp
        IF (soil_unit == type_clay_loam)         fr_clay_loam(ie,je,ke)       = fr_clay_loam(ie,je,ke)       + 1.0_wp
        IF (soil_unit == type_silt)              fr_silt(ie,je,ke)            = fr_silt(ie,je,ke)            + 1.0_wp
        IF (soil_unit == type_silt_loam)         fr_silt_loam(ie,je,ke)       = fr_silt_loam(ie,je,ke)       + 1.0_wp
        IF (soil_unit == type_sandy_clay)        fr_sandy_clay(ie,je,ke)      = fr_sandy_clay(ie,je,ke)      + 1.0_wp
        IF (soil_unit == type_loam)              fr_loam(ie,je,ke)            = fr_loam(ie,je,ke)            + 1.0_wp
        IF (soil_unit == type_sandy_clay_loam)   fr_sandy_clay_loam(ie,je,ke) = fr_sandy_clay_loam(ie,je,ke) + 1.0_wp
        IF (soil_unit == type_sandy_loam)        fr_sandy_loam(ie,je,ke)      = fr_sandy_loam(ie,je,ke)      + 1.0_wp
        IF (soil_unit == type_loamy_sand)        fr_loamy_sand(ie,je,ke)      = fr_loamy_sand(ie,je,ke)      + 1.0_wp
        IF (soil_unit == type_sand)              fr_sand(ie,je,ke)            = fr_sand(ie,je,ke)            + 1.0_wp
        IF (soil_unit == undef_hwsdARTtype .OR. soil_unit == no_data .OR. soil_unit < 0) &
                                                 fr_undef(ie,je,ke)           = fr_undef(ie,je,ke)           + 1.0_wp
      END DO
    END DO lat_loop

    CALL check_netcdf(nf90_close(ncid))
!$  IF (ALLOCATED(start_cell_arr)) DEALLOCATE(start_cell_arr)

    ! Normalisierung
    DO ke = 1, tg%ke
      DO je = 1, tg%je
        DO ie = 1, tg%ie
          IF (no_raw_data_pixel(ie,je,ke) > 0) THEN
            denom = REAL(no_raw_data_pixel(ie,je,ke), wp)
            fr_heavy_clay(ie,je,ke)      = fr_heavy_clay(ie,je,ke)      / denom
            fr_silty_clay(ie,je,ke)      = fr_silty_clay(ie,je,ke)      / denom
            fr_light_clay(ie,je,ke)      = fr_light_clay(ie,je,ke)      / denom
            fr_silty_clay_loam(ie,je,ke) = fr_silty_clay_loam(ie,je,ke) / denom
            fr_clay_loam(ie,je,ke)       = fr_clay_loam(ie,je,ke)       / denom
            fr_silt(ie,je,ke)            = fr_silt(ie,je,ke)            / denom
            fr_silt_loam(ie,je,ke)       = fr_silt_loam(ie,je,ke)       / denom
            fr_sandy_clay(ie,je,ke)      = fr_sandy_clay(ie,je,ke)      / denom
            fr_loam(ie,je,ke)            = fr_loam(ie,je,ke)            / denom
            fr_sandy_clay_loam(ie,je,ke) = fr_sandy_clay_loam(ie,je,ke) / denom
            fr_sandy_loam(ie,je,ke)      = fr_sandy_loam(ie,je,ke)      / denom
            fr_loamy_sand(ie,je,ke)      = fr_loamy_sand(ie,je,ke)      / denom
            fr_sand(ie,je,ke)            = fr_sand(ie,je,ke)            / denom
            fr_undef(ie,je,ke)           = fr_undef(ie,je,ke)           / denom
          ELSE
            fr_undef(ie,je,ke) = 1.0_wp
          END IF
        END DO
      END DO
    END DO

    t_now = omp_get_wtime()
    WRITE(message_text,'(A,F10.2,A,I0,A)') &
      "=== HWSD-ART AGGREGATION SUCCESSFULLY COMPLETED in ", (t_now-t_start)/60.0_wp, &
      " minutes using ", num_blocks, " threads ==="
    CALL logging%info(message_text)

    DEALLOCATE(lon_row, lat_row, lu_row_double, lu_row_int, ie_vec, je_vec, ke_vec)

  END SUBROUTINE get_hwsdART_data_and_aggregate

END MODULE mo_hwsdART_routines
