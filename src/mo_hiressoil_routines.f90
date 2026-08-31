!================================================================================
! mo_hiressoil_routines.f90
! FINAL – läuft mit 1 bis 64+ Threads, GCC 13.2 + -fcheck=all, kein Absturz mehr!
!================================================================================
MODULE mo_hiressoil_routines
  USE mo_logging
  USE mo_kind, ONLY: wp, i4
  USE netcdf, ONLY: nf90_open, nf90_close, nf90_inquire, nf90_inquire_dimension, &
                    nf90_inquire_variable, nf90_inq_varid, nf90_get_var, &
                    NF90_NOWRITE, NF90_DOUBLE, NF90_INT
  USE mo_io_utilities, ONLY: check_netcdf
  USE mo_io_units, ONLY: filename_max
  USE mo_GRID_structures, ONLY: target_grid_def, igrid_icon
  USE mo_utilities_extpar, ONLY: free_un
!  USE mo_hiressoil_data, ONLY: undef_hiressoiltype, no_data, &
!    type_clay_heavy, type_silty_clay, type_clay_light, type_silty_clay_loam, &
!    type_clay_loam, type_silt, type_silt_loam, type_sandy_clay, type_loam, &
!    type_sandy_clay_loam, type_sandy_loam, type_loamy_sand, type_sand
  USE mo_hiressoil_tg_fields, ONLY: hiressoil_mean, hiressoil_min, &
                                   hiressoil_max, hiressoil_var
  
  USE mo_target_grid_data, ONLY: no_raw_data_pixel, search_res
  USE mo_search_target_grid, ONLY: find_nearest_target_grid_element
  USE omp_lib
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: get_hiressoil_data_and_aggregate
  PUBLIC :: get_dimension_hiressoil_data          ! Korrekt geschrieben
  PUBLIC :: read_namelists_extpar_hiressoil       ! Korrekt geschrieben
  PUBLIC :: read_hiressoil_namelist_entry         ! Korrekt geschrieben
  PUBLIC :: read_hiressoil_namelist_all         ! Korrekt geschrieben
  PUBLIC :: read_hiressoil_namelist_by_index
  PUBLIC :: split_hiressoil_namelist
CONTAINS


SUBROUTINE read_hiressoil_namelist_all(namelist_file,n_entries)
  USE mo_kind, ONLY: i4
  USE mo_io_units, ONLY: filename_max
  USE mo_logging
  IMPLICIT NONE

  CHARACTER(len=*), INTENT(IN)  :: namelist_file
  INTEGER(KIND=i4),  INTENT(OUT) :: n_entries

  INTEGER(KIND=i4) :: unit_nml, ios
  CHARACTER(len=filename_max) :: raw_data_hiressoil_path
  CHARACTER(len=filename_max) :: raw_data_hiressoil_filename
  CHARACTER(len=filename_max) :: raw_data_hiressoil_varname
  CHARACTER(len=filename_max) :: hiressoil_output_file

  NAMELIST /hiressoil_nml/ raw_data_hiressoil_path, &
                           raw_data_hiressoil_filename, &
                           raw_data_hiressoil_varname, &
                           hiressoil_output_file

  n_entries = 0

  OPEN(NEWUNIT=unit_nml, FILE=TRIM(namelist_file), STATUS='OLD', &
       ACTION='READ', IOSTAT=ios)
  IF (ios /= 0) THEN
    CALL logging%error('Cannot open namelist: '//TRIM(namelist_file))
    STOP 1
  END IF

  DO
    raw_data_hiressoil_path     = ''
    raw_data_hiressoil_filename = ''
    raw_data_hiressoil_varname  = ''
    hiressoil_output_file       = ''

    READ(unit_nml, NML=hiressoil_nml, IOSTAT=ios)
    IF (ios < 0) EXIT
    IF (ios > 0) CYCLE

    IF (LEN_TRIM(hiressoil_output_file) > 0) n_entries = n_entries + 1
  END DO

  CLOSE(unit_nml)

  WRITE(message_text,'(A,I0,A)') &
    'Namelist '//TRIM(namelist_file)//' contains ', n_entries, ' hiressoil entries.'
  CALL logging%info(message_text)

END SUBROUTINE read_hiressoil_namelist_all

!==============================================================================
SUBROUTINE split_hiressoil_namelist(input_namelist)
!==============================================================================
  USE mo_kind,      ONLY: i4
  USE mo_io_units,  ONLY: filename_max
  USE mo_logging

  IMPLICIT NONE

  CHARACTER(len=*), INTENT(IN)  :: input_namelist
  INTEGER(KIND=i4)  :: n_entries

  INTEGER(KIND=i4) :: unit_in, unit_out, ios, entry_id, ipos
  CHARACTER(len=512) :: line, trimmed
  CHARACTER(len=filename_max) :: out_filename
  CHARACTER(len=filename_max) :: path_line, file_line, var_line, output_line
  CHARACTER(len=filename_max) :: tmp, suffix
  CHARACTER(len=16) :: tmp_str
  
  n_entries = 0
  entry_id  = 0

  OPEN(NEWUNIT=unit_in, FILE=TRIM(input_namelist), STATUS='OLD', ACTION='READ', IOSTAT=ios)
  IF (ios /= 0) THEN
    CALL logging%error('Cannot open INPUT_hiressoil')
    STOP 1
  END IF

  DO
    path_line   = ''
    file_line   = ''
    var_line   = ''
    output_line = ''

    ! Suche nach einem vollständigen Eintrag (3 Zeilen)
    DO
      READ(unit_in, '(A)', IOSTAT=ios) line
      IF (ios /= 0) EXIT

      trimmed = ADJUSTL(line)

      IF (INDEX(trimmed, 'raw_data_hiressoil_path') > 0) THEN
        path_line = TRIM(line)
      ELSE IF (INDEX(trimmed, 'raw_data_hiressoil_filename') > 0) THEN
        file_line = TRIM(line)
      ELSE IF (INDEX(trimmed, 'raw_data_hiressoil_varname') > 0) THEN
        var_line = TRIM(line)
     ELSE IF (INDEX(trimmed, 'hiressoil_output_file') > 0) THEN
        output_line = TRIM(line)
      END IF

      ! Wenn wir alle 4 Zeilen haben → neuen Eintrag schreiben
      IF (LEN_TRIM(path_line) > 0 .AND. LEN_TRIM(file_line) > 0 .AND. LEN_TRIM(var_line) > 0 .AND. LEN_TRIM(output_line) > 0) THEN
        EXIT
      END IF
    END DO

    IF (LEN_TRIM(output_line) == 0) EXIT   ! Kein weiterer Eintrag

    entry_id = entry_id + 1
    n_entries = n_entries + 1

    ! ---------------------------------------------------------------
    ! Suffix aus hiressoil_output_file extrahieren
    ! Beispiel: 'KSAT_extpar_ICON_hiressoil.nc' → 'KSAT'
    ! ---------------------------------------------------------------
    suffix = ''
    tmp    = ADJUSTL(output_line)

    ipos = INDEX(tmp, '=')
    IF (ipos > 0) THEN
      tmp = ADJUSTL(tmp(ipos+1:))

      IF (tmp(1:1) == "'" .OR. tmp(1:1) == '"') tmp = tmp(2:)

      ipos = INDEX(tmp, "'")
      IF (ipos > 0) tmp = tmp(1:ipos-1)

      ipos = INDEX(tmp, '"')
      IF (ipos > 0) tmp = tmp(1:ipos-1)

      ipos = INDEX(tmp, ',')
      IF (ipos > 0) tmp = tmp(1:ipos-1)

      tmp = ADJUSTL(TRIM(tmp))

      ipos = INDEX(tmp, '_')
      IF (ipos > 1) THEN
        suffix = tmp(1:ipos-1)
      ELSE
        ipos = INDEX(tmp, '.')
        IF (ipos > 1) THEN
          suffix = tmp(1:ipos-1)
        ELSE
          suffix = TRIM(tmp)
        END IF
      END IF
    END IF

    IF (LEN_TRIM(suffix) > 0) THEN
      WRITE(out_filename, '(A,A)') 'INPUT_HHS_', TRIM(suffix)
    ELSE
      WRITE(out_filename, '(A,I0)') 'INPUT_hiressoil_', entry_id
    END IF

    ! --- Wichtig: OPEN mit IOSTAT prüfen ---
    OPEN(NEWUNIT=unit_out, FILE=TRIM(out_filename), STATUS='REPLACE', &
         ACTION='WRITE', IOSTAT=ios)
    IF (ios /= 0) THEN
      CALL logging%error('Cannot create '//TRIM(out_filename)// &
                         '  (IOSTAT='//TRIM(ADJUSTL(tmp_str))//')')
      ! optional: Fallback-Name
      WRITE(out_filename, '(A,I0)') 'INPUT_hiressoil_', entry_id
      OPEN(NEWUNIT=unit_out, FILE=TRIM(out_filename), STATUS='REPLACE', &
           ACTION='WRITE', IOSTAT=ios)
      IF (ios /= 0) THEN
        CALL logging%error('Also failed to create fallback file')
        STOP 1
      END IF
    END IF

    WRITE(unit_out, '(A)') '&hiressoil_nml'
    WRITE(unit_out, '(A)') TRIM(path_line)
    WRITE(unit_out, '(A)') TRIM(file_line)
    WRITE(unit_out, '(A)') TRIM(var_line)
    WRITE(unit_out, '(A)') TRIM(output_line)
    WRITE(unit_out, '(A)') '/'
    CLOSE(unit_out)

    CALL logging%info('Created: '//TRIM(out_filename))
 END DO
 CLOSE(unit_in)

!  CALL logging%info('Successfully split into '//TRIM(ADJUSTL(int2string(n_entries)))// &
!                    ' single-entry namelist files.')

END SUBROUTINE split_hiressoil_namelist



!==============================================================================
SUBROUTINE read_hiressoil_namelist_by_index(namelist_file, idx, &
                                            raw_data_path,      &
                                            raw_data_filename,  &
                                            raw_data_varname , &
                                            output_file)
!==============================================================================
  USE mo_kind,      ONLY: i4
  USE mo_io_units,  ONLY: filename_max
  USE mo_logging

  IMPLICIT NONE

  CHARACTER(len=*), INTENT(IN)  :: namelist_file
  INTEGER(KIND=i4), INTENT(IN)  :: idx
  CHARACTER(len=*), INTENT(OUT) :: raw_data_path
  CHARACTER(len=*), INTENT(OUT) :: raw_data_filename
  CHARACTER(len=*), INTENT(OUT) :: raw_data_varname
  CHARACTER(len=*), INTENT(OUT) :: output_file

  
  INTEGER(KIND=i4) :: unit_nml, ios, current_idx
 
  CHARACTER(len=filename_max) :: raw_data_hiressoil_path     = './'
  CHARACTER(len=filename_max) :: raw_data_hiressoil_filename = ''
  CHARACTER(len=filename_max) :: raw_data_hiressoil_varname = ''
  CHARACTER(len=filename_max) :: hiressoil_output_file       = ''
  
  NAMELIST /hiressoil_nml/ raw_data_hiressoil_path,     &
                           raw_data_hiressoil_filename,  &
                           raw_data_hiressoil_varname , &
                           hiressoil_output_file

  raw_data_path     = ''
  raw_data_filename = ''
  raw_data_varname = ''
  output_file       = ''

  OPEN(NEWUNIT=unit_nml, FILE=TRIM(namelist_file), STATUS='OLD', ACTION='READ', IOSTAT=ios)
  IF (ios /= 0) THEN
    CALL logging%error('Cannot open namelist')
    STOP 1
  END IF

  current_idx = 0

  DO
    raw_data_hiressoil_path     = './'
    raw_data_hiressoil_filename = ''
    raw_data_hiressoil_varname = ''
    hiressoil_output_file       = ''

    READ(unit_nml, NML=hiressoil_nml, IOSTAT=ios)
    IF (ios < 0) EXIT
    IF (ios > 0) CYCLE

    current_idx = current_idx + 1

    IF (current_idx == idx) THEN
      raw_data_path     = TRIM(ADJUSTL(raw_data_hiressoil_path))
      raw_data_filename = TRIM(ADJUSTL(raw_data_hiressoil_filename))
      raw_data_varname  = TRIM(ADJUSTL(raw_data_hiressoil_varname))
      output_file       = TRIM(ADJUSTL(hiressoil_output_file))
      CLOSE(unit_nml)
      RETURN
    END IF
  END DO

  CLOSE(unit_nml)

END SUBROUTINE read_hiressoil_namelist_by_index

!==============================================================================
SUBROUTINE read_hiressoil_namelist_entry(namelist_file,     &
                                         raw_data_path,     &
                                         raw_data_filename, &
                                         raw_data_varname, &
                                         output_file,       &
                                         entry_number)
!==============================================================================
! Liest den n-ten Eintrag (1-basiert)
!==============================================================================
  USE mo_kind,      ONLY: i4
  USE mo_io_units,  ONLY: filename_max
  USE mo_logging

  IMPLICIT NONE

  CHARACTER(len=*), INTENT(IN)    :: namelist_file
  CHARACTER(len=*), INTENT(OUT)   :: raw_data_path
  CHARACTER(len=*), INTENT(OUT)   :: raw_data_filename
  CHARACTER(len=*), INTENT(OUT)   :: raw_data_varname
  CHARACTER(len=*), INTENT(OUT)   :: output_file
  INTEGER(KIND=i4), INTENT(IN)    :: entry_number

  INTEGER(KIND=i4) :: ierr, unit_nml, current_entry
  INTEGER(KIND=i4) :: ios

  CHARACTER(len=filename_max) :: raw_data_hiressoil_path     = './'
  CHARACTER(len=filename_max) :: raw_data_hiressoil_filename = ''
  CHARACTER(len=filename_max) :: raw_data_hiressoil_varname           = ''
  CHARACTER(len=filename_max) :: hiressoil_output_file       = ''

  NAMELIST /hiressoil_nml/ raw_data_hiressoil_path,     &
                           raw_data_hiressoil_filename,  &
                           raw_data_hiressoil_varname , &
                           hiressoil_output_file

  raw_data_path     = ''
  raw_data_filename = ''
  raw_data_varname = ''
  output_file       = ''

  OPEN(NEWUNIT=unit_nml, FILE=TRIM(namelist_file), STATUS='OLD', ACTION='READ', IOSTAT=ierr)
  IF (ierr /= 0) THEN
    CALL logging%error('Cannot open namelist')
    STOP 1
  END IF

  current_entry = 0

  DO
    raw_data_hiressoil_path     = './'
    raw_data_hiressoil_filename = ''
    raw_data_hiressoil_varname           = ''
    hiressoil_output_file       = ''

    READ(unit_nml, NML=hiressoil_nml, IOSTAT=ios)

    IF (ios < 0) EXIT          ! Ende der Datei
    IF (ios > 0) THEN
      CALL logging%warning('Skipping malformed namelist entry')
      CYCLE
    END IF

    current_entry = current_entry + 1

    IF (current_entry == entry_number) THEN
      raw_data_path     = TRIM(ADJUSTL(raw_data_hiressoil_path))
      raw_data_filename = TRIM(ADJUSTL(raw_data_hiressoil_filename))
      raw_data_varname  = TRIM(ADJUSTL(raw_data_hiressoil_varname))
      output_file       = TRIM(ADJUSTL(hiressoil_output_file))
      CLOSE(unit_nml)
      RETURN
    END IF
  END DO

  CLOSE(unit_nml)

END SUBROUTINE read_hiressoil_namelist_entry
  

  
  SUBROUTINE read_namelists_extpar_hiressoil(namelist_file, raw_data_hiressoil_path, &
                                           raw_data_hiressoil_filename, raw_data_hiressoil_varname, hiressoil_output_file)
    CHARACTER(len=filename_max), INTENT(IN)  :: namelist_file
    CHARACTER(len=filename_max), INTENT(OUT) :: raw_data_hiressoil_path
    CHARACTER(len=filename_max), INTENT(OUT) :: raw_data_hiressoil_filename
    CHARACTER(len=filename_max), INTENT(OUT) :: raw_data_hiressoil_varname
    CHARACTER(len=filename_max), INTENT(OUT) :: hiressoil_output_file
    NAMELIST /hiressoil_nml/ raw_data_hiressoil_path, raw_data_hiressoil_filename, raw_data_hiressoil_varname, hiressoil_output_file
    INTEGER :: nuin, ierr
    nuin = free_un()
    OPEN(nuin, FILE=TRIM(namelist_file), IOSTAT=ierr)
    IF (ierr /= 0) CALL logging%error('Cannot open namelist file', __FILE__, __LINE__)
    READ(nuin, NML=hiressoil_nml, IOSTAT=ierr)
    IF (ierr /= 0) CALL logging%error('Error reading hiressoil_nml', __FILE__, __LINE__)
    CLOSE(nuin)
  END SUBROUTINE read_namelists_extpar_hiressoil


  SUBROUTINE get_dimension_hiressoil_data(path_hiressoil_file, nlon_hiressoil, nlat_hiressoil)
    CHARACTER(len=*), INTENT(in)  :: path_hiressoil_file
    INTEGER(i4),      INTENT(out) :: nlon_hiressoil, nlat_hiressoil
    INTEGER :: ncid, ndims, dimid, length
    CHARACTER(len=80) :: dimname
    CALL check_netcdf(nf90_open(TRIM(path_hiressoil_file), NF90_NOWRITE, ncid))
    CALL check_netcdf(nf90_inquire(ncid, nDimensions=ndims))
    nlon_hiressoil = 0; nlat_hiressoil = 0
    DO dimid = 1, ndims
      CALL check_netcdf(nf90_inquire_dimension(ncid, dimid, dimname, length))
      IF (TRIM(dimname) == 'lon'  .OR. TRIM(dimname) == 'longitude') nlon_hiressoil = length
      IF (TRIM(dimname) == 'lat'  .OR. TRIM(dimname) == 'latitude')  nlat_hiressoil = length
    END DO
    CALL check_netcdf(nf90_close(ncid))
  END SUBROUTINE get_dimension_hiressoil_data


  SUBROUTINE get_hiressoil_data_and_aggregate(path_hiressoil_file, varname_hiressoil_file, tg)
    CHARACTER(len=*),      INTENT(in) :: path_hiressoil_file
    CHARACTER(len=*),      INTENT(in) :: varname_hiressoil_file
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


    REAL(wp), ALLOCATABLE :: sum_soil   (:,:,:)
    REAL(wp), ALLOCATABLE :: min_soil   (:,:,:)
    REAL(wp), ALLOCATABLE :: max_soil   (:,:,:)
    REAL(wp), ALLOCATABLE :: sumsq_soil (:,:,:)   ! für Varianz-Berechnung
    REAL(wp), ALLOCATABLE :: mean_soil   (:,:,:)
    REAL(wp), ALLOCATABLE :: var_soil   (:,:,:)
    REAL(wp) :: mean_val, soil_unit_real 
    
    !$  INTEGER :: thread_id

    REAL(wp) :: bound_north, bound_south, bound_west, bound_east
    REAL(wp) :: denom


    
    CALL get_dimension_hiressoil_data(path_hiressoil_file, nlon, nlat)

    WRITE(message_text,'(A,I0,A,I0,A)') &
      "=== HIRESSOIL Aggregation: ", nlon, " x ", nlat, " → target grid (FULLY PARALLEL) ==="
    CALL logging%info(message_text)

    ALLOCATE(lon_row(nlon), lat_row(nlat), lu_row_double(nlon), lu_row_int(nlon))
    ALLOCATE(ie_vec(nlon), je_vec(nlon), ke_vec(nlon))

    ALLOCATE ( sum_soil   (tg%ie,tg%je,tg%ke))
    ALLOCATE ( min_soil(tg%ie,tg%je,tg%ke))
    ALLOCATE ( max_soil(tg%ie,tg%je,tg%ke))
    ALLOCATE ( sumsq_soil(tg%ie,tg%je,tg%ke))
    ALLOCATE ( mean_soil (tg%ie,tg%je,tg%ke))
    ALLOCATE ( var_soil  (tg%ie,tg%je,tg%ke))
    
    CALL check_netcdf(nf90_open(TRIM(path_hiressoil_file), NF90_NOWRITE, ncid))
    CALL check_netcdf(nf90_inq_varid(ncid, "lon", varid_lon))
    CALL check_netcdf(nf90_get_var(ncid, varid_lon, lon_row))
    CALL check_netcdf(nf90_inq_varid(ncid, "lat", varid_lat))
    CALL check_netcdf(nf90_get_var(ncid, varid_lat, lat_row))
    CALL check_netcdf(nf90_inq_varid(ncid,varname_hiressoil_file, varid_lu))
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
    hiressoil_mean = 0.0_wp; hiressoil_min = 0.0_wp; hiressoil_max = 0.0_wp
    hiressoil_var = 0.0_wp


    sum_soil   = 0.0_wp
    min_soil   = HUGE(1.0_wp)          ! sehr großer Wert
    max_soil   = -HUGE(1.0_wp)         ! sehr kleiner Wert
    sumsq_soil = 0.0_wp
    
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


        soil_unit = lu_row_int(ir)
        soil_unit_real = REAL(lu_row_int(ir), wp) / 10000._wp
        
        IF (lu_row_int(ir) >= 0) no_raw_data_pixel(ie,je,ke) = no_raw_data_pixel(ie,je,ke) + 1
        
        IF (lu_row_int(ir) >= 0) THEN
        sum_soil(ie,je,ke) = sum_soil(ie,je,ke) + soil_unit_real
        sumsq_soil(ie,je,ke) = sumsq_soil(ie,je,ke) + soil_unit_real**2
        
        IF (soil_unit_real < min_soil(ie,je,ke)) THEN
           min_soil(ie,je,ke) = soil_unit_real
        END IF
        
        IF (soil_unit_real > max_soil(ie,je,ke)) THEN
           max_soil(ie,je,ke) = soil_unit_real
        END IF
        END IF
        
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

             mean_val = sum_soil(ie,je,ke) / denom
             ! Mittelwert
             mean_soil(ie,je,ke) = mean_val   ! ← neues Feld, das du definieren musst

             ! Minimum & Maximum (falls kein Wert → bleibt HUGE / -HUGE)
             ! (optional: auf -9999 oder ähnliches setzen, wenn gewünscht)

             ! Varianz (Stichprobenvarianz – ddof=1)
             IF (no_raw_data_pixel(ie,je,ke) >= 2) THEN
                var_soil(ie,je,ke) = (sumsq_soil(ie,je,ke) - &
                               (sum_soil(ie,je,ke)**2 / denom) ) / (denom - 1.0_wp)
             ELSE
                var_soil(ie,je,ke) = 0.0_wp
             END IF


          ELSE

             ! Keine gültigen Rohdaten → Fill-Value setzen
             mean_soil(ie,je,ke) = -999._wp
             min_soil (ie,je,ke) = -999._wp
             max_soil (ie,je,ke) = -999._wp
             var_soil (ie,je,ke) = -999._wp
          END IF

            hiressoil_mean(ie,je,ke)     = mean_soil(ie,je,ke)
            hiressoil_min(ie,je,ke)      = min_soil (ie,je,ke)
            hiressoil_max(ie,je,ke)      = max_soil (ie,je,ke)
            hiressoil_var(ie,je,ke)      = var_soil (ie,je,ke)
          
       END DO
      END DO
    END DO

    t_now = omp_get_wtime()
    WRITE(message_text,'(A,F10.2,A,I0,A)') &
      "=== HIRESSOIL AGGREGATION SUCCESSFULLY COMPLETED in ", (t_now-t_start)/60.0_wp, &
      " minutes using ", num_blocks, " threads ==="
    CALL logging%info(message_text)

    DEALLOCATE(lon_row, lat_row, lu_row_double, lu_row_int, ie_vec, je_vec, ke_vec)
    DEALLOCATE ( sum_soil)
    DEALLOCATE ( min_soil)
    DEALLOCATE ( max_soil)
    DEALLOCATE ( sumsq_soil)
    DEALLOCATE ( mean_soil) 
    DEALLOCATE ( var_soil ) 
  END SUBROUTINE get_hiressoil_data_and_aggregate

END MODULE mo_hiressoil_routines
