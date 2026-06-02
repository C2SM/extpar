PROGRAM extpar_hiressoil_to_buffer

USE mo_logging
USE info_extpar,       ONLY: info_print
USE mo_io_units,       ONLY: filename_max
USE mo_kind,           ONLY: wp, i4
USE mo_target_grid_data, ONLY: tg
USE mo_grid_structures, ONLY: igrid_icon, igrid_cosmo
USE mo_cosmo_grid,     ONLY: cosmo_grid
USE mo_icon_grid_data, ONLY: icon_grid

USE mo_hiressoil_routines, ONLY: read_hiressoil_namelist_entry, &
                                 read_namelists_extpar_hiressoil, &
                                 split_hiressoil_namelist, &
                                 read_hiressoil_namelist_by_index, &
                                 read_hiressoil_namelist_all,   &
                                 get_hiressoil_data_and_aggregate
!USE mo_hiressoil_data,        ONLY: define_hiressoil_type
USE mo_hiressoil_tg_fields,   ONLY: allocate_hiressoil_target_fields, &
                                    deallocate_hiressoil_target_fields
USE mo_hiressoil_output_nc,   ONLY: write_netcdf_hiressoil_icon_grid, &
                                    write_netcdf_hiressoil_cosmo_grid
USE mo_target_grid_routines, ONLY: init_target_grid

IMPLICIT NONE

!========================================================================
INTEGER(KIND=i4), PARAMETER :: max_files = 50

CHARACTER(len=filename_max) :: raw_data_hiressoil_path
CHARACTER(len=filename_max) :: raw_data_hiressoil_filename
CHARACTER(len=filename_max) :: hiressoil_output_file
CHARACTER(len=filename_max) :: namelist_single
CHARACTER(len=filename_max) :: input_file_full
REAL(KIND=wp)    :: undefined = -999.0_wp
INTEGER(KIND=i4) :: i
INTEGER(KIND=i4) :: n_entries
CHARACTER(len=8) :: tmp_str   ! für Zahl → String Umwandlung

!========================================================================
CALL initialize_logging("extpar_hiressoil_to_buffer.log")
CALL info_print()

CALL logging%info('')
CALL logging%info('============= start hiressoil_to_buffer =============')
CALL logging%info('')

CALL init_target_grid('INPUT_grid_org')

CALL read_hiressoil_namelist_all('INPUT_hiressoil', n_entries)

!=== Namelist aufteilen ===
CALL split_hiressoil_namelist('INPUT_hiressoil', n_entries)


IF (n_entries == 0) THEN
  CALL logging%error('No entries found in INPUT_hiressoil')
  STOP 1
END IF


!========================================================================
! Schleife: Mehrere Einträge im selben Namelist-Block lesen
!========================================================================
DO i = 1, n_entries
 WRITE(namelist_single, '(A,I0)') 'INPUT_hiressoil_', i
  CALL read_namelists_extpar_hiressoil(namelist_single,       &
                                     raw_data_hiressoil_path,  &
                                     raw_data_hiressoil_filename, &
                                     hiressoil_output_file)


   
  input_file_full = TRIM(raw_data_hiressoil_path) // TRIM(raw_data_hiressoil_filename)

  ! Zahl in String umwandeln
  WRITE(tmp_str,'(I0)') i
  
  CALL logging%info('')
  CALL logging%info('--- Datei '//TRIM(ADJUSTL(tmp_str))//' ---')
  CALL logging%info('Input : '//TRIM(input_file_full))
  CALL logging%info('Output: '//TRIM(hiressoil_output_file))

  !Allocate Fields
  CALL allocate_hiressoil_target_fields(tg)
  ! Verarbeitung
  CALL get_hiressoil_data_and_aggregate(input_file_full, tg)

  ! Ausgabe
  SELECT CASE(tg%igrid_type)
    CASE(igrid_icon)
      CALL write_netcdf_hiressoil_icon_grid(TRIM(hiressoil_output_file), &
                                            icon_grid, tg, undefined)
    CASE(igrid_cosmo)
      CALL write_netcdf_hiressoil_cosmo_grid(TRIM(hiressoil_output_file), &
                                             cosmo_grid, tg, undefined)
  END SELECT

  CALL logging%info('→ Erfolgreich geschrieben: '//TRIM(hiressoil_output_file))

  ! === Deallocate nach der Verarbeitung ===
  CALL deallocate_hiressoil_target_fields()

END DO

WRITE(tmp_str,'(I0)') n_entries
CALL logging%info('')
CALL logging%info('============= hiressoil_to_buffer beendet ('// &
                  TRIM(ADJUSTL(tmp_str))//' Dateien) =============')

END PROGRAM extpar_hiressoil_to_buffer
