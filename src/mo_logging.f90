MODULE mo_logging

  USE mo_io_units,          ONLY: filename_max
#ifdef NAGFOR
  USE f90_unix, ONLY: exit
#endif

  IMPLICIT NONE


  !Integer for debugging levels
  INTEGER, PARAMETER, PUBLIC :: verbose = 2 ! verbosity of extpar, add to namelist
  INTEGER, PARAMETER, PUBLIC :: idbg_low  = 1 ! low debug output
  INTEGER, PARAMETER, PUBLIC :: idbg_high = 2 ! high debug output
  INTEGER, PARAMETER :: closed = -1
  
  TYPE, PUBLIC :: logger
    CHARACTER(len=:), ALLOCATABLE :: logfile    
    INTEGER                       :: fileunit
  CONTAINS
    PROCEDURE :: message            => logger_message
    PROCEDURE :: info               => logger_info
    PROCEDURE :: warning            => logger_warning
    PROCEDURE :: error              => logger_error
  END TYPE logger

  TYPE(logger), PUBLIC :: logging

  PUBLIC :: initialize_logging

  CHARACTER(len=filename_max) :: message_text = ""

  PUBLIC :: message_text
  
CONTAINS

  FUNCTION constructor(logfile) RESULT(this)
    TYPE(logger) :: this
    CHARACTER(len=*), INTENT(in)  :: logfile
    INTEGER :: flag
      this%logfile = logfile
      this%fileunit= free_unit_number()
      OPEN(newunit=this%fileunit,file=this%logfile,action='write',asynchronous='yes',iostat=flag,status='replace')
  END FUNCTION constructor

  SUBROUTINE initialize_logging(logfile)
    CHARACTER(len=*), INTENT(in)  :: logfile
      logging = constructor(logfile)
  END SUBROUTINE initialize_logging

  SUBROUTINE logger_message(this, message)
    CLASS(logger), INTENT(in)    :: this
    CHARACTER(len=*), INTENT(in) :: message
      WRITE(this%fileunit,*) trim(message)
  END SUBROUTINE logger_message

  SUBROUTINE logger_info(this, info)
    CLASS(logger), INTENT(in)    :: this
    CHARACTER(len=*), INTENT(in) :: info
      CALL this%message(info)
  END SUBROUTINE logger_info
  
  SUBROUTINE logger_warning(this, message, file, line)
    CLASS(logger), INTENT(in)    :: this
    CHARACTER(len=*), INTENT(in) :: message
    CHARACTER(len=*), INTENT(in), OPTIONAL :: file
    INTEGER,          INTENT(in), OPTIONAL :: line
    CHARACTER(len=filename_max) :: warning
      WRITE(warning, *)'***WARNING: ', TRIM(message)
      CALL this%message(warning)
  END SUBROUTINE logger_warning

  SUBROUTINE logger_error(this, error_message, file, line, rc)
    CLASS(logger), INTENT(in)    :: this
    CHARACTER(len=*), INTENT(in) :: error_message
    CHARACTER(len=*), INTENT(in), OPTIONAL :: file
    INTEGER,          INTENT(in), OPTIONAL :: line
    INTEGER,          INTENT(in), OPTIONAL :: rc    

    CHARACTER(len=filename_max) :: pfile
    CHARACTER(len=filename_max) :: error_info
    INTEGER :: pline
    INTEGER :: prc
    
    IF (PRESENT(file)) THEN
      pfile = file
    ELSE
      pfile=""
    ENDIF

    IF (PRESENT(line)) THEN
      pline = line
    ENDIF

    IF (PRESENT(rc)) THEN
      prc = rc
    ELSE
      prc = 1
    ENDIF
    
    WRITE(error_info,*) TRIM(pfile),' at line ', pline,' : ', TRIM(error_message)

    CALL this%message('*********************************************')
    CALL this%message('')
    CALL this%message('Abort generation of external parameters:')
    CALL this%message(error_info)
    CALL this%message('ABORT')
    CALL this%message('')
    CALL this%message('*********************************************')
    CALL exit(prc)
    
  END SUBROUTINE logger_error

  FUNCTION current_time()
    CHARACTER(len=19) :: current_time
    INTEGER :: time_vals(8)
    CHARACTER(len=48), PARAMETER :: time_format = '(i4,"-",i2.2,"-",i2.2,1x,i2.2,":",i2.2,":",i2.2)'
      CALL date_and_time(values=time_vals)
      WRITE(current_time,time_format) time_vals(1), time_vals(2), time_vals(3), time_vals(5), time_vals(6), time_vals(7)
  END FUNCTION current_time

  !> Function to get free FORTRAN unit number
  INTEGER FUNCTION free_unit_number()

    integer :: unit_no
    logical :: is_open

    free_unit_number = 1 ! start with unit 1
    is_open=.true.
    un_search: DO unit_no=1,999 
      INQUIRE (UNIT=unit_no, OPENED=is_open)
      IF (.NOT. is_open ) then
        free_unit_number = unit_no
        exit un_search
      END IF
    END DO un_search

    !jj_tmp: figure out a better way with mo_logging and mo_utitlities_extpar
!    IF (is_open) CALL abort_extpar('No free FORTRAN unit!')

  END FUNCTION free_unit_number

END MODULE mo_logging
