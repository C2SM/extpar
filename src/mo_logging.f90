MODULE mo_logging

  IMPLICIT NONE

  !Integer for debugging levels
  INTEGER, PARAMETER, PUBLIC :: verbose = 0 ! verbosity of extpar, add to namelist
  INTEGER, PARAMETER, PUBLIC :: idbg_low  = 1 ! low debug output
  INTEGER, PARAMETER, PUBLIC :: idbg_high = 2 ! high debug output
  INTEGER, PARAMETER :: closed = -1
  
  TYPE, PUBLIC :: logger
    CHARACTER(len=:), ALLOCATABLE :: logfile    
    INTEGER                       :: fileunit
  END TYPE logger

  TYPE(logger), PUBLIC :: logging

  PUBLIC :: initialize_logging

  CHARACTER(len=132) :: message_text = ""

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
