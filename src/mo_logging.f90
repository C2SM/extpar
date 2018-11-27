MODULE mo_logging

  USE iso_fortran_env, ONLY: error_unit, output_unit

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: notset   =  0
  INTEGER, PARAMETER, PUBLIC :: debug    = 10
  INTEGER, PARAMETER, PUBLIC :: info     = 20
  INTEGER, PARAMETER, PUBLIC :: warning  = 30
  INTEGER, PARAMETER, PUBLIC :: error    = 40
  INTEGER, PARAMETER, PUBLIC :: critical = 50

  INTEGER, PARAMETER :: default_logfile_level = warning
  INTEGER, PARAMETER :: default_stderr_level = error
  INTEGER, PARAMETER :: default_stdout_level = info

  INTEGER, PARAMETER :: closed = -1
  
  CHARACTER(len=26), PARAMETER :: logging_format = "('[',a,']','[',a,']',1x,a)"
  
  TYPE, PUBLIC :: logger
    PRIVATE
    CHARACTER(len=:), ALLOCATABLE :: logfile    
    INTEGER                       :: stdout = output_unit
    INTEGER                       :: stderr = error_unit
    INTEGER                       :: stderr_level  = notset
    INTEGER                       :: stdout_level  = notset
    INTEGER                       :: logfile_level = notset
    INTEGER                       :: fileunit = closed
  CONTAINS
    PROCEDURE :: debug    => logger_debug
    PROCEDURE :: info     => logger_info
    PROCEDURE :: warning  => logger_warning
    PROCEDURE :: error    => logger_error
    PROCEDURE :: critical => logger_critical
    PROCEDURE :: message  => logger_message
    PROCEDURE :: is_open  => logger_is_open  
  END TYPE logger

  INTERFACE logger
    MODULE PROCEDURE :: constructor
  END INTERFACE logger

  TYPE(logger), PUBLIC :: logging

  PUBLIC :: initialize_logging

  CHARACTER(len=132) :: message_text = ""

  PUBLIC :: message_text
  
CONTAINS

  FUNCTION constructor(logfile, stderr_level, stdout_level, logfile_level) RESULT(this)
    TYPE(logger) :: this
    CHARACTER(len=*), INTENT(in)  :: logfile
    INTEGER, INTENT(in), OPTIONAL :: stderr_level
    INTEGER, INTENT(in), OPTIONAL :: stdout_level
    INTEGER, INTENT(in), OPTIONAL :: logfile_level
    INTEGER :: flag
    this%logfile = logfile
    OPEN(newunit=this%fileunit,file=this%logfile,action='write',asynchronous='yes',iostat=flag,status='replace')
    IF (flag /= 0) ERROR STOP 'Error opening log file.'
    IF (PRESENT(stderr_level)) THEN
      this%stderr_level = stderr_level
    ELSE
      this%stderr_level = default_stderr_level
    END IF
    IF (PRESENT(stdout_level)) THEN
      this%stdout_level = stdout_level
    ELSE
      this%stdout_level = default_stdout_level
    END IF
    IF (PRESENT(logfile_level)) THEN
      this%logfile_level = logfile_level
    ELSE
      this%logfile_level = default_logfile_level
    END IF
  END FUNCTION constructor

  SUBROUTINE initialize_logging(logfile, stderr_level, stdout_level, logfile_level)
    CHARACTER(len=*), INTENT(in)  :: logfile
    INTEGER, INTENT(in), OPTIONAL :: stderr_level
    INTEGER, INTENT(in), OPTIONAL :: stdout_level
    INTEGER, INTENT(in), OPTIONAL :: logfile_level
    logging = logger(logfile, stderr_level, stdout_level, logfile_level)
  END SUBROUTINE initialize_logging
  
  FUNCTION current_time()
    CHARACTER(len=19) :: current_time
    INTEGER :: time_vals(8)
    CHARACTER(len=48), PARAMETER :: time_format = '(i4,"-",i2.2,"-",i2.2,1x,i2.2,":",i2.2,":",i2.2)'
    CALL date_and_time(values=time_vals)
    WRITE(current_time,time_format) time_vals(1), time_vals(2), time_vals(3), time_vals(5), time_vals(6), time_vals(7)
  END FUNCTION current_time

  SUBROUTINE logger_message(this, source, level, message)
    CLASS(logger), INTENT(in)    :: this
    CHARACTER(len=*), INTENT(in) :: source
    INTEGER, INTENT(in)          :: level
    CHARACTER(len=*), INTENT(in) :: message
    IF (level >= this%stderr_level) THEN
      WRITE(this%stderr,logging_format) current_time(), trim(source), trim(message)
    ELSE IF (level >= this%stdout_level) THEN
      WRITE(this%stdout,logging_format) current_time(), trim(source), trim(message)
    END IF
    IF (level >= this%logfile_level) THEN
      WRITE(this%fileunit,logging_format) current_time(), trim(source), trim(message)
    END IF
  END SUBROUTINE logger_message
  
  SUBROUTINE logger_debug(this, message, file, line)
    CLASS(logger), INTENT(in)    :: this
    CHARACTER(len=*), INTENT(in) :: message
    CHARACTER(len=*), INTENT(in) :: file
    INTEGER, INTENT(in) :: line
    CHARACTER(len=1024) :: source
    WRITE(source,'(a,":",i0)') TRIM(file), line
    CALL this%message(source, debug, message)
  END SUBROUTINE logger_debug

  SUBROUTINE logger_info(this, message, file, line)
    CLASS(logger), INTENT(in)    :: this
    CHARACTER(len=*), INTENT(in) :: message
    CHARACTER(len=*), INTENT(in) :: file
    INTEGER, INTENT(in) :: line
    CHARACTER(len=1024) :: source
    WRITE(source,'(a,":",i0)') TRIM(file), line
    CALL this%message(source, info, message)
  END SUBROUTINE logger_info
  
  SUBROUTINE logger_warning(this, message, file, line)
    CLASS(logger), INTENT(in)    :: this
    CHARACTER(len=*), INTENT(in) :: message
    CHARACTER(len=*), INTENT(in) :: file
    INTEGER, INTENT(in) :: line
    CHARACTER(len=1024) :: source
    WRITE(source,'(a,":",i0)') TRIM(file), line
    CALL this%message(source, warning, message)
  END SUBROUTINE logger_warning
  
  SUBROUTINE logger_error(this, message, file, line)
    CLASS(logger), INTENT(in)    :: this
    CHARACTER(len=*), INTENT(in) :: message
    CHARACTER(len=*), INTENT(in) :: file
    INTEGER, INTENT(in) :: line
    CHARACTER(len=1024) :: source
    WRITE(source,'(a,":",i0)') TRIM(file), line
    CALL this%message(source, error, message)
  END SUBROUTINE logger_error
  
  SUBROUTINE logger_critical(this, message, file, line)
    CLASS(logger), INTENT(in)    :: this
    CHARACTER(len=*), INTENT(in) :: message
    CHARACTER(len=*), INTENT(in) :: file
    INTEGER, INTENT(in) :: line
    CHARACTER(len=1024) :: source
    WRITE(source,'(a,":",i0)') TRIM(file), line
    CALL this%message(trim(source), critical, message)
  END SUBROUTINE logger_critical

  PURE FUNCTION logger_is_open(this)
    CLASS(logger), INTENT(in) :: this
    LOGICAL :: logger_is_open
    logger_is_open = (this%fileunit /= closed)
  END FUNCTION logger_is_open

END MODULE mo_logging
