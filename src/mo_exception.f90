!+ Fortran module with routines for exception handling
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module with routines for exception handling
MODULE mo_exception

  USE mo_io_units, ONLY: nerr, nlog 

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: message_text
  PUBLIC :: message, finish
  PUBLIC :: em_none, em_info, em_warn

  PUBLIC :: open_log, close_log

  INTEGER, PARAMETER :: em_none = 0 
  INTEGER, PARAMETER :: em_info = 1
  INTEGER, PARAMETER :: em_warn = 2

  CHARACTER(1024) :: message_text = ''

  LOGICAL :: l_log = .FALSE.

CONTAINS

  SUBROUTINE finish (name, text, exit_no)


    CHARACTER(len=*) :: name
    CHARACTER(len=*), OPTIONAL :: text
    INTEGER, OPTIONAL :: exit_no
    INTEGER           :: iexit


    WRITE (nerr,'(/,80("*"),/)')
    IF (l_log) WRITE (nlog,'(/,80("*"),/)')

    IF (PRESENT(exit_no)) THEN
       iexit = exit_no
    ELSE
       iexit = 1
    END IF

    IF (PRESENT(text)) THEN
      WRITE (nerr,'(1x,a,a,a)') TRIM(name), ': ', TRIM(text)
      IF (l_log) WRITE (nlog,'(1x,a,a,a)') TRIM(name), ': ', TRIM(text)
    ELSE
      WRITE (nerr,'(1x,a,a)') TRIM(name), ': '
      IF (l_log) WRITE (nlog,'(1x,a,a)') TRIM(name), ': '
    ENDIF


    WRITE (nerr,'(/,80("-"),/,/)')
    IF (l_log) WRITE (nlog,'(/,80("-"),/,/)')


    WRITE (nerr,'(/,80("*"),/)')

       STOP 'mo_exception: finish ..'

  END SUBROUTINE finish

  SUBROUTINE message (name, text, out, level, all_print)

    CHARACTER (*) :: name, text
    INTEGER, INTENT(in), OPTIONAL :: out
    INTEGER, INTENT(in), OPTIONAL :: level
    LOGICAL, INTENT(in), OPTIONAL :: all_print

    INTEGER :: iout
    INTEGER :: ilevel
    LOGICAL :: lprint

    name=name
    text=text

    IF (PRESENT(all_print)) THEN
      lprint = all_print
    ELSE
      lprint = .FALSE.
    ENDIF

    IF (PRESENT(out)) THEN
      iout = out
    ELSE
      iout = nerr
    END IF

    IF (PRESENT(level)) THEN
      ilevel = level
    ELSE
      ilevel = em_none
    END IF

    SELECT CASE (ilevel)
    CASE (em_none)
    CASE (em_info)
    CASE (em_warn)
    END SELECT


  END SUBROUTINE message

  SUBROUTINE open_log (logfile_name)

    CHARACTER(len=*), INTENT(in) :: logfile_name
    LOGICAL                      :: l_opened

    INQUIRE (UNIT=nlog,OPENED=l_opened)

    IF (l_opened) THEN
      WRITE (message_text,'(a)') 'log file unit has been used already.'
      CALL message ('open_log', message_text)
      WRITE (message_text,'(a)') 'Close unit and reopen for log file.'
      CALL message ('open_log', message_text, level=em_warn)
      CLOSE (nlog)
    ENDIF

    OPEN (nlog,file=TRIM(logfile_name))
  
    l_log = .TRUE.

  END SUBROUTINE open_log

  SUBROUTINE close_log
    LOGICAL :: l_opened
   
    INQUIRE (UNIT=nlog,OPENED=l_opened)
    IF (l_opened) THEN
      CLOSE (nlog)
    ENDIF

    l_log = .FALSE.

  END SUBROUTINE close_log

END MODULE mo_exception
