!+ Sets the standard I/O units
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_14        2014-07-18 Juergen Helmert
!  Combined COSMO Release
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!>
!! Sets the standard I/O units.
!! 
!! 
!! @par Revision History
!!  Initial version  by Luis Kornblueh (2004).
!!  Modification by Thomas Heinze (2005-08-29):
!!   - no more comments in ifdef declarations
!!   - included VERSION CONTROL
!!  Modification by Thomas Heinze (2006-02-21):
!!  - renamed m_modules to mo_modules
!!  Modification by Hui Wan (2007-02):
!!  - parameters <i>input1</i> and <i>input3</i> deleted.
!!    (They specified the input namelist files of the grid generator,
!!    but were not used at all in the shallow water model.)
!!  - 'ini_file' renamed as 'namelist_file' and set to 'NAMELIST_ICON'
!!  - 'root_dir' defined, which is used for setting the default path to the
!!     files providing grid/patch/topography data.
!!  Modification by Thomas Heinze, DWD, (2007-03-05):
!!  - reintroduced nlog
!!  Modification by Thomas Heinze, DWD, (2007-05-10):
!!  - moved unit1, unit2 and output1 from mo_io_graph here to avoid
!!    cyclic dependencies in Makefile
!!  Modifications by Luis Kornblueh, MPI-M, 2008-04-21:
!!  - going back to original design
!!  - added some predefined units
!!  - added function to find next free unit
!! 
MODULE mo_io_units
  IMPLICIT NONE

  
  PUBLIC

! This paramter is taken from /usr/include/stdio.h (ANSI C standard). If problems
! with filename length appear, check the before mentioned file.

  INTEGER, PARAMETER :: filename_max = 1024
  
! Standard I/O-units

  INTEGER, PARAMETER :: nerr  = 0     ! error output
  INTEGER, PARAMETER :: nlog  = 1     ! standard log file unit
  INTEGER, PARAMETER :: nnml  = 2     ! standard namelist file unit
  INTEGER, PARAMETER :: nstat = 3     ! standard statistics file unit
  INTEGER, PARAMETER :: ngmt  = 4     ! standard GMT output file unit
  INTEGER, PARAMETER :: nin   = 5     ! standard input
  INTEGER, PARAMETER :: nout  = 6     ! standard output  

  INTEGER, PARAMETER :: nan   = -1    ! unit given back, when nothing 
                                      ! in the allowed range is available

  INTEGER :: nnml_output ! unit of the ASCII output that contains the 
                         ! namelist variables and their actual values.
CONTAINS

  FUNCTION find_next_free_unit(istart,istop) RESULT(iunit)
    INTEGER :: iunit
    INTEGER, INTENT(in) :: istart, istop
    INTEGER :: kstart, kstop
    LOGICAL :: found, opened
    INTEGER :: i
    CHARACTER(len=32) :: info
    
    found = .FALSE.
    kstart = istart
    kstop  = istop
    IF (kstart < 10) kstart = 10
    IF (kstop <= kstart) kstop = kstart+10
    DO i = kstart, kstop
      INQUIRE(unit=i, opened=opened)
      IF (.NOT. opened) THEN
        iunit = i
        found = .TRUE.
        EXIT
      END IF
    END DO
    
    IF (.NOT. found) THEN
      WRITE(info,'(a,i0,a,i0,a)') &
           'No unit in range <', kstart, ':', kstop, '> free.'
      WRITE(nerr,'(a,a,a)') 'find_next_free_unit', ': ',TRIM(info)
      iunit = nan
    END IF
    
  END FUNCTION find_next_free_unit

END MODULE mo_io_units

!EOC













