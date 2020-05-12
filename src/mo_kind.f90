!+ Module determines kinds for different precisions.
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================!>
!!  Module determines kinds for different precisions.
!! 
!!  Module determines kinds for different precisions
!!  Number model from which the SELECTED_*\\_KIND are requested: <br>
!!
!! @f{tabular}{{r@{\hspace*{3em}}c@{\hspace*{3em}}c}
!!                     &4 byte REAL     &8 byte REAL        \\\
!!        CRAY:        &-               &precision =   13   \\\
!!                     &                &exponent  = 2465   \\\
!!        IEEE:        &precision = 6   &precision =   15   \\\
!!                     &exponent  = 37  &exponent  =  307
!! @f}
!! \\medskip
!! 
!!  Most likely this are the only possible models.
!! 
!! @par Revision History
!!  Working precision and comments are added by Luis Kornblueh (2001)
!!  Modification by Luis Kornblueh (2010-02-16):
!!  Working precision selection is parameterized by values, which can be as 
!!  well used for selection of MPI data types associated with the respective 
!!  Fortran datatypes.
!! 
MODULE mo_kind

  IMPLICIT NONE

  ! floating point section 

  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6,37)  
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)
  INTEGER, PARAMETER :: qp = SELECTED_REAL_KIND(33,4931)

  INTEGER, PARAMETER :: wp = dp

  ! integer section

  INTEGER, PARAMETER :: i1 = SELECTED_INT_KIND(2)
  INTEGER, PARAMETER :: i2 = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)

END MODULE mo_kind
