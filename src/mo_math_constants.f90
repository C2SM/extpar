!+ Module determines main mathematical constants for ICON
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
!> Module determines main mathematical constants for ICON
!!  prototype. Values are taken from glibc 2.2.5: /usr/include/math.h
!!  These constants
!!  are provided to more significant digits than is necessary for a 64-bit
!!  double precision number; they may be used for other purposes where the
!!  extra precision is necessary or useful.
!! 
!! @par Revision History
!!  Developed  by Luis Kornblueh (2004)
!!  Modified to ProTeX-style by  Luca Bonaventura and Thomas Heinze (2004).
!!  Modified according to style guide by Thomas Heinze (2005-06-24):
!!   - module renamed from mo_math to mo_math_constants
!!   - eps moved from mo_physical_constants
!!   - pid180i renamed to rad2deg
!!  Including some more constants from math.h by Thomas Heinze (2005-07-18)
!!  pid5 renamed to pi_5 by Thomas Heinze (2005-07-26)
!!  Modification by Thomas Heinze (2006-02-21):
!!  - renamed m_modules to mo_modules
!!  Modification by Thomas Heinze (2006-05-18):
!!  - introduced dbl_eps
!! 
MODULE mo_math_constants
USE mo_kind,      ONLY:  wp

IMPLICIT NONE

! !VERSION CONTROL:
CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id: mo_math_constants.f90,v 1.11 2013-04-16 11:09:20 for0adm Exp $'

PUBLIC

! !DEFINED PARAMETERS:

! 
! ! Mathematical constants (in brackets original C names)
!
! ! euler     (mo_E       )  -- e
! ! log2e     (mo_LOG2E   )  -- log2(e)
! ! log10e    (mo_LOG10E  )  -- log10(e)
! ! ln2       (mo_LN2     )  -- ln(2)
! ! ln10      (mo_LN10    )  -- ln(10)
! ! pi        (mo_PI      )  -- pi
! ! pi_2      (mo_PI_2    )  -- pi/2
! ! pi_4      (mo_PI_4    )  -- pi/4
! ! rpi       (mo_1_PI    )  -- 1/pi
! ! rpi_2     (mo_2_PI    )  -- 2/pi
! ! rsqrtpi_2 (mo_2_SQRTPI)  -- 2/(sqrt(pi))
! ! sqrt2     (mo_SQRT2   )  -- sqrt(2)
! ! sqrt1_2   (mo_SQRT1_2 )  -- 1/sqrt(2)
! !
! 
!

REAL (wp), PARAMETER ::  euler     = 2.71828182845904523536028747135266250_wp
REAL (wp), PARAMETER ::  log2e     = 1.44269504088896340735992468100189214_wp
REAL (wp), PARAMETER ::  log10e    = 0.434294481903251827651128918916605082_wp
REAL (wp), PARAMETER ::  ln2       = 0.693147180559945309417232121458176568_wp
REAL (wp), PARAMETER ::  ln10      = 2.30258509299404568401799145468436421_wp
REAL (wp), PARAMETER ::  pi        = 3.14159265358979323846264338327950288_wp
REAL (wp), PARAMETER ::  pi_2      = 1.57079632679489661923132169163975144_wp
REAL (wp), PARAMETER ::  pi_4      = 0.785398163397448309615660845819875721_wp
REAL (wp), PARAMETER ::  rpi       = 0.318309886183790671537767526745028724_wp
REAL (wp), PARAMETER ::  rpi_2     = 0.636619772367581343075535053490057448_wp
REAL (wp), PARAMETER ::  rsqrtpi_2 = 1.12837916709551257389615890312154517_wp
REAL (wp), PARAMETER ::  sqrt2     = 1.41421356237309504880168872420969808_wp
REAL (wp), PARAMETER ::  sqrt1_2   = 0.707106781186547524400844362104849039_wp

!
! ! some more useful constants
! ! pi_5    -- half angle of pentagon
! ! rad2deg -- conversion factor from radians to degree
! ! deg2rad -- conversion factor from degree to radians
! ! eps     -- residual bound for solvers
!

REAL (wp), PARAMETER ::  pi_5      = pi*0.2_wp     
REAL (wp), PARAMETER ::  rad2deg   = 180.0_wp/pi 
REAL (wp), PARAMETER ::  deg2rad   = pi/180.0_wp
REAL (wp), PARAMETER ::  eps       = 1.e-8_wp      
REAL (wp), PARAMETER ::  dbl_eps   = EPSILON(1._wp)      

!
! ! phi0 is  the latitude of the lowest major triangle corner
! ! and the latitude of the major hexagonal faces centers
! ! phi0 = 0.5_wp*pi - 2._wp*acos(1.0_wp/(2._wp*sin(pi/5._wp)))  
!

REAL (wp), PARAMETER :: phi0  = 0.46364760900080614903_wp 


!EOP  
!--------------------------------------------------------------------  
!BOC  

END MODULE mo_math_constants

!EOC 
