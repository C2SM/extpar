!+ Fortran module for utility routines
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_4         2011/04/21 Anne Roches
! two new SR extend_field2D, horizontal_filtering
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module for utility routines
!>
!> This module is based on the lm_f90 module utilities.f90 by Ulrich Schaettler
!> with some additional routines by Hermann Asensio. The module provides service utilities.
!!
!! Description:
!!  This module provides service utilities for the model. All routines are 
!!   written in a manner that also other models can use it. That means:
!!     - no routine uses other modules, except the declarations for the 
!!       KIND-type parameter; the data access is by parameter list only
!!     - no routine allocates dynamic memory; work space needed is
!!       provided via the parameter list
!!     - no derived data types are used
!!
!!     - free_un:
!!       Function to get free FORTRAN unit number. 
!!     
!!     - abort_extpar:
!!       Subroutine to abort generation of external parameter program
!!
!!     - phirot2phi:
!!       Converts phi from the rotated system to phi in the real
!!       geographical system.
!!
!!     - phi2phirot:
!!       Converts phi from the real geographical system to phi
!!       in the rotated system.
!!
!!     - rlarot2rla:
!!       Converts lambda from the rotated system to lambda in the real
!!       geographical system.
!!
!!     - rla2rlarot:
!!       Converts lambda from the real geographical system to lambda 
!!       in the rotated system.
!!
!!     - uvrot2uv:
!!       Converts the wind components u and v from the rotated system
!!       to the real geographical system.
!!
!!     - uvrot2uv_vec:
!!       the same as above, but for a whole 2D field (in vectorized form).
!!
!!     - uv2uvrot:
!!       Converts the wind components u and v from the real geographical
!!       system to the rotated system.
!!
!!     - uv2uvrot_vec:
!!       the same as above, but for a whole 2D field (in vectorized form).
!!
!!     - uv2df:
!!       Converts the wind components u and v to wind direction and speed.
!!
!!     - uv2df_vec:
!!       the same as above, but for a whole 2D field (in vectorized form).
!!
!roa>
!!     - extend_field2D:
!!       extend a 2D field with additional rows and columns.
!!
!!     - horizontal_filtering:
!!       filter horizontally a field using a local filter
!!       (weights determined by J. Foerstner, DWD).
!!
!!
!roa<

!!
MODULE  mo_utilities_extpar

  USE mo_kind,     ONLY: wp, i4, i4
  USE mo_logging
  USE mo_io_units, ONLY: filename_max
  
#ifdef NAGFOR
  USE f90_unix, ONLY: exit
#endif
    
  IMPLICIT NONE

  PUBLIC

CONTAINS

  SUBROUTINE check_input_file(filename, file, line)
    CHARACTER(len=*), INTENT(in) :: filename
    CHARACTER(len=*), INTENT(in) :: file
    INTEGER, INTENT(in) :: line
    LOGICAL :: exists = .FALSE.
    INQUIRE(file=TRIM(filename), exist=exists)
    IF (exists) THEN
      WRITE(logging%fileunit,*)TRIM(filename)//' ... exists'
    ELSE
      WRITE(logging%fileunit,*)TRIM(filename)//' ... no such file'
      CALL abort_extpar('Missing input file ...')
    ENDIF
    
  END SUBROUTINE check_input_file
    
  !> Function to get free FORTRAN unit number
  !!
  !! Method: use INQUIRE for unit number
  INTEGER FUNCTION free_un()

    integer :: unit_no
    logical :: is_open

    free_un = 1 ! start with unit 1
    is_open=.true.
    un_search: DO unit_no=1,999 
      INQUIRE (UNIT=unit_no, OPENED=is_open)
      IF (.NOT. is_open ) then
        free_un = unit_no
        exit un_search
      END IF
    END DO un_search

    IF (is_open) CALL abort_extpar('No free FORTRAN unit!')

  END FUNCTION free_un


  !> Subroutine to abort generation of external parameter program
  SUBROUTINE abort_extpar(error_message, file, line, rc)
    CHARACTER(len=*), INTENT(in) :: error_message
    CHARACTER(len=*), INTENT(in), OPTIONAL :: file
    INTEGER,          INTENT(in), OPTIONAL :: line
    INTEGER,          INTENT(in), OPTIONAL :: rc    

    CHARACTER(len=filename_max) :: pfile
    INTEGER :: pline
    INTEGER :: prc
    
    IF (PRESENT(file)) THEN
      pfile = file
    ELSE
      pfile = __FILE__
    ENDIF

    IF (PRESENT(line)) THEN
      pline = line
    ELSE
      pline = __LINE__
    ENDIF

    IF (PRESENT(rc)) THEN
      prc = rc
    ELSE
      prc = 1
    ENDIF
    
    WRITE(logging%fileunit,*)'*********************************************'
    WRITE(logging%fileunit,*)''
    WRITE(logging%fileunit,*)'Abort generation of external parameters:'
    WRITE(logging%fileunit,*)TRIM(error_message)
    WRITE(logging%fileunit,*)'ABORT'
    WRITE(logging%fileunit,*)''
    WRITE(logging%fileunit,*)'*********************************************'
    CALL exit(prc)
    
  END SUBROUTINE abort_extpar

  !> Subroutine to get coordinates of rotated soutpole from coordinates of rotated northpol
  SUBROUTINE get_rot_spol_coor(npollon,npollat,spollon,spollat)

    REAL (KIND=WP), INTENT(IN)  :: npollon !< longitude of the rotated north pole
    REAL (KIND=WP), INTENT(IN)  :: npollat !< latitude of the rotated north pole
    REAL (KIND=WP), INTENT(OUT) :: spollon !< longitude of the rotated south pole
    REAL (KIND=WP), INTENT(OUT) :: spollat !< latitude of the rotated south pole

    spollat = -1. * npollat ! latitude coordinate of rotated southpole

    spollon = npollon + 180. ! longitude coordinate of rotated soutpole
    IF (spollon > 180.00) THEN ! shift the longitude range
      spollon = spollon - 360.
    ENDIF

  END SUBROUTINE get_rot_spol_coor

  !> Function for rotation of geographical coordinates
  !>
  !! Description:
  !!   This function converts phi from one rotated system to phi in another
  !!   system. The other system
  !!   can also be a rotated one, where polgam is the angle between the two
  !!   north poles.
  !!   If polgam =0 , the other system is the real geographical
  !!   system.
  !!
  !! Method:
  !!   Transformation formulas for converting between these two systems.
  !!
  FUNCTION  phirot2phi ( phirot, rlarot, polphi, polgam )

    REAL (KIND=wp), INTENT (IN)      ::    polphi !< latitude of the rotated north pole
    REAL (KIND=wp), INTENT (IN)      ::    phirot !< latitude in the rotated system
    REAL (KIND=wp), INTENT (IN)      ::    rlarot !< longitude in the rotated system

    REAL (KIND=wp), INTENT (IN)      ::        &
         polgam       !< angle between the north poles of the systems

    REAL (KIND=wp)                   ::        &
         phirot2phi  !< latitude in the geographical system

    ! Local variables
    REAL (KIND=wp)                   ::        &
         zsinpol, zcospol, zphis, zrlas, zarg, zgam

    REAL (KIND=wp), PARAMETER        ::        &
         zrpi18 = 57.29577951308232_wp,                  & !
         zpir18 = 0.017453292519943295_wp

    !------------------------------------------------------------------------------

    ! Begin function phirot2phi

    zsinpol     = SIN (zpir18 * polphi)
    zcospol     = COS (zpir18 * polphi)

    zphis       = zpir18 * phirot
    IF (rlarot > 180.0_wp) THEN
      zrlas = rlarot - 360.0_wp
    ELSE
      zrlas = rlarot
    ENDIF
    zrlas       = zpir18 * zrlas

    IF (polgam /= 0.0_wp) THEN
      zgam  = zpir18 * polgam
      zarg  = zsinpol*SIN (zphis) +                                           &
           zcospol*COS(zphis) * ( COS(zrlas)*COS(zgam) - SIN(zgam)*SIN(zrlas) )
    ELSE
      zarg  = zcospol * COS (zphis) * COS (zrlas) + zsinpol * SIN (zphis)
    ENDIF

    phirot2phi  = zrpi18 * ASIN (zarg)

  END FUNCTION phirot2phi

  !==============================================================================
  !==============================================================================

  !------------------------------------------------------------------------------

  !> This routine converts phi from the real geographical system to phi in the rotated system
  FUNCTION  phi2phirot ( phi, rla, polphi, pollam )

    !------------------------------------------------------------------------------
    ! Description:
    !   This routine converts phi from the real geographical system to phi
    !   in the rotated system.
    !
    ! Method:
    !   Transformation formulas for converting between these two systems.
    !
    !------------------------------------------------------------------------------
    ! Parameter list:

    REAL (KIND=wp), INTENT (IN)      ::    polphi !< latitude of the rotated north pole
    REAL (KIND=wp), INTENT (IN)      ::    pollam !< longitude of the rotated north pole
    REAL (KIND=wp), INTENT (IN)      ::    phi !< latitude in the geographical system
    REAL (KIND=wp), INTENT (IN)      ::    rla !< longitude in the geographical system

    REAL (KIND=wp)                   ::        &
         phi2phirot !< longitude in the rotated system

    ! Local variables
    REAL (KIND=wp)                       ::    &
         zsinpol, zcospol, zlampol, zphi, zrla, zarg1, zarg2, zrla1

    REAL (KIND=wp), PARAMETER            ::    &
         zrpi18 = 57.29577951308232_wp,                  & !
         zpir18 = 0.017453292519943295_wp

    !------------------------------------------------------------------------------

    ! Begin function phi2phirot

    zsinpol  = SIN (zpir18 * polphi)
    zcospol  = COS (zpir18 * polphi)
    zlampol  =      zpir18 * pollam
    zphi     =      zpir18 * phi
    IF (rla > 180.0_wp) THEN
      zrla1  = rla - 360.0_wp
    ELSE
      zrla1  = rla
    ENDIF
    zrla     = zpir18 * zrla1

    zarg1    = SIN (zphi) * zsinpol
    zarg2    = COS (zphi) * zcospol * COS (zrla - zlampol)

    phi2phirot = zrpi18 * ASIN (zarg1 + zarg2)

  END FUNCTION phi2phirot

  !==============================================================================
  !==============================================================================

  !------------------------------------------------------------------------------
  !> This function converts lambda from one rotated system to lambda in another.
  !!
  !! Description:
  !!   This function converts lambda from one rotated system to lambda in another
  !!   system. The other system
  !!   can also be a rotated one, where polgam is the angle between the two
  !!   north poles.
  !!   If polgam =0 , the other system is the real geographical
  !!   system.
  FUNCTION  rlarot2rla (phirot, rlarot, polphi, pollam, polgam)

    !
    ! Method:
    !   Transformation formulas for converting between these two systems.
    !
    ! Modules used:    NONE
    !
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    !
    ! Declarations:
    !
    !------------------------------------------------------------------------------

    ! Parameter list:
    REAL (KIND=wp), INTENT (IN)      ::   polphi !< latitude of the rotated north pole
    REAL (KIND=wp), INTENT (IN)      ::   pollam !< longitude of the rotated north pole
    REAL (KIND=wp), INTENT (IN)      ::   phirot !< latitude in the rotated system
    REAL (KIND=wp), INTENT (IN)      ::   rlarot !< longitude in the rotated system

    REAL (KIND=wp), INTENT (IN)      ::        &
         polgam      !< angle between the north poles of the systems

    REAL (KIND=wp)                   ::        &
         rlarot2rla  !< latitude in the geographical system

    ! Local variables
    REAL (KIND=wp)                   ::        &
         zsinpol, zcospol, zlampol, zphis, zrlas, zarg1, zarg2, zgam

    REAL (KIND=wp), PARAMETER        ::        &
         zrpi18 = 57.29577951308232_wp,                  & !
         zpir18 = 0.017453292519943295_wp

    !------------------------------------------------------------------------------

    ! Begin function rlarot2rla

    zsinpol = SIN (zpir18 * polphi)
    zcospol = COS (zpir18 * polphi)

    zlampol = zpir18 * pollam
    zphis   = zpir18 * phirot
    IF (rlarot > 180.0_wp) THEN
      zrlas = rlarot - 360.0_wp
    ELSE
      zrlas = rlarot
    ENDIF
    zrlas   = zpir18 * zrlas

    IF (polgam /= 0.0_wp) THEN
      zgam    = zpir18 * polgam
      zarg1   = SIN (zlampol) *                                                &
           (- zsinpol*COS(zphis) * (COS(zrlas)*COS(zgam) - SIN(zrlas)*SIN(zgam))  &
           + zcospol * SIN(zphis))                                               &
           - COS (zlampol)*COS(zphis) * (SIN(zrlas)*COS(zgam) + COS(zrlas)*SIN(zgam))

      zarg2   = COS (zlampol) *                                                &
           (- zsinpol*COS(zphis) * (COS(zrlas)*COS(zgam) - SIN(zrlas)*SIN(zgam))  &
           + zcospol * SIN(zphis))                                               &
           + SIN (zlampol)*COS(zphis) * (SIN(zrlas)*COS(zgam) + COS(zrlas)*SIN(zgam))
    ELSE
      zarg1   = SIN (zlampol) * (-zsinpol * COS(zrlas) * COS(zphis)  +    &
           zcospol *              SIN(zphis)) -    &
           COS (zlampol) *             SIN(zrlas) * COS(zphis)
      zarg2   = COS (zlampol) * (-zsinpol * COS(zrlas) * COS(zphis)  +    &
           zcospol *              SIN(zphis)) +   &
           SIN (zlampol) *             SIN(zrlas) * COS(zphis)
    ENDIF

    IF (zarg2 == 0.0) zarg2 = 1.0E-20_wp

    rlarot2rla = zrpi18 * ATAN2(zarg1,zarg2)

  END FUNCTION rlarot2rla

  !==============================================================================
  !==============================================================================

  !------------------------------------------------------------------------------
  !>  This routine converts lambda from the real geographical system to lambda in the rotated system.
  FUNCTION  rla2rlarot ( phi, rla, polphi, pollam, polgam )

    !------------------------------------------------------------------------------
    !
    ! Description:
    !   This routine converts lambda from the real geographical system to lambda 
    !   in the rotated system.
    !
    ! Method:
    !   Transformation formulas for converting between these two systems.
    !
    !------------------------------------------------------------------------------
    !
    ! Parameter list:
    REAL (KIND=wp), INTENT (IN)      ::   polphi !< latitude of the rotated north pole
    REAL (KIND=wp), INTENT (IN)      ::   pollam !< longitude of the rotated north pole
    REAL (KIND=wp), INTENT (IN)      ::   phi    !< latitude in geographical system
    REAL (KIND=wp), INTENT (IN)      ::   rla        !< longitude in geographical system

    REAL (KIND=wp), INTENT (IN)      ::        &
         polgam      !< angle between the north poles of the systems

    REAL (KIND=wp)                   ::        &
         rla2rlarot !< latitude in the the rotated system

    ! Local variables
    REAL (KIND=wp)                       ::    &
         zsinpol, zcospol, zlampol, zphi, zrla, zarg1, zarg2, zrla1

    REAL (KIND=wp), PARAMETER            ::    &
         zrpi18 = 57.29577951308232_wp,                  & !
         zpir18 = 0.017453292519943295_wp

    !------------------------------------------------------------------------------

    ! Begin function rla2rlarot

    zsinpol  = SIN (zpir18 * polphi)
    zcospol  = COS (zpir18 * polphi)
    zlampol  =      zpir18 * pollam
    zphi     =      zpir18 * phi
    IF (rla > 180.0_wp) THEN
      zrla1  = rla - 360.0_wp
    ELSE
      zrla1  = rla
    ENDIF
    zrla     = zpir18 * zrla1

    zarg1    = - SIN (zrla-zlampol) * COS(zphi)
    zarg2    = - zsinpol * COS(zphi) * COS(zrla-zlampol) + zcospol * SIN(zphi)

    IF (zarg2 == 0.0) zarg2 = 1.0E-20_wp

    rla2rlarot = zrpi18 * ATAN2 (zarg1,zarg2)

    IF (polgam /= 0.0_wp ) THEN
      rla2rlarot = polgam + rla2rlarot
      IF (rla2rlarot > 180._wp) rla2rlarot = rla2rlarot -360._wp
    ENDIF

  END FUNCTION rla2rlarot

  !==============================================================================
  !==============================================================================

  !------------------------------------------------------------------------------
  !>  This routine converts the wind components u and v from the rotated system to the real geographical system.
  SUBROUTINE uvrot2uv (urot, vrot, rlat, rlon, pollat, pollon, u, v)


    !------------------------------------------------------------------------------
    !
    ! Description:
    !   This routine converts the wind components u and v from the rotated system
    !   to the real geographical system.
    !
    ! Method:
    !   Transformation formulas for converting between these two systems.
    !
    !------------------------------------------------------------------------------

    ! Parameter list:
    REAL (KIND=wp), INTENT (IN)          ::  urot      !< wind components in the rotated grid
    REAL (KIND=wp), INTENT (IN)          ::  vrot      !< wind components in the rotated grid
    REAL (KIND=wp), INTENT (IN)          ::  rlat      !< latitude in the true geographical system
    REAL (KIND=wp), INTENT (IN)          ::  rlon      !< longitude in the true geographical system
    REAL (KIND=wp), INTENT (IN)          ::  pollat    !< latitude  of the north pole of the rotated grid
    REAL (KIND=wp), INTENT (IN)          ::  pollon    !<  longitude of the north pole of the rotated grid

    REAL (KIND=wp), INTENT (OUT)         ::   u        !< wind components in the true geographical system
    REAL (KIND=wp), INTENT (OUT)         ::   v        !< wind components in the true geographical system

    ! Local variables

    REAL (KIND=wp)                       ::    &
         zsinpol, zcospol, zlonp, zlat, zarg1, zarg2, znorm

    REAL (KIND=wp)                       ::    &
         zpir18 = 0.0174532925_wp

    !------------------------------------------------------------------------------
    ! Begin subroutine uvrot2uv
    !------------------------------------------------------------------------------

    ! Converting from degree to radians
    zsinpol = SIN(pollat * zpir18)
    zcospol = COS(pollat * zpir18)
    zlonp   = (pollon-rlon) * zpir18
    zlat    =         rlat  * zpir18

    zarg1   = zcospol*SIN(zlonp)
    zarg2   = zsinpol*COS(zlat) - zcospol*SIN(zlat)*COS(zlonp)
    znorm   = 1.0/SQRT(zarg1**2 + zarg2**2)

    ! Convert the u- and v-components
    u       =   urot*zarg2*znorm + vrot*zarg1*znorm
    v       = - urot*zarg1*znorm + vrot*zarg2*znorm

  END SUBROUTINE uvrot2uv

  !==============================================================================
  !==============================================================================

  !------------------------------------------------------------------------------
  !> This routine converts the wind components u and v from the rotated system to the real geographical system.
  !> This is the vectorized form.
  SUBROUTINE uvrot2uv_vec(u, v, rlat, rlon, pollat, pollon, idim, jdim)

    !------------------------------------------------------------------------------
    !
    ! Description:
    !   This routine converts the wind components u and v from the rotated
    !   system to the real geographical system. This is the vectorized form
    !   of the routine above, i.e. the computation is for a whole 2D field.
    !
    ! Method:
    !   Transformation formulas for converting between these two systems.
    !
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Parameter list:
    INTEGER, INTENT(IN)     ::    &
         idim        !< dimensions of the field
    INTEGER, INTENT(IN)     ::    &
         jdim        !< dimensions of the field


    REAL (KIND=wp), INTENT (INOUT)       ::    &
         u  (idim,jdim) !< wind components 
    REAL (KIND=wp), INTENT (INOUT)       ::    &
         v  (idim,jdim)    !< wind components 

    REAL (KIND=wp), INTENT (IN)          ::  rlat(idim,jdim)      !< latitude in the true geographical system
    REAL (KIND=wp), INTENT (IN)          ::  rlon (idim,jdim)     !< longitude in the true geographical system
    REAL (KIND=wp), INTENT (IN)          ::  pollat    !< latitude  of the north pole of the rotated grid
    REAL (KIND=wp), INTENT (IN)          ::  pollon    !<  longitude of the north pole of the rotated grid


    ! Local variables
    REAL (KIND=wp)                       ::    &
         zsinpol, zcospol, zlonp, zlat, zarg1, zarg2, znorm, zugeo, zvgeo

    INTEGER ::    i, j
    REAL (KIND=wp)                       ::    &
         zpir18 = 0.0174532925_wp

    !------------------------------------------------------------------------------
    ! Begin Subroutine uvrot2uv_vec
    !------------------------------------------------------------------------------

    ! Converting from degree to radians
    zsinpol = SIN(pollat * zpir18)
    zcospol = COS(pollat * zpir18)

    DO j = 1, jdim
      DO i = 1, idim

        zlonp   = (pollon-rlon(i,j)) * zpir18
        zlat    =         rlat(i,j)  * zpir18

        zarg1   = zcospol*SIN(zlonp)
        zarg2   = zsinpol*COS(zlat) - zcospol*SIN(zlat)*COS(zlonp)
        znorm   = 1.0/SQRT(zarg1**2 + zarg2**2)

        ! Convert the u- and v-components
        zugeo   =  u(i,j)*zarg2*znorm + v(i,j)*zarg1*znorm
        zvgeo   = -u(i,j)*zarg1*znorm + v(i,j)*zarg2*znorm
        u(i,j) = zugeo
        v(i,j) = zvgeo

      ENDDO
    ENDDO

  END SUBROUTINE uvrot2uv_vec

  !==============================================================================
  !==============================================================================

  !------------------------------------------------------------------------------
  !> This routine converts the wind components u and v from the real geographical system to the rotated system.
  SUBROUTINE uv2uvrot(u, v, rlat, rlon, pollat, pollon, urot, vrot)

    !------------------------------------------------------------------------------
    !
    ! Description:
    !   This routine converts the wind components u and v from the real
    !   geographical system to the rotated system.
    !
    ! Method:
    !   Transformation formulas for converting between these two systems.
    !
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Parameter list:
    REAL (KIND=wp), INTENT (IN)          ::    &
         u  !< wind components in the true geographical system
    REAL (KIND=wp), INTENT (IN)          ::    &
         v    !< wind components in the true geographical system
    REAL (KIND=wp), INTENT (IN)          ::    &
         rlat !< coordinates in the true geographical system
    REAL (KIND=wp), INTENT (IN)          ::    &
         rlon !< coordinates in the true geographical system
    REAL (KIND=wp), INTENT (IN)          ::    &
         pollat !< latitude of the north pole of the rotated grid
    REAL (KIND=wp), INTENT (IN)          ::    &
         pollon    !< longitude of the north pole of the rotated grid

    REAL (KIND=wp), INTENT (OUT)         ::    &
         urot !< wind components in the rotated grid 
    REAL (KIND=wp), INTENT (OUT)         ::    &
         vrot        !< wind components in the rotated grid             

    ! Local variables

    REAL (KIND=wp)                       ::    &
         zsinpol, zcospol, zlonp, zlat, zarg1, zarg2, znorm

    REAL (KIND=wp)                       ::    &
         zpir18 = 0.0174532925_wp

    !------------------------------------------------------------------------------
    ! Begin Subroutine uv2uvrot
    !------------------------------------------------------------------------------

    zsinpol = SIN(pollat * zpir18)
    zcospol = COS(pollat * zpir18)
    zlonp   = (pollon-rlon) * zpir18
    zlat    =         rlat  * zpir18

    zarg1   = zcospol*SIN(zlonp)
    zarg2   = zsinpol*COS(zlat) - zcospol*SIN(zlat)*COS(zlonp)
    znorm   = 1.0_wp/SQRT( zarg1**2 + zarg2**2 )

    ! Transform the u and v wind components
    urot   =  u*zarg2*znorm - v*zarg1*znorm
    vrot   =  u*zarg1*znorm + v*zarg2*znorm

  END SUBROUTINE uv2uvrot

  !==============================================================================
  !==============================================================================

  !------------------------------------------------------------------------------
  !>  This routine converts the wind components u and v from the real  geographical system to the rotated system.
  !>  This is the vectorized form.
  SUBROUTINE uv2uvrot_vec(u, v, rlat, rlon, pollat, pollon, idim, jdim)

    !------------------------------------------------------------------------------
    !
    ! Description:
    !   This routine converts the wind components u and v from the real
    !   geographical system to the rotated system. This is the vectorized form
    !   of the routine above, i.e. the computation is for a whole 2D field.
    !
    ! Method:
    !   Transformation formulas for converting between these two systems.
    !
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    ! Parameter list:
    INTEGER, INTENT(IN)     ::    &
         idim        !< dimensions of the field
    INTEGER, INTENT(IN)     ::    &
         jdim        !< dimensions of the field

    REAL (KIND=wp), INTENT (INOUT)       ::    &
         u  (idim,jdim) !< wind components 
    REAL (KIND=wp), INTENT (INOUT)       ::    &
         v  (idim,jdim)  !< wind components 

    REAL (KIND=wp), INTENT (IN)          ::  rlat(idim,jdim)      !< latitude in the true geographical system
    REAL (KIND=wp), INTENT (IN)          ::  rlon (idim,jdim)     !< longitude in the true geographical system
    REAL (KIND=wp), INTENT (IN)          ::  pollat    !< latitude  of the north pole of the rotated grid
    REAL (KIND=wp), INTENT (IN)          ::  pollon    !<  longitude of the north pole of the rotated grid

    ! Local variables
    REAL (KIND=wp)                       ::    &
         zsinpol, zcospol, zlonp, zlat, zarg1, zarg2, znorm, zurot, zvrot

    INTEGER                 ::    i, j
    REAL (KIND=wp)                       ::    &
         zpir18 = 0.0174532925_wp

    !------------------------------------------------------------------------------
    ! Begin Subroutine uv2uvrot_vec
    !------------------------------------------------------------------------------

    zsinpol = SIN ( pollat * zpir18 )
    zcospol = COS ( pollat * zpir18 )

    DO j = 1, jdim
      DO i = 1, idim

        zlonp   = ( pollon - rlon(i,j) ) * zpir18
        zlat    =            rlat(i,j)   * zpir18

        zarg1   = zcospol*SIN(zlonp)
        zarg2   = zsinpol*COS(zlat) - zcospol*SIN(zlat)*COS(zlonp)
        znorm   = 1.0_wp/SQRT( zarg1**2 + zarg2**2 )

        ! Transform the u and v wind components
        zurot =  u(i,j)*zarg2*znorm - v(i,j)*zarg1*znorm
        zvrot =  u(i,j)*zarg1*znorm + v(i,j)*zarg2*znorm
        u(i,j) = zurot
        v(i,j) = zvrot

      ENDDO
    ENDDO

  END SUBROUTINE uv2uvrot_vec

  !==============================================================================
  !==============================================================================

  !------------------------------------------------------------------------------
  !> This routine computes wind speed and wind direction from the wind components.
  SUBROUTINE uv2df (u, v, d, f)

    !------------------------------------------------------------------------------
    !
    ! Description:
    !   This routine computes wind speed and wind direction from the wind
    !   components.
    !
    ! Method:
    !   Straightforward.
    !
    !------------------------------------------------------------------------------
    !
    ! Parameter list:
    REAL (KIND=wp), INTENT (IN)          ::    &
         u   !< wind components in the true geographical system

    REAL (KIND=wp), INTENT (IN)          ::    &
         v           !< wind components in the true geographical system

    REAL (KIND=wp), INTENT (OUT)         ::    &
         f   !< wind speed 

    REAL (KIND=wp), INTENT (OUT)         ::    &
         d   !< wind direction

    ! Local variables

    REAL (KIND=wp)                       ::    &

         zrpi18 = 57.2957795_wp,       & ! conversion from radians to degrees
         zsmall = 0.001_wp

    !------------------------------------------------------------------------------
    ! Begin Subroutine uv2df
    !------------------------------------------------------------------------------

    IF (ABS(u) > zsmall) THEN
      f  =  SQRT( u*u + v*v )
      d  =  v / u
      d  =  180.0_wp + SIGN( 90.0_wp , u ) - ATAN( d ) *zrpi18
    ELSEIF (ABS(v) > zsmall) THEN
      f  =  ABS( v )
      d  =  270.0_wp - SIGN( 90.0_wp , v )
    ELSE
      f  =    0.0_wp
      d  =    0.0_wp
    ENDIF

  END SUBROUTINE uv2df

  !==============================================================================
  !==============================================================================

  !------------------------------------------------------------------------------
  !> This routine computes wind speed and wind direction from the wind components. This is the vectorized form. 
  SUBROUTINE uv2df_vec (u, v, d, f, idim, jdim)

    !------------------------------------------------------------------------------
    !
    ! Description:
    !   This routine computes wind speed and wind direction from the wind
    !   components.  This is the vectorized form of the routine above,
    !   i.e. the computation is for a whole 2D field.
    !
    ! Method:
    !   Straightforward.
    !
    !------------------------------------------------------------------------------
    !
    ! Parameter list:
    INTEGER, INTENT(IN)     ::    &
         idim        !< dimensions of the field
    INTEGER, INTENT(IN)     ::    &
         jdim        !< dimensions of the field


    REAL (KIND=wp), INTENT (IN)          ::    &
         u  (idim,jdim)  !< wind components in the true geographical system
    REAL (KIND=wp), INTENT (IN)          ::    &
         v  (idim,jdim)  !< wind components in the true geographical system


    REAL (KIND=wp), INTENT (OUT)         ::    &
         f  (idim,jdim)  !< wind speed
    REAL (KIND=wp), INTENT (OUT)         ::    &
         d  (idim,jdim)    !< wind direction

    ! Local variables

    INTEGER                 ::    i, j
    REAL (KIND=wp)                       ::    &
         zrpi18 = 57.2957795_wp,       & ! conversion from radians to degrees
         zsmall = 0.001_wp

    !------------------------------------------------------------------------------
    ! Begin Subroutine uv2df_vec
    !------------------------------------------------------------------------------

    DO j = 1, jdim
      DO i = 1, idim

        IF (ABS(u(i,j)) > zsmall) THEN
          f (i,j)  =  SQRT( u(i,j)*u(i,j) + v(i,j)*v(i,j) )
          d (i,j)  =  180.0_wp + SIGN( 90.0_wp , u(i,j) )               &
               - ATAN( v(i,j) / u(i,j) ) *zrpi18
        ELSEIF (ABS(v(i,j)) > zsmall) THEN
          f (i,j)  =  ABS( v(i,j) )
          d (i,j)  =  270.0_wp - SIGN( 90.0_wp , v(i,j) )
        ELSE
          f (i,j)  =    0.0_wp
          d (i,j)  =    0.0_wp
        ENDIF

      ENDDO
    ENDDO

  END SUBROUTINE uv2df_vec

  !==============================================================================
  !==============================================================================
  !roa>
  !==============================================================================
  !------------------------------------------------------------------------------
  !> This subroutine extends a given 2D field with additional rows and columns

  SUBROUTINE extend_field2D ( field_in,  ie_in,  je_in,                          &
       field_out, ie_out, je_out,                           &
       nextlines)

    !------------------------------------------------------------------------------
    !
    ! Description:
    !   This routine extends field_in with nextlines additional rows and columns
    !   around. At the borders of the total domain, all values are set to the
    !   values of the neighboring row/column of field_in. 
    !
    ! Method:
    !
    !------------------------------------------------------------------------------

    ! Subroutine arguments:
    ! ---------------------
    INTEGER (KIND=i4), INTENT (IN) ::  &
         ie_in,  je_in,         & ! horizontal dimensions of field_in
         ie_out, je_out           ! dimensions of field_out

    INTEGER (KIND=i4), INTENT (IN) ::  &  
         nextlines                ! number of additional rows/columns at each side
    ! for the extension

    REAL    (KIND=wp  ), INTENT (IN) ::  &
         field_in (ie_in, je_in)

    REAL    (KIND=wp   ), INTENT(OUT) ::  &
         field_out(ie_out,je_out)

    ! Local scalars:
    ! -------------
    INTEGER (KIND=i4) ::  i, j

    !------------------------------------------------------------------------------

    field_out = 0.0

    ! fill domain existing in both field_in and field_out
    DO j = 1, je_in
      DO i = 1, ie_in
        field_out(i+nextlines,j+nextlines) = field_in(i,j)
      END DO
    END DO

    ! set values for the boundary lines in the extended frame
    ! western border
    DO j = 1, je_in
      DO i = 1, nextlines
        field_out(i,j+nextlines) = field_in(1,j)
      ENDDO
    ENDDO

    ! eastern border
    DO j = 1, je_in
      DO i = ie_in+1, ie_in+nextlines
        field_out(i+nextlines,j+nextlines) = field_in(ie_in,j)
      ENDDO
    ENDDO


    ! southern border
    DO j = 1, nextlines
      DO i = 1, ie_in
        field_out(i+nextlines,j) = field_in(i,1)
      ENDDO
    ENDDO

    ! northern border
    DO j = je_in+1, je_in+nextlines
      DO i = 1, ie_in
        field_out(i+nextlines,j+nextlines) = field_in(i,je_in)
      ENDDO
    ENDDO


    ! southwest corner
    DO j = 1, nextlines
      DO i = 1, nextlines
        field_out(i,j) = field_in(1,1)
      ENDDO
    ENDDO

    ! northwest corner
    DO j = je_in+1, je_in+nextlines
      DO i = 1, nextlines
        field_out(i,j+nextlines) = field_in(1,je_in)
      ENDDO
    ENDDO

    ! southeast corner
    DO j = 1, nextlines
      DO i = ie_in+1, ie_in+nextlines
        field_out(i+nextlines,j) = field_in(ie_in,1)
      ENDDO
    ENDDO

    ! northeast corner
    DO j = je_in+1, je_in+nextlines
      DO i = ie_in+1, ie_in+nextlines
        field_out(i+nextlines,j+nextlines) = field_in(ie_in,je_in)
      ENDDO
    ENDDO



  END SUBROUTINE extend_field2D

  !==============================================================================
  !==============================================================================
  !------------------------------------------------------------------------------
  !> This subroutine performs horizontal fitlering of a field, using a local
  !  filter. The stencil width can be 17-, 13-, 9-, or 3-points and the weights
  !  have been determined by J. Foerstner (DWD).

  SUBROUTINE horizontal_filtering( field_flt, ie_in, je_in,                 &
       nflt_width, ncutoff,                   &
       hfx_mask, hfy_mask )

    !------------------------------------------------------------------------------
    !
    ! Description:
    !
    ! Method:
    !
    !------------------------------------------------------------------------------

    ! Subroutine arguments:
    ! ---------------------
    INTEGER (KIND=i4), INTENT(IN) :: &
         ie_in, je_in            ! Dimensions of the field to be filtered

    INTEGER (KIND=i4), INTENT(IN) :: &
         nflt_width,           & ! width of field extension for filtering
         ncutoff                 ! filter-value for cutoff

    REAL    (KIND=wp   ), INTENT(INOUT) ::  &
         field_flt(ie_in, je_in)

    LOGICAL, INTENT(in), OPTIONAL ::  &
         hfx_mask(ie_in, je_in), hfy_mask(ie_in, je_in)

    ! Local scalars:
    ! -------------
    INTEGER (KIND=i4)  ::  &
         ilow, iup,           & !
         jlow, jup,           & !
         i, j,                & !  Loop indices
         iend,        &
         jend

    INTEGER (KIND=i4) ::  &
         l, nfw_m_nb

    ! Local (automatic) arrays:
    ! -------------------------
    REAL    (KIND=wp   ) ::  &
         field_tmp (ie_in, je_in), &
         zfwnp(-nflt_width:nflt_width),   & ! filter weights for n-point filter
         zfw3p(-1:1)                        ! filter weights for 3-point filter

    !------------------------------------------------------------------------------

    nfw_m_nb = nflt_width 

    iend   = ie_in - 2*nfw_m_nb 
    jend   = je_in - 2*nfw_m_nb 

    ! filter weights for n-point filter
    IF (ncutoff == 3 .AND. nflt_width == 4) THEN
      ! --> dfilt4
      ! filter weights for 9-point filter (approx. cutoff = 3)
      zfwnp = (/ -0.390625E-02,     &
           +0.3125E-01,       &
           -0.109375,         &
           +0.21875,          &
           +0.7265625,        &
           +0.21875,          &
           -0.109375,         &
           +0.3125E-01,       &
           -0.390625E-02 /)
    ELSEIF (ncutoff == 3 .AND. nflt_width == 8) THEN
      ! --> dfilt8
      ! filter weights for 17-point filter (approx. cutoff = 3)
      zfwnp = (/ -0.15259E-04,      &
           +0.2441406E-03,    &
           -0.18310546E-02,   &
           +0.85449218E-02,   &
           -0.27770996E-01,   &
           +0.666503906E-01,  &
           -0.1221923828,     &
           +0.1745605469,     &
           +0.8036193848,     &
           +0.1745605469,     &
           -0.1221923828,     &
           +0.666503906E-01,  &
           -0.27770996E-01,   &
           +0.85449218E-02,   &
           -0.18310546E-02,   &
           +0.2441406E-03,    &
           -0.15259E-04 /)
    ELSEIF (ncutoff == 4 .AND. nflt_width == 4) THEN
      ! filter weights for 9-point filter (approx. cutoff = 4)
      zfwnp = (/ +0.1171875E-01,    &
           -0.3125E-01,       &
           -0.46875E-01,      &
           +0.28125,          &
           +0.5703125,        &
           +0.28125,          &
           -0.46875E-01,      &
           -0.3125E-01,       &
           +0.1171875E-01 /)
    ELSEIF (ncutoff == 5 .AND. nflt_width == 6) THEN
      ! filter weights for 13-point filter (approx. cutoff = 5)
      zfwnp = (/ +0.44023278E-02,   &
           +0.13175894E-01,   &
           -0.477203075E-01,  &
           -0.435555245E-01,  &
           +0.94700467E-01,   &
           +0.2888298641,     &
           +0.3803345582,     &
           +0.2888298641,     &
           +0.94700467E-01,   &
           -0.435555245E-01,  &
           -0.477203075E-01,  &
           +0.13175894E-01,   &
           +0.44023278E-02 /)
    ELSEIF (ncutoff == 6 .AND. nflt_width == 4) THEN
      ! filter weights for 9-point filter (approx. cutoff = 6)
      zfwnp = (/ -0.4694126E-01,    &
           -0.50095541E-02,   &
           +0.13528415,       &
           +0.25500955,       &
           +0.32331423,       &
           +0.25500955,       &
           +0.13528415,       &
           -0.50095541E-02,   &
           -0.4694126E-01 /)
    ELSEIF (ncutoff == 8 .AND. nflt_width == 6) THEN
      ! filter weights for 13-point filter (approx. cutoff = 8)
      zfwnp = (/ -0.16638111E-01,   &
           -0.30753028E-01,   &
           -0.17361869E-02,   &
           +0.65428931E-01,   &
           +0.14784805,       &
           +0.2153241,        &
           +0.2410525,        &
           +0.2153241,        &
           +0.14784805,       &
           +0.65428931E-01,   &
           -0.17361869E-02,   &
           -0.30753028E-01,   &
           -0.16638111E-01 /)
    ELSE
      PRINT *, ' ERROR *** Wrong cutoff value for filtering        or *** '
      PRINT *, ' ERROR *** wrong value for filter/field extension.    *** '
    ENDIF

    ! filter weights for 3-point filter (approx. cutoff = 4)
    zfw3p = (/ 0.25, 0.5, 0.25 /)

    ! west
    ilow = 1 + 2*nflt_width

    ! east
    iup = iend 

    ! south
    jlow = 1 + 2*nflt_width

    ! north
    jup = jend 


    ! init working array
    field_tmp (:,:) = field_flt(:,:)


    IF ( PRESENT( hfx_mask ) ) THEN

      ! apply n-point-filter in x-direction
      DO j = 1, je_in
        DO i = ilow, iup
          IF ( hfx_mask(i,j) ) THEN
            field_tmp(i,j) = 0.0
          ENDIF
        ENDDO
        DO l = -nflt_width, nflt_width
          DO i = ilow, iup
            IF ( hfx_mask(i,j) ) THEN
              field_tmp(i,j) = field_tmp(i,j)               &
                   + zfwnp(l)*field_flt(i+l,j)
            END IF
          END DO
        END DO
      END DO

      ! apply 3-point-filter in x-direction at west boundary
      DO j = 1, je_in
        DO i = nfw_m_nb+1, ilow-1
          IF ( hfx_mask(i,j) ) THEN
            field_tmp(i,j) = 0.0
          ENDIF
        ENDDO
        DO l = -1, 1
          DO i = nfw_m_nb+1, ilow-1
            IF ( hfx_mask(i,j) ) THEN
              field_tmp(i,j) = field_tmp(i,j)             &
                   + zfw3p(l)*field_flt(i+l,j)
            END IF
          END DO
        END DO
      END DO

      ! apply 3-point-filter in x-direction at east boundary
      DO j = 1, je_in
        DO i = iup+1, ie_in-nfw_m_nb
          IF ( hfx_mask(i,j) ) THEN
            field_tmp(i,j) = 0.0
          ENDIF
        ENDDO
        DO l = -1, 1
          DO i = iup+1, ie_in-nfw_m_nb
            IF ( hfx_mask(i,j) ) THEN
              field_tmp(i,j) = field_tmp(i,j)             &
                   + zfw3p(l)*field_flt(i+l,j)
            END IF
          END DO
        END DO
      END DO

    ELSE

      !
      ! apply n-point-filter in x-direction
      !
      DO j = 1, je_in
        DO i = ilow, iup
          field_tmp(i,j) = 0.0
        END DO
        DO l = -nflt_width, nflt_width
          DO i = ilow, iup
            field_tmp(i,j) = field_tmp(i,j)                 &
                 + zfwnp(l)*field_flt(i+l,j)
          END DO
        END DO
      END DO

      ! apply 3-point-filter in x-direction at west boundary
      DO j = 1, je_in
        DO i = nfw_m_nb+1, ilow-1
          field_tmp(i,j) = 0.0
        END DO
        DO l = -1, 1
          DO i = nfw_m_nb+1, ilow-1
            field_tmp(i,j) = field_tmp(i,j)               &
                 + zfw3p(l)*field_flt(i+l,j)
          END DO
        END DO
      END DO

      ! apply 3-point-filter in x-direction at east boundary
      DO j = 1, je_in
        DO i = iup+1, ie_in-nfw_m_nb
          field_tmp(i,j) = 0.0
        END DO
        DO l = -1, 1
          DO i = iup+1, ie_in-nfw_m_nb
            field_tmp(i,j) = field_tmp(i,j)               &
                 + zfw3p(l)*field_flt(i+l,j)
          END DO
        END DO
      END DO

    END IF


    IF ( PRESENT( hfy_mask ) ) THEN

      ! apply n-point-filter in y-direction
      DO j = jlow, jup
        DO i = 1, ie_in
          IF ( hfy_mask(i,j) ) THEN
            field_flt(i,j) = 0.0
          ELSE
            field_flt(i,j) = field_tmp(i,j)
          ENDIF
        ENDDO
        DO l = -nflt_width, nflt_width
          DO i = 1, ie_in
            IF ( hfy_mask(i,j) ) THEN
              field_flt(i,j) = field_flt(i,j) + zfwnp(l)*field_tmp(i,j+l)
            END IF
          END DO
        END DO
      END DO

      ! apply 3-point-filter in y-direction at south boundary
      DO j = nfw_m_nb+1, jlow-1
        DO i = 1, ie_in
          IF ( hfy_mask(i,j) ) THEN
            field_flt(i,j) = 0.0
          ELSE
            field_flt(i,j) = field_tmp(i,j)
          ENDIF
        ENDDO
        DO l = -1, 1
          DO i = 1, ie_in
            IF ( hfy_mask(i,j) ) THEN
              field_flt(i,j) = field_flt(i,j)+zfw3p(l)*field_tmp(i,j+l)
            END IF
          END DO
        END DO
      END DO

      ! apply 3-point-filter in y-direction at north boundary
      DO j = jup+1, je_in-nfw_m_nb
        DO i = 1, ie_in
          IF ( hfy_mask(i,j) ) THEN
            field_flt(i,j) = 0.0
          ELSE
            field_flt(i,j) = field_tmp(i,j)
          ENDIF
        ENDDO
        DO l = -1, 1
          DO i = 1, ie_in
            IF ( hfy_mask(i,j) ) THEN
              field_flt(i,j) = field_flt(i,j)+zfw3p(l)*field_tmp(i,j+l)
            END IF
          END DO
        END DO
      END DO

    ELSE

      !
      ! apply n-point-filter in y-direction
      !
      DO j = jlow, jup
        DO i = 1, ie_in
          field_flt(i,j) = 0.0
        ENDDO
        DO l = -nflt_width, nflt_width
          DO i = 1, ie_in
            field_flt(i,j) = field_flt(i,j)+zfwnp(l)*field_tmp(i,j+l)
          END DO
        END DO
      END DO

      ! apply 3-point-filter in y-direction at south boundary
      DO j = nfw_m_nb+1, jlow-1
        DO i = 1, ie_in
          field_flt(i,j) = 0.0
        ENDDO
        DO l = -1, 1
          DO i = 1, ie_in
            field_flt(i,j) = field_flt(i,j)+zfw3p(l)*field_tmp(i,j+l)
          END DO
        END DO
      END DO

      ! apply 3-point-filter in y-direction at north boundary
      DO j = jup+1, je_in-nfw_m_nb
        DO i = 1, ie_in
          field_flt(i,j) = 0.0
        ENDDO
        DO l = -1, 1
          DO i = 1, ie_in
            field_flt(i,j) = field_flt(i,j)+zfw3p(l)*field_tmp(i,j+l)
          END DO
        END DO
      END DO

    END IF

  END SUBROUTINE horizontal_filtering

  !==============================================================================
  !roa<
  !==============================================================================

END MODULE mo_utilities_extpar
