!+  Fortran module to find grid element index in GME grid
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
!
! Code Description:
! Language: Fortran 95.
!====================================================================!>
!! Fortran module to find grid element index in GME grid
!!
!> Note: this program is used to generate external parameters for the GME. 
!>       It originates from old code for the GME grid from
!>       B. Ritter, D. Majewski and other DWD staff members.
MODULE mo_search_gme_grid

  USE mo_kind,            ONLY: wp, i4, i8
  USE mo_io_units,        ONLY: filename_max
  USE mo_grid_structures, ONLY: gme_triangular_grid
  USE mo_grid_structures, ONLY: target_grid_def
  USE mo_grid_structures, ONLY: igrid_gme

  USE mo_utilities_extpar, ONLY: abort_extpar

  USE mo_gme_grid, ONLY: factorni


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pp_ll2gp

  CONTAINS
    !> find nearest grid point in GME grid
    SUBROUTINE pp_ll2gp(xn,slon,slat,nip1,zx,zy,zz,spd_t,kd,kj1,kj2,sp,ldebug)
    ! Description:
    !     *pp_ll2gp* finds the nearest grid point on the unit sphere of
    !                the triangular grid for a point of specified
    !                longitude (slon) and latitude (slat)
    ! Note: this program is used to generate external parameters for the GME. 
    !       It originates from old code from
    !       B. Ritter, A. Mueller, J. Helmert and other DWD staff members.
    ! Current Code Owner: DWD, Hermann Asensio
    !  phone:  +49  69  8062 2709
    !  fax:    +49  69  8062 3721
    !  email:  Hermann.Asensio@dwd.de 
    !
    ! Code Description:
    ! Language: Fortran 90.
    !=======================================================================
    ! Declarations:
    !
    IMPLICIT NONE

    !     Input:
    !
    INTEGER                  :: nip1  ! grid mesh dimension
    REAL (KIND=wp)           :: xn(nip1,nip1,3,10)   !< cartesian coordinates of all nodes
    !                                   ATTENTION: In the GME "xn" is
    !                                   defined as xn(0:ni, 1:ni+1,3,10)
    REAL (KIND=wp)           :: slon  ! longitude (degrees) of point     
    REAL (KIND=wp)           :: slat  ! latitude  (degrees) of point  

    REAL (KIND=wp)           :: spd_t ! threshold value for scalar product 
                                   ! for initiation of additional fine
                                   ! scale search with routine ZU_FUSS
    !     Output:
    !
    REAL (KIND=wp)            :: zx,zy,zz
    REAL (KIND=wp)            :: sp         ! scalar product between point and
                                  ! nearest GME nodal point
    INTEGER                  ::  kd    ! diamond containing point
    INTEGER                  ::  kj1,kj2 
    !                                  nodal indices of nearest grid point
    !                                  on entry, kj1 and kj2 are first guess values

    INTEGER                  :: itest, icalls, n, ni, n2, n3, ierr, &
                                n3p1, n2pn3, jd, kpa, j1top, j2top,  &
                                j1ll, j2ll, j1lr, j2lr, jj1top,  &
                                jj2top, jj1ll, jj2ll, jj1lr, jj2lr, & 
                                ji, j1, j2, j1_ini, j2_ini

    INTEGER                  ::  mfirst = 0 ! utility variable
    INTEGER                  ::  m_error! error return code for search
    INTEGER                  ::  method = 2 ! choice of search method
                                        ! 1 = subdivision of triangles
                                        ! 2 = subdivison of diamond

    REAL (KIND=wp)                     :: api, d2r,zlon, zlat, smax, zxpt, &
                                zypt, zzpt, zxat, zyat, zzat, spp, & 
                                spa, sptop, spll, splr, sp_ini

    LOGICAL ldebug

    SAVE api,d2r,n2,n3
    ! DATA ldebug/.false./
    !
    !=======================================================================
    !
    !     Presettings
    !
    n     = nip1
    ni    = nip1-1

    if (mfirst.eq.0) then
    api   = 2.*asin(1.)
    d2r   = api/180.
    CALL  factorni(ni,ldebug,n2,n3,ierr)
    mfirst = 1
    end if

    zlon  = slon*d2r       
    zlat  = slat*d2r       
    zx    = cos(zlon)*cos(zlat)
    zy    = sin(zlon)*cos(zlat)
    zz    =           sin(zlat)
    !!$#ifdef DEBUG
    !!$print *,'=== pp_ll2gp.f in DEBUG MODE ==='
    !!$print *,'=== searching for zlon: ',slon    
    !!$print *,'===               zlat: ',slat    
    !!$print *,'===               zx  : ',zx      
    !!$print *,'===               zy  : ',zy      
    !!$print *,'===               zz  : ',zz      
    !!$#endif


    n3p1  = n3+1
    n2pn3 = n2+n3

    smax  = -999.
    kd    = 0

    !!$#ifdef DEBUG
    !!$print *,'=== diamond search loop      '
    !!$print *,'=== nip1=',nip1
    !!$print *,'=== n   =',n
    !!$#endif
    DO jd = 1,10
      zxpt  = xn(1,1,1,jd)+xn(1,n,1,jd)+xn(n,1,1,jd)    ! Zentrum des pol-
      zypt  = xn(1,1,2,jd)+xn(1,n,2,jd)+xn(n,1,2,jd)    ! waertigen Icosa-
      zzpt  = xn(1,1,3,jd)+xn(1,n,3,jd)+xn(n,1,3,jd)    ! dreiecks
      zxat  = xn(n,n,1,jd)+xn(1,n,1,jd)+xn(n,1,1,jd)    ! Zentrum des anti-
      zyat  = xn(n,n,2,jd)+xn(1,n,2,jd)+xn(n,1,2,jd)    ! polwaertigen
      zzat  = xn(n,n,3,jd)+xn(1,n,3,jd)+xn(n,1,3,jd)    ! Icosaederdreiecks
      spp   = zx*zxpt + zy*zypt + zz*zzpt   
      spa   = zx*zxat + zy*zyat + zz*zzat   
      !!$#ifdef DEBUG
      !!$print *,'=== diamond : ',jd
      !!$print *,'=== spp     : ',spp
      !!$print *,'=== spa     : ',spa
      !!$print *,'=== smax    : ',smax
      !!$print *,'=== xn(1,1,1,jd): ',xn(1,1,1,jd)
      !!$print *,'=== xn(1,1,2,jd): ',xn(1,1,2,jd)
      !!$print *,'=== xn(1,1,3,jd): ',xn(1,1,3,jd)
      !!$!     print *,'=== zxpt    : ',zxpt
      !!$!     print *,'=== zypt    : ',zypt
      !!$!     print *,'=== zzpt    : ',zzpt
      !!$!     print *,'=== zxat    : ',zxat
      !!$!     print *,'=== zyat    : ',zyat
      !!$!     print *,'=== zzat    : ',zzat
      !!$#endif

      IF (spp > smax ) THEN
        kd    = jd
        kpa   =  1
        smax  = spp
      END IF
      IF (spa > smax ) THEN
        kd    = jd
        kpa   = -1
        smax  = spa
      END IF
    END DO

    IF (method == 1) then  ! search via subdivision of triangles
    !!$#ifdef DEBUG
    !!$print *,'=== search method no. 1'          
    !!$#endif

    IF (kpa == 1) THEN     ! Point in poleward triangle   

    !!$#ifdef DEBUG
    !!$print *,'=== point belongs to diamond:',kd
    !!$print *,'=== in poleward triangle     '
    !!$#endif
    !     Binaersuche nach Unterdreieck der naechsten Stufe
      j1top = 1     ! Zeilenindex der polwaertigen Dreiecksspitze
      j2top = 1     ! Spaltenindex der polwaertigen Dreiecksspitze
      j1ll  = n     ! Zeilenindex der linken unteren Dreiecksspitze
      j2ll  = 1     ! Spaltenindex der linken unteren Dreiecksspitze
      j1lr  = 1     ! Zeilenindex der rechten unteren Dreiecksspitze
      j2lr  = n     ! Spaltenindex der rechten unteren Dreiecksspitze
    !
    ELSE 
    !
    !!$#ifdef DEBUG
    !!$print *,'=== point belongs to diamond:',kd
    !!$print *,'=== in antipoleward triangle     '
    !!$#endif
      j1top = n     ! Zeilenindex der antipolwaertigen Dreiecksspitze
      j2top = n     ! Spaltenindex der antipolwaertigen Dreiecksspitze
      j1ll  = 1     ! Zeilenindex der linken unteren Dreiecksspitze
      j2ll  = n     ! Spaltenindex der linken unteren Dreiecksspitze
      j1lr  = n     ! Zeilenindex der rechten unteren Dreiecksspitze
      j2lr  = 1     ! Spaltenindex der rechten unteren Dreiecksspitze
    END IF

    IF (n3 == 1) THEN
    !!$#ifdef DEBUG
    !!$print *,'=== calling sub_t9 ===='            
    !!$print *,'=== old indices  ===='            
    !!$print *,'=== j1top: ',j1top,'  j2top:',j2top
    !!$print *,'=== j1ll : ',j1ll ,'  j2ll :',j2ll 
    !!$print *,'=== j1lr : ',j1lr ,'  j2lr :',j2lr 
    !!$#endif
          
      CALL sub_t9(xn,n,kd,zx,zy,zz, &
                  j1top ,j2top  , &
                  j1ll  ,j2ll   , &
                  j1lr  ,j2lr   ,ldebug, &
                  jj1top,jj2top , &
                  jj1ll ,jj2ll  , &
                  jj1lr ,jj2lr )
         
    !     Swap: sub triangle corners ---> new main triangle corners
      j1top = jj1top
      j1ll  = jj1ll
      j1lr  = jj1lr
      j2top = jj2top
      j2ll  = jj2ll
      j2lr  = jj2lr
    !!$#ifdef DEBUG
    !!$      print *,'=== new indices  ===='            
    !!$      print *,'=== j1top: ',j1top,'  j2top:',j2top
    !!$      print *,'=== j1ll : ',j1ll ,'  j2ll :',j2ll 
    !!$      print *,'=== j1lr : ',j1lr ,'  j2lr :',j2lr 
    !!$#endif

    END IF

    !!$#ifdef DEBUG
    !!$      print *,'=== loop for calling sub_t4 ===='            
    !!$      print *,'=== old indices  ===='            
    !!$      print *,'=== j1top: ',j1top,'  j2top:',j2top
    !!$      print *,'=== j1ll : ',j1ll ,'  j2ll :',j2ll 
    !!$      print *,'=== j1lr : ',j1lr ,'  j2lr :',j2lr 
    !!$#endif
    DO ji=n3p1,n2pn3     ! Anzahl der Unterteilungsniveaus
    !!$#ifdef DEBUG
    !!$      print *,'=== loop index :',ji
    !!$#endif
          
    CALL sub_t4(xn,n,kd,zx,zy,zz, &
                        j1top ,j2top  , &
                        j1ll  ,j2ll   , &
                        j1lr  ,j2lr   ,ldebug, &
                        jj1top,jj2top , &
                        jj1ll ,jj2ll  , &
                        jj1lr ,jj2lr )
         
    !     Swap: sub triangle corners ---> new main triangle corners
    j1top = jj1top
    j1ll  = jj1ll
    j1lr  = jj1lr
    j2top = jj2top
    j2ll  = jj2ll
    j2lr  = jj2lr
    !!$#ifdef DEBUG
    !!$print *,'=== new indices  ===='            
    !!$print *,'=== j1top: ',j1top,'  j2top:',j2top
    !!$print *,'=== j1ll : ',j1ll ,'  j2ll :',j2ll 
    !!$print *,'=== j1lr : ',j1lr ,'  j2lr :',j2lr 
    !!$#endif
    ENDDO 

    !     At end of subdivision sequence find nearest grid point:

    sptop = zx*xn(j1top,j2top,1,kd)+zy*xn(j1top,j2top,2,kd)+zz*xn(j1top,j2top,3,kd)
    spll  = zx*xn(j1ll ,j2ll ,1,kd)+zy*xn(j1ll ,j2ll ,2,kd)+zz*xn(j1ll ,j2ll ,3,kd)
    splr  = zx*xn(j1lr ,j2lr ,1,kd)+zy*xn(j1lr ,j2lr ,2,kd)+zz*xn(j1lr ,j2lr ,3,kd)
    !!$#ifdef DEBUG
    !!$      print *,'=== end of subdivision loop  ===='            
    !!$      print *,'=== sptop: ',sptop
    !!$      print *,'=== spll : ',spll 
    !!$      print *,'=== splr : ',splr 
    !!$      print *,'=== j1top: ',j1top,'  j2top:',j2top
    !!$      print *,'=== j1ll : ',j1ll ,'  j2ll :',j2ll 
    !!$      print *,'=== j1lr : ',j1lr ,'  j2lr :',j2lr 
    !!$#endif

    IF (sptop >= spll .AND. sptop >= splr) THEN
      j1 = j1top
      j2 = j2top
      sp = sptop
    ELSE IF (spll >= splr) THEN
      j1 = j1ll
      j2 = j2ll
      sp = spll 
    ELSE
      j1 = j1lr
      j2 = j2lr
      sp = splr 
    END IF
    !!$#ifdef DEBUG
    !!$      print *,'=== nearest grid point ===='            
    !!$      print *,'=== sp: ',sp
    !!$      print *,'=== j1: ',j1,'  j2:',j2
    !!$#endif
    !
    !=======================================================================
    !
    IF (sp.lt.spd_t) THEN
    !!$#ifdef DEBUG
    !!$      print *,'=== final walk ===='            
    !!$      print *,'=== spd_t: ',spd_t
    !!$#endif
      j1_ini = j1
      j2_ini = j2
      sp_ini = sp
      CALL final_walk (xn,n,kd,zx,zy,zz,j1_ini,j2_ini,sp_ini,ldebug,j1,j2,sp)
    !!$#ifdef DEBUG
    !!$      print *,'=== after final walk   ===='            
    !!$      print *,'=== sp: ',sp
    !!$      print *,'=== j1: ',j1,'  j2:',j2
    !!$#endif
    END IF

    ELSE       ! use binary search method

    !!$#ifdef DEBUG
    !!$      print *,'=== search method no. 2'          
    !!$#endif
    j1 = kj1+1   ! first guess value for first index
    j2 = kj2     ! first guess value for sec.  index
    m_error = 0  ! preset error return code for search
    CALL  R_TO_G(zx,zy,zz,            & ! cartesian coordinates of pixel
      &            xn  ,              & ! target grid cartesian coordinates
      &            n   ,kd,  j1  ,j2 ,& ! target grid resolution and point
      &            m_error,ldebug)      ! error code
    if (m_error.ne.0) then
          print *,'ERROR error ERROR error ERROR'
          print *,'Search of target grid point failed'
          stop 'search'
    end if
    sp = zx*xn(j1,j2,1,kd)+zy*xn(j1,j2,2,kd)+zz*xn(j1,j2,3,kd)

    END IF     ! choice of search method

    kj1 = j1-1  ! indices in calling routine     
    kj2 = j2    ! indices in calling routine
    !
    !=======================================================================
    !
    END SUBROUTINE pp_ll2gp

    !> determines the grid point indices of a GME nodal point which is closest to a specified point on the unit sphere
    SUBROUTINE final_walk (xn,n,kd,zx,zy,zz,j1_ini,j2_ini,sp_ini,ldebug,jj1_1,jj2_1,sp_1)
    ! Description:
    !
    !     *final_walk*  determines the grid point indices of a GME nodal point
    !     which is closest to a specified point on the unit sphere
    !
    !
    ! Note: this program is used to generate external parameters for the GME. 
    !       It originates from old code from
    !       B. Ritter, A. Mueller, J. Helmert and other DWD staff members.
    !       In a future version all code related to external parameters should be 
    !       ported into a new version using Fortran Modules.
    ! Current Code Owner: DWD, Hermann Asensio
    !  phone:  +49  69  8062 2709
    !  fax:    +49  69  8062 3721
    !  email:  Hermann.Asensio@dwd.de 
    !
    ! Code Description:
    ! Language: Fortran 90.
    !=======================================================================

    !
    IMPLICIT NONE

    INTEGER :: n, kd, j1_ini, j2_ini, jj1_1, jj2_1,  &
               j1_o, j2_o, jstart, jloop, jdir, j1_n, j2_n
    REAL (KIND=wp)  :: sp_ini, sp_1, zx, zy, zz, sp_o, sp_n
    REAL (KIND=wp)  :: xn(n,n,3,10)   ! cartesian coordinates of all grid points

    INTEGER :: j1n(11),j2n(11) ! index increments for relevant neigbours
                               ! of any GME nodal point
                               ! as each node has only six neighbours,
                               ! elements 7 to 11 are just repetitions
                               ! of elements 1 to 5 to allow each neigbour
                               ! (from 1 to 6) as starting point
                               !  in a search loop, cf. Fig.1:
                            !
                            !      5---------4
                            !     / \       / \            \       /
                            !    /   \     /   \            \     /
                            !   /     \   /     \            \   /
                            !  /       \ /       \            \ /
                            ! 6---------*---------3          j1,j2
                            !  \       / \       /            / \
                            !   \     /   \     /            /   \
                            !    \   /     \   /            /     \
                            !     \ /       \ /            /       \
                            !      1---------2           j1+1     j2+1


    REAL (KIND=wp)  :: SPD  (6)       ! scalar product of neigbour nodes with reference point
    LOGICAL ldebug
    !
    !=======================================================================
    !
    !     Spoke No.  1  2  3  4  5  6    1  2  3  4  5
    !     ==============================================
    data j1n / 1, 0,-1,-1, 0, 1,   1, 0,-1,-1, 0 / 
    data j2n / 0, 1, 1, 0,-1,-1,   0, 1, 1, 0,-1 /
    !
    !=======================================================================
    !
    j1_o  = j1_ini                   ! 1.index of first guess
    j2_o  = j2_ini                   ! 2.index of first guess
    sp_o  = sp_ini                   ! first guess scalar product
    jstart  = 1                        ! initial search direction

    !     Walking loop                                      
    !     ============

    DO 2000 jloop=1,100000

    !     Directional loop
    !     ----------------

     DO 1000 jdir=jstart,jstart+5
     j1_n = j1_o + j1n(jdir)
     j2_n = j2_o + j2n(jdir)

     IF (j1_n >= 1 .AND. j1_n <= n .AND.j2_n >= 1 .and. j2_n <= n) THEN
     sp_n=zx*xn(j1_n,j2_n,1,kd)+zy*xn(j1_n,j2_n,2,kd)+zz*xn(j1_n,j2_n,3,kd)

        IF (sp_n > sp_o) THEN      ! improvement !
          j1_o = j1_n
          j2_o = j2_n
          sp_o = sp_n
          jstart = jdir-((jdir-1)/6) * 6  ! try same direction first
          GO TO 1001
        END IF                         ! improvement
      END IF                           ! domain boundaries
     1000   CONTINUE                         ! loop over directions

    !      If this part of the code is reached, no improvement could be
    !      achieved in any direction, i.e. the fit is optimal
    !      ------------------------------------------------------------

      GO TO 2001

     1001   CONTINUE            ! branch adress for new, improved point

     2000 CONTINUE              ! end of walking loop           

    !     If this part of the code is reached, the walk was not completed
    !     within the walking loop 
    !     ----------------------------------------------------------------

    STOP '  Error in *final_walk*, walk not completed'

     2001 CONTINUE            ! best fit branch address      


    !     Store final result in dummy arguments
    !     -------------------------------------
    !
    jj1_1 = j1_o       
    jj2_1 = j2_o       
    sp_1  = sp_o       

    !=======================================================================
    !
    END SUBROUTINE final_walk

    !
    !>  determines grid point indices of a subtriangle
    SUBROUTINE sub_t9(xn,n,kd,zx,zy,zz, &
                           j1top,j2top, &
                           j1ll ,j2ll , &
                           j1lr ,j2lr ,ldebug, &
                           jj1top,jj2top, &
                           jj1ll ,jj2ll , &
                           jj1lr ,jj2lr )
    ! Description:
    !     *sub_t9* determines the grid point indices of a subtriangle
    !              which contains a specified point on the unit sphere
    !              The sub triangle is one of nine which partition a
    !              known main triangle, defined by the three pairs of
    !              grid point indices of it's corners (top, lower left,
    !              lower right). The main triangle must contain the
    !              the point which is searched for!
    !
    !
    !
    ! Note: this program is used to generate external parameters for the GME. 
    !       It originates from old code from
    !       B. Ritter, A. Mueller, J. Helmert and other DWD staff members.
    !       In a future version all code related to external parameters should be 
    !       ported into a new version using Fortran Modules.
    ! Current Code Owner: DWD, Hermann Asensio
    !  phone:  +49  69  8062 2709
    !  fax:    +49  69  8062 3721
    !  email:  Hermann.Asensio@dwd.de 
    !


    ! ---------- ---------- ----
    ! !VERSION!  !DATE!     B. Ritter
    !  Initial release
    !
    !

    INTEGER   :: j1(7),j2(7) ! indices of new subtriangle corner points
    INTEGER   :: n, kd, j1top, j2top, j1ll, j2ll, &
                 j1lr, j2lr, jj1top, jj2top, jj1ll, &
                 jj2ll, jj1lr, jj2lr, jd, ide
    REAL (KIND=wp)      ::       &
         xn(n,n,3,10), &  ! cartesian coordinates of all grid points
         xtc(9),       &  ! x-coordinate of subtriangle centers 
         ytc(9),       &  ! y-coordinate of subtriangle centers  
         ztc(9)           ! z-coordinate of subtriangle centers 

    REAL (KIND=wp)                     :: zx, zy, zz, spmin, rnorm, sp

    logical ldebug
    !
    !=======================================================================
    !
    !     Structure of main triangle and the points defining the subtriangles

    !
    !     a) polward triangle             j1top,j2top
    !                                         / \
    !                                        /   \
    !                                       /     \
    !                                      /       \
    !                                     /         \
    !                                    /    D1     \
    !                                   /             \
    !                                  /               \
    !                                 /                 \
    !                           j1(1),j2(1)---------j1(6),j2(6)
    !                               / \                 / \
    !                              /   \               /   \
    !                             /     \     D7      /     \
    !                            /       \           /       \
    !                           /         \         /         \
    !                          /    D2     \       /     D6    \
    !                         /             \     /             \
    !                        /               \   /               \
    !                       /                 \ /                 \
    !                 j1(2),j2(2)---------j1(7),j2(7)---------j1(5),j2(5)
    !                     / \                 / \                 / \
    !                    /   \               /   \               /   \
    !                   /     \     D8      /     \      D9     /     \
    !                  /       \           /       \           /       \
    !                 /         \         /         \         /         \
    !                /    D3     \       /    D4     \       /    D5     \
    !               /             \     /             \     /             \
    !              /               \   /               \   /               \
    !             /                 \ /                 \ /                 \
    !        j1ll,j2ll----------j1(3),j2(3)---------j1(4),j2(4)----------j1lr,j2lr
    !
    !
    !     b) antipolward triangle: analogous to case a), but top is oriented
    !                              away from pole and lower left refers to
    !                              'eastern' corner, lower right refers to
    !                              'western' corner
    !
    !=======================================================================
    !
    !     Subindices
    j1(1) = j1top +   (j1ll-j1top)/3
    j2(1) = j2top +   (j2ll-j2top)/3
    j1(2) = j1top + 2*(j1ll-j1top)/3
    j2(2) = j2top + 2*(j2ll-j2top)/3
    j1(3) = j1top + 2*(j1ll-j1top)/3
    j2(3) = j2top +   (j2lr-j2top)/3
    j1(4) = j1top +   (j1ll-j1top)/3
    j2(4) = j2top + 2*(j2lr-j2top)/3
    j1(5) = j1top + 2*(j1lr-j1top)/3
    j2(5) = j2top + 2*(j2lr-j2top)/3
    j1(6) = j1top +   (j1lr-j1top)/3
    j2(6) = j2top +   (j2lr-j2top)/3
    j1(7) = (j1(3)+j1(6))/2          
    j2(7) = (j2(1)+j2(4))/2          

    !     Centres of the 9 sub-triangles

    xtc(1)  = xn(j1top,j2top,1,kd)+xn(j1(1),j2(1),1,kd)+xn(j1(6),j2(6),1,kd)
    ytc(1)  = xn(j1top,j2top,2,kd)+xn(j1(1),j2(1),2,kd)+xn(j1(6),j2(6),2,kd)
    ztc(1)  = xn(j1top,j2top,3,kd)+xn(j1(1),j2(1),3,kd)+xn(j1(6),j2(6),3,kd)

    xtc(2)  = xn(j1(1),j2(1),1,kd)+xn(j1(2),j2(2),1,kd) +xn(j1(7),j2(7),1,kd)
    ytc(2)  = xn(j1(1),j2(1),2,kd)+xn(j1(2),j2(2),2,kd) +xn(j1(7),j2(7),2,kd)
    ztc(2)  = xn(j1(1),j2(1),3,kd)+xn(j1(2),j2(2),3,kd) +xn(j1(7),j2(7),3,kd)

    xtc(3)  = xn(j1(3),j2(3),1,kd)+xn(j1(2),j2(2),1,kd)+xn(j1ll ,j2ll ,1,kd)
    ytc(3)  = xn(j1(3),j2(3),2,kd)+xn(j1(2),j2(2),2,kd)+xn(j1ll ,j2ll ,2,kd)
    ztc(3)  = xn(j1(3),j2(3),3,kd)+xn(j1(2),j2(2),3,kd)+xn(j1ll ,j2ll ,3,kd)

    xtc(4)  = xn(j1(7),j2(7),1,kd)+xn(j1(3),j2(3),1,kd)+xn(j1(4),j2(4),1,kd)
    ytc(4)  = xn(j1(7),j2(7),2,kd)+xn(j1(3),j2(3),2,kd)+xn(j1(4),j2(4),2,kd)
    ztc(4)  = xn(j1(7),j2(7),3,kd)+xn(j1(3),j2(3),3,kd)+xn(j1(4),j2(4),3,kd)

    xtc(5)  = xn(j1(5),j2(5),1,kd)+xn(j1(4),j2(4),1,kd)+xn(j1lr ,j2lr ,1,kd)
    ytc(5)  = xn(j1(5),j2(5),2,kd)+xn(j1(4),j2(4),2,kd)+xn(j1lr ,j2lr ,2,kd)
    ztc(5)  = xn(j1(5),j2(5),3,kd)+xn(j1(4),j2(4),3,kd)+xn(j1lr ,j2lr ,3,kd)

    xtc(6)  = xn(j1(6),j2(6),1,kd)+xn(j1(7),j2(7),1,kd)+xn(j1(5),j2(5),1,kd)
    ytc(6)  = xn(j1(6),j2(6),2,kd)+xn(j1(7),j2(7),2,kd)+xn(j1(5),j2(5),2,kd)
    ztc(6)  = xn(j1(6),j2(6),3,kd)+xn(j1(7),j2(7),3,kd)+xn(j1(5),j2(5),3,kd)

    xtc(7)  = xn(j1(1),j2(1),1,kd)+xn(j1(7),j2(7),1,kd)+xn(j1(6),j2(6),1,kd)
    ytc(7)  = xn(j1(1),j2(1),2,kd)+xn(j1(7),j2(7),2,kd)+xn(j1(6),j2(6),2,kd)
    ztc(7)  = xn(j1(1),j2(1),3,kd)+xn(j1(7),j2(7),3,kd)+xn(j1(6),j2(6),3,kd)

    xtc(8)  = xn(j1(2),j2(2),1,kd)+xn(j1(3),j2(3),1,kd)+xn(j1(7),j2(7),1,kd)
    ytc(8)  = xn(j1(2),j2(2),2,kd)+xn(j1(3),j2(3),2,kd)+xn(j1(7),j2(7),2,kd)
    ztc(8)  = xn(j1(2),j2(2),3,kd)+xn(j1(3),j2(3),3,kd)+xn(j1(7),j2(7),3,kd)

    xtc(9)  = xn(j1(7),j2(7),1,kd)+xn(j1(4),j2(4),1,kd)+xn(j1(5),j2(5),1,kd)
    ytc(9)  = xn(j1(7),j2(7),2,kd)+xn(j1(4),j2(4),2,kd)+xn(j1(5),j2(5),2,kd)
    ztc(9)  = xn(j1(7),j2(7),3,kd)+xn(j1(4),j2(4),3,kd)+xn(j1(5),j2(5),3,kd)

    spmin = 0.0

    DO jd = 1,9
    rnorm = 1./sqrt(xtc(jd)**2+ytc(jd)**2+ztc(jd)**2)
    xtc(jd) = rnorm*xtc(jd)     
    ytc(jd) = rnorm*ytc(jd)     
    ztc(jd) = rnorm*ztc(jd)     
    sp      = xtc(jd)*zx+ytc(jd)*zy+ztc(jd)*zz
    IF (sp > spmin) THEN
      ide     = jd
      spmin   = sp
    END IF
    END DO
          
    IF (ide == 1) THEN      ! Point in D1
      jj1top = j1top
      jj2top = j2top
      jj1ll  = j1(1)
      jj2ll  = j2(1)
      jj1lr  = j1(6)
      jj2lr  = j2(6)
    ELSE IF (ide == 2) THEN ! Point in D2
      jj1top = j1(1)
      jj2top = j2(1)
      jj1ll  = j1(2)
      jj2ll  = j2(2)
      jj1lr  = j1(7)
      jj2lr  = j2(7)
    ELSE IF (ide == 3) THEN ! Point in D3
      jj1top = j1(2)
      jj2top = j2(2)
      jj1ll  = j1ll
      jj2ll  = j2ll
      jj1lr  = j1(3)
      jj2lr  = j2(3)
    ELSE IF (ide == 4) THEN ! Point in D4
      jj1top = j1(7)
      jj2top = j2(7)
      jj1ll  = j1(3)
      jj2ll  = j2(3)
      jj1lr  = j1(4)
      jj2lr  = j2(4)
    ELSE IF (ide == 5) THEN ! Point in D5
      jj1top = j1(5)
      jj2top = j2(5)
      jj1ll  = j1(4)
      jj2ll  = j2(4)
      jj1lr  = j1lr  
      jj2lr  = j2lr  
    ELSE IF (ide == 6) THEN ! Point in D6
      jj1top = j1(6)
      jj2top = j2(6)
      jj1ll  = j1(7)
      jj2ll  = j2(7)
      jj1lr  = j1(5)
      jj2lr  = j2(5)
    ELSE IF (ide == 7) THEN ! Point in D7
      jj1top = j1(7)
      jj2top = j2(7)
      jj1ll  = j1(6)
      jj2ll  = j2(6)
      jj1lr  = j1(1)
      jj2lr  = j2(1)
    ELSE IF (ide == 8) THEN ! Point in D8
      jj1top = j1(3)
      jj2top = j2(3)
      jj1ll  = j1(7)
      jj2ll  = j2(7)
      jj1lr  = j1(2)
      jj2lr  = j2(2)
    ELSE IF (ide == 9) THEN ! Point in D9
      jj1top = j1(4)
      jj2top = j2(4)
      jj1ll  = j1(5)
      jj2ll  = j2(5)
      jj1lr  = j1(7)
      jj2lr  = j2(7)
    ELSE                    ! ERROR
      PRINT *,' Error in SUB_t9 ide=',ide
      STOP 'Error in SUB_t9'
    END IF
    !
    !=======================================================================
    END SUBROUTINE sub_t9

    !> determines grid point indices  of a subtriangle
    SUBROUTINE sub_t4(xn,n,kd,zx,zy,zz, &
                           j1top,j2top, &
                           j1ll ,j2ll , &
                           j1lr ,j2lr ,ldebug, &
                           jj1top,jj2top, &
                           jj1ll ,jj2ll , &
                           jj1lr ,jj2lr )
    ! Description:
    !     *sub_t4* determines the grid point indices of a subtriangle
    !              which contains a specified point on the unit sphere
    !              The sub triangle is one of four which partition a
    !              known main triangle, defined by the three pairs of
    !              grid point indices of it's corners (top, lower left,
    !              lower right). The main triangle must contain the
    !              the point which is searched for!

    !
    !
    ! Note: this program is used to generate external parameters for the GME. 
    !       It originates from old code from
    !       B. Ritter, A. Mueller, J. Helmert and other DWD staff members.
    !       In a future version all code related to external parameters should be 
    !       ported into a new version using Fortran Modules.
    ! Current Code Owner: DWD, Hermann Asensio
    !  phone:  +49  69  8062 2709
    !  fax:    +49  69  8062 3721
    !  email:  Hermann.Asensio@dwd.de 
    !
    ! ---------- ---------- ----
    ! !VERSION!  !DATE!     B. Ritter
    !  Initial release
    !
    !
    IMPLICIT NONE

    INTEGER          :: n, kd, j1top, j2top, j1ll, j2ll, &
                        j1lr, j2lr, jj1top, jj2top, jj1ll, &
                        jj2ll, jj1lr, jj2lr, jd, ide
    REAL (KIND=wp)                     :: &
         xn(n,n,3,10)   ! cartesian coordinates of all grid points
    INTEGER                   :: &
            j1(3),j2(3) ! indices of new subtriangle corner points
    REAL (KIND=wp)                     :: &
         xtc(4),  &     ! x-coordinate of subtriangle centers 
         ytc(4),  &     ! y-coordinate of subtriangle centers  
         ztc(4)         ! z-coordinate of subtriangle centers 
    LOGICAL ldebug

    REAL (KIND=wp)            :: zx, zy, zz, spmin, rnorm, sp

    !=======================================================================
    !     Structure of main triangle and the points defining the subtriangles


    !     a) polward triangle             j1top,j2top
    !                                         / \
    !                                        /   \
    !                                       /     \
    !                                      /       \
    !                                     /         \
    !                                    /    D1     \
    !                                   /             \
    !                                  /               \
    !                                 /                 \
    !                       j1(1),j2(1)-----------------j1(3),j2(3)
    !                               / \                 / \
    !                              /   \               /   \
    !                             /     \    D4       /     \
    !                            /       \           /       \
    !                           /         \         /         \
    !                          /    D2     \       /     D3    \
    !                         /             \     /             \
    !                        /               \   /               \
    !                       /                 \ /                 \
    !                  j1ll,j2ll----------j1(2),j2(2)----------j1lr,j2lr
    !
    !
    !     b) antipolward triangle: analogous to case a), but top is oriented
    !                              away from pole and lower left refers to
    !                              'eastern' corner, lower right refers to
    !                              'western' corner
    !
    !
    !=======================================================================
    !
    !     Subindices
    j1(1) = j1top + 0.5*(j1ll-j1top)
    j2(1) = j2top + 0.5*(j2ll-j2top) 
    j1(2) = j1top + 0.5*(j1ll-j1top)
    j2(2) = j2top + 0.5*(j2lr-j2top) 
    j1(3) = j1top + 0.5*(j1lr-j1top)
    j2(3) = j2top + 0.5*(j2lr-j2top) 
    !
    !     Centres of the 4 sub-triangles
    !
    xtc(1)  = xn(j1top,j2top,1,kd)+xn(j1(1),j2(1),1,kd)+xn(j1(3),j2(3),1,kd)
    ytc(1)  = xn(j1top,j2top,2,kd)+xn(j1(1),j2(1),2,kd)+xn(j1(3),j2(3),2,kd)
    ztc(1)  = xn(j1top,j2top,3,kd)+xn(j1(1),j2(1),3,kd)+xn(j1(3),j2(3),3,kd)
    !
    xtc(2)  = xn(j1(1),j2(1),1,kd)+xn(j1ll ,j2ll ,1,kd)+xn(j1(2),j2(2),1,kd)
    ytc(2)  = xn(j1(1),j2(1),2,kd)+xn(j1ll ,j2ll ,2,kd)+xn(j1(2),j2(2),2,kd)
    ztc(2)  = xn(j1(1),j2(1),3,kd)+xn(j1ll ,j2ll ,3,kd)+xn(j1(2),j2(2),3,kd)

    xtc(3)  = xn(j1(3),j2(3),1,kd)+xn(j1(2),j2(2),1,kd)+xn(j1lr ,j2lr ,1,kd)
    ytc(3)  = xn(j1(3),j2(3),2,kd)+xn(j1(2),j2(2),2,kd)+xn(j1lr ,j2lr ,2,kd)
    ztc(3)  = xn(j1(3),j2(3),3,kd)+xn(j1(2),j2(2),3,kd)+xn(j1lr ,j2lr ,3,kd)

    xtc(4)  = xn(j1(2),j2(2),1,kd)+xn(j1(1),j2(1),1,kd)+xn(j1(3),j2(3),1,kd)
    ytc(4)  = xn(j1(2),j2(2),2,kd)+xn(j1(1),j2(1),2,kd)+xn(j1(3),j2(3),2,kd)
    ztc(4)  = xn(j1(2),j2(2),3,kd)+xn(j1(1),j2(1),3,kd)+xn(j1(3),j2(3),3,kd)

    spmin = 0.0

    DO jd = 1,4
    rnorm = 1./sqrt(xtc(jd)**2+ytc(jd)**2+ztc(jd)**2)
    xtc(jd) = rnorm*xtc(jd)     
    ytc(jd) = rnorm*ytc(jd)     
    ztc(jd) = rnorm*ztc(jd)     
    sp      = xtc(jd)*zx+ytc(jd)*zy+ztc(jd)*zz
    IF (sp > spmin) THEN
      ide     = jd
      spmin   = sp
    END IF
    END DO
          
    IF (ide == 1) THEN                       ! Point in D1
      jj1top = j1top
      jj2top = j2top
      jj1ll  = j1(1)
      jj2ll  = j2(1)
      jj1lr  = j1(3)
      jj2lr  = j2(3)
    ELSE IF (ide == 2) THEN  ! Point in D2
      jj1top = j1(1)
      jj2top = j2(1)
      jj1ll  = j1ll 
      jj2ll  = j2ll 
      jj1lr  = j1(2)
      jj2lr  = j2(2)
    ELSE IF (ide == 3) THEN ! Point in D3
      jj1top = j1(3)
      jj2top = j2(3)
      jj1ll  = j1(2)
      jj2ll  = j2(2)
      jj1lr  = j1lr
      jj2lr  = j2lr
    ELSE                                     ! Point in D4
      jj1top = j1(2)
      jj2top = j2(2)
      jj1ll  = j1(3)
      jj2ll  = j2(3)
      jj1lr  = j1(1)
      jj2lr  = j2(1)
    END IF
    !
    !=======================================================================
    END SUBROUTINE sub_t4

    !> associates pixels with cartesian coordinates px,py,pz on the unit sphere with the 'nearest' grid point of a target grid
    SUBROUTINE R_TO_G (px,py,pz,          &  ! cartesian coordinates of pixel
      &                 xn  ,               &  ! target grid cartesian coordinates
      &                 n   ,kd,  k1  ,k2  ,&  ! target grid resolution and point
      &                 k_error,ldebug)        ! error code and debug switch   
    ! Description:
    !
    !     ====================================================================
    !     R_TO_G  associates pixels with cartesian coordinates px,py,pz on the
    !             unit sphere with the 'nearest' grid point of a target grid
    !             specified by its cartesian coordinates
    !     ====================================================================
    !
    !
    ! Note: this program is used to generate external parameters for the GME. 
    !       It originates from old code from
    !       B. Ritter, A. Mueller, J. Helmert and other DWD staff members.
    !       In a future version all code related to external parameters should be 
    !       ported into a new version using Fortran Modules.
    ! Current Code Owner: DWD, Hermann Asensio
    !  phone:  +49  69  8062 2709
    !  fax:    +49  69  8062 3721
    !  email:  Hermann.Asensio@dwd.de 

    !     Input arrays and variables:

          INTEGER, INTENT(IN)            :: n            ! ni+1 horizontal resolution
          REAL (KIND=wp), INTENT(IN) :: px,py,pz           ! pixel coordinates
          REAL (KIND=wp), INTENT(IN) :: xn(n,n,3,10)       ! cartesian coordinates of target grid

          INTEGER, INTENT(IN)            :: kd           ! diamond containg the target grid point
                                                         ! (computed in calling program!)
          INTEGER, INTENT(INOUT) :: k1,k2           ! nearest target grid point
          INTEGER, INTENT(INOUT) :: k_error         ! error code (0, if successful)
          LOGICAL, INTENT(INOUT) :: ldebug          ! debug switch


    !     Local arrays and variables

          REAL (KIND=wp)   ::  zsp,zsp_max          ! scalar products
          REAL (KIND=wp)   :: zs1,zs2,zs3          ! scalar products
    !     real zsp,zsp_max          ! scalar products
          INTEGER :: jnb               ! loop variable ( 8 neighbours )

          integer jn1,jn2           ! 'neighbour' indices in target grid
          integer nbstep            ! step size for shell search
          INTEGER :: iter              ! iteration count


          logical lstartf           ! indicator for fine/coarse mode at
                                    ! start of module 

          integer mspoke1(8)  ! offset of internal neighbours from point in 1.index
          integer mspoke2(8)  ! offset of internal neighbours from point in 2.index
          data mspoke1 / 0, 1, 1, 1, 0,-1,-1,-1/
          data mspoke2 / 1, 1, 0,-1,-1,-1, 0, 1/


          lstartf=.true.            ! indicate, that search shall start in
                                    ! fine scale mode, i.e. step size 1

    !      if (k_error > 0) ldebug = .true.
          if (ldebug) then
          print *,'=========================================='
          print *,'========= ldebug ========================='
          print *,'R_TO_GR_TO_GR_TO_GR_TO_GR_TO_GR_TO_GR_TO_G'
          print *,'n  = ',n 
          print *,'k1 = ',k1
          print *,'k2 = ',k2
          print *,'kd = ',kd
          print *,'px = ',px
          print *,'py = ',py
          print *,'pz = ',pz
          print *,'x  = ',xn(k1,k2,1,kd)
          print *,'y  = ',xn(k1,k2,2,kd)
          print *,'z  = ',xn(k1,k2,3,kd)
          end if
    !     scalar product for first guess point
          zs1     = px*xn(k1,k2,1,kd)
          zs2     = py*xn(k1,k2,2,kd)
          zs3     = pz*xn(k1,k2,3,kd)
          zsp_max = zs1+zs2+zs3
          if (ldebug) then
          print *,'=========================================='
          print *,'initial zsp_max = ',zsp_max
          end if            

       10 continue              ! branch adress for complete rerun of module,
                                ! if 'fine scale search' was initially used
                                ! and sucessfull

          nbstep=n              ! initiate step width for 'coarse' start
          if (lstartf) nbstep=2 ! initiate step width for 'fine'   start

           if (ldebug) then
           print *,'============================================'
           print *,'after 10 continue: current step size=',nbstep
           end if            

      100 continue              ! branch adress for halving of step width
          iter  = 0
          nbstep=max(1,nbstep/2)

           if (ldebug) then
           print *,'============================================'
           print *,'after 100 continue: current step size=',nbstep
           end if            

      200 continue               ! branch adress for continued search
          iter  = iter + 1       ! counter for iterations

           if (ldebug) then
           print *,'============================================'
           print *,'after 200 continue: iteration count =',iter  
           end if            


    !     if (ldebug) then
    !     print *,'=========================================='
    !     print *,'nbstep,iter     = ',nbstep,iter    
    !     end if            

          if (iter.gt.1000000) then ! abort, if search does not converge
          print *,'=============================================='
          print *,'= R_TO_G can not find a suitable target grid ='
          print *,'=        point for this pixel                ='
          print *,'=============================================='
          print *,'= current target node k1,k2 = ',k1,k2
          print *,'= location                x = ',xn(k1,k2,1,kd)
          print *,'=                         y = ',xn(k1,k2,2,kd)
          print *,'=                         z = ',xn(k1,k2,3,kd)
          print *,'= pixel                  px = ',px 
          print *,'=                        py = ',py 
          print *,'=                        pz = ',pz 
          print *,'= best scalar product spmax = ',zsp_max 
          print *,'=============================================='
          k_error = 1
          return
          end if


          do jnb=1,8
          jn1= k1 + mspoke1(jnb)*nbstep
          jn1= max(1,min(n,jn1))
          jn2= k2 + mspoke2(jnb)*nbstep
          jn2= max(1,min(n,jn2))
          zs1     = px*xn(jn1,jn2,1,kd)
          zs2     = py*xn(jn1,jn2,2,kd)
          zs3     = pz*xn(jn1,jn2,3,kd)
          zsp     = zs2+zs3+zs1

           if (ldebug) then
           print *,'================================================'
           print *,'within neighbour loop of current central point:',k1,k2,kd
           print *,'neighbour indices                             :',jn1,jn2,' neighbour no.',jnb
           print *,'scalar product with search location           :',zsp                         
           print *,'current maximum scalar product                :',zsp_max                         
           print *,'x-neighbour:',xn(jn1,jn2,1,kd),' x-searched   :',px                          
           print *,'y-neighbour:',xn(jn1,jn2,2,kd),' y-searched   :',py                          
           print *,'z-neighbour:',xn(jn1,jn2,3,kd),' z-searched   :',pz                          
           end if            

           if (zsp > zsp_max) then  ! 'neighbour' is closer to pixel
           zsp_max = zsp
           k1      = jn1
           k2      = jn2
            if (ldebug) then 
            print *,'     improved location indices found:',k1,k2,' new maximum scalar product:',zsp_max
            endif

            if (lstartf.and.iter >  1) then ! increase search radius, because
            lstartf = .false.               ! 'fine' scale search caused 
            go to 10                        ! several index changes, when
            end if                          ! starting from first guess

           go to 200                 ! continue search for closer points
           end if

          end do

      if (nbstep > 1) go to 100 ! reduce step width, if not yet one

    !     best point found

    k_error = 0

    END SUBROUTINE R_TO_G

!=======================================================================
END MODULE mo_search_gme_grid

