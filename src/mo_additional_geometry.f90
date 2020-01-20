!+ Fortran module with additional geometry routines
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2010/12/21 Hermann Asensio
!  Initial release
! V1_2         2011/03/25 Hermann Asensio
!  add point in polygon test for convex polygons (like grid elements)
!  new subroutine point_in_grid_element
!
! Code Description:
! Language: Fortran 2003.
!=======================================================================
!> Fortran module with additional geometry routines
!!
!!
!! - point_in_polygon_sp
!! - point_in_grid_element
!!
!! @author Hermann Asensio, DWD
!!
!!
MODULE mo_additional_geometry

  USE mo_kind,            ONLY: wp
  USE mo_math_constants,  ONLY: pi_2, dbl_eps
  USE mo_base_geometry,   ONLY: geographical_coordinates
  USE mo_base_geometry,   ONLY: cartesian_coordinates

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cc2gc
  PUBLIC :: gc2cc
  PUBLIC :: vector_product
  PUBLIC :: vecpro
  PUBLIC :: cos_arc_length
  PUBLIC :: arc_length
  PUBLIC :: scal_pro
  PUBLIC :: chord_d
  PUBLIC :: inter_section

  PUBLIC :: point_in_grid_element

  PUBLIC :: point_in_polygon_sp

  CONTAINS

  !> Check if a point is inside a grid element (convex polygon) for ICON grid
  !! The coordinates ar given in carteseian coordinates
  SUBROUTINE point_in_grid_element(point,nvertices,cc_vertices,inflag)
    TYPE(cartesian_coordinates), INTENT(IN) :: point        !< cartesian coordinates of the "test point"
    INTEGER, INTENT(IN)                     :: nvertices    !< number of vertices of the grid element
    TYPE(cartesian_coordinates), INTENT(IN) :: cc_vertices(1:nvertices) !< cartesian coordinates of the vertices

    INTEGER, INTENT(OUT)                    :: inflag    !< flag for polygon test,
                                                         !! 0 if point is "outside",
                                                         !! 1 if point is "inside",
                                                         !! 2 if point is on the edge,
                                                         !! -1 for an error
    ! local variables


    TYPE(cartesian_coordinates) :: p1
    TYPE(cartesian_coordinates) :: v0,v1,v2
    INTEGER :: nv !< counter
    INTEGER :: ps !< position of dimension of maximal value of normal vector
    REAL (KIND=wp) :: sn(nvertices)
    TYPE(cartesian_coordinates) :: rn(nvertices)

    inflag = -1 ! set to "error"

    v1%x = cc_vertices(2)%x - cc_vertices(1)%x
    v2%x = cc_vertices(nvertices)%x - cc_vertices(1)%x
    p1 = vector_product(v1,v2)
    ps = MAXLOC(ABS(p1%x),DIM=1) ! this is the dimension (x1, x2 or x3) with the maximum value for a 
                                 ! vector rectangular to two edges of the grid element use this dimension "ps"
                                 ! as projection surface of the coordinates for a point in polygon test
   ! instead of the vector product and the test for equal sign in the "projection" surface
   ! one could perform a 3D "point in tetrahedron" test with a triple product and test for equal sign
    DO nv=1,nvertices-1
      v0%x = point%x - cc_vertices(nv)%x
      v1%x = cc_vertices(nv+1)%x - cc_vertices(nv)%x
      rn(nv) = vector_product(v1,v0)
      sn(nv) = SIGN(1.0_wp,rn(nv)%x(ps))
    ENDDO
    nv = nvertices ! last vertex
    v0%x = point%x - cc_vertices(nvertices)%x
    v1%x = cc_vertices(1)%x - cc_vertices(nvertices)%x
    rn(nv) = vector_product(v1,v0)
    sn(nv) = SIGN(1.0_wp,rn(nv)%x(ps))

    ! Test if point is on the edge of grid element
    DO nv=1,nvertices
      IF (ABS(rn(nv)%x(ps))== 0.0_wp) THEN ! point is on the edge of grid element
        inflag = 2
        RETURN
      ENDIF
    ENDDO

    ! Test if point is inside the edge with an "point in convex polygon test"
    DO nv=2,nvertices
      IF ((sn(nv) == sn(1))) THEN ! the sign of sn(nv) and sn(1) are the same if the point is within the convex polygon
        inflag = 1
      ELSE
        inflag = 0
        RETURN ! exit if the signs are different, which means that the point is outside the polygon
      ENDIF
    ENDDO


  END SUBROUTINE

  !> Check if a point is inside a polyon on a sphere
  !!
  !! This subroutine perfoms a point in polygon test on a sphere,
  !! where strictly speaking "inside" and "outside" is not really well defined.
  !! The polygon is assumed to be "simply connected", i.e. the edges do not intesect itself
  !! (no "8") and that there are no wholes. A test point defines the "outside" area.
  !! The intesections of the segment "test point" to the "target point" 
  !! (arc on the great circle which is defined by the two points)
  !! with the edges of the polygon are counted. 
  !! For an even number of interesctions (0, 2, 4, ...) the test point is in the "outside" area, inflag = 0
  !! for an odd number if intersections (1, 3, 5, ...) the test point is in the "inside" area. inflag = 1
  !! The input polygon is given with an array of the coordinates of the vertices,
  !! with dimension of polygon(nvertices+1), and with polygon(1) = polygon(nvertices+1).
  SUBROUTINE point_in_polygon_sp(point,nvertices, polygon, out_point, inflag)
    TYPE(cartesian_coordinates), INTENT(IN) :: point        !< coordinates of the "test point"
    INTEGER, INTENT(IN)                     :: nvertices    !< number of vertices of the polygon
    TYPE(cartesian_coordinates), INTENT(IN) :: polygon(1:nvertices+1) 
    !< coordinates of the "simple connected" polygon, and with  polygon(1) = polygon(nvertices+1)
    TYPE(cartesian_coordinates), INTENT(IN) :: out_point !< coordinates of a point which defines the "outside" area

    INTEGER, INTENT(OUT)                    :: inflag    
    !< flag for polygon test, 0 if point is "outside", 1 if point is "inside", -1 for an error

    ! local variables
    INTEGER :: crossn  !< number of intersections
    
    INTEGER :: i       !< counter

    TYPE(cartesian_coordinates)  :: isp !< intersection point

    LOGICAL :: isp_on_testarc ! flag if intersection point is within test arc 
    LOGICAL :: isp_on_polyedge! flag if intersection point is within polygon edge


    REAL(wp) :: c_arcl_pt_op      ! cos arc length (distance measure) of point - out_point arc
    REAL(wp) :: c_arcl_pt_isp     ! cos arc length (distance measure) of point - intersection point arc
    REAL(wp) :: c_arcl_poly_edge  ! cos arc length (distance measure) of polygon edge
    REAL(wp) :: c_arcl_poly1_ips  ! cos arc length (distance measure) of point polygon(i)- intersection point arc

    TYPE(cartesian_coordinates)  :: normal_pt_op     !  normal vector 
    TYPE(cartesian_coordinates)  :: normal_pt_isp    !  normal vector
    TYPE(cartesian_coordinates)  :: normal_poly_edge !  normal vector
    TYPE(cartesian_coordinates)  :: normal_poly1_ips !  normal vector
    
    REAL(KIND=wp) :: s_rot ! help variable for sense of rotation



    inflag = 0 ! set to a default value
    crossn = 0 ! start with crossn = 0

   allvertices: DO i=1,nvertices

     ! first deal with some pathological cases

     IF ((out_point%x(1) == polygon(i)%x(1)).AND.    &
         (out_point%x(2) == polygon(i)%x(2)).AND.    &
         (out_point%x(3) == polygon(i)%x(3)) ) THEN   !  -> outpoint not outside :-( -> inflag = -1 (error) or 
                                                      !  0 (outside) and exit subroutine
       inflag = -1 ! Error
       RETURN
     ENDIF

     
     IF ((out_point%x(1) == polygon(i+1)%x(1)).AND.    &
         (out_point%x(2) == polygon(i+1)%x(2)).AND.    &
         (out_point%x(3) == polygon(i+1)%x(3)) ) THEN !  -> outpoint not outside :-( -> inflag = -1 (error) or 
                                                      !  0 (outside) and exit subroutine
       inflag = -1 ! Error
       RETURN
     ENDIF

     IF ((point%x(1) == polygon(i)%x(1)).AND.    &
         (point%x(2) == polygon(i)%x(2)).AND.    &
         (point%x(3) == polygon(i)%x(3)) ) THEN   !  target point on polygon vertex -> inside

       inflag = 1 
       RETURN
     ENDIF


     IF ((point%x(1) == polygon(i+1)%x(1)).AND.    &
         (point%x(2) == polygon(i+1)%x(2)).AND.    &
         (point%x(3) == polygon(i+1)%x(3)) ) THEN !  target point on polygon vertex -> inside

       inflag = 1 
       RETURN
     ENDIF


     ! check if the test arc point - out_point interesects with polygon edge
     !----------------------------------------------------------------------------------------------------------------


     ! calculate some distance measures and normal vectors 
     ! calculate cos arc length (distance measure)
     c_arcl_pt_op     = cos_arc_length(point, out_point)
     normal_pt_op     = vector_product(point, out_point)  ! normal vector for plain point - (0,0,0) (earth centre) - out-point
     
     c_arcl_poly_edge = cos_arc_length(polygon(i), polygon(i+1))
     normal_poly_edge = vector_product(polygon(i), polygon(i+1)) ! normal vector for plain (polygon(i) - (0,0,0) 
                                                                 ! (earth centre) - (polygon(i+1)
     !----------------------------------------------------------------------------------------------------------------
     
     
     !  calculate the intesection points of the two great circles (given by point-test_point and polygon edge)
     !----------------------------------------------------------------------------------------------------------------
     isp =  inter_section (point, out_point, polygon(i), polygon(i+1))
     !isp2%x = -1. * isp%x ! on the sphere there are always two intersection points of great circles, 
                           ! here the sphere has the radius 1
                           ! the function "inter_section" returns the intersection point 
                           ! which is nearer to the first point (p0) in the argument list
                           ! of FUNCTION inter_section (p0, p1, v0, v1)
     !----------------------------------------------------------------------------------------------------------------
     
     ! pathological case
     !----------------------------------------------------------------------------------------------------------------
     IF (DOT_PRODUCT(isp%x, isp%x) == 0)  THEN ! all points are on one great cirlce, no intersection point
        ! all four points are on one plane in this case
        
        ! check if the target point is on polygon edge
         isp_on_polyedge = .false. ! default set to "false", i.e. the target point is not within the segment
     
         normal_poly1_ips = vector_product(polygon(i), point)
         s_rot = DOT_PRODUCT(normal_poly_edge%x,normal_poly1_ips%x) 
         IF (s_rot > 0 ) THEN  ! the sense of roation is not switched
          ! calculate cos arc length (distance measure)
           c_arcl_poly1_ips = cos_arc_length(polygon(i), point)
           IF (c_arcl_poly1_ips  >  c_arcl_poly_edge) THEN 
             !so the target point is within segment
               isp_on_polyedge = .true.
               inflag = 1 ! if the target point is on the polygon edge, define this as "inside" and exit subroutine
               RETURN
           ENDIF
         ENDIF


        ! check if the out_point is on polygon edge

        isp_on_polyedge = .false. ! default set to "false", i.e. the target point is not within the segment
     
         normal_poly1_ips = vector_product(polygon(i), out_point)
         s_rot = DOT_PRODUCT(normal_poly_edge%x,normal_poly1_ips%x) 
         IF (s_rot > 0 ) THEN  ! the sense of roation is not switched
          ! calculate cos arc length (distance measure)
           c_arcl_poly1_ips = cos_arc_length(polygon(i), out_point)
           IF (c_arcl_poly1_ips  >  c_arcl_poly_edge) THEN 
             !so the out point is within segment, ERROR, the out_point should not be on the polygon edge
               isp_on_polyedge = .false.
               inflag = -1 ! Error
               RETURN
           ENDIF
         ENDIF

         ! if the polygon edge is within test arc, do not count for crossing numbers
         ! crossn = crossn + 0   ! don't count for crossing number

         ! if the polygon edge is outside test arc, do not count for crossing numbers
         ! crossn = crossn + 0   ! don't count for crossing number

           CYCLE allvertices ! Continue with next polygon edge

     ENDIF  ! all points are on one great cirlce
     !----------------------------------------------------------------------------------------------------------------
     !----------------------------------------------------------------------------------------------------------------

     

     
     !----------------------------------------------------------------------------------------------------------------
     ! test if first intersection point is a crossing point
     ! test if isp point is within test arc point - out_point (check with sense of rotation and arc length)

     isp_on_testarc = .false. ! default set to "false", i.e. the interesection point is not within the segment
     ! calculate cos arc length (distance measure)
     c_arcl_pt_isp = cos_arc_length(point, isp)  ! if point==isp, cos_arc_length(point, isp) = 1
     IF (c_arcl_pt_isp == 1._wp) THEN  ! target point is the same as intersection point
        !so the intersection point is within segment
           isp_on_testarc = .true.
     ELSE
       normal_pt_isp = vector_product(point, isp)  
       s_rot = DOT_PRODUCT(normal_pt_op%x,normal_pt_isp%x) 
       IF (s_rot > 0 ) THEN  ! the sense of roation is not switched 
         IF (c_arcl_pt_isp > c_arcl_pt_op) THEN ! 
           !so the intersection point is within segment
             isp_on_testarc = .true.
         ENDIF
       ENDIF
     ENDIF
     !----------------------------------------------------------------------------------------------------------------

     ! test if isp point is within polygon edge
     isp_on_polyedge = .false. ! default set to "false", i.e. the interesection point is not within the segment

     ! calculate cos arc length (distance measure)
     c_arcl_poly1_ips = cos_arc_length(polygon(i), isp)

     IF  (c_arcl_poly1_ips == 1._wp)  THEN ! .OR. (c_arcl_poly1_ips==c_arcl_poly_edge) ) THEN 
                                           ! the intesection point is on a polygon vertex
       normal_poly1_ips = vector_product(isp, polygon(i+1))
       s_rot = DOT_PRODUCT(normal_pt_op%x,normal_poly1_ips%x) ! check if the normal vectors of the test arc and 
                                                              ! the polygon edge have the same direction
       IF (s_rot > 0 ) THEN ! define as "valid" intersection
          isp_on_polyedge = .true.
       ENDIF
       
     ELSEIF (c_arcl_poly1_ips == c_arcl_poly_edge) THEN !  the intesection point is on a polygon vertex
       normal_poly1_ips = vector_product(isp, polygon(i))
       s_rot = DOT_PRODUCT(normal_pt_op%x,normal_poly1_ips%x) ! check if the normal vectors of the test arc and 
                                                              ! the polygon edge have the same direction
       IF (s_rot > 0 ) THEN ! define as "valid" intersection
          isp_on_polyedge = .true.
       ENDIF

     ELSE
       normal_poly1_ips = vector_product(polygon(i), isp)
       s_rot = DOT_PRODUCT(normal_poly_edge%x,normal_poly1_ips%x) 
       IF (s_rot > 0 ) THEN  ! the sense of roation is not switched
         IF (c_arcl_poly1_ips  >  c_arcl_poly_edge) THEN 
           !so the intersection point is within segment
             isp_on_polyedge = .true.
             ! (c_arcl_poly1_ips==c_arcl_poly_edge)
         ENDIF
       ENDIF
     ENDIF

     IF (isp_on_testarc.AND.isp_on_polyedge) THEN
       crossn = crossn +1
     !----------------------------------------------------------------------------------------------------------------
     !----------------------------------------------------------------------------------------------------------------
     ELSE ! test if the second intersection point is a crossing point
       isp%x = -1. * isp%x ! the function "inter_section" returns the intersection point 
                           ! wich is nearer to the first point (p0) in the argument list
                           ! of FUNCTION inter_section (p0, p1, v0, v1) 
                           ! if the first arc (po-p1) has an angle > 90 degrees, 
                           ! the second intersection point of the great circles might be a crossing point
     !----------------------------------------------------------------------------------------------------------------
     !----------------------------------------------------------------------------------------------------------------
     !----------------------------------------------------------------------------------------------------------------
     !----------------------------------------------------------------------------------------------------------------
     !----------------------------------------------------------------------------------------------------------------
     !----------------------------------------------------------------------------------------------------------------
     !----------------------------------------------------------------------------------------------------------------
     !----------------------------------------------------------------------------------------------------------------


      ! test if isp point is within test arc point - out_point (check with sense of rotation and arc length)

       isp_on_testarc = .false. ! default set to "false", i.e. the interesection point is not within the segment
       ! calculate cos arc length (distance measure)
       c_arcl_pt_isp = cos_arc_length(point, isp)  ! if point==isp, cos_arc_length(point, isp) = 1
       IF (c_arcl_pt_isp == 1._wp) THEN  ! target point is the same as intersection point
          !so the intersection point is within segment
             isp_on_testarc = .true.
       ELSE
         normal_pt_isp = vector_product(point, isp)  
         s_rot = DOT_PRODUCT(normal_pt_op%x,normal_pt_isp%x) 
         IF (s_rot > 0 ) THEN  ! the sense of roation is not switched 
           IF (c_arcl_pt_isp > c_arcl_pt_op) THEN ! 
             !so the intersection point is within segment
               isp_on_testarc = .true.
           ENDIF
         ENDIF
       ENDIF
       !----------------------------------------------------------------------------------------------------------------

       ! test if isp point is within polygon edge
       isp_on_polyedge = .false. ! default set to "false", i.e. the interesection point is not within the segment

       ! calculate cos arc length (distance measure)
       c_arcl_poly1_ips = cos_arc_length(polygon(i), isp)

       IF  (c_arcl_poly1_ips == 1._wp)  THEN ! .OR. (c_arcl_poly1_ips==c_arcl_poly_edge) ) THEN 
                                             ! the intesection point is on a polygon vertex
         normal_poly1_ips = vector_product(isp, polygon(i+1))
         s_rot = DOT_PRODUCT(normal_pt_op%x,normal_poly1_ips%x) ! check if the normal vectors of the test arc and 
                                                                ! the polygon edge have the same direction
         IF (s_rot > 0 ) THEN ! define as "valid" intersection
            isp_on_polyedge = .true.
         ENDIF
         
       ELSEIF (c_arcl_poly1_ips == c_arcl_poly_edge) THEN !  the intesection point is on a polygon vertex
         normal_poly1_ips = vector_product(isp, polygon(i))
         s_rot = DOT_PRODUCT(normal_pt_op%x,normal_poly1_ips%x) ! check if the normal vectors of the test arc and 
                                                                ! the polygon edge have the same direction
         IF (s_rot > 0 ) THEN ! define as "valid" intersection
            isp_on_polyedge = .true.
         ENDIF

       ELSE
         normal_poly1_ips = vector_product(polygon(i), isp)
         s_rot = DOT_PRODUCT(normal_poly_edge%x,normal_poly1_ips%x) 
         IF (s_rot > 0 ) THEN  ! the sense of roation is not switched
           IF (c_arcl_poly1_ips  >  c_arcl_poly_edge) THEN 
             !so the intersection point is within segment
               isp_on_polyedge = .true.
               ! (c_arcl_poly1_ips==c_arcl_poly_edge)
           ENDIF
         ENDIF
       ENDIF

        IF (isp_on_testarc.AND.isp_on_polyedge) THEN
          crossn = crossn +1
        ENDIF

     ENDIF

   !----------------------------------------------------------------------------------------------------------------
   !----------------------------------------------------------------------------------------------------------------
   ENDDO allvertices

   if (MODULO(crossn,2).eq.1) then ! check for even or odd number of intersections
     inflag = 1
   else
     inflag = 0
   endif



  END  SUBROUTINE point_in_polygon_sp




!> Converts cartesian coordinates to geographical.
  ELEMENTAL FUNCTION cc2gc(x) RESULT (position)
!!
!! REVISION HISTORY:  
!! Developed  by Luis Kornblueh  (2004).
!! Completely new version by Thomas Heinze (2006-07-20)

! !INPUT PARAMETERS:  
    TYPE(cartesian_coordinates), INTENT(IN) :: x

! !RETURN VALUE:  
    TYPE(geographical_coordinates)          :: position

! !LOCAL VARIABLES:  
    REAL(wp)                                :: z_x, z_y, z_z, z_r

    z_x = x%x(1)
    z_y = x%x(2)
    z_z = x%x(3)

    z_r = z_x * z_x + z_y * z_y
    z_r = SQRT(z_r)
    
    IF (ABS(z_r) < dbl_eps) THEN    ! one of the poles

      IF (z_z > 0.0_wp) THEN
        position%lat = pi_2
      ELSE
        position%lat = -1._wp * pi_2
      END IF
      position%lon = 0._wp

    ELSE  

      position%lat = ATAN2 ( z_z, z_r) 

      IF (ABS(z_x) < dbl_eps) THEN    ! z_x == 0 ?

        IF (z_y >= 0.0_wp) THEN
          position%lon = pi_2
        ELSE
          position%lon = -1._wp * pi_2
        END IF

      ELSE
       position%lon = ATAN2( z_y, z_x)
      END IF

    END IF

  END FUNCTION  cc2gc

!> Converts longitude and latitude to cartesian coordinates.
  ELEMENTAL FUNCTION gc2cc(position) RESULT(x)
!! Converts longitude and latitude to cartesian coordinates.
!!
!! REVISION HISTORY:  
!! Developed  by Luis Kornblueh  (2004).
 
    TYPE(geographical_coordinates), INTENT(IN) :: position
    TYPE(cartesian_coordinates) :: x

    REAL (wp) :: z_cln, z_sln, z_clt, z_slt

    z_sln = SIN(position%lon)
    z_cln = COS(position%lon)
    z_slt = SIN(position%lat)
    z_clt = COS(position%lat)

    x%x(1) = z_cln*z_clt
    x%x(2) = z_sln*z_clt
    x%x(3) = z_slt

  END FUNCTION gc2cc


 !>  Computes vector product of x0,x1
  ELEMENTAL FUNCTION vector_product (x0, x1) RESULT(x2)
!!
!! REVISION HISTORY:  
!! Developed  by Luis Kornblueh  (2004).
!-----------------------------------------------------------------------  
    TYPE(cartesian_coordinates), INTENT(IN) :: x0, x1
    TYPE(cartesian_coordinates) :: x2

    x2%x(1) = x0%x(2)*x1%x(3) - x0%x(3)*x1%x(2)
    x2%x(2) = x0%x(3)*x1%x(1) - x0%x(1)*x1%x(3)
    x2%x(3) = x0%x(1)*x1%x(2) - x0%x(2)*x1%x(1)

  END FUNCTION vector_product


!>  Computes vector product of 2 unit vectors
ELEMENTAL FUNCTION  vecpro(geog1,geog2) RESULT (xx)
   
TYPE(geographical_coordinates), INTENT(IN) :: geog2, geog1

TYPE(cartesian_coordinates)  :: xx
REAL (wp)  :: long1, lat1,long2,lat2
REAL (wp) :: sn1,sn2,sn3,sn4,cs1,cs2

long1=geog1%lon
lat1=geog1%lat
long2=geog2%lon
lat2=geog2%lat
   
sn1=SIN(lat1-lat2)
sn2=SIN(lat1+lat2)
sn3=SIN(0.5_wp*(long1-long2))
sn4=SIN(0.5_wp*(long1+long2))
cs1=COS(0.5_wp*(long1-long2))
cs2=COS(0.5_wp*(long1+long2))
   
xx%x(1)=sn1*sn4*cs1-sn2*cs2*sn3
xx%x(2)=sn1*cs2*cs1+sn2*sn4*sn3
xx%x(3)=COS(lat1)*COS(lat2)*SIN(long1-long2)
    
END FUNCTION vecpro

!> Computes length of geodesic arc with endpoints x0,x1.
  ELEMENTAL FUNCTION arc_length (x0, x1) RESULT(arc)
!! DESCRIPTION:
!! Computes length of geodesic arc with endpoints x0,x1.
!!
!! !REVISION HISTORY:  
!! Developed by Th.Heinze (2006-09-19).
!! Previous version by Luis Kornblueh (2004) discarded.

    TYPE(cartesian_coordinates), INTENT(IN) :: x0, x1
    REAL(wp) :: arc
    REAL(wp) :: z_d0,  z_d1,  z_cc

    z_d0 = SQRT(DOT_PRODUCT(x0%x,x0%x))
    z_d1 = SQRT(DOT_PRODUCT(x1%x,x1%x))

    z_cc = DOT_PRODUCT(x0%x, x1%x)/(z_d0*z_d1)

! in case we get numerically incorrect solutions 
    IF (z_cc > 1._wp )  z_cc =  1._wp
    IF (z_cc < -1._wp ) z_cc = -1._wp

    arc = ACOS(z_cc)
 
  END FUNCTION arc_length


!>  Computes the cosine of the length of geodesic arc with endpoints x0,x1.
  ELEMENTAL FUNCTION cos_arc_length (x0, x1) RESULT(carc)
!! !DESCRIPTION:
!! Computes the cosine of the length of geodesic arc with endpoints x0,x1.
!!
!! !REVISION HISTORY:  
!! Developed by Almut Gassmann (2007-03-13).

    TYPE(cartesian_coordinates), INTENT(IN) :: x0, x1
    REAL(wp) :: carc
    REAL(wp) :: z_d0,  z_d1,  z_cc

    z_d0 = SQRT(DOT_PRODUCT(x0%x,x0%x))
    z_d1 = SQRT(DOT_PRODUCT(x1%x,x1%x))

    z_cc = DOT_PRODUCT(x0%x, x1%x)/(z_d0*z_d1)

! in case we get numerically incorrect solutions 
    IF (z_cc > 1._wp )  z_cc =  1._wp
    IF (z_cc < -1._wp ) z_cc = -1._wp

    carc = z_cc
 
  END FUNCTION cos_arc_length

    !> Computes the scalar product of to points in cartesian coordinates
  ELEMENTAL FUNCTION scal_pro(x0, x1) RESULT(carc)
! !REVISION HISTORY:  
! Developed by Hermann Asensio (2009).
    TYPE(cartesian_coordinates), INTENT(IN) :: x0, x1
    REAL(wp) :: carc

    carc = x0%x(1)*x1%x(1) + x0%x(2)*x1%x(2) + x0%x(3)*x1%x(3) 
 
  END FUNCTION scal_pro

     !> Computes a distance measure of to points in cartesian coordinates
     !! by the sqaure of the chord length (Sehne)
  ELEMENTAL FUNCTION chord_d(x0, x1) RESULT(dist)
! !REVISION HISTORY:  
! Developed by Hermann Asensio (2009).
    TYPE(cartesian_coordinates), INTENT(IN) :: x0, x1
    REAL(wp) :: dist

    dist =   (x0%x(1)-x1%x(1)) * (x0%x(1)-x1%x(1)) &
           + (x0%x(2)-x1%x(2)) * (x0%x(2)-x1%x(2)) &
           + (x0%x(3)-x1%x(3)) * (x0%x(3)-x1%x(3))
 
  END FUNCTION chord_d




!> Finds the intersection of two great circles.
 ELEMENTAL FUNCTION inter_section (p0, p1, v0, v1) RESULT(p)
!! DESCRIPTION:
!! Finds the intersection of two great circles.
!
!! !REVISION HISTORY:  
!! Developed  by Luis Kornblueh  (2004).
!! Bug fix: avoid divison by zero, Hermann Asensio

    

    TYPE(cartesian_coordinates), INTENT(IN) :: p0, p1
    TYPE(cartesian_coordinates), INTENT(IN) :: v0, v1
    TYPE(cartesian_coordinates) :: p

    TYPE(cartesian_coordinates) :: en, cn
    
    REAL(wp) :: z_pn

    ! normal to plane where edge lies

    en = vector_product (v0, v1)

    ! normal tp plane where triangle centers lies

    cn = vector_product (p0, p1)
 
    p = vector_product (en, cn)

    ! rescale on sphere surface

    z_pn = SQRT(DOT_PRODUCT(p%x, p%x))
    
    IF (z_pn /= 0.) then
        p%x = p%x/z_pn
    ENDIF

    IF (DOT_PRODUCT(p0%x, p%x) < 0.0_wp) THEN
      p%x = -p%x
    END IF
    
  END FUNCTION inter_section

END MODULE mo_additional_geometry


