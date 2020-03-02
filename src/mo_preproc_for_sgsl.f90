!-----------------------------------------------------------------------------------!
! Read in topography from NetCDF file and calculate maximum gradient of the         !
! topography. Required as input field for the runoff calculations of terra_ml.      !
! Input data is the "altitude" parameter from the GLOBE dataset. The resolution     !
! of the data set is 30''.                                                          !
!                                                                                   !

! based on code from Linda Schlemmer, IAC ETH Zurich, July 2016                                        !
! Daniel Luethi, IAC ETH Zurich, August 2016                                        !
!-----------------------------------------------------------------------------------!
MODULE mo_preproc_for_sgsl
  
  USE mo_logging

  USE mo_kind,                  ONLY: wp, i4

  USE mo_grid_structures,       ONLY: reg_lonlat_grid !< Definition of Data Type to describe a regular (lonlat) grid

  PUBLIC :: preprocess_globe_for_sgsl

  CONTAINS

  SUBROUTINE preprocess_globe_for_sgsl(h_block,&
       &                               sl_block, &
       &                               undef_topo, &
       &                               ta_grid)

    IMPLICIT NONE     

    TYPE(reg_lonlat_grid), INTENT(IN) :: ta_grid

    INTEGER(KINd=i4), INTENT(IN)      :: undef_topo

    INTEGER(KIND=i4), INTENT(IN)      :: h_block(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg)

    INTEGER(KIND=i4), INTENT(OUT)     :: sl_block(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg)

    ! local variables
    INTEGER(KIND=i4)                  :: nx, ny, i, j, x_inner, y_inner, errorcode, &
                                         ! local working array for h_block
         &                               zh_block(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg)


    REAL(KIND=wp)                     :: dx, dy, len, oolen, oolenx, ooleny, grad(9), zlats, crlat, &
         &                               dx0, dx2, zlats0, zlats2, crlat0, crlat2, len0, len2, &
         &                               oolen0, oolen2


    REAL(KIND=wp),ALLOCATABLE         :: h_block_inner(:,:), &
         &                               lat(:), &
         &                               lon(:)

    REAL(KIND=wp), PARAMETER          :: r_earth  =  6371.229E3, & ! radius of the earth
         &                               pi       =  4.0 * ATAN (1.0), &
         &                               degrad   =   pi / 180.0, &
         &                               dlat     =  30./3600., & ! resolution
         &                               dlon     =  30./3600., & ! resolution
         &                               eps      =  1.E-9, &
         &                               add_offset = 0., &
         &                               scale_factor = 0.001, &
         &                               r_scfct = 1. / scale_factor

    x_inner=ta_grid%nlon_reg-2
    y_inner=ta_grid%nlat_reg-2

    ALLOCATE (h_block_inner(x_inner, y_inner ), STAT=errorcode)
    IF(errorcode/=0) CALL logging%error('Cant allocate h_block_inner',__FILE__,__LINE__)

    ALLOCATE (lat(0:y_inner+1 ), STAT=errorcode)
    IF(errorcode/=0) CALL logging%error('Cant allocate lat',__FILE__,__LINE__)

    ALLOCATE (lon(0:x_inner+1 ), STAT=errorcode)
    IF(errorcode/=0) CALL logging%error('Cant allocate lat',__FILE__,__LINE__)

    zh_block = h_block

    h_block_inner = h_block(2:x_inner +1, 2:y_inner +1)

    ! determine lat/lon
    DO i = 0, y_inner +1 
      lat(i)= ta_grid%start_lat_reg + (i)*ta_grid%dlat_reg
    ENDDO

    DO i = 1, x_inner +1 
      lon(i)= ta_grid%start_lon_reg + (i)*ta_grid%dlon_reg
    ENDDO

    ! Calculations
    WHERE (ABS(zh_block-undef_topo) <= eps) zh_block = 0.0

    dy = r_earth * dlat * degrad
    ooleny = 1./dy
    DO   j = 2,y_inner 
       zlats0 = lat(j-1)
       zlats  = lat(j)
       zlats2 = lat(j+1)
       crlat0 = COS ( zlats0  * degrad )
       crlat  = COS ( zlats  * degrad )
       crlat2 = COS ( zlats2  * degrad )
       dx0    = dlon * r_earth * degrad * crlat0
       dx     = dlon * r_earth * degrad * crlat
       dx2    = dlon * r_earth * degrad * crlat2
       len0   = sqrt(dx0**2+dy**2)
       len2   = sqrt(dx2**2+dy**2)
       oolen0 = 1./len0
       oolen2 = 1./len2
       oolenx = 1./dx
       DO i = 2, x_inner
         IF (abs(h_block_inner(i,j)-undef_topo) <= eps) THEN
           grad(1)      = oolenx * (zh_block(i,j)-zh_block(i-1,j  ))
           grad(2)      = oolenx * (zh_block(i,j)-zh_block(i+1,j  ))
           grad(3)      = ooleny * (zh_block(i,j)-zh_block(i  ,j-1))
           grad(4)      = ooleny * (zh_block(i,j)-zh_block(i  ,j+1))
           grad(5)      = oolen0 * (zh_block(i,j)-zh_block(i-1,j-1))
           grad(6)      = oolen2 * (zh_block(i,j)-zh_block(i-1,j+1))
           grad(7)      = oolen0 * (zh_block(i,j)-zh_block(i+1,j-1))
           grad(8)      = oolen2 * (zh_block(i,j)-zh_block(i+1,j+1))
           grad(9)      = 0.0
           sl_block(i,j)   = NINT((MAXVAL(grad)-add_offset) * r_scfct)
         ELSE
          sl_block(i,j)   = undef_topo
         ENDIF
       END DO
    END DO

    DEALLOCATE (h_block_inner, STAT=errorcode)
    IF(errorcode/=0) CALL logging%error('Cant deallocate h_block_inner',__FILE__,__LINE__)

    DEALLOCATE (lat, STAT=errorcode)
    IF(errorcode/=0) CALL logging%error('Cant deallocate lat',__FILE__,__LINE__)

    DEALLOCATE (lon, STAT=errorcode)
    IF(errorcode/=0) CALL logging%error('Cant deallocate lon',__FILE__,__LINE__)
 
  END SUBROUTINE preprocess_globe_for_sgsl

END MODULE mo_preproc_for_sgsl
