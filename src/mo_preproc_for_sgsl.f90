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
  
  USE netcdf
  USE mo_logging

  USE mo_kind,                  ONLY: wp, i4

  USE mo_grid_structures,       ONLY: reg_lonlat_grid !< Definition of Data Type to describe a regular (lonlat) grid
  USE mo_sgsl_data,             ONLY: max_tiles
  USE mo_utilities_extpar,      ONLY: free_un

  PUBLIC :: preprocess_globe_for_sgsl, &
       &    topo_grad_globe

  CONTAINS

  SUBROUTINE prepare_preproc (preproc_namelist, &
       &                      output_files)

    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN)     :: preproc_namelist
    
    CHARACTER (LEN=1024), INTENT(IN)  :: output_files(1:max_tiles)

    CHARACTER (LEN=1024)              :: raw_data_orography_path        !< path to raw data

    CHARACTER (LEN=1024)              :: topo_files(1:max_tiles)        !< filenames globe raw data
    
    INTEGER (KIND=i4)                 :: ntiles_column, &      !< number of tile columns
         &                               ntiles_row, &         !< number of tile rows
         &                               itopo_type, idx, &
         &                               ierr, nuin

    LOGICAL                           :: lsso_param, &
         &                               lsubtract_mean_slope

    !> namelist with information on orography data input
    NAMELIST /orography_raw_data/ itopo_type, lsso_param, lsubtract_mean_slope, &
         &                        raw_data_orography_path, ntiles_column, ntiles_row, topo_files

    nuin = free_un()  ! function free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(preproc_namelist), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(preproc_namelist)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=orography_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist orography_io_extpar',__FILE__, __LINE__) 
    ENDIF
  
    DO idx= 1, max_tiles
      IF (LEN(TRIM(topo_files(idx)))<=50) THEN 
        CALL topo_grad_globe(topo_files(idx), output_files(idx))
      ENDIF
    END DO

  END SUBROUTINE prepare_preproc

  SUBROUTINE preprocess_globe_for_sgsl(h_block,&
       &                               sl_block, &
       &                               undef_topo, &
       &                               ta_grid)

    IMPLICIT NONE     

    TYPE(reg_lonlat_grid), INTENT(IN) :: ta_grid

    REAL(KIND=wp), INTENT(IN)         :: undef_topo

    INTEGER(KIND=i4), INTENT(IN)      :: h_block(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg)

    REAL(KIND=wp), INTENT(OUT)        :: sl_block(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg)

    ! local variables
    INTEGER(KIND=i4)                  ::  i, j, x_inner, y_inner, errorcode, &
                                         ! local working array for h_block
         &                               zh_block(1:ta_grid%nlon_reg,1:ta_grid%nlat_reg)


    REAL(KIND=wp)                     :: dx, dy, oolenx, ooleny, grad(9), zlats, crlat, &
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
         &                               eps      =  1.E-9

    x_inner=ta_grid%nlon_reg-2
    y_inner=ta_grid%nlat_reg-2

    ALLOCATE (h_block_inner(x_inner, y_inner ), STAT=errorcode)
    IF(errorcode/=0) CALL logging%error('Cant allocate h_block_inner',__FILE__,__LINE__)

    ALLOCATE (lat(0:y_inner+1 ), STAT=errorcode)
    IF(errorcode/=0) CALL logging%error('Cant allocate lat',__FILE__,__LINE__)

    ALLOCATE (lon(0:x_inner+1 ), STAT=errorcode)
    IF(errorcode/=0) CALL logging%error('Cant allocate lat',__FILE__,__LINE__)

    zh_block = h_block

 !   WRITE(logging%fileunit,*)MAXVAL(zh_block), MINVAL(zh_block)

    h_block_inner = h_block(2:x_inner +1, 2:y_inner +1)

    ! determine lat/lon
    DO i = 0, y_inner +1 
      lat(i)= ta_grid%start_lat_reg + (i)*ta_grid%dlat_reg
    ENDDO

    DO i = 0, x_inner +1 
      lon(i)= ta_grid%start_lon_reg + (i)*ta_grid%dlon_reg
    ENDDO

    ! Calculations
    WHERE (ABS(zh_block-undef_topo) <= eps) zh_block = 0

 !   WRITE(logging%fileunit,*)MAXVAL(zh_block), MINVAL(zh_block)

    dy = r_earth * dlat * degrad
 !   WRITE(logging%fileunit,*) dy, lat(0), lon(0)
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
         IF (abs(h_block_inner(i,j)-undef_topo) > eps) THEN
           grad(1)      = oolenx * (zh_block(i,j)-zh_block(i-1,j  ))
           grad(2)      = oolenx * (zh_block(i,j)-zh_block(i+1,j  ))
           grad(3)      = ooleny * (zh_block(i,j)-zh_block(i  ,j-1))
           grad(4)      = ooleny * (zh_block(i,j)-zh_block(i  ,j+1))
           grad(5)      = oolen0 * (zh_block(i,j)-zh_block(i-1,j-1))
           grad(6)      = oolen2 * (zh_block(i,j)-zh_block(i-1,j+1))
           grad(7)      = oolen0 * (zh_block(i,j)-zh_block(i+1,j-1))
           grad(8)      = oolen2 * (zh_block(i,j)-zh_block(i+1,j+1))
           grad(9)      = 0.0
           sl_block(i,j)   = MAXVAL(grad)
          ! IF (MOD(j, 100) == 0)THEN
          !    WRITE(logging%fileunit,*)'start output'
          !    WRITE(logging%fileunit,*) crlat0, crlat,crlat2, dx0, dx,len0,len2,oolen0,'linebreak', &
          !        & oolen2, oolenx, MAXVAL(grad)
          !    WRITE(logging%fileunit,*)'end output'
          !    WRITE(logging%fileunit,*)sl_block(i,j)
          !  ENDIF

         ELSE
          sl_block(i,j)   = undef_topo
         ENDIF
       END DO
    END DO

    WRITE(logging%fileunit,*)'MAXVAL(sl_block: ',MAXVAL(sl_block),'MINVAL(sl_block): ', MINVAL(sl_block)

    DEALLOCATE (h_block_inner, STAT=errorcode)
    IF(errorcode/=0) CALL logging%error('Cant deallocate h_block_inner',__FILE__,__LINE__)

    DEALLOCATE (lat, STAT=errorcode)
    IF(errorcode/=0) CALL logging%error('Cant deallocate lat',__FILE__,__LINE__)

    DEALLOCATE (lon, STAT=errorcode)
    IF(errorcode/=0) CALL logging%error('Cant deallocate lon',__FILE__,__LINE__)
 
  END SUBROUTINE preprocess_globe_for_sgsl

  SUBROUTINE topo_grad_globe( infile, outfile)


    IMPLICIT NONE     

    INTEGER nx, ny, nt, nxp2, nyp2

    CHARACTER(LEN=80), INTENT(IN)  :: infile,outfile
    INTEGER, PARAMETER :: short = SELECTED_INT_KIND(4)
    INTEGER, PARAMETER :: double = SELECTED_REAL_KIND(15)

    ! input fields
    REAL, ALLOCATABLE :: &

         hsurf  (:,:), hsurf_inner(:,:)

    ! output fields
    INTEGER(short), ALLOCATABLE :: &

         s_oro  (:,:)

    ! grid
    REAL(double), ALLOCATABLE :: &

         lat   (:),    & !
         lon   (:)

    !* netCDF id
    INTEGER  ncid, ncido, status
    !* dimension ids
    INTEGER londim, latdim
    !* variable ids
    INTEGER lonid, latid
    INTEGER varid, outid, mapid

    INTEGER i, j, nargs

    CHARACTER(LEN=100) name, char, comment

    REAL(double) :: dx, dy, len, oolen, oolenx, ooleny, grad(9), zlats, crlat 
    REAL(double) :: dx0, dx2, zlats0, zlats2, crlat0, crlat2, len0, len2
    REAL(double) :: oolen0, oolen2
    INTEGER(short) :: mdv

    REAL(double), PARAMETER :: r_earth  =  6371.229E3 ! radius of the earth
    REAL(double), PARAMETER :: pi       =  4.0 * ATAN (1.0)
    REAL(double), PARAMETER :: degrad   =   pi / 180.0
    REAL(double), PARAMETER :: dlat     =  30./3600. ! resolution
    REAL(double), PARAMETER :: dlon     =  30./3600. ! resolution
    REAL(double), PARAMETER :: eps      =  1.E-9
    REAL(double), PARAMETER :: add_offset = 0.
    REAL(double), PARAMETER :: scale_factor = 0.001
    REAL(double), PARAMETER :: r_scfct = 1. / scale_factor

    ! Read user input
    
    PRINT*,'in: ', TRIM(infile), 'out: ',TRIM(outfile)

    ! Open file
    
    status = nf90_open(infile, nf90_nowrite, ncid)
    IF (status .NE. NF90_NOERR) THEN
       PRINT *, NF90_STRERROR(status)
       STOP
    ENDIF
    PRINT *,'file infile opened'
    ! inquire dimensions

    status=nf90_inq_dimid(ncid,"lon",lonid)
    status=nf90_inquire_dimension(ncid,lonid,len=nxp2)
    status=nf90_inq_dimid(ncid,"lat",latid)
    status=nf90_inquire_dimension(ncid,latid,len=nyp2)
    PRINT *,'dimensions read:', nxp2, nyp2

    ! allocate fields
    nx = nxp2 - 2
    ny = nyp2 - 2

    ALLOCATE(hsurf(0:nx+1,0:ny+1))
    ALLOCATE(s_oro(nx,ny),hsurf_inner(nx,ny))
    PRINT *,'2d fields allocated'

    ALLOCATE(lon(0:nx+1))
    ALLOCATE(lat(0:ny+1))
    PRINT *,'1d fields allocated'

  !  s_oro(1:nx,1:ny)    = 0.0
  !  print *,'s_oro initialized'
  !  hsurf(0:nx+1,0:ny+1)    = 0.0
  !  print *,'hsurf initialized'

    lat(:)       = 0.0
    lon(:)       = 0.0
    PRINT *,'arrays initialized'

    ! Read in variables
    
    !print*,'READING VARIABLES'

    status = nf90_inq_varid(ncid,"lon", lonid)
    IF (status .NE. NF90_NOERR) PRINT *, NF90_STRERROR(status)
    status = nf90_get_var(ncid,lonid,lon)
    IF (status .NE. NF90_NOERR) PRINT *, NF90_STRERROR(status)
    PRINT *,'lon read'
   
    status = nf90_inq_varid(ncid,"lat", latid)
    IF (status .NE. NF90_NOERR) PRINT *, NF90_STRERROR(status)
    status = nf90_get_var(ncid,latid,lat)
    IF (status .NE. NF90_NOERR) PRINT *, NF90_STRERROR(status)
    PRINT *,'lat read'
    
    status = nf90_inq_varid(ncid,"altitude", varid)
    IF (status .NE. NF90_NOERR) PRINT *, NF90_STRERROR(status)
    status = nf90_get_var(ncid,varid,hsurf)
    PRINT *,'hsurf read'
    IF (status .NE. NF90_NOERR) PRINT *, NF90_STRERROR(status)
    status = nf90_get_att(ncid, varid,'_FillValue',mdv)
    CALL check_err(status)
    status = nf90_get_att(ncid, varid,'comment',comment)
    CALL check_err(status)

    PRINT *,'mdv = ',mdv
    
    ! Calculations
    hsurf_inner(:,:) = hsurf(1:nx,1:ny)
    WHERE (ABS(hsurf-mdv) <= eps) hsurf=0.0

    dy = r_earth * dlat * degrad
    ooleny = 1./dy
    DO   j = 1, ny
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
       DO i = 1, nx
         IF (abs(hsurf_inner(i,j)-mdv).gt.eps) THEN
           grad(1)      = oolenx * (hsurf(i,j)-hsurf(i-1,j  ))
           grad(2)      = oolenx * (hsurf(i,j)-hsurf(i+1,j  ))
           grad(3)      = ooleny * (hsurf(i,j)-hsurf(i  ,j-1))
           grad(4)      = ooleny * (hsurf(i,j)-hsurf(i  ,j+1))
           grad(5)      = oolen0 * (hsurf(i,j)-hsurf(i-1,j-1))
           grad(6)      = oolen2 * (hsurf(i,j)-hsurf(i-1,j+1))
           grad(7)      = oolen0 * (hsurf(i,j)-hsurf(i+1,j-1))
           grad(8)      = oolen2 * (hsurf(i,j)-hsurf(i+1,j+1))
           grad(9)      = 0.0
           s_oro(i,j)   = NINT((MAXVAL(grad)-add_offset) * r_scfct)
         ELSE
           s_oro(i,j)   = mdv
         ENDIF
       END DO

    END DO
   
  !  print*,'CREATE NEW NetCDF FILE'

    !* enter define mode
    status = nf90_create (outfile,  NF90_NETCDF4, ncido)
    IF (STATUS .NE. NF90_NOERR) PRINT *, NF90_STRERROR(STATUS)


    !* define dimensions
    status = nf90_def_dim(ncido, 'lon', nx, londim)
    IF (STATUS .NE. NF90_NOERR) PRINT *, NF90_STRERROR(STATUS)
    status = nf90_def_dim(ncido, 'lat', ny, latdim)
    IF (STATUS .NE. NF90_NOERR) PRINT *, NF90_STRERROR(STATUS)

    !* define variables

    status = nf90_def_var(ncido, 'lon', NF90_DOUBLE,(/londim/), lonid)
    CALL check_err(status)

    status = nf90_def_var(ncido, 'lat', NF90_DOUBLE,(/latdim/), latid)
    CALL check_err(status)

    status = nf90_def_var(ncido, 'S_ORO', NF90_SHORT,(/londim,latdim/),outid)
    CALL check_err(status)

    status = nf90_def_var(ncido, 'regular_grid', NF90_CHAR,mapid)
    CALL check_err(status)

  !  print*,'ATTRIBUTES'
    status = nf90_inq_varid(ncid,"lat", varid)
    status = nf90_get_att(ncid, varid,'long_name',name)
    CALL check_err(status)
    status = nf90_put_att(ncido, latid,'long_name',TRIM(name))
    CALL check_err(status)
    status = nf90_get_att(ncid, varid,'units',name)
    CALL check_err(status)
    status = nf90_put_att(ncido, latid,'units',TRIM(name))
    CALL check_err(status)

    status = nf90_inq_varid(ncid,"lon", varid)
    status = nf90_get_att(ncid, varid,'long_name',name)
    CALL check_err(status)
    status = nf90_put_att(ncido, lonid,'long_name',TRIM(name))
    CALL check_err(status)
    status = nf90_get_att(ncid, varid,'units',name)
    CALL check_err(status)
    status = nf90_put_att(ncido, lonid,'units',TRIM(name))
    CALL check_err(status)

    status = nf90_put_att(ncido, outid,'standard_name','topography gradient')
    CALL check_err(status)
    status = nf90_put_att(ncido, outid,'long_name','maximum local gradient of surface height')
    CALL check_err(status)
    status = nf90_put_att(ncido, outid,'scale_factor',scale_factor)
    CALL check_err(status)
    status = nf90_put_att(ncido, outid,'add_offset',add_offset)
    CALL check_err(status)
    status = nf90_put_att(ncido, outid,'units','')
    CALL check_err(status)

    status = nf90_put_att(ncido, outid,'grid_mapping','regular_grid')
    CALL check_err(status)
    status = nf90_put_att(ncido, outid,'_FillValue',mdv)
    CALL check_err(status)
    status = nf90_put_att(ncido, outid,'comment',trim(comment))
    CALL check_err(status)

    status = nf90_put_att(ncido, mapid,'grid_mapping_name','latitude_longitude')
    CALL check_err(status)

  !* leave define mode
    status = NF90_ENDDEF(ncido)
    IF (status .NE. NF90_NOERR) PRINT *, NF90_STRERROR(status)

  !* store variables
    STATUS = NF90_PUT_VAR(ncido, lonid, lon(1:nx))
    IF (STATUS .NE. NF90_NOERR) PRINT *, NF90_STRERROR(STATUS)

    STATUS = NF90_PUT_VAR(ncido, latid, lat(1:ny))
    IF (STATUS .NE. NF90_NOERR) PRINT *, NF90_STRERROR(STATUS)
    
    STATUS = NF90_PUT_VAR(ncido, outid,s_oro)
    IF (STATUS .NE. NF90_NOERR) PRINT *, NF90_STRERROR(STATUS)

    DEALLOCATE(hsurf)
    DEALLOCATE(hsurf_inner)
    DEALLOCATE(s_oro)

    DEALLOCATE(lon)
    DEALLOCATE(lat)

  END SUBROUTINE topo_grad_globe

  SUBROUTINE check_err(iret)

    IMPLICIT NONE
    INTEGER iret
    IF (iret .NE. NF90_NOERR) THEN
       PRINT *, nf90_strerror(iret)
       STOP
    ENDIF
  
  END SUBROUTINE check_err

END MODULE mo_preproc_for_sgsl
