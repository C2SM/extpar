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
  USE mo_topo_data,             ONLY: max_tiles
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

    CHARACTER (LEN=1024)              :: topo_files(1:max_tiles), &         !< filenames globe raw data
         &                               sgsl_files(1:max_tiles)
    
    INTEGER (KIND=i4)                 :: ntiles_column, &      !< number of tile columns
         &                               ntiles_row, &         !< number of tile rows
         &                               itopo_type, idx, &
         &                               ierr, nuin

    LOGICAL                           :: lsso_param, &
         &                               lsubtract_mean_slope

    !> namelist with information on orography data input
    NAMELIST /orography_raw_data/ itopo_type, lsso_param, lsubtract_mean_slope, &
         &                        raw_data_orography_path, ntiles_column, ntiles_row, &
         &                        sgsl_files,topo_files

    nuin = free_un()  ! function free_un returns free Fortran unit number
    OPEN(nuin,FILE=TRIM(preproc_namelist), IOSTAT=ierr)
    IF (ierr /= 0) THEN
      WRITE(message_text,*)'Cannot open ', TRIM(preproc_namelist)
      CALL logging%error(message_text,__FILE__, __LINE__) 
    ENDIF

    READ(nuin, NML=orography_raw_data, IOSTAT=ierr)
    IF (ierr /= 0) THEN
      CALL logging%error('Cannot read in namelist orography_raw_data',__FILE__, __LINE__) 
    ENDIF
  
    DO idx= 1, max_tiles
      IF (LEN(TRIM(topo_files(idx)))<=50) THEN 
        CALL topo_grad_globe(topo_files(idx), output_files(idx))
      ENDIF
    END DO

  END SUBROUTINE prepare_preproc

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
    ALLOCATE(s_oro(0:nx+1,0:ny+1),hsurf_inner(nx,ny))
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
    
    PRINT*, latid, lonid
    PRINT*, MAXVAL(lon), MINVAL (lon)
    PRINT*, MAXVAL(lat), MINVAL (lat)
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
   
    s_oro(0,:)  = mdv
    s_oro(nx+1,:)  = mdv
    s_oro(:,0)  = mdv
    s_oro(:,ny+1:)  = mdv
  !  print*,'CREATE NEW NetCDF FILE'

    !* enter define mode
    status = nf90_create (outfile,  NF90_NETCDF4, ncido)
    IF (STATUS .NE. NF90_NOERR) PRINT *, NF90_STRERROR(STATUS)


    !* define dimensions
    status = nf90_def_dim(ncido, 'lon', nxp2, londim)
    IF (STATUS .NE. NF90_NOERR) PRINT *, NF90_STRERROR(STATUS)
    status = nf90_def_dim(ncido, 'lat', nyp2, latdim)
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
    STATUS = NF90_PUT_VAR(ncido, lonid, lon)
    IF (STATUS .NE. NF90_NOERR) PRINT *, NF90_STRERROR(STATUS)

    STATUS = NF90_PUT_VAR(ncido, latid, lat)
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
