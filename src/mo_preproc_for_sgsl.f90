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

  PUBLIC :: preproc_orography, &
       &    topo_grad_globe


  CONTAINS

  ! wrapper function for the preprocessing of raw orography data
  ! and calls the right subroutine for itopo_type (GLOBE or ASTER)
  SUBROUTINE preproc_orography (raw_data_orography_path, &
       &                        topo_files, &
       &                        sgsl_files, &
       &                        itopo_type, &
       &                        ntiles_row, &
       &                        ntiles_column)

    IMPLICIT NONE

    CHARACTER (LEN=1024),INTENT(IN)   :: raw_data_orography_path, &
         &                               topo_files(1:max_tiles), & 
         &                               sgsl_files(1:max_tiles)
    
    INTEGER (KIND=i4),INTENT(IN)      :: ntiles_column, &      !< number of tile columns
         &                               ntiles_row, &         !< number of tile rows
         &                               itopo_type

    !local variables
    INTEGER(KIND=i4)                  :: ntiles_tot, idx
    
    CALL logging%info('SGSL: Enter routine: preproc_orography')

    IF (itopo_type == 2) THEN
      CALL logging%warning('SGSL: Only tested and validated for GLOBE data')
    ENDIF

    ntiles_tot = ntiles_row * ntiles_column

    CALL logging%info('Start processing topo tiles...')


    DO idx= 1, ntiles_tot 
      CALL topo_subgrid_slope( TRIM(raw_data_orography_path) // TRIM(topo_files(idx)), &
           &                   sgsl_files(idx), &
           &                   itopo_type)
    END DO

    CALL logging%info('                       ...done')

    CALL logging%info('SGSL: Exit routine: preproc_orography')

  END SUBROUTINE preproc_orography

  ! preprocessing of raw topography (GLOBE or ASTER) to determine subgrid-slope (SGSL)
  SUBROUTINE topo_subgrid_slope( infile, outfile, itopo_type)

    IMPLICIT NONE


    CHARACTER(LEN=*), INTENT(IN)  :: infile, outfile

    REAL(KIND=wp), ALLOCATABLE    :: hsurf(:,:), &
         &                           hsurf_inner(:,:), &
         &                           lat(:), &
         &                           lon(:)

    INTEGER(KIND=i4),INTENT(IN)   :: itopo_type

    ! output fields
    INTEGER(KIND=i4),ALLOCATABLE  :: s_oro  (:,:)

    !* netCDF id
    INTEGER(KIND=i4)              :: ncid, ncido, status, &
         &                           londim, latdim, &
         &                           lonid, latid, &
         &                           varid, outid, mapid, &
         &                           i, j, &
         &                           mdv, &
         &                           nx, ny, nxp2, nyp2

    CHARACTER(LEN=100)            :: name, comment

    REAL(KIND=wp)                 :: dx, dy, oolenx, ooleny, grad(9), zlats, crlat, &
         &                           dx0, dx2, zlats0, zlats2, crlat0, crlat2, len0, len2, &
         &                           oolen0, oolen2

    REAL(KIND=wp), PARAMETER      :: r_earth        =  6371.229E3, & ! radius of the earth
         &                           pi             =  4.0 * ATAN (1.0), &
         &                           degrad         =   pi / 180.0, &
         &                           dlat_globe     =  30./3600. , &! resolution
         &                           dlon_globe     =  30./3600. , &! resolution
         &                           dlat_aster     =  30./3600. , &! resolution
         &                           dlon_aster     =  30./3600. , &! resolution
         &                           eps            =  1.E-9, &
         &                           add_offset     = 0., &
         &                           scale_factor   = 0.001, &
         &                           r_scfct        = 1. / scale_factor


    WRITE(message_text,*) TRIM(infile), ' --> ', TRIM(outfile)
    CALL logging%info(message_text)

    ! Open file
    status = nf90_open(infile, nf90_nowrite, ncid)
    CALL check_err(status, __FILE__, __LINE__)
    
    ! inquire dimensions
    status=nf90_inq_dimid(ncid,"lon",lonid)
    CALL check_err(status, __FILE__, __LINE__)
    status=nf90_inquire_dimension(ncid,lonid,len=nxp2)
    CALL check_err(status, __FILE__, __LINE__)
    status=nf90_inq_dimid(ncid,"lat",latid)
    CALL check_err(status, __FILE__, __LINE__)
    status=nf90_inquire_dimension(ncid,latid,len=nyp2)
    CALL check_err(status, __FILE__, __LINE__)

    ! allocate fields
    nx = nxp2 - 2
    ny = nyp2 - 2

    ALLOCATE(hsurf(0:nx+1,0:ny+1))
    ALLOCATE(s_oro(0:nx+1,0:ny+1),hsurf_inner(nx,ny))

    ALLOCATE(lon(0:nx+1))
    ALLOCATE(lat(0:ny+1))

    lat(:)       = 0.0
    lon(:)       = 0.0

    ! Read in variables
    status = nf90_inq_varid(ncid,"lon", lonid)
    CALL check_err(status, __FILE__, __LINE__)
    status = nf90_get_var(ncid,lonid,lon)
    CALL check_err(status, __FILE__, __LINE__)
   
    status = nf90_inq_varid(ncid,"lat", latid)
    CALL check_err(status, __FILE__, __LINE__)
    status = nf90_get_var(ncid,latid,lat)
    CALL check_err(status, __FILE__, __LINE__)
    
    IF (itopo_type == 1)THEN
      status = nf90_inq_varid(ncid,"altitude", varid)
    ELSE
      status = nf90_inq_varid(ncid,"Z", varid)
    ENDIF
    CALL check_err(status, __FILE__, __LINE__)

    status = nf90_get_var(ncid,varid,hsurf)
    CALL check_err(status, __FILE__, __LINE__)

    status = nf90_get_att(ncid, varid,'_FillValue',mdv)
    CALL check_err(status, __FILE__, __LINE__)

    IF (itopo_type == 1) THEN
      status = nf90_get_att(ncid, varid,'comment',comment)
      CALL check_err(status, __FILE__, __LINE__)
    ENDIF

    ! Calculations
    hsurf_inner(:,:) = hsurf(1:nx,1:ny)
    WHERE (ABS(hsurf-mdv) <= eps) hsurf=0.0

    IF (itopo_type == 1) THEN !GLOBE
      dy = r_earth * dlat_globe * degrad
      ooleny = 1./dy
      DO   j = 1, ny
        zlats0 = lat(j-1)
        zlats  = lat(j)
        zlats2 = lat(j+1)
        crlat0 = COS ( zlats0  * degrad )
        crlat  = COS ( zlats  * degrad )
        crlat2 = COS ( zlats2  * degrad )
        dx0    = dlon_globe * r_earth * degrad * crlat0
        dx     = dlon_globe * r_earth * degrad * crlat
        dx2    = dlon_globe * r_earth * degrad * crlat2
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

    ELSE !ASTER
      dy = r_earth * dlat_aster * degrad
      ooleny = 1./dy
      DO   j = 1, ny
        zlats  = lat(j)
        crlat  = COS ( zlats  * degrad )
        dx     = dlon_aster * r_earth * degrad * crlat
        len0   = sqrt(dx**2+dy**2)
        oolen0  = 1./len0
        oolenx = 1./dx
        DO i = 1, nx
          IF (abs(hsurf_inner(i,j)-mdv).gt.eps) THEN
            grad(1)      = oolenx * (hsurf(i,j)-hsurf(i-1,j  ))
            grad(2)      = oolenx * (hsurf(i,j)-hsurf(i+1,j  ))
            grad(3)      = ooleny * (hsurf(i,j)-hsurf(i  ,j-1))
            grad(4)      = ooleny * (hsurf(i,j)-hsurf(i  ,j+1))
            grad(5)      = oolen0  * (hsurf(i,j)-hsurf(i-1,j-1))
            grad(6)      = oolen0  * (hsurf(i,j)-hsurf(i-1,j+1))
            grad(7)      = oolen0  * (hsurf(i,j)-hsurf(i+1,j-1))
            grad(8)      = oolen0  * (hsurf(i,j)-hsurf(i+1,j+1))
            grad(9)      = 0.0
            s_oro(i,j)   = NINT((MAXVAL(grad)-add_offset) * r_scfct)
          ELSE
            s_oro(i,j)   = mdv
          ENDIF
        END DO
      END DO
      
    ENDIF !itopo_type
   
    ! fill edges with undef values
    s_oro(0,:)     = mdv
    s_oro(nx+1,:)  = mdv
    s_oro(:,0)     = mdv
    s_oro(:,ny+1:) = mdv

    ! enter define mode
    status = nf90_create (outfile,  NF90_NETCDF4, ncido)
    CALL check_err(status, __FILE__, __LINE__)

    ! define dimensions
    status = nf90_def_dim(ncido, 'lon', nxp2, londim)
    CALL check_err(status, __FILE__, __LINE__)

    status = nf90_def_dim(ncido, 'lat', nyp2, latdim)
    CALL check_err(status, __FILE__, __LINE__)

    ! define variables
    status = nf90_def_var(ncido, 'lon', NF90_DOUBLE,(/londim/), lonid)
    CALL check_err(status, __FILE__, __LINE__)

    status = nf90_def_var(ncido, 'lat', NF90_DOUBLE,(/latdim/), latid)
    CALL check_err(status, __FILE__, __LINE__)

    status = nf90_def_var(ncido, 'S_ORO', NF90_INT,(/londim,latdim/),outid)
    CALL check_err(status, __FILE__, __LINE__)

    status = nf90_def_var(ncido, 'regular_grid', NF90_CHAR,mapid)
    CALL check_err(status, __FILE__, __LINE__)

    status = nf90_inq_varid(ncid,"lat", varid)
    CALL check_err(status, __FILE__, __LINE__)
    status = nf90_get_att(ncid, varid,'long_name',name)
    CALL check_err(status, __FILE__, __LINE__)
    status = nf90_put_att(ncido, latid,'long_name',TRIM(name))
    CALL check_err(status, __FILE__, __LINE__)
    status = nf90_get_att(ncid, varid,'units',name)
    CALL check_err(status, __FILE__, __LINE__)
    status = nf90_put_att(ncido, latid,'units',TRIM(name))
    CALL check_err(status, __FILE__, __LINE__)

    status = nf90_inq_varid(ncid,"lon", varid)
    CALL check_err(status, __FILE__, __LINE__)
    status = nf90_get_att(ncid, varid,'long_name',name)
    CALL check_err(status, __FILE__, __LINE__)
    status = nf90_put_att(ncido, lonid,'long_name',TRIM(name))
    CALL check_err(status, __FILE__, __LINE__)
    status = nf90_get_att(ncid, varid,'units',name)
    CALL check_err(status, __FILE__, __LINE__)
    status = nf90_put_att(ncido, lonid,'units',TRIM(name))
    CALL check_err(status, __FILE__, __LINE__)

    status = nf90_put_att(ncido, outid,'standard_name','topography gradient')
    CALL check_err(status, __FILE__, __LINE__)
    status = nf90_put_att(ncido, outid,'long_name','maximum local gradient of surface height')
    CALL check_err(status, __FILE__, __LINE__)
    status = nf90_put_att(ncido, outid,'scale_factor',scale_factor)
    CALL check_err(status, __FILE__, __LINE__)
    status = nf90_put_att(ncido, outid,'add_offset',add_offset)
    CALL check_err(status, __FILE__, __LINE__)
    status = nf90_put_att(ncido, outid,'units','')
    CALL check_err(status, __FILE__, __LINE__)

    status = nf90_put_att(ncido, outid,'grid_mapping','regular_grid')
    CALL check_err(status, __FILE__, __LINE__)
    status = nf90_put_att(ncido, outid,'_FillValue',mdv)
    CALL check_err(status, __FILE__, __LINE__)

    IF (itopo_type == 1) THEN
      status = nf90_put_att(ncido, outid,'comment',trim(comment))
    ELSE
      status = nf90_put_att(ncido, outid,'comment','ASTER tile: '//infile)
    ENDIF
    CALL check_err(status, __FILE__, __LINE__)

    status = nf90_put_att(ncido, mapid,'grid_mapping_name','latitude_longitude')
    CALL check_err(status, __FILE__, __LINE__)

    ! leave define mode
    status = NF90_ENDDEF(ncido)
    CALL check_err(status, __FILE__, __LINE__)

    ! store variables
    STATUS = NF90_PUT_VAR(ncido, lonid, lon)
    CALL check_err(status, __FILE__, __LINE__)

    STATUS = NF90_PUT_VAR(ncido, latid, lat)
    CALL check_err(status, __FILE__, __LINE__)
    
    STATUS = NF90_PUT_VAR(ncido, outid,s_oro)
    CALL check_err(status, __FILE__, __LINE__)

    DEALLOCATE(hsurf)
    DEALLOCATE(hsurf_inner)
    DEALLOCATE(s_oro)

    DEALLOCATE(lon)
    DEALLOCATE(lat)

  END SUBROUTINE topo_subgrid_slope

  SUBROUTINE check_err(iret,FILE,LINE)

    IMPLICIT NONE
  
    INTEGER(KIND=i4) iret, line
    CHARACTER(LEN=*), INTENT(IN)  :: file
    
    IF (iret .NE. NF90_NOERR) THEN
      CALL logging%error(nf90_strerror(iret), file, line)
    ENDIF
  
  END SUBROUTINE check_err

END MODULE mo_preproc_for_sgsl
