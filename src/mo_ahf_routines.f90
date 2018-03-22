!+ Fortran module with AHF data handling routines
!
!
! Description:
! Fortran module with AHF data handling routines
!
! Current Code Owner: DWD, <Name of person responsible for this code>
!    <smail, phone, fax and email>
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V2_4         2015-05-21 Hendrik Wouters
!  Initial release
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================
!> Fortran module with AHF data handling routines
!> \author Hermann Asensio
!>
MODULE mo_ahf_routines

USE mo_kind, ONLY: wp, &
                   sp, &
                   i8, &
                   i4

USE netcdf,      ONLY :   &
  nf90_open,              &
  nf90_close,             &
  nf90_inquire,           &
  nf90_inquire_dimension, &
  nf90_inquire_variable,  &
  nf90_inq_attname,       &
  nf90_inquire_attribute, &
  nf90_get_att,           &
  nf90_inquire_dimension, &
  nf90_inq_varid,          &
  nf90_get_var,            &
  nf90_noerr,              &
  nf90_strerror

USE netcdf,      ONLY:     &
  nf90_create,             &
  nf90_def_dim,            &
  nf90_def_var,            &
  nf90_enddef,             &
  nf90_redef,              &
  nf90_put_att,            &
  nf90_put_var

 
USE netcdf,      ONLY :   &
  NF90_CHAR,               &
  NF90_DOUBLE,             &
  NF90_FLOAT,              &
  NF90_INT,                &
  NF90_BYTE,               &
  NF90_SHORT


USE netcdf,      ONLY :   &
  NF90_GLOBAL,             &
  NF90_UNLIMITED,          &
  NF90_CLOBBER,            &
  NF90_NOWRITE


!> abort_extpar defined in MODULE utilities_extpar
USE mo_utilities_extpar, ONLY: abort_extpar


USE mo_io_utilities,     ONLY: check_netcdf
USE mo_io_units,         ONLY: filename_max


USE mo_grid_structures,  ONLY: reg_lonlat_grid



IMPLICIT NONE

PRIVATE
!
PUBLIC :: open_netcdf_AHF_data, &
          close_netcdf_AHF_data, &
          read_namelists_extpar_ahf, &
          read_ahf_data_input_namelist, &
          get_dimension_AHF_data, &
          get_AHF_data_coordinates, &
          get_one_row_AHF_data, &
          get_block_AHF_data, &
          get_pixel_AHF_data


CONTAINS

!---------------------------------------------------------------------------
!> subroutine to read namelist for AHF data settings for EXTPAR 
SUBROUTINE read_namelists_extpar_ahf(namelist_file, &
                                      iahf_type,    & !_br 14.04.16
                                      raw_data_ahf_path, &
                                      raw_data_ahf_filename, &
                                      ahf_buffer_file, &
                                      ahf_output_file)

  USE mo_utilities_extpar, ONLY: free_un ! function to get free unit number

  
  CHARACTER (len=filename_max), INTENT(IN) :: namelist_file !< filename with namelists for for EXTPAR settings

  ! AHF
  CHARACTER (len=filename_max) :: raw_data_ahf_path        !< path to raw data
  CHARACTER (len=filename_max) :: raw_data_ahf_filename !< filename AHF raw data

  CHARACTER (len=filename_max) :: ahf_buffer_file !< name for AHF buffer file
  CHARACTER (len=filename_max) :: ahf_output_file !< name for AHF output file

  INTEGER (KIND=i4)            :: iahf_type  !< ID of dataset used !_br 14.04.16

  !> namelist with filenames for AHF data input
  NAMELIST /ahf_raw_data/ raw_data_ahf_path, raw_data_ahf_filename, iahf_type !_br 14.04.16
  !> namelist with filenames for AHF data output
  NAMELIST /ahf_io_extpar/ ahf_buffer_file, ahf_output_file
  
  INTEGER           :: nuin !< unit number
  INTEGER (KIND=i4) :: ierr !< error flag

  nuin = free_un()  ! functioin free_un returns free Fortran unit number
  OPEN(nuin,FILE=TRIM(namelist_file), IOSTAT=ierr)
  
  READ(nuin, NML=ahf_raw_data, IOSTAT=ierr)
  READ(nuin, NML=ahf_io_extpar, IOSTAT=ierr)
  
  CLOSE(nuin)
  
END SUBROUTINE read_namelists_extpar_ahf
!---------------------------------------------------------------------------

    !> open netcdf-file and get netcdf unit file number
       SUBROUTINE open_netcdf_AHF_data(path_ahf_file, &
                                        ncid)
         CHARACTER (len=*), INTENT(in) :: path_ahf_file         !< filename with path for AHF raw data
         INTEGER, INTENT(out) :: ncid                             !< netcdf unit file number

        !! open netcdf file 
            call check_netcdf( nf90_open(TRIM(path_ahf_file),NF90_NOWRITE, ncid))

       END SUBROUTINE open_netcdf_AHF_data

        !> close netcdf-file 
       SUBROUTINE close_netcdf_AHF_data(ncid)
         INTEGER, INTENT(in) :: ncid                             !< netcdf unit file number

        !! close netcdf file 
            call check_netcdf( nf90_close( ncid))

       END SUBROUTINE close_netcdf_AHF_data


        !> read namelist with settings for COSMO target grid
        !> \author Hermann Asensio
       SUBROUTINE read_ahf_data_input_namelist(input_namelist_file,        &
                                                 raw_data_path,             &
                                                 raw_data_ahf_filename,    &
                                                 outputgrid_ahf_filename)

         USE mo_utilities_extpar, ONLY: free_un ! function to get free unit number

         USE mo_io_units,          ONLY: filename_max


           CHARACTER (len=*), INTENT(in) :: input_namelist_file !< file with input namelist 
           CHARACTER (len=filename_max), INTENT(out) :: raw_data_path        !< path to raw data
           CHARACTER (len=filename_max), INTENT(out) :: raw_data_ahf_filename !< filename ahf raw data
           CHARACTER (len=filename_max), INTENT(out) :: outputgrid_ahf_filename !< output filename

           !>Define the namelist group
           NAMELIST /ahf_data_input/ raw_data_path, raw_data_ahf_filename, outputgrid_ahf_filename

           
           INTEGER (KIND=i4) :: ierr !< error flag
           INTEGER                  :: nuin !< unit number

              nuin = free_un()  ! functioin free_un returns free Fortran unit number
              open(nuin,FILE=input_namelist_file, IOSTAT=ierr)
              !print *, ierr
              read(nuin, NML=ahf_data_input, IOSTAT=ierr)
              !print *, ierr

              close(nuin)

       END SUBROUTINE read_ahf_data_input_namelist


            !> inquire dimension information for AHF raw data 
       SUBROUTINE get_dimension_AHF_data(ncid, &
                                          nlon_ahf, &
                                          nlat_ahf)

        INTEGER, INTENT(in) :: ncid                             !< netcdf unit file number

        INTEGER (KIND=i4), INTENT(out) :: nlon_ahf !< number of grid elements in zonal direction for AHF data
        INTEGER (KIND=i4), INTENT(out) :: nlat_ahf !< number of grid elements in meridional direction for AHF data

        !local variables
        INTEGER :: ndimension                       !< number of dimensions in netcdf file
        INTEGER :: nVars                            !< number of variables in netcdf file
        INTEGER :: nGlobalAtts                      !< number of gloabal Attributes in netcdf file
        INTEGER :: unlimdimid                       !< id of unlimited dimension (e.g. time) in netcdf file

        INTEGER :: dimid                            !< id of dimension
        CHARACTER (len=80) :: dimname               !< name of dimensiona
        INTEGER :: length                           !< length of dimension



                !! look for numbers of dimensions, Variable, Attributes, and the dimid for the one possible unlimited dimension ( &
!& probably time)
                !; nf90_inquire input: ncid; nf90_inquire output: ndimension, nVars, nGlobalAtts,unlimdimid
                call check_netcdf (nf90_inquire(ncid,ndimension, nVars, nGlobalAtts,unlimdimid))
                !print *,'ncid,ndimension, nVars, nGlobalAtts,unlimdimid',ncid,ndimension, nVars, nGlobalAtts,unlimdimid

                !; the dimid in netcdf-files is counted from 1 to ndimension
                !; look for the name and length of the dimension with f90_inquire_dimension
                !; nf90_inquire_dimension input: ncid, dimid; nf90_inquire_dimension output: name, length
                do dimid=1,ndimension
                    print *,'dimension loop dimid ',dimid
                    call check_netcdf( nf90_inquire_dimension(ncid,dimid, dimname, length) )
                            print*, 'ncid,dimid, dimname, length',ncid,dimid, trim(dimname), length
                             if ( trim(dimname) == 'lon') nlon_ahf=length          ! here I know that the name of zonal dimensio &
!& n is 'lon'
                             if ( trim(dimname) == 'lat') nlat_ahf=length          ! here I know that the name of meridional dim &
!& ension is 'lat'
                enddo


       END SUBROUTINE get_dimension_AHF_data


       
        !> read coordinates for AHF raw data from netcdf file
       SUBROUTINE get_AHF_data_coordinates(ncid,          &
                                            nlon_ahf,     &
                                            nlat_ahf,     &  
                                            startlon_ahf, &
                                            startlat_ahf, &
                                            dlon_ahf,     &
                                            dlat_ahf,     &
                                            lon,           &
                                            lat)

        INTEGER, INTENT(in) :: ncid                             !< netcdf unit file number
        INTEGER (KIND=i4), INTENT(in) :: nlon_ahf !< number of grid elements in zonal direction for AHF data
        INTEGER (KIND=i4), INTENT(in) :: nlat_ahf !< number of grid elements in meridional direction for AHF data

        REAL (KIND=wp), INTENT(out) :: startlon_ahf !< longitude of lower left grid element for AHF data 
        REAL (KIND=wp), INTENT(out) :: startlat_ahf !< latitude of lower left grid element for AHF data

        REAL (KIND=wp), INTENT(out) :: dlon_ahf !< grid point distance in zonal direction (in degrees) for AHF data
        REAL (KIND=wp), INTENT(out) :: dlat_ahf !< grid point distance in meridional direction (in degrees) for AHF data

        REAL (KIND=wp), INTENT(out) :: lon(1:nlon_ahf)      !< longitude of ahf raw data in geographical system
        REAL (KIND=wp), INTENT(out) :: lat(1:nlat_ahf)      !< latitude of ahf raw date in geographical system
        


        !local variables
        INTEGER :: ndimension                       !< number of dimensions in netcdf file
        INTEGER :: nVars                            !< number of variables in netcdf file
        INTEGER :: nGlobalAtts                      !< number of gloabal Attributes in netcdf file
        INTEGER :: unlimdimid                       !< id of unlimited dimension (e.g. time) in netcdf file

        INTEGER :: dimid                            !< id of dimension
        CHARACTER (len=80) :: dimname               !< name of dimension
        INTEGER :: length                           !< length of dimension

        INTEGER :: varid                            !< id of variable
        CHARACTER (len=80) :: varname               !< name of variable
        INTEGER :: xtype                            !< netcdf type of variable/attribute
        INTEGER :: ndim                             !< number of dimensions of variable
        INTEGER, ALLOCATABLE :: var_dimids(:)       !< id of variable dimensions, vector, maximal dimension ndimension
        INTEGER :: nAtts                            !< number of attributes for a netcdf variable
        
        INTEGER :: errorcode                        !< error status variable
        
        INTEGER :: attnum                           !< counter for attribute number
        CHARACTER (len=80) :: attname               !< name of attribute
        CHARACTER (len=80) :: attributetext         !< attribute text

        !! look for numbers of dimensions, Variable, Attributes, and the dimid of a possible unlimited dimension (probably time)
        !; nf90_inquire input: ncid; nf90_inquire output: ndimension, nVars, nGlobalAtts,unlimdimid
        CALL check_netcdf (nf90_inquire(ncid,ndimension, nVars, nGlobalAtts,unlimdimid))
        !print *,'ncid,ndimension, nVars, nGlobalAtts,unlimdimid',ncid,ndimension, nVars, nGlobalAtts,unlimdimid
        ALLOCATE (var_dimids(ndimension), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array var_dimids')
        var_dimids = 0
        
        !; the varid in netcdf files is counted from 1 to nVars
        !; look for variable names, number of dimension, var_dimids etc with nf90_inquire_variable
        !; nf90_inquire_variable input: ncid, varid; nf90_inquire_variable output name xtype, ndim, var_dimids, nAtts
        !print *,'nVars ',nVars
        variables: DO varid=1,nVars
          !print *,'variable loop, varid ',varid
          CALL check_netcdf(nf90_inquire_variable(ncid,varid,varname,xtype, ndim, var_dimids, nAtts))
          !print *,'------------------------------------'
          !print *,'ncid,varid,vrname',ncid,varid,trim(varname)
          !print *, 'xtype',xtype
          !print *, 'ndim', ndim
          !print *, 'var_dimids', var_dimids
          !print *, 'nAtts', nAtts
          !print *,'------------------------------------'
          
          IF (nAtts.gt.0) THEN
            ! get Attribute name
            DO attnum=1,nAtts
              ! input nf90_inq_attname: ncid, varid, attnum, output nf90_inq_attname name
              CALL check_netcdf ( nf90_inq_attname(ncid, varid, attnum, attname ))
              ! print *,'ncid,varid,attnumm, attname',ncid,varid,attnum,trim(attname)
              ! input nf90_inquire_attribute: ncid, varid, name; nf90_inquire_attribute output: xtype, length, attnum
              CALL check_netcdf(nf90_inquire_attribute(ncid, varid, attname, xtype, length))
              !print *, 'ncid, varid, attname, xtype, length, attnum',ncid, varid, trim(attname),' ', xtype, length, attnum
              ! nf90_get_att input: ncid, varid, name; nf90_get_att ou put: attributetext
              ! note, attributetext should be the right type of variable
              getattribute: SELECT CASE (xtype)
              CASE (NF90_CHAR) ! for character attributes
                CALL check_netcdf(nf90_get_att(ncid, varid, attname, attributetext))
                ! print *,'get attribute: ncid, varid, attname, attributetext',ncid, varid, trim(attname),' ', trim(attributetext)
              END SELECT getattribute
            END DO ! done with attributes
          ENDIF
          
          IF (trim(varname) == 'lon') THEN             ! here I know that the variable with longitude coordinates is called 'lon'

            !print *,'------------------------------------'
            !print *,'ncid,varid,vrname',ncid,varid,trim(varname)
            !print *, 'xtype',xtype
            !print *, 'ndim', ndim
            !print *, 'var_dimids', var_dimids
            !print *, 'nAtts', nAtts
            !print *,'------------------------------------'
            
            CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(1), dimname, len=length) )
            if (length /= nlon_ahf)  CALL abort_extpar('nlon_ahf is not equal data dimension')
            CALL check_netcdf(nf90_get_var(ncid,varid,lon,start=(/1/),count=(/length/))) ! read from netcdf file into lon(:)
          ENDIF

          IF (trim(varname) == 'lat') THEN              ! here I know that the variable with latitude coordinates is called 'lat'

            !print *,'------------------------------------'
            !print *,'ncid,varid,vrname',ncid,varid,trim(varname)
            !print *, 'xtype',xtype
            !print *, 'ndim', ndim
            !print *, 'var_dimids', var_dimids
            !print *, 'nAtts', nAtts
            !print *,'------------------------------------'
            
            CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(1), dimname, len=length) )
            if (length /= nlat_ahf)  CALL abort_extpar('nlat_ahf is not equal data dimension')
            !print *,'dimt ',dimt
            CALL check_netcdf(nf90_get_var(ncid,varid,lat,start=(/1/),count=(/length/))) ! read from netcdf file into lat(:)
          ENDIF

        ENDDO variables

        startlon_ahf = lon(1) ! longitude of the upper left grid element
        startlat_ahf = lat(1) ! latitude of the upper left grid element

        dlon_ahf = (lon(nlon_ahf)-lon(1))/float(nlon_ahf-1) ! dlon_ahf in degrees
        dlat_ahf = (lat(1)-lat(nlat_ahf))/float(nlat_ahf-1) ! dlat_ahf in degrees
        ! dlon_ahf = 360._wp/float(nlon_ahf) ! dlon_ahf in degrees
        ! dlat_ahf = 180._wp/float(nlat_ahf) ! dlat_ahf in degrees
        
        
        
        
      END SUBROUTINE get_AHF_data_coordinates

       
   !> get one row of AHF raw data from netcdf file (along zonal direction)
       SUBROUTINE get_one_row_AHF_data(ncid,                &
                                        nlon_ahf,           & 
                                        nlat_ahf,           &
                                        row_index,           &
                                        ahf_raw_data_lonrow)

        INTEGER, INTENT(in) :: ncid                             !< netcdf unit file number
        INTEGER (KIND=i4), INTENT(in) :: nlon_ahf       !< number of grid elements in zonal direction for AHF data
        INTEGER (KIND=i4), INTENT(in) :: nlat_ahf       !< number of grid elements in meridional direction for AHF data

        INTEGER (KIND=i4), INTENT(in) :: row_index       !< the index of the data row to read in

        REAL (KIND=wp), INTENT(out) :: ahf_raw_data_lonrow(1:nlon_ahf)      !< longitude of ahf raw data in geographical syste &
!& m
        

        !local variables
        INTEGER :: ndimension                       !< number of dimensions in netcdf file
        INTEGER :: nVars                            !< number of variables in netcdf file
        INTEGER :: nGlobalAtts                      !< number of gloabal Attributes in netcdf file
        INTEGER :: unlimdimid                       !< id of unlimited dimension (e.g. time) in netcdf file

        INTEGER :: dimid                            !< id of dimension
        CHARACTER (len=80) :: dimname               !< name of dimension
        INTEGER :: length                           !< length of dimension
        INTEGER :: dim_lon                           !< length of dimension lon
        INTEGER :: dim_lat                           !< length of dimension lat

        

        INTEGER :: varid                            !< id of variable
        CHARACTER (len=80) :: varname               !< name of variable
        INTEGER :: xtype                            !< netcdf type of variable/attribute
        INTEGER :: ndim                             !< number of dimensions of variable
        INTEGER, ALLOCATABLE :: var_dimids(:)       !< id of variable dimensions, vector, maximal dimension ndimension
        INTEGER :: nAtts                            !< number of attributes for a netcdf variable
        
        INTEGER :: errorcode                        !< error status variable
        
        INTEGER :: attnum                           !< counter for attribute number
        CHARACTER (len=80) :: attname               !< name of attribute
        CHARACTER (len=80) :: attributetext         !< attribute text

        !! look for numbers of dimensions, Variable, Attributes, and the dimid of a possible unlimited dimension (probably time)
        !; nf90_inquire input: ncid; nf90_inquire output: ndimension, nVars, nGlobalAtts,unlimdimid
        CALL check_netcdf (nf90_inquire(ncid,ndimension, nVars, nGlobalAtts,unlimdimid))
        !print *,'ncid,ndimension, nVars, nGlobalAtts,unlimdimid',ncid,ndimension, nVars, nGlobalAtts,unlimdimid
        ALLOCATE (var_dimids(ndimension), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array var_dimids')
        var_dimids = 0
        
        !; the varid in netcdf files is counted from 1 to nVars
        !; look for variable names, number of dimension, var_dimids etc with nf90_inquire_variable
        !; nf90_inquire_variable input: ncid, varid; nf90_inquire_variable output name xtype, ndim, var_dimids, nAtts
        !print *,'nVars ',nVars
        variables: DO varid=1,nVars
          !print *,'variable loop, varid ',varid
          CALL check_netcdf(nf90_inquire_variable(ncid,varid,varname,xtype, ndim, var_dimids, nAtts))
          !print *,'------------------------------------'
          !print *,'ncid,varid,vrname',ncid,varid,trim(varname)
          !print *, 'xtype',xtype
          !print *, 'ndim', ndim
          !print *, 'var_dimids', var_dimids
          !print *, 'nAtts', nAtts
          !print *,'------------------------------------'

          IF (nAtts.gt.0) THEN
            ! get Attribute name
            DO attnum=1,nAtts
              ! input nf90_inq_attname: ncid, varid, attnum, output nf90_inq_attname name
              CALL check_netcdf ( nf90_inq_attname(ncid, varid, attnum, attname ))
              ! print *,'ncid,varid,attnumm, attname',ncid,varid,attnum,trim(attname)
              ! input nf90_inquire_attribute: ncid, varid, name; nf90_inquire_attribute output: xtype, length, attnum
              CALL check_netcdf(nf90_inquire_attribute(ncid, varid, attname, xtype, length))
              !print *, 'ncid, varid, attname, xtype, length, attnum',ncid, varid, trim(attname),' ', xtype, length, attnum
              ! nf90_get_att input: ncid, varid, name; nf90_get_att ou put: attributetext
              ! note, attributetext should be the right type of variable
              getattribute: SELECT CASE (xtype)
              CASE (NF90_CHAR) ! for character attributes
                CALL check_netcdf(nf90_get_att(ncid, varid, attname, attributetext))
                ! print *,'get attribute: ncid, varid, attname, attributetext',ncid, varid, trim(attname),' ', trim(attributetext)
              END SELECT getattribute
            END DO ! done with attributes
          ENDIF
          
          
          IF (trim(varname) == 'AHF') THEN            ! here I know that the variable with latitude coordinates is called 'AHF'

            !print *,'------------------------------------'
            !print *,'ncid,varid,vrname',ncid,varid,trim(varname)
            !print *, 'xtype',xtype
            !print *, 'ndim', ndim
            !print *, 'var_dimids', var_dimids
            !print *, 'nAtts', nAtts
            !print *,'------------------------------------'
            
            CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(1), dimname, len=dim_lon) )
            ! print *, 'var_dimids(1), dimname, dim_lon: ', var_dimids(1), TRIM(dimname), dim_lon
            if (dim_lon /= nlon_ahf)   CALL abort_extpar('nlon_ahf is not equal data dimension')
            
            CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(2), dimname, len=dim_lat) )
            ! print *, 'var_dimids(2), dimname, dim_lat: ', var_dimids(2), TRIM(dimname), dim_lat
            if (dim_lat /= nlat_ahf)   CALL abort_extpar('nlat_ahf is not equal data dimension')

            !print *,'dimt ',dimt
            CALL check_netcdf(nf90_get_var(ncid,varid,ahf_raw_data_lonrow,                          &
                 start=(/1,row_index/),count=(/nlon_ahf,1/))) ! read from netcdf file into ahf_raw_data_lonrow(:)
          ENDIF

        ENDDO variables

       END SUBROUTINE get_one_row_AHF_data

       
   !> get a block of AHF raw data from netcdf file (given startrow, endrow, startcolumn, endcolumn)
       SUBROUTINE get_block_AHF_data(ncid,                &
                                      nlon_ahf,           & 
                                      nlat_ahf,           &
                                      startcolumn_index,   &
                                      endcolumn_index,     &
                                      startrow_index,      &
                                      endrow_index,        &
                                      ncolumns,            &
                                      nrows,               &
                                      ahf_data_block)

        INTEGER, INTENT(in) :: ncid                             !< netcdf unit file number
        INTEGER (KIND=i4), INTENT(in) :: nlon_ahf       !< number of grid elements in zonal direction for AHF data
        INTEGER (KIND=i4), INTENT(in) :: nlat_ahf       !< number of grid elements in meridional direction for AHF data
        INTEGER (KIND=i4), INTENT(in) :: startcolumn_index    !< the index of the startcolumn of data to read in
        INTEGER (KIND=i4), INTENT(in) :: endcolumn_index      !< the index of the endcolumn of data to read in
        INTEGER (KIND=i4), INTENT(in) :: startrow_index       !< the index of the startrow of data to read in
        INTEGER (KIND=i4), INTENT(in) :: endrow_index         !< the index of the endrow of data to read in
        INTEGER (KIND=i4), INTENT(in) ::   ncolumns           !< number of columns of data block
        INTEGER (KIND=i4), INTENT(in) ::   nrows              !< number of rows of data block
        

        REAL (KIND=wp), INTENT(out) :: ahf_data_block(1:ncolumns,1:nrows)      !< longitude of ahf raw data in geographical sys &
!& tem
        

        !local variables
        INTEGER :: ndimension                       !< number of dimensions in netcdf file
        INTEGER :: nVars                            !< number of variables in netcdf file
        INTEGER :: nGlobalAtts                      !< number of gloabal Attributes in netcdf file
        INTEGER :: unlimdimid                       !< id of unlimited dimension (e.g. time) in netcdf file

        INTEGER :: dimid                            !< id of dimension
        CHARACTER (len=80) :: dimname               !< name of dimension
        INTEGER :: length                           !< length of dimension
        INTEGER :: dim_lon                           !< length of dimension lon
        INTEGER :: dim_lat                           !< length of dimension lat

        

        INTEGER :: varid                            !< id of variable
        CHARACTER (len=80) :: varname               !< name of variable
        INTEGER :: xtype                            !< netcdf type of variable/attribute
        INTEGER :: ndim                             !< number of dimensions of variable
        INTEGER, ALLOCATABLE :: var_dimids(:)       !< id of variable dimensions, vector, maximal dimension ndimension
        INTEGER :: nAtts                            !< number of attributes for a netcdf variable
        
        INTEGER :: errorcode                        !< error status variable
        
        INTEGER :: attnum                           !< counter for attribute number
        CHARACTER (len=80) :: attname               !< name of attribute
        CHARACTER (len=80) :: attributetext         !< attribute text

        
        IF ( (startcolumn_index < 1) .or. (startcolumn_index > nlon_ahf)) then
          CALL abort_extpar('startcolumn_index out of range')
        ENDIF
        IF ( (endcolumn_index < 1) .or. (endcolumn_index > nlon_ahf)) then
          CALL abort_extpar('endcolumn_index out of range')
        ENDIF
        IF ( (startrow_index < 1) .or. (startrow_index > nlat_ahf)) then
          CALL abort_extpar('startrow_index out of range')
        ENDIF
        
        IF ( (endrow_index < 1) .or. (endrow_index > nlat_ahf)) then
          CALL abort_extpar('endrow_index out of range')
        ENDIF
        
        !ncolumns = endcolumn_index - startcolumn_index + 1
        !nrows    = endrow_index - startrow_index + 1

        
        !! look for numbers of dimensions, Variable, Attributes, and the dimid of a possible unlimited dimension (probably time)
        !; nf90_inquire input: ncid; nf90_inquire output: ndimension, nVars, nGlobalAtts,unlimdimid
        CALL check_netcdf (nf90_inquire(ncid,ndimension, nVars, nGlobalAtts,unlimdimid))
        !print *,'ncid,ndimension, nVars, nGlobalAtts,unlimdimid',ncid,ndimension, nVars, nGlobalAtts,unlimdimid
        ALLOCATE (var_dimids(ndimension), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array var_dimids')
        var_dimids = 0

        !; the varid in netcdf files is counted from 1 to nVars
        !; look for variable names, number of dimension, var_dimids etc with nf90_inquire_variable
        !; nf90_inquire_variable input: ncid, varid; nf90_inquire_variable output name xtype, ndim, var_dimids, nAtts
        !print *,'nVars ',nVars
        variables: DO varid=1,nVars
          !print *,'variable loop, varid ',varid
          CALL check_netcdf(nf90_inquire_variable(ncid,varid,varname,xtype, ndim, var_dimids, nAtts))
          !print *,'------------------------------------'
          !print *,'ncid,varid,vrname',ncid,varid,trim(varname)
          !print *, 'xtype',xtype
          !print *, 'ndim', ndim
          !print *, 'var_dimids', var_dimids
          !print *, 'nAtts', nAtts
          !print *,'------------------------------------'

          IF (nAtts.gt.0) THEN
            ! get Attribute name
            DO attnum=1,nAtts
              ! input nf90_inq_attname: ncid, varid, attnum, output nf90_inq_attname name
              CALL check_netcdf ( nf90_inq_attname(ncid, varid, attnum, attname ))
              ! print *,'ncid,varid,attnumm, attname',ncid,varid,attnum,trim(attname)
              ! input nf90_inquire_attribute: ncid, varid, name; nf90_inquire_attribute output: xtype, length, attnum
              CALL check_netcdf(nf90_inquire_attribute(ncid, varid, attname, xtype, length))
              !print *, 'ncid, varid, attname, xtype, length, attnum',ncid, varid, trim(attname),' ', xtype, length, attnum
              ! nf90_get_att input: ncid, varid, name; nf90_get_att ou put: attributetext
              ! note, attributetext should be the right type of variable
              getattribute: SELECT CASE (xtype)
              CASE (NF90_CHAR) ! for character attributes
                CALL check_netcdf(nf90_get_att(ncid, varid, attname, attributetext))
                ! print *,'get attribute: ncid, varid, attname, attributetext',ncid, varid, trim(attname),' ', trim(attributetext)
              END SELECT getattribute
            END DO ! done with attributes
          ENDIF

          IF (trim(varname) == 'AHF') THEN           ! here I know that the variable with latitude coordinates is called 'AHF'

            !print *,'------------------------------------'
            !print *,'ncid,varid,vrname',ncid,varid,trim(varname)
            !print *, 'xtype',xtype
            !print *, 'ndim', ndim
            !print *, 'var_dimids', var_dimids
            !print *, 'nAtts', nAtts
            !print *,'------------------------------------'
            
            CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(1), dimname, len=dim_lon) )
            ! print *, 'var_dimids(1), dimname, dim_lon: ', var_dimids(1), TRIM(dimname), dim_lon
            if (dim_lon /= nlon_ahf)   CALL abort_extpar('nlon_ahf is not equal data dimension')
            
            CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(2), dimname, len=dim_lat) )
            ! print *, 'var_dimids(2), dimname, dim_lat: ', var_dimids(2), TRIM(dimname), dim_lat
            if (dim_lat /= nlat_ahf)   CALL abort_extpar('nlat_ahf is not equal data dimension')

            !print *,'dimt ',dimt
            ! read from netcdf file in to ahf_raw_data_lonrow(:)
            CALL check_netcdf(nf90_get_var(ncid,varid,ahf_data_block,                          &
                 start=(/startcolumn_index,startrow_index/),count=(/ncolumns,nrows/))) 
          ENDIF

        ENDDO variables



       END SUBROUTINE get_block_AHF_data

       !> get a pixel of AHF raw data from netcdf file (given grid index)
       SUBROUTINE get_pixel_AHF_data(ncid,                &
                                      nlon_ahf,           & 
                                      nlat_ahf,           &
                                      column_index,        &
                                      row_index,           &
                                      ahf_pixel_data)

        INTEGER, INTENT(in) :: ncid                             !< netcdf unit file number
        INTEGER (KIND=i4), INTENT(in) :: nlon_ahf       !< number of grid elements in zonal direction for AHF data
        INTEGER (KIND=i4), INTENT(in) :: nlat_ahf       !< number of grid elements in meridional direction for AHF data
        INTEGER (KIND=i4), INTENT(in) :: column_index    !< the index of the column of data to read in
        INTEGER (KIND=i4), INTENT(in) :: row_index       !< the index of the trow of data to read in

        REAL (KIND=wp), INTENT(out) :: ahf_pixel_data      !< value of ahf raw data pixel
        

        !local variables
        INTEGER :: ndimension                       !< number of dimensions in netcdf file
        INTEGER :: nVars                            !< number of variables in netcdf file
        INTEGER :: nGlobalAtts                      !< number of gloabal Attributes in netcdf file
        INTEGER :: unlimdimid                       !< id of unlimited dimension (e.g. time) in netcdf file

        INTEGER :: dimid                            !< id of dimension
        CHARACTER (len=80) :: dimname               !< name of dimension
        INTEGER :: length                           !< length of dimension
        INTEGER :: dim_lon                          !< length of dimension lon
        INTEGER :: dim_lat                          !< length of dimension lat

        

        INTEGER :: varid                            !< id of variable
        CHARACTER (len=80) :: varname               !< name of variable
        INTEGER :: xtype                            !< netcdf type of variable/attribute
        INTEGER :: ndim                             !< number of dimensions of variable
        INTEGER, ALLOCATABLE :: var_dimids(:)       !< id of variable dimensions, vector, maximal dimension ndimension
        INTEGER :: nAtts                            !< number of attributes for a netcdf variable
        
        INTEGER :: errorcode                        !< error status variable
        
        INTEGER :: attnum                           !< counter for attribute number
        CHARACTER (len=80) :: attname               !< name of attribute
        CHARACTER (len=80) :: attributetext         !< attribute text

        
        IF ( (column_index < 1) .or. (column_index > nlon_ahf)) then
          CALL abort_extpar('column_index out of range')
        ENDIF
        IF ( (row_index < 1) .or. (row_index > nlat_ahf)) then
          CALL abort_extpar('row_index out of range')
        ENDIF
        
        !! look for numbers of dimensions, Variable, Attributes, and the dimid of a possible unlimited dimension (probably time)
        !; nf90_inquire input: ncid; nf90_inquire output: ndimension, nVars, nGlobalAtts,unlimdimid
        CALL check_netcdf (nf90_inquire(ncid,ndimension, nVars, nGlobalAtts,unlimdimid))
        !print *,'ncid,ndimension, nVars, nGlobalAtts,unlimdimid',ncid,ndimension, nVars, nGlobalAtts,unlimdimid
        ALLOCATE (var_dimids(ndimension), STAT=errorcode)
        IF(errorcode.NE.0) CALL abort_extpar('Cant allocate the array var_dimids')
        var_dimids = 0
        

        !; the varid in netcdf files is counted from 1 to nVars
        !; look for variable names, number of dimension, var_dimids etc with nf90_inquire_variable
        !; nf90_inquire_variable input: ncid, varid; nf90_inquire_variable output name xtype, ndim, var_dimids, nAtts
        !print *,'nVars ',nVars
        variables: DO varid=1,nVars
          !print *,'variable loop, varid ',varid
          CALL check_netcdf(nf90_inquire_variable(ncid,varid,varname,xtype, ndim, var_dimids, nAtts))
          !print *,'------------------------------------'
          !print *,'ncid,varid,vrname',ncid,varid,trim(varname)
          !print *, 'xtype',xtype
          !print *, 'ndim', ndim
          !print *, 'var_dimids', var_dimids
          !print *, 'nAtts', nAtts
          !print *,'------------------------------------'
          
          IF (trim(varname) == 'AHF') THEN   ! here I know that the variable with latitude coordinates is called 'AHF'

            !print *,'------------------------------------'
            !print *,'ncid,varid,vrname',ncid,varid,trim(varname)
            !print *, 'xtype',xtype
            !print *, 'ndim', ndim
            !print *, 'var_dimids', var_dimids
            !print *, 'nAtts', nAtts
            !print *,'------------------------------------'
            
            CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(1), dimname, len=dim_lon) )
            ! print *, 'var_dimids(1), dimname, dim_lon: ', var_dimids(1), TRIM(dimname), dim_lon
            if (dim_lon /= nlon_ahf)   CALL abort_extpar('nlon_ahf is not equal data dimension')
            
            CALL check_netcdf( nf90_inquire_dimension(ncid,var_dimids(2), dimname, len=dim_lat) )
            ! print *, 'var_dimids(2), dimname, dim_lat: ', var_dimids(2), TRIM(dimname), dim_lat
            if (dim_lat /= nlat_ahf)   CALL abort_extpar('nlat_ahf is not equal data dimension')

            !print *,'dimt ',dimt
            CALL check_netcdf(nf90_get_var(ncid,varid,ahf_pixel_data,                      &
                 start=(/column_index,row_index/))) ! read from netcdf file into ahf_pixel_data
          ENDIF
          

        ENDDO variables


        
      END SUBROUTINE get_pixel_AHF_data



END MODULE mo_ahf_routines
