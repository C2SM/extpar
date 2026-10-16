MODULE mo_hiressoil_tg_fields
  USE mo_logging
  USE mo_kind, ONLY: wp,i4
  USE mo_grid_structures, ONLY: target_grid_def

  IMPLICIT NONE 

  PRIVATE

  PUBLIC :: hiressoil_mean, &
       hiressoil_min, &
       hiressoil_max, &
       hiressoil_var

 PUBLIC :: allocate_hiressoil_target_fields, &
      deallocate_hiressoil_target_fields
 
  REAL(wp), ALLOCATABLE :: hiressoil_mean(:,:,:)
  REAL(wp), ALLOCATABLE :: hiressoil_min(:,:,:)
  REAL(wp), ALLOCATABLE :: hiressoil_max(:,:,:)
  REAL(wp), ALLOCATABLE :: hiressoil_var(:,:,:)
  

  CONTAINS


  !> allocate fields for GLOBE target data 
  SUBROUTINE allocate_hiressoil_target_fields(tg)

    TYPE(target_grid_def), INTENT(IN) :: tg  !< structure with target grid description
    INTEGER :: errorcode !< error status variable

    ALLOCATE (hiressoil_mean(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array hiressoil_mean',__FILE__,__LINE__)
    hiressoil_mean = 0.0_wp

    ALLOCATE (hiressoil_min(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array hiressoil_min',__FILE__,__LINE__)
    hiressoil_min = 0.0_wp
    
    ALLOCATE (hiressoil_max(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array hiressoil_max',__FILE__,__LINE__)
    hiressoil_max = 0.0_wp

    ALLOCATE (hiressoil_var(1:tg%ie,1:tg%je,1:tg%ke), STAT=errorcode)
    IF(errorcode.NE.0) CALL logging%error('Cant allocate the array hiressoil_var',__FILE__,__LINE__)
    hiressoil_var = 0.0_wp
    
  END SUBROUTINE allocate_hiressoil_target_fields

  !> Deallocate fields after processing one variable
  SUBROUTINE deallocate_hiressoil_target_fields()
    IF (ALLOCATED(hiressoil_mean)) THEN
      DEALLOCATE(hiressoil_mean)
    END IF
    IF (ALLOCATED(hiressoil_min)) THEN
      DEALLOCATE(hiressoil_min)
    END IF
    IF (ALLOCATED(hiressoil_max)) THEN
      DEALLOCATE(hiressoil_max)
    END IF
    IF (ALLOCATED(hiressoil_var)) THEN
      DEALLOCATE(hiressoil_var)
    END IF

    CALL logging%info('hiressoil target fields deallocated')
  END SUBROUTINE deallocate_hiressoil_target_fields
  
END MODULE mo_hiressoil_tg_fields

