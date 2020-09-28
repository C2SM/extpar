module mo_array_cache

  use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_sizeof, c_associated, c_null_char
  use mo_util_mmap_cache,          only: allocate_cache, deallocate_cache  
  use mo_kind,                     only: i4, i8, sp, dp
  
  implicit none
  
  private

  public :: variable
  public :: allocate_cached
  public :: deallocate_cached  

  integer, parameter :: of_real_type    = 0
  integer, parameter :: of_integer_type = 1  
  
  type variable
    integer :: hash
    character(len=64) :: name
    integer :: sizes(5)
    type(c_ptr) :: ptr_for_deallocation
    integer :: real_or_integer_type
    integer :: kind_type

! preprocessing flag needed because of a compiler bug
! remove once bug is fixed
#if 0
  contains
    procedure :: write_variable_cache_formatted
    generic   :: write(formatted) => write_variable_cache_formatted
#endif
  end type variable

  !! allocate vector in blocks of this size (> 0)
  integer, parameter :: chunk_size = 256 
  type(variable), allocatable ::  cached_variables(:)
  integer :: number_of_cached_variables = 0
  integer :: index_of_cached_variable = 0
  
  interface allocate_cached
    module procedure allocate_3d_sp_cached
    module procedure allocate_4d_sp_cached
    module procedure allocate_5d_sp_cached    
    module procedure allocate_3d_dp_cached
    module procedure allocate_4d_dp_cached
    module procedure allocate_5d_dp_cached    
    module procedure allocate_3d_i4_cached
    module procedure allocate_4d_i4_cached
    module procedure allocate_5d_i4_cached    
  end interface allocate_cached
  
contains

  subroutine allocate_3d_sp_cached(variable_name, variable_pointer, variable_size)

    character(len=*), intent(in) :: variable_name
    real(sp), pointer :: variable_pointer(:,:,:) 
    integer, intent(in) :: variable_size(3)

    integer :: size_to_save(5)

    type(c_ptr) :: addr
    type(variable) :: variable_to_cache
    
    if (any(variable_size == 0)) then
      allocate(variable_pointer(variable_size(1), variable_size(2), variable_size(3)))
    else
      addr = allocate_cache(trim(variable_name)//c_null_char, int(product(variable_size)*c_sizeof(0.0_sp),i8))      
      call c_f_pointer(addr, variable_pointer, shape=variable_size)

      size_to_save(:) = 0
      size_to_save(1:3) = variable_size

      index_of_cached_variable = index_of_cached_variable + 1      
      variable_to_cache = variable(index_of_cached_variable, &
           &                       variable_name,            &
           &                       size_to_save,             &
           &                       addr,                     &
           &                       of_real_type,             &
           &                       sp)
      call add_to_variable_cache(cached_variables, variable_to_cache, number_of_cached_variables)
    endif

  end subroutine allocate_3d_sp_cached

  subroutine allocate_3d_dp_cached(variable_name, variable_pointer, variable_size)

    character(len=*), intent(in) :: variable_name
    real(dp), pointer :: variable_pointer(:,:,:) 
    integer, intent(in) :: variable_size(3)

    integer :: size_to_save(5)
    
    type(c_ptr) :: addr
    type(variable) :: variable_to_cache
    
    if (any(variable_size == 0)) then
      allocate(variable_pointer(variable_size(1), variable_size(2), variable_size(3)))
    else
      addr = allocate_cache(trim(variable_name)//c_null_char, int(product(variable_size)*c_sizeof(0.0_dp),i8))
      call c_f_pointer(addr, variable_pointer, shape=variable_size)

      size_to_save(:) = 0
      size_to_save(1:3) = variable_size
      
      index_of_cached_variable = index_of_cached_variable + 1      
      variable_to_cache = variable(index_of_cached_variable, &
           &                       variable_name,            &
           &                       size_to_save,             &
           &                       addr,                     &
           &                       of_real_type,             &
           &                       dp)
      call add_to_variable_cache(cached_variables, variable_to_cache, number_of_cached_variables)
    endif

  end subroutine allocate_3d_dp_cached

    subroutine allocate_3d_i4_cached(variable_name, variable_pointer, variable_size)

    character(len=*), intent(in) :: variable_name
    integer(i4), pointer :: variable_pointer(:,:,:) 
    integer, intent(in) :: variable_size(3)

    integer :: size_to_save(5)
    
    type(c_ptr) :: addr
    type(variable) :: variable_to_cache
    
    if (any(variable_size == 0)) then
      allocate(variable_pointer(variable_size(1), variable_size(2), variable_size(3)))
    else
      addr = allocate_cache(trim(variable_name)//c_null_char, int(product(variable_size)*c_sizeof(0_i4),i8))
      call c_f_pointer(addr, variable_pointer, shape=variable_size)

      size_to_save(:) = 0
      size_to_save(1:3) = variable_size
      
      index_of_cached_variable = index_of_cached_variable + 1      
      variable_to_cache = variable(index_of_cached_variable, &
           &                       variable_name,            &
           &                       size_to_save,             &
           &                       addr,                     &
           &                       of_integer_type,          &
           &                       i4)
      call add_to_variable_cache(cached_variables, variable_to_cache, number_of_cached_variables)
    endif

  end subroutine allocate_3d_i4_cached

  subroutine allocate_4d_sp_cached(variable_name, variable_pointer, variable_size)

    character(len=*), intent(in) :: variable_name
    real(sp), pointer :: variable_pointer(:,:,:,:) 
    integer, intent(in) :: variable_size(4)

    integer :: size_to_save(5)
    
    type(c_ptr) :: addr
    type(variable) :: variable_to_cache
    
    if (any(variable_size == 0)) then
      allocate(variable_pointer(variable_size(1), variable_size(2), variable_size(3), variable_size(4)))
    else
      addr = allocate_cache(trim(variable_name)//c_null_char, int(product(variable_size)*c_sizeof(0.0_sp),i8))
      call c_f_pointer(addr, variable_pointer, shape=variable_size)
      
      size_to_save(:) = 0
      size_to_save(1:4) = variable_size
      
      index_of_cached_variable = index_of_cached_variable + 1      
      variable_to_cache = variable(index_of_cached_variable, &
           &                       variable_name,            &
           &                       size_to_save,             &
           &                       addr,                     &
           &                       of_real_type,             &
           &                       sp)
      call add_to_variable_cache(cached_variables, variable_to_cache, number_of_cached_variables)
    endif
    
  end subroutine allocate_4d_sp_cached

  subroutine allocate_4d_dp_cached(variable_name, variable_pointer, variable_size)

    character(len=*), intent(in) :: variable_name
    real(dp), pointer :: variable_pointer(:,:,:,:) 
    integer, intent(in) :: variable_size(4)

    integer :: size_to_save(5)
    
    type(c_ptr) :: addr
    type(variable) :: variable_to_cache
    
    if (any(variable_size == 0)) then
      allocate(variable_pointer(variable_size(1), variable_size(2), variable_size(3), variable_size(4)))
    else
      addr = allocate_cache(trim(variable_name)//c_null_char, int(product(variable_size)*c_sizeof(0.0_dp),i8))
      call c_f_pointer(addr, variable_pointer, shape=variable_size)
      
      size_to_save(:) = 0
      size_to_save(1:4) = variable_size
      
      index_of_cached_variable = index_of_cached_variable + 1      
      variable_to_cache = variable(index_of_cached_variable, &
           &                       variable_name,            &
           &                       size_to_save,             &
           &                       addr,                     &
           &                       of_real_type,             &
           &                       dp)
      call add_to_variable_cache(cached_variables, variable_to_cache, number_of_cached_variables)
    endif
    
  end subroutine allocate_4d_dp_cached

    subroutine allocate_4d_i4_cached(variable_name, variable_pointer, variable_size)

    character(len=*), intent(in) :: variable_name
    integer(i4), pointer :: variable_pointer(:,:,:,:) 
    integer, intent(in) :: variable_size(4)

    integer :: size_to_save(5)
    
    type(c_ptr) :: addr
    type(variable) :: variable_to_cache
    
    if (any(variable_size == 0)) then
      allocate(variable_pointer(variable_size(1), variable_size(2), variable_size(3), variable_size(4)))
    else
      addr = allocate_cache(trim(variable_name)//c_null_char, int(product(variable_size)*c_sizeof(0_i4),i8))
      call c_f_pointer(addr, variable_pointer, shape=variable_size)
      
      size_to_save(:) = 0
      size_to_save(1:4) = variable_size
      
      index_of_cached_variable = index_of_cached_variable + 1      
      variable_to_cache = variable(index_of_cached_variable, &
           &                       variable_name,            &
           &                       size_to_save,             &
           &                       addr,                     &
           &                       of_integer_type,          &
           &                       i4)
      call add_to_variable_cache(cached_variables, variable_to_cache, number_of_cached_variables)
    endif
    
  end subroutine allocate_4d_i4_cached

  subroutine allocate_5d_sp_cached(variable_name, variable_pointer, variable_size)

    character(len=*), intent(in) :: variable_name
    real(sp), pointer :: variable_pointer(:,:,:,:,:) 
    integer, intent(in) :: variable_size(5)

    integer :: size_to_save(5)
    
    type(c_ptr) :: addr
    type(variable) :: variable_to_cache
    
    if (any(variable_size == 0)) then
      allocate(variable_pointer(variable_size(1), &
           &                    variable_size(2), &
           &                    variable_size(3), &
           &                    variable_size(4), &
           &                    variable_size(5)))
    else
      addr = allocate_cache(trim(variable_name)//c_null_char, int(product(variable_size)*c_sizeof(0.0_sp),i8))
      call c_f_pointer(addr, variable_pointer, shape=variable_size)
      
      size_to_save(:) = variable_size
      
      index_of_cached_variable = index_of_cached_variable + 1      
      variable_to_cache = variable(index_of_cached_variable, &
           &                       variable_name,            &
           &                       size_to_save,             &
           &                       addr,                     &
           &                       of_real_type,             &
           &                       sp)
      call add_to_variable_cache(cached_variables, variable_to_cache, number_of_cached_variables)
    endif
    
  end subroutine allocate_5d_sp_cached

  subroutine allocate_5d_dp_cached(variable_name, variable_pointer, variable_size)

    character(len=*), intent(in) :: variable_name
    real(dp), pointer :: variable_pointer(:,:,:,:,:) 
    integer, intent(in) :: variable_size(5)

    integer :: size_to_save(5)
    
    type(c_ptr) :: addr
    type(variable) :: variable_to_cache
    
    if (any(variable_size == 0)) then
      allocate(variable_pointer(variable_size(1), &
           &                    variable_size(2), &
           &                    variable_size(3), &
           &                    variable_size(4), &
           &                    variable_size(5)))
    else
      addr = allocate_cache(trim(variable_name)//c_null_char, int(product(variable_size)*c_sizeof(0.0_dp),i8))
      call c_f_pointer(addr, variable_pointer, shape=variable_size)
      
      size_to_save(:) = variable_size
      
      index_of_cached_variable = index_of_cached_variable + 1      
      variable_to_cache = variable(index_of_cached_variable, &
           &                       variable_name,            &
           &                       size_to_save,             &
           &                       addr,                     &
           &                       of_real_type,             &
           &                       dp)
      call add_to_variable_cache(cached_variables, variable_to_cache, number_of_cached_variables)
    endif
    
  end subroutine allocate_5d_dp_cached

  subroutine allocate_5d_i4_cached(variable_name, variable_pointer, variable_size)

    character(len=*), intent(in) :: variable_name
    integer(i4), pointer :: variable_pointer(:,:,:,:,:) 
    integer, intent(in) :: variable_size(5)

    integer :: size_to_save(5)
    
    type(c_ptr) :: addr
    type(variable) :: variable_to_cache
    
    if (any(variable_size == 0)) then
      allocate(variable_pointer(variable_size(1), &
           &                    variable_size(2), &
           &                    variable_size(3), &
           &                    variable_size(4), &
           &                    variable_size(5)))
    else
      addr = allocate_cache(trim(variable_name)//c_null_char, int(product(variable_size)*c_sizeof(0_i4),i8))
      call c_f_pointer(addr, variable_pointer, shape=variable_size)
      
      size_to_save(:) = variable_size
      
      index_of_cached_variable = index_of_cached_variable + 1      
      variable_to_cache = variable(index_of_cached_variable, &
           &                       variable_name,            &
           &                       size_to_save,             &
           &                       addr,                     &
           &                       of_integer_type,          &
           &                       i4)
      call add_to_variable_cache(cached_variables, variable_to_cache, number_of_cached_variables)
    endif
    
  end subroutine allocate_5d_i4_cached

  subroutine deallocate_cached(variable_name)
    character(len=*), intent(in) :: variable_name
    integer :: i
    integer(i8) :: v_sizeof

    v_sizeof = -1

    do i = 1, number_of_cached_variables

      select case (cached_variables(i)%real_or_integer_type)
      case(of_real_type)
        select case (cached_variables(i)%kind_type)
        case(sp)
          v_sizeof = c_sizeof(0.0_sp)
        case(dp)
          v_sizeof = c_sizeof(0.0_dp)
        end select
      case(of_integer_type)
        select case (cached_variables(i)%kind_type)
        case(i4)
          v_sizeof = c_sizeof(0_i4)
        end select
      end select

      if (variable_name == cached_variables(i)%name) then
        if (.not. c_associated(cached_variables(i)%ptr_for_deallocation)) then
          call deallocate_cache(cached_variables(i)%ptr_for_deallocation, &
               &                int(size(cached_variables(i)%sizes)*v_sizeof,i8))
          exit
        endif
      endif
    enddo
    
  end subroutine deallocate_cached

  subroutine add_to_variable_cache(vector, val, n, finished)
  
    implicit none
  
    type(variable), allocatable, intent(inout) :: vector(:)
    type(variable), intent(in) :: val
    !! index of last element added to vector
    !! must be initialized to size(vector)
    !! (or 0 if not allocated) before first call
    integer, intent(inout) :: n  

    !! set to true to return vector resized to its final, correct size (n)
    logical, intent(in), optional :: finished 
  
    type(variable), allocatable :: tmp(:)
    logical :: final_resize

    final_resize = .false.
    if (present(finished)) then
      if (finished) final_resize = .true.
    endif

    if (allocated(vector)) then
      if (n == size(vector)) then
        allocate(tmp(size(vector)+chunk_size))
        tmp(1:size(vector)) = vector
        call move_alloc(tmp,vector)
      end if
      n = n + 1
    else
      allocate(vector(chunk_size))
      n = 1
    end if
    
    vector(n) = val

    if (final_resize) then
      if (allocated(tmp)) deallocate(tmp)
      allocate(tmp(n))
      tmp = vector(1:n)
      call move_alloc(tmp,vector)
    end if
    
  end subroutine add_to_variable_cache

! preprocessing flag needed because of a compiler bug
! remove once bug is fixed
#if 0  
  subroutine write_variable_cache_formatted(val, i_unit, i_iotype, v_list, i_iostat, i_iomsg)
    class(variable), intent(in)     :: val
    integer, intent(in)             :: i_unit
    character(len=*), intent(in)    :: i_iotype
    integer, intent(in)             :: v_list(:)
    integer, intent(out)            :: i_iostat
    character(len=*), intent(inout) :: i_iomsg

    i_iostat = 0

    write (i_unit, '(i4.4)',     iostat=i_iostat) val%hash
    write (i_unit, '(1x,a)',     iostat=i_iostat) trim(val%name)
    write (i_unit, '(5(1x,i4))', iostat=i_iostat) val%sizes
    write (i_unit, '(i0)',       iostat=i_iostat) val%ptr_for_deallocation
    write (i_unit, '(i0)',       iostat=i_iostat) val%real_or_integer_type
    write (i_unit, '(i0)',       iostat=i_iostat) val%kind_type
        
  end subroutine write_variable_cache_formatted
#endif
  
end module mo_array_cache
