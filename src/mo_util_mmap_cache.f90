module mo_util_mmap_cache

  use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_size_t, c_null_char

  implicit none

  public
  
  !DO NOT DELETE!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !C signature corresponding to the iso_c_binding implementation:
  !
  ! C prototype: void *allocate_cache(char *variable, size_t length);
  !
  !     type(c_ptr) :: addr
  !     real(wp), pointer :: array(:)
  !     addr = allocate_cache(c_char_"<variable name>"//c_null_char, length)
  !     call c_f_pointer(addr, array)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  interface
    function allocate_cache(variable, length) result(addr) bind(c)
      import :: c_ptr, c_char, c_size_t
      type(c_ptr) :: addr
      character(c_char), intent(in) :: variable(*)
      integer(c_size_t), value, intent(in) :: length
    end function allocate_cache
  end interface
 
  ! C prototype: void deallocate_cache(void *addr, size_t length);

  interface
    subroutine deallocate_cache(addr, length) bind(c)
      import :: c_ptr, c_char, c_size_t
      type(c_ptr), value, intent(in) :: addr
      integer(c_size_t), value, intent(in) :: length
    end subroutine deallocate_cache
  end interface

  
end module mo_util_mmap_cache

! #if 0 
! program test_cache
!   use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_loc, &
!        &                                 c_size_t, c_sizeof
!   use mo_util_mmap_cache
!   implicit none

!   type(c_ptr) :: addr
!   real, pointer :: sic(:)

!   addr = allocate_cache("sic", 3*c_sizeof(0.0))
!   call c_f_pointer(addr, sic, [3])

!   sic = [ 0.3, 0.0, 0.1 ]

!   call deallocate_cache(c_loc(sic), 3*c_sizeof(0.0))
  
! end program test_cache
! #endif
