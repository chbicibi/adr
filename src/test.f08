subroutine add(a, b)
    implicit none
    integer(8), intent(in) :: a
    integer(8), intent(inout):: b
    b = a + b
end subroutine

subroutine add_array(a, b, N)
    implicit none
    real(8), intent(in) :: a(100)
    real(8), intent(inout):: b(100)
    integer(8), intent(in) :: N
    print *, a(1)
    b(1:N) = a(1:N) + b(1:N)
end subroutine

module test_module
implicit none

contains
    subroutine add(a, b)
        implicit none
        integer(8), intent(in) :: a
        integer(8), intent(inout):: b
        b = a + b
    end subroutine
end module test_module

#define MAXSIZE 1024
subroutine add_array_with_lim(a, b, N)
    implicit none
    integer(8),intent(in) :: N
    real(8), dimension(0:MAXSIZE), intent(in) :: a
    real(8), dimension(0:MAXSIZE), intent(inout):: b
#ifndef NDEBUG
    if (N > MAXSIZE) then
        print *,"MAXSIZE is too small."
    endif
#endif
    b(0:N-1) = b(0:N-1) + a(0:N-1)
end subroutine

! コンパイルコマンド
! gfortran -shared -cpp -o [出力ファイル].dll [*ソースファイル]
