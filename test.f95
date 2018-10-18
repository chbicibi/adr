program Test
  implicit none

  real(8), allocatable :: array(:)
  integer :: i

  array = [(rand(), i = 1, 1000000)]

  print "(5f10.5)", array(1:10)
  array =  qsort(array)
  print *
  print "(5f10.5)", array(1:1000:100)

  contains

  recursive function qsort(src) result (res)
    real(8), allocatable :: res(:)
    real(8), intent(in)  :: src(:)
    integer              :: len
    real(8)              :: pivot

    len = size(src)

    if (len <= 1) then
      res = src
    else
      pivot = src(len / 2)
      res   = [qsort(pack(src, src <  pivot)), &
                     pack(src, src == pivot),  &
               qsort(pack(src, src >  pivot))]
    end if
  end function qsort

  real(8) function rand(range, center)
    real(8), intent(in), optional :: range, center
    real(8)                       :: r = 1d0

    call random_number(rand)

    if (present(range)) then
      rand = rand * range
      r        = range
    end if

    if (present(center)) rand = rand - 0.5d0 * r + center
  end function rand
end program Test
