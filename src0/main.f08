program Main
  use Mod
  use Mod2

  implicit none

  integer     :: a, sub1
  type(type)  :: t
  type(type2) :: t2

  a    = 1
  sub1 = 0
  t = type(10)
  t2%v = a * 3
  print *, a, t, t2, sub1

  call t%new
  print *, a, t, t2, sub1
end program Main
