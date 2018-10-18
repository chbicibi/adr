module Base_Func
  implicit none

  ! private
  ! public :: Hash, SetRandomSeed, Random, Sort, Reverse, Shuffle, CheckDup
  ! public :: Sgn, Factorial, Pick, EulerAngle, Dot, Cross, Norm, Unit, Mul, Angle, Rotate
  ! public :: PI, PI2

  type :: TIAry
    integer, allocatable :: elm(:)
  end type TIAry

  type :: TRAry
    real(8), allocatable :: elm(:)
  end type TRAry

  type :: TRAry2
    real(8), allocatable :: elm(:, :)
  end type TRAry2

  type :: Hash
    integer :: key   = 0
    real(8) :: value = 0d0
  end type Hash

  interface Sgn
    procedure :: Sgn_i, Sgn_d
  end interface Sgn

  interface Random
    procedure :: Random_i, Random_d, Random_d2
  end interface Random

  interface Sort
    procedure :: Sort_i, Sort_d, Sort_h
  end interface Sort

  interface Reverse
    procedure :: Reverse_i, Reverse_d, Reverse_h
  end interface Reverse

  interface Shuffle
    procedure :: Shuffle_i, Shuffle_a
  end interface Shuffle

  interface Angle
    procedure :: Angle_acos, Angle2, Angle3
  end interface Angle

  interface Pick
    procedure :: Pick_i, Pick_d
  end interface Pick

  real(8), parameter :: PI  = acos(-1d0)
  real(8), parameter :: PI2 = 2d0 * PI

  contains

  function Str(n)
    character(:), allocatable :: Str
    integer,      intent(in)  :: n
    integer                   :: i

    i = int(log10((max(abs(n), 1) + 1d-1))) + (3 - Sgn(n)) / 2

    allocate(character(i) :: Str)
    ! allocate(character(int(abs(n) / 10) + (3 - Sgn(n)) / 2) :: Str)
    write(Str, "(i0)") n
  end function Str

  integer function StoI(str_in)
    character(*), intent(in) :: Str_in

    read(Str_in, *) StoI
  end function StoI

  integer function Sgn_i(a)
    integer, intent(in) :: a

    if (a < 0) then
      Sgn_i = -1
    else if (a == 0) then
      Sgn_i = 0
    else
      Sgn_i = 1
    end if
  end function Sgn_i

  integer function Sgn_d(a)
    real(8), intent(in) :: a

    if (a < 0) then
      Sgn_d = -1
    else
      Sgn_d = 1
    end if
  end function Sgn_d

  integer function Factorial(n)
    integer, intent(in) :: n
    integer             :: i

    Factorial = 1
    do i = 2, n
      Factorial = Factorial * i
    end do
  end function Factorial

  ! ----------------------------------------------------------------------------
  ! vector / matrix
  ! ----------------------------------------------------------------------------
  function Cross(a, b)
    real(8)             :: Cross(3)
    real(8), intent(in) :: a(3), b(3)
    integer             :: i

    Cross = [(a(mod(i, 3) + 1) * b(mod(i + 1, 3) + 1) - a(mod(i + 1, 3) + 1) * b(mod(i, 3) + 1), i = 1, 3)]
  end function Cross

  real(8) function Norm(a)
    real(8), intent(in) :: a(:)

    ! Norm = sqrt(sum(max(a ** 2, 0d0)))
    Norm = sqrt(sum(a ** 2))
  end function Norm

  function Unit(a)
    real(8), allocatable :: Unit(:)
    real(8), intent(in)  :: a(:)

    Unit = a / Norm(a)
  end function Unit

  function Rotate(v, axis, theta)
    real(8)             :: Rotate(3)
    real(8), intent(in) :: v(3), axis(3), theta
    real(8)             :: i(3, 3), r(3, 3), m(3, 3)

    i      = reshape([1d0, 0d0, 0d0, 0d0, 1d0, 0d0, 0d0, 0d0, 1d0], shape(i))
    r      = reshape([0d0, axis(3), -axis(2), -axis(3), 0d0, axis(1), axis(2), -axis(1), 0d0], shape(r))
    m      = i + sin(theta) * r + (1d0 - cos(theta)) * matmul(r, r)
    Rotate = matmul(m, v)
  end function Rotate

  function EulerAngle(theta, axis)
    real(8)              :: EulerAngle(3, 3)
    real(8), intent(in)  :: theta(3)
    integer, intent(in)  :: axis(3)
    real(8)              :: s(3), c(3), r(3, 3)
    integer              :: i

    s = sin(theta)
    c = cos(theta)

    do i = 1, 3
      select case (axis(i))
      case (1)
        r = reshape([1d0, 0d0, 0d0, 0d0, c(i), -s(i), 0d0, s(i), c(i)], shape(r))
      case (2)
        r = reshape([c(i), 0d0, s(i), 0d0, 1d0, 0d0, -s(i), 0d0, c(i)], shape(r))
      case (3)
        r = reshape([c(i), -s(i), 0d0, s(i), c(i), 0d0, 0d0, 0d0, 1d0], shape(r))
      end select

      if (i == 1) then
        EulerAngle = r
      else
        EulerAngle = matmul(EulerAngle, r)
      end if
    end do
  end function EulerAngle

  ! ----------------------------------------------------------------------------
  !
  ! ----------------------------------------------------------------------------
  function CheckDup(array1, array2)
    logical             :: CheckDup
    integer, intent(in) :: array1(:), array2(:)
    integer             :: i

    CheckDup = .false.

    do i = 1, size(array2)
      if (any(array1 == array2(i))) then
        CheckDup = .true.
        return
      end if
    end do
  end function CheckDup

  ! ----------------------------------------------------------------------------
  ! random
  ! ----------------------------------------------------------------------------
  real(8) function Random_d(range, center)
    real(8), intent(in), optional :: range, center
    real(8)                       :: r = 1d0

    call random_number(Random_d)

    if (present(range)) then
      Random_d = Random_d * range
      r        = range
    end if

    if (present(center)) Random_d = Random_d - 0.5d0 * r + center
  end function Random_d

  integer function Random_i(max)
    integer, intent(in) :: max

    Random_i = int(Random() * max) + 1
  end function Random_i

  function Random_d2(range, center)
    real(8), allocatable          :: Random_d2(:, :)
    real(8), intent(in)           :: range(:, :)
    real(8), intent(in), optional :: center
    integer                       :: p, q, i, j

    p = ubound(range, 1)
    q = ubound(range, 2)

    Random_d2 = reshape([([(Random(range(i, j), center), i = 1, p)], j = 1, q)], [p, q])
  end function Random_d2

  ! ----------------------------------------------------------------------------
  ! quicksort
  ! ----------------------------------------------------------------------------
  recursive function Sort_i(array_in) result (array_out)
    integer, allocatable          :: array_out(:)
    integer, intent(in)           :: array_in(:)
    integer                       :: size_array, pivot

    size_array = size(array_in)

    if (size_array <= 1) then
      array_out = array_in
    else
      pivot     = array_in(size_array / 2)
      array_out = [Sort_i(pack(array_in, array_in <  pivot)), &
                          pack(array_in, array_in == pivot),  &
                   Sort_i(pack(array_in, array_in >  pivot))]
    end if
  end function Sort_i

  recursive function Sort_d(array_in) result (array_out)
    real(8), allocatable :: array_out(:)
    real(8), intent(in)  :: array_in(:)
    integer              :: size_array
    real(8)              :: pivot

    size_array = size(array_in)

    if (size_array <= 1) then
      array_out = array_in
    else
      pivot     = array_in(size_array / 2)
      array_out = [Sort_d(pack(array_in, array_in <  pivot)), &
                          pack(array_in, array_in == pivot),  &
                   Sort_d(pack(array_in, array_in >  pivot))]
    end if
  end function Sort_d

  recursive function Sort_h(hash_in) result (hash_out)
    type(Hash), allocatable :: hash_out(:)
    type(Hash), intent(in)  :: hash_in(:)
    integer                     :: size_hash
    real(8)                     :: pivot

    size_hash = size(hash_in)

    if (size_hash <= 1) then
      hash_out = hash_in
    else
      pivot    = hash_in(size_hash / 2)%value
      hash_out = [Sort_h(pack(hash_in, hash_in%value <  pivot)), &
                         pack(hash_in, hash_in%value == pivot),  &
                  Sort_h(pack(hash_in, hash_in%value >  pivot))]
    end if
  end function Sort_h

  ! ----------------------------------------------------------------------------
  ! reverse
  ! ----------------------------------------------------------------------------
  function Reverse_i(array)
    integer, allocatable :: Reverse_i(:)
    integer, intent(in)  :: array(:)
    integer              :: i

    Reverse_i = [(array(i), i = size(array), 1, -1)]
  end function Reverse_i

  function Reverse_d(array)
    real(8), allocatable :: Reverse_d(:)
    real(8), intent(in)  :: array(:)
    integer              :: i

    Reverse_d = [(array(i), i = size(array), 1, -1)]
  end function Reverse_d

  function Reverse_h(hash_in)
    type(Hash), allocatable :: Reverse_h(:)
    type(Hash), intent(in)  :: hash_in(:)
    integer                 :: i

    Reverse_h = [(hash_in(i), i = size(hash_in), 1, -1)]
  end function Reverse_h

  ! ----------------------------------------------------------------------------
  ! shuffle
  ! ----------------------------------------------------------------------------
  function Shuffle_i(range, size_in)
    integer, allocatable          :: Shuffle_i(:)
    integer, intent(in)           :: range
    integer, intent(in), optional :: size_in
    integer                       :: size_array, array(range), pos, i

    if (present(size_in) .and. size_in < range) then
      size_array = size_in
    else
      size_array = range
    end if

    Shuffle_i = [(0, i = 1, size_array)]
    array     = [(i, i = 1, range)]

    do i = 1, size_array
      pos          = Random(range - i + 1)
      Shuffle_i(i) = array(pos)
      array(pos)   = array(range - i + 1)
    end do
  end function Shuffle_i

  function Shuffle_a(array, size_in)
    integer, allocatable           :: Shuffle_a(:)
    integer, intent(inout)         :: array(:)
    integer, intent(in),  optional :: size_in
    integer                        :: range, size_array, pos, i

    range = size(array)

    if (present(size_in) .and. size_in < range) then
      size_array = size_in
    else
      size_array = range
    end if

    Shuffle_a  = [(0, i = 1, size_array)]

    do i = 1, size_array
      pos            = Random(range - i + 1)
      Shuffle_a(i)   = array(pos)
      array(pos)     = array(range - i + 1)
    end do
  end function Shuffle_a

  ! ----------------------------------------------------------------------------
  integer function Pick_i(a, p)
    integer, intent(in) :: a(:)
    integer, intent(in) :: p

    Pick_i = a(p)
  end function Pick_i

  real(8) function Pick_d(a, p)
    real(8), intent(in) :: a(:)
    integer, intent(in) :: p

    Pick_d = a(p)
  end function Pick_d

  ! ----------------------------------------------------------------------------
  real(8) function Deg(a)
    real(8), intent(in) :: a

    Deg = a * 360d0 / PI2
  end function Deg

  real(8) function Rad(a)
    real(8), intent(in) :: a

    Rad = a * PI2 / 360d0
  end function Rad

  real(8) function Angle_acos(cos_, sign_)
    real(8), intent(in) :: cos_, sign_

    Angle_acos = modulo(Sgn(sign_) * acos(cos_), PI2)
  end function Angle_acos

  real(8) function Angle2(a, b)
    real(8), intent(in) :: a(3), b(3)

    Angle2 = acos(min(max(dot_product(a, b) / (Norm(a) * Norm(b)), -1d0), 1d0))
  end function Angle2

  real(8) function Angle3(a, b, c)
    real(8), intent(in) :: a(3), b(3), c(3)

    Angle3 = modulo(Sgn(dot_product(Cross(a, b), c)) * Angle(a, b), PI2)
  end function Angle3
end module Base_Func
