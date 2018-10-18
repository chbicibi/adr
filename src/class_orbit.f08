module Class_Orbit
  use Solver_Orbit
  use Solver_Debri
  use IO_Output

  implicit none

  private
  public :: COrbit
  public :: GetDebriOrbit

  type, abstract :: OrbitFunc
    contains
    procedure, nopass :: Debri => GetDebriOrbit
  end type OrbitFunc

  type, extends(TOrbit) :: COrbit
    contains
    procedure :: Show   => ShowElem
    procedure :: Output => CallOutputOrbit
  end type COrbit

  contains

  type(COrbit) function GetDebriOrbit(n, date)
    integer, intent(in) :: n
    real(8), intent(in) :: date

    GetDebriOrbit = COrbit(DEBRIS(n)%orbit%CalcOrbit(date))
  end function GetDebriOrbit

  subroutine ShowElem(self)
    class(COrbit), intent(in) :: self
    print "(a, f11.5)", "epoch: ", self%epoch
    print "(a, f11.5)", "sma:   ", self%sma
    print "(a, f11.5)", "slr:   ", self%slr
    print "(a, f11.5)", "inc:   ", self%inc
    print "(a, f11.5)", "ran:   ", self%ran
    print "(a, f11.5)", "ecc:   ", self%ecc
    print "(a, f11.5)", "ap:    ", self%ap
    print "(a, f11.5)", "ma:    ", self%ma
    print "(a, f11.5)", "mm:    ", self%mm
    print *
  end subroutine ShowElem

  subroutine CallOutputOrbit(self, date, dt)
    class(COrbit), intent(in) :: self
    real(8),       intent(in) :: date, dt

    call OutputOrbit(self, date, dt)
  end subroutine CallOutputOrbit
end module Class_Orbit
