module orbit_optimizer
  use util
  use orbit_func
  use orbit_base
  use problem
  use soga
  use nsga2
  implicit none

  private
  public :: optimizer_main
  ! public :: Sgn, EulerAngle, cross_product, Unit, calc_angle, Rotate, mjd
  ! public :: PI, PI2, RADIANS, DEGREES

  type, extends(TProblem) :: TOP
    type(TOrbit) :: dep, arr, trf
    real(8) :: start, theta, offset

    contains

    generic :: initialize => initialize_op

    procedure :: initialize_op
    procedure :: call_func_so_con
    procedure :: call_func_mo
  end type TOP

  contains

  subroutine initialize_op(this)
    implicit none
    class(TOP), intent(inout) :: this
    real(8) :: start
    integer :: i

    start = mjd(2018, 0d0)
    this%start = start
    i = 10
    call this%dep%initialize(epc=start, inc=99d0, ran=0d0, ecc=0d0, ap=0d0, ma=0d0, mm=86400d0/calc_period(7178d0))
    call this%arr%initialize(epc=start, inc=99d0, ran=30d0, ecc=0d0, ap=0d0, ma=3d0 * i, mm=86400d0/calc_period(8078d0))
  end subroutine initialize_op

  subroutine call_func_so_con(this, variables, objective, feasible)
    implicit none
    class(TOP), intent(inout) :: this
    real(8), intent(in) :: variables(:)
    real(8), intent(out) :: objective
    logical, intent(out) :: feasible
    real(8) :: range, dd, dt, dv

    range = 6000d0
    dd = this%start + range * variables(1) * SEC_DAY
    dt = range * variables(2)

    call lambert(dd, dt, this%dep, this%arr, dv)

    objective = dv
    feasible = .true.
  end subroutine call_func_so_con

  subroutine call_func_mo(this, variables, objectives)
    implicit none
    class(TOP), intent(inout) :: this
    real(8), intent(in) :: variables(:)
    real(8), intent(out), allocatable :: objectives(:)
    real(8) :: range, dd, dt, dv

    range = 6000d0
    dd = this%start + range * variables(1) * SEC_DAY
    dt = range * variables(2)

    call lambert(dd, dt, this%dep, this%arr, dv)

    objectives = [dv]
  end subroutine call_func_mo

  subroutine optimizer_main
    implicit none
    type(TNSGA2) :: optimizer
    type(TOP) :: problem
    real(8), allocatable :: best(:)
    integer :: num_var, pop_size, num_step

    num_var = 2
    pop_size = 500
    num_step = 500

    call optimizer%initialize(nx=num_var, N=pop_size, m=1, selection="T", crossover="BLX", mutation="PM")

    call optimizer%set_fitness_type("VALUE")
    optimizer%elite_preservation = .true.
    optimizer%sharing = .true.
    optimizer%history_preservation = .false.

    call problem%initialize
    call optimizer%set_problem(problem)

    call optimizer%set_scaling_function(scaling)
    call optimizer%prepare_calculation
    call optimizer%run(num_step)
    print *, "A"
    call optimizer%best(best)

    call flight(problem%dep, problem%arr, problem%start, best(1), best(2))

    contains

    function scaling(variables) result(res)
      implicit none
      real(8), intent(in) :: variables(:)
      real(8), allocatable :: res(:)
      real(8) :: low(2), upp(2)

      low = [0d0, 0d0]
      upp = [6000d0, 6000d0]
      res = low + (upp - low) * variables
    end function scaling
  end subroutine optimizer_main

  subroutine flight(dep, arr, st, of, dt)
    implicit none
    type(TOrbit), intent(inout) :: dep, arr
    type(TOrbit) :: trf
    real(8), intent(in) :: st, of, dt
    ! real(8) :: axis(3)!, vel(3), pos(3)
    real(8) :: dd, dv!, theta, offset

    dd = st + of * SEC_DAY

    call lambert(dd, dt, dep, arr, dv, trf)

    print "(3(af15.8))", "best of=", of, " dt= ", dt, " dv=", dv

    call dep%save(dd, dt, "orbit1.csv", .true.)
    call arr%save(dd, dt, "orbit2.csv", .true.)
    call trf%save(dd, dt, "orbit3.csv", .true.)
  end subroutine flight
end module orbit_optimizer
