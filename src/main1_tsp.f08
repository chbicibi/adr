module class_main1
  use util
  use problem
  use individual
  use sogai
  implicit none

  private
  public :: TProb1, TIndiv1, TOpt1

  type, extends(TProblem) :: TProb1
    integer :: num_city
    real(8), allocatable :: loc(:, :)
    real(8), allocatable :: table(:, :)

    contains

    generic :: initialize => initialize_tsp
    generic :: call => call_indiv1

    procedure :: initialize_tsp
    procedure :: call_indiv1
    procedure :: print_route
  end type TProb1

  type, extends(TIndivSI) :: TIndiv1
    contains

    procedure :: evaluate => eval_prob1
  end type TIndiv1

  type, extends(TSOGAI) :: TOpt1
    contains
  end type TOpt1

  integer, parameter :: ROUND = 1

  contains


  ! ============================================================================
  ! Problem Class
  ! ============================================================================

  subroutine initialize_tsp(this, n)
    implicit none
    class(TProb1), intent(inout) :: this
    integer, intent(in) :: n
    integer :: i, j

    this%num_city = n
    this%loc = reshape(random_array(2 * n), [2, n])
    allocate(this%table(n, n))

    do j = 1, n
      this%table(j, j) = huge(0d0)
      do i = j + 1, n
        this%table(i, j) = distance(this%loc(:, i), this%loc(:, j))
        this%table(j, i) = this%table(i, j)
      end do
    end do
  end subroutine initialize_tsp

  subroutine call_indiv1(this, indiv)
    implicit none
    class(TProb1), intent(inout) :: this
    class(TIndiv1), intent(inout) :: indiv
    real(8) :: dist
    integer :: s, i

    dist = 0d0
    s = size(indiv%ivariables)

    do i = 1, s - 1 + ROUND
      dist = dist + this%table(indiv%ivariables(i), indiv%ivariables(mod(i, s) + 1))
    end do

    allocate(indiv%objectives(1), source=dist)
    indiv%feasible = .true.
  end subroutine call_indiv1

  subroutine print_route(this, variables, file)
    implicit none
    class(TProb1), intent(inout) :: this
    integer, intent(in) :: variables(:)
    character(*), intent(in) :: file
    integer :: unit, s, p, i

    s = size(variables)
    open(newunit=unit, file=file)
    write(unit, "(a)") "i,x,y"
    do i = 1, s + ROUND
      p = mod(i - 1, s) + 1
      write(unit, "(i0,2(','g0))") variables(p), this%loc(:, variables(p))
    end do
    close(unit)
  end subroutine print_route


  ! ============================================================================
  ! Individual Class
  ! ============================================================================

  subroutine eval_prob1(this, problem)
    implicit none
    class(TIndiv1), intent(inout) :: this
    class(TProblem), intent(inout) :: problem

    if (this%evaluated) return

    select type(problem)
    class is (TProb1)
      call problem%call(this)
    end select

    this%evaluated = .true.
  end subroutine eval_prob1


  ! ============================================================================
  ! Optimizer Class
  ! ============================================================================
end module class_main1

module mod_main1
  use util
  use class_main1
  implicit none

  contains

  subroutine main1
    implicit none
    type(TProb1) :: tsp
    type(TIndiv1) :: indiv
    type(TOpt1) :: optimizer
    integer :: nx, N, steps, start_time
    character(32) :: ts, tc, tm
    integer, allocatable :: best(:)

    call set_randomseed

    nx = 30
    N = 100
    ts = "ROULETTE"
    tc = "ORDER"
    tm = "SWAP"
    steps = 1000

    call tsp%initialize(nx)
    call optimizer%initialize(nx=nx, N=N, selection=ts, crossover=tc, mutation=tm)
    call optimizer%set_problem(tsp)
    call optimizer%set_prototype(indiv)
    call optimizer%prepare_calculation
    optimizer%elite_preservation = .true.
    optimizer%dup_rejection = .true.

    ! call optimizer%best(best)
    ! call tsp%print_route(best, "tsp_init.csv")

    call system_clock(start_time)
    call optimizer%run(steps)
    call elapsed_time(start_time)

    call optimizer%best(best)
    call tsp%print_route(best, "tsp_best.csv")

    call optimizer%save_result("tsp_result.csv")
    call optimizer%save_history("tsp_history.csv", elite="only")
  end subroutine main1
end module mod_main1

program main
  use mod_main1
  implicit none

  call main1
end program main
