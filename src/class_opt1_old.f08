module class_opt1
  use interface
  use util
  use individual
  use problem
  use ga_unit
  use soga
  use nsga2
  use moead
  use orbit_func
  use orbit_base
  use orbit_debri
  use class_prob1
  implicit none

  private
  public :: TOpt1

  type, extends(TSOGA) :: TOpt1
    integer :: num_var1
    class(TCrossover), allocatable :: crossover1
    class(TMutation), allocatable :: mutation1

    contains

    procedure :: initialize_ex => initialize_opt1
    procedure :: init_indiv => init_indiv_opt1
    procedure :: reproduce
    procedure :: reproduce_m => reproduce
    procedure :: logger
  end type TOpt1

  contains


  ! ============================================================================
  ! Optimizer Class
  ! ============================================================================

  subroutine initialize_opt1(this, nx, crossover, mutation)
    implicit none
    class(TOpt1), intent(inout) :: this
    integer, intent(in) :: nx
    character(*), intent(in) :: crossover, mutation

    this%num_var1 = nx

    allocate(TCrossover::this%crossover1)
    allocate(TMutation::this%mutation1)

    call this%crossover1%initialize(crossover, 1d0, 0.5d0)
    call this%mutation1%initialize(mutation, 0.1d0, 20d0)
  end subroutine initialize_opt1

  subroutine init_indiv_opt1(this, indiv)
    implicit none
    class(TOpt1), intent(inout) :: this
    class(TIndiv), intent(inout) :: indiv

    select type(indiv)
    class is (TIndiv1)
      call indiv%initialize(this%num_var, this%num_var1)
    end select
  end subroutine init_indiv_opt1

  subroutine reproduce(this, index, children)
    implicit none
    class(TOpt1), intent(in) :: this
    integer, intent(in) :: index(:)
    class(TIndiv), intent(out), allocatable :: children(:)
    integer, allocatable :: parents_i(:, :), children_i(:, :)
    real(8), allocatable :: parents_d(:, :), children_d(:, :)
    integer :: p, i

    allocate(children(2), source=this%prototype)

    ! *** ivar ***
    parents_i = reshape([(this%population(index(i))%indiv%ivariables, i = 1, 2)], &
                        [this%num_var, 2])

    call this%crossover%call(parents_i, children_i)
    ! allocate(children(size(children_i, dim=2)), source=this%prototype)

    do i = 1, size(children)
      ! call children(i)%set_variables(children_i(:, i))
      ! call this%mutation%call(children(i)%ivariables)
      call children(i)%set_variables(this%population(index(i))%indiv%ivariables)
    end do

    ! *** dvar ***
    parents_d = reshape([(this%population(index(i))%indiv%dvariables, i = 1, 2)], &
                        [this%num_var1, 2])

    call this%crossover1%call(parents_d, children_d)
    ! allocate(children(size(children_d, dim=2)), source=this%prototype)

    do i = 1, size(children)
      p = min(i, size(children_d, dim=2))
      call children(i)%set_variables(children_d(:, p))
      call children(i)%clamp_variables(lower=0d0, upper=1d0)
      call this%mutation1%call(children(i)%dvariables)
      call children(i)%clamp_variables(lower=0d0, upper=1d0)
    end do
  end subroutine reproduce

  subroutine reproduce_s(this, index, child)
    implicit none
    class(TOpt1), intent(in) :: this
    integer, intent(in) :: index(:)
    class(TIndiv), intent(out), allocatable :: child
    class(TIndiv), allocatable :: children(:)

    call this%reproduce_m(index, children)
    child = children(1)
  end subroutine reproduce_s

  subroutine logger(this, n)
    implicit none
    class(TOpt1), intent(in) :: this
    integer, intent(in) :: n
    real(8) :: best

    if (n == 0) then
      print "(2a)", "Start: ", this%problem_type
    else if (n == -1) then
      print "(2a/)", "End: ", this%problem_type
    else if (mod(n, 1) == 0) then
      call this%best_obj(best)
      print "(a i0 a es15.8)", "Progress: ", n, " steps finished best=", best
    end if
  end subroutine logger
end module class_opt1
