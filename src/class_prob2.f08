module class_prob2
  ! 再帰型GAの外側問題クラス
  use interface
  use util
  use individual
  use problem
  use orbit_func
  use orbit_base
  use orbit_debri
  use class_prob1
  use class_opt1
  implicit none

  private
  public :: TIndiv2, TProb2

  type, extends(TIndivI) :: TIndiv2
    contains

    procedure :: initialize_n => initialize_indiv1
    procedure :: print_u => print_indiv1
  end type TIndiv2

  type, extends(TProblem) :: TProb2
    type(TOpt1), allocatable :: opt_table(:, :)
    type(TDebri), allocatable :: debris(:)
    real(8), allocatable :: table(:, :)

    contains

    generic :: initialize => initialize_prob

    procedure :: initialize_prob => initialize_prob1
    procedure :: call_indiv => call_prob1
    procedure :: set_debris
    procedure :: optimize
    procedure :: save_table, load_table
  end type TProb2

  integer :: NUM_DEBRIS

  contains


  ! ============================================================================
  ! Problem Class
  ! ============================================================================

  subroutine initialize_prob1(this, n, order, opt, scaling)
    implicit none
    class(TProb2), intent(inout) :: this
    integer, intent(in) :: n, order(0:)
    type(TOpt1), intent(in) :: opt
    procedure(func_1d_1d) :: scaling
    type(TProb1) :: inner_problem
    integer :: i, j

    NUM_DEBRIS = n
    call this%set_debris(order)

    allocate(this%opt_table(0:n, 0:n))
    allocate(this%table(0:n, 0:n))

    do j = 0, n
      do i = 0, n
        if (i == j) cycle
        this%opt_table(i, j) = opt
        call inner_problem%initialize(this%debris, i, j)
        call this%opt_table(i, j)%set_problem(inner_problem)
        call this%opt_table(i, j)%set_scaling_function(scaling)
      end do
    end do
  end subroutine initialize_prob1

  subroutine call_prob1(this, indiv)
    implicit none
    class(TProb2), intent(inout) :: this
    class(TIndiv), intent(inout) :: indiv
    real(8) :: dv_total, rcs
    integer :: s, p, n, i

    dv_total = 0d0
    rcs = 0d0
    s = size(indiv%ivariables)

    do i = 1, s
      if (i == 1) then
        p = 0
      else
        p = indiv%ivariables(i-1)
      end if
      n = indiv%ivariables(i)

      dv_total = dv_total + this%table(p, n)
      rcs = rcs + this%debris(n)%rcs
    end do

    allocate(indiv%objectives(2), source=[dv_total, 1/rcs])
    indiv%feasible = .true.
  end subroutine call_prob1

  subroutine set_debris(this, order)
    implicit none
    class(TProb2), intent(inout) :: this
    integer, intent(in) :: order(0:)

    allocate(this%debris(0:ubound(order, dim=1)), source=DEBRIS(order))
  end subroutine set_debris

  subroutine optimize(this, step)
    implicit none
    class(TProb2), intent(inout) :: this
    integer, intent(in) :: step
    real(8) :: best
    integer :: i, j, n

    n = ubound(this%table, dim=1)

    do j = 0, n
      if (j > 0 .and. mod(j, 10) == 0) print *, "optimize", j
      do i = 0, n
        if (i == j) cycle
        call this%opt_table(i, j)%prepare_calculation
        call this%opt_table(i, j)%run(step)
      end do
    end do

    do j = 0, n
      do i = 0, n
        if (i == j) cycle
        call this%opt_table(i, j)%best_obj(best)
        ! print *, "best", i, j, best
        this%table(i, j) = best
      end do
    end do
  end subroutine optimize

  subroutine save_table(this, filename)
    implicit none
    class(TProb2), intent(inout) :: this
    character(*), intent(in) :: filename
    real(8), allocatable :: dvar(:)
    integer :: unit, i, j, n

    n = ubound(this%table, dim=1)
    open(newunit=unit, file=filename)
      write(unit, "(a)") "i,j,dv"
      do j = 0, n
        do i = 0, n
          if (i == j) cycle
          call this%opt_table(i, j)%best(dvar)
          write(unit, "(2(i0',')3(es15.8','))") i, j, this%table(i, j), dvar
        end do
      end do
    close(unit)
  end subroutine save_table

  subroutine load_table(this, filename)
    implicit none
    class(TProb2), intent(inout) :: this
    character(*), intent(in) :: filename
    integer :: unit, i, j, n, idum

    n = ubound(this%table, dim=1)
    open(newunit=unit, file=filename)
      read(unit, *)
      do j = 0, n
        do i = 0, n
          if (i == j) cycle
          read(unit, *) idum, idum, this%table(i, j)
        end do
      end do
    close(unit)
  end subroutine load_table


  ! ============================================================================
  ! Individual Class
  ! ============================================================================

  subroutine initialize_indiv1(this, nvar)
    implicit none
    class(TIndiv2), intent(inout) :: this
    integer, intent(in) :: nvar

    ! this%dvariables = random_array(nvar1)
    this%ivariables = shuffle(NUM_DEBRIS, nvar)
    ! this%ivariables = integers(nvar)
    ! this%ivariables = [1995, 1681]
  end subroutine initialize_indiv1

  logical function is_same_indiv1(this, other) result(res)
    implicit none
    class(TIndiv2), intent(in) :: this
    class(TIndiv), intent(in) :: other

    ! res = all(this%ivariables == other%ivariables) .and. &
    !       all(abs(this%dvariables - other%dvariables) < 1d-4)
    res = all(this%ivariables == other%ivariables)
    ! res = all(this%dvariables == other%dvariables)
  end function is_same_indiv1

  subroutine print_header_indiv1(this, unit)
    implicit none
    class(TIndiv2), intent(in) :: this
    integer, intent(in) :: unit
    type(string), allocatable :: s(:)
    integer :: l(3)

    l(1) = size(this%objectives)
    l(2) = size(this%ivariables)
    l(3) = size(this%dvariables)
    allocate(s(sum(l)), source=[serstr("obj", integers(l(1))),  &
                                serstr("ivar", integers(l(2))), &
                                serstr("dvar", integers(l(3)))])
    write(unit, "(a)", advance='no') join(s, ",")
  end subroutine print_header_indiv1

  subroutine print_indiv1(this, unit)!, proc)
    implicit none
    class(TIndiv2), intent(in) :: this
    integer, intent(in) :: unit
    ! procedure(func_1d_1d) :: proc
    character(:), allocatable :: format
    real(8) :: obj(2)

    obj = this%objectives
    obj(2) = 1 / obj(2)

    format = "(" // nformat(size(this%objectives), "es15.8','") &
                 // nformat(size(this%ivariables), "i0", ",") // ")"!     &
                 ! // nformat(size(this%dvariables), "es15.8", ",") // ")"
    write(unit, format, advance='no') obj, this%ivariables!, proc(this%dvariables)
  end subroutine print_indiv1
end module class_prob2
