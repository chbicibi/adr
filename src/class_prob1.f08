module class_prob1
  ! 再帰型GAの内側問題クラス
  use interface
  use util
  use individual
  use problem
  use orbit_func
  use orbit_base
  use orbit_debri
  implicit none

  private
  public :: TIndiv1, TProb1

  type, extends(TIndivS) :: TIndiv1
    contains

    ! procedure :: initialize_n => initialize_indiv1
    ! procedure :: is_same => is_same_indiv1
  end type TIndiv1

  type, extends(TProblem) :: TProb1
    type(TOrbit) :: dep, arr
    real(8) :: start

    contains

    generic :: initialize => initialize_prob

    procedure :: initialize_prob => initialize_prob1
    procedure :: call_indiv => call_prob1
  end type TProb1

  contains


  ! ============================================================================
  ! Problem Class
  ! ============================================================================

  subroutine initialize_prob1(this, debris, id, ia)
    implicit none
    class(TProb1), intent(inout) :: this
    type(TDebri), intent(in) :: debris(0:)
    integer, intent(in) :: id, ia

    this%start = mjd(2018, 0d0)
    this%dep = debris(id)%orbit
    this%arr = debris(ia)%orbit
  end subroutine initialize_prob1

  subroutine call_prob1(this, indiv)
    implicit none
    class(TProb1), intent(inout) :: this
    class(TIndiv), intent(inout) :: indiv
    real(8) :: var(2)
    real(8) :: start, duration, dv

    var = this%scaling_func(indiv%dvariables)
    start = this%start + var(1)
    duration = var(2)

    call lambert(start, duration, this%dep, this%arr, dv)

    allocate(indiv%objectives(1), source=dv)
    indiv%feasible = .true.
  end subroutine call_prob1


  ! ============================================================================
  ! Individual Class
  ! ============================================================================

  ! subroutine initialize_indiv1(this, nvar)
  !   implicit none
  !   class(TIndiv1), intent(inout) :: this
  !   integer, intent(in) :: nvar

  !   this%dvariables = random_array(nvar)
  !   ! this%ivariables = shuffle(NUM_DEBRIS, nvar)
  !   ! this%ivariables = integers(nvar)
  !   ! this%ivariables = [1995, 1681]
  ! end subroutine initialize_indiv1

  ! logical function is_same_indiv1(this, other) result(res)
  !   implicit none
  !   class(TIndiv1), intent(in) :: this
  !   class(TIndiv), intent(in) :: other

  !   res = all(abs(this%dvariables - other%dvariables) < 1d-4)
  ! end function is_same_indiv1

  ! subroutine print_header_indiv1(this, unit)
  !   implicit none
  !   class(TIndiv1), intent(in) :: this
  !   integer, intent(in) :: unit
  !   type(string), allocatable :: s(:)
  !   integer :: l(3)

  !   l(1) = size(this%objectives)
  !   l(2) = size(this%ivariables)
  !   l(3) = size(this%dvariables)
  !   allocate(s(sum(l)), source=[serstr("obj", integers(l(1))), &
  !                               serstr("ivar", integers(l(2))), &
  !                               serstr("dvar", integers(l(3)))])
  !   write(unit, "(a)", advance='no') join(s, ",")
  ! end subroutine print_header_indiv1

  ! subroutine print_indiv1(this, unit, proc)
  !   implicit none
  !   class(TIndiv1), intent(in) :: this
  !   integer, intent(in) :: unit
  !   procedure(func_1d_1d) :: proc
  !   character(:), allocatable :: format

  !   format = "(" // nformat(size(this%objectives), "es15.8','") &
  !                // nformat(size(this%ivariables), "i0','")     &
  !                // nformat(size(this%dvariables), "es15.8", ",") // ")"
  !   write(unit, format, advance='no') this%objectives, this%ivariables, proc(this%dvariables)
  ! end subroutine print_indiv1
end module class_prob1
