module ModUnitsSA
  use ModOperations
  use ModObject
  use ModIndividual
  use ModGroup

  implicit none

  private
  public :: TypeDataSA, SearchSA

  type :: TypeDataSA
    ! ------------------------------------------------------
    ! parameters
    ! ------------------------------------------------------
    integer              :: num_object_all, num_object_select, limit_step, convergence
    real(8)              :: cool, cost, temperature0
    logical              :: output

    ! ------------------------------------------------------
    ! status
    ! ------------------------------------------------------
    integer              :: generation
    real(8)              :: temperature
    type(TypeIndividual) :: individual_new, individual_best
  end type TypeDataSA

  contains

  subroutine SearchSA(dataSA)
    type(TypeDataSA), intent(inout) :: dataSA
    type(TypeIndividual)            :: individual_new
    real(8)                         :: cost(2), cost_best, probability

    individual_new  = Neightbour(dataSA)

    cost(1)         = dataSA%individual_new%cost_total
    cost(2)         = individual_new%cost_total
    cost_best       = dataSA%individual_best%cost_total

    cost            = cost * dataSA%cost

    probability     = ProbabilitySA(cost, dataSA%temperature)

    if (Random() < probability) dataSA%individual_new  = individual_new
    if (cost(2)  < cost_best)   dataSA%individual_best = individual_new
  end subroutine SearchSA

  ! ============================================================================
  ! private
  ! ============================================================================
  real(8) function ProbabilitySA(energy, temperature)
    real(8), intent(in) :: energy(2), temperature

    if (energy(2) <= energy(1)) then
      ProbabilitySA = 1.
    else
      ProbabilitySA = EXP((energy(1) - energy(2)) / temperature)
    end if
  end function ProbabilitySA

  ! ----------------------------------------------------------------------------

  type(TypeIndividual) function Neightbour(dataSA)
    type(TypeDataSA), intent(in) :: dataSA

    Neightbour       = CreateIndividual(dataSA%generation, dataSA%num_object_select, 710)
    Neightbour%order = Swap2Opt(dataSA%individual_new%order)
    call EvalIndividual(Neightbour)
  end function Neightbour

  function Swap2Opt(order)
    integer, allocatable :: Swap2Opt(:)
    integer, intent(in)  :: order(:)
    integer              :: pos(2)

    Swap2Opt                = order
    pos                     = Shuffle(SIZE(order) - 1, 2, 1) + iRandom(2) - 1
    Swap2Opt(pos(1):pos(2)) = Reverse(order(pos(1):pos(2)))
  end function Swap2Opt
end module ModUnitsSA

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

module ModMainSA
  use ModIndividual
  use ModGroup
  use ModOutput
  use ModUnitsSA

  implicit none

  private
  public :: OptimizeSA

  contains

  subroutine OptimizeSA(dataSA, individual_in)
    type(TypeDataSA),     intent(inout)        :: dataSA
    type(TypeIndividual), intent(in), optional :: individual_in
    type(TypeIndividual)                       :: individual_initial
    real(8)                                    :: cool
    real(8)                                    :: temperature0
    integer                                    :: num_object, limit_step, convergence, counter
    integer                                    :: step_output_window, step_output_file, step_output_check, i
    logical                                    :: output

    ! --------------------------------------------------------------------------
    ! read parameters
    ! --------------------------------------------------------------------------
    num_object   = dataSA%num_object_select
    cool         = dataSA%cool
    temperature0 = dataSA%temperature0
    limit_step   = dataSA%limit_step
    convergence  = dataSA%convergence
    output       = dataSA%output

    ! --------------------------------------------------------------------------
    ! set output interval
    ! --------------------------------------------------------------------------
    step_output_window = 1000
    step_output_file   = 100
    step_output_check  = 10

    ! --------------------------------------------------------------------------
    ! set initial individuals
    ! --------------------------------------------------------------------------
    if (PRESENT(individual_in)) then
      individual_initial = individual_in
    else
      individual_initial = CreateIndividual(0, num_object, 710)
    end if

    dataSA%individual_new  = individual_initial
    dataSA%individual_best = individual_initial

    ! --------------------------------------------------------------------------
    ! output initial data
    ! --------------------------------------------------------------------------
    if (output) then
      call OutputRoute(individual_initial, 0)
      call OutputCost(individual_initial, 0)
    end if

    counter = 0

    do i = 1, limit_step
      dataSA%generation = i

      ! ------------------------------------------------------------------------
      ! set temperature
      ! ------------------------------------------------------------------------
      if (MOD(i, step_output_check) == 0) then
        dataSA%temperature = temperature0 * cool ** (DBLE(i) / 100000.)
      else
        dataSA%temperature = 0.
      end if

      ! ------------------------------------------------------------------------
      ! main
      ! ------------------------------------------------------------------------
      call SearchSA(dataSA)

      ! ------------------------------------------------------------------------
      ! output
      ! ------------------------------------------------------------------------
      if (MOD(i, step_output_window) == 0) then
        print '(2(i10, f10.4, "|"), 2(f12.4, "|"))', i, dataSA%individual_new%cost_total, &
          dataSA%individual_best%generation, dataSA%individual_best%cost_total,           &
          dataSA%individual_best%cost_total - dataSA%individual_new%cost_total,           &
          dataSA%temperature
      end if

      if (output .and. MOD(i, step_output_file) == 0) call OutputCost(dataSA%individual_best, i)

      if (i == dataSA%individual_best%generation) then
        if (output .and. MOD(counter, step_output_file) == 0) call OutputRoute(dataSA%individual_best, i)
        counter = counter + 1
      end if

      ! ------------------------------------------------------------------------
      ! end determination
      ! ------------------------------------------------------------------------
      if (MOD(i, step_output_check) == 0 .and. dataSA%temperature < 1. &
          .and. i - dataSA%individual_best%generation > convergence) then
        if (output .and. MOD(counter, step_output_file) /= 1) call OutputRoute(dataSA%individual_best, i)
        exit
      end if
    end do

    print *, 'update', counter
  end subroutine OptimizeSA
end module ModMainSA
