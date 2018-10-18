module pareto_hyouka_mod
  use globals
  use base_func
  implicit none

  private
  public :: RANK_TYPE, ORDER_LIST_COPY, FIT_COPY, RANK_PARETO
  public :: pareto_hyouka

  logical :: is_show = .false.

  integer :: RANK_TYPE = 2
  integer :: ORDER_LIST_COPY(n, tar)
  real(8) :: FIT_COPY(n)
  integer :: RANK_PARETO(n)

contains
  subroutine pareto_hyouka(order_list, fit, elite_group, local_group, elite_fit, &
    local_fit, total_RCS, elite_RCS, local_RCS, c, elite_flight_time, total_flight_time)

    integer, intent(in) :: order_list(n, tar)
    real(8), intent(in) :: fit(n)
    real(8), intent(in) :: total_RCS(n)
    real(8), intent(in) :: total_flight_time(n)
    integer, intent(out) :: c
    integer, intent(out), allocatable :: elite_group(:, :)                      ! 上位5個体
    integer, intent(out), allocatable :: local_group(:, :)                      ! 上位5個体を取り除いた集団
    real(8), intent(out), allocatable :: elite_fit(:)                           ! 評価値の上位5個体
    real(8), intent(out), allocatable :: local_fit(:)                           ! 上位5個体以外の評価値
    real(8), intent(out), allocatable :: elite_RCS(:)
    real(8), intent(out), allocatable :: local_RCS(:)
    real(8), intent(out), allocatable :: elite_flight_time(:)
    integer, allocatable :: elite_num(:)
    integer, allocatable :: local_num(:)
    real(8) :: object(2, n)
    integer :: rank(n)
    logical :: mask(n)
    integer :: i

    object(1, :) = fit
    object(2, :) = 1d0 / total_RCS

    call rank_pareto_NS(object, rank)

    mask = rank == 1
    call reject_duplicate(order_list, mask)

    c = count(mask)
    elite_num = pack([(i, i = 1, n)], mask)
    local_num = pack([(i, i = 1, n)], .not. mask)

    elite_group = order_list(elite_num, :)
    elite_fit = fit(elite_num)
    elite_RCS = total_RCS(elite_num)
    elite_flight_time = total_flight_time(elite_num)

    local_group = order_list(local_num, :)
    local_fit = fit(local_num)
    local_RCS = total_RCS(local_num)

    ! for external process
    ORDER_LIST_COPY = order_list
    FIT_COPY = fit
    RANK_PARETO = rank

    if (is_show) then
      print "(" // Str(n) // "(i3))", Sort(RANK_PARETO)
    end if
  end subroutine pareto_hyouka

  subroutine reject_duplicate(order_list, mask)
    integer, intent(in) :: order_list(:, :)
    logical, intent(inout) :: mask(:)
    integer :: i, j

    do i = 1, n
      if (.not. mask(i)) cycle

      where([(i < j .and.                                            &
              all(order_list(i, :) == order_list(j, :)), j = 1, n)]) &
        mask = .false.
    end do
  end subroutine reject_duplicate

  subroutine rank_pareto_NS(object, rank_result)
    real(8), intent(in)  :: object(:, :)
    integer, intent(out) :: rank_result(:)
    integer, allocatable :: num_dominate(:)
    logical, allocatable :: non_dominated(:)
    integer              :: num_indiv, rank, i, j

    num_indiv = size(object)

    num_dominate = [(count([(all(object(:, i) >= object(:, j)) .and.         &
                             any(object(:, i) /= object(:, j)), j = 1, n)]), &
                     i = 1, n)]

    rank_result = 0

    rank = 1
    do while (any(num_dominate > 0))
      non_dominated = rank_result == 0 .and. num_dominate == 0

      where (non_dominated) rank_result = rank

      do i = 1, n
        if (non_dominated(i)) then
          where ([(all(object(:, i) <= object(:, j)) .and. &
                   any(object(:, i) /= object(:, j)), j = 1, n)])
            num_dominate = num_dominate - 1
          end where
        end if
      end do

      rank = rank + 1
    end do

    where (rank_result == 0 .and. num_dominate == 0) rank_result = rank
  end subroutine rank_pareto_NS
end module pareto_hyouka_mod
