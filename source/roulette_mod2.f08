module roulette_mod                                         ! 選択モジュール
  use globals
  use random_mod
  use pareto_hyouka_mod

  implicit none

  private
  public :: roulette_par
  public :: select_wait

  real(8) :: select_wait = 0.1

contains
  subroutine roulette_par(local_group, probability, select_order, local_fit, select_fit, c)
    integer, intent(in) :: c
    real(8), intent(in) :: probability(:)
    integer, intent(in) :: local_group(:, :)
    real(8), intent(in) :: local_fit(:)
    integer, intent(out) :: select_order(nselect, tar)
    real(8), intent(out) :: select_fit(nselect)
    real(8) :: prob(n)
    logical :: rest(n)
    integer :: selected, i

    if (RANK_TYPE /= 2) then
      print *, "error: 'pareto_hyouka_mod2.f08' is not included"
      stop
    end if

    prob = select_wait * (1.0d0 - select_wait) ** (1 - RANK_PARETO)
    rest = .true.

    ! call reject_duplicate(ORDER_LIST_COPY, rest)

    do i = 1, nselect
      call select_1(prob, rest, selected)

      select_order(i, :) = ORDER_LIST_COPY(selected, :)
      select_fit(i) = FIT_COPY(selected)
    end do
  end subroutine roulette_par

  subroutine reject_duplicate(order_list, rest)
    integer, intent(in) :: order_list(:, :)
    logical, intent(inout) :: rest(:)
    integer :: i, j

    do i = 1, n
      if (.not. rest(i)) cycle

      where([(i < j .and.                                            &
              all(order_list(i, :) == order_list(j, :)), j = 1, n)]) &
        rest = .false.
    end do
  end subroutine reject_duplicate

  subroutine select_1(prob, rest, selected)
    real(8), intent(in) :: prob(:)
    logical, intent(inout) :: rest(:)
    integer, intent(out) :: selected
    real(8) :: r
    integer :: i

    call random_single(r)

    r = r * sum(prob, mask=rest)

    do i = 1, n
      if (.not. rest(i)) cycle

      r = r - prob(i)
      if (r < 0) then
        rest(i) = .false.
        selected = i
        return
      end if
    end do

    print *, "error: incorrect roulette algorithm (in subroutine 'select_1')"
    stop
  end subroutine select_1
end module roulette_mod
