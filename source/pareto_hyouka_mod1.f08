module pareto_hyouka_mod
  implicit none

  integer :: RANK_TYPE = 1

contains
  subroutine pareto_hyouka(order_list, fit, elite_group, local_group, elite_fit, &
    local_fit, total_RCS, elite_RCS, local_RCS, c, elite_flight_time, total_flight_time)
    use globals

    real(8), intent(in) :: total_RCS(n)
    integer, intent(in) :: order_list(n, tar)
    real(8), intent(in) :: fit(n)
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
    integer :: rank(n)
    integer :: count_num(n)
    integer :: ni, a
    integer :: i, j

    !-----------------------------------パレート解を抽出する
    c = 0
    do i = 1, n
      ni = 1

      do j = 1, n
        if (fit(i) > fit(j) .and. total_RCS(i) < total_RCS(j)) then
          ni = ni + 1
        end if
      end do

      rank(i) = ni
      if (rank(i) == 1) then
        c = c + 1
        count_num(c) = i
      end if

      if (c == 20) then                                     ! 解の数を制限(仮)
        exit
      end if
    end do

    !-----------------------------------パレート解の番号をコピー
    elite_num = count_num(1:c)
    local_num = pack([(i, i = 1, n)], [(all(elite_num /= i), i = 1, n)])

    a = 0
    do i = 1, n
      if (all(elite_num(1:c) /= i)) then
        a = a + 1
        local_num(a) = i
      end if
    end do

    !-----------------------------------パレート保存戦略
    elite_group = order_list(elite_num, :)
    elite_fit = fit(elite_num)
    elite_RCS = total_RCS(elite_num)
    elite_flight_time = total_flight_time(elite_num)

    !-----------------------------------その他
    local_group = order_list(local_num, :)
    local_fit = fit(local_num)
    local_RCS = total_RCS(local_num)
  end subroutine pareto_hyouka
end module pareto_hyouka_mod
