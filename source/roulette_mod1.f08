module roulette_mod                                         ! 選択モジュール
  implicit none

contains
  subroutine roulette_par(local_group, probability, select_order, local_fit, select_fit, c)
    use globals
    use random_mod

    integer, intent(in) :: c
    real(8), intent(in) :: probability(:)
    integer, intent(in) :: local_group(:, :)
    real(8), intent(in) :: local_fit(:)
    integer, intent(out) :: select_order(nselect, tar)
    real(8), intent(out) :: select_fit(nselect)
    real(8), allocatable :: kakuritu(:)
    integer ns(nselect)
    integer min_num                                         ! デブリ番号一時格納用
    real(8) min_prob                                        ! ルーレット点数一時格納用
    integer i, j
    real(8) r

    !-----------------------------------c個体を選び出す
    ! probabilityは選択される確率で，値が小さいほど選ばれやすい
    ! probabilityにランダムな値rを乗じてkakurituとし，これが小さいものから順に選ばれる
    allocate(kakuritu(n - c))

    do i = 1, n - c
      call random_single(r)
      kakuritu(i) = probability(i) * r
    end do

    do i = 1, nselect
      min_num = 0
      min_prob = 1.0d0
      do j = 1, n - c
        if (kakuritu(j) < min_prob) then
          min_num = j
          min_prob = kakuritu(j)
        end if
      end do
      ns(i) = min_num
      kakuritu(min_num) = 1.0d0
    end do

    deallocate(kakuritu)

    !---------------------------------------------------------------------------
    select_order(1:nselect, 1:tar) = local_group(ns(1:nselect), 1:tar)          ! select_orderに格納
    select_fit(1:nselect) = local_fit(ns(1:nselect))

  end subroutine roulette_par
end module roulette_mod
