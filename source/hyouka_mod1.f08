module hyouka_mod
  implicit none
  !-----------------------------------各種パラメータ
  ! real(8), parameter :: g = 8.1717302d12                    ! 重力パラメータ
  real(8), parameter :: GE = 3.986d5                        ! 地心重力定数
  real(8), parameter :: GS = 1.327d11                       ! 太陽重力定数
  real(8), parameter :: AU = 149597870.7d0                  ! 天文単位
  real(8), parameter :: PI = 3.141592653d0                  ! 円周率

  real(8), parameter :: r_time = 100.0d0                    ! 時間スケール係数(地球->太陽)
  real(8) r_length                                          ! 距離スケール係数(地球->太陽)
  real(8) r_vel                                             ! 速度スケール係数(地球->太陽)

contains
  subroutine read_debridata
    use globals

    integer i
    integer geny(db)                                        ! 元期(年)
    real(8) gend(db)                                        ! 元期(日)
    character(len = 14) gens(db)
    character(len = 1) dum(db)
    character(len = 6) dum1(db), dum2(db)
    character(len = 16) dum3(db), kari
    integer line(db), catalogs(db)
    integer yearg, monthg, dayg
    real y, m
    integer mjd

    !-----------------------------------係数
    r_length = (GS / GE * r_time ** 2.0) ** (1.0 / 3.0)
    r_vel = (GS / GE / r_time) ** (1.0 / 3.0)

    !-----------------------------------ファイル読み込み
    open(10, file = 'debri_elements.txt')
    do i = 1, db
      read(10, *)
      read(10, *) dum(i), dum1(i), dum2(i), gens(i)
      read(10, *) line(i), catalogs(i), incl(i), omega(i), ecc(i), argper(i), ma(i), dum3(i)
    end do
    close(10)

    open(11, file = 'RCS_list.txt')
    do i = 1, db
      read(11, *) RCS(i)
    end do
    close(11)

    !-----------------------------------文字列の変換
    do i = 1, db
      kari = dum3(i)
      read(kari(1:11), *) mm(i)
      read(kari(12:16), "(a)") syukai_num(i) ! 2017/10/12
      read(gens(i)(1:2), *) geny(i)
      read(gens(i)(3:14), *) gend(i)
    end do

    !-----------------------------------離心率
    ecc(1:db) = ecc(1:db) / 10000000d0

    !-----------------------------------起動長半径計算(太陽中心スケール)
    smma(1:db) = r_length * (GE * (86400.0 / mm(1:db) / (2.0 * PI)) ** 2) ** (1.0 / 3.0)

    !-----------------------------------修正ユリウス暦計算/元期修正
    ! yearg = 2015
    ! monthg = 1
    ! dayg = 1

    ! ! フリーゲルの公式
    ! if (monthg >= 3) then
    !   y = yearg
    !   m = monthg
    ! else
    !   y = yearg - 1
    !   m = monthg + 12
    ! end if
    ! mjd = floor(365.25 * y) + floor(y / 400.0) - floor(y / 100.0) + floor(30.59 * (m - 2)) + dayg - 678912

    ! gen(1:db) = int((geny(1:db) - 15) * 365 + gend(1:db)) * r_time + mjd      ! 太陽中心スケールに変換
    gen(1:db) = 57023                                       ! 元期固定(テスト用)

  end subroutine read_debridata

  subroutine hyouka(order_list, fit, total_RCS, total_flight_time, k)
    use globals
    use random_mod

    integer, intent(in) :: order_list(n, tar)
    real(8), intent(out) :: fit(n)
    real(8), intent(out) :: total_RCS(n)
    real(8), intent(out) :: total_flight_time(n)

    integer i, j, k, num
    real(8) a, b
    real(8) aa, ab, ac, ad, ae
    real(8) delta_V_d
    real(8) delta_V_a
    real(8) total_delta_V

    !-----------------------------------遷移時間
    real(8) flight_time(n, tar)                             ! 各遷移時間
    ! real(8) total_flight_time(n)                          ! 合計遷移時間
    real(8) tot_time

    !-----------------------------------各個体の軌道要素
    integer select_gen(n, tar)
    real(8) select_incl(n, tar)
    real(8) select_omega(n, tar)
    real(8) select_ecc(n, tar)
    real(8) select_argper(n, tar)
    real(8) select_ma(n, tar)
    real(8) select_mm(n, tar)
    real(8) select_smma(n, tar)
    real(8) sele_syukai_num(n, tar)
    real(8) select_RCS(n, tar)

    !-----------------------------------出発軌道要素
    integer, parameter :: gen1 = 57023                      ! 2015/1/1 00:00:00
    real(8), parameter :: incl1 = 30.0d0
    real(8), parameter :: omega1 = 0.0d0
    real(8), parameter :: ecc1 = 0.0d0
    real(8), parameter :: argper1 = 0.0d0
    real(8), parameter :: ma1 = 0.0d0
    real(8), parameter :: long_a = 6578.0d0

    !-----------------------------------出発時刻
    integer, parameter :: year = 2015
    integer, parameter :: month = 1
    integer, parameter :: day = 1
    integer, parameter :: hour = 0
    integer, parameter :: minute = 0
    integer, parameter :: second = 0

    !---------------------------------------------------------------------------
    integer :: flag = 0                                     ! ADR_GA__trj.outファイル出力の有無
    character(80) str                                       ! ARDR_GA.outファイル読み込み用
    real(8) TOF                                             ! 伝播期間
    integer :: ref = 1                                      ! 軌道解参照[11]

    !-----------------------------------それぞれselectに代入
    do i = 1, n
      do j = 1, tar
        num = order_list(i, j)
        select_gen(i, j) = gen(num)
        select_incl(i, j) = incl(num)
        select_omega(i, j) = omega(num)
        select_ecc(i, j) = ecc(num)
        select_argper(i, j) = argper(num)
        select_ma(i, j) = ma(num)
        select_mm(i, j) = mm(num)
        select_smma(i, j) = smma(num)
        sele_syukai_num(i, j) = syukai_num(num)
        select_RCS(i, j) = RCS(num)
      end do
    end do

    !---------------------------------------------------------------------------
    total_RCS(1:n) = 0

    do i = 1, n
      aa = 0
      ab = 0
      ac = 0
      ad = 0
      ae = 0
      total_delta_V = 0
      tot_time = 0

      do j = 1, tar
        !-------------------------------遷移時間設定
        flight_time(i, j) = 1.0 / 24.0                      ! 遷移時間固定(テスト用)
        TOF = flight_time(i, j) * r_time                    ! 太陽中心スケールに変更

        !-------------------------------ELEMENTS_JPLに出力
        open(20, file = '../dat/ELEMENTS_JPL.in', status = 'old')
        read(20, *)
        read(20, *)
        write(20, '(a, i25, f11.7, f11.8, f10.5, f10.5, f10.5, f12.6, i18, $)') &
       & 'debri', select_gen(i, j), select_smma(i, j) / AU, select_ecc(i, j), &
       & select_incl(i, j), select_argper(i, j), select_omega(i, j), select_ma(i, j), ref
        close(20)

        !-------------------------------ADR_GA_Para.inに出力
        open(30, file = '../dat/ADR_GA_Para.in', status = 'old')
        read(30, *)
        read(30, *)
        if (j == 1) then
          write(30, '(i4, a, i2, a, i2, i3.2, a, i2.2, a, i2.2, &
         & f12.3, i12, f18.3, f12.8, f12.5, f12.5, f12.5, f12.5, i5, $)') &
         & year, '/', month, '/', day, hour, ':', minute, ':', second, &
         & TOF, gen1, long_a * r_length, ecc1, omega1, incl1, argper1, ma1, flag
        else
          write(30, '(i4, a, i2, a, i2, i3.2, a, i2.2, a, i2.2, &
         & f12.3, i12, f18.3, f12.8, f12.5, f12.5, f12.5, f12.5, i5, $)') &
         & year, '/', month, '/', day, hour, ':', minute, ':', second, &
         & TOF, select_gen(i, j - 1), select_smma(i, j - 1), select_ecc(i, j - 1), select_omega(i, j - 1), &
         & select_incl(i, j - 1), select_argper(i, j - 1), select_ma(i, j - 1), flag
        end if
        close(30)

        !-------------------------------軌道計算
        call system('ADR_GA.exe > NUL')

        !-------------------------------------------------------------------------
        open(40, file = '../dat/ARDR_GA.out')
        read(40, *)
        read(40, *)
        read(40, '(a)') str
        close(40)
        str(1:index(str, ':') + 5) = ''                     ! ファイル読み込み用処理
        read(str, *) TOF, delta_V_d, delta_V_a

        !-------------------------------------------------------------------------
        delta_V_d = delta_V_d / r_vel                       ! 地球中心スケールに変更
        delta_V_a = delta_V_a / r_vel

        !-------------------------------合計ΔV算出
        total_delta_V = total_delta_V + delta_V_d + delta_V_a

       !  write(*, '(2(i3, a), i1, i5, a8, f10.5, 2a5, f12.5, a3)') k, ': ', i, '-', j, order_list(i, j), &
       ! & 'delV=', delta_V_d + delta_V_a, 'km/s', &
       ! & 'a=', select_smma(i, j) / r_length - 6378.0, 'km'

        !-------------------------------不良個体処理(仮)
        if (delta_V_d + delta_V_a < 1.0d-5) then
          write(*, *) 'error: ', order_list(i, j)
          total_delta_V = 99.0
          ! exit
        end if

        !-------------------------------RCS算出
        aa = aa + select_RCS(i, j)

        !-------------------------------遷移時間算出
        ! if (j == 1) then
        !   ab = (ma1 - ecc1 * sin(ma1))
        ! else
        !   ab = (select_ma(i, j - 1) - select_ecc(i, j - 1) * sin(select_ma(i, j - 1)))
        ! end if

        ! ac = (select_ma(i, j) - select_ecc(i, j) * sin(select_ma(i, j)))
        ! ad = 2 * sele_syukai_num(i, j) / 50 * PI
        ! a = GE / (select_mm(i, j) ** 2)
        ! b = (a / g) ** (1.0 / 2.0)
        ! ae = b * (ad + ab - ac)
        ! flight_time(i, j) = ae                              ! 各遷移時間
        tot_time = tot_time + flight_time(i, j)             ! 合計遷移時間

        if (j < tar) then
          ! 時刻を加算する処理を書く
        end if
      end do

      fit(i) = total_delta_V
      total_RCS(i) = aa
      total_flight_time(i) = tot_time

      write(*, '(2(i3, a), a9, f10.5, a)') k, ': ', i, ' total', 'delV=', total_delta_V, ' km/s'
    end do

  end subroutine hyouka
end module hyouka_mod
