module globals                                              ! グローバル変数mod
  integer, parameter :: n = 40                              ! n:個体数(20) <- 40
  integer, parameter :: db = 100                            ! デブリ数(100個)
  integer, parameter :: tar = 5                             ! 回収目標数(5個)
  integer, parameter :: nselect = 10                        ! 選択数(偶数)
  real(8), parameter :: crate = 0.9d0                       ! 交叉確率 <- 1.0
  real(8), parameter :: mttn_rate = 0.01d0                  ! 変異確率 <- 0.0
  integer, parameter :: T = 100                             ! 世代数

  !-------------------------------------デブリデータ
  integer gen(db)                                           ! 元期(整数, 太陽中心スケール)
  real(8)  incl(db)                                         ! 傾斜角[10]
  real(8)  omega(db)                                        ! 昇交点経度[10]
  real(8)  ecc(db)                                          ! 離心率[11]
  real(8)  argper(db)                                       ! 近日点引数[10]
  real(8)  ma(db)                                           ! 平均近点離角[12]
  real(8)  mm(db)                                           ! 平均運動
  real(8)  smma(db)                                         ! 起動長半径(太陽中心スケール)
  integer  syukai_num(db)                                   ! 周回数
  real(8)  RCS(db)                                          ! レーダー反射断面積
end module globals

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

module random_mod                                           ! 乱数を作るモジュール
  implicit none

contains
  subroutine set_seed                                       ! 初回seed設定
    integer, allocatable :: seed(:)
    integer nrand, clock

    call random_seed(size = nrand)
    allocate(seed(nrand))
    if (.false.) then
      seed(1:nrand) = 123456789                             ! seed固定(テスト用)
    else
      call system_clock(count = clock)
      seed(1:nrand) = clock
    end if

    call random_seed(put = seed)
    deallocate(seed)
  end subroutine set_seed

  subroutine random_single(r)
    real(8), intent(out) :: r

    call random_number(r)

  end subroutine random_single
end module random_mod

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

module debri_order_mod
  implicit none

contains
  subroutine debri_order_table(order_list)
    use globals
    use random_mod

    integer, intent(out) :: order_list(n, tar)
    integer i, j, k, num
    real(8) r

    do i = 1, n
      !---------------------------------順列の生成
      do j = 1, tar
        call random_single(r)
        num = int(real(db) * r) + 1

        if (j > 1) then
          do
            if (any(order_list(i, 1:j-1) == num)) then
              call random_single(r)
              num = int(real(db) * r) + 1
            else
              exit
            end if
          end do
        end if

        order_list(i, j) = num
      end do
    end do

  end subroutine debri_order_table
end module debri_order_mod

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

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
      read(kari(12:16), *) syukai_num(i)
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

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

module pareto_hyouka_mod
  implicit none

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
    integer i, j
    integer ni, a
    integer rank(n)
    integer, allocatable :: count_num(:)
    integer, allocatable :: elite_num(:)
    integer, allocatable :: local_num(:)

    allocate(count_num(n))
    count_num(1:n) = 0

    !-----------------------------------パレート解を抽出する
    c = 0
    do i = 1, n
      exit
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

    !-----------------------------------割り付け
    allocate(elite_num(c))
    allocate(elite_group(c, tar))
    allocate(elite_fit(c))
    allocate(elite_RCS(c))
    allocate(local_num(n - c))
    allocate(local_group(n - c, tar))
    allocate(local_fit(n - c))
    allocate(local_RCS(n - c))
    allocate(elite_flight_time(c))

    !-----------------------------------パレート解の番号をコピー
    elite_num(1:c) = count_num(1:c)
    deallocate(count_num)

    a = 0
    do i = 1, n
      if (all(elite_num(1:c) /= i)) then
        a = a + 1
        local_num(a) = i
      end if
    end do

    !-----------------------------------パレート保存戦略
    elite_group(1:c, 1:tar) = order_list(elite_num(1:c), 1:tar)
    elite_fit(1:c) = fit(elite_num(1:c))
    elite_RCS(1:c) = total_RCS(elite_num(1:c))
    elite_flight_time(1:c) = total_flight_time(elite_num(1:c))

    !-----------------------------------その他
    local_group(1:n-c, 1:tar) = order_list(local_num(1:n-c), 1:tar)
    local_fit(1:n-c) = fit(local_num(1:n-c))
    local_RCS(1:n-c) = total_RCS(local_num(1:n-c))

    deallocate(elite_num)
    deallocate(local_num)

  end subroutine pareto_hyouka
end module pareto_hyouka_mod

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

module cumlation_mod                                        ! 確率を計算するmod
  implicit none

contains
  subroutine cumlation(probability, local_fit, c)
    use globals

    integer, intent(in) :: c
    real(8), intent(in) :: local_fit(:)
    real(8), intent(out), allocatable :: probability(:)
    integer i
    real(8) tot

    !---------------------------------------------------------------------------
    tot = sum(local_fit(1:n-c))                             ! 適応度totを計算→ここでは距離の合計が適応度totになる

    !---------------------------------------------------------------------------
    allocate(probability(n - c))
    probability(1:n-c) = local_fit(1:n-c) / tot             ! 選択確率を計算

  end subroutine cumlation
end module cumlation_mod

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

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

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

module crossover_mod                                        ! 交叉モジュール
  implicit none

contains
  subroutine nc_ket(nc, ns)                                 ! nc決定サブルーチン（配列形状決定）
    use globals
    use random_mod

    integer, intent(out) :: nc
    integer, intent(out) :: ns(nselect)
    integer empty(nselect)                                  ! ns以外
    integer i
    integer a
    real(8) r

    ns(1:nselect) = 0
    empty(1:nselect) = 0

    !-----------------------------------親の選択
    do
      nc = 0                                                ! 0クリア
      do i = 1, nselect
        call random_single(r)                               ! 乱数生成
        if (r < crate) then                                 ! 発生した乱数がcrateより小さいかチェック
          nc = nc + 1
          ns(nc) = i                                        ! 小さければ生成された乱数に対応する染色体番号がns()に蓄えられる
        else
          empty(i - nc) = i
        end if                                              ! ここでncはカウンター
      end do
      if (nc == 0) then
        cycle
      end if
      exit
    end do

    !---------------------------------------------------------------------------
    if (mod(nc, 2) == 1) then                               ! 変数ncが偶数か判定
      call random_single(r)                                 ! 奇数ならば乱数を新たに作成しもう一つの染色体番号を決定
      a = int(real(nselect - nc) * r) + 1

      nc = nc + 1
      ns(nc) = empty(a)
    end if

  end subroutine nc_ket

  !-----------------------------------------------------------------------------
  subroutine crossover(select_order, select_fit, newgen, nc, ns)
    use globals
    use random_mod

    integer, intent(in) :: select_order(nselect, tar)
    real(8), intent(in) :: select_fit(nselect)
    integer, intent(in) :: nc
    integer, intent(in) :: ns(nselect)
    integer, allocatable, intent(out) :: newgen(:, :)
    integer, parameter :: d = tar / 2
    integer i, j, k
    integer c, lis
    real(8) r
    integer, allocatable :: cross_par(:, :)
    integer na(d)
    integer newch(nc, tar)
    integer nuki1(d), nuki2(d)
    integer aa, pattern
    integer cross_point

    allocate(newgen(nc, tar))

    !-----------------------------------cross_parに親候補を格納
    allocate(cross_par(nc, tar))

    cross_par(1:nc, 1:tar) = select_order(ns(1:nc), 1:tar)

    !-----------------------------------交叉操作
    aa = 1
    do k = 1, nc / 2
      !---------------------------------交叉パターン決定
      ! 2つの親が1つでも同じデブリ番号を持っていたら順序交叉(pattern = 1)
      ! それ以外は1点交叉(pattern = 2)
      pattern = 2
      do i = 1, tar
        if (any(cross_par(aa, 1:tar) == cross_par(aa + 1, i))) then
          pattern = 1
          exit
        end if
      end do
      print '(a, i1)', 'pattern=', pattern

      if (pattern == 1) then
        !-------------------------------交叉操作 パターン1(順序交叉)
        na(1:d) = 0                                         ! 0クリア
        nuki1(1:d) = 0

        do i = 1, d                                         ! ランダムにd個の整数を作る
          call random_single(r)
          lis = int(real(tar) * r) + 1
          if (i > 1) then
            do
              if (any(na(1:i-1) == lis)) then
                call random_single(r)
                lis = int(real(tar) * r) + 1
              else
                exit
              end if
            end do
          end if
          na(i) = lis
        end do

        nuki1(1:d) = cross_par(aa + 1, na(1:d))             ! 親p1からd個の整数に対応する位置の数をnuki1()に格納

        !-------------------------------親p2からnuki1()の値を除いた(ここでは0にする)部分を遺伝させる
        do i = 1, tar
          newch(aa, i) = cross_par(aa, i)
          do j = 1, d
            if (cross_par(aa, i) == nuki1(j)) then
              newch(aa, i) = 0
              exit
            end if
          end do
        end do

        !-------------------------------nuki1を順番に0のところに入れていく
        c = 1
        do i = 1, tar
          if (newch(aa, i) == 0) then
            newch(aa, i) = nuki1(c)
            c = c + 1
          end if
        end do

        !-----------------------------------------------------------------------
        na(1:d) = 0                                         ! 0クリア
        nuki2(1:d) = 0

        do i = 1, d                                         ! ランダムにd個の整数を作る
          call random_single(r)
          lis = int(real(tar) * r) + 1
          if (i > 1) then
            do
              if (any(na(1:i-1) == lis)) then
                call random_single(r)
                lis = int(real(tar) * r) + 1
              else
                exit
              end if
            end do
          end if
          na(i) = lis
        end do

        nuki2(1:d) = cross_par(aa, na(1:d))                 ! 親p1から5個の整数に対応する位置の数をnuki1()に格納

        !-------------------------------親p2からnuki1()の値を除いた(ここでは0にする)部分を遺伝させる
        do i = 1, tar
          newch(aa + 1, i) = cross_par(aa + 1, i)
          do j = 1, d
            if (cross_par(aa + 1, i) == nuki2(j)) then
              newch(aa + 1, i) = 0
              exit
            end if
          end do
        end do

        c = 1
        do i = 1, tar
          if (newch(aa + 1, i) == 0) then
            newch(aa + 1, i) = nuki2(c)
            c = c + 1
          end if
        end do
        !-------------------------------交叉操作 パターン1 ここまで
      else
        !-------------------------------交叉操作 パターン2(1点交叉)
        cross_point = 0                                     ! 0クリア

        call random_single(r)
        cross_point = int(real(tar - 1) * r) + 1

        ! 子1生成
        newch(aa, 1:cross_point) = cross_par(aa, 1:cross_point)
        newch(aa, cross_point+1:tar) = cross_par(aa + 1, cross_point+1:tar)

        ! 子2生成
        newch(aa + 1, 1:cross_point) = cross_par(aa + 1, 1:cross_point)
        newch(aa + 1, cross_point+1:tar) = cross_par(aa, cross_point+1:tar)
        !-------------------------------交叉操作 パターン2 ここまで
      end if

      ! 子をnewgenに格納
      newgen(aa:aa+1, 1:tar) = newch(aa:aa+1, 1:tar)

      aa = aa + 2
    end do

  end subroutine crossover
end module crossover_mod

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

module mttn_mod
  implicit none

contains
  subroutine mttn(nc, newgen, mttn_group)
    use globals
    use random_mod

    integer, intent(in) :: nc
    integer, intent(in) :: newgen(nc, tar)
    integer, intent(out), allocatable :: mttn_group(:, :)
    integer i, j
    integer a, b, c1, c2
    real(8) r
    integer as(nc, tar)

    allocate(mttn_group(nc, tar))

    as(1:nc, 1:tar) = newgen(1:nc, 1:tar)

    do i = 1, nc
      call random_single(r)
      if (r < mttn_rate) then
        call random_single(r)
        a = int(real(tar) * r) + 1
        call random_single(r)
        b = int(real(tar - 1) * r) + 1
        if (b == a) then
          b = b + 1
        end if

        c1 = as(i, a)
        c2 = as(i, b)

        as(i, a) = c2
        as(i, b) = c1
      end if
    end do

    mttn_group(1:nc, 1:tar) = as(1:nc, 1:tar)

  end subroutine mttn
end module mttn_mod

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

module add_mod                                              ! 集団追加
  implicit none

contains
  subroutine addition(add_order)
    use globals
    use random_mod

    integer, intent(out) :: add_order(tar)
    integer i, j
    integer num
    real(8) r

    do i = 1, tar
      call random_single(r)
      num = int(real(db) * r) + 1
      if (i > 1) then
        do
          if (any(add_order(1:i-1) == num)) then
            call random_single(r)
            num = int(real(db) * r) + 1
          else
            exit
          end if
        end do
      end if
      add_order(i) = num
    end do

  end subroutine addition
end module add_mod

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

module output_mod
  implicit none

contains
  subroutine output_result(elite_fit, elite_RCS, c, k, elite_flight_time, elite_group)
    use globals

    real(8), intent(in) :: elite_fit(:)
    real(8), intent(in) :: elite_RCS(:)
    integer, intent(in) :: elite_group(:, :)
    integer, intent(in) :: c
    integer, intent(in) :: k
    real(8), intent(in) :: elite_flight_time(:)
    integer i
    character(len = 44) fmt
    write(fmt, '(a, i3, a)') '(f12.4, ",", f9.4, ",", f12.6, ', tar + 1, '(",", i3))'

    if (k == 1) then
      open(50, file = '../result/delta_V_result.csv', status = 'old')
      read(50, '()')
      read(50, '()')
    else
      open(50, file = '../result/delta_V_result.csv', position = 'append')
    end if

    do i = 1, c
      write(50, fmt) elite_fit(i), elite_RCS(i), elite_flight_time(i), k, elite_group(i, 1:tar)
    end do
    close(50)

  end subroutine output_result
end module output_mod

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

program main
  use globals
  use random_mod
  use debri_order_mod
  use hyouka_mod
  use pareto_hyouka_mod
  use cumlation_mod
  use roulette_mod
  use crossover_mod
  use mttn_mod
  use add_mod
  use output_mod

  implicit none

  integer i, j, k
  integer a
  integer c
  real(8) time_begin, time_end                             ! 時間計測のための変数
  integer order_list(n, tar)                               ! 初期集団
  integer new_order_list(n, tar)                           ! 次世代集団
  real(8) fit(n)                                           ! 適応度
  real(8) total_RCS(n)                                     ! 合計RCS
  real(8) total_flight_time(n)                             ! 合計飛行時間
  real(8), allocatable :: elite_flight_time(:)
  real(8), allocatable :: elite_RCS(:)                     ! 上位5RCS
  real(8), allocatable :: local_RCS(:)
  integer, allocatable :: elite_group(:, :)                ! 上位5個体
  integer, allocatable :: local_group(:, :)                ! 上位5個体を取り除いたもの
  real(8), allocatable :: elite_fit(:)                     ! 評価値の上位5個体
  real(8), allocatable :: local_fit(:)                     ! 上位5個体以外の評価値
  real(8), allocatable :: probability(:)                   ! 確率
  integer select_order(nselect, tar)                       ! 選択配列
  real(8) select_fit(nselect)
  integer nc
  integer ns(nselect)
  integer, allocatable :: newgen(:, :)
  integer, allocatable :: mttn_group(:, :)
  integer add_order(tar)
  integer min_root(tar)
  character(len = 16) fmt
  write(fmt, '(a, i3, a)') '(f8.4, a, ', tar, 'i5)'

  !-------------------------------------時間計測（開始）
  call cpu_time(time_begin)

  !-------------------------------------乱数シード設定
  call set_seed

  !-------------------------------------デブリデータ読み込み
  call read_debridata

  !-------------------------------------初期集団生成
  call debri_order_table(order_list)

  ! open(100, file = 'test.txt', status = 'replace')
  ! close(100)

  do k = 1, T
    ! open(101, file = 'test.txt', position = 'append')
    ! do i = 1, n
    !   write(101, '(i0, 5(a, i3))') &
    !  & k, ', ', order_list(i, 1), ', ', order_list(i, 2), ', ', order_list(i, 3), ', ', &
    !  & order_list(i, 4), ', ', order_list(i, 5)
    ! end do
    ! close(101)

    !-------------------------------------評価
    call hyouka(order_list, fit, total_RCS, total_flight_time, k)

    call pareto_hyouka(order_list, fit, elite_group, local_group, elite_fit, &
   & local_fit, total_RCS, elite_RCS, local_RCS, c, elite_flight_time, total_flight_time)

    write(*, '(a, i3)') 'elite:', c
    ! do i = 1, c
    !   write(*, fmt) elite_fit(i), ' km/s', elite_group(i, 1:tar)
    ! end do
    do i = 1, 5
      write(*, fmt) local_fit(i), ' km/s', local_group(i, 1:tar)
    end do

    !-------------------------------------出力
    ! call output_result(elite_fit, elite_RCS, c, k, elite_flight_time, elite_group)
    call output_result(local_fit, local_RCS, 5, k, total_flight_time, local_group)

    !-------------------------------------終了判定
    if (k == T) then
      write(*, '(a)') ''
      write(*, '(a)') '---------------------------------'
      write(*, '(a23)') 'Calculate End'
      write(*, '(a11, i6, a)') 'Total:' , k , ' Generation'
      write(*, '(a)') '---------------------------------'
      exit
    end if

    !-------------------------------------選択
    call cumlation(probability, local_fit, c)
    call roulette_par(local_group, probability, select_order, local_fit, select_fit, c)

    !-------------------------------------交叉
    call nc_ket(nc, ns)
    write(*, '(a, i3)') 'cross:', nc

    call crossover(select_order, select_fit, newgen, nc, ns)

    !-------------------------------------変異
    call mttn(nc, newgen, mttn_group)

    !-------------------------------------次世代へ受け継ぐ
    new_order_list(1:c, 1:tar) = elite_group(1:c, 1:tar)
    new_order_list(c+1:c+nc, 1:tar) = mttn_group(1:nc, 1:tar)

    !-------------------------------------集団追加
    if (c + nc < n) then
      do i = c + nc + 1, n
        call addition(add_order)
        new_order_list(i, 1:tar) = add_order(1:tar)
      end do
    end if

    order_list(1:n, 1:tar) = new_order_list(1:n, 1:tar)

    write(*, '(a)') ''
    write(*, '(a)') '---------------------------------'
    write(*, '(a11, i6, a)') 'Finish' , k , ' Generation'
    write(*, '(a)') '---------------------------------'
    deallocate(newgen)                                      !配列の初期化
    deallocate(mttn_group)
    deallocate(elite_group)
    deallocate(elite_fit)
    deallocate(elite_RCS)
    deallocate(local_group)
    deallocate(local_fit)
    deallocate(local_RCS)
    deallocate(elite_flight_time)
    deallocate(probability)
  end do

  !-------------------------------------時間計測（終了）
  call cpu_time(time_end)

  write(*, '(a, f10.5)') 'Calculate Time:' , time_end - time_begin

end program main
