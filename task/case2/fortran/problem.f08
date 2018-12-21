module problem
use iso_c_binding
use util
use orbit_func
use orbit_base
use orbit_debri
implicit none

real(8) :: base_time
integer, allocatable :: ORDER(:)

contains

! ##############################################################################
! # IO
! ##############################################################################

! 初期化(兼fortranの動作確認用)
subroutine initialize() bind(c, name="initialize")
    implicit none

    ! 基準時刻を設定
    base_time = mjd(2018, 0d0)
    print *, 'initialize(fortran)', base_time
end subroutine initialize


! デブリデータ初期化
subroutine init_debri(tle_file, rcs_file, n_debris) bind(c, name="init_debri")
    implicit none
    character(1, c_char), intent(in) :: tle_file(*)
    character(1, c_char), intent(in) :: rcs_file(*)
    integer(c_int), intent(out) :: n_debris
    integer :: i

    ! ファイルからデブリデータを読み込み
    call load_debri(-1, fstring(tle_file), fstring(rcs_file))
    ! base_time = DEBRIS(1)%orbit%epc

    ! デブリ位置を基準時刻のものに初期化
    do i = 0, NUM_DEBRIS_ALL
      call DEBRIS(i)%orbit%move(base_time)
    end do

    ! デブリデータを昇交点赤系昇順に並べ替え
    call select_debri(NUM_DEBRIS_ALL, ORDER)

    ! カタログ上のデブリデータ数を返す
    n_debris = NUM_DEBRIS_ALL
end subroutine init_debri


! 目的関数(ΔV)計算
subroutine call_problem(iargs, n1, dargs, n2, result) bind(c, name="call_problem")
    ! ### 引数 ###
    ! iargs: デブリ番号の列 # [例] (1, 2, 3, 4) => 1 -> 2, 2 -> 3, 3 -> 4
    ! dargs(2n-1): 待機時間 [day]
    ! dargs(2n): 遷移時間 [s]
    ! result(1): ΔV [km/s]
    ! result(2): RCS [m^2]

    implicit none
    integer(c_int), intent(in), value :: n1, n2
    integer(c_int), intent(in) :: iargs(n1) ! 32bit整数型配列引数
    real(c_double), intent(in) :: dargs(n2) ! 64bit実数型配列引数
    real(c_double), intent(out) :: result(2)
    type(TOrbit) :: orbit1, orbit2, orbit3
    real(8) :: total_delta_v, total_delta_t
    real(8) :: date_start, del_v, del_t
    real(8) :: rho, dnu, theta(3), ncycle(3)
    integer :: n_iargs, n_dargs

    ! 引数チェック
    n_iargs = 2
    n_dargs = 7
    if (n1 /= n_iargs .or. n2 /= n_dargs) then
        print *, "ArgumentError: invalid argment size"
        call exit(1)
    end if

    ! 計算結果格納変数
    total_delta_v = 0d0
    total_delta_t = 0d0

    ! 計算用軌道設定
    orbit1 = DEBRIS(iargs(1))%orbit
    orbit3 = DEBRIS(iargs(2))%orbit

    ! パラメータ設定
    date_start = dargs(1)
    rho = dargs(2)
    dnu = dargs(3)
    theta = [dargs(4:5), 1d0]
    ncycle = [dargs(6:7), 0d0]

    ! 軌道遷移
    call orbit_transfer(orbit1, orbit3, date_start, rho, theta, ncycle, del_v, del_t, orbit2)
    total_delta_v = total_delta_v + del_v
    total_delta_t = total_delta_t + del_t

    ! 対象に接触
    date_start = date_start + total_delta_t
    call orbit_approach(orbit2, orbit3, date_start, dnu, del_v, del_t)
    total_delta_v = total_delta_v + del_v
    total_delta_t = total_delta_t + del_t

    ! 結果格納
    result(1) = total_delta_v
    result(2) = total_delta_t
end subroutine call_problem


! ##############################################################################
! # Helpers
! ##############################################################################

! デブリデータを昇交点赤系昇順に並べ替え・一定個数選択
subroutine select_debri(n, order)
    implicit none
    integer, intent(in) :: n
    integer, intent(out), allocatable :: order(:)
    real(8), allocatable :: data(:)

    allocate(data(NUM_DEBRIS_ALL), source=DEBRIS(1:)%orbit%raan)
    allocate(order(0:n), source=[0, sort(data, n)])
end subroutine select_debri


! 軌道遷移・ΔV計算
subroutine orbit_transfer(orbit_dep, orbit_arr, date_s, rho, theta, ncycle, del_v, del_t, orbit_res)
    ! ### 引数 ###
    ! orbit_dep: 出発軌道
    ! orbit_arr: 到達軌道
    ! date_s: 出発日時
    ! rho: 出発軌道と遷移軌道の長半径比 (a_tra / a_dep)
    ! theta: 軌道面変更角の比
    ! ncycle: 時間調整用周回数

    ! ### 戻り値 ###
    ! del_v: 合計速度増分 [km/s]
    ! del_t: 合計遷移時間 [s]
    ! orbit_res: 遷移終了軌道

    implicit none
    type(TOrbit), intent(in) :: orbit_dep, orbit_arr
    real(8), intent(in)  :: date_s, rho, theta(3), ncycle(3)
    real(8), intent(out) :: del_v, del_t
    type(TOrbit), intent(out) :: orbit_res
    type(TOrbit) :: orbit1, orbit2, orbit_t
    real(8) :: date, date_d, date_a, dtheta, axis(3), axis_t(3), offset
    real(8) :: sma_t, pos_d(3), pos_a(3), vel_d(3), vel_t(3)
    integer :: i

    ! 計算結果格納変数
    del_v = 0d0
    del_t = 0d0

    ! 計算用時刻設定
    date = base_time + date_s

    ! 計算用軌道設定
    orbit1 = orbit_dep
    orbit2 = orbit_arr

    ! 衛星を基準時刻地点に移動
    call orbit1%move(date)
    call orbit2%move(date)

    do i = 1, 3
        ! 出発軌道と遷移軌道の軌道面交差軸
        axis = unit(cross_product(orbit1%axis_w, orbit2%axis_w)) * (-1) ** (i - 1)

        ! 軌道面変更角
        dtheta = sgn(dot_product(cross_product(orbit1%axis_w, orbit2%axis_w), axis)) &
            * calc_angle(orbit1%axis_w, orbit2%axis_w)

        ! 衛星が遷移開始地点に到達するまでの時間 [s]
        offset = calc_angle(convert_anomaly(calc_angle(orbit1%axis_p, axis, orbit1%axis_w), &
            orbit1%eccentricity, from='T', to='M') - orbit1%mean_anomaly) / orbit1%mean_motion

        ! 遷移開始日時
        date_d = date + offset * SEC2DAY

        ! 衛星を遷移開始地点に移動
        call orbit1%move(date_d)

        ! 遷移開始位置ベクトル
        pos_d = orbit1%get_position()

        ! 遷移開始速度ベクトル
        vel_d = orbit1%get_velocity()

        select case (i)
        case (1)
            ! 遷移軌道の長半径
            sma_t = 0.5d0 * norm(pos_d) * rho + 1d-10
        case (2)
            ! 軌道遷移時間 [s]
            offset = calc_angle(convert_anomaly(calc_angle(orbit2%axis_p, axis, orbit2%axis_w), &
                orbit2%eccentricity, from='T', to='M') - orbit1%mean_anomaly) / orbit1%mean_motion

            ! 目標地点到達時刻
            date_a = date + offset * SEC2DAY

            ! 目標地点位置ベクトル
            call orbit2%move(date_a)

            ! 目標地点位置ベクトル
            pos_a = orbit2%get_position()

            ! 遷移軌道の長半径
            sma_t = 0.5d0 * (norm(pos_d) + norm(pos_a)) + 1d-10
            ! sma_t = (norm(pos_d) + orbit2%semimajor_axis - 30d0) / 2d0
        case (3)
            ! 遷移軌道の長半径
            sma_t = norm(pos_d) + 1d-10
        end select

        ! 遷移軌道面直交軸
        axis_t = rotate(orbit1%axis_w, axis, dtheta * theta(i))

        ! 遷移軌道
        orbit_t = get_orbit(pos_d, sma_t, axis_t, date_d)

        ! 速度変更直後の速度ベクトル
        vel_t = orbit_t%get_velocity()

        ! 速度増分
        del_v = del_v + norm(vel_t - vel_d)

        ! 時刻更新
        date = date_d + orbit_t%period * ncycle(i) * SEC2DAY

        ! 軌道更新
        orbit1 = orbit_t
        call orbit1%move(date)
        call orbit2%move(date)
    end do

    ! 最終軌道更新
    orbit_res = orbit1
end subroutine orbit_transfer


! 位相差を利用して目標に接近
subroutine orbit_approach(orbit_dep, orbit_arr, date_s, dnu, del_v, del_t)
    ! ### 引数 ###
    ! orbit_dep: 出発軌道
    ! orbit_arr: 到達軌道
    ! date_s: 出発日時
    ! dnu: 遷移軌道の中心角比

    ! ### 戻り値 ###
    ! del_v: 合計速度増分 [km/s]
    ! del_t: 合計遷移時間 [s]

    implicit none
    type(TOrbit), intent(in) :: orbit_dep, orbit_arr
    real(8), intent(in)  :: date_s, dnu
    real(8), intent(out) :: del_v, del_t
    type(TOrbit) :: orbit1, orbit2
    real(8) :: date, date_d, offset
    real(8) :: tof, pos_d(3), pos_a(3), vel_d(3), vel_a(3), v1(3), v2(3)

    ! 計算結果格納変数
    del_v = 0d0
    del_t = 0d0

    ! 計算用時刻設定
    date = base_time + date_s

    ! 計算用軌道設定
    orbit1 = orbit_dep
    orbit2 = orbit_arr

    ! 衛星を基準時刻地点に移動
    call orbit1%move(date)
    call orbit2%move(date)

    ! 軌道遷移時間を設定
    tof = calc_period(0.5d0 * (orbit1%semimajor_axis + orbit2%semimajor_axis)) * dnu

    ! 衛星を接触地点(仮)に移動
    call orbit1%move(date + tof * SEC2DAY)
    call orbit2%move(date + tof * SEC2DAY)

    ! 接触地点(仮)と実際の接触地点の差を計算 [s]
    offset = calc_angle(sgn(orbit2%mean_motion - orbit1%mean_motion) &
           * (convert_anomaly(calc_angle(orbit2%axis_p, orbit1%get_position(), orbit2%axis_w), &
              orbit2%eccentricity, from='T', to='M') - orbit2%mean_anomaly))                  &
              / abs(orbit2%mean_motion - orbit1%mean_motion)

    ! 実際の遷移開始時刻を計算
    date_d = date + offset * SEC2DAY + 1d-1

    ! 遷移開始地点の位置と速度を計算
    call orbit1%move(date_d)
    pos_d = orbit1%get_position()
    vel_d = orbit1%get_velocity()

    ! 遷移終了地点の位置と速度を計算
    call orbit2%move(date_d + tof * SEC2DAY)
    pos_a = orbit2%get_position()
    vel_a = orbit2%get_velocity()

    ! 速度増分と遷移時間を計算
    call lambert(date_d, tof, pos_d, pos_a, orbit2%axis_w, v1, v2, orbit1)
    del_v = del_v + norm(v1 - vel_d) + norm(vel_a - v2)
    del_t = date_d - date_s + tof * SEC2DAY
end subroutine orbit_approach


! ##############################################################################
! # Utilities
! ##############################################################################

function fstring(string)
    ! c言語の文字列ポインタをFortranの文字列配列に変換する
    implicit none
    character(1, c_char), intent(in) :: string(*)
    character(:, c_char), allocatable :: fstring
    integer :: len, i

    len = 1
    do while (string(len) /= c_null_char)
        len = len + 1
    end do
    len = len - 1
    allocate(character(len, c_char) :: fstring)
    do i = 1, len
      fstring(i:i) = string(i)
    end do
end function fstring

end module problem
