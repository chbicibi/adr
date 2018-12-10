module orbit_base
  use util
  use orbit_func
  implicit none

  private
  public :: TOrbit
  public :: MU, J2, SEC_DAY, AXIS_I, AXIS_J, AXIS_K
  public :: GetOrbit, calc_period, convert_anomaly, hohmann, lambert
  public :: intersection

  type :: TOrbit
    real(8) :: epc       = -1d0
    real(8) :: inc       = -1d0
    real(8) :: ran       = -1d0
    real(8) :: ecc       = -1d0
    real(8) :: ap        = -1d0
    real(8) :: ma        = -1d0
    real(8) :: mm        = -1d0

    integer :: type      = 1
    real(8) :: sma       = -1d0
    real(8) :: slr       = -1d0
    real(8) :: prd       = -1d0
    real(8) :: axis_p(3) = 0d0
    real(8) :: axis_q(3) = 0d0
    real(8) :: axis_w(3) = 0d0
    real(8) :: amom(3)   = 0d0 ! angular momentum, h

    real(8) :: dt        = -1d0
    real(8) :: ea        = -1d0 ! eccentric Anomaly, E /
    real(8) :: ta        = -1d0 ! true anomaly, f

    real(8) :: pos(3)    = 0d0 ! position, r
    real(8) :: vel(3)    = 0d0 ! velocity, v

    real(8) :: dran      = 0d0
    real(8) :: dap       = 0d0

    contains

    generic :: initialize => initialize_elm, initialize_tle

    procedure :: initialize_elm, initialize_tle
    procedure :: CalcOrbit
    procedure :: get_position
    procedure :: get_velocity
    procedure :: move
    procedure :: save
  end type TOrbit

  interface GetOrbit
    procedure :: CreateOrbit, PosToOrbit, PeriToOrbit
  end interface GetOrbit

  interface lambert
    procedure :: lambert_orbit, lambert_core
  end interface lambert

  real(8), parameter :: MU        = 3.9860044d5
  real(8), parameter :: J2        = 4.4041964d4 ! => J2 * re ** 2
  real(8), parameter :: SEC_DAY   = 1d0 / 86400d0
  real(8), parameter :: AXIS_I(3) = [1d0, 0d0, 0d0]
  real(8), parameter :: AXIS_J(3) = [0d0, 1d0, 0d0]
  real(8), parameter :: AXIS_K(3) = [0d0, 0d0, 1d0]

  contains

  subroutine initialize_elm(this, epc, inc, ran, ecc, ap, ma, mm)
    implicit none
    class(TOrbit), intent(inout) :: this
    real(8), intent(in) :: epc, inc, ran, ecc, ap, ma, mm
    real(8) :: matrix(3, 3)

    this%epc    = epc
    this%inc    = inc * RADIANS
    this%ran    = ran * RADIANS
    this%ecc    = ecc
    this%ap     = ap * RADIANS
    this%ma     = ma * RADIANS
    this%mm     = mm * PI2 * SEC_DAY
    this%sma    = (MU / this%mm ** 2) ** (1d0 / 3d0)
    this%slr    = this%sma * (1d0 - ecc ** 2)
    this%prd    = 86400d0 / mm

    matrix      = EulerAngle([this%ap, this%inc, this%ran], [3, 1, 3])
    this%axis_p = matrix(1, :)
    this%axis_q = matrix(2, :)
    this%axis_w = matrix(3, :)

    this%dran   = -1.5d0 * J2 * this%mm / this%slr ** 2 * cos(this%inc)
    this%dap    =  1.5d0 * J2 * this%mm / this%slr ** 2 * (2d0 - 2.5d0 * sin(this%inc) ** 2)
  end subroutine initialize_elm

  subroutine initialize_tle(this, line1, line2)
    implicit none
    class(TOrbit), intent(inout) :: this
    character(69), intent(inout) :: line1, line2
    real(8) :: day, inc, ran, ecc, ap, ma, mm
    integer :: year

    read(line1(19:20), *) year
    read(line1(21:32), *) day

    read(line2( 9:16), *) inc
    read(line2(18:25), *) ran
    line2(25:26) = "0."
    read(line2(25:33), *) ecc
    read(line2(35:42), *) ap
    read(line2(44:51), *) ma
    read(line2(53:63), *) mm

    if (year < 57) then
      year = year + 2000
    else
      year = year + 1900
    end if

    call this%initialize(mjd(year, day), inc, ran, ecc, ap, ma, mm)
  end subroutine initialize_tle


  ! ============================================================================
  ! ============================================================================

  !-----------------------------------------------------------------------------
  ! Constructor
  !-----------------------------------------------------------------------------
  type(TOrbit) function CreateOrbit(epc, inc, ran, ecc, ap, ma, mm) result (orbit)
    implicit none
    real(8), intent(in) :: epc, inc, ran, ecc, ap, ma, mm
    real(8)             :: matrix(3, 3)

    orbit%epc    = epc
    orbit%inc    = inc * RADIANS
    orbit%ran    = ran * RADIANS
    orbit%ecc    = ecc
    orbit%ap     = ap   * RADIANS
    orbit%ma     = ma   * RADIANS
    orbit%mm     = mm   * PI2 * SEC_DAY
    orbit%sma    = (MU / orbit%mm ** 2) ** (1d0 / 3d0)
    orbit%slr    = orbit%sma * (1d0 - ecc ** 2)
    orbit%prd    = 86400d0 / mm

    matrix       = EulerAngle([orbit%ap, orbit%inc, orbit%ran], [3, 1, 3])
    orbit%axis_p = matrix(1, :)
    orbit%axis_q = matrix(2, :)
    orbit%axis_w = matrix(3, :)

    orbit%dran   = -1.5d0 * J2 * orbit%mm / orbit%slr ** 2 * cos(orbit%inc)
    orbit%dap    =  1.5d0 * J2 * orbit%mm / orbit%slr ** 2 * (2d0 - 2.5d0 * sin(orbit%inc) ** 2)
  end function CreateOrbit

  type(TOrbit) function PosToOrbit(pos, vel, date) result (orbit)
    implicit none
    real(8), intent(in) :: pos(3), vel(3), date
    real(8)             :: r, sma, ecc, esin, ecos, ea
    real(8)             :: h(3), axis_p(3), axis_q(3), axis_w(3), axis_n(3)
    integer             :: type

    r      = norm(pos)
    sma    = 1d0 / (2d0 / r - dot_product(vel, vel) / MU)
    h      = cross_product(pos, vel)
    axis_p = Unit(-MU * pos / r - cross_product(h, vel))
    axis_w = Unit(h)
    axis_q = cross_product(axis_w, axis_p)
    axis_n = Unit(cross_product(AXIS_K, axis_w))
    ecc    = sqrt(1d0 - dot_product(h, h) / (MU * sma))
    type   = Sgn(1d0 - ecc)
    esin   = dot_product(pos, vel) / sqrt(mu * sma)

    if (ecc > 0) then
      if (type > 0) then
        ecos = 1d0 - r / sma
        ea   = calc_angle(ecos / ecc, esin)
      else
        ea   = asinh(esin / ecc)
      end if
    else
      ea     = 0d0
    end if

    orbit%type   = type
    orbit%epc    = date
    orbit%inc    = calc_angle(AXIS_K, axis_w)
    orbit%ran    = calc_angle(AXIS_I, axis_n, AXIS_K)
    orbit%ecc    = ecc
    orbit%ap     = calc_angle(axis_n, axis_p, axis_w)
    orbit%ma     = type * (ea - esin)
    orbit%mm     = PI2 / calc_period(sma)
    orbit%sma    = sma
    orbit%slr    = sma * (1d0 - ecc ** 2)
    orbit%prd    = PI2 / orbit%mm
    orbit%axis_p = axis_p
    orbit%axis_q = axis_q
    orbit%axis_w = axis_w
    ! orbit%amom   = h
    ! orbit%dt     = orbit%ma / orbit%mm
    ! orbit%ea     = ea
    orbit%pos    = pos
    orbit%vel    = vel

    orbit%dran   = -1.5d0 * J2 * orbit%mm / orbit%slr ** 2 * cos(orbit%inc)
    orbit%dap    =  1.5d0 * J2 * orbit%mm / orbit%slr ** 2 * (2d0 - 2.5d0 * sin(orbit%inc) ** 2)
  end function PosToOrbit

  type(TOrbit) function PeriToOrbit(pos, sma, axis_w, date) result (orbit)
    implicit none
    real(8), intent(in) :: pos(3), sma, axis_w(3), date
    real(8)  :: r, ecc, axis_p(3), axis_q(3), axis_n(3)
    integer  :: dir

    r      = norm(pos)
    dir    = Sgn(sma - r)
    axis_p = dir * Unit(pos)
    axis_q = cross_product(axis_w, axis_p)
    axis_n = Unit(cross_product(AXIS_K, axis_w))
    ecc    = dir * (1d0 - r / sma)

    orbit%type   = Sgn(1d0 - ecc)
    orbit%epc    = date
    orbit%ecc    = ecc
    orbit%ma     = 0.5d0 * PI * (1 - dir)
    orbit%mm     = PI2 / calc_period(sma)
    orbit%sma    = sma
    orbit%slr    = sma * (1d0 - ecc ** 2)
    orbit%prd = PI2 / orbit%mm
    orbit%axis_p = axis_p
    orbit%axis_q = axis_q
    orbit%axis_w = axis_w
    orbit%inc    = calc_angle(AXIS_K, axis_w)
    orbit%ran    = calc_angle(AXIS_I, axis_n, AXIS_K)
    orbit%ap     = calc_angle(axis_n, axis_p, axis_w)

    orbit%dran   = -1.5d0 * J2 * orbit%mm / orbit%slr ** 2 * cos(orbit%inc)
    orbit%dap    =  1.5d0 * J2 * orbit%mm / orbit%slr ** 2 * (2d0 - 2.5d0 * sin(orbit%inc) ** 2)
  end function PeriToOrbit

  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  type(TOrbit) function CalcOrbit(this, date) result (orbit)
    implicit none
    class(TOrbit), intent(in) :: this
    real(8), intent(in) :: date
    real(8) :: offset, matrix(3, 3)

    orbit%type   = this%type
    orbit%inc    = this%inc
    orbit%ecc    = this%ecc
    orbit%mm     = this%mm
    orbit%sma    = this%sma
    orbit%slr    = this%slr
    orbit%prd    = this%prd
    orbit%dran   = this%dran
    orbit%dap    = this%dap

    offset       = 86400d0 * (date - this%epc)
    orbit%epc    = date
    orbit%ran    = modulo(this%ran + this%dran * offset, PI2)
    orbit%ap     = modulo(this%ap   + this%dap   * offset, PI2)
    orbit%ma     = modulo(this%ma   + this%mm    * offset, PI2)
    orbit%ea     = convert_anomaly(orbit%ma, this%ecc, from='M', to='E')
    matrix       = EulerAngle([orbit%ap, orbit%inc, orbit%ran], [3, 1, 3])
    orbit%axis_p = matrix(1, :)
    orbit%axis_q = matrix(2, :)
    orbit%axis_w = matrix(3, :)
    ! orbit%dt    = ma / this%mm
  end function CalcOrbit

  function get_position(this) result (pos)
    implicit none
    class(TOrbit), intent(inout) :: this
    real(8) :: pos(3), s, c

    if (this%ea < 0) this%ea = convert_anomaly(this%ma, this%ecc, from='M', to='E')

    if (this%type > 0) then
      s = sin(this%ea)
      c = cos(this%ea)
    else
      s = sinh(this%ea)
      c = cosh(this%ea)
    end if
    pos = this%sma * (this%type * (c - this%ecc) * this%axis_p &
                 + sqrt(1d0 - this%ecc ** 2) * s * this%axis_q)
  end function get_position

  function get_velocity(this) result (vel)
    implicit none
    class(TOrbit), intent(inout) :: this
    real(8) :: vel(3), r, s, c

    if (this%ea < 0) this%ea = convert_anomaly(this%ma, this%ecc, from='M', to='E')

    if (ALL(this%pos == 0d0)) then
      r = norm(this%get_position())
    else
      r = norm(this%pos)
    end if
    if (this%type > 0) then
      s = sin(this%ea)
      c = cos(this%ea)
    else
      s = sinh(this%ea)
      c = cosh(this%ea)
    end if
    vel = -sqrt(MU) / r * sqrt(this%sma) * (s * this%axis_p &
             -  sqrt(1d0 - this%ecc ** 2) * c * this%axis_q)
  end function get_velocity

  subroutine move(this, arg, mode)
    implicit none
    class(TOrbit), intent(inout) :: this
    real(8), intent(in) :: arg
    integer, intent(in) :: mode
    real(8) :: epc, dt, matrix(3, 3)

    select case (mode)
    case (0) ! 時刻 (day)
      epc = arg
      dt = (arg - this%epc) * 86400d0
    case (1) ! 変動時間 (day)
      epc = this%epc + arg
      dt = arg * 86400d0
    case (2) !
      epc = this%epc + this%prd * arg / PI2 / 86400d0
      dt = this%prd * arg / PI2
    case (-1)
      epc = this%epc
      dt = 0d0
    case default
      return
    end select

    this%epc    = epc
    this%ran    = calc_angle(this%ran + this%dran * dt)
    this%ap     = calc_angle(this%ap  + this%dap  * dt)
    this%ma     = calc_angle(this%ma  + this%mm   * dt)
    this%ea     = convert_anomaly(this%ma, this%ecc, from='M', to='E')
    matrix      = EulerAngle([this%ap, this%inc, this%ran], [3, 1, 3])
    this%axis_p = matrix(1, :)
    this%axis_q = matrix(2, :)
    this%axis_w = matrix(3, :)
    ! orbit%dt    = ma / this%mm
  end subroutine move

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

  real(8) function convert_anomaly(a, e, from, to) result (z)
    implicit none
    real(8), intent(in) :: a, e
    character, intent(in) :: from, to
    real(8) :: z_next, f, dfdz, c
    integer :: i

    z = a
    select case (from)
    case ('M')
      do i = 1, 100
        if (e < 1) then
          f = z   - e * sin(z) - a
          dfdz = 1d0 - e * cos(z)
        else
          f = e * sinh(z) - z - a
          dfdz = e * cosh(z) - 1d0
        end if
        z_next = z - (f / dfdz)

        if (abs(z_next - z) < 1d-6) exit
        z = z_next
      end do
    case ('T')
      c = cos(z)
      z = calc_angle((c + e) / (1d0 + e * c), PI - z)
    end select

    select case (to)
    case ('M')
      z = z - e * sin(z)
    case ('T')
      c = cos(z)
      z = calc_angle((c - e) / (1d0 - e * c), PI - z)
    end select
  end function convert_anomaly

  real(8) function calc_period(sma) result(period)
    implicit none
    real(8), intent(in) :: sma

    period = PI2 * sqrt(sma ** 3 / MU)
  end function calc_period


  ! ============================================================================
  ! IO
  ! ============================================================================

  subroutine save(this, date, dt, file, init)
    class(TOrbit), intent(in) :: this
    real(8), intent(in) :: date, dt
    character(*), intent(in) :: file
    logical, intent(in), optional :: init
    type(TOrbit) :: orbit
    ! real(8), allocatable :: data(:, :)
    real(8) :: pos(3), ea1, ea2, dea
    integer :: unit, n, i

    orbit = this%CalcOrbit(date)

    ea1 = convert_anomaly(orbit%ma + orbit%mm * min(dt, 0d0), orbit%ecc, from='M', to='E')
    ea2 = convert_anomaly(orbit%ma + orbit%mm * max(dt, 0d0), orbit%ecc, from='M', to='E')
    dea = ea2 - ea1

    n = int(100 * abs(dt) / orbit%prd)
    ! allocate(data(3, 0:n))
    ! print *, n
    ! print *, orbit%prd
    ! stop

    if (present(init) .and. init) then
      open(newunit=unit, file=file, status="replace")
    else
      open(newunit=unit, file=file, position="append")
    end if
      do i = 0, n
        orbit%ea = calc_angle(ea1 + dea * i / dble(n))
        ! orbit%ma = convert_anomaly(orbit%ea, orbit%ecc, from='E', to='M')
        ! call orbit%move(0d0, -1)
        pos = orbit%get_position()
        write(unit, "(2(f15.5',')f15.5)") pos
      end do
      write(unit, *)
    close(unit)
  end subroutine save


  ! ============================================================================
  ! orbital transfer
  ! ============================================================================

  subroutine hohmann(date, dep, arr, trf, date_d, dt, dv)
    implicit none
    real(8), intent(in) :: date
    type(TOrbit), intent(inout) :: dep, arr
    type(TOrbit), intent(out) :: trf
    real(8), intent(out) :: date_d, dt, dv
    type(TOrbit) :: orbit1, orbit2, orbit_t
    real(8) :: theta, offset, date_a, sma_t, ratio
    real(8) :: axis(3), axis_t(3), pos_d(3), vel_d(3), vel_a(3)
    integer :: i

    i = 1
    ratio = 0.5d0

    orbit1 = dep
    orbit2 = arr

    axis = Unit(cross_product(orbit1%axis_w, orbit2%axis_w)) * (-1) ** (i - 1)
    theta = Sgn(dot_product(cross_product(orbit1%axis_w, orbit2%axis_w), axis)) &
                  * calc_angle(orbit1%axis_w, orbit2%axis_w)

    offset = calc_angle(convert_anomaly(calc_angle(orbit1%axis_p, axis, orbit1%axis_w), &
                                        orbit1%ecc, from='T', to='M') - orbit1%ma) / orbit1%mm

    date_d = date + offset * SEC_DAY
    call orbit1%move(date_d, 0)
    pos_d = orbit1%get_position()
    vel_d = orbit1%get_velocity()

    date_a  = date + offset * SEC_DAY
    call orbit2%move(date_a, 0)
    vel_a = orbit2%get_velocity()

    sma_t = norm(pos_d)
    axis_t = Rotate(orbit1%axis_w, axis, ratio * theta)
    orbit_t = GetOrbit(pos_d, sma_t, axis_t, date_d)

    trf = orbit_t
    dt = 0.5d0 * trf%prd

    dv = norm(vel_d - orbit_t%get_velocity())
    call orbit_t%move(date_a, 0)
    dv = dv + norm(vel_a - orbit_t%get_velocity())
  end subroutine hohmann

  subroutine intersection(dep, arr, axis, theta, offset)
    implicit none
    type(TOrbit), intent(in) :: dep, arr
    real(8), intent(out) :: axis(3), theta, offset
    integer :: i

    i = 1

    axis    = Unit(cross_product(dep%axis_w, arr%axis_w)) * (-1) ** (i - 1)
    theta   = Sgn(dot_product(cross_product(dep%axis_w, arr%axis_w), axis)) &
                  * calc_angle(dep%axis_w, arr%axis_w)
    offset  = calc_angle(convert_anomaly(calc_angle(dep%axis_p, axis, dep%axis_w), &
                                         dep%ecc, from='T', to='M') - dep%ma) / dep%mm
  end subroutine intersection

  subroutine lambert_orbit(date, delta_t, dep, arr, delta_v, trf)
    implicit none
    real(8), intent(in) :: date, delta_t
    type(TOrbit), intent(inout) :: dep, arr
    real(8), intent(out) :: delta_v
    type(TOrbit), intent(out), optional :: trf
    real(8) :: pos1(3), pos2(3), vel1(3), vel2(3), axis(3)

    call dep%move(date, 0)
    call arr%move(date + delta_t * SEC_DAY, 0)
    pos1 = dep%get_position()
    pos2 = arr%get_position()
    axis = arr%axis_w
    call lambert(date, delta_t, pos1, pos2, axis, vel1, vel2, trf)

    delta_v = norm(vel1 - dep%get_velocity()) + norm(vel2 - arr%get_velocity())
  end subroutine lambert_orbit

  subroutine lambert_intermediate(date, delta_t, dep, arr, mid, delta_v, trf)
    implicit none
    real(8), intent(in) :: date, delta_t, mid(:, :)
    type(TOrbit), intent(inout) :: dep, arr
    real(8), intent(out) :: delta_v
    type(TOrbit), intent(out), optional :: trf
    real(8) :: pos1(3), pos2(3), vel1(3), vel2(3), axis(3)
    real(8), allocatable :: dv(:)
    integer :: s, i

    ! pos_s = dep%get_position()
    ! pos_d = arr%get_position()
    ! delta_nu = calc_angle(pos1, pos2, axis)

    s = size(mid, dim=2)
    allocate(dv(s + 1))

    do i = 1, s
      call dep%move(date, 0)
      call arr%move(date + delta_t * SEC_DAY, 0)
      pos1 = dep%get_position()
      pos2 = arr%get_position()
      axis = arr%axis_w
      call lambert(date, delta_t, pos1, pos2, axis, vel1, vel2, trf)

      dv(i) = norm(vel1 - dep%get_velocity()) + norm(vel2 - arr%get_velocity())
    end do
  end subroutine lambert_intermediate

  subroutine lambert_core(date, delta_t, pos1, pos2, axis, vel1, vel2, trf)
    implicit none
    real(8), intent(in) :: date, delta_t, pos1(3), pos2(3), axis(3)
    real(8), intent(out) :: vel1(3), vel2(3)
    type(TOrbit), intent(out), optional :: trf
    real(8) :: r1, r2, delta_nu, cos_dnu, sin_dnu, c, s, dtm, dtp
    real(8) :: a, am, t, tm, dt, ddtda, range
    real(8) :: alpha, beta, betam, gamma, delta, f, gi, gd
    real(8) :: cos_nu1, sin_nu1, p, e, ta1, ea1, ma1, axis_n(3)
    integer :: sector, type, round, roundn, i
    logical :: flag

    r1       = norm(pos1)
    r2       = norm(pos2)
    delta_nu = calc_angle(pos1, pos2, axis)
    sin_dnu  = sin(delta_nu)
    cos_dnu  = cos(delta_nu)
    sector   = Sgn(PI - delta_nu)

    c        = norm(pos2 - pos1)
    am       = 0.25d0 * (r1 + r2 + c)
    s        = 2d0 * am
    tm       = calc_period(am)
    betam    = 2d0 * asin(sqrt(1d0 - c / s))

    dtm      = tm * (0 + 0.5d0 - sector * (betam - sin(betam)) / PI2)
    dtp      = sqrt(2d0 / MU) * (s ** 1.5d0 - sector * (s - c) ** 1.5d0) / 3d0

    type     = Sgn(delta_t - dtp)
    round    = Sgn(delta_t - dtm)
    roundn   = (round + 1) / 2 + 0

    a        = 1.1d0 * am
    alpha    = 0d0
    beta     = 0d0
    gamma    = 0d0
    delta    = 0d0
    flag     = .true.
    range    = 0.5d0 * (a - am)

    do i = 1, 1000
      t = calc_period(a)

      if (type > 0) then
        alpha = 2d0 * asin(sqrt(s       / (2d0 * a)))
        beta  = 2d0 * asin(sqrt((s - c) / (2d0 * a)))
        dt    = t * (roundn - (round * (alpha - sin(alpha)) + sector * (beta - sin(beta))) / PI2)
        ! ddtda = 1.5d0 * (roundn * t + dt) / a &
        !       + (round * s ** 2 / sin(alpha) + sector * (s - c) ** 2 / sin(beta)) / sqrt(MU * a ** 3)
        if (round * (delta_t - dt) > 0) then
          if (flag) then
            a     = 2d0 * a
            range = 0.5d0 * a
          else
            a     = a + range
            range = 0.5d0 * range
          end if
        else
          if (flag) flag = .false.
          a     = a - range
          range = 0.5d0 * range
        end if
      else
        gamma = 2d0 * asinh(sqrt(s       / (2d0 * a)))
        delta = 2d0 * asinh(sqrt((s - c) / (2d0 * a)))
        dt    = t / PI2 * (sinh(gamma) - gamma - sector * (sinh(delta) - delta))
        ddtda = 1.5d0 * dt / a &
              - (s ** 2 / sinh(gamma) - sector * (s - c) ** 2 / sinh(delta)) / sqrt(MU * a ** 3)
        a     = max(a + (delta_t - dt) / ddtda, 0.1d0 * a)
      end if

      if (abs(1d0 - dt / delta_t) < 1d-5) exit
    end do

    if (type > 0) then
      p = 4d0 * a * (s - r1) * (s - r2) / c ** 2 * sin(0.5d0 * (alpha - round * sector * beta)) ** 2
      e = sqrt(1d0 - p / a)
    else
      p = 4d0 * a * (s - r1) * (s - r2) / c ** 2 * sinh(0.5d0 * (gamma + sector * delta)) ** 2
      e = sqrt(1d0 + p / a)
    end if

    cos_nu1 = (p - r1) / (e * r1)
    sin_nu1 = (cos_nu1 * cos_dnu - (p - r2) / (e * r2)) / sin_dnu

    ta1     = calc_angle(cos_nu1, sin_nu1)
    ea1     = convert_anomaly(ta1, e, from='T', to='E')
    ma1     = convert_anomaly(ea1, e, from='E', to='M')

    gi      = sqrt(MU * p) / (r1 * r2 * sin_dnu)
    f       = 1d0 - r2 * (1d0 - cos_dnu) / p
    gd      = 1d0 - r1 * (1d0 - cos_dnu) / p

    vel1    = gi * (pos2 - f * pos1)
    vel2    = gi * (gd * pos2 - pos1)

    if (present(trf)) then
      trf%type   = type
      trf%epc    = date
      trf%ecc    = e
      trf%sma    = a
      trf%slr    = p
      trf%mm     = PI2 / calc_period(a)
      trf%prd    = PI2 / trf%mm
      trf%axis_w = sector * Unit(cross_product(pos1, pos2))
      trf%axis_p = Unit(Rotate(pos1, trf%axis_w, -ta1))
      trf%axis_q = cross_product(trf%axis_w, trf%axis_p)
      axis_n     = Unit(cross_product(AXIS_K, trf%axis_w))
      trf%inc    = calc_angle(AXIS_K, trf%axis_w)
      trf%ran    = calc_angle(AXIS_I, axis_n, AXIS_K)
      trf%ap     = calc_angle(axis_n, trf%axis_p, trf%axis_w)
      trf%ta     = ta1
      trf%ea     = ea1
      trf%ma     = ma1
      trf%pos    = pos1
      trf%vel    = vel1
      trf%dran   = -1.5d0 * J2 * trf%mm / p ** 2 * cos(trf%inc)
      trf%dap    =  1.5d0 * J2 * trf%mm / p ** 2 * (2d0 - 2.5d0 * sin(trf%inc) ** 2)
    end if
  end subroutine lambert_core
end module orbit_base
