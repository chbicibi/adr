module Solver_Orbit
  use Base_Func

  implicit none

  private
  public :: TOrbit
  public :: MU, J2, axis_i, axis_j, axis_k
  public :: GetOrbit, CalcLambert, CalcPeriod, MAtoEA, TAtoMA

  type :: TOrbit
    real(8) :: epoch     = -1d0
    real(8) :: incl      = -1d0
    real(8) :: raan      = -1d0
    real(8) :: ecc       = -1d0
    real(8) :: ap        = -1d0
    real(8) :: ma        = -1d0
    real(8) :: mm        = -1d0

    integer :: type      = 1
    real(8) :: sma       = -1d0
    real(8) :: slr       = -1d0
    real(8) :: period    = -1d0
    real(8) :: axis_p(3) = 0d0
    real(8) :: axis_q(3) = 0d0
    real(8) :: axis_w(3) = 0d0
    real(8) :: amom(3)   = 0d0 ! angular momentum, h

    real(8) :: dt        = -1d0
    real(8) :: ea        = -1d0 ! eccentric Anomaly, E /
    real(8) :: ta        = -1d0 ! true anomaly, f

    real(8) :: pos(3)    = 0d0 ! position, r
    real(8) :: vel(3)    = 0d0 ! velocity, v

    real(8) :: draan     = 0d0
    real(8) :: dap       = 0d0

    contains
    procedure :: CalcOrbit
    procedure :: GetPosition
    procedure :: GetVelocity
  end type TOrbit

  interface GetOrbit
    procedure :: CreateOrbit, PosToOrbit, PeriToOrbit
  end interface GetOrbit

  real(8), parameter :: MU        = 3.9860044d5
  real(8), parameter :: J2        = 4.4041964d4! == J2 * re ** 2
  real(8), parameter :: axis_i(3) = [1d0, 0d0, 0d0]
  real(8), parameter :: axis_j(3) = [0d0, 1d0, 0d0]
  real(8), parameter :: axis_k(3) = [0d0, 0d0, 1d0]

  contains
  subroutine CalcLambert(date, delta_t, pos1, pos2, axis, vel1, vel2, orbit)
    real(8),         intent(in)  :: date, delta_t, pos1(3), pos2(3), axis(3)
    real(8),         intent(out) :: vel1(3), vel2(3)
    type(TOrbit), intent(out) :: orbit
    real(8)                      :: r1, r2, delta_nu, cos_dnu, sin_dnu, c, s, dtm, dtp
    real(8)                      :: a, am, t, tm, dt, ddtda, range
    real(8)                      :: alpha, beta, betam, gamma, delta, f, gi, gd
    real(8)                      :: cos_nu1, sin_nu1, p, e, ta1, ea1, ma1, axis_n(3)
    integer                      :: sector, type, round, roundn, i
    logical                      :: flag

    r1       = Norm(pos1)
    r2       = Norm(pos2)
    delta_nu = Angle(pos1, pos2, axis)
    sin_dnu  = sin(delta_nu)
    cos_dnu  = cos(delta_nu)
    sector   = Sgn(PI - delta_nu)

    c        = Norm(pos2 - pos1)
    am       = 0.25d0 * (r1 + r2 + c)
    s        = 2d0 * am
    tm       = CalcPeriod(am)
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
    range    = (a - am) / 2d0

    do i = 1, 1000
      t = CalcPeriod(a)

      if (type > 0) then
        alpha = 2d0 * asin(sqrt(s       / (2d0 * a)))
        beta  = 2d0 * asin(sqrt((s - c) / (2d0 * a)))
        dt    = t * (roundn - (round * (alpha - sin(alpha)) + sector * (beta - sin(beta))) / PI2)
        ! ddtda = 1.5d0 * (roundn * t + dt) / a &
        !       + (round * s ** 2 / sin(alpha) + sector * (s - c) ** 2 / sin(beta)) / sqrt(MU * a ** 3)
        if (round * (delta_t - dt) > 0) then
          if (flag) then
            a     = a * 2d0
            range = a / 2d0
          else
            a     = a + range
            range = range / 2d0
          end if
        else
          if (flag) flag = .false.
          a     = a - range
          range = range / 2d0
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
      p = 4d0 * a * (s - r1) * (s - r2) / c ** 2 * sin((alpha - round * sector * beta) / 2d0) ** 2
      e = sqrt(1d0 - p / a)
    else
      p = 4d0 * a * (s - r1) * (s - r2) / c ** 2 * sinh((gamma + sector * delta) / 2d0) ** 2
      e = sqrt(1d0 + p / a)
    end if

    cos_nu1 = (p - r1) / (e * r1)
    sin_nu1 = (cos_nu1 * cos_dnu - (p - r2) / (e * r2)) / sin_dnu

    ta1     = Angle(cos_nu1, sin_nu1)
    ea1     = TAtoEA(ta1, e)
    ma1     = EAtoMA(ea1, e)

    gi      = sqrt(MU * p) / (r1 * r2 * sin_dnu)
    f       = 1d0 - r2 * (1d0 - cos_dnu) / p
    gd      = 1d0 - r1 * (1d0 - cos_dnu) / p

    vel1    = gi * (pos2 - f * pos1)
    vel2    = gi * (gd * pos2 - pos1)

    orbit%type   = type
    orbit%epoch  = date
    orbit%ecc    = e
    orbit%sma    = a
    orbit%slr    = p
    orbit%mm     = PI2 / CalcPeriod(a)
    orbit%period = PI2 / orbit%mm
    orbit%axis_w = sector * Unit(Cross(pos1, pos2))
    orbit%axis_p = Unit(Rotate(pos1, orbit%axis_w, -ta1))
    orbit%axis_q = Cross(orbit%axis_w, orbit%axis_p)
    axis_n       = Unit(Cross(axis_k, orbit%axis_w))
    orbit%incl   = Angle(axis_k, orbit%axis_w)
    orbit%raan   = Angle(axis_i, axis_n, axis_k)
    orbit%ap     = Angle(axis_n, orbit%axis_p, orbit%axis_w)
    orbit%ta     = ta1
    orbit%ea     = ea1
    orbit%ma     = ma1
    orbit%pos    = pos1
    orbit%vel    = vel1

    orbit%draan = -1.5d0 * J2 * orbit%mm / p ** 2 * cos(orbit%incl)
    orbit%dap   =  1.5d0 * J2 * orbit%mm / p ** 2 * (2d0 - 2.5d0 * sin(orbit%incl) ** 2)
  end subroutine CalcLambert

  !-----------------------------------------------------------------------------
  ! Constructor
  !-----------------------------------------------------------------------------
  type(TOrbit) function CreateOrbit(epoch, incl, raan, ecc, ap, ma, mm) result (orbit)
    real(8), intent(in) :: epoch, incl, raan, ecc, ap, ma, mm
    real(8)             :: matrix(3, 3)

    orbit%epoch  = epoch
    orbit%incl   = incl * PI2 / 360d0
    orbit%raan   = raan * PI2 / 360d0
    orbit%ecc    = ecc
    orbit%ap     = ap   * PI2 / 360d0
    orbit%ma     = ma   * PI2 / 360d0
    orbit%mm     = mm   * PI2 / 86400d0
    orbit%sma    = (MU / orbit%mm ** 2) ** (1d0 / 3d0)
    orbit%slr    = orbit%sma * (1d0 - ecc ** 2)
    orbit%period = 86400d0 / mm

    matrix       = EulerAngle([orbit%ap, orbit%incl, orbit%raan], [3, 1, 3])
    orbit%axis_p = matrix(1, :)
    orbit%axis_q = matrix(2, :)
    orbit%axis_w = matrix(3, :)

    orbit%draan  = -1.5d0 * J2 * orbit%mm / orbit%slr ** 2 * cos(orbit%incl)
    orbit%dap    =  1.5d0 * J2 * orbit%mm / orbit%slr ** 2 * (2d0 - 2.5d0 * sin(orbit%incl) ** 2)
  end function CreateOrbit

  type(TOrbit) function PosToOrbit(pos, vel, date) result (orbit)
    real(8), intent(in) :: pos(3), vel(3), date
    real(8)             :: r, sma, ecc, esin, ecos, ea
    real(8)             :: h(3), axis_p(3), axis_q(3), axis_w(3), axis_n(3)
    integer             :: type

    r      = Norm(pos)
    sma    = 1d0 / (2d0 / r - dot_product(vel, vel) / MU)
    h      = Cross(pos, vel)
    axis_p = Unit(-MU * pos / r - Cross(h, vel))
    axis_w = Unit(h)
    axis_q = Cross(axis_w, axis_p)
    axis_n = Unit(Cross(axis_k, axis_w))
    ecc    = sqrt(1d0 - dot_product(h, h) / (MU * sma))
    type   = Sgn(1d0 - ecc)
    esin   = dot_product(pos, vel) / sqrt(mu * sma)

    if (ecc > 0) then
      if (type > 0) then
        ecos = 1d0 - r / sma
        ea   = Angle(ecos / ecc, esin)
      else
        ea   = asinh(esin / ecc)
      end if
    else
      ea     = 0d0
    end if

    orbit%type   = type
    orbit%epoch  = date
    orbit%incl   = Angle(axis_k, axis_w)
    orbit%raan   = Angle(axis_i, axis_n, axis_k)
    orbit%ecc    = ecc
    orbit%ap     = Angle(axis_n, axis_p, axis_w)
    orbit%ma     = type * (ea - esin)
    orbit%mm     = PI2 / CalcPeriod(sma)
    orbit%sma    = sma
    orbit%slr    = sma * (1d0 - ecc ** 2)
    orbit%period = PI2 / orbit%mm
    orbit%axis_p = axis_p
    orbit%axis_q = axis_q
    orbit%axis_w = axis_w
    ! orbit%amom   = h
    ! orbit%dt     = orbit%ma / orbit%mm
    ! orbit%ea     = ea
    orbit%pos    = pos
    orbit%vel    = vel

    orbit%draan = -1.5d0 * J2 * orbit%mm / orbit%slr ** 2 * cos(orbit%incl)
    orbit%dap   =  1.5d0 * J2 * orbit%mm / orbit%slr ** 2 * (2d0 - 2.5d0 * sin(orbit%incl) ** 2)
  end function PosToOrbit

  type(TOrbit) function PeriToOrbit(pos, sma, axis_w, date) result (orbit)
    real(8), intent(in) :: pos(3), sma, axis_w(3), date
    real(8)             :: r, ecc, axis_p(3), axis_q(3), axis_n(3)
    integer             :: dir

    r      = Norm(pos)
    dir    = Sgn(sma - r)
    axis_p = dir * Unit(pos)
    axis_q = Cross(axis_w, axis_p)
    axis_n = Unit(Cross(axis_k, axis_w))
    ecc    = dir * (1d0 - r / sma)

    orbit%type   = Sgn(1d0 - ecc)
    orbit%epoch  = date
    orbit%ecc    = ecc
    orbit%ma     = PI * (1 - dir) / 2
    orbit%mm     = PI2 / CalcPeriod(sma)
    orbit%sma    = sma
    orbit%slr    = sma * (1d0 - ecc ** 2)
    orbit%period = PI2 / orbit%mm
    orbit%axis_p = axis_p
    orbit%axis_q = axis_q
    orbit%axis_w = axis_w
    orbit%incl   = Angle(axis_k, axis_w)
    orbit%raan   = Angle(axis_i, axis_n, axis_k)
    orbit%ap     = Angle(axis_n, axis_p, axis_w)

    orbit%draan = -1.5d0 * J2 * orbit%mm / orbit%slr ** 2 * cos(orbit%incl)
    orbit%dap   =  1.5d0 * J2 * orbit%mm / orbit%slr ** 2 * (2d0 - 2.5d0 * sin(orbit%incl) ** 2)
  end function PeriToOrbit

  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  type(TOrbit) function CalcOrbit(self, date) result (orbit)
    class(TOrbit), intent(in) :: self
    real(8),          intent(in) :: date
    real(8)                      :: offset, matrix(3, 3)

    orbit%type   = self%type
    orbit%incl   = self%incl
    orbit%ecc    = self%ecc
    orbit%mm     = self%mm
    orbit%sma    = self%sma
    orbit%slr    = self%slr
    orbit%period = self%period
    orbit%draan  = self%draan
    orbit%dap    = self%dap

    offset       = 86400d0 * (date - self%epoch)
    orbit%epoch  = date
    orbit%raan   = modulo(self%raan + self%draan * offset, PI2)
    orbit%ap     = modulo(self%ap   + self%dap   * offset, PI2)
    orbit%ma     = modulo(self%ma   + self%mm    * offset, PI2)
    orbit%ea     = MAtoEA(orbit%ma, self%ecc)
    matrix       = EulerAngle([orbit%ap, orbit%incl, orbit%raan], [3, 1, 3])
    orbit%axis_p = matrix(1, :)
    orbit%axis_q = matrix(2, :)
    orbit%axis_w = matrix(3, :)
    ! orbit%dt    = ma / self%mm
  end function CalcOrbit

  function GetPosition(self) result (pos)
    real(8)                         :: pos(3)
    class(TOrbit), intent(inout) :: self
    real(8)                         :: s, c

    if (self%ea < 0) self%ea = MAtoEA(self%ma, self%ecc)

    if (self%type > 0) then
      s = sin(self%ea)
      c = cos(self%ea)
    else
      s = sinh(self%ea)
      c = cosh(self%ea)
    end if
    pos = self%sma * (self%type * (c - self%ecc) * self%axis_p &
                 + sqrt(1d0 - self%ecc ** 2) * s * self%axis_q)
  end function GetPosition

  function GetVelocity(self) result (vel)
    real(8)                         :: vel(3)
    class(TOrbit), intent(inout) :: self
    real(8)                         :: r, s, c

    if (self%ea < 0) self%ea = MAtoEA(self%ma, self%ecc)

    if (ALL(self%pos == 0d0)) then
      r = Norm(GetPosition(self))
    else
      r = Norm(self%pos)
    end if
    if (self%type > 0) then
      s = sin(self%ea)
      c = cos(self%ea)
    else
      s = sinh(self%ea)
      c = cosh(self%ea)
    end if
    vel = -sqrt(MU) / r * sqrt(self%sma) * (s * self%axis_p &
             -  sqrt(1d0 - self%ecc ** 2) * c * self%axis_q)
  end function GetVelocity

  !-----------------------------------------------------------------------------
  ! convert angle
  !-----------------------------------------------------------------------------
  real(8) function MAtoEA(ma, ecc)
    real(8), intent(in) :: ma, ecc
    real(8)             :: z, z_next, f, dfdz
    integer             :: i

    z = ma
    do i = 1, 1000
      if (ecc < 1) then
        f    = z   - ecc * sin(z) - ma
        dfdz = 1d0 - ecc * cos(z)
      else
        f    = ecc * sinh(z) - z - ma
        dfdz = ecc * cosh(z) - 1d0
      end if
      z_next = z - (f / dfdz)

      if (abs(z_next - z) < 1d-6) exit
      z = z_next
    end do
    MAtoEA = z_next
  end function MAtoEA

  real(8) function MAtoTA(ea, ecc)
    real(8), intent(in) :: ea, ecc

    MAtoTA = EAtoMA(MAtoEA(ea, ecc), ecc)
  end function MAtoTA

  real(8) function EAtoMA(ea, ecc)
    real(8), intent(in) :: ea, ecc

    EAtoMA = ea - ecc * sin(ea)
  end function EAtoMA

  real(8) function EAtoTA(ea, ecc)
    real(8), intent(in) :: ea, ecc
    real(8)             :: c

    c      = cos(ea)
    EAtoTA = Angle((c - ecc) / (1d0 - ecc * c), PI - ea)
  end function EAtoTA

  real(8) function TAtoEA(ta, ecc)
    real(8), intent(in) :: ta, ecc
    real(8)             :: c

    c      = cos(ta)
    TAtoEA = Angle((c + ecc) / (1d0 + ecc * c), PI - ta)
  end function TAtoEA

  real(8) function TAtoMA(ta, ecc)
    real(8), intent(in) :: ta, ecc

    TAtoMA = EAtoMA(TAtoEA(ta, ecc), ecc)
  end function TAtoMA

  real(8) function CalcPeriod(sma)
    real(8), intent(in) :: sma

    CalcPeriod = PI2 * sqrt(sma ** 3 / MU)
  end function CalcPeriod
end module Solver_Orbit
