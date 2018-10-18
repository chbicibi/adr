module XMod
  use Base_Func
  use Solver_Orbit
  use Solver_Debri
  use IO_Output

  implicit none

  contains

  subroutine Transfer0(debri_num, args, op, del_v, del_t)
    integer,     intent(in)  :: debri_num(2), op
    real(8),     intent(in)  :: args(7)
    real(8),     intent(out) :: del_v, del_t
    type(TOrbit)             :: orbit1, orbit2, orbit_d, orbit_t, orbit_a, orbit1_a, orbit2_a
    real(8)                  :: date_s, rho, dnu, t(3), n(3), date, date_d, date_a, theta, axis(3), axis_t(3), offset
    real(8)                  :: sma_t, tof, pos_d(3), pos_a(3), vel_d(3), vel_t(3), vel_a(3), v1(3), v2(3)
    integer                  :: i

    date_s = args(1)
    rho    = args(2)
    dnu    = args(3)
    t      = [args(4:5), 1d0]
    n      = [args(6:7), 0d0]

    del_v  = 0d0
    del_t  = 0d0

    if (debri_num(1) == debri_num(2)) return

    date     = date_s
    orbit1   = DEBRIS(debri_num(1))%orbit%CalcOrbit(date)
    orbit2   = DEBRIS(debri_num(2))%orbit%CalcOrbit(date)

    do i = 1, 3
      axis    = Unit(Cross(orbit1%axis_w, orbit2%axis_w)) * (-1) ** (i - 1)
      theta   = Sgn(dot_product(Cross(orbit1%axis_w, orbit2%axis_w), axis)) * Angle(orbit1%axis_w, orbit2%axis_w)

      offset  = modulo(TAtoMA(Angle(orbit1%axis_p, axis, orbit1%axis_w), orbit1%ecc) - orbit1%ma, PI2) / orbit1%mm

      date_d  = date + offset / 86400d0
      orbit_d = orbit1%CalcOrbit(date_d)
      pos_d   = orbit_d%GetPosition()
      vel_d   = orbit_d%GetVelocity()

      select case (i)
      case (1)
        sma_t   = Norm(pos_d) * (1d0 + rho) / 2d0 + 1d-10
      case (2)
        offset  = modulo(TAtoMA(Angle(orbit2%axis_p, axis, orbit2%axis_w), orbit2%ecc) - orbit1%ma, PI2) / orbit1%mm
        date_a  = date + offset / 86400d0
        orbit_a = orbit2%CalcOrbit(date_a)
        pos_a   = orbit_a%GetPosition()
        ! sma_t   = (Norm(pos_d) + Norm(pos_a)) / 2d0 + 1d-10
        sma_t   = (Norm(pos_d) + orbit_a%sma + 0d0) / 2d0
      case (3)
        sma_t   = Norm(pos_d) + 1d-10
      end select

      axis_t    = Rotate(orbit_d%axis_w, axis, theta * t(i))
      orbit_t   = GeTOrbit(pos_d, sma_t, axis_t, date_d)

      vel_t     = orbit_t%GetVelocity()

      ! print *, orbit_t%GetVelocity()
      ! print *, Norm(orbit_t%GetVelocity())
      ! print *, Norm(vel_d)
      ! print *, Angle(vel_t, vel_d)

      del_v     = del_v + Norm(vel_t - vel_d)

      if (op == 1) then
        print "(3f10.5)", vel_d
        print "(3f10.5)", vel_t
        print "(3f10.5)", vel_t - vel_d
        print "(f10.5/)", Norm(vel_t - vel_d)
        ! if (i < 2) call OutputOrbit(orbit_d, date, offset)
        ! if (i < 3) call OutputOrbit(orbit_t, date_d, orbit_t%prd / 2d0)
      end if

      date      = date_d + orbit_t%prd * n(i) / 86400d0
      orbit1    = orbit_t%CalcOrbit(date)
      orbit2    = orbit2%CalcOrbit(date)
    end do

    del_t    = date_d - date_s
    return

    tof      = CalcPeriod((orbit1%sma + orbit2%sma) / 2d0) * dnu
    orbit1_a = orbit1%CalcOrbit(date + tof / 86400d0)
    orbit2_a = orbit2%CalcOrbit(date + tof / 86400d0)
    offset   = modulo(Sgn(orbit2_a%mm - orbit1_a%mm) * (TAtoMA(Angle(orbit2_a%axis_p, orbit1_a%GetPosition(), orbit2_a%axis_w), &
                orbit2_a%ecc) - orbit2_a%ma), PI2) / abs(orbit2_a%mm - orbit1_a%mm)
    date_d   = date + offset / 86400d0 + 1d-1
    orbit_d  = orbit1%CalcOrbit(date_d)
    pos_d    = orbit_d%GetPosition()
    vel_d    = orbit_d%GetVelocity()

    orbit_a  = orbit2%CalcOrbit(date_d + tof / 86400d0)
    pos_a    = orbit_a%GetPosition()
    vel_a    = orbit_a%GetVelocity()

    call CalcLambert(date_d, tof, pos_d, pos_a, orbit_a%axis_w, v1, v2, orbit1)
    del_v    = del_v + Norm(v1 - vel_d)
    del_v    = del_v + Norm(vel_a - v2)

    del_t    = date_d - date_s + tof / 86400d0

    if (op == 1) then
      print "(3f10.5)", v1
      print "(3f10.5)", vel_d
      print "(3f10.5)", vel_d - v1
      print "(f10.5/)", Norm(vel_d - v1)
      print "(3f10.5)", vel_t
      print "(3f10.5)", v2
      print "(3f10.5)", v2 - vel_a
      print "(f10.5/)", Norm(vel_a - v2)
      ! call OutputOrbit(orbit1_a, date, modulo(offset, orbit1_a%prd))
      ! call OutputOrbit(orbit1, date_d, tof)
      ! call OutputOrbit(orbit2, date_d, tof)
    end if
    ! print "(2f15.5)", Norm(vel_a - v2) + Norm(vel_d - v1), offset / 86400d0
    ! del_v = Norm(vel_a - v2) + Norm(vel_d - v1)
    ! del_t = offset / 86400d0
    ! print *, Angle(pos_d, pos_a) * 360d0 / PI2
    ! print *, Angle(orbit_d%axis_w, orbit_a%axis_w) * 360d0 / PI2
    ! print *, Angle(Cross(pos_a, pos_d), orbit_a%axis_w) * 360d0 / PI2
    ! print "(3f15.5)", Cross(Unit(pos_d), Unit(pos_a))
  end subroutine Transfer0

  subroutine Transfer2(orbit1, orbit2, args, op, del_v, del_t)
    real(8),     intent(in)  :: args(7)
    type(TOrbit),intent(inout)  :: orbit1, orbit2
    integer,     intent(in)  :: op
    real(8),     intent(out) :: del_v, del_t
    type(TOrbit)             :: orbit_d, orbit_t, orbit_a, orbit1_a, orbit2_a
    real(8)                  :: date_s, rho, dnu, t(3), n(3), date, date_d, date_a, theta, axis(3), axis_t(3), offset
    real(8)                  :: sma_t, tof, pos_d(3), pos_a(3), vel_d(3), vel_t(3), vel_a(3), v1(3), v2(3)
    integer                  :: i

    date_s = args(1)
    rho    = args(2)
    dnu    = args(3)
    t      = [args(4:5), 1d0]
    n      = [args(6:7), 0d0]

    del_v  = 0d0
    del_t  = 0d0

    date     = date_s
    ! orbit1   = DEBRIS(debri_num(1))%orbit%CalcOrbit(date)
    ! orbit2   = DEBRIS(debri_num(2))%orbit%CalcOrbit(date)

    do i = 1, 3
      axis    = Unit(Cross(orbit1%axis_w, orbit2%axis_w)) * (-1) ** (i - 1)
      theta   = Sgn(dot_product(Cross(orbit1%axis_w, orbit2%axis_w), axis)) * Angle(orbit1%axis_w, orbit2%axis_w)

      offset  = modulo(TAtoMA(Angle(orbit1%axis_p, axis, orbit1%axis_w), orbit1%ecc) - orbit1%ma, PI2) / orbit1%mm

      date_d  = date + offset / 86400d0
      orbit_d = orbit1%CalcOrbit(date_d)
      pos_d   = orbit_d%GetPosition()
      vel_d   = orbit_d%GetVelocity()

      select case (i)
      case (1)
        sma_t   = Norm(pos_d) * (1d0 + rho) / 2d0 + 1d-10
      case (2)
        offset  = modulo(TAtoMA(Angle(orbit2%axis_p, axis, orbit2%axis_w), orbit2%ecc) - orbit1%ma, PI2) / orbit1%mm
        date_a  = date + offset / 86400d0
        orbit_a = orbit2%CalcOrbit(date_a)
        pos_a   = orbit_a%GetPosition()
        ! sma_t   = (Norm(pos_d) + Norm(pos_a)) / 2d0 + 1d-10
        sma_t   = (Norm(pos_d) + orbit_a%sma + 30d0) / 2d0
      case (3)
        sma_t   = Norm(pos_d) + 1d-10
      end select
      axis_t    = Rotate(orbit_d%axis_w, axis, theta * t(i))
      orbit_t   = GeTOrbit(pos_d, sma_t, axis_t, date_d)

      vel_t     = orbit_t%GetVelocity()
      del_v     = del_v + Norm(vel_t - vel_d)

      if (op == 1) then
        print "(3f10.5)", vel_d
        print "(3f10.5)", vel_t
        print "(3f10.5)", vel_t - vel_d
        print "(f10.5/)", Norm(vel_t - vel_d)
        ! if (i < 2) call OutputOrbit(orbit_d, date, offset)
        ! if (i < 3) call OutputOrbit(orbit_t, date_d, orbit_t%prd / 2d0)
      end if

      date      = date_d + orbit_t%prd * n(i) / 86400d0
      orbit1    = orbit_t%CalcOrbit(date)
      orbit2    = orbit2%CalcOrbit(date)
    end do

    tof      = CalcPeriod((orbit1%sma + orbit2%sma) / 2d0) * dnu
    orbit1_a = orbit1%CalcOrbit(date + tof / 86400d0)
    orbit2_a = orbit2%CalcOrbit(date + tof / 86400d0)
    offset   = modulo(Sgn(orbit2_a%mm - orbit1_a%mm) * (TAtoMA(Angle(orbit2_a%axis_p, orbit1_a%GetPosition(), orbit2_a%axis_w), &
                orbit2_a%ecc) - orbit2_a%ma), PI2) / abs(orbit2_a%mm - orbit1_a%mm)
    date_d   = date + offset / 86400d0 + 1d-1
    orbit_d  = orbit1%CalcOrbit(date_d)
    pos_d    = orbit_d%GetPosition()
    vel_d    = orbit_d%GetVelocity()

    orbit_a  = orbit2%CalcOrbit(date_d + tof / 86400d0)
    pos_a    = orbit_a%GetPosition()
    vel_a    = orbit_a%GetVelocity()

    call CalcLambert(date_d, tof, pos_d, pos_a, orbit_a%axis_w, v1, v2, orbit1)
    del_v    = del_v + Norm(v1 - vel_d)
    del_v    = del_v + Norm(vel_a - v2)

    del_t    = date_d - date_s + tof / 86400d0

    if (op == 1) then
      print "(3f10.5)", v1
      print "(3f10.5)", vel_d
      print "(3f10.5)", vel_d - v1
      print "(f10.5/)", Norm(vel_d - v1)
      print "(3f10.5)", vel_t
      print "(3f10.5)", v2
      print "(3f10.5)", v2 - vel_a
      print "(f10.5/)", Norm(vel_a - v2)
      ! call OutputOrbit(orbit1_a, date, modulo(offset, orbit1_a%prd))
      ! call OutputOrbit(orbit1, date_d, tof)
      ! call OutputOrbit(orbit2, date_d, tof)
    end if
    ! print "(2f15.5)", Norm(vel_a - v2) + Norm(vel_d - v1), offset / 86400d0
    del_v = Norm(vel_a - v2) + Norm(vel_d - v1)
    del_t = offset / 86400d0
    ! print *, Angle(pos_d, pos_a) * 360d0 / PI2
    ! print *, Angle(orbit_d%axis_w, orbit_a%axis_w) * 360d0 / PI2
    ! print *, Angle(Cross(pos_a, pos_d), orbit_a%axis_w) * 360d0 / PI2
    ! print "(3f15.5)", Cross(Unit(pos_d), Unit(pos_a))
  end subroutine Transfer2

  subroutine CalcLambert2(date, pos1, pos2, axis, vel1, vel2, dt, orbit)
    real(8),         intent(in)  :: date, pos1(3), pos2(3), axis(3)
    real(8),         intent(out) :: vel1(3), vel2(3), dt
    type(TOrbit), intent(out) :: orbit
    real(8)                      :: r1, r2, delta_nu, cos_dnu, sin_dnu, c, s, dtm, dtp
    real(8)                      :: a, am, t, tm, ddtda, range
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

    type     = 1
    round    = 1
    roundn   = (round + 1) / 2 + 0

    a        = 1.0d0 * am
    alpha    = 0d0
    beta     = 0d0
    gamma    = 0d0
    delta    = 0d0
    flag     = .true.

    t = CalcPeriod(a)
    alpha = 2d0 * asin(sqrt(s       / (2d0 * a)))
    beta  = 2d0 * asin(sqrt((s - c) / (2d0 * a)))
    dt    = t * (roundn - (round * (alpha - sin(alpha)) + sector * (beta - sin(beta))) / PI2)

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
    ! ea1     = TAtoEA(ta1, e)
    ! ma1     = EAtoMA(ea1, e)

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
    orbit%prd = PI2 / orbit%mm
    orbit%axis_w = sector * Unit(Cross(pos1, pos2))
    orbit%axis_p = Unit(Rotate(pos1, orbit%axis_w, -ta1))
    orbit%axis_q = Cross(orbit%axis_w, orbit%axis_p)
    axis_n       = Unit(Cross(axis_k, orbit%axis_w))
    orbit%inc   = Angle(axis_k, orbit%axis_w)
    orbit%ran   = Angle(axis_i, axis_n, axis_k)
    orbit%ap     = Angle(axis_n, orbit%axis_p, orbit%axis_w)
    orbit%ta     = ta1
    ! orbit%ea     = ea1
    ! orbit%ma     = ma1
    orbit%pos    = pos1
    orbit%vel    = vel1

    orbit%dran = -1.5d0 * J2 * orbit%mm / p ** 2 * cos(orbit%inc)
    orbit%dap   =  1.5d0 * J2 * orbit%mm / p ** 2 * (2d0 - 2.5d0 * sin(orbit%inc) ** 2)
  end subroutine CalcLambert2

end module XMod
