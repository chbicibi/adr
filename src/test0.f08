module test
  use Base_Func
  use IO_Output
  use Solver_Orbit
  use Solver_Debri

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
      theta   = Sgn(DOT_PRODUCT(Cross(orbit1%axis_w, orbit2%axis_w), axis)) * Angle(orbit1%axis_w, orbit2%axis_w)

      offset  = MODULO(TAtoMA(Angle(orbit1%axis_p, axis, orbit1%axis_w), orbit1%ecc) - orbit1%ma, PI2) / orbit1%mm

      date_d  = date + offset / 86400d0
      orbit_d = orbit1%CalcOrbit(date_d)
      pos_d   = orbit_d%GetPosition()
      vel_d   = orbit_d%GetVelocity()

      select case (i)
      case (1)
        sma_t   = Norm(pos_d) * (1d0 + rho) / 2d0 + 1d-10
      case (2)
        offset  = MODULO(TAtoMA(Angle(orbit2%axis_p, axis, orbit2%axis_w), orbit2%ecc) - orbit1%ma, PI2) / orbit1%mm
        date_a  = date + offset / 86400d0
        orbit_a = orbit2%CalcOrbit(date_a)
        pos_a   = orbit_a%GetPosition()
        ! sma_t   = (Norm(pos_d) + Norm(pos_a)) / 2d0 + 1d-10
        sma_t   = (Norm(pos_d) + orbit_a%sma + 30d0) / 2d0
      case (3)
        sma_t   = Norm(pos_d) + 1d-10
      end select
      print *, "sma", sma
      axis_t    = Rotate(orbit_d%axis_w, axis, theta * t(i))
      orbit_t   = GetOrbit(pos_d, sma_t, axis_t, date_d)

      vel_t     = orbit_t%GetVelocity()
      del_v     = del_v + Norm(vel_t - vel_d)

      if (op == 1) then
        print "(3f10.5)", vel_d
        print "(3f10.5)", vel_t
        print "(3f10.5)", vel_t - vel_d
        print "(f10.5/)", Norm(vel_t - vel_d)
        if (i < 2) call OutputOrbit(orbit_d, date, offset)
        if (i < 3) call OutputOrbit(orbit_t, date_d, orbit_t%period / 2d0)
      end if

      date      = date_d + orbit_t%period * n(i) / 86400d0
      orbit1    = orbit_t%CalcOrbit(date)
      orbit2    = orbit2%CalcOrbit(date)
    end do

    tof      = CalcPeriod((orbit1%sma + orbit2%sma) / 2d0) * dnu
    orbit1_a = orbit1%CalcOrbit(date + tof / 86400d0)
    orbit2_a = orbit2%CalcOrbit(date + tof / 86400d0)
    offset   = MODULO(Sgn(orbit2_a%mm - orbit1_a%mm) * (TAtoMA(Angle(orbit2_a%axis_p, orbit1_a%GetPosition(), orbit2_a%axis_w), &
                orbit2_a%ecc) - orbit2_a%ma), PI2) / ABS(orbit2_a%mm - orbit1_a%mm)
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
      call OutputOrbit(orbit1_a, date, MODULO(offset, orbit1_a%period))
      call OutputOrbit(orbit1, date_d, tof)
      call OutputOrbit(orbit2, date_d, tof)
    end if
    ! print "(2f15.5)", Norm(vel_a - v2) + Norm(vel_d - v1), offset / 86400d0
    del_v = Norm(vel_a - v2) + Norm(vel_d - v1)
    del_t = offset / 86400d0
    ! print *, Angle(pos_d, pos_a) * 360d0 / PI2
    ! print *, Angle(orbit_d%axis_w, orbit_a%axis_w) * 360d0 / PI2
    ! print *, Angle(Cross(pos_a, pos_d), orbit_a%axis_w) * 360d0 / PI2
    ! print "(3f15.5)", Cross(Unit(pos_d), Unit(pos_a))

  end subroutine Transfer0

  subroutine Test_
    type(TOrbit) :: orbit
    real(8)      :: raand, apd
    integer      :: i

    orbit = DEBRIS(1)%orbit

    do i = 1, 100
      raand = -3 * PI * J2 / (orbit%slr ** 2 * orbit%period) * COS(orbit%incl)
      apd   = -3 * PI * J2 / (orbit%slr ** 2 * orbit%period) * (2d0 - 2.5d0 * SIN(orbit%incl) ** 2)
      print *, raand * 86400d0 * 360d0 / PI2
    end do
  end subroutine Test_

end module test
