program XTest
  use Base
  use Base_Func
  use Solver_Orbit
  use Solver_Debri
  use IO_Output
  ! use IO_Path

  implicit none

  integer      :: i, j, n, imax, jmax

  real(8)      :: date, delt, st
  real(8)      :: pos1(3), pos2(3), axis(3)
  real(8)      :: vel1(3), vel2(3)
  type(TOrbit) :: orbit, orbit1, orbit2

  real(8)      :: del_v, del_t, v, t, v2, t2, dt

  real(8)      :: out(100)

  call Initialize
  call InitializeFile

  ! dt   = 0.001d0
  date = MJD(2015, 0d0)

  orbit1 = GetOrbit(date, 0d0, 90d0, 0.01d0, 0d0, 0d0, 86400d0 / CalcPeriod(6378 + 4000d0))
  orbit2 = GetOrbit(date, 0d0, 90d0, 0.3d0, 60d0, 45d0, 86400d0 / CalcPeriod(6378 + 4000d0))

  call orbit1%Show
  call orbit2%Show

  ! st   = date + dt * dble(0)
  ! delt =        dt * dble(i + 1) / 50

  ! call orbit1%Move(st, 0)
  ! call orbit2%Move(st + delt, 0)
  ! pos1 = orbit1%GetPosition()
  ! pos2 = orbit2%GetPosition()
  ! call CalcLambert(st, delt * 86400d0, pos1, pos2, axis, vel1, vel2, orbit)
  ! del_v = Norm(vel1 - orbit1%GetVelocity()) + Norm(vel2 - orbit2%GetVelocity())

  ! print *, i, del_v

  call orbit1%Output(date, orbit1%prd)
  call orbit2%Output(date, orbit2%prd)

  call OutputPoint(orbit1%GetPosition())
  call OutputPoint(orbit2%GetPosition())
  call OutputArrow([0d0, 0d0, 0d0, orbit1%GetPosition()])
  call OutputArrow([0d0, 0d0, 0d0, orbit2%GetPosition()])
  call orbit2%Move(-orbit2%ma, 2)
  call OutputArrow([0d0, 0d0, 0d0, orbit2%GetPosition()])
  call orbit1%Move(-PI2 / 4, 2)
  call OutputArrow([0d0, 0d0, 0d0, orbit1%GetPosition()])
  call OutputArrow([0d0, 0d0, 0d0, Unit(Cross(orbit2%axis_w, orbit2%GetPosition())) * Norm(orbit2%GetPosition())])
  call OutputArrow([0d0, 0d0, 0d0, orbit2%axis_w * Norm(orbit2%GetPosition())])
  ! call OutputOrbit(orbit,  st, delt * 86400)
  ! call unlink("plot\orbit" // Str(i) // ".dat")
  ! call rename(FILEPATH_ORBIT, "plot\orbit" // Str(i) // ".dat")

  call Terminate

  contains
end program XTest
