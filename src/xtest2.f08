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
  type(TOrbit) :: orbit, orbit1, orbit2, orbit3, orbit4, orbit5

  real(8)      :: del_v, del_t, v, t, v2, t2, dt

  real(8)      :: out(100)

  call Initialize

  ! dt   = 0.001d0
  date = MJD(2015, 0d0)

  ! 0-2
  orbit1 = GetOrbit(date, 100d0, 0d0, 10d-6, 0d0, 0d0, 86400d0 / CalcPeriod(6378 + 1000d0))
  orbit2 = GetOrbit(date, 100d0, 0d0, 10d-6, 0d0, 0d0, 86400d0 / CalcPeriod(6378 + 3000d0))
  orbit3 = GetOrbit(date, 100d0, 0d0, 10d-6, 0d0, 0d0, 86400d0 / CalcPeriod(6378 + 5000d0))
  orbit4 = GetOrbit(date, 100d0, 0d0, 10d-6, 0d0, 0d0, 86400d0 / CalcPeriod(6378 + 7000d0))
  orbit5 = GetOrbit(date, 100d0, 0d0, 10d-6, 0d0, 0d0, 86400d0 / CalcPeriod(6378 + 9000d0))
  ! 3-5
  orbit1 = GetOrbit(date, 100d0, 0d0, 10d-6,                0d0, 0d0, 86400d0 / CalcPeriod(7378 / (1 - 1d-6)))
  orbit2 = GetOrbit(date, 100d0, 0d0, 1 - 7378d0 /  9379d0, 0d0, 0d0, 86400d0 / CalcPeriod(6378 + 3000d0))
  orbit3 = GetOrbit(date, 100d0, 0d0, 1 - 7378d0 / 11379d0, 0d0, 0d0, 86400d0 / CalcPeriod(6378 + 5000d0))
  orbit4 = GetOrbit(date, 100d0, 0d0, 1 - 7378d0 / 13379d0, 0d0, 0d0, 86400d0 / CalcPeriod(6378 + 7000d0))
  orbit5 = GetOrbit(date, 100d0, 0d0, 1 - 7378d0 / 15379d0, 0d0, 0d0, 86400d0 / CalcPeriod(6378 + 9000d0))
  ! 6-8
  orbit1 = GetOrbit(date, 100d0, 0d0, 1d-6, 0d0, 0d0, 86400d0 / CalcPeriod(6378 + 1000d0))
  orbit2 = GetOrbit(date, 45d0,  0d0, 1d-6, 0d0, 0d0, 86400d0 / CalcPeriod(6378 + 1000d0))
  orbit3 = GetOrbit(date, 90d0,  0d0, 1d-6, 0d0, 0d0, 86400d0 / CalcPeriod(6378 + 1000d0))
  orbit4 = GetOrbit(date, 135d0, 0d0, 1d-6, 0d0, 0d0, 86400d0 / CalcPeriod(6378 + 1000d0))
  orbit5 = GetOrbit(date, 175d0, 0d0, 1d-6, 0d0, 0d0, 86400d0 / CalcPeriod(6378 + 1000d0))

  ! call orbit1%Show
  ! call orbit2%Show

  ! do i = 0, 100
  !   call InitializeFile
  !   do j = 1, 100
  !     call DEBRIS(j)%orbit%Output(date + i, DEBRIS(j)%orbit%prd)
  !     call OutputPoint(DEBRIS(j)%orbit%axis_p * DEBRIS(j)%orbit%sma * (1 - DEBRIS(j)%orbit%ecc))
  !     call DEBRIS(j)%orbit%Move(1d0, 1)
  !   end do
  !   call unlink("plot\orbit" // Str(i) // ".dat")
  !   call rename(FILEPATH_ORBIT, "plot\orbit" // Str(i) // ".dat")
  ! end do

  do i = 0, 359
    call InitializeFile
    call orbit1%Output(date + i, orbit1%prd)
    call orbit2%Output(date + i, orbit2%prd)
    call orbit3%Output(date + i, orbit3%prd)
    call orbit4%Output(date + i, orbit4%prd)
    call orbit5%Output(date + i, orbit5%prd)
    call unlink("plot\orbit" // Str(i) // ".dat")
    call rename(FILEPATH_ORBIT, "plot\orbit" // Str(i) // ".dat")

    call OutputPoint(orbit1%axis_p * orbit1%sma * (1 - orbit1%ecc))
    call OutputPoint(orbit2%axis_p * orbit2%sma * (1 - orbit2%ecc))
    call OutputPoint(orbit3%axis_p * orbit3%sma * (1 - orbit3%ecc))
    call OutputPoint(orbit4%axis_p * orbit4%sma * (1 - orbit4%ecc))
    call OutputPoint(orbit5%axis_p * orbit5%sma * (1 - orbit5%ecc))
    call unlink("plot\point" // Str(i) // ".dat")
    call rename(FILEPATH_POINT, "plot\point" // Str(i) // ".dat")

    if (i >= 359) exit

    call orbit1%Move(1d0, 1)
    call orbit2%Move(1d0, 1)
    call orbit3%Move(1d0, 1)
    call orbit4%Move(1d0, 1)
    call orbit5%Move(1d0, 1)
  end do

  call Terminate

  contains
end program XTest
