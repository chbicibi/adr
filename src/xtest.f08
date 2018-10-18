program XTest
  use Base
  use Base_Func
  use Solver_Orbit
  use Solver_Debri
  use IO_Output
  use IO_Path
  use XMod

  implicit none

  integer :: i, j, n, imax, jmax, num

  real(8) :: date, delt, st
  real(8) :: pos1(3), pos2(3), axis(3)
  real(8) :: vel1(3), vel2(3), vel1a(3), vel2a(3)
  type(TOrbit) :: orbit, orbit1, orbit2
  real(8) :: del_v, del_t, v, t, v2, t2, dt, dtx, dty, a_ave, a_sd
  real(8) :: out(100)
  type(Hash) :: ran(100)

  character(:), allocatable :: odir, ofile

  call Initialize

  call Transfer0([1, 2], [MJD(2015, 0d0), 1.2d0, (180d0 - 2d1) / 360d0, 0.5d0, 0.5d0, 10d0, 10d0], 0, del_v, del_t)
  print "(2f15.5)", del_v, del_t


  do i = 2, 10
    call Transfer0([i, i + 1], [MJD(2015, 0d0), 1.2d0, (180d0 - 2d1) / 360d0, 0.5d0, 0.5d0, 10d0, 10d0], 0, del_v, del_t)
    print "(2f15.5)", del_v, del_t
  end do
  stop

  ran = Sort([(Hash(i, DEBRIS(i)%orbit%ran), i = 1, 100)])

  ! a_ave = sum(DEBRIS%orbit%sma) / 100
  ! a_sd  = sqrt(sum((DEBRIS%orbit%sma - a_ave) ** 2) / 100)
  ! print *, a_ave, a_sd

  date = MJD(2015, 0d0)

  do i = 1, 1
    orbit1 = DEBRIS(1)%orbit
    orbit2 = DEBRIS(2)%orbit
    call orbit1%Move(date, 0)
    call orbit2%Move(date + 20d0 / 864d0, 0)
    pos1 = orbit1%GetPosition()
    pos2 = orbit2%GetPosition()
    vel1a = orbit1%GetVelocity()
    vel2a = orbit2%GetVelocity()
    call CalcLambert(date, 2000d0, pos1, pos2, orbit1%axis_w, vel1, vel2, orbit)
  end do
  print "(4f10.5)", vel1a - vel1, Norm(vel1a - vel1)
  print "(4f10.5)", vel2a - vel2, Norm(vel2a - vel2)
  print *, Norm(vel1a - vel1) + Norm(vel2a - vel2)

  call Terminate

  ! orbit1 = GetOrbit(date, 100d0, 0d0,  1d-6, 0d0, 0d0, 86400d0 / CalcPeriod(a_ave))
  ! orbit2 = GetOrbit(date, 100d0, 10d0, 1d-6, 0d0, 0d0, 86400d0 / CalcPeriod(a_ave + a_sd))
  orbit1 = DEBRIS(ran(1)%key)%orbit
  orbit2 = DEBRIS(ran(4)%key)%orbit
  ! orbit2%sma = orbit1%sma + 50
  ! call orbit2%Move(date, 0)

  date = date + 60
  axis = orbit1%axis_w + orbit2%axis_w
  ! dt   = orbit1%prd / 86400d0
  dtx  = 0.1d0
  dty  = 0.001d0

  print "(f7.1)", orbit1%sma
  print "(f7.1)", orbit2%sma
  print "(f8.3)", Deg(Angle(orbit1%axis_w, orbit2%axis_w))
  ! stop
  ! call orbit1%Show
  ! call orbit2%Show

  imax = 1000
  jmax = 200

  open(10, file = "dv.dat", status = "replace")
    do i = 0, imax
      st = date + dtx * dble(i)
      do j = 1, jmax
        delt = dty * dble(j)

        call orbit1%Move(st, 0)
        call orbit2%Move(st + delt, 0)
        pos1 = orbit1%GetPosition()
        pos2 = orbit2%GetPosition()

        call CalcLambert(st, delt * 86400d0, pos1, pos2, axis, vel1, vel2, orbit)
        del_v = Norm(vel1 - orbit1%GetVelocity()) + Norm(vel2 - orbit2%GetVelocity())

        ! write(10, "(3f15.5)") dble(i) / n, dble(j) / n, del_v
        write(10, "(3f15.5)") st - date, delt, del_v
      end do
      write(10, *)
    end do
  close(10)
  call Terminate

  num  = 4
  ! open(10, file = "dvp.dat", status = "replace")
    do n = 0, 99
      call InitializeFile
      ! i = int(n / 100d0 * 0.5d0  / dtx)
      ! j = int(n / 100d0 * 0.05d0 / dty) + 1
      ! i = n
      ! j = nint(i * 0.054 + 38) + 1
      ! j = nint(i * (65 - 55) / 100d0 + 55) + 1
      i = 0
      j = n * 2 + 1

      st   = date + dtx * dble(i)
      delt =        dty * dble(j)

      call orbit1%Move(st, 0)
      call orbit2%Move(st + delt, 0)
      pos1 = orbit1%GetPosition()
      pos2 = orbit2%GetPosition()

      call CalcLambert(st, delt * 86400d0, pos1, pos2, axis, vel1, vel2, orbit)
      del_v = Norm(vel1 - orbit1%GetVelocity()) + Norm(vel2 - orbit2%GetVelocity())

      print *, n, i, j, del_v

      call orbit1%Output(st, delt * 86400d0)
      call orbit2%Output(st, delt * 86400d0)
      call orbit%Output( st, delt * 86400)

      odir  = "plot\data" // Str(num)
      ofile = odir // "\orbit" // Str(n) // ".dat"
      call system("mkdir -p " // odir // " > NUL 2>&1")
      call unlink(ofile)
      call rename(FILEPATH_ORBIT, ofile)
    end do
  ! close(10)
  call Terminate

  do i = 1, -1
    date = date + dt
    call orbit1%Move(dt, 1)
    call orbit2%Move(dt, 1)
    pos1 = orbit1%GetPosition()
    pos2 = orbit2%GetPosition()

    call Transfer2(orbit1, orbit2, [date, 1d0, 0.5d0, 0d0, 0d0, 0d0, 0d0], 0, del_v, del_t)
    out(i) = del_v
  end do
  ! print "(6f10.5)", out

  call Terminate

  contains
end program XTest
