program main
  use orbit_func
  use orbit_base
  use orbit_debri
  use orbit_optimizer
  implicit none

  ! call test33
  call optimizer_main

  contains

  subroutine test1
    implicit none
    type(TOrbit) :: orbit1
    real(8) :: vel(3), pos(3)

    call orbit1%initialize(epc=mjd(2017, 0d0), inc=99d0, ran=0d0, ecc=0d0, ap=0d0, ma=0d0, mm=86400d0/calc_period(6578d0))

    vel = orbit1%get_velocity()
    pos = orbit1%get_position()

    print *, vel
    print *, pos
  end subroutine test1

  subroutine test2
    implicit none
    type(TOrbit) :: orbit1
    real(8) :: vel(3), pos(3)

    call load_debri(num_debris=-1, file_tle="../dat/debri_elements.txt", file_rcs="../dat/RCS_list.txt")
    orbit1 = DEBRIS(1)%orbit

    vel = orbit1%get_velocity()
    pos = orbit1%get_position()

    print *, orbit1%epc
    print *, vel
    print *, pos
  end subroutine test2

  subroutine test31
    implicit none
    type(TOrbit) :: dep, arr, trf
    ! real(8) :: mid(3), vel(3), pos(3)
    real(8) :: start, dt, dv

    print *, "test31"

    start = mjd(2018, 0d0)
    call dep%initialize(epc=start, inc=99d0, ran=0d0, ecc=0d0, ap=0d0, ma=0d0, mm=86400d0/calc_period(6578d0))
    call arr%initialize(epc=start, inc=99d0, ran=30d0, ecc=0d0, ap=0d0, ma=10d0, mm=86400d0/calc_period(6678d0))

    dt = dep%prd * 0.5d0

    ! mid = []

    call lambert(start, dt, dep, arr, dv, trf)

    print *, dv
    call dep%save(start, dt, "orbit1.csv", .true.)
    call arr%save(start, dt, "orbit2.csv", .true.)
    call trf%save(start, dt, "orbit3.csv", .true.)
  end subroutine test31

  subroutine test32
    implicit none
    type(TOrbit) :: dep, arr, trf
    ! real(8) :: mid(3), vel(3), pos(3)
    real(8) :: start, dd, dt, dv

    print *, "test32"

    start = mjd(2018, 0d0)
    call dep%initialize(epc=start, inc=99d0, ran=0d0, ecc=0d0, ap=0d0, ma=0d0, mm=86400d0/calc_period(6578d0))
    call arr%initialize(epc=start, inc=99d0, ran=30d0, ecc=0d0, ap=0d0, ma=0d0, mm=86400d0/calc_period(6678d0))

    call hohmann(start, dep, arr, trf, dd, dt, dv)

    print *, dv

    call dep%save(dd, dt, "orbit1.csv", .true.)
    call arr%save(dd, dt, "orbit2.csv", .true.)
    call trf%save(dd, dt, "orbit3.csv", .true.)
  end subroutine test32

  subroutine test33
    implicit none
    type(TOrbit) :: dep, arr, trf
    real(8) :: axis(3)!, vel(3), pos(3)
    real(8) :: start, dd, dt, dv, theta, offset

    print *, "test33"

    start = mjd(2018, 0d0)

    call dep%initialize(epc=start, inc=99d0, ran=0d0, ecc=0d0, ap=0d0, ma=0d0, mm=86400d0/calc_period(6578d0))
    call arr%initialize(epc=start, inc=99d0, ran=30d0, ecc=0d0, ap=0d0, ma=0d0, mm=86400d0/calc_period(6678d0))

    call intersection(dep, arr, axis, theta, offset)

    dt = dep%prd * 0.5d0
    dd = start + offset * SEC_DAY

    call dep%move(dd, 0)
    call lambert(dd, dt, dep, arr, dv, trf)

    print "('dv='f15.8)", dv

    call dep%save(dd, dt, "orbit1.csv", .true.)
    call arr%save(dd, dt, "orbit2.csv", .true.)
    call trf%save(dd, dt, "orbit3.csv", .true.)
  end subroutine test33

  subroutine test4
    implicit none
    type(TOrbit) :: dep, arr, trf
    real(8) :: axis(3)!, vel(3), pos(3)
    real(8) :: start, dd, dt, dv, theta, offset
    integer :: i

    print *, "test4"

    start = mjd(2018, 0d0)

    do i = 1, 100
      call dep%initialize(epc=start, inc=99d0, ran=0d0, ecc=0d0, ap=0d0, ma=0d0, mm=86400d0/calc_period(7178d0))
      call arr%initialize(epc=start, inc=99d0, ran=30d0, ecc=0d0, ap=0d0, ma=3d0 * i, mm=86400d0/calc_period(7378d0))

      call intersection(dep, arr, axis, theta, offset)

      dt = dep%prd * 0.5d0
      dd = start + offset * SEC_DAY

      call dep%move(dd, 0)
      call lambert(dd, dt, dep, arr, dv, trf)

      print "(i3f15.8)", i, dv
      ! if (dv > 20) exit
    end do

    call dep%save(dd, dt, "orbit1.csv", .true.)
    call arr%save(dd, dt, "orbit2.csv", .true.)
    call trf%save(dd, dt, "orbit3.csv", .true.)
  end subroutine test4

end program main
