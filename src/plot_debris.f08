module mod_main
  use util
  use orbit_func
  use orbit_base
  use orbit_debri
  implicit none

  private
  public :: main1

  contains

  subroutine main1
    implicit none
    type(TOrbit) :: orbit
    integer, allocatable :: order(:)
    real(8) :: start, dt
    integer :: num_debris, i

    call set_randomseed
    call init_debri
    num_debris = 100
    call select_debri(num_debris, order)

    start = mjd(2018, 0d0)

    do i = 1, num_debris
      orbit = DEBRIS(order(i))%orbit
      dt = orbit%prd
      call orbit%save(start, dt, "debris.csv", .false.)
    end do
  end subroutine main1

  subroutine flight(debris, start, ivar, dvar)
    implicit none
    type(TDebri), intent(in) :: debris(0:)
    real(8), intent(in) :: start, dvar(:)
    integer, intent(in) :: ivar(:)
    type(TOrbit) :: dep, arr, trf
    real(8) :: date, dt, dv

    dep = debris(ivar(1))%orbit
    arr = debris(ivar(2))%orbit

    date = start + dvar(1)
    dt = dvar(2)

    call lambert(date, dt, dep, arr, dv, trf)

    print "(a,2i3,a,2es12.5)", "best ivar=", ivar, " dvar=", dvar
    print "(3(aes12.5))", "best date=", date, " dt=", dt, " dv=", dv

    print *, "*** output orbit ***"
    call dep%save(date, dt, "orbit1.csv", .true.)
    call arr%save(date, dt, "orbit2.csv", .true.)
    call trf%save(date, dt, "orbit3.csv", .true.)
  end subroutine flight

  subroutine init_debri
    implicit none
    character(:), allocatable :: file_tle, file_rcs
    real(8) :: base_time
    integer :: i

    file_tle = "../dat/debri_elements.txt"
    file_rcs = "../dat/RCS_list.txt"

    call load_debri(-1, file_tle, file_rcs)
    ! base_time = DEBRIS(1)%orbit%epc
    base_time = mjd(2018, 0d0)

    do i = 0, NUM_DEBRIS_ALL
      call DEBRIS(i)%orbit%move(base_time, 0)
    end do
  end subroutine init_debri

  subroutine select_debri(n, order)
    implicit none
    integer, intent(in) :: n
    integer, intent(out), allocatable :: order(:)
    real(8), allocatable :: data(:)

    allocate(data(NUM_DEBRIS_ALL), source=DEBRIS(1:)%orbit%ran)
    allocate(order(0:n), source=[0, sort(data, n)])
  end subroutine select_debri
end module mod_main

program main
  use util
  use mod_main
  implicit none

  integer :: argc
  type(string), allocatable :: argv(:)

  call get_argv(argc, argv)
  if (argc == 0) then
    print *, "no arguments"
    stop
  end if

  select case (argv(1)%s)
  case ("1")
    call main1
  case default
    print *, "unknown argument: ", argv(1)%s
  end select

  contains

  subroutine get_argv(argc, argv)
    implicit none
    integer, intent(out) :: argc
    type(string), intent(out), allocatable :: argv(:)
    character(:), allocatable :: arg
    integer :: i, length, status

    argc = command_argument_count()
    allocate(argv(0:argc))

    do i = 0, argc
      call get_command_argument(i, length=length, status=status)
      if (status == 0) then
        allocate (character(length) :: arg)
        call get_command_argument(i, arg, status=status)
        if (status == 0) argv(i) = string(arg)
        deallocate (arg)
      end if
      if (status /= 0)  print *, 'Error', status, 'on argument', i
    end do
  end subroutine get_argv
end program main
