module mod_main
  use util
  use orbit_func
  use orbit_base
  use orbit_debri
  use class_prob1
  use class_opt1
  use class_opt2
  use class_opt3
  implicit none

  private
  public :: main1, main2, main3

  logical, allocatable :: LOC_ODD(:)

  contains

  subroutine main1
    ! SOGA
    implicit none
    type(TProb1) :: problem
    type(TIndiv1) :: indiv
    type(TOpt1) :: optimizer
    integer :: nx, nx1, N, steps
    character(32) :: ts, tc, tm, tc1, tm1
    ! integer, allocatable :: best(:)
    real(8), allocatable :: dbest(:)
    integer, allocatable :: ibest(:), order(:)
    integer :: start_time, num_debris

    call set_randomseed
    call init_debri

    num_debris = 30

    N = 10000
    steps = 100
    ! **** ivar ***
    nx = 10
    ts = "T"
    tc = "ORDER"
    tm = "SWAP"
    ! **** dvar ***
    nx1 = nx * 2
    tc1 = "BLX"
    tm1 = "PM"

    print *, "*** initialize problem ***"
    call set_odd(nx1)
    call select_debri(num_debris, order)
    call problem%initialize(num_debris, order)

    print *, "*** initialize optimizer ***"
    call optimizer%initialize(nx=nx, N=N, selection=ts, crossover=tc, mutation=tm)
    call optimizer%initialize_ex(nx=nx1, crossover=tc1, mutation=tm1)
    call optimizer%set_problem(problem)
    call optimizer%set_prototype(indiv)
    call optimizer%set_scaling_function(scaling)
    optimizer%elite_preservation = .true.
    optimizer%dup_rejection = .false.

    ! call optimizer%best(best)
    ! call problem%print_route(best, "init.csv")

    print *, "*** optimize ***"
    call system_clock(start_time)
    call optimizer%prepare_calculation
    call optimizer%run(steps)
    call elapsed_time(start_time)

    print *, "*** output result ***"
    call optimizer%save_result("soga_result.csv")
    call optimizer%save_history("soga_history.csv", elite="only")
    ! call optimizer%save_history("soga_history_all.csv", elite="all")

    call optimizer%best(ibest)
    call optimizer%best(dbest)
    ! stop
    ! call flight(problem%debris, problem%start, ibest, dbest)
  end subroutine main1

  subroutine main2
    ! NSGA2
    implicit none
    type(TProb1) :: problem
    type(TIndiv1) :: indiv
    type(TOpt2) :: optimizer
    integer :: nx, nx1, N, steps
    character(32) :: ts, tc, tm, tc1, tm1
    ! integer, allocatable :: best(:)
    real(8), allocatable :: dbest(:)
    integer, allocatable :: ibest(:), order(:)
    integer :: start_time, num_debris

    call set_randomseed
    call init_debri

    num_debris = 100

    N = 100
    steps = 1000
    ! **** ivar ***
    nx = 10
    ts = "ROULETTE"
    tc = "ORDER"
    tm = "SWAP"
    ! **** dvar ***
    nx1 = nx * 2
    tc1 = "BLX"
    tm1 = "PM"

    print *, "*** initialize problem ***"
    call set_odd(nx1)
    call select_debri(num_debris, order)
    call problem%initialize(num_debris, order)

    print *, "*** initialize optimizer ***"
    call optimizer%initialize(nx=nx, N=N, m=2, selection=ts, crossover=tc, mutation=tm)
    call optimizer%initialize_ex(nx=nx1, crossover=tc1, mutation=tm1)
    call optimizer%set_problem(problem)
    call optimizer%set_prototype(indiv)
    call optimizer%set_scaling_function(scaling)
    call optimizer%prepare_calculation
    optimizer%elite_preservation = .true.
    optimizer%dup_rejection = .true.

    print *, "*** optimize ***"
    call system_clock(start_time)
    call optimizer%run(steps)
    call elapsed_time(start_time)

    print *, "*** output result ***"
    call optimizer%save_result("nsga2_result.csv")
    call optimizer%save_history("nsga2_history.csv", elite="only")
    call optimizer%save_history("nsga2_history_all.csv", elite="all")

    call optimizer%best(dbest)
    call optimizer%best(ibest)
    ! stop
    ! call problem%print_route(best, "best.csv")
    ! call flight(problem%dep, problem%arr, problem%start, dbest(1), dbest(2))
    ! call flight(DEBRIS(ibest(1))%orbit, DEBRIS(ibest(2))%orbit, problem%start, dbest)
  end subroutine main2

  subroutine main3
    ! MOEAD
    implicit none
    type(TProb1) :: problem
    type(TIndiv1) :: indiv
    type(TOpt3) :: optimizer
    integer :: nx, nx1, N, steps
    character(32) :: ts, tc, tm, tc1, tm1
    ! integer, allocatable :: best(:)
    real(8), allocatable :: dbest(:)
    integer, allocatable :: ibest(:), order(:)
    integer :: start_time, num_debris

    call set_randomseed
    call init_debri

    num_debris = 1000

    N = 100
    steps = 100
    ! **** ivar ***
    nx = 10
    ts = "ROULETTE"
    tc = "ORDER"
    tm = "SWAP"
    ! **** dvar ***
    nx1 = nx * 2
    tc1 = "BLX"
    tm1 = "PM"

    print *, "initialize problem"

    call select_debri(num_debris, order)
    call problem%initialize(num_debris, order)

    print *, "initialize optimizer"

    call optimizer%initialize(nx=nx, N=N, m=2, T=5, g="CH", crossover=tc, mutation=tm)
    call optimizer%initialize_ex(nx=nx1, crossover=tc1, mutation=tm1)
    call optimizer%set_problem(problem)
    call optimizer%set_prototype(indiv)
    call optimizer%set_scaling_function(scaling)
    call optimizer%prepare_calculation
    ! optimizer%elite_preservation = .true.
    ! optimizer%dup_rejection = .true.

    ! call optimizer%best(best)
    ! call problem%print_route(best, "init.csv")

    print *, "optimize"

    call system_clock(start_time)
    call optimizer%run(steps)
    call elapsed_time(start_time)

    call optimizer%save_result("moead_result.csv")
    call optimizer%save_history("moead_history.csv")

    call optimizer%best(dbest)
    call optimizer%best(ibest)
    ! stop
    ! call problem%print_route(best, "best.csv")
    ! call flight(problem%dep, problem%arr, problem%start, dbest(1), dbest(2))
    ! call flight(DEBRIS(ibest(1))%orbit, DEBRIS(ibest(2))%orbit, problem%start, dbest)
  end subroutine main3

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

  function scaling(variables) result(res)
    implicit none
    real(8), intent(in) :: variables(:)
    real(8), allocatable :: res(:)

    allocate(res(size(variables)))
    where (LOC_ODD)      ! Departure(day)
      res = variables * 0
    else where           ! Transfer(s)
      res = variables * 6000
    end where
  end function scaling

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

  subroutine set_odd(n)
    implicit none
    integer, intent(in) :: n
    integer :: i

    allocate(LOC_ODD(n), source=.false.)
    do i = 1, n, 2
      LOC_ODD(i) = .true.
    end do
  end subroutine set_odd
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
