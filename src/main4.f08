module mod_main
  use util
  use orbit_func
  use orbit_base
  use orbit_debri
  use class_prob1
  use class_prob2
  use class_opt1
  use class_opt2
  implicit none

  private
  public :: main_test

  logical, allocatable :: LOC_ODD(:)

  contains

  subroutine main_test(opt)
    implicit none
    integer, intent(in) :: opt
    type(TOpt1) :: inner_optimizer
    type(TOpt2) :: outer_optimizer
    type(TProb2) :: outer_problem
    type(TIndiv2) :: outer_indiv
    integer, allocatable :: order(:)
    real(8) :: dv_best
    integer :: num_debris, nx, N, steps, start_time!, unit, i
    character(32) :: ts, tc, tm

    num_debris = 100
    ! **** ivar ***
    N = 1000
    steps = 100
    nx = 10
    ts = "T"
    tc = "ORDER"
    tm = "SWAP"

    call set_randomseed
    call set_odd(1000)
    call init_debri

    call select_debri(num_debris, order)
    call get_optimizer(inner_optimizer)

    print *, "*** initialize outer_problem ***"
    call outer_problem%initialize(num_debris, order, inner_optimizer, scaling)

    if (btest(opt, 0)) then
      call system_clock(start_time)
      call outer_problem%optimize(50)
      call outer_problem%save_table("table.csv")
      call elapsed_time(start_time)
    else
      call outer_problem%load_table("table.csv")
    end if

    print *, "*** initialize outer_optimizer ***"
    call outer_optimizer%initialize(nx=nx, m=2, N=N, selection=ts, crossover=tc, mutation=tm) ! NSGA2
    ! call outer_optimizer%initialize(nx=nx, m=2, N=N, T=5, g="CH", selection=ts, crossover=tc, mutation=tm) ! MOEAD
    call outer_optimizer%set_problem(outer_problem)
    call outer_optimizer%set_prototype(outer_indiv)
    ! outer_optimizer%elite_preservation = .true.

    print *, "*** optimize ***"
    call system_clock(start_time)
    call outer_optimizer%prepare_calculation
    call outer_optimizer%run(steps)
    call elapsed_time(start_time)

    print *, "*** output result ***"
    call outer_optimizer%best_obj(dv_best)
    print "(aes15.8)", "best=", dv_best

    ! *** NSGA2 ***
    call outer_optimizer%save_result("res3_nsga2_result.csv")
    call outer_optimizer%save_history("res3_nsga2_history.csv", elite="only")
    call outer_optimizer%save_history("res3_nsga2_history_all.csv", elite="all")

    ! *** MOEAD ***
    ! call outer_optimizer%save_result("res3_moead_result.csv")
    ! call outer_optimizer%save_history("res3_moead_history.csv")
  end subroutine main_test

  subroutine get_optimizer(inner_optimizer)
    implicit none
    type(TOpt1), intent(inout) :: inner_optimizer
    type(TIndiv1) :: inner_indiv
    integer :: nx, N
    character(32) :: ts, tc, tm

    N = 100
    ! **** dvar ***
    nx = 2
    ts = "T"
    tc = "BLX"
    tm = "PM"

    print *, "*** initialize inner_optimizer ***"
    call inner_optimizer%initialize(nx=nx, N=N, selection=ts, crossover=tc, mutation=tm)
    call inner_optimizer%set_prototype(inner_indiv)
    inner_optimizer%elite_preservation = .true.
    inner_optimizer%dup_rejection = .false.
    inner_optimizer%show_log = .false.
  end subroutine get_optimizer

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
    integer :: s

    s = size(variables)

    allocate(res(s))
    where (LOC_ODD(1:s)) ! Departure(day)
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

  integer :: argc, opt
  type(string), allocatable :: argv(:)
  integer, parameter :: INIT_TABLE = 0

  call get_argv(argc, argv)
  if (argc == 0) then
    print *, "no arguments"
    stop
  end if

  opt = INIT_TABLE

  select case (argv(1)%s)
  case ("1")
    call main_test(opt)
  case ("2")
    print *, btest(2, 0)
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
