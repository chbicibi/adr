module class_main
  use interface
  use util
  use problem
  use individual
  use ga_unit
  use soga
  use nsga2
  use moead
  use orbit_func
  use orbit_base
  use orbit_debri
  implicit none

  private
  public :: TProb1, TIndiv1, TOpt1
  public :: DEBRI_ORDER, NUM_DEBRI

  type, extends(TProblem) :: TProb1
    type(TOrbit) :: dep, arr, trf
    real(8) :: start, theta, offset
    real(8), allocatable :: table(:, :)
    integer, allocatable :: order(:)
    integer :: num_debri

    contains

    generic :: initialize => initialize_prob2
    generic :: call => call_prob2

    procedure :: initialize_prob2
    procedure :: call_prob2
    ! procedure :: print_route
  end type TProb1


  type, extends(TIndiv) :: TIndiv1
    contains

    generic :: initialize => initialize_n2

    procedure :: initialize_n2 => initialize_indiv1
    procedure :: evaluate => eval_indiv1
    procedure :: is_same => is_same_indiv1
    ! procedure :: print_u => print_d
    procedure :: print_up => print_indiv1
    procedure :: print_header_u => print_header_indiv1
  end type TIndiv1


  ! type, extends(TSOGA) :: TOpt1
  !   integer :: num_var1
  !   class(TCrossover), allocatable :: crossover1
  !   class(TMutation), allocatable :: mutation1

  !   contains

  !   procedure :: initialize_ex => initialize_opt1
  !   procedure :: init_indiv => init_indiv_opt1
  !   procedure :: reproduce
  !   procedure :: reproduce_m => reproduce
  ! end type TOpt1


  type, extends(TNSGA2) :: TOpt1
    integer :: num_var1
    class(TCrossover), allocatable :: crossover1
    class(TMutation), allocatable :: mutation1

    contains

    procedure :: initialize_ex => initialize_opt1
    procedure :: init_indiv => init_indiv_opt1
    procedure :: reproduce
    procedure :: reproduce_m => reproduce
  end type TOpt1

  ! type, extends(TMOEAD) :: TOpt1
  !   integer :: num_var1
  !   class(TCrossover), allocatable :: crossover1
  !   class(TMutation), allocatable :: mutation1

  !   contains

  !   procedure :: initialize_ex => initialize_opt1
  !   procedure :: init_indiv => init_indiv_opt1
  !   procedure :: reproduce => reproduce_s
  !   procedure :: reproduce_m => reproduce
  ! end type TOpt1

  ! integer, parameter :: ROUND = 1
  integer, allocatable :: DEBRI_ORDER(:)
  integer :: NUM_DEBRI

  contains


  ! ============================================================================
  ! Problem Class
  ! ============================================================================

  subroutine initialize_prob1(this, n)
    implicit none
    class(TProb1), intent(inout) :: this
    integer, intent(in) :: n
    real(8) :: start

    start = mjd(2018, 1d0 * n)
    this%start = start
    call this%dep%initialize(epc=start, inc=99d0, ran=0d0, ecc=0d0, ap=0d0, ma=0d0, mm=86400d0/calc_period(7178d0))
    call this%arr%initialize(epc=start, inc=99d0, ran=30d0, ecc=0d0, ap=0d0, ma=30d0, mm=86400d0/calc_period(8078d0))
  end subroutine initialize_prob1

  subroutine initialize_prob2(this, n)
    implicit none
    class(TProb1), intent(inout) :: this
    integer, intent(in) :: n
    type(TOrbit) :: dep, arr
    real(8) :: date, dt, dv
    integer :: i, j

    this%num_debri = n
    this%start = mjd(2018, 0d0)

    allocate(this%table(0:n, 0:n))

    do j = 0, n
      do i = 0, n
        if (i == j) then
          this%table(i, j) = huge(0d0)
        else
          dep = DEBRIS(DEBRI_ORDER(i))%orbit
          arr = DEBRIS(DEBRI_ORDER(j))%orbit
          date = this%start
          dt = 3600d0

          call lambert(date, dt, dep, arr, dv)

          this%table(i, j) = dv
        end if
      end do
    end do
  end subroutine initialize_prob2

  subroutine call_prob1(this, indiv)
    implicit none
    class(TProb1), intent(inout) :: this
    class(TIndiv1), intent(inout) :: indiv
    type(TOrbit) :: dep, arr
    ! real(8) :: start, theta, offset
    integer :: idx(2)
    real(8) :: var(2)
    real(8) :: date, delta_t, delta_v

    idx = indiv%ivariables(1:2)
    dep = DEBRIS(idx(1))%orbit
    arr = DEBRIS(idx(2))%orbit

    var = this%scaling_func(indiv%dvariables)
    date = this%start + var(1) * SEC_DAY
    delta_t = var(2)

    call lambert(date, delta_t, dep, arr, delta_v)

    allocate(indiv%objectives(1), source=delta_v)
    indiv%feasible = .true.
  end subroutine call_prob1

  subroutine call_prob2(this, indiv)
    implicit none
    class(TProb1), intent(inout) :: this
    class(TIndiv1), intent(inout) :: indiv
    real(8) :: dist, rcs
    integer :: s, p, i

    dist = 0d0
    rcs = 0d0
    s = size(indiv%ivariables)

    do i = 1, s
      if (i == 1) then
        p = 0
      else
        p = indiv%ivariables(i-1)
      end if
      dist = dist + this%table(p, indiv%ivariables(i))
      rcs = rcs + DEBRIS(DEBRI_ORDER(indiv%ivariables(i)))%rcs
    end do

    allocate(indiv%objectives(2), source=[dist, 1d0/rcs])
    indiv%feasible = .true.
  end subroutine call_prob2


  ! ============================================================================
  ! Individual Class
  ! ============================================================================

  subroutine initialize_indiv1(this, nvar, nvar1)
    implicit none
    class(TIndiv1), intent(inout) :: this
    integer, intent(in) :: nvar, nvar1

    this%dvariables = random_array(nvar1)
    this%ivariables = shuffle(NUM_DEBRI, nvar)
    ! this%ivariables = [1995, 1681]
  end subroutine initialize_indiv1

  subroutine eval_indiv1(this, problem)
    implicit none
    class(TIndiv1), intent(inout) :: this
    class(TProblem), intent(inout) :: problem

    if (this%evaluated) return

    select type(problem)
    class is (TProb1)
      call problem%call(this)
    end select

    this%evaluated = .true.
  end subroutine eval_indiv1

  logical function is_same_indiv1(this, other) result(res)
    implicit none
    class(TIndiv1), intent(in) :: this
    class(TIndiv), intent(in) :: other

    ! res = all(this%dvariables == other%dvariables) .and. &
    !       all(this%ivariables == other%ivariables)
    ! res = all(this%ivariables == other%ivariables)
    res = all(this%dvariables == other%dvariables)
  end function is_same_indiv1

  subroutine print_header_indiv1(this, unit)
    implicit none
    class(TIndiv1), intent(in) :: this
    integer, intent(in) :: unit
    type(string), allocatable :: s(:)
    integer :: l(3)

    l(1) = size(this%objectives)
    l(2) = size(this%ivariables)
    l(3) = size(this%dvariables)
    allocate(s(sum(l)), source=[serstr("obj", integers(l(1))), &
                                serstr("ivar", integers(l(2))), &
                                serstr("dvar", integers(l(3)))])
    write(unit, "(a)", advance='no') join(s, ",")
  end subroutine print_header_indiv1

  subroutine print_indiv1(this, unit, proc)
    implicit none
    class(TIndiv1), intent(in) :: this
    integer, intent(in) :: unit
    procedure(func_1d_1d) :: proc
    character(:), allocatable :: format
    integer :: l

    l = size(this%objectives)
    if (l == 1) then
      format = "(es15.8," &
                   // str(size(this%ivariables)) // "(','i0)" &
                   // str(size(this%dvariables)) // "(','es15.8))"
    else
      format = "(es15.8," // str(l-1) // "(','es15.8)" &
                   // str(size(this%ivariables)) // "(','i0)" &
                   // str(size(this%dvariables)) // "(','es15.8))"
    end if
    write(unit, format, advance='no') this%objectives, this%ivariables, proc(this%dvariables)
  end subroutine print_indiv1


  ! ============================================================================
  ! Optimizer Class
  ! ============================================================================

  subroutine initialize_opt1(this, nx, crossover, mutation)
    implicit none
    class(TOpt1), intent(inout) :: this
    integer, intent(in) :: nx
    character(*), intent(in) :: crossover, mutation

    this%num_var1 = nx

    allocate(TCrossover::this%crossover1)
    allocate(TMutation::this%mutation1)

    call this%crossover1%initialize(crossover, 1d0, 0.5d0)
    call this%mutation1%initialize(mutation, 0.1d0, 20d0)
  end subroutine initialize_opt1

  subroutine init_indiv_opt1(this, indiv)
    implicit none
    class(TOpt1), intent(inout) :: this
    class(TIndiv), intent(inout) :: indiv

    select type(indiv)
    class is (TIndiv1)
      call indiv%initialize(this%num_var, this%num_var1)
    end select
  end subroutine init_indiv_opt1

  subroutine reproduce(this, index, children)
    implicit none
    class(TOpt1), intent(in) :: this
    integer, intent(in) :: index(:)
    class(TIndiv), intent(out), allocatable :: children(:)
    integer, allocatable :: parents_i(:, :), children_i(:, :)
    real(8), allocatable :: parents_d(:, :), children_d(:, :)
    integer :: p, i

    allocate(children(2), source=this%prototype)

    ! *** ivar ***
    parents_i = reshape([(this%population(index(i))%indiv%ivariables, i = 1, 2)], &
                        [this%num_var, 2])

    call this%crossover%call(parents_i, children_i)
    ! allocate(children(size(children_i, dim=2)), source=this%prototype)

    do i = 1, size(children)
      call children(i)%set_variables(children_i(:, i))
      call this%mutation%call(children(i)%ivariables)
      ! call children(i)%set_variables(this%population(index(i))%indiv%ivariables)
    end do

    ! *** dvar ***
    parents_d = reshape([(this%population(index(i))%indiv%dvariables, i = 1, 2)], &
                        [this%num_var1, 2])

    call this%crossover1%call(parents_d, children_d)
    ! allocate(children(size(children_d, dim=2)), source=this%prototype)

    do i = 1, size(children)
      p = min(i, size(children_d, dim=2))
      call children(i)%set_variables(children_d(:, p))
      call children(i)%clamp_variables(lower=0d0, upper=1d0)
      call this%mutation1%call(children(i)%dvariables)
      call children(i)%clamp_variables(lower=0d0, upper=1d0)
    end do
  end subroutine reproduce

  subroutine reproduce_s(this, index, child)
    implicit none
    class(TOpt1), intent(in) :: this
    integer, intent(in) :: index(:)
    class(TIndiv), intent(out), allocatable :: child
    class(TIndiv), allocatable :: children(:)

    call this%reproduce_m(index, children)
    child = children(1)
  end subroutine reproduce_s
end module class_main

module mod_main
  use util
  use class_main
  use orbit_func
  use orbit_base
  use orbit_debri
  implicit none

  private
  public :: main1, main2, main3

  contains

  ! subroutine main1
  !   ! SOGA
  !   implicit none
  !   type(TProb1) :: problem
  !   type(TIndiv1) :: indiv
  !   type(TOpt1) :: optimizer
  !   integer :: nx, nx1, N, steps, start_time
  !   character(32) :: ts, tc, tm, tc1, tm1
  !   ! integer, allocatable :: best(:)
  !   real(8), allocatable :: dbest(:)
  !   integer, allocatable :: ibest(:)

  !   call set_randomseed
  !   call init_debri

  !   NUM_DEBRI = 200

  !   N = 100
  !   steps = 100
  !   ! **** ivar ***
  !   nx = 2
  !   ts = "ROULETTE"
  !   tc = "ORDER"
  !   tm = "SWAP"
  !   ! **** dvar ***
  !   nx1 = 2
  !   tc1 = "BLX"
  !   tm1 = "PM"

  !   call select_debri(NUM_DEBRI, DEBRI_ORDER)
  !   call problem%initialize(NUM_DEBRI)

  !   call optimizer%initialize(nx=nx, N=N, selection=ts, crossover=tc, mutation=tm)
  !   call optimizer%initialize_ex(nx=nx1, crossover=tc1, mutation=tm1)
  !   call optimizer%set_problem(problem)
  !   call optimizer%set_prototype(indiv)
  !   call optimizer%set_scaling_function(scaling)
  !   call optimizer%prepare_calculation
  !   optimizer%elite_preservation = .true.
  !   optimizer%dup_rejection = .true.

  !   ! call optimizer%best(best)
  !   ! call problem%print_route(best, "init.csv")

  !   call system_clock(start_time)
  !   call optimizer%run(steps)
  !   call elapsed_time(start_time)

  !   call optimizer%save_result("soga_result.csv")
  !   call optimizer%save_history("soga_history.csv", elite="only")
  !   ! call optimizer%save_history("soga_history_all.csv", elite="all")

  !   call optimizer%best(dbest)
  !   call optimizer%best(ibest)
  !   ! stop
  !   ! call problem%print_route(best, "best.csv")
  !   ! call flight(problem%dep, problem%arr, problem%start, dbest(1), dbest(2))
  !   call flight(DEBRIS(ibest(1))%orbit, DEBRIS(ibest(2))%orbit, problem%start, dbest)
  ! end subroutine main1

  subroutine main1
    ! NSGA2
    implicit none
    type(TProb1) :: problem
    type(TIndiv1) :: indiv
    type(TOpt1) :: optimizer
    integer :: nx, nx1, N, steps, start_time
    character(32) :: ts, tc, tm, tc1, tm1
    ! integer, allocatable :: best(:)
    real(8), allocatable :: dbest(:)
    integer, allocatable :: ibest(:)

    call set_randomseed
    call init_debri

    NUM_DEBRI = 1000

    N = 100
    steps = 100
    ! **** ivar ***
    nx = 10
    ts = "ROULETTE"
    tc = "ORDER"
    tm = "SWAP"
    ! **** dvar ***
    nx1 = 2
    tc1 = "BLX"
    tm1 = "PM"

    call select_debri(NUM_DEBRI, DEBRI_ORDER)
    call problem%initialize(NUM_DEBRI)

    call optimizer%initialize(nx=nx, N=N, m=2, selection=ts, crossover=tc, mutation=tm)
    call optimizer%initialize_ex(nx=nx1, crossover=tc1, mutation=tm1)
    call optimizer%set_problem(problem)
    call optimizer%set_prototype(indiv)
    call optimizer%set_scaling_function(scaling)
    call optimizer%prepare_calculation
    optimizer%elite_preservation = .true.
    optimizer%dup_rejection = .true.

    ! call optimizer%best(best)
    ! call problem%print_route(best, "init.csv")

    call system_clock(start_time)
    call optimizer%run(steps)
    call elapsed_time(start_time)

    call optimizer%save_result("nsga2_result.csv")
    call optimizer%save_history("nsga2_history.csv", elite="only")
    call optimizer%save_history("nsga2_history_all.csv", elite="all")

    call optimizer%best(dbest)
    call optimizer%best(ibest)
    ! stop
    ! call problem%print_route(best, "best.csv")
    ! call flight(problem%dep, problem%arr, problem%start, dbest(1), dbest(2))
    ! call flight(DEBRIS(ibest(1))%orbit, DEBRIS(ibest(2))%orbit, problem%start, dbest)
  end subroutine main1

  subroutine main1
    ! MOEAD
    implicit none
    type(TProb1) :: problem
    type(TIndiv1) :: indiv
    type(TOpt1) :: optimizer
    integer :: nx, nx1, N, steps, start_time
    character(32) :: ts, tc, tm, tc1, tm1
    ! integer, allocatable :: best(:)
    real(8), allocatable :: dbest(:)
    integer, allocatable :: ibest(:)

    call set_randomseed
    call init_debri

    NUM_DEBRI = 1000

    N = 100
    steps = 100
    ! **** ivar ***
    nx = 10
    ts = "ROULETTE"
    tc = "ORDER"
    tm = "SWAP"
    ! **** dvar ***
    nx1 = 2
    tc1 = "BLX"
    tm1 = "PM"

    call select_debri(NUM_DEBRI, DEBRI_ORDER)
    call problem%initialize(NUM_DEBRI)

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
  end subroutine main1

  subroutine flight(dep, arr, start, dvar)
    implicit none
    type(TOrbit), intent(inout) :: dep, arr
    real(8), intent(in) :: start, dvar(:)
    type(TOrbit) :: trf
    real(8) :: date, dt, dv

    date = start + dvar(1) * SEC_DAY
    dt = dvar(2)

    call lambert(date, dt, dep, arr, dv, trf)

    print "(a,3es15.5)", "best dvar=", dvar
    print "(3(aes15.5))", "best date=", date, " dt= ", dt, " dv=", dv

    call dep%save(date, dt, "orbit1.csv", .true.)
    call arr%save(date, dt, "orbit2.csv", .true.)
    call trf%save(date, dt, "orbit3.csv", .true.)
  end subroutine flight

  function scaling(variables) result(res)
    implicit none
    real(8), intent(in) :: variables(:)
    real(8), allocatable :: res(:)
    real(8) :: low(2), upp(2)

    low = [0d0, 0d0]
    upp = [20 * 86400d0, 6000d0]
    res = low + (upp - low) * variables
  end function scaling

  subroutine init_debri
    implicit none
    character(:), allocatable :: file_tle, file_rcs
    real(8) :: base_time
    integer :: i

    file_tle = "../dat/debri_elements.txt"
    file_rcs = "../dat/RCS_list.txt"

    call load_debri(-1, file_tle, file_rcs)
    base_time = DEBRIS(1)%orbit%epc

    do i = 1, NUM_DEBRIS_ALL
      call DEBRIS(i)%orbit%move(base_time, 0)
    end do
  end subroutine init_debri

  subroutine select_debri(n, order)
    implicit none
    integer, intent(in) :: n
    integer, intent(out), allocatable :: order(:)
    real(8), allocatable :: data(:)
    integer, allocatable :: ord(:)

    allocate(data(NUM_DEBRIS_ALL), source=DEBRIS(1:)%orbit%ran)
    ord = sort(data)
    allocate(order(0:n), source=[0, ord(1:n)])
  end subroutine select_debri

  subroutine main2
    implicit none
    real(8), allocatable :: data(:)
    integer, allocatable :: order(:)
    integer :: unit, i, j, idx(2)

    call init_debri
    allocate(data(NUM_DEBRIS_ALL), source=DEBRIS(1:)%orbit%dran)
    order = sort(data)
    idx(1) = order(1)
    idx(2) = order(NUM_DEBRIS_ALL)

    do i = 0, 100
      call DEBRIS(idx(1))%orbit%move(mjd(2018, 1d0 * i), 0)
      call DEBRIS(idx(2))%orbit%move(mjd(2018, 1d0 * i), 0)

      print *, i, calc_angle(DEBRIS(idx(1))%orbit%axis_w, DEBRIS(idx(2))%orbit%axis_w) * DEGREES
    end do

    stop

    open(newunit=unit, file="debri_data.csv")
      write(unit, "(a)") "i,sma,inc,ran,ecc,ap,ma,mm,dran,dap"
      do i = 1, NUM_DEBRIS_ALL
        j = order(i)
        write(unit, "(i0,9(','es15.8))") j, DEBRIS(j)%orbit%sma, &
                                            DEBRIS(j)%orbit%inc * DEGREES, &
                                            DEBRIS(j)%orbit%ran * DEGREES, &
                                            DEBRIS(j)%orbit%ecc, &
                                            DEBRIS(j)%orbit%ap * DEGREES, &
                                            DEBRIS(j)%orbit%ma * DEGREES, &
                                            DEBRIS(j)%orbit%mm, &
                                            DEBRIS(j)%orbit%dran * DEGREES * 86400, &
                                            DEBRIS(j)%orbit%dap * DEGREES * 86400
      end do
    close(unit)
  end subroutine main2

  subroutine main3
    ! ΔVプロット
    implicit none
    type(TOrbit) :: dep, arr
    real(8), allocatable :: data(:)
    integer, allocatable :: order(:)
    real(8) :: date, dt, dv, dvs(10001, 21)
    integer :: unit, i, j, idx(2)

    call init_debri
    allocate(data(NUM_DEBRIS_ALL), source=DEBRIS(1:)%orbit%dran)
    order = sort(data)
    idx(1) = order(1)
    idx(2) = order(NUM_DEBRIS_ALL)

    ! call dep%initialize(epc=mjd(2018, 0d0), inc=99d0, ran=0d0, ecc=0d0, ap=0d0, ma=0d0, mm=86400d0/calc_period(7178d0))
    ! call arr%initialize(epc=mjd(2018, 0d0), inc=99d0, ran=30d0, ecc=0d0, ap=0d0, ma=30d0, mm=86400d0/calc_period(8078d0))

    dep = DEBRIS(idx(1))%orbit
    arr = DEBRIS(idx(2))%orbit

    do j = 0, 20
      if (mod(j, 10) == 0 .and. j > 0) print *, j
      do i = 0, 10000
        date = mjd(2018, 200d0 * i * SEC_DAY) + 10
        dt = 600d0 + 200d0 * j
        ! call dep%move(mjd(2018, 0.1d0 * i), 0)
        ! call arr%move(mjd(2018, 0.1d0 * i), 0)
        call lambert(date, dt, dep, arr, dv)
        dvs(i+1, j+1) = dv
      end do
    end do

    print *, minval(dvs), minloc(dvs)

    open(newunit=unit, file="dv_data1.csv")
      write(unit, "(a)") "date,dt,dv"
      ! do j = 0, 50
      !   do i = 0, 5000
      !     date = 100d0 * i * SEC_DAY
      !     dt = 1000d0 + 100d0 * j
      !     dv = dvs(i+1, j+1)
      !     write(unit, "(2(es10.3',')es10.3)") date, dt, dv
      !   end do
      ! end do
      do i = 0, 10000
        do j = 0, 20
          date = 200d0 * i * SEC_DAY
          dt = 600d0 + 200d0 * j
          dv = dvs(i+1, j+1)
          write(unit, "(2(es10.3' ')es10.3)") date, dt, dv
        end do
        write(unit, *)
      end do
    close(unit)
  end subroutine main3
end module mod_main

program main
  use mod_main
  implicit none

  call main1
end program main
