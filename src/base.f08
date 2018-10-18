module Base
  use Solver_Debri
  use IO_Input
  use IO_Output

  implicit none

  private
  public :: Initialize, Terminate

  type :: Time
    integer :: begin_total, end_total, rate, max, diff
    real(8) :: begin_cpu,   end_cpu
  end type Time

  type(Time) :: progtime

  contains

  subroutine Initialize
    call system_clock(progtime%begin_total)
    call cpu_time(progtime%begin_cpu)
    call ReadInputFile
    call CheckError
    call SetRandomSeed(SEED_RAND)
    ! call InitializeFile

    if (METHOD_TYPE < 10) then
      call ReadDebri
    end if

    print "(a)", "[START]"
  end subroutine Initialize

  subroutine Terminate
    print *
    print "(a)", "[FINISH]"

    call system_clock(progtime%end_total, progtime%rate, progtime%max)
    call cpu_time(progtime%end_cpu)

    if (progtime%end_total < progtime%begin_total) then
      progtime%diff = (progtime%max - progtime%begin_total) + progtime%end_total + 1
    else
      progtime%diff = progtime%end_total - progtime%begin_total
    endif

    print "(a, f13.5)", "Calculate Time (TOTAL):", progtime%diff / DBLE(progtime%rate)
    print "(a, f13.5)", "Calculate Time (CPU)  :", progtime%end_cpu - progtime%begin_cpu
    stop
  end subroutine Terminate

  subroutine SetRandomSeed(fixed_num)
    integer,      intent(in), optional :: fixed_num
    integer,      allocatable          :: seed(:)
    character(:), allocatable          :: arg
    integer                            :: length, nrand, clock

    call random_seed(size = nrand)
    allocate(seed(nrand))

    if (present(fixed_num) .and. fixed_num > 0) then
      seed = fixed_num
    else if (command_argument_count() > 0) then
      call get_command_argument(1, length = length)
      allocate(character(length) :: arg)
      call get_command_argument(1, arg)
      seed = StoI(arg)
      print *, arg
    else
      call system_clock(count = clock)
      seed = clock
    end if

    call random_seed(put = seed)
  end subroutine SetRandomSeed
end module Base
