program Main
  use Base
  use Base_Indiv
  use GA_main
  use GA_base
  use IO_Input

  implicit none

  select case (METHOD_TYPE)
  case (0)
    step_output_window = 1000
    step_output_file   = 1000
    lim_all            = 10 ** 7
    call Main_Debri_0
  case (1)
    step_output_window = 100
    step_output_file   = 100
    lim_all            = 10 ** 7
    call Main_Debri_1
  case (2)
    step_output_window = 1
    step_output_file   = 1
    lim_all            = 10 ** 7
    call Main_Debri_2
  end select

  call Terminate

  contains

  subroutine Main_Debri_0
    type(TGAhandle), pointer :: handle
    real(8)                  :: result(3)

    allocate(handle)
    handle%mode           = [0, 3, 1, 3]
    handle%num_target_all = INPUT_NUM_OBJECT_ALL
    handle%num_chrom      = INPUT_NUM_OBJECT_SELECT
    handle%num_indiv      = INPUT_GA_NUM_SAMPLE
    handle%num_elite      = INPUT_GA_NUM_ELITE
    handle%num_cross      = INPUT_GA_NUM_CROSS
    handle%rate_mutate    = INPUT_GA_RATE_MUTATE
    handle%lim_gener      = INPUT_GA_LIMIT_STEP
    handle%convergence    = INPUT_GA_CONVERGENCE

    if (handle%mode(1) == 0) handle%num_elite = handle%num_indiv ! for pareto

    call GAMain0(handle)

    result = abs(handle%indiv_best%object)
    print "(a)", repeat("-", 80)
    print "(a, 3f15.5)", "Elite: ", result
    print "(a, i10)",    "Total: ", Counter_Eval
  end subroutine Main_Debri_0

  subroutine Main_Debri_1
    type(TGAhandle), pointer :: handle, subhandle
    real(8)                  :: result(3)

    allocate(handle, subhandle)

    handle%mode           = [0, 0, 1, 3]
    handle%num_target_all = INPUT_NUM_OBJECT_ALL
    handle%num_chrom      = INPUT_NUM_OBJECT_SELECT
    handle%num_indiv      = 30
    handle%num_elite      = INPUT_GA_NUM_ELITE
    handle%num_cross      = INPUT_GA_NUM_CROSS
    handle%rate_mutate    = INPUT_GA_RATE_MUTATE
    handle%lim_gener      = INPUT_GA_LIMIT_STEP
    handle%convergence    = INPUT_GA_CONVERGENCE
    handle%interval_migration =  100

    subhandle             = handle
    subhandle%mode        = [0, 3, 1, 0]
    subhandle%num_indiv   = INPUT_GA_NUM_SAMPLE

    if (subhandle%mode(1) == 0) subhandle%num_elite = subhandle%num_indiv ! for pareto
    handle%subhandle      => subhandle

    call GAMain1(handle)

    result = abs(handle%indiv_best%object)
    print "(a)", repeat("-", 80)
    print "(a, 3f15.5)", "Elite: ", result
    print "(a, i10)",    "Total: ", Counter_Eval
  end subroutine Main_Debri_1

  subroutine Main_Debri_2
    type(TGAhandle), pointer :: handle, subhandle
    real(8)                     :: result(3)

    allocate(handle, subhandle)

    handle%mode           = [0, 1, 1, 3]
    handle%num_target_all = INPUT_NUM_OBJECT_ALL
    handle%num_chrom      = INPUT_NUM_OBJECT_SELECT
    handle%num_indiv      = INPUT_GA_NUM_SAMPLE
    handle%num_elite      = INPUT_GA_NUM_ELITE
    handle%num_cross      = INPUT_GA_NUM_CROSS
    handle%rate_mutate    = INPUT_GA_RATE_MUTATE
    handle%lim_gener      = INPUT_GA_LIMIT_STEP
    handle%convergence    = INPUT_GA_CONVERGENCE

    subhandle             = handle
    subhandle%mode        = [2, 2, 1, 0]
    subhandle%num_indiv   = INPUT_GA_NUM_SAMPLE_INNER
    subhandle%num_elite   = 5
    subhandle%lim_gener   = 200
    subhandle%convergence = 30
    handle%subhandle      => subhandle

    if (handle%mode(1) == 0) handle%num_elite = handle%num_indiv ! for pareto

    call GAMain2(handle)

    result = abs(handle%indiv_best%object)
    print "(a)", repeat("-", 80)
    print "(a, 3f10.5)", "Elite: ", result
    print "(a, i10)",    "Total: ", Counter_Eval
  end subroutine Main_Debri_2
end program Main
