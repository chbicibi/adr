module IO_Input
  use IO_Path

  implicit none

  ! --------------------------------------------------------
  ! SEED_RAND: seed of random (-1 -> auto)
  ! --------------------------------------------------------
  integer :: SEED_RAND

  ! --------------------------------------------------------
  ! METHOD_TYPE:
  !   1 -> HC / SA
  !   2 -> GA
  !   3 -> GA (create initial individuals in HC)
  ! --------------------------------------------------------
  integer :: METHOD_TYPE

  ! --------------------------------------------------------
  ! object
  ! --------------------------------------------------------
  integer :: INPUT_NUM_OBJECT_ALL
  integer :: INPUT_NUM_OBJECT_SELECT

  ! --------------------------------------------------------
  ! coordinate plane
  ! --------------------------------------------------------
  real(8) :: INPUT_MAP_WIDTH      = 2d0
  real(8) :: INPUT_MAP_HEIGHT     = 2d0

  ! --------------------------------------------------------
  ! parameters of SA
  ! --------------------------------------------------------
  real(8) :: INPUT_SA_TEMPERATURE = 100d0
  real(8) :: INPUT_SA_COOL        = 0.9d0
  real(8) :: INPUT_SA_COST        = 1000d0
  integer :: INPUT_SA_CONVERGENCE = 1000
  integer :: INPUT_SA_LIMIT_STEP  = 100000

  ! --------------------------------------------------------
  ! parameters of GA
  ! --------------------------------------------------------
  integer :: INPUT_GA_NUM_SAMPLE
  integer :: INPUT_GA_NUM_SAMPLE_INNER
  integer :: INPUT_GA_NUM_ELITE
  integer :: INPUT_GA_CONVERGENCE
  integer :: INPUT_GA_LIMIT_STEP

  integer :: INPUT_GA_NUM_CROSS   ! number of change-cities in ordered crossover
  real(8) :: INPUT_GA_RATE_MUTATE

  contains

  subroutine ReadInputFile
    open(10, file = FILEPATH_INPUT)
      ! --------------------------------------------------------
      ! SEED_RAND
      ! --------------------------------------------------------
      read(10, *) SEED_RAND

      ! --------------------------------------------------------
      ! METHOD_TYPE:
      ! --------------------------------------------------------
      read(10, *) METHOD_TYPE

      ! --------------------------------------------------------
      ! object
      ! --------------------------------------------------------
      read(10, *) INPUT_NUM_OBJECT_ALL
      read(10, *) INPUT_NUM_OBJECT_SELECT

      ! --------------------------------------------------------
      ! parameters of GA
      ! --------------------------------------------------------
      read(10, *) INPUT_GA_NUM_SAMPLE
      read(10, *) INPUT_GA_NUM_SAMPLE_INNER
      read(10, *) INPUT_GA_NUM_ELITE
      read(10, *) INPUT_GA_CONVERGENCE
      read(10, *) INPUT_GA_LIMIT_STEP

      read(10, *) INPUT_GA_NUM_CROSS
      read(10, *) INPUT_GA_RATE_MUTATE
    close(10)
  end subroutine ReadInputFile

  subroutine CheckError
    if (INPUT_NUM_OBJECT_SELECT > INPUT_NUM_OBJECT_ALL) then
      print "(a)", "Error: Input 'INPUT_NUM_OBJECT_SELECT' &
        & must be smaller than or equal to 'INPUT_NUM_OBJECT_ALL'"
      stop
    end if

    if (INPUT_GA_NUM_CROSS >= INPUT_NUM_OBJECT_SELECT)then
      print "(a)", "Error: Input 'INPUT_GA_NUM_CROSS' &
        & must be smaller than 'INPUT_NUM_OBJECT_SELECT'"
        print "(a)", INPUT_GA_NUM_CROSS
      stop
    end if
    if (INPUT_GA_NUM_ELITE >= INPUT_GA_NUM_SAMPLE)then
      print "(a)", "Error: Input 'INPUT_GA_NUM_ELITE' &
        & must be smaller than 'INPUT_GA_NUM_SAMPLE'"
      stop
    end if
  end subroutine CheckError
end module IO_Input
