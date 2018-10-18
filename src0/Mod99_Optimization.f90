module ModMultiStageGeneticAlgorithm
  use ModOperations
  use ModIndividual
  use ModOutput
  implicit none
  private
  public :: TypeParamsMSGA, MultiStageGeneticAlgorithm

  type :: TypeNode
    ! self data
    real(8)                 :: priority
    integer                 :: height, position, swap(2)
    ! other node
    type(TypeNode), pointer :: node_parent, node_next, node_child1, node_child2
  end type TypeNode

  type :: TypeChain
    integer                  :: num_node, num_node_set, num_node_select, index
    type(TypeChain), pointer :: chain_next
    type(TypeNode),  pointer :: node_front
  end type TypeChain

  type :: TypeParamsMSGA
    integer                      :: num_individual, num_city, limit_step, convergence
    real(8)                      :: rate_mutate

    integer                      :: generation
    type(TypeNode)               :: node
    type(TypeChain), allocatable :: chain(:)
    type(TypeGroup)              :: group
    type(TypeIndividual)         :: individual, individual_best
  end type TypeParamsMSGA

  contains
  subroutine MultiStageGeneticAlgorithm(paramsMSGA, individual_in)
    type(TypeParamsMSGA), intent(inout)        :: paramsMSGA
    type(TypeIndividual), intent(in), optional :: individual_in
    type(TypeChain),      allocatable          :: chain(:)
    type(TypeNode),       pointer              :: node, node_front
    type(TypeIndividual)                       :: individual_initial
    ! real(8) :: cool, temperature
    integer,              allocatable          :: route(:)
    integer                                    :: num_city, limit_step, convergence, counter
    integer                                    :: i, j

    num_city    = paramsMSGA%num_city
    ! cool = paramsSA%cool
    ! temperature = paramsSA%temperature
    limit_step  = paramsMSGA%limit_step
    convergence = paramsMSGA%convergence

    allocate(chain(10))
    do i = 1, 10
      chain(i)%index        = i
      chain(i)%num_node     = 0
      chain(i)%num_node_set = 1024
      nullify(chain(i)%node_front)
    end do

    if (PRESENT(individual_in)) then
      individual_initial = individual_in
    else
      individual_initial = CreateIndividual(num_city, 0, 2)
    end if

    paramsMSGA%individual      = individual_initial
    paramsMSGA%individual_best = individual_initial
    counter = 0
    ! call OutputRoute(individual_initial, 0)
    ! call OutputCost(individual_initial, 0)

    do i = 1, limit_step
      paramsMSGA%generation = i

      do j = 1, 10
        call AddNode(chain, j, num_city)
        call EvalChain(chain(j), paramsMSGA%individual)
        call SelectNode(chain(j))
        print *, chain(j)%num_node
      end do

      ! route = paramsMSGA%individual%route
      ! print '(10i4)', route

      stop

      node => chain(1)%node_front
      do while (associated(node))
        print *, node%priority
        node => node%node_next
      end do


      ! if (MOD(i, 10000) == 0) print '(i10, f10.5, i12, f10.5)',             &
      !   i, paramsMSGA%individual%distance, paramsMSGA%individual_best%generation, &
      !   paramsMSGA%individual_best%distance

      ! if (i == paramsMSGA%individual_best%generation) then
      !   if (MOD(counter, 100) == 0) call OutputRoute(paramsMSGA%individual_best, i)
      !   counter = counter + 1
      ! end if
      ! if (MOD(i, 10000) == 0) call OutputCost(paramsMSGA%individual_best, i)

      ! if (MOD(i, 10) == 0 .and. paramsMSGA%temperature < 1. .and. &
      !     i - paramsMSGA%individual_best%generation > convergence) then
      !   if (MOD(counter, 100) /= 1) call OutputRoute(paramsMSGA%individual_best, i)
      !   exit
      ! end if

    end do
    print *, counter
  end subroutine MultiStageGeneticAlgorithm

  ! private --------------------------------------------------------------------
  subroutine AddNode(chain, height, num_city)
    type(TypeChain), intent(inout) :: chain(:)
    integer,         intent(in)    :: height, num_city
    type(TypeNode),  pointer       :: node_front, node_new
    integer                        :: num_node, i

    node_front => chain(height)%node_front
    num_node   = chain(height)%num_node_set

    do i = 1, num_node
      node_new           => CreateNode(chain, height, num_city)
      node_new%node_next => node_front
      node_front         => node_new
    end do
    chain(height)%node_front => node_front
    chain(height)%num_node   = chain(height)%num_node + num_node
  end subroutine AddNode

  subroutine EvalChain(chain, individual)
    type(TypeChain),      intent(in) :: chain
    type(TypeIndividual), intent(in) :: individual
    type(TypeNode),       pointer    :: node

    node => chain%node_front
    do while (associated(node))
      call EvalNode(node, individual)
      node => node%node_next
    end do
  end subroutine EvalChain

  subroutine EvalNode(node, individual)
    type(TypeNode),       intent(inout), pointer :: node
    type(TypeIndividual), intent(in)             :: individual
    integer,              allocatable            :: route(:)
    real(8)                                      :: distance_0, distance_1, distance_2
    real(8)                                      :: delta_1, delta_2
    integer                                      :: height, i

    height     = node%height
    route      = individual%route
    distance_0 = CalcRoute(route)

    if (height == 1) then
      distance_1    = CalcRoute(AlterRoute(route, node, 0))

      delta_1       = ABS(distance_1 - distance_0)
      node%priority = 0.5 ** delta_1
    else
      if (associated(node%node_child1) .and. associated(node%node_child2)) then
        distance_1    = CalcRoute(AlterRoute(route, node, 0))

        distance_2    = CalcRoute(AlterRoute(route, node, 1))

        delta_1       = ABS(distance_1 - distance_0)
        delta_2       = ABS(distance_2 - distance_0)

        node%priority = (delta_2 - delta_1) ** 2
      else
        node%priority = 0.
      end if
    end if
  end subroutine EvalNode

  subroutine SelectNode(chain)
    type(TypeChain), intent(inout) :: chain
    type(TypeNode),  pointer       :: node, node_previous
    integer                        :: height, n, m
    real(8)                        :: sum, average

    height = chain%index

    n   = 0
    sum = 0

    node => chain%node_front
    do while (associated(node))
      if (ABS(node%priority) > 10e-5) then
        n    = n + 1
        sum  = sum + node%priority
      end if
      node => node%node_next
    end do

    average = sum / n

    node => chain%node_front
    nullify(node_previous)
    m = 0

    do while (associated(node))
      if (node%priority < average) then
        if (associated(node_previous)) then
          node_previous%node_next => node%node_next
        else
          chain%node_front => node%node_next
        end if
      else
        m = m + 1
      end if
      node_previous => node
      node          => node%node_next
    end do
    chain%num_node = m
  end subroutine SelectNode

  function CreateNode(chain, height, num_city)
    type(TypeNode),  pointer       :: CreateNode
    type(TypeChain), intent(inout) :: chain(:)
    integer,         intent(in)    :: height, num_city
    type(TypeNode),  pointer       :: node_child, node_child1, node_child2
    integer                        :: num_node, list_node(2), swap(2), i

    allocate(CreateNode)
    nullify(CreateNode%node_next)

    if (height == 1) then
      nullify(node_child1, node_child2)
      swap = Shuffle(num_city, 2, 1)
    else
      num_node   = chain(height - 1)%num_node
      list_node  = Shuffle(num_node, 2)
      node_child => chain(height - 1)%node_front

      do i = 1, num_node
        if (i == list_node(1)) node_child1 => node_child
        if (i == list_node(2)) node_child2 => node_child
        node_child => node_child%node_next
      end do
      swap = 0
    end if

    CreateNode%height      = height
    CreateNode%node_child1 => node_child1
    CreateNode%node_child2 => node_child2
    CreateNode%swap        = swap
  end function CreateNode

  recursive function AlterRoute(route_in, node, mode) result (route_out)
    integer,        allocatable :: route_out(:)
    integer,        intent(in)  :: route_in(:), mode
    type(TypeNode), intent(in)  :: node
    integer                     :: swap(2)
    ! mode -> 0: swap, 1: swap(reverse), 2: connect, 3: connect(reverse)

    route_out = route_in

    if (node%height == 1) then
      swap = node%swap
      where (route_out == swap(2)) route_out = 0
      where (route_out == swap(1)) route_out = swap(2)
      where (route_out == 0)       route_out = swap(1)
    else
      if (mode == 0) then
        if (associated(node%node_child1)) route_out = AlterRoute(route_out, node%node_child1, 0)
        if (associated(node%node_child2)) route_out = AlterRoute(route_out, node%node_child2, 0)
      else
        if (associated(node%node_child2)) route_out = AlterRoute(route_out, node%node_child2, 0)
        if (associated(node%node_child1)) route_out = AlterRoute(route_out, node%node_child1, 0)
      end if
    end if
  end function AlterRoute
end module ModMultiStageGeneticAlgorithm
