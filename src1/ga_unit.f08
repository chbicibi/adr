module GA_Unit
  use Base_Func
  use Base_Indiv
  use Base_Group
  use GA_Base
  use IO_Output

  implicit none

  private
  public :: TGAHandle
  public :: lim_all, step_output_window, step_output_file

  type :: TGAHandle
    ! ------------------------------------------------------
    ! parameters
    ! ------------------------------------------------------
    integer                  :: mode(4)            =  0

    integer                  :: num_indiv          =  0
    integer                  :: num_chrom          =  0
    integer                  :: lim_gener          =  0
    integer                  :: convergence        =  0
    integer                  :: num_target_all     =  0

    integer                  :: num_elite          =  0
    integer                  :: lim_create         =  5

    integer                  :: num_cross          =  0
    real(8)                  :: rate_mutate        =  0d0

    integer                  :: interval_migration =  0
    real(8)                  :: rate_migration     =  0d0

    ! ------------------------------------------------------
    ! status
    ! ------------------------------------------------------
    integer                  :: step               =  0
    integer                  :: gener              =  0
    integer                  :: status             =  0
    integer                  :: counter_update     =  0
    logical                  :: update             =  .false.

    ! ------------------------------------------------------
    ! variables
    ! ------------------------------------------------------
    type(TGroup),    pointer :: group              => null()
    type(TIndiv)             :: indiv_best

    ! ------------------------------------------------------
    ! controller
    ! ------------------------------------------------------
    integer                  :: serial              =  0
    type(TGAHandle), pointer :: subhandle           => null()
    type(TGAHandle), pointer :: subhandles(:)       => null()

    contains
    procedure                :: Init
    procedure                :: Operate
    procedure                :: Eval
    procedure                :: Select, Crossover, Mutate, Migration
  end type TGAHandle

  integer :: step_output_window = 1
  integer :: step_output_file   = 1

  integer :: lim_all = 100000

  contains

  subroutine Init(self, eval)
    class(TGAHandle), intent(inout)        :: self
    logical,          intent(in), optional :: eval
    type(TGroup),     pointer              :: group
    integer                                :: mode(4), num_indiv, num_chrom, i
    integer,          allocatable          :: list(:)

    mode      =  self%mode
    num_indiv =  self%num_indiv
    num_chrom =  self%num_chrom
    group     => self%group

    if (.not. associated(group)) then
      allocate(group)
      call group%Init(0, num_indiv, num_chrom)

      if ((.not. present(eval) .or. eval) .and. mode(2) < 3) then
        call group%indiv_ref%Init(0, num_chrom)
      end if
    end if

    self%status         =  1
    self%gener          =  0
    self%counter_update =  0
    self%update         =  .false.
    self%group          => group

    if (.not. present(eval) .or. eval) then
      call group%Eval(mode(3))
      group%indiv_best%gener = 0
      self%indiv_best        = group%indiv_best
    else
      return
    end if

    if (btest(mode(4), 1)) then
      select case (mode(1))
      case (0)
        list = Reverse(PACK(group%score%key, group%indivs(group%score%key)%rank_pareto == 1))
        do i = 1, group%num_elite
          call OutputDebri(group%indivs(list(i)), 0)
        end do
      case (1:)
        call OutputDebri(self%indiv_best, 0)
      end select
      call LineFeed
    end if
  end subroutine Init

  subroutine Operate(self)
    class(TGAHandle), intent(inout) :: self
    type(TIndiv)                    :: indiv_new
    type(TIndiv),     allocatable   :: indivs_new(:)
    type(TGroup),     pointer       :: group
    integer,          allocatable   :: pos(:)
    integer                         :: mode(4), gener, num_indiv, num_chrom
    integer                         :: num_elite, lim_create, parents(2), lim, i

    mode       =  self%mode
    gener      =  self%gener + 1
    num_chrom  =  self%num_chrom
    num_indiv  =  self%num_indiv
    num_elite  =  self%num_elite
    lim_create =  self%lim_create
    group      => self%group

    allocate(indivs_new(num_indiv))

    if (.not. group%evald) then
      call group%Eval(mode(3))
      self%indiv_best = group%indiv_best
    end if

    self%gener     = gener
    group%gener    = gener
    group%indivs%c = 0

    if (num_elite > 0) then
      select case (mode(1))
      case (0) ! pareto
        if (num_elite >= group%num_elite) then
          num_elite               = group%num_elite
          indivs_new(1:num_elite) = PACK(group%indivs, group%indivs%rank_pareto == 1)
        else
          ! num_elite                    = MAX(num_elite, 2)
          pos                     = Shuffle(group%num_elite - 1) + 1
          indivs_new(1)           = group%indivs(group%score(1)%key)
          indivs_new(2:num_elite) = PACK(group%indivs(group%score(pos)%key), &
                                         group%indivs(group%score(pos)%key)%rank_pareto == 1)
          ! indivs_new(2) = group%indivs(group%score(num_elite)%key)
        end if
      case (1:)
        if (num_elite >= num_indiv) then
          num_elite = 0
        else
          indivs_new(1:num_elite) = group%indivs(group%score(1:num_elite)%key)
        end if
      end select
    end if

    lim = lim_create
    i   = num_elite + 1

    do while (i <= num_indiv)
      call indiv_new%Init(gener, num_chrom)

      if (lim > 0) then
        call self%Select(parents)
        call self%Crossover(indiv_new, parents)
        call self%Mutate(indiv_new)
      ! else
      !   print "(a)", "ERROR (create)"
      !   indivs_new(i) = indiv_new
      !   lim              = lim_create
      !   i                  = i + 1
      end if
      if (i > 1 .and. indiv_new%CheckDup(indivs_new(1:i - 1), mode(2))) then
        lim = lim - 1
      else
        indivs_new(i) = indiv_new
        lim           = lim_create
        i             = i + 1
        ! group%indivs(parents)%c = group%indivs(parents)%c + 1
      end if
    end do

    call group%Set(indivs_new)
  end subroutine Operate

  subroutine Operate2(self)
    class(TGAHandle), intent(inout) :: self
    type(TIndiv)                    :: indiv_new
    type(TIndiv),     allocatable   :: indivs_new(:)
    type(TGroup),     pointer       :: group, subgroup
    integer,          allocatable   :: pos(:)
    integer                         :: mode(4), gener, num_indiv, num_chrom
    integer                         :: num_elite, lim_create, parents(2), lim, i

    mode       =  self%mode
    gener      =  self%gener + 1
    num_chrom  =  self%num_chrom
    num_indiv  =  self%num_indiv
    num_elite  =  self%num_elite
    lim_create =  self%lim_create
    group      => self%group

    self%gener     = gener
    group%gener    = gener

    allocate(subgroup)
    allocate(indivs_new(30))

    parents         = Shuffle(num_indiv)
    indivs_new(1:2) = group%indivs(parents)

    lim = lim_create
    i   = 1

    do while (i <= 30)
      call indiv_new%Init(gener, num_chrom)

      if (lim > 0) then
        call self%Crossover(indiv_new, parents)
        call self%Mutate(indiv_new)
      end if
      if (i > 1 .and. indiv_new%CheckDup(indivs_new(1:i - 1), mode(2))) then
        lim = lim - 1
      else
        indivs_new(i) = indiv_new
        lim           = lim_create
        i             = i + 1
      end if
    end do

    call subgroup%Set(indivs_new)
    call subgroup%Eval(mode(2))
  end subroutine Operate2

  subroutine Eval(self)
    class(TGAHandle), intent(inout) :: self
    type(TGroup),     pointer       :: group
    integer                         :: mode(4), gener
    integer,          allocatable   :: list(:)
    logical                         :: suspend
    integer                         :: i

    mode  =  self%mode
    gener =  self%gener
    group => self%group

    if (.not. group%evald) then
      call group%Eval(mode(3))
    end if

    if (self%counter_update == 0 .or. &
        group%indiv_best%object(mode(3)) < self%indiv_best%object(mode(3))) then
      self%indiv_best     = group%indiv_best
      self%update         = .true.
      self%counter_update = self%counter_update + 1
    else
      self%update         = .false.
    end if

    suspend = .false.

    if (gener >= self%lim_gener .or. &
        gener -  self%indiv_best%gener >= self%convergence) then
      suspend     = .true.
      self%status = 2
    end if

    if (btest(mode(4), 0)) then
      if (mod(gener, step_output_window) == 0) then
        print "(2(a, i6, f8.3), a, i3, a, i0)", "gene:", gener, group%indiv_best%object(1),         &
        "   best:", self%indiv_best%gener, self%indiv_best%object(1), "   elite:", group%num_elite, &
        "   total: ", Counter_Eval
      end if
    end if

    if (btest(mode(4), 1)) then
      if (mod(gener, step_output_file) == 0) then
        select case (mode(1))
        case (0)
          list = PACK(group%score%key, group%indivs(group%score%key)%rank_pareto == 1)
          do i = 1, group%num_elite
            call OutputDebri(group%indivs(list(i)), gener)
          end do
          call OutputHV(group%hv, gener)
        case (1:)
          ! call OutputDebri(indiv_best, gener)
          ! do i = 1, group%num_indiv ! output all
          do i = 1, group%num_elite
            call OutputDebri(group%indivs(group%score(i)%key), gener)
          end do
        end select
        call LineFeed
        call OutputHV(group%hv, gener)
      end if
    end if
  end subroutine Eval

  ! ----------------------------------------------------------------------------
  ! call operator selection
  ! ----------------------------------------------------------------------------
  subroutine Select(self, parents)
    class(TGAHandle), intent(in)  :: self
    integer,          intent(out) :: parents(:)
    type(THash),      allocatable :: score(:)
    type(TGroup),     pointer     :: group
    integer                       :: mode

    mode  =  self%mode(1)
    group => self%group
    score =  PACK(group%score, group%indivs(group%score%key)%c < 5)

    if (group%num_elite >= 2 .and. mode == 0) then
      parents = SelectPateto(group%num_indiv, group%indivs%rank_pareto)
    else
      select case (mode)
      case (0, 1)
        parents = SelectSingle_Roulette(score)
      case (2:)
        ! parents = SelectSingle_Roulette(group%score)
        parents = SelectSingle_Tournament(score)
      end select
    end if
  end subroutine Select

  ! ----------------------------------------------------------------------------
  ! call operator crossover
  ! ----------------------------------------------------------------------------
  subroutine Crossover(self, indiv, parents)
    class(TGAHandle), intent(in)    :: self
    type(TIndiv),     intent(inout) :: indiv
    integer,          intent(in)    :: parents(:)

    type(TGroup),     pointer       :: group
    type(TIAry),      allocatable   :: order_parent(:), order_child(:)
    type(TRAry2),     allocatable   :: var_parent(:), var_child(:)
    integer                         :: mode, num_parents, i

    mode        =  self%mode(2)
    num_parents =  size(parents)
    group       => self%group

    ! --------------------------------------------------------------------------
    ! integer variable 1 (order)
    ! --------------------------------------------------------------------------
    if (btest(mode, 0)) then
      allocate(order_parent(num_parents))
      do i = 1, num_parents
        order_parent(i) = TIAry(group%indivs(parents(i))%order)
      end do

      if (CheckDup(order_parent(1)%elm, order_parent(2)%elm)) then
        order_child = Crossover_Order(order_parent, self%num_cross)
      else
        order_child = Crossover_SinglePoint(order_parent)
      end if
      indiv%order = order_child(1)%elm
    else if (allocated(group%indiv_ref%order)) then
      indiv%order = group%indiv_ref%order
    end if

    ! --------------------------------------------------------------------------
    ! integer variable 2 (var)
    ! --------------------------------------------------------------------------
    if (btest(mode, 1)) then
      allocate(var_parent(num_parents))
      do i = 1, num_parents
        var_parent(i) = TRAry2(group%indivs(parents(i))%var)
      end do

      var_child = Crossover_BLXalpha(var_parent, 0.1d0)
      indiv%var = var_child(1)%elm
    else if (allocated(group%indiv_ref%var)) then
      indiv%var = group%indiv_ref%var
    end if
  end subroutine Crossover

  ! ----------------------------------------------------------------------------
  ! call operator mutation
  ! ----------------------------------------------------------------------------
  subroutine Mutate(self, indiv)
    class(TGAHandle), intent(in)    :: self
    type(TIndiv),     intent(inout) :: indiv
    integer                         :: mode

    mode = self%mode(2)

    ! --------------------------------------------------------------------------
    ! integer variable 1 (order)
    ! --------------------------------------------------------------------------
    if (btest(mode, 0)) then
      if (Random() < self%rate_mutate) then
        if (Random() < 0.5d0) then
          indiv%order = Mutate_Swap(indiv%order)
        else
          indiv%order = Mutate_Change(indiv%order, self%num_target_all)
        end if
      end if
    end if

    ! --------------------------------------------------------------------------
    ! integer variable 2 (var)
    ! --------------------------------------------------------------------------
    if (btest(mode, 1)) then
      if (Random() < self%rate_mutate) then
        indiv%var = Mutate_Neighbor(indiv%var)
      end if
    end if
  end subroutine Mutate

  ! ----------------------------------------------------------------------------
  ! migration
  ! ----------------------------------------------------------------------------
  subroutine Migration(self)
    class(TGAHandle), intent(inout) :: self
    type(TGAHandle),  pointer       :: subhandles(:)
    type(TIndiv)                    :: indiv
    integer,          allocatable   :: list(:)
    integer                         :: num_island, p1, p2, i

    num_island        =  self%num_indiv
    subhandles        => self%subhandles

    list  = Shuffle(num_island)
    p1    = Random(subhandles(list(1))%num_indiv)
    indiv = subhandles(list(1))%group%indivs(p1)

    do i = 1, num_island
      if (i < num_island) then
        p2 = Random(subhandles(list(i + 1))%num_indiv)
        subhandles(list(i))%group%indivs(p1) = subhandles(list(i + 1))%group%indivs(p2)
        p1 = p2
      else
        subhandles(list(i))%group%indivs(p1) = indiv
      end if
    end do

    do i = 1, num_island
      call subhandles(i)%group%Eval(self%mode(3))
    end do
  end subroutine Migration
end module GA_Unit
