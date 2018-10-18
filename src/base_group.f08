module Base_Group
  use Base_Func
  use Base_Indiv
  use Solver_Debri

  implicit none

  private
  public :: BaseGroup

  type :: BaseGroup
    logical                   :: evald           = .false.

    integer                   :: num_chrom       = 0
    integer                   :: num_indiv       = 0
    integer                   :: num_indiv_valid = 0
    integer                   :: num_elite       = 0
    type(Hash),   allocatable :: score(:)
    real(8)                   :: object_mean(3)  = 0d0
    real(8)                   :: hv              = 0d0

    type(BaseIndiv)              :: indiv_ref
    type(BaseIndiv)              :: indiv_best
    type(BaseIndiv), allocatable :: indivs(:)

    contains
    procedure                 :: Init
    procedure                 :: Set
    procedure                 :: Eval
    procedure                 :: CalcMean
    procedure                 :: RankPareto
    procedure                 :: RankPareto_NS
    procedure                 :: CalcHV
    final                :: Dest
    final                :: Dest_array
  end type BaseGroup

  contains

  ! ----------------------------------------------------------------------------
  ! Constructor
  ! ----------------------------------------------------------------------------
  subroutine Init(self, gener, num_indiv, num_chrom)
    class(BaseGroup), intent(inout) :: self
    integer,       intent(in)    :: gener, num_indiv, num_chrom
    integer                      :: i

    if (allocated(self%indivs)) deallocate(self%indivs)
    if (allocated(self%score))  deallocate(self%score)

    self%num_chrom = num_chrom
    self%num_indiv = num_indiv
    self%evald     = .false.

    if (num_indiv > 0) then
      allocate(self%indivs(num_indiv))
      do i = 1, num_indiv
        call self%indivs(i)%Init(gener, num_chrom)
      end do
    end if
  end subroutine Init

  subroutine Set(self, indivs)
    class(BaseGroup), intent(inout)              :: self
    type(BaseIndiv),  intent(inout), allocatable :: indivs(:)
    integer                                   :: i

    if (allocated(self%score)) deallocate(self%score)

    self%num_chrom = indivs(1)%num_chrom
    self%num_indiv = size(indivs)
    self%evald     = .false.

    call move_alloc(indivs, self%indivs)
  end subroutine Set

  ! ----------------------------------------------------------------------------
  !
  ! ----------------------------------------------------------------------------
  subroutine Eval(self, obj)
    class(BaseGroup), intent(inout) :: self
    integer,       intent(in)    :: obj
    integer,       allocatable   :: list(:)
    integer                      :: num_indiv, num_indiv_valid, i, n
    character(40)                :: num0, num1

    num_indiv      = size(self%indivs)
    self%num_indiv = num_indiv

    write(num0, "(i0)") self%num_chrom
    write(num1, "(i0)") self%num_chrom + 1

    ! --------------------------------------------------------------------------
    ! do i = 1, num_indiv
    !   if (self%indivs(i)%evald) cycle
    !   call self%indivs(i)%Eval
    ! end do

    n = count(.not. self%indivs%evald)
    if (n > 0) then
      list = pack([(i, i = 1, num_indiv)], .not. self%indivs%evald)
      !$omp parallel do num_threads(n)
      do i = 1, n
        call self%indivs(list(i))%Eval(Counter_Eval + i)
      end do
      !$omp end parallel do
      Counter_Eval = Counter_Eval + n
    end if
    ! --------------------------------------------------------------------------

    num_indiv_valid = count(self%indivs%valid)
    if (num_indiv_valid <= 1) then
      print *, "Error: all indivs are invalid"
      stop
    end if
    self%num_indiv_valid = num_indiv_valid

    if (allocated(self%score)) deallocate(self%score)
    self%score      = Sort(pack([(Hash(i, self%indivs(i)%object(obj)), i = 1, num_indiv)], self%indivs%valid))
    self%indiv_best = self%indivs(self%score(1)%key)

    call self%CalcMean
    call self%RankPareto_NS
    call self%CalcHV
    self%evald = .true.
  end subroutine Eval

  subroutine CalcMean(self)
    class(BaseGroup), intent(inout) :: self
    integer                      :: i

    self%object_mean =  [(sum(self%indivs%object(i), mask = self%indivs%valid) / self%num_indiv_valid, i = 1, 3)]
  end subroutine CalcMean

  subroutine RankPareto(self)
    class(BaseGroup), intent(inout) :: self
    type(BaseIndiv),  allocatable   :: indivs(:)
    integer                      :: num_indiv, i, j, k(2)

    call move_alloc(self%indivs, indivs)
    num_indiv =  self%num_indiv
    k         =  [1, 3]

    where (indivs%valid)
      indivs%rank = [(count([(indivs(j)%valid      .and. &
          all(indivs(i)%object(k) >= indivs(j)%object(k)) .and. &
          any(indivs(i)%object(k) /= indivs(j)%object(k)), j = 1, num_indiv)]), i = 1, num_indiv)] + 1
    else where
      indivs%rank = -1
    end where

    self%num_elite = count(indivs%rank == 1)
    call move_alloc(indivs, self%indivs)
  end subroutine RankPareto

  subroutine RankPareto_NS(self)
    class(BaseGroup), intent(inout) :: self
    type(BaseIndiv),  allocatable   :: indivs(:)
    integer,       allocatable   :: num_dominate(:)
    logical,       allocatable   :: non_dominated(:)
    integer,       allocatable   :: k(:)
    integer                      :: num_indiv, rank, i, j

    call move_alloc(self%indivs, indivs)
    num_indiv =  self%num_indiv
    k         =  [1, 3]

    num_dominate = [(count([(    indivs(j)%valid      .and. &
      all(indivs(i)%object(k) >= indivs(j)%object(k)) .and. &
      any(indivs(i)%object(k) /= indivs(j)%object(k)), j = 1, num_indiv)]), i = 1, num_indiv)]

    where (indivs%valid)
      indivs%rank = 0
    else where
      indivs%rank = -1
      num_dominate       = -1
    end where

    rank = 1
    do while (any(num_dominate > 0))
      non_dominated = indivs%rank == 0 .and. num_dominate == 0

      where (non_dominated) indivs%rank = rank

      do i = 1, num_indiv
        if (non_dominated(i)) then
          where ([(                    indivs(j)%valid      .and. &
            all(indivs(i)%object(k) <= indivs(j)%object(k)) .and. &
            any(indivs(i)%object(k) /= indivs(j)%object(k)), j = 1, num_indiv)])
            num_dominate = num_dominate - 1
          end where
        end if
      end do

      rank = rank + 1
    end do

    where (indivs%rank == 0 .and. num_dominate == 0) indivs%rank = rank

    self%num_elite = count(indivs%rank == 1)
    call move_alloc(indivs, self%indivs)
  end subroutine RankPareto_NS

  subroutine CalcHV(self)
    class(BaseGroup), intent(inout) :: self
    integer,       allocatable   :: key(:), list(:)
    integer                      :: num_indiv, num_elite, i

    num_indiv =  self%num_indiv
    num_elite =  self%num_elite
    key       =  self%score%key
    list      =  pack(key, self%indivs(key)%rank == 1)

    self%hv   = sum([(((self%indivs(list(i))%object(1) - self%indivs(list(i + 1))%object(1)) &
                      * self%indivs(list(i))%object(3)), i = 1, num_elite - 1)])        &
          - max(100d0 - self%indivs(list(num_elite))%object(1), 0d0) * self%indivs(list(num_elite))%object(3)
  end subroutine CalcHV

  subroutine Dest(self)
    type(BaseGroup) :: self
    ! print *, "Dest", self%index
    ! if (associated(self)) then
      if (allocated(self%indivs)) deallocate(self%indivs)
      if (allocated(self%score))  deallocate(self%score)
      ! deallocate(self)
    ! end if
    ! stop
  end subroutine Dest

  subroutine Dest_array(self)
    type(BaseGroup) :: self(:)
    integer         :: i

    ! if (associated(self)) then
      do i = 1, size(self)
        ! print *, "Dest_a", self(i)%index
        if (allocated(self(i)%indivs)) deallocate(self(i)%indivs)
        if (allocated(self(i)%score))  deallocate(self(i)%score)
      end do
      ! deallocate(self)
    ! end if
  end subroutine Dest_array
end module Base_Group
