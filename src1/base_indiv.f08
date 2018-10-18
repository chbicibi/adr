module Base_Indiv
  use Base_Func
  use Solver_Orbit
  use Solver_Debri

  implicit none
  private
  public :: TIndiv
  public :: Counter_Eval
  ! public :: Dest, Dest_array

  type :: TIndiv
    integer              :: index       =  0
    integer              :: gener       =  0
    integer              :: num_chrom   =  0
    integer, allocatable :: order(:)
    real(8), allocatable :: var(:, :)
    real(8)              :: object(3)   =  0d0
    integer              :: rank_pareto =  0
    logical              :: valid       =  .false.
    logical              :: evald       =  .false.

    integer              :: c           =  0

    contains
    procedure            :: Init
    procedure            :: Eval
    procedure            :: CheckDup    => CheckDupIndiv
    procedure            :: Clone
    final                :: Dest
    final                :: Dest_array
  end type TIndiv

  integer                :: Counter_Eval =  0

  contains

  ! ----------------------------------------------------------------------------
  ! Constructor
  ! ----------------------------------------------------------------------------
  subroutine Init(self, gener, num_chrom)
    class(TIndiv), intent(inout) :: self
    integer,       intent(in)    :: gener, num_chrom
    integer                      :: i

    if (allocated(self%order)) deallocate(self%order)
    if (allocated(self%var))   deallocate(self%var)
    self%evald       = .false.

    self%gener       = gener
    self%num_chrom   = num_chrom
    self%order       = Shuffle(NUM_OBJECT_all, num_chrom)
    allocate(self%var(7, num_chrom))

    self%var(1, :)   = [(Random(10d0),               i = 1, num_chrom)]
    self%var(2, :)   = [(Random(0.3d0, 0d0) + 1.1d0, i = 1, num_chrom)]
    self%var(3, :)   = [(Random(0.3d0, 0d0) + 0.5d0, i = 1, num_chrom)]
    self%var(4:5, :) = reshape([(Random(1d0),        i = 1, num_chrom * 2)], [2, num_chrom])
    self%var(6:7, :) = reshape([(Random(50d0),       i = 1, num_chrom * 2)], [2, num_chrom])
    self%var(1, 1)   = self%var(1, 1) * 10d0
  end subroutine Init

  ! ----------------------------------------------------------------------------
  !
  ! ----------------------------------------------------------------------------
  subroutine Eval(self, index)
    class(TIndiv), intent(inout) :: self
    integer,       intent(in)    :: index
    real(8)                      :: date_s, del_v, del_t
    integer                      :: i

    if (self%evald) return

    self%valid = .true.
    self%evald = .true.

    self%object(2) =  sum(self%var(1, 2:))
    self%object(3) = -sum(DEBRIS(self%order)%rcs)

    do i = 1, self%num_chrom
      if (i == 1) then
        date_s = self%var(1, 1) + MJD(2015, 0d0)
        call Transfer([0, self%order(1)],  date_s, self%var(2:, i), del_v, del_t)
      else
        date_s = date_s + del_t + self%var(1, i)
        call Transfer(self%order(i - 1:i), date_s, self%var(2:, i), del_v, del_t)
      end if

      self%object(1) = self%object(1) + del_v
      self%object(2) = self%object(2) + del_t

      ! if (self%object(1) > 150d0 .or. self%object(2) > 365d0) then
      !   self%object = 0d0
      !   self%valid  = .false.
      !   exit
      ! end if
    end do

    self%index   = index
    ! self%index   = Counter_Eval
    ! Counter_Eval = Counter_Eval + 1
  end subroutine Eval

  logical function CheckDupIndiv(self, indivs, mode) result (res)
    class(TIndiv), intent(in) :: self
    type(TIndiv),  intent(in) :: indivs(:)
    integer,       intent(in) :: mode
    logical                   :: dup_order, dup_var
    integer                   :: i

    res = .false.

    do i = 1, size(indivs)
      if (btest(mode, 0)) then
        dup_order = all(self%order == indivs(i)%order)
      else
        dup_order = .true.
      end if

      if (btest(mode, 1)) then
        dup_var = all(abs(self%var - indivs(i)%var) < 1e-2)
      else
        dup_var = .true.
      end if

      res = dup_order .and. dup_var

      if (res) then
        return
      end if
    end do
  end function CheckDupIndiv

  function Clone(self) result (indiv)
    type(TIndiv),  pointer    :: indiv
    class(TIndiv), intent(in) :: self

    allocate(indiv)
    indiv = self
  end function Clone

  subroutine Dest(self)
    type(TIndiv) :: self
    ! print *, "Dest", self%index
    ! if (associated(self)) then
      if (allocated(self%order)) deallocate(self%order)
      if (allocated(self%var))   deallocate(self%var)
      ! deallocate(self)
    ! end if
    ! stop
  end subroutine Dest

  subroutine Dest_array(self)
    type(TIndiv) :: self(:)
    integer               :: i

    ! if (associated(self)) then
      do i = 1, size(self)
        ! print *, "Dest_a", self(i)%index
        if (allocated(self(i)%order)) deallocate(self(i)%order)
        if (allocated(self(i)%var))   deallocate(self(i)%var)
      end do
      ! deallocate(self)
    ! end if
  end subroutine Dest_array
end module Base_Indiv
