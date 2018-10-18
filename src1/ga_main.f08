module GA_Main
  use GA_Base
  use GA_Unit
  use Base_Indiv
  use Base_Group

  implicit none

  ! private
  ! public :: GAMain1, GAMain2

  contains

  subroutine GAMain0(handle)
    type(TGAHandle), intent(inout), pointer :: handle

    if (handle%status == 0) call handle%Init

    do
      if (ACCESS("stop", "") == 0) exit
      call handle%Operate
      call handle%Eval
      ! if (handle%status > 1 .or. ACCESS("stop", "") == 0) exit
      if (lim_all > 0 .and. Counter_Eval > lim_all) exit
    end do
  end subroutine GAMain0

  subroutine GAMain1(handle)
    type(TGAHandle), intent(inout), pointer :: handle
    type(TGAHandle),                pointer :: subhandles(:)
    integer                                 :: num_indiv, interval_migration, i

    do
      if (ACCESS("stop", "") == 0) exit

      num_indiv  =  handle%num_indiv
      interval_migration = handle%interval_migration

      ! print *, handle%step

      select case (handle%step)
      case (0)
        call handle%Init(.false.)

        allocate(subhandles(num_indiv))

        do i = 1, num_indiv
          subhandles(i) = handle%subhandle
          call subhandles(i)%Init
        end do
        handle%subhandles => subhandles

        handle%step = 1
      case (1)
        handle%gener = handle%gener + 1
        do i = 1, num_indiv
          call handle%subhandles(i)%Operate
          call handle%subhandles(i)%Eval
        end do

        if (MOD(handle%gener, interval_migration) == 0) handle%step = 2
      case (2)
        handle%group%evald = .false.
        do i = 1, num_indiv
          handle%group%indivs(i) = subhandles(i)%indiv_best
        end do
        call handle%Eval
        call handle%Migration

        if (lim_all > 0 .and. Counter_Eval > lim_all) exit
        handle%step = 1
      end select
    end do
  end subroutine GAMain1

  subroutine GAMain2(handle)
    type(TGAHandle), intent(inout), pointer :: handle
    type(TGAHandle),                pointer :: subhandles(:)
    type(TIndiv),               allocatable :: indivs(:)
    type(TGroup),                   pointer :: group, group_sub
    real(8)                                 :: min, mean
    real(8),                    allocatable :: list(:)
    integer                                 :: gener, num_chrom, num_indiv, num_indiv_sub
    integer                                 :: i, j, c, p

    do
      if (ACCESS("stop", "") == 0) exit

      gener      =  handle%gener
      num_chrom  =  handle%num_chrom
      num_indiv  =  handle%num_indiv
      group      => handle%group
      subhandles => handle%subhandles

      select case (handle%step)
      case (0) ! init
        call handle%Init(.false.)

        allocate(subhandles(num_indiv))
        do i = 1, num_indiv
          subhandles(i) = handle%subhandle
        end do
        handle%subhandles => subhandles

        handle%step = 1
      case (1) ! allocate
        do i = 1, num_indiv
          if (associated(subhandles(i)%group)) deallocate(subhandles(i)%group)
          allocate(subhandles(i)%group)
          group_sub => subhandles(i)%group
          call group_sub%Init(0, 0, num_chrom)

          group_sub%indiv_ref       = group%indivs(i)
          group_sub%indiv_ref%gener = 0

          num_indiv_sub = subhandles(i)%num_indiv

          allocate(indivs(num_indiv_sub))
          do j = 1, num_indiv_sub
            if (j == 1 .and. group_sub%indiv_ref%rank_pareto == 1) then
              indivs(j) = group_sub%indiv_ref
            else
              call indivs(j)%Init(0, num_chrom)
            end if
          end do
          group_sub%num_indiv = num_indiv_sub
          call move_alloc(indivs, group_sub%indivs)
          call subhandles(i)%Init
        end do

        handle%step = 2
      case (2) ! RGA
        ! c = int(sqrt(dble(num_indiv * gener + 1)))
        ! c = int(exp(gener * 1d-1))
        ! c = int((1d0 + gener * 1d-1) ** 2) + 0
        c = gener * num_indiv / 10

        do i = 1, c
          if (all(subhandles%status > 1))               exit
          if (lim_all > 0 .and. Counter_Eval > lim_all) exit

          if (any(subhandles%gener == 0)) then
            p       = Pick(minloc(subhandles%gener, subhandles%status == 1), 1)
          else
            allocate(list(num_indiv))
            list(:) = [(subhandles(j)%group%score(1)%value, j = 1, num_indiv)]
            min     = minval(list)
            mean    = sum(list) / num_indiv
            ! p    = Pick(minloc(list / min - sqrt(log(dble(i)) / dble(2 * subhandles%gener)), &
            !   subhandles%status == 1), 1)
            ! p    = Pick(maxloc(2 / (1 + exp(list / min - 1)) + sqrt(log(dble(i)) / dble(2 * subhandles%gener)), &
            !   subhandles%status == 1), 1)
            p       = Pick(maxloc(-list / num_chrom + sqrt(log(dble(i)) / dble(2 * subhandles%gener)), &
              subhandles%status == 1), 1)
            ! p    = Pick(maxloc(-list + subhandles%gener * 1d-2, subhandles%status == 1), 1)
            deallocate(list)
          end if

          call subhandles(p)%Operate
          call subhandles(p)%Eval
        end do

        print "(i0)", c
        print "(25(20i4))", subhandles%gener

        handle%step = 3
      case (3) ! IGA
        group%evald = .false.

        allocate(indivs(num_indiv))
        do i = 1, num_indiv
          indivs(i) = subhandles(i)%indiv_best
        end do
        call move_alloc(indivs, group%indivs)

        call handle%Eval
        call handle%Operate

        if (lim_all < 1 .and. handle%status > 1) exit
        if (lim_all > 0 .and. Counter_Eval > lim_all) exit
        handle%step = 1
      end select
    end do
  end subroutine GAMain2

  recursive subroutine GAMain3(handle)
    type(TGAHandle), intent(inout), pointer :: handle
    ! type(TGAHandle),                pointer :: subhandle
    ! type(TGroup),                   pointer :: group

    ! subhandle => handle%subhandle

    select case (handle%step)
    case (0)
      call handle%Init
      handle%mode(2) = 1
      handle%step    = 1
      call GAMain3(handle)
    case (1)
      select case (handle%mode(2))
      case (1)
        handle%mode(2) = 2
      case (2)
        handle%mode(2) = 4
      case (4)
        handle%mode(2) = 1
      end select

      handle%step = 2
      call GAMain3(handle)
    case (2)
      call handle%Operate
      call handle%Eval

      if (handle%status == 1) then
        handle%step = 1
        call GAMain2(handle)
      end if
    end select
  end subroutine GAMain3
end module GA_Main
