module globals                                              ! �O���[�o���ϐ�mod
  integer, parameter :: n = 100                             ! n:�̐�(20) <- 40
  integer, parameter :: db = 2613                           ! �f�u����(100��)
  integer, parameter :: tar = 5                             ! ����ڕW��(5��)
  integer, parameter :: nselect = 80                        ! �I��(����)
  real(8), parameter :: crate = 0.9d0                       ! �����m�� <- 1.0
  real(8), parameter :: mttn_rate = 0.01d0                  ! �ψيm�� <- 0.0
  integer, parameter :: T = 500                             ! ���㐔

  !-------------------------------------�f�u���f�[�^
  integer gen(db)                                           ! ����(����, ���z���S�X�P�[��)
  real(8)  incl(db)                                         ! �X�Ίp[10]
  real(8)  omega(db)                                        ! ����_�o�x[10]
  real(8)  ecc(db)                                          ! ���S��[11]
  real(8)  argper(db)                                       ! �ߓ��_����[10]
  real(8)  ma(db)                                           ! ���ϋߓ_���p[12]
  real(8)  mm(db)                                           ! ���ω^��
  real(8)  smma(db)                                         ! �N�������a(���z���S�X�P�[��)
  integer  syukai_num(db)                                   ! ����
  real(8)  RCS(db)                                          ! ���[�_�[���˒f�ʐ�
end module globals

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

module random_mod                                           ! ��������郂�W���[��
  implicit none

contains
  subroutine set_seed                                       ! ����seed�ݒ�
    integer, allocatable :: seed(:)
    integer nrand, clock

    call random_seed(size = nrand)
    allocate(seed(nrand))
    if (.false.) then
      seed(1:nrand) = 123456789                             ! seed�Œ�(�e�X�g�p)
    else
      call system_clock(count = clock)
      seed(1:nrand) = clock
    end if

    call random_seed(put = seed)
    deallocate(seed)
  end subroutine set_seed

  subroutine random_single(r)
    real(8), intent(out) :: r

    call random_number(r)

  end subroutine random_single
end module random_mod

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

module debri_order_mod
  implicit none

contains
  subroutine debri_order_table(order_list)
    use globals
    use random_mod

    integer, intent(out) :: order_list(n, tar)
    integer i, j, k, num
    real(8) r

    do i = 1, n
      !---------------------------------����̐���
      do j = 1, tar
        call random_single(r)
        num = int(real(db) * r) + 1

        if (j > 1) then
          do
            if (any(order_list(i, 1:j-1) == num)) then
              call random_single(r)
              num = int(real(db) * r) + 1
            else
              exit
            end if
          end do
        end if

        order_list(i, j) = num
      end do
    end do

  end subroutine debri_order_table
end module debri_order_mod

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

! original version (using external solver)
! include "hyouka_mod1.f08"

! new version (using internal solver)
include "hyouka_mod2.f08"

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

! new version (only selecting elite individuals)
! include "pareto_hyouka_mod1.f08"

! new version (Non-dominated Sorting)
include "pareto_hyouka_mod2.f08"

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

module cumlation_mod                                        ! �m�����v�Z����mod
  use pareto_hyouka_mod

  implicit none

contains
  subroutine cumlation(probability, local_fit, c)
    use globals

    integer, intent(in) :: c
    real(8), intent(in) :: local_fit(:)
    real(8), intent(out), allocatable :: probability(:)
    integer i
    real(8) tot

    !---------------------------------------------------------------------------
    tot = sum(local_fit(1:n-c))                             ! �K���xtot���v�Z�������ł͋����̍��v���K���xtot�ɂȂ�

    !---------------------------------------------------------------------------
    allocate(probability(n - c))
    probability(1:n-c) = local_fit(1:n-c) / tot             ! �I���m�����v�Z

  end subroutine cumlation
end module cumlation_mod

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

! new version
! include "roulette_mod1.f08"

! new version (it require "pareto_hyouka_mod2.f08")
include "roulette_mod2.f08"

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

module crossover_mod                                        ! �������W���[��
  implicit none

contains
  subroutine nc_ket(nc, ns)                                 ! nc����T�u���[�`���i�z��`�󌈒�j
    use globals
    use random_mod

    integer, intent(out) :: nc
    integer, intent(out) :: ns(nselect)
    integer empty(nselect)                                  ! ns�ȊO
    integer i
    integer a
    real(8) r

    ns(1:nselect) = 0
    empty(1:nselect) = 0

    !-----------------------------------�e�̑I��
    do
      nc = 0                                                ! 0�N���A
      do i = 1, nselect
        call random_single(r)                               ! ��������
        if (r < crate) then                                 ! ��������������crate��菬�������`�F�b�N
          nc = nc + 1
          ns(nc) = i                                        ! ��������ΐ������ꂽ�����ɑΉ�������F�̔ԍ���ns()�ɒ~������
        else
          empty(i - nc) = i
        end if                                              ! ������nc�̓J�E���^�[
      end do
      if (nc == 0) then
        cycle
      end if
      exit
    end do

    !---------------------------------------------------------------------------
    if (mod(nc, 2) == 1) then                               ! �ϐ�nc������������
      call random_single(r)                                 ! ��Ȃ�Η�����V���ɍ쐬��������̐��F�̔ԍ�������
      a = int(real(nselect - nc) * r) + 1

      nc = nc + 1
      ns(nc) = empty(a)
    end if

  end subroutine nc_ket

  !-----------------------------------------------------------------------------
  subroutine crossover(select_order, select_fit, newgen, nc, ns)
    use globals
    use random_mod

    integer, intent(in) :: select_order(nselect, tar)
    real(8), intent(in) :: select_fit(nselect)
    integer, intent(in) :: nc
    integer, intent(in) :: ns(nselect)
    integer, intent(out), allocatable :: newgen(:, :)
    integer, parameter :: d = 3
    integer i, j, k
    integer c, lis
    real(8) r
    integer, allocatable :: cross_par(:, :)
    integer na(d)
    integer newch(nc, tar)
    integer nuki1(d), nuki2(d)
    integer aa, pattern
    integer cross_point

    allocate(newgen(nc, tar))

    !-----------------------------------cross_par�ɐe�����i�[
    allocate(cross_par(nc, tar))

    cross_par(1:nc, 1:tar) = select_order(ns(1:nc), 1:tar)

    !-----------------------------------��������
    aa = 1
    do k = 1, nc / 2
      !---------------------------------�����p�^�[������
      ! 2�̐e��1�ł������f�u���ԍ��������Ă����珇������(pattern = 1)
      ! ����ȊO��1�_����(pattern = 2)
      pattern = 2
      do i = 1, tar
        if (any(cross_par(aa, 1:tar) == cross_par(aa + 1, i))) then
          pattern = 1
          exit
        end if
      end do

      if (pattern == 1) then
        !-------------------------------�������� �p�^�[��1(��������)
        na(1:d) = 0                                         ! 0�N���A
        nuki1(1:d) = 0

        do i = 1, d                                         ! �����_����d�̐��������
          call random_single(r)
          lis = int(real(tar) * r) + 1
          if (i > 1) then
            do
              if (any(na(1:i-1) == lis)) then
                call random_single(r)
                lis = int(real(tar) * r) + 1
              else
                exit
              end if
            end do
          end if
          na(i) = lis
        end do

        nuki1(1:d) = cross_par(aa + 1, na(1:d))             ! �ep1����d�̐����ɑΉ�����ʒu�̐���nuki1()�Ɋi�[

        !-------------------------------�ep2����nuki1()�̒l��������(�����ł�0�ɂ���)��������`������
        do i = 1, tar
          newch(aa, i) = cross_par(aa, i)
          do j = 1, d
            if (cross_par(aa, i) == nuki1(j)) then
              newch(aa, i) = 0
              exit
            end if
          end do
        end do

        !-------------------------------nuki1�����Ԃ�0�̂Ƃ���ɓ���Ă���
        c = 1
        do i = 1, tar
          if (newch(aa, i) == 0) then
            newch(aa, i) = nuki1(c)
            c = c + 1
          end if
        end do

        !-----------------------------------------------------------------------
        na(1:d) = 0                                         ! 0�N���A
        nuki2(1:d) = 0

        do i = 1, d                                         ! �����_����d�̐��������
          call random_single(r)
          lis = int(real(tar) * r) + 1
          if (i > 1) then
            do
              if (any(na(1:i-1) == lis)) then
                call random_single(r)
                lis = int(real(tar) * r) + 1
              else
                exit
              end if
            end do
          end if
          na(i) = lis
        end do

        nuki2(1:d) = cross_par(aa, na(1:d))                 ! �ep1����5�̐����ɑΉ�����ʒu�̐���nuki1()�Ɋi�[

        !-------------------------------�ep2����nuki1()�̒l��������(�����ł�0�ɂ���)��������`������
        do i = 1, tar
          newch(aa + 1, i) = cross_par(aa + 1, i)
          do j = 1, d
            if (cross_par(aa + 1, i) == nuki2(j)) then
              newch(aa + 1, i) = 0
              exit
            end if
          end do
        end do

        c = 1
        do i = 1, tar
          if (newch(aa + 1, i) == 0) then
            newch(aa + 1, i) = nuki2(c)
            c = c + 1
          end if
        end do
        !-------------------------------�������� �p�^�[��1 �����܂�
      else
        !-------------------------------�������� �p�^�[��2(1�_����)
        cross_point = 0                                     ! 0�N���A

        call random_single(r)
        cross_point = int(real(tar - 1) * r) + 1

        ! �q1����
        newch(aa, 1:cross_point) = cross_par(aa, 1:cross_point)
        newch(aa, cross_point+1:tar) = cross_par(aa + 1, cross_point+1:tar)

        ! �q2����
        newch(aa + 1, 1:cross_point) = cross_par(aa + 1, 1:cross_point)
        newch(aa + 1, cross_point+1:tar) = cross_par(aa, cross_point+1:tar)
        !-------------------------------�������� �p�^�[��2 �����܂�
      end if

      ! �q��newgen�Ɋi�[
      newgen(aa:aa+1, 1:tar) = newch(aa:aa+1, 1:tar)

      aa = aa + 2
    end do

  end subroutine crossover
end module crossover_mod

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

module mttn_mod
  implicit none

contains
  subroutine mttn(nc, newgen, mttn_group)
    use globals
    use random_mod

    integer, intent(in) :: nc
    integer, intent(in) :: newgen(nc, tar)
    integer, intent(out), allocatable :: mttn_group(:, :)
    integer i, j
    integer a, b, c1, c2
    real(8) r
    integer as(nc, tar)

    allocate(mttn_group(nc, tar))

    as(1:nc, 1:tar) = newgen(1:nc, 1:tar)

    do i = 1, nc
      call random_single(r)
      if (r < mttn_rate) then
        call random_single(r)
        a = int(real(tar) * r) + 1
        call random_single(r)
        b = int(real(tar - 1) * r) + 1
        if (b == a) then
          b = b + 1
        end if

        c1 = as(i, a)
        c2 = as(i, b)

        as(i, a) = c2
        as(i, b) = c1
      end if
    end do

    mttn_group(1:nc, 1:tar) = as(1:nc, 1:tar)

  end subroutine mttn
end module mttn_mod

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

module add_mod                                              ! �W�c�ǉ�
  implicit none

contains
  subroutine addition(add_order)
    use globals
    use random_mod

    integer, intent(out) :: add_order(tar)
    integer i, j
    integer num
    real(8) r

    do i = 1, tar
      call random_single(r)
      num = int(real(db) * r) + 1
      if (i > 1) then
        do
          if (any(add_order(1:i-1) == num)) then
            call random_single(r)
            num = int(real(db) * r) + 1
          else
            exit
          end if
        end do
      end if
      add_order(i) = num
    end do

  end subroutine addition
end module add_mod

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

module output_mod
  implicit none

contains
  subroutine output_result(elite_fit, elite_RCS, c, k, elite_flight_time, elite_group)
    use globals

    real(8), intent(in) :: elite_fit(:)
    real(8), intent(in) :: elite_RCS(:)
    integer, intent(in) :: elite_group(:, :)
    integer, intent(in) :: c
    integer, intent(in) :: k
    real(8), intent(in) :: elite_flight_time(:)

    integer i

    if (k == 1) then
      open(50, file = '../result/delta_V_result.csv', status = 'old')
      read(50, '()')
      read(50, '()')
    else
      open(50, file = '../result/delta_V_result.csv', position = 'append')
    end if

    do i = 1, c
      write(50, '(f12.4, a, f6.4, a, f12.6, 6(a, i5))') &
     & elite_fit(i), ', ', elite_RCS(i), ', ', elite_flight_time(i), ', ', k, &
     & ', ', elite_group(i, 1), ', ', elite_group(i, 2), ', ', elite_group(i, 3), ', ', &
     & elite_group(i, 4), ', ', elite_group(i, 5)
    end do
    close(50)

  end subroutine output_result
end module output_mod

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

program main
  use globals
  use random_mod
  use debri_order_mod
  use hyouka_mod
  use pareto_hyouka_mod
  use cumlation_mod
  use roulette_mod
  use crossover_mod
  use mttn_mod
  use add_mod
  use output_mod

  implicit none

  integer i, j, k
  integer a
  integer c
  real(8) time_begin, time_end                             ! ���Ԍv���̂��߂̕ϐ�
  integer order_list(n, tar)                               ! �����W�c
  integer new_order_list(n, tar)                           ! ������W�c
  real(8) fit(n)                                           ! �K���x
  real(8) total_RCS(n)                                     ! ���vRCS
  real(8) total_flight_time(n)                             ! ���v��s����
  real(8), allocatable :: elite_flight_time(:)
  real(8), allocatable :: elite_RCS(:)                     ! ���5RCS
  real(8), allocatable :: local_RCS(:)
  integer, allocatable :: elite_group(:, :)                ! ���5��
  integer, allocatable :: local_group(:, :)                ! ���5�̂���菜��������
  real(8), allocatable :: elite_fit(:)                     ! �]���l�̏��5��
  real(8), allocatable :: local_fit(:)                     ! ���5�̈ȊO�̕]���l
  real(8), allocatable :: probability(:)                   ! �m��
  integer select_order(nselect, tar)                       ! �I��z��
  real(8) select_fit(nselect)
  integer nc
  integer ns(nselect)
  integer, allocatable :: newgen(:, :)
  integer, allocatable :: mttn_group(:, :)
  integer add_order(tar)
  integer min_root(tar)

  !-------------------------------------���Ԍv���i�J�n�j
  call cpu_time(time_begin)

  !-------------------------------------�����V�[�h�ݒ�
  call set_seed

  !-------------------------------------�f�u���f�[�^�ǂݍ���
  call read_debridata

  !-------------------------------------�����W�c����
  call debri_order_table(order_list)

  do k = 1, T
    !-------------------------------------�]��
    call hyouka(order_list, fit, total_RCS, total_flight_time, k)

    call pareto_hyouka(order_list, fit, elite_group, local_group, elite_fit, &
                       local_fit, total_RCS, elite_RCS, local_RCS, c,        &
                       elite_flight_time, total_flight_time)

    do i = 1, c
      write(*, '(2(f8.4, a), 5i5)') elite_fit(i), ' [km/s]', &
                                    elite_RCS(i), ' [m2]', elite_group(i, 1:5)
    end do

    !-------------------------------------�o��
    call output_result(elite_fit, elite_RCS, c, k, elite_flight_time, elite_group)

    !-------------------------------------�I������
    if (k == T) then
      write(*, '(a)') ''
      write(*, '(a)') '---------------------------------'
      write(*, '(a23)') 'Calculate End'
      write(*, '(a11, i6, a)') 'Total:' , k , ' Generation'
      write(*, '(a)') '---------------------------------'
      exit
    end if

    !-------------------------------------�I��
    call cumlation(probability, local_fit, c)
    call roulette_par(local_group, probability, select_order, local_fit, select_fit, c)

    !-------------------------------------����
    call nc_ket(nc, ns)

    if (c + nc > n) then ! 2017/11/4
      nc = n - c
      nc = nc - mod(nc, 2)
    end if

    print *, 'elite:', c, 'cross:', nc

    call crossover(select_order, select_fit, newgen, nc, ns)

    !-------------------------------------�ψ�
    call mttn(nc, newgen, mttn_group)

    !-------------------------------------������֎󂯌p��
    new_order_list(1:c, 1:tar) = elite_group(1:c, 1:tar)
    new_order_list(c+1:c+nc, 1:tar) = mttn_group(1:nc, 1:tar)

    !-------------------------------------�W�c�ǉ�
    if (c + nc < n) then
      do i = c + nc + 1, n
        call addition(add_order)
        new_order_list(i, 1:tar) = add_order(1:tar)
      end do
    end if

    order_list(1:n, 1:tar) = new_order_list(1:n, 1:tar)

    write(*, '(a)') ''
    write(*, '(a)') '---------------------------------'
    write(*, '(a11, i6, a)') 'Finish' , k , ' Generation'
    write(*, '(a)') '---------------------------------'
    deallocate(newgen)                                      !�z��̏�����
    deallocate(mttn_group)
    deallocate(probability)
  end do

  !-------------------------------------���Ԍv���i�I���j
  call cpu_time(time_end)

  write(*, '(a, f10.5, a)') 'Calculate Time:' , time_end - time_begin, ' s'

end program main
