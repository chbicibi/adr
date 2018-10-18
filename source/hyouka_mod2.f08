include "lib/base_func.f08"
include "lib/solver_orbit.f08"
include "lib/solver_debri.f08"

module hyouka_mod
  use globals
  use Base_Func
  use Solver_Orbit
  use Solver_Debri

  implicit none

  integer :: NUM_OBJECT_SELECT
  logical :: is_show = .false.

contains

  subroutine read_debridata
    FILEPATH_DEBRI = "debri_elements.txt"
    FILEPATH_RCS   = "RCS_list.txt"
    NUM_OBJECT_ALL = db
    NUM_OBJECT_SELECT = tar

    call ReadDebri

    ! 出発軌道
    DEBRIS(0)%orbit = GetOrbit(MJD(2017, 0d0), 99d0, 0d0, 1d-2, 0d0, 0d0, 86400d0 / CalcPeriod(6578d0))
  end subroutine read_debridata

  subroutine hyouka(order_list, fit, total_RCS, total_flight_time, k)
    integer, intent(in) :: order_list(n, tar)
    real(8), intent(out) :: fit(n)
    real(8), intent(out) :: total_RCS(n)
    real(8), intent(out) :: total_flight_time(n)
    integer :: i, k

    do i = 1, n
      call hyouka_lambert(order_list(i, :), fit(i), total_RCS(i), total_flight_time(i))

      !-------------------------------不良個体処理(仮)
      if (fit(i) < 1.0d-5) then
        print *, "error: ", order_list(i, :)
        fit(i) = 99.0
      end if

      if (is_show) then
        print "(2(i5, a), 3(f10.5, a))", k, ": ", i, " total delV=", fit(i), " [km/s] RCS=" &
                                       , total_RCS(i), " [m2] delT=", total_flight_time(i), " day"
      end if
    end do
  end subroutine hyouka

  subroutine hyouka_lambert(order, delta_V, RCS, flight_time)
    integer, intent(in) :: order(:)
    real(8), intent(out) :: delta_V
    real(8), intent(out) :: RCS
    real(8), intent(out) :: flight_time
    real(8), allocatable :: periods(:)
    type(TOrbit), allocatable :: orbits(:)
    type(TOrbit) :: orbit, orbit1, orbit2
    real(8) :: date_start, period_trans
    real(8) :: pos1(3), pos2(3), axis(3)
    real(8) :: vel1(3), vel2(3), vel1_res(3), vel2_res(3), dv
    integer :: i

    orbits = DEBRIS([0, order])%orbit
    periods = [(3600d0, i = 1, tar)]

    RCS = sum(DEBRIS(order)%rcs)
    delta_V = 0d0
    flight_time = 0d0

    date_start = orbits(1)%epoch

    do i = 1, tar
      orbit1 = orbits(i)
      orbit2 = orbits(i + 1)

      period_trans = periods(i)

      call orbit1%Move(date_start, 0)
      call orbit2%Move(date_start + period_trans / 86400d0, 0)

      pos1 = orbit1%GetPosition()
      pos2 = orbit2%GetPosition()
      vel1 = orbit1%GetVelocity()
      vel2 = orbit2%GetVelocity()

      call CalcLambert(date_start, period_trans, pos1, pos2, orbit1%axis_w, vel1_res, vel2_res, orbit)

      dv = Norm(vel1_res - vel1) + Norm(vel2_res - vel2)

      delta_V = delta_V + dv
      flight_time = flight_time + period_trans / 86400d0

      date_start = date_start + period_trans / 86400d0
    end do
  end subroutine hyouka_lambert
end module hyouka_mod
