module IO_Output
  use Base_Func
  ! use Base_Indiv
  ! use Solver_Orbit
  ! use Solver_Debri
  use IO_Path

  implicit none

  type :: IO
    character, allocatable :: path

    contains
    procedure, nopass :: func
  end type IO

  interface
    function func(a)
      real(8) :: a
    end function func
  end interface

  contains

  subroutine InitializeFile
    ! open(10, file = FILEPATH_ROUTE, status = "replace")
    ! close(10)

    ! open(11, file = FILEPATH_COST, status = "replace")
    !   write(11, "(a, a10, a13, 2a15)") "#", "step", "step_create", "distance", "value"
    ! close(11)

    open(11, file = FILEPATH_ORBIT, status = "replace")
      write(11, "(a)") "#x,y,z"
      write(11, "(3i15)") 0, 0, 0
      write(11, *)
    close(11)

    open(12, file = FILEPATH_POINT, status = "replace")
      ! write(12, "(a)") "#x,y,z"
    close(12)

    open(13, file = FILEPATH_ARROW, status = "replace")
      ! write(13, "(a)") "#x,y,z,u,v,w"
    close(13)

    ! open(12, file = FILEPATH_oDEBRI, status = "replace")
    !   write(12, "(a)") "#step,delta_V,delta_T,totalRCS,rank,gener,serial,total,params"
    !   ! write(12, "(a)") "# [km/s],[m^2],[-]"
    ! close(12)

    ! open(13, file = FILEPATH_oHV, status = "replace")
    !   write(13, "(a)") "#step,eval,hv,hv(logscale)"
    !   ! write(13, "(a)") "# [km/s],[m^2],[-]"
    ! close(13)

    ! open(14, file = FILEPATH_oPARAM, status = "replace")
    ! close(14)
  end subroutine InitializeFile

  ! ----------------------------------------------------------------------------

  ! subroutine OutputOrbit(orbit_in, date, dt)
  !   type(TOrbit), intent(in) :: orbit_in
  !   real(8),      intent(in) :: date, dt
  !   type(TOrbit)             :: orbit
  !   real(8)                  :: ea1, ea2, dea
  !   integer                  :: i, n

  !   orbit = orbit_in%CalcOrbit(date)

  !   ea1 = MAtoEA(orbit%ma, orbit%ecc)
  !   ea2 = MAtoEA(orbit%ma + orbit%mm * dt, orbit%ecc)
  !   dea = ea2 - ea1

  !   n = int(100 * dt / orbit%prd) + 1

  !   open(10, file = FILEPATH_ORBIT, position = "append")
  !     do i = 0, n
  !       orbit%ea = ea1 + dea * i / dble(n)
  !       ! orbit = GetOrbit(orbit_in, date + dt * i / 30d0 / 86400d0)
  !       write(10, "(3f15.5)") orbit%GetPosition()
  !     end do
  !     write(10, *)
  !   close(10)
  ! end subroutine OutputOrbit

  subroutine OutputOrbit(data)
    real(8), intent(in) :: data(:, :)

    open(10, file = FILEPATH_ORBIT, position = "append")
      write(10, "(3f15.5)") data
      write(10, *)
    close(10)
  end subroutine OutputOrbit

  subroutine OutputPoint(data)
    real(8), intent(in) :: data(:)

    open(10, file = FILEPATH_POINT, position = "append")
      write(10, "(3f15.5)") data
    close(10)
  end subroutine OutputPoint

  subroutine OutputArrow(data)
    real(8), intent(in) :: data(:)

    open(10, file = FILEPATH_ARROW, position = "append")
      write(10, "(6f15.5)") data
    close(10)
  end subroutine OutputArrow

  ! ----------------------------------------------------------------------------

  ! subroutine OutputCost(indiv, gener)
  !   type(BaseIndiv), intent(in) :: indiv
  !   integer,              intent(in) :: gener
  !   integer                          :: i

  !   open(10, file = FILEPATH_COST, position = "append")
  !     write(10, "(i11, i13, 2f15.5)") gener, indiv%gener, &
  !       indiv%cost_total, indiv%value_total
  !   close(10)
  ! end subroutine OutputCost

  ! ----------------------------------------------------------------------------

  ! subroutine OutputDebri(indiv, step)
  !   type(BaseIndiv), intent(in)  :: indiv
  !   integer,         intent(in)  :: step
  !   character(:),    allocatable :: format
  !   character(8)                 :: step_str
  !   integer                      :: i
  !   ! return

  !   format   = "(a, 3(f11.5, ','), 2(i6, ','), 2(i10, ','), ' |,', " // Str(indiv%num_chrom) // "(i4, ','))"
  !   step_str = Str(step) // ","

  !   open(10, file = FILEPATH_oDEBRI, position = "append")
  !     write(10, format) step_str, abs(indiv%object), indiv%rank, &
  !                       indiv%gener, indiv%index, Counter_Eval, indiv%order
  !   close(10)

  !   open(11, file = FILEPATH_oPARAM, position = "append")
  !     write(11, "(a, 3f11.5, 5i4)") step_str, abs(indiv%object)
  !     write(11, "(5(i3, x, 7f11.5/))") (indiv%order(i), indiv%var(:, i), i = 1, indiv%num_chrom)
  !   close(11)
  ! end subroutine OutputDebri

  subroutine LineFeed
    open(10, file = FILEPATH_oDEBRI, position = "append")
      write(10, *)
    close(10)
  end subroutine LineFeed

  ! subroutine OutputHV(hv, step)
  !   real(8), intent(in) :: hv
  !   integer, intent(in) :: step

  !   open(10, file = FILEPATH_oHV, position = "append")
  !     write(10, "(2(i10, ','), 2(f12.3, ','))") step, Counter_Eval, hv, log10(hv)
  !   close(10)
  ! end subroutine OutputHV
end module IO_Output
