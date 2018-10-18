
  subroutine CalcLambert(pos1, pos2, delta_t)
    real(8), intent(in) :: pos1(3), pos2(3), delta_t
    real(8)             :: r1, r2, delta_nu, cos_delta_nu, a, y, z, z_next
    real(8)             :: delta_t_, ddel_t_dz, s, c, dsdz, dcdz, chi
    real(8)             :: f, g, gd, sma, slr, ecc, v1(3), v2(3)
    integer             :: i, j

    r1           = Norm(pos1)
    r2           = Norm(pos2)
    cos_delta_nu = DotProduct(pos1, pos2) / (r1 * r2)
    delta_nu     = ACOS(cos_delta_nu)

    a            = -1 * SQRT(r1 * r2 * (1d0 + cos_delta_nu))
    z            = 0d0

    do j = 1, 10
      print *, j, z
      s            = SUM((/(       (-z) ** (i - 1) / Factorial(2 * i + 1), i = 1, 6)/))
      c            = SUM((/(       (-z) ** (i - 1) / Factorial(2 * i    ), i = 1, 6)/))
      dsdz         = SUM((/((-i) * (-z) ** (i - 1) / Factorial(2 * i + 3), i = 1, 5)/))
      dcdz         = SUM((/((-i) * (-z) ** (i - 1) / Factorial(2 * i + 2), i = 1, 5)/))
      y            = r1 + r2 - a * SGN(PI ** 2 - z ** 2) * SQRT(2d0 - z * c)
      if (y > 0) then
        chi          = SQRT(y / c)
        delta_t_     = (a * SQRT(y) + chi ** 3 * s) / rMU
        ddel_t_dz    = (0.125d0 * a * (a / chi + 3d0 * SQRT(y) * s / c) + chi ** 3 * (dsdz - 1.5d0 * s / c * dcdz)) / rMU
        z_next       = z + (delta_t - delta_t_) / ddel_t_dz

        print "(3f15.5)", z, s, c
        print "(10f8.3)", (/(       (-z) ** (i - 1), i = 1, 6)/)
        print "(10f8.3)", (/(       (-z) ** (i - 1) / Factorial(2 * i + 1), i = 1, 6)/)
        print "(10f8.3)", (/(       (-z) ** (i - 1) / Factorial(2 * i    ), i = 1, 6)/)
        print *

        if (ABS(delta_t - delta_t_) < 1d-6) exit
        z = z_next
      else
        ! temp
        print "(a)", "check"

        print "(f15.5)", z, delta_t_
        print *
        z = z - (delta_t - delta_t_) / (2d0 * ddel_t_dz)
      end if
    end do

    f  = 1d0 - y / r1
    g  = a * SQRT(y / MU)
    gd = 1d0 - y / r2
    v1 = (pos2 - f * pos1) / g
    v2 = (gd * pos2 - pos1) / g

    sma = chi ** 2 / z_next
    slr = r1 * r2 * (1d0 - cos_delta_nu) / y
    ecc = SQRT(1d0 - slr / sma)

    print "(3f15.5)", v1, v2
    print "(3f15.5)", sma, slr,ecc
  end subroutine CalcLambert