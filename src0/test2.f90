program cputime
  implicit none

  call RankPareto

  contains

  subroutine RankPareto
    logical :: valid(10)
    integer :: cost(10), value(10), rank(10), i

    cost  = (/1, 2, 2, 3, 2, 2, 4, 3, 4, 3/)
    value = (/1, 1, 2, 2, 3, 3, 2, 3, 3, 4/)

    valid = .true.

    where (valid)
      rank = (/(COUNT(valid .and. (cost(i) >= cost .and. value(i) >= value) .and. (cost(i) /= cost &
       .or. value(i) /= value)), i = 1, 10)/) + 1
    else where
      rank = -1
    end where

    print "(3i10)", (cost(i), value(i), rank(i), i = 1, 10)
    print *
  end subroutine RankPareto

  subroutine RankPareto_NS
    logical :: valid(10), non_dominated(10)
    integer :: cost(10), value(10), rank(10), num_dominate(10), i, j

    cost  = (/1, 2, 2, 3, 2, 2, 4, 3, 4, 3/)
    value = (/1, 1, 2, 2, 3, 3, 2, 3, 3, 4/)

    valid = .true.
    ! valid(1) = .false.

    where (valid)
      rank = 0
    else where
      rank = -1
    end where

    num_dominate = (/(COUNT(valid .and. (cost(i) >= cost .and. value(i) >= value) .and. (cost(i) /= cost &
     .or. value(i) /= value)), i = 1, 10)/)

    print "(3i10)", (cost(i), value(i), num_dominate(i), i = 1, 10)
    print *

    i = 1
    do while (ANY(num_dominate > 0))
      non_dominated = rank == 0 .and. num_dominate == 0

      where (non_dominated) rank = i

      do j = 1, 10
        if (non_dominated(j)) then
          where (       valid .and. ( &
            value(j) <= value .and.   &
            cost(j)  <= cost) .and. ( &
            value(j) /= value .or.    &
            cost(j)  /= cost))
            num_dominate = num_dominate - 1
          end where
        end if
      end do

      print "(4i10)", (cost(j), value(j), rank(j), num_dominate(j), j = 1, 10)
      print *

      if (i > 10) stop

      i = i + 1
    end do

    where (rank == 0 .and. num_dominate == 0) rank = i

    print "(3i3)", (cost(i), value(i), rank(i), i = 1, 10)

  end subroutine RankPareto_NS
end program cputime
