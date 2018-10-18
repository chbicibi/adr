submodule (IO_Output) IO_Output2
contains
  module function point_dist(a, b) result(distance)
    type(point), intent(in) :: a, b
    real :: distance
    distance = sqrt((a%x - b%x)**2 + (a%y - b%y)**2)
  end function point_dist
end submodule IO_Output2
