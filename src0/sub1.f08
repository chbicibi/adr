module Mod
  implicit none

  public

  type :: type
    integer  :: v = 0
    contains
    procedure :: new => create
  end type type

  interface type
    procedure :: const
  end interface type

  contains

  type(type) function const(num) result (this)
    integer :: num

    this%v = num
  end function const

  subroutine create(self)
    class(type) :: self

    self%v = 50
  end subroutine create
end module Mod
