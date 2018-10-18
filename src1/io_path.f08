module IO_Path
  implicit none

  ! character(50), parameter :: DIRPATH
  character(50), parameter :: FILEPATH_INPUT  = "input.txt"

  character(50), parameter :: FILEPATH_COST   = "ost.dat"
  character(50), parameter :: FILEPATH_oDEBRI = "result.csv"
  ! character(50), parameter :: FILEPATH_oDEBRI = "L:/result.csv"
  character(50), parameter :: FILEPATH_oHV    = "hypervolume.csv"
  character(50), parameter :: FILEPATH_oPARAM = "params.csv"

  character(50), parameter :: FILEPATH_DEBRI  = "dat/debri_elements.dat"
  character(50), parameter :: FILEPATH_RCS    = "dat/RCS_list.dat"

  character(50), parameter :: FILEPATH_ORBIT  = "orbit.txt"
end module IO_Path
