module constants
  use, intrinsic :: iso_fortran_env
  implicit none

  integer, parameter :: cr = real32
  integer, parameter :: dp = real64
  integer, parameter :: MAX_STRING_LEN = 512
  integer, parameter :: IORD = 2
  real(kind=dp), parameter :: TRBDNDW = 0.3
  real(kind=dp), parameter :: APARM = 30.
  integer, parameter :: PASSES = 2

end module constants

module config
  use constants

end module config