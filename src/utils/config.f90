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
  ! number of entries in window_chi output file
  integer, parameter :: N_MEASUREMENT = 5
  integer, parameter :: NCHI = 3*(N_MEASUREMENT-1) + 8
  real(kind=dp), parameter :: PI = 3.141592653589793_dp
  integer, parameter :: IMEAS_RF = 11

end module constants

module config
  use constants

end module config