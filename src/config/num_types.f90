module num_types
  use, intrinsic :: iso_fortran_env
  use, intrinsic :: iso_c_binding
  private
  integer, public, parameter :: i2 = INT16
  integer, public, parameter :: i4 = INT32
  integer, public, parameter :: i8 = INT64
  integer, public, parameter :: sp = REAL32
  integer, public, parameter :: dp = REAL64
  integer, public, parameter :: qp = REAL128
  !> Global precision used in computations
  integer, public, parameter :: rp = dp
  integer, public, parameter :: c_rp = c_double
  integer, public, parameter :: xp = dp
  integer, public, parameter :: c_xp = c_double
end module num_types


