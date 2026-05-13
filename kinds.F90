!     
! File:   kinds.F90
! Author: zyy
!
!

module kinds

implicit none
public
integer, parameter :: i4b = selected_int_kind(9)
integer, parameter :: DP = selected_real_kind(14, 200)

end module kinds
