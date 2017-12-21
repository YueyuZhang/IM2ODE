!     
! File:   kinds.F90
! Author: zyy
!
! Created on 2013年4月19日, 上午8:58
!

module kinds

implicit none
public
integer, parameter :: i4b = selected_int_kind(9)
integer, parameter :: DP = selected_real_kind(14, 200)

end module kinds
