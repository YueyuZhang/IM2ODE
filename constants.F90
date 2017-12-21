!     
! File:   constants.F90
! Author: zyy
!
! Created on 2013年4月19日, 上午9:37
!

module constants
    use kinds
    implicit none
    public
    save
    
    real(dp), parameter :: pi = 3.14159265358979323846_DP
    
    integer(i4b), parameter :: max_struct = 250
    integer(i4b), parameter :: max_type = 15
    integer(i4b), parameter :: max_atom = 150
    
end module constants
   

