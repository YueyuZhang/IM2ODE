!     
! File:   differencial_evolution.F90
! Author: zyy
!
! Created on 2013年4月23日, 下午9:47
!

module differencial_evolution
    contains
    subroutine run_de_lammps(step)
        use parameters, only : population, pstruct, pool
        use init_struct, only : find_dis
        implicit none
        integer(i4b), intent(in) :: step
        integer(i4b) :: i,j,k
        real(dp) :: tmp, cr = 0.9, f = 0.5
        integer(i4b) :: a,b,c, a1, a2
        real(dp) :: p1(3), p2(3), p3(3), y(3), lat(3,3), dis, atom_dis
        logical :: f1
        do i = 1, population
            do
                call random_number(tmp)
                a = integer(tmp * population + 1)
                if(a /= i) exit
            end do
            do
                call random_number(tmp)
                b = integer(tmp * population + 1)
                if(b /= i) exit
            end do
            do
                call random_number(tmp)
                c = integer(tmp * population + 1)
                if(c /= i) exit
            end do
            n = pstruct(i) % natom
            lat = pstruct(i) % lat
            do j = 1, n
                p1(:) = pool(a) % pos(:, j)
                p2(:) = pool(b) % pos(:, j)
                p3(:) = pool(c) % pos(:, j)
                y(:) = p1(:) + f * (p2(:) - p3(:))
                do k = 1, 3
                    y(k) = y(k) - lat(k, k) * real(floor(y(k) / lat(k, k)))
                end do
                f1 = .true.
                do k = 1, n
                    if(k == j) continue
                    find_dis(lat, pstruct(i) % pos(:, k), y(:), dis)
                    a1 = pstruct(i) % ptype(k)
                    a2 = pstruct(i) % ptype(j)
                    if(dis < atom_dis(a1, a2)) then
                        f1 = .false.
                    end if
                end do
                if(f1) then
                    pstruct(i) % pos(:, j) = y(:)
                end if
            end do
        end do
    end subroutine run_de
end module differencial_evolution

