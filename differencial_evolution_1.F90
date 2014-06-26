!     
! File:   differencial_evolution.F90
! Author: zyy
!
! Created on 2013年4月23日, 下午9:47
!

module differencial_evolution
    contains
    subroutine run_de_lammps(step)
        use kinds
        use lattice_mod, only : init_lat_nosym
        use init_struct, only : fill_atom
        use parameters, only : population, pstruct, pool, atom_dis, de_ratio
        use sort_atom_idx, only : kmatch
        use init_struct
        implicit none
        integer(i4b), intent(in) :: step
        integer(i4b) :: i,j,k
        real(dp) :: tmp, cr = 0.9, f = 0.5
        integer(i4b) :: a,b,c, a1, a2, n
        real(dp) :: p1(3), p2(3), p3(3), y(3), lat(3,3), dis
        logical :: f1
        do i = 1, population
            if(i > floor(population * de_ratio)) then
                call init_lat_nosym(i)
                call fill_atom(i)
                cycle
            end if
            do
                call random_number(tmp)
                a = floor(tmp * population * de_ratio + 1.0)
                if(a /= i) exit
            end do
            do
                call random_number(tmp)
                b = floor(tmp * population * de_ratio + 1.0)
                if(b /= i .and. b /= a) exit
            end do
            do
                call random_number(tmp)
                c = floor(tmp * population * de_ratio + 1.0)
                if(c /= i .and. c /= a .and. c /= b) exit
            end do
            call kmatch(pstruct(i), pool(a))
            call kmatch(pstruct(i), pool(b))
            call kmatch(pstruct(i), pool(c))
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
                    if(k == j) cycle
                    call find_dis(lat, pstruct(i) % pos(:, k), y(:), dis)
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
    end subroutine run_de_lammps
    
    subroutine run_de_vasp(step)
        use kinds
        use parameters, only : pstruct, pool, population, de_ratio
        use parameters, only : Q2D, vacuum_layer, cluster
        use sort, only : sort_results
        use init_struct_all, only : init_struct_tot
        use sort_atom_idx, only : kmatch
        implicit none
        integer(i4b), intent(in) :: step
        integer(i4b) :: i,j,k, ratio_cut, a,b
        real(dp) :: tmp, y(3), p1(3), p2(3), p3(3), pbest(3)
        real(dp) :: lambd = 0.2, F = 0.2
        real(dp) :: c1, c2, c3, ratio, max_z, min_z
        logical :: flag
        write(*, *) "start de_operation"
        do i = 1, population
            write(*, *) "energy now : old", pstruct(i) % energy, pool(i) % energy
            if (pstruct(i) % energy < pool(i) % energy) then
                write(*, *) "renew: ", i
                pool(i) = pstruct(i)
            end if
        end do
        call sort_results()
        ratio_cut = floor(de_ratio * population)
        
        do i = 1, population
            write(*, *) "de_struct: ", i
            write(*, *) "atom_type: ", (pstruct(i) % ptype(k), k = 1, pstruct(i) % natom)
            if(i > ratio_cut) then
                call init_struct_tot(i)
                cycle
            end if
            do
                call random_number(tmp)
                a = floor(tmp * ratio_cut + 1)
                if(a /= i) exit
            end do
            do
                call random_number(tmp)
                b = floor(tmp * ratio_cut + 1)
                if(b /= i .and. b /= a) exit
            end do
            !call kmatch(pool(1), pstruct(i))
            !call kmatch(pool(1), pool(a))
            !call kmatch(pool(1), pool(b))
            do j = 1, pstruct(i) % natom
                pbest(:) = pool(1) % pos(:, j)
                p1(:) = pool(i) % pos(:, j)
                p2(:) = pool(a) % pos(:, j)
                p3(:) = pool(b) % pos(:, j)
                y(:) = lambd * pbest(:) + (1-lambd) * p1(:) + F * (p2(:) - p3(:))
                y(:) = y(:) - floor(y(:))
                pstruct(i) % pos(:, j) = y(:)
            end do
            call check(i, flag)
            if (flag .eqv. .false.) then
                call init_struct_tot(i)                
            end if
            ! for quasi-2D case
            ! to keep the vacuum_layer after DE operations
            if(Q2D) then
                min_z = 1.0
                max_z = 0.0
                do j = 1, pstruct(i) % natom
                    if(min_z > pstruct(i) % pos(3, j)) min_z = pstruct(i) % pos(3, j)
                    if(max_z < pstruct(i) % pos(3, j)) max_z = pstruct(i) % pos(3, j)
                end do
                c1 = (1.0 - (max_z - min_z)) * pstruct(i) % lat(3, 3)
                c2 = (max_z - min_z) * pstruct(i) % lat(3, 3) + vacuum_layer
                ratio = c1 / c2
                do j = 1, pstruct(i) % natom
                    pstruct(i) % pos(3, j) = pstruct(i) % pos(3, j) * ratio
                end do
                pstruct(i) % lat(3, 3) = c2
            end if
        end do
    end subroutine run_de_vasp
    
    subroutine check(tag, flag)
        use kinds
        use parameters, only : pstruct, atom_dis
        use init_struct, only : dir2car, find_dis
        implicit none
        integer(i4b), intent(in) :: tag
        logical, intent(out) :: flag
        integer(i4b) :: i,j,k,a1,a2,n
        real(dp) :: dir0(3), dir1(3), dir2(3), car1(3), car2(3), lattice(3,3)
        real(dp) :: distance
        integer(i4b) :: dpos(3, 27), dnum
        n = pstruct(tag) % natom
        lattice = pstruct(tag) % lat
        
        dnum = 0
        do i = -1, 1
            do j = -1, 1
                do k = -1, 1
                    dnum = dnum + 1
                    dpos(1,dnum) = i
                    dpos(2,dnum) = j
                    dpos(3,dnum) = k
                end do
            end do
        end do
        
        
        flag = .true.
        do i = 1, n
            do j = 1, i - 1
                dir0 = pstruct(tag) % pos(:, i)
                dir2 = pstruct(tag) % pos(:, j)
                do k = 1, 27
                    dir1(:) = dpos(:, k) + dir0(:)
                    call dir2car(lattice, dir1, car1)
                    call dir2car(lattice, dir2, car2)
                    call find_dis(lattice, car1, car2, distance)
                    a1 = pstruct(tag) % ptype(i)
                    a2 = pstruct(tag) % ptype(j)
                    if(distance < atom_dis(a1, a2) * 0.9) then
                        flag = .false.
                        return
                    end if
                end do
            end do
        end do
    end subroutine check
    
    subroutine run_mode_vasp(step)
        use kinds
        use parameters, only : pstruct, pool, population, de_ratio
        use sort, only : sort_results_mode
        use init_struct_spg, only : init_struct_sym
        implicit none
        integer(i4b), intent(in) :: step
        integer(i4b) :: i,j,k,ratio_cut, a, b, c, nf
        real(dp) :: tmp, y(3), p1(3), p2(3), p3(3), pbest(3)
        real(dp) :: lambd = 0.1, F = 0.2
        logical :: flag
        do i = 1, population
            if(pstruct(i) % energy < pool(i) % energy .and. pstruct(i) % hardness < pool(i) % hardness) then
                pool(i) = pstruct(i)
            end if
            call random_number(tmp)
            if(pstruct(i) % energy < pool(i) % energy .or. pstruct(i) % hardness < pool(i) % hardness) then
                if(tmp < 0.5) then
                    pool(i) = pstruct(i)
                end if
            end if
        end do
        write(*, *) "start sort"
        call sort_results_mode()
        write(*, *) "end sort"
        ratio_cut = floor(de_ratio * population)
        nf = 0
        do i = 1, population
            if(pool(i) % frontier == 1) then
                nf = nf + 1
            end if
        end do
        write(*, *) "number of struct in 1st front: ", nf
        do i = 1, population
            write(*, *) "de on: ", i
            if(pool(i) % energy > 100 .or. pool(i) % hardness > 100) then
                write(*, *) "de new struct1: ", i
                call init_struct_sym(i)
                cycle
            end if
            if(i > ratio_cut) then
                write(*, *) "de new struct2: ", i
                call init_struct_sym(i)
                cycle
            end if
            write(*, *) "de operation: ", i
            do
                call random_number(tmp)
                a = floor(tmp * ratio_cut + 1)
                if(a /= i) exit
            end do
            do
                call random_number(tmp)
                b = floor(tmp * ratio_cut + 1)
                if(b /= i .and. b /= a) exit
            end do
            if(pool(i) % frontier /= 1) then
                call random_number(tmp)
                c = floor(tmp * nf + 1)
            end if
            if(pool(i) % frontier == 1) then
                do j = 1, pstruct(i) % natom
                    p1(:) = pool(i) % pos(:, j)
                    p2(:) = pool(a) % pos(:, j)
                    p3(:) = pool(b) % pos(:, j)
                    y(:) = p1(:) + F * (p2(:) - p3(:))
                    y(:) = y(:) - floor(y(:))
                    pstruct(i) % pos(:, j) = y(:)
                end do
            else
                do j = 1, pstruct(i) % natom
                    pbest(:) = pool(c) % pos(:, j)
                    p1(:) = pool(i) % pos(:, j)
                    p2(:) = pool(a) % pos(:, j)
                    p3(:) = pool(b) % pos(:, j)
                    y(:) = lambd * pbest(:) + (1 - lambd) * p1(:) + F * (p2(:) - p3(:))
                    y(:) = y(:) - floor(y(:))
                    pstruct(i) % pos(:, j) = y(:)
                end do
            end if
            call check(i,flag)
            if(flag .eqv. .false.) then
                call init_struct_sym(i)
            end if
        end do
    end subroutine run_mode_vasp
end module differencial_evolution

