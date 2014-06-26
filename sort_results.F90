!     
! File:   sort_results.F90
! Author: zyy
!
! Created on 2013年4月24日, 下午2:00
!
module sort
    contains
    subroutine sort_results()
        use kinds
        use parameters, only : struct_info, pstruct, pool, population, volumn
        use lattice_mod, only : get_volumn
        implicit none
        type(struct_info) :: tmp
        real(dp) :: min_e
        real(dp) :: lattice(3,3)
        integer(i4b) :: i,j,k
        do i = 1, population
            min_e = pool(i) % energy
            k = i
            do j = i + 1, population
                if(pool(j) % energy < min_e) then
                    min_e = pool(j) % energy
                    k = j
                end if
            end do
            if(k /= i) then
                tmp = pool(i)
                pool(i) = pool(k)
                pool(k) = tmp
                tmp = pstruct(i)
                pstruct(i) = pstruct(k)
                pstruct(k) = tmp
            end if
        end do
        do i = 1, population
            write(1224, "(1X, I5, A8, F10.3)"), i, "energy= ", pool(i) % energy
        end do
        lattice = pool(1) % lat
        call get_volumn(lattice, volumn)
        write(1224, "(1X, A8, F10.3)"), "volumn = ", volumn
    end subroutine sort_results
    
    subroutine sort_results_mode()
        use kinds
        use constants
        use init_struct_all, only : init_struct_tot
        use parameters, only : struct_info, pstruct, pool, population
        use parameters, only : HSE, energy_cut, gap_cut, LDA_population, HSE_population
        implicit none
        type(struct_info) :: tmp
        type(struct_info) :: pstruct_HSE(max_struct)
        integer(i4b) :: flag(max_struct), fcnt
        integer(i4b) :: i,j,k
        real(dp) :: fts, inf
        logical :: f1
        
        inf = 1000
        call sort_results()
        flag = 0
        do i = 1, population
            pool(i) % frontier = 0
        end do
        fcnt = 0
        do
            fcnt = fcnt + 1
            if(sum(flag) == population) exit
            do i = 1, population
                if(pool(i) % frontier == 0) then
                    pool(i) % frontier = fcnt
                    fts = pool(i) % hardness
                    flag(i) = 1
                    exit
                end if
            end do
            do i = 1, population
                if(pool(i) % frontier == 0) then
                    if(pool(i) % hardness < fts) then
                        fts = pool(i) % hardness
                        pool(i) % frontier = fcnt
                        flag(i) = 1
                    end if
                end if
            end do
        end do
        do i = 1, population
            fcnt = pool(i) % frontier
            k = i
            do j = i + 1, population
                if(pool(j) % frontier < fcnt) then
                    fcnt = pool(j) % frontier
                    k = j
                end if
            end do
            if(k /= i) then
                tmp = pool(i)
                pool(i) = pool(k)
                pool(k) = tmp
                tmp = pstruct(i)
                pstruct(i) = pstruct(k)
                pstruct(k) = pstruct(i)
            end if
        end do
        do i = 1, population
            write(1224, "(1X, 2I5, A20, 2F10.3)"), i, pool(i) % frontier, "pool_fitness= ", pool(i) % energy, pool(i) % hardness
            write(1224, "(1X, 2I5, A20, 2F10.3)"), i, pool(i) % frontier, "pstruct_fitness= ", &
            & pstruct(i) % energy, pstruct(i) % hardness
            write(*, "(1X, 2I5, A20, 2F10.3)"), i, pool(i) % frontier, "pool_fitness= ", pool(i) % energy, pool(i) % hardness
            write(*, "(1X, 2I5, A20, 2F10.3)"), i, pool(i) % frontier, "pstruct_fitness= ", &
            & pstruct(i) % energy, pstruct(i) % hardness
        end do
        if(HSE) then
            do i = 1, HSE_population
                f1 = .true.
                do j = 1, LDA_population
                    if((pstruct(j) % energy / pstruct(j) % natom) < energy_cut .and. pstruct(j) % hardness < gap_cut) then
                        write(*, *), "HSE select: ", j
                        write(*, *), pstruct(j) % energy, pstruct(j) % hardness
                        f1 = .false.
                        pstruct_HSE(i) = pstruct(j)
                        pstruct(j) % energy = inf
                        pstruct(j) % hardness = inf
                        write(*, *), pstruct(j) % energy, pstruct(j) % hardness
                        exit
                    end if
                end do
                if(f1) then
                    write(*, *), "In HSE, generate new structure: ", i
                    call init_struct_tot(i)
                    pstruct_HSE(i) = pstruct(i)
                    pstruct(i) % energy = inf
                    pstruct(i) % hardness = inf
                end if
            end do
            do i = 1, HSE_population
                pstruct(i) = pstruct_HSE(i)
            end do
        end if
    end subroutine sort_results_mode
end module sort


module sort_atom_idx
    use kinds
    use constants, only : max_atom
    implicit none
    real(dp) :: INF = 100000000.0
    real(dp) :: eps = 1e-12
    integer(i4b), dimension(max_atom) :: mx, my
    real(dp), dimension(max_atom) :: lx, ly
    logical, dimension(max_atom) :: sx, sy
    integer(i4b) :: nx, ny
    real(dp) :: g(max_atom, max_atom)
    contains
    subroutine kmatch(ptl1, ptl2)
        use parameters, only : struct_info
        implicit none
        type(struct_info), intent(in) :: ptl1
        type(struct_info), intent(inout) :: ptl2
        type(struct_info) :: tmp_ptl
        integer(i4b) :: i,j,u
        real(dp) :: ex, dis
        nx = ptl1 % natom
        ny = ptl2 % natom
        ! build graph
        do i = 1, nx
            do j = 1, ny
                call re_dis(ptl1 % pos(:, i), ptl2 % pos(:, j), dis)
                g(i, j) = - dis
            end do
        end do
        ! KM
        ly = 0.0
        mx = -1
        my = -1
        do i = 1, nx
            lx(i) = - INF
            do j = 1, ny
                lx(i) = max(lx(i), g(i, j))
            end do
        end do
        do u = 1, nx
            if(mx(u) == -1) then
                sx = .false.
                sy = .false.
                do
                    if(path(u) .eqv. .true.) exit
                    if(path(u) .eqv. .false.) then
                        ex = INF
                        do i = 1, nx
                            if(sx(i) .eqv. .true.) then
                                do j = 1, ny
                                    if(sy(j) .eqv. .false.) then
                                        if(lx(i) + ly(j) - g(i, j) < ex) then
                                            ex = lx(i) + ly(j) - g(i, j)
                                        end if
                                    end if
                                end do
                            end if
                        end do
                        do i = 1, nx
                            if(sx(i) .eqv. .true.) then
                                lx(i) = lx(i) - ex
                                sx(i) = .false.
                            end if
                        end do
                        do j = 1, ny
                            if(sy(j) .eqv. .true.) then
                                ly(j) = ly(j) - ex
                                sy(j) = .false.
                            end if
                        end do
                    end if
                end do
            end if
        end do
        !write(*, *) (mx(i), i = 1, nx)
        tmp_ptl = ptl2
        do i = 1, nx
            tmp_ptl % pos(:, i) = ptl2 % pos(:, mx(i))
        end do
    end subroutine kmatch
    
    recursive function path(u) result (ans)
        implicit none
        integer(i4b), intent(in) :: u
        logical :: ans
        integer(i4b) :: v
        sx(u) = .true.
        do v = 1, ny
            if((g(u, v) < lx(u) + ly(v) + eps) .and. (g(u, v) > lx(u) + ly(v) - eps) .and. (sy(v) .eqv. .false.)) then
                sy(v) = .true.
                if(my(v) == -1 .or. (path(my(v)) .eqv. .true.)) then
                    mx(u) = v
                    my(v) = u
                    ans = .true.
                    return
                else
                    ans = .false.
                end if
            else
                ans = .false.
            end if
        end do
    end function path
    
    subroutine re_dis(pos1, pos2, dist)
        use kinds
        implicit none
        real(dp), intent(in) :: pos1(3), pos2(3)
        real(dp), intent(out) :: dist
        real(dp) :: cnt
        integer(i4b) :: i,j
        cnt = 0.0
        do i = 1, 3
            cnt = cnt + (pos1(i) - pos2(i)) ** 2
        end do
        dist = sqrt(cnt)
    end subroutine re_dis
end module sort_atom_idx

