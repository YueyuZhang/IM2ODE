!     
! File:   init_struct.F90
! Author: zyy
!
! Created on 2013年4月21日, 上午1:25
!

module init_struct
    contains
    subroutine init_struct_nosym()
        use kinds
        use parameters, only : struct_info
        use parameters, only : pstruct, population, num_species, num_ele, atom_dis
        implicit none
        integer(i4b) :: i,j,k,i1
        do i = 1, population
            pstruct(i) % ntyp = num_species
            pstruct(i) % nelement = num_ele
            pstruct(i) % natom = sum(num_ele)
            i1 = 0
            do j = 1, num_species
                do k = 1, num_ele(j)
                    i1 = i1 + 1
                    pstruct(i) % ptype(i1) = j
                end do
            end do
            call fill_atom(i)
        end do
    end subroutine init_struct_nosym
    
    subroutine fill_atom(i)
        use kinds
        use parameters, only : pstruct, atom_dis
        implicit none
        integer(i4b), intent(in) :: i
        integer(i4b) :: a1, a2, i1, j, k
        real(dp) :: pos_car(3), pos_dir(3), pos_tmp(3), lattice(3,3), distance
        logical :: flag
        lattice = pstruct(i) % lat
        do j = 1, pstruct(i) % natom
            flag = .true.
            i1 = 1
            do
                i1 = i1 + 1
                do k = 1, 3
                    call random_number(pos_dir(k))
                end do
                pos_car = pos_dir
                call dir2car(lattice, pos_dir, pos_car)
                do k = 1, j - 1
                    pos_tmp = pstruct(i) % pos(:, k)
                    call find_dis(lattice, pos_car, pos_tmp, distance)
                    a1 = pstruct(i) % ptype(k)
                    a2 = pstruct(i) % ptype(j)
                    if(distance < atom_dis(a1,a2)) then
                        flag = .false.
                    end if
                end do
                if(flag .or. i1 > 1000) exit
            end do
            pstruct(i) % pos(:,j) = pos_car
        end do
    end subroutine fill_atom
    
    subroutine dir2car(lat, pos_old, pos_new)
        use kinds
        implicit none
        real(dp), intent(in) :: lat(3, 3)
        real(dp), intent(in) :: pos_old(3)
        real(dp), intent(out) :: pos_new(3)
        real(dp) :: posx(3), cnt
        integer(i4b) :: i,j,k
        posx = pos_old
        do i = 1, 3
            cnt = 0.0
            do j = 1, 3
                cnt = cnt + pos_old(j) * lat(j, i)
            end do
            posx(i) = cnt
        end do
        pos_new = posx
    end subroutine dir2car
                
    subroutine find_dis(lat, pos1, pos2, dist)
        use kinds
        implicit none
        real(dp), intent(in) :: lat(3,3), pos1(3), pos2(3)
        real(dp), intent(out) :: dist
        real(dp) :: cnt
        integer(i4b) :: i,j
        cnt = 0.0
        do i = 1, 3
            cnt = cnt + (pos1(i) - pos2(i)) ** 2
        end do
        dist = sqrt(cnt)
    end subroutine find_dis
end module init_struct
        
module init_struct_spg
    use kinds
    integer(i4b) :: n_tot, length(32), ans(32), res(32, 100000), cnt, sum_atom, spg_index
    
    contains
    
    subroutine fill_atom_spg(tag, flag)
        use kinds
        use parameters, only : pstruct, atom_dis
        use init_struct
        use spacegroup_data, only : spg_opr, spg_ctl, spg_x, spg_y, spg_z, spg_add, spg_pro
        implicit none
        integer(i4b), intent(in) :: tag
        logical, intent(out) :: flag
        integer(i4b) :: i,j,k,p,q,cho_res(100000), res_idx, iatom, patom, pt, ptyp, nf2
        logical :: f1, f2, f3
        real(dp) :: a,b,c,x(192), y(192), z(192)
        real(dp) :: pos_car(3), pos_dir(3), pos_tmp(3), pos_tmp_car(3), lattice(3,3), distance
        integer(i4b) :: type1, type2, natom
        integer(i4b) :: dx(27), dy(27), dz(27), dnum
        dnum = 0
        do i = -1, 1
            do j = -1, 1
                do k = -1, 1
                    dnum = dnum + 1
                    dx(dnum) = i
                    dy(dnum) = j
                    dz(dnum) = k
                end do
            end do
        end do
        
        spg_index = pstruct(tag) % spg_idx
        lattice = pstruct(tag) % lat
        ptyp = pstruct(tag) % ntyp
        patom = 0
        do pt = 1, ptyp
            natom = pstruct(tag) % nelement(pt)
            n_tot = spg_opr(spg_index)
            length(:) = spg_ctl(:, spg_index)
            res = 0
            ans = 0
            cnt = 0
            sum_atom = pstruct(tag) % nelement(pt)
            call find(1,0)
            write(*, *)"spg = ", spg_index
            if(cnt == 0) then
                write(*, *) "no permutation"
                flag = .false.
                return
            end if
            !do i = 1, cnt
            !    write(*, *) (res(j, i), j = 1, n_tot)
            !end do
            cho_res = 0
            flag = .true.
            f1 = .false.
            do
                if(sum(cho_res) == cnt) then
                    flag = .false.
                    return
                end if
                ! find permutation
                do
                    call random_number(a)
                    res_idx = floor(cnt * a + 1)
                    if(cho_res(res_idx) == 0) then
                        cho_res(res_idx) = 1
                        f1 = .true.
                    end if
                    if(f1) exit
                    if(sum(cho_res) > floor(cnt * 0.7)) then
                        do i = 1, cnt
                            if(cho_res(i) == 0) then
                                res_idx = i
                                cho_res(res_idx) = 1
                                f1 = .true.
                            end if
                        end do
                    end if
                    if(f1) exit
                end do
                !
                ans(:) = res(:, res_idx)
                iatom = 0
                do i = 1, n_tot
                    do j = 1, ans(i)
                        f2 = .true.
                        nf2 = 0
                        do
                            f2 = .true.
                            nf2 = nf2 + 1
                            call random_number(a)
                            call random_number(b)
                            call random_number(c)
                            do k = 1, spg_ctl(i, spg_index)
                                x(k) = spg_x(1,k,i,spg_index)*a + spg_x(2,k,i,spg_index)*b + spg_x(3,k,i,spg_index)*c
                                y(k) = spg_y(1,k,i,spg_index)*a + spg_y(2,k,i,spg_index)*b + spg_y(3,k,i,spg_index)*c
                                z(k) = spg_z(1,k,i,spg_index)*a + spg_z(2,k,i,spg_index)*b + spg_z(3,k,i,spg_index)*c
                                !x(k) = spg_pro(1,k,i,spg_index) * x(k)
                                x(k) = spg_add(1,k,i,spg_index) + x(k)
                                !y(k) = spg_pro(2,k,i,spg_index) * y(k)
                                y(k) = spg_add(2,k,i,spg_index) + y(k)
                                !z(k) = spg_pro(3,k,i,spg_index) * z(k)
                                z(k) = spg_add(3,k,i,spg_index) + z(k)
                                x(k) = x(k) - floor(x(k))
                                y(k) = y(k) - floor(y(k))
                                z(k) = z(k) - floor(z(k))
                            end do
                            do k = 1, spg_ctl(i, spg_index)
                                do p = 1, patom + iatom
                                    pos_tmp(:) = pstruct(tag) % pos(:, p)
                                    do q = 1, 27
                                        pos_dir(1) = x(k) + dx(q)
                                        pos_dir(2) = y(k) + dy(q)
                                        pos_dir(3) = z(k) + dz(q)
                                        call dir2car(lattice, pos_tmp, pos_tmp_car)
                                        call dir2car(lattice, pos_dir, pos_car)
                                        call find_dis(lattice, pos_car, pos_tmp_car, distance)
                                        type1 = pstruct(tag) % ptype(p)
                                        type2 = pstruct(tag) % ptype(patom + iatom + k)
                                        if(distance < atom_dis(type1, type2)) then
                                            f2 = .false.
                                        end if
                                    end do
                                end do
                                do p = 1, k-1
                                    pos_tmp(1) = x(p)
                                    pos_tmp(2) = y(p)
                                    pos_tmp(3) = z(p)
                                    do q = 1, 27
                                        pos_dir(1) = x(k) + dx(q)
                                        pos_dir(2) = y(k) + dy(q)
                                        pos_dir(3) = z(k) + dz(q)
                                        call dir2car(lattice, pos_tmp, pos_tmp_car)
                                        call dir2car(lattice, pos_dir, pos_car)
                                        call find_dis(lattice, pos_car, pos_tmp_car, distance)
                                        type1 = pstruct(tag) % ptype(patom + iatom + p)
                                        type2 = pstruct(tag) % ptype(patom + iatom + k)
                                        if(distance < atom_dis(type1, type2)) then
                                            f2 = .false.
                                        end if
                                    end do
                                end do
                            end do
                            if(f2) then
                                do k = 1, spg_ctl(i, spg_index)
                                    pstruct(tag) % pos(1, patom + iatom + k) = x(k)
                                    pstruct(tag) % pos(2, patom + iatom + k) = y(k)
                                    pstruct(tag) % pos(3, patom + iatom + k) = z(k)
                                end do
                                iatom = iatom + spg_ctl(i, spg_index)
                            end if
                            if(f2 .or. nf2 > 5000) exit
                        end do
                    end do
                end do
                if(iatom == natom) then
                    exit
                end if
            end do
            patom = patom + iatom
            if(patom == pstruct(tag) % natom) then
                flag = .true.
                return
            end if
        end do
    end subroutine fill_atom_spg
    
    recursive subroutine find(k, add)
        use spacegroup_data, only : spg_pro
        implicit none
        integer(i4b), intent(in) :: k, add
        integer(i4b) :: i,j,m, i1
        real(dp) :: tmp(32)
        real(dp) :: t1, eps = 1e-8
        if(cnt > 100000) return
        if(k > n_tot + 1) return
        if(add > sum_atom) return
        if(k == n_tot + 1 .and. add == sum_atom) then
            cnt = cnt + 1
            res(:, cnt) = ans(:)
            return
        end if
        j = floor(real(sum_atom - add) / real(length(k))) + 1
        m = length(k)
        t1 = 0.0
        do i = 1, m
            do i1 = 1, 32
                tmp(i1) = abs(spg_pro(i1,i,k,spg_index))
            end do
            t1 = t1 + sum(tmp)
        end do
        if(t1 < eps) then
            j = 1
        end if
        do i = 0, j
            ans(k) = i
            call find(k + 1, add + length(k) * i)
        end do
    end subroutine find
    
    subroutine init_struct_sym(i)
        use kinds
        use parameters, only: struct_info
        use parameters, only: pstruct, population, num_species, num_ele, atom_dis
        use parameters, only: Q2D, vacuum_layer
        use lattice_mod, only: init_lat_sym
        implicit none
        integer(i4b), intent(in) :: i
        integer(i4b) :: j,k,i1
        real(dp) :: c1, c2, ratio
        logical :: f1
        !do i = 1, population
        write(*, *) "pop = ", i
        pstruct(i) % ntyp = num_species
        pstruct(i) % nelement = num_ele
        pstruct(i) % natom = sum(num_ele)
        i1 = 0
        do j = 1, num_species
            do k = 1, num_ele(j)
                i1 = i1 + 1
                pstruct(i) % ptype(i1) = j
            end do
        end do
        f1 = .true.
        do
            call init_lat_sym(i)
            !write(*, *) "end init lat"
            call fill_atom_spg(i, f1)
            !write(*, *) "end fill atom"
            if(f1) exit
        end do
        !end do
        if(Q2D) then
            c1 = pstruct(i) % lat(3,3)
            c2 = vacuum_layer + c1
            ratio = c1 / c2
            do j = 1, pstruct(i) % natom
                pstruct(i) % pos(3, j) = pstruct(i) % pos(3, j) * ratio
            end do
            pstruct(i) % lat(3, 3) = c2
        end if
    end subroutine init_struct_sym
    
end module init_struct_spg

