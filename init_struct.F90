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
        use constants, only : pi
        use parameters, only : pstruct, atom_dis
        use parameters, only : cluster, init_radius, cluster_ctr_x, cluster_ctr_y, cluster_ctr_z
        use parameters, only : model_ball, model_shell, model_plate
        use parameters, only : shell_radius_in, shell_radius_out
        use parameters, only : shell_ctr_x, shell_ctr_y, shell_ctr_z
        use parameters, only : plate_radius, plate_height
        use parameters, only : plate_ctr_x, plate_ctr_y, plate_ctr_z
        implicit none
        integer(i4b), intent(in) :: i
        integer(i4b) :: a1, a2, i1, j, k
        real(dp) :: pos_car(3), pos_dir(3), pos_tmp(3), pos_tmp_car(3), lattice(3,3), pos_bak(3)
        real(dp) :: distance, distance_tot, distance_max
        real(dp) :: r, theta, phi, height
        real(dp) :: tmp
        logical :: flag
        distance_max = 0.0
        lattice = pstruct(i) % lat
        call random_number(tmp)
        do j = 1, pstruct(i) % natom
            flag = .true.
            i1 = 1
            distance_max = 0.0
            do
                flag = .true.
                i1 = i1 + 1
                do k = 1, 3
                    call random_number(pos_dir(k))
                end do
                if(cluster) then
                    if(tmp < model_ball) then
                        call random_number(r)
                        call random_number(theta)
                        call random_number(phi)
                        r = r * init_radius
                        theta = theta * 2.0 * pi
                        phi = phi * pi
                        pos_car(1) = r * sin(phi) * cos(theta)
                        pos_car(2) = r * sin(phi) * sin(theta)
                        pos_car(3) = r * cos(phi)
                        call car2dir(lattice, pos_car, pos_dir)
                        pos_dir(1) = pos_dir(1) + cluster_ctr_x
                        pos_dir(2) = pos_dir(2) + cluster_ctr_y
                        pos_dir(3) = pos_dir(3) + cluster_ctr_z
                    else if (tmp < model_ball + model_shell) then
                        call random_number(r)
                        call random_number(theta)
                        call random_number(phi)
                        r = r * (shell_radius_out - shell_radius_in) + shell_radius_in
                        theta = theta * 2.0 * pi
                        phi = phi * pi
                        pos_car(1) = r * sin(phi) * cos(theta)
                        pos_car(2) = r * sin(phi) * sin(theta)
                        pos_car(3) = r * cos(phi)
                        call car2dir(lattice, pos_car, pos_dir)
                        pos_dir(1) = pos_dir(1) + shell_ctr_x
                        pos_dir(2) = pos_dir(2) + shell_ctr_y
                        pos_dir(3) = pos_dir(3) + shell_ctr_z
                    else
                        call random_number(r)
                        call random_number(theta)
                        call random_number(height)
                        r = r * plate_radius
                        theta = theta * 2.0 * pi
                        height = (height - 0.5) * plate_height
                        pos_car(1) = r * sin(theta)
                        pos_car(2) = r * cos(theta)
                        pos_car(3) = height
                        call car2dir(lattice, pos_car, pos_dir)
                        pos_dir(1) = pos_dir(1) + plate_ctr_x
                        pos_dir(2) = pos_dir(2) + plate_ctr_y
                        pos_dir(3) = pos_dir(3) + plate_ctr_z
                    end if
                end if
                pos_car = pos_dir
                call dir2car(lattice, pos_dir, pos_car)
                distance_tot = 1e5
                do k = 1, j - 1
                    pos_tmp = pstruct(i) % pos(:, k)
                    call dir2car(lattice, pos_tmp, pos_tmp_car)
                    call find_dis(lattice, pos_car, pos_tmp_car, distance)
                    a1 = pstruct(i) % ptype(k)
                    a2 = pstruct(i) % ptype(j)
                    if(distance < atom_dis(a1,a2)) then
                        flag = .false.
                        if(.not. flag .and. distance < distance_tot) distance_tot = distance
                    end if    
                end do
                if(.not. flag .and. distance_tot > distance_max) then
                    distance_max = distance_tot
                    pos_bak = pos_dir
                end if
                if(flag .or. i1 > 50000) then
                    !if(i1 > 5000) write(*, *) "too small radius"
                    if(i1 > 50000) then
                        pos_dir = pos_bak
                    end if
                    exit
                end if
            end do
            pstruct(i) % pos(:,j) = pos_dir
        end do
    end subroutine fill_atom
    
    subroutine fill_atom_cluster_substrate(i)
        use kinds
        use constants, only : pi
        use parameters, only : pstruct, atom_dis
        use parameters, only : cluster, init_radius, cluster_ctr_x, cluster_ctr_y, cluster_ctr_z
        use parameters, only : cluster_substrate, substrate
        use parameters, only : model_ball, model_shell, model_plate
        use parameters, only : shell_radius_in, shell_radius_out
        use parameters, only : shell_ctr_x, shell_ctr_y, shell_ctr_z
        use parameters, only : plate_radius, plate_height
        use parameters, only : plate_ctr_x, plate_ctr_y, plate_ctr_z
        use parameters, only : selective_dynamics
        use parameters, only : find_defect
        use parameters, only : model_surface, surface_height
        implicit none
        integer(i4b), intent(in) :: i
        integer(i4b) :: a1, a2, i1, j, k
        integer(i4b) :: dx(27), dy(27), dz(27), dnum, q
        real(dp) :: pos_car(3), pos_dir(3), pos_tmp(3), pos_tmp_car(3), lattice(3,3), pos_bak(3)
        real(dp) :: distance, distance_tot, distance_max, top_substrate
        real(dp) :: r, theta, phi, height
        real(dp) :: tmp
        logical :: flag
        distance_max = 0.0
        top_substrate = 0.0
        lattice = pstruct(i) % lat
        write(*, *) "fill substrate"
        if(selective_dynamics) pstruct(i) % SelectiveDynamics = substrate % SelectiveDynamics
        do j = 1, pstruct(i) % natom
            if(substrate % ptype(j) /= -1) then
                pstruct(i) % pos(:, j) = substrate % pos(:, j)
                write(*, *) pstruct(i) % pos(:, j)
                if(substrate % pos(3, j) > top_substrate) then
                    top_substrate = substrate % pos(3, j)
                end if
            end if
        end do
        
        dnum = 0
        do k = -1, 1
            do j = -1, 1
                dnum = dnum + 1
                dx(dnum) = k
                dy(dnum) = j
                dz(dnum) = 0
            end do
        end do
        
        call random_number(tmp)
        do j = 1, pstruct(i) % natom
            if(substrate % ptype(j) == -1) then
                flag = .true.
                i1 = 1
                distance_max = 0.0
                do
                    flag = .true.
                    i1 = i1 + 1
                    
                    if(model_surface) then
                        ! surface model
                        call random_number(r)
                        call random_number(theta)
                        call random_number(phi)
                        pos_dir(1) = r
                        pos_dir(2) = theta
                        pos_dir(3) = phi * surface_height / lattice(3, 3) + top_substrate
                    else if(tmp < model_ball) then
                        ! ball model
                        call random_number(r)
                        call random_number(theta)
                        call random_number(phi)
                        r = r * init_radius
                        theta = theta * 2.0 * pi
                        phi = phi * pi
                        pos_car(1) = r * sin(phi) * cos(theta)
                        pos_car(2) = r * sin(phi) * sin(theta)
                        pos_car(3) = r * cos(phi)
                        call car2dir(lattice, pos_car, pos_dir)
                        pos_dir(1) = pos_dir(1) + cluster_ctr_x
                        pos_dir(2) = pos_dir(2) + cluster_ctr_y
                        pos_dir(3) = pos_dir(3) + cluster_ctr_z
                    else if (tmp < model_ball + model_shell) then
                        ! shell model
                        call random_number(r)
                        call random_number(theta)
                        call random_number(phi)
                        r = r * (shell_radius_out - shell_radius_in) + shell_radius_in
                        theta = theta * 2.0 * pi
                        phi = phi * pi
                        pos_car(1) = r * sin(phi) * cos(theta)
                        pos_car(2) = r * sin(phi) * sin(theta)
                        pos_car(3) = r * cos(phi)
                        call car2dir(lattice, pos_car, pos_dir)
                        pos_dir(1) = pos_dir(1) + shell_ctr_x
                        pos_dir(2) = pos_dir(2) + shell_ctr_y
                        pos_dir(3) = pos_dir(3) + shell_ctr_z
                    else
                        ! plate model
                        call random_number(r)
                        call random_number(theta)
                        call random_number(height)
                        r = r * plate_radius
                        theta = theta * 2.0 * pi
                        height = (height - 0.5) * plate_height
                        pos_car(1) = r * sin(theta)
                        pos_car(2) = r * cos(theta)
                        pos_car(3) = height
                        call car2dir(lattice, pos_car, pos_dir)
                        pos_dir(1) = pos_dir(1) + plate_ctr_x
                        pos_dir(2) = pos_dir(2) + plate_ctr_y
                        pos_dir(3) = pos_dir(3) + plate_ctr_z
                    end if
                    
                    if(pos_dir(3) < top_substrate .and. (.not. find_defect)) then
                        flag = .false.
                    end if
                    pos_car = pos_dir
                    call dir2car(lattice, pos_dir, pos_car)
                    distance_tot = 1e5
                    do k = 1, j - 1
                        pos_tmp = pstruct(i) % pos(:, k)
                        call dir2car(lattice, pos_tmp, pos_tmp_car)
                        call find_dis(lattice, pos_car, pos_tmp_car, distance)
                        a1 = pstruct(i) % ptype(k)
                        a2 = pstruct(i) % ptype(j)
                        if(distance < atom_dis(a1,a2)) then
                            flag = .false.
                            if(.not. flag .and. distance < distance_tot) distance_tot = distance
                        end if
                        if(model_surface) then
                            do q = 1, dnum
                                pos_tmp(1) = pstruct(i) % pos(1, k) + dx(q)
                                pos_tmp(2) = pstruct(i) % pos(2, k) + dy(q)
                                pos_tmp(3) = pstruct(i) % pos(3, k) + dz(q)
                                call dir2car(lattice, pos_tmp, pos_tmp_car)
                                call find_dis(lattice, pos_car, pos_tmp_car, distance)
                                if(distance < atom_dis(a1, a2)) then
                                    flag = .false.
                                    if(.not. flag .and. distance < distance_tot) distance_tot = distance
                                end if
                            end do
                        end if
                    end do
                    do k = j + 1, pstruct(i) % natom
                        if(substrate % ptype(k) /= -1) then
                            pos_tmp = substrate % pos(:, k)
                            call dir2car(lattice, pos_tmp, pos_tmp_car)
                            call find_dis(lattice, pos_car, pos_tmp_car, distance)
                            a1 = pstruct(i) % ptype(k)
                            a2 = pstruct(i) % ptype(j)
                            if(distance < atom_dis(a1,a2)) then
                                flag = .false.
                                if(.not. flag .and. distance < distance_tot) distance_tot = distance
                            end if
                        end if
                    end do
                    if(.not. flag .and. distance_tot > distance_max) then
                        distance_max = distance_tot
                        pos_bak = pos_dir
                    end if
                    if(flag .or. i1 > 50000) then
                        !if(i1 > 5000) write(*, *) "too small radius"
                        if(i1 > 50000) then
                            pos_dir = pos_bak
                        end if
                        exit
                    end if
                end do
                pstruct(i) % pos(:,j) = pos_dir
            end if
        end do
    end subroutine fill_atom_cluster_substrate
            
    
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
    
    subroutine car2dir(lat, car, dir)
        use kinds
        implicit none
        real(dp), intent(in) :: lat(3, 3), car(3)
        real(dp), intent(out) :: dir(3)
        integer(i4b) :: n, i, j, k
        real(dp) :: maxp, tmp, aa(3, 3)
        n = 3
        dir(:) = car(:)
        aa = lat
        do k = 1, n
          maxp = lat(k, k)
          do j = k + 1, n
            aa(k, j) = aa(k, j) / maxp
            do i = k + 1, n
              aa(i, j) = aa(i, j) - aa(i, k) * aa(k, j)
            end do
          end do
          dir(k) = dir(k) / maxp
          do i = k + 1, n
            dir(i) = dir(i) - dir(k) * aa(i, k)
          end do
        end do
        do i = n, 1, -1
          do j = i + 1, n
            dir(i) = dir(i) - aa(i, j) * dir(j)
          end do
        end do
    end subroutine car2dir
                
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
    
    subroutine init_struct_cluster(i)
        use kinds
        use parameters, only: struct_info
        use parameters, only: pstruct, population,num_species, num_ele, atom_dis
        use parameters, only: cluster_substrate
        use lattice_mod, only: init_lat_sym
        implicit none
        integer(i4b), intent(in) :: i
        integer(i4b) :: j,k,i1
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
        call init_lat_sym(i)
        write(*, *) pstruct(i) % lat
        if(cluster_substrate) then
            call fill_atom_cluster_substrate(i)
        else
            call fill_atom(i)
        end if
        do j = 1, pstruct(i) % natom
            write(*, *) pstruct(i) % pos(:, j)
        end do
        write(*, *) "atom_type: ", (pstruct(i) % ptype(k), k = 1, pstruct(i) % natom)
    end subroutine init_struct_cluster
    
end module init_struct
        
module init_struct_spg
    use kinds
    integer(i4b) :: n_tot, length(32), ans(32), res(32, 100000), cnt, sum_atom, spg_index
    
    contains
    
    subroutine fill_atom_spg(tag, flag)
        use kinds
        use parameters, only : pstruct, atom_dis
        use parameters, only : Q2D
        use init_struct
        use spacegroup_data, only : spg_opr, spg_ctl, spg_x, spg_y, spg_z, spg_add, spg_pro
        implicit none
        integer(i4b), intent(in) :: tag
        logical, intent(out) :: flag
        integer(i4b) :: i,j,k,p,q,cho_res(100000), res_idx, iatom, patom, pt, ptyp, nf2
        logical :: f1, f2, f3
        real(dp) :: a,b,c,x(192), y(192), z(192)
        real(dp) :: pos_car(3), pos_dir(3), pos_tmp(3), pos_tmp_car(3), lattice(3,3), distance
        integer(i4b) :: type1, type2, natom, spg_old
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
        
        if(Q2D) then
            dnum = 0
            do i = -1, 1
                do j = -1, 1
                    dnum = dnum + 1
                    dx(dnum) = i
                    dy(dnum) = j
                    dz(dnum) = 0
                end do
            end do
        end if
           
        spg_index = pstruct(tag) % spg_idx
        spg_old = pstruct(tag) % spg_idx
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
            if(spg_index .ne. spg_old) then
                write(*, *) "nspg = ", spg_index, spg_old
                spg_index = spg_old
            end if
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
                                    do q = 1, dnum
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
                                    do q = 1, dnum
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
                            if(f2 .or. nf2 > 100) exit
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
        if(cnt > 10000) return
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
        ! generate structures
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
                pstruct(i) % pos(3, j) = pstruct(i) % pos(3, j) * ratio + 0.5
            end do
            pstruct(i) % lat(3, 3) = c2
        end if
    end subroutine init_struct_sym
    
end module init_struct_spg





module init_struct_all
    contains
    subroutine init_struct_tot(i)
        use kinds
        use parameters, only : cluster, cluster_substrate
        use init_struct_spg, only : init_struct_sym
        use init_struct, only : init_struct_cluster
        implicit none
        integer(i4b), intent(in) :: i
        if(cluster) then
            call init_struct_cluster(i)
        else
            call init_struct_sym(i)
        end if
    end subroutine init_struct_tot
end module init_struct_all
