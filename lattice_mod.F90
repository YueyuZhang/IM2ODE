!     
! File:   lattice_mod.F90
! Author: zyy
!
! Created on 2013年4月20日, 下午11:57
!

module lattice_mod
    contains
    subroutine init_lat_nosym(i)
        use kinds
        use constants, only : pi
        use parameters, only : population, pstruct, volumn
        use parameters, only : fix_lat, fix_a, fix_b, fix_c, fix_alpha, fix_beta, fix_gama
        implicit none
        integer(i4b), intent(in) :: i
        real(dp) :: a,b,c, alph, beta, gama, ratio, vtmp
        real(dp) :: lat_tmp(3, 3), vec(6)
        integer(i4b) :: j
        logical :: f1
        !do i = 1, population
        f1 = .false.
        do
            call random_number(a)
            call random_number(b)
            call random_number(c)
            !call random_number(alph)
            alph = 1.0
            alph = alph * pi / 2.0
            !call random_number(beta)
            beta = 1.0
            beta = beta * pi / 2.0
            !call random_number(gama)
            gama = 1.0
            gama = gama * pi / 2.0
            f1 = .true.
            if(a / b < 0.2 .or. a / b > 5.0) f1 = .false.
            if(a / c < 0.2 .or. a / c > 5.0) f1 = .false.
            if(b / c < 0.2 .or. b / c > 5.0) f1 = .false.
            if(alph < 10.0 / 180.0 * pi) f1 = .false.
            if(beta < 10.0 / 180.0 * pi) f1 = .false.
            if(gama < 10.0 / 180.0 * pi) f1 = .false.
            if(f1) exit
        end do
        if(fix_lat) then
            a = fix_a
            b = fix_b
            c = fix_c
            alph = fix_alpha / 90.0 * pi / 2.0
            beta = fix_beta / 90.0 * pi / 2.0
            gama = fix_gama / 90.0 * pi / 2.0
        end if
        vec(1) = a
        vec(2) = b
        vec(3) = c
        vec(4) = alph
        vec(5) = beta
        vec(6) = gama
        call vec2mat(vec, lat_tmp)
        call get_volumn(lat_tmp, vtmp)
        ratio = (volumn / vtmp) ** 0.33333333
        vec(1) = a * ratio
        vec(2) = b * ratio
        vec(3) = c * ratio
        call vec2mat(vec, lat_tmp)
        pstruct(i)%lat = lat_tmp
        !end do
    end subroutine init_lat_nosym
    
    subroutine vec2mat(vecx, matx)
        use kinds
        implicit none
        real(dp), intent(in) :: vecx(6)
        real(dp), intent(out) :: matx(3,3)
        matx = 0.0
        matx(1, 1) = vecx(1)
        matx(2, 1) = vecx(2) * cos(vecx(4))
        matx(2, 2) = vecx(2) * sin(vecx(4))
        matx(3, 1) = vecx(3) * cos(vecx(5))
        matx(3, 2) = vecx(3) * cos(vecx(6)) * sin(vecx(4)) - ((vecx(3) * cos(vecx(5))&
        &-vecx(3) * cos(vecx(6)) * cos(vecx(4))) / tan(vecx(4)))
        matx(3, 3) = sqrt(vecx(3) * vecx(3) - matx(3,1) * matx(3,1) - matx(3,2) * matx(3,2))
    end subroutine vec2mat
    
    subroutine get_volumn(matx, vol)
        use kinds
        implicit none
        real(dp), intent(in) :: matx(3,3)
        real(dp), intent(out) :: vol
        vol = 0.0
        vol = (matx(1,2)*matx(2,3)-matx(1,3)*matx(2,2))*matx(3,1)
        vol = vol + (matx(2,1)*matx(1,3)-matx(2,3)*matx(1,1))*matx(3,2)
        vol = vol + (matx(1,1)*matx(2,2)-matx(1,2)*matx(2,1))*matx(3,3)
    end subroutine get_volumn
    
    subroutine init_lat_sym(i)
        use kinds
        use constants, only : pi
        use parameters, only : population, pstruct, volumn
        use parameters, only : spg_index, spacegroup_log
        use parameters, only : spg_front, spg_rear
        use parameters, only : fix_lat, fix_a, fix_b, fix_c, fix_alpha, fix_beta, fix_gama
        implicit none
        integer(i4b), intent(in) :: i
        real(dp) :: a,b,c, alph, beta, gama, ratio, vtmp
        real(dp) :: lat_tmp(3, 3), vec(6)
        integer(i4b) :: j, n_spg
        logical :: f1
        !do i = 1, population
        f1 = .false.
        n_spg = spg_rear - spg_front + 1
        do
            call random_number(a)
            spg_index = floor(n_spg * a + spg_front)
            if(spacegroup_log(spg_index) == 0) then
                spacegroup_log(spg_index) = 1
                f1 = .true.
            end if
            if(n_spg < 100) then
                f1 = .true.
            end if
            if(f1) exit
            if(sum(spacegroup_log) > n_spg * 0.75)then
                do j = spg_front, floor(n_spg * a + spg_front)
                    if(spacegroup_log(j) == 0) then
                        spg_index = j
                        spacegroup_log(j) = 1
                        f1 = .true.
                    end if
                end do
            end if
            if(f1) exit
            if(sum(spacegroup_log) > floor(n_spg * 0.9)) then
                spacegroup_log = 0
            end if
        end do
        !if(spg_index > 230 .or. spg_index < 1) spg_index = 1
        if(spg_index > spg_rear .or. spg_index < spg_front) spg_index = spg_front
        pstruct(i) % spg_idx = spg_index
        !pstruct(i) % spg_idx = 3
        !write(*, *) "end init spg_INDEX"
        do
            if(fix_lat) exit
            call random_number(a)
            call random_number(b)
            call random_number(c)
            call random_number(alph)            
            call random_number(beta)            
            call random_number(gama)
            alph = alph * 1.8
            beta = beta * 1.8
            gama = gama * 1.8
            if(spg_index >= 3 .and. spg_index <= 15) then
                alph = 1.0
                gama = 1.0
            else if(spg_index >= 16 .and. spg_index <= 74) then
                alph = 1.0
                beta = 1.0
                gama = 1.0
            else if(spg_index >= 75 .and. spg_index <= 142) then
                b = a
                alph = 1.0
                beta = 1.0
                gama = 1.0
            else if(spg_index >= 143 .and. spg_index <= 167) then
                b = a
                c = a
                beta = alph
                gama = alph
            else if(spg_index >= 168 .and. spg_index <= 194) then
                b = a
                alph = 1.0
                beta = 1.0
                gama = 1.3333333
            else
                b = a
                c = a
                alph = 1.0
                beta = 1.0
                gama = 1.0
            end if
            alph = alph * pi / 2.0
            beta = beta * pi / 2.0
            gama = gama * pi / 2.0
            f1 = .true.
            if(a / b < 0.2 .or. a / b > 5.0) f1 = .false.
            if(a / c < 0.2 .or. a / c > 5.0) f1 = .false.
            if(b / c < 0.2 .or. b / c > 5.0) f1 = .false.
            if(alph < 10.0 / 180.0 * pi) f1 = .false.
            if(beta < 10.0 / 180.0 * pi) f1 = .false.
            if(gama < 10.0 / 180.0 * pi) f1 = .false.
            if(f1) exit
        end do
        if(fix_lat) then
            a = fix_a
            b = fix_b
            c = fix_c
            alph = fix_alpha / 90.0 * pi / 2.0
            beta = fix_beta / 90.0 * pi / 2.0
            gama = fix_gama / 90.0 * pi / 2.0
        end if
        !write(*, *) "end init LAT"
        vec(1) = a
        vec(2) = b
        vec(3) = c
        vec(4) = alph
        vec(5) = beta
        vec(6) = gama
        call vec2mat(vec, lat_tmp)
        call get_volumn(lat_tmp, vtmp)
        ratio = (volumn / vtmp) ** 0.33333333
        if(fix_lat) ratio = 1.0
        vec(1) = a * ratio
        vec(2) = b * ratio
        vec(3) = c * ratio
        call vec2mat(vec, lat_tmp)
        pstruct(i)%lat = lat_tmp
    end subroutine init_lat_sym
end module lattice_mod
    
            
