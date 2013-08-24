!     
! File:   tools.F90
! Author: zyy
!
! Created on 2013年4月28日, 上午8:35
!

module de_tools
    implicit none
    contains
    
    function itoa(i)
        use kinds
        implicit none
        integer(i4b), intent(in) :: i
        character(len=10) :: itoa
        if(i < 10) then
            write(itoa, "(I1)"), i
        else if(i < 100) then
            write(itoa, "(I2)"), i
        else if(i < 1000) then
            write(itoa, "(I3)"), i
        else if(i < 10000) then
            write(itoa, "(I4)"), i
        else 
            write(itoa, "(I5)"), i
        end if
        return
    end function itoa
    
    subroutine write_vasp(tag)
        use kinds
        use parameters, only : pstruct, sys_name, name_element
        implicit none
        integer(i4b), intent(in) :: tag
        character(len=10) :: fname
        integer(i4b) :: i,j,k
        fname = "POSCAR"//trim(itoa(tag))
        open(4111, file=fname)
        write(4111, *) sys_name
        write(4111,"(A6)")"1.0000"
        do i = 1, 3
            write(4111,"(3F10.6)")(pstruct(tag)%lat(i,j),j=1,3)
        end do
        write(4111,*) (name_element(i)//" ",i=1,pstruct(tag)%ntyp)
        write(4111, *) (pstruct(tag)%nelement(i), i=1,pstruct(tag)%ntyp)
        write(4111,"(A6)")"Direct"
        do i = 1, pstruct(tag)%natom
            write(4111, "(3F10.6)")(pstruct(tag)%pos(j, i),j=1,3)
        end do
        close(4111)
    end subroutine write_vasp
    
    
    subroutine write_input_all(step)
        use kinds
        use parameters, only : pstruct, population, sys_name, name_element
        implicit none
        integer(i4b), intent(in) :: step
        integer(i4b) :: tag,i,j
        character(len = 20) :: fname
        fname = "results/de_ini_"//trim(itoa(step))
        open(4313, file=fname)
        do tag = 1, population
            !write(4313, *) pstruct(tag) % energy / pstruct(tag) % natom
            write(4313, *) sys_name
            write(4313,"(A6)")"1.0000"
            do i = 1, 3
                write(4313,"(3F10.6)")(pstruct(tag)%lat(i,j),j=1,3)
            end do
            write(4313,*) (name_element(i)//" ",i=1,pstruct(tag)%ntyp)
            write(4313, *) (pstruct(tag)%nelement(i), i=1,pstruct(tag)%ntyp)
            write(4313,"(A6)")"Direct"
            do i = 1, pstruct(tag)%natom
                write(4313, "(3F10.6)")(pstruct(tag)%pos(j, i),j=1,3)
            end do
        end do
        close(4313)
    end subroutine write_input_all
    
    subroutine write_vasp_all(step)
        use kinds
        use parameters, only : pstruct, population, sys_name, name_element
        use parameters, only : mode
        implicit none
        integer(i4b), intent(in) :: step
        integer(i4b) :: tag,i,j
        character(len = 20) :: fname
        fname = "results/de_opt_"//trim(itoa(step))
        open(4313, file=fname)
        do tag = 1, population
            write(4313, *) pstruct(tag) % energy / pstruct(tag) % natom
            if(mode) then
                write(4313, *) pstruct(tag) % hardness
            end if
            write(4313, *) sys_name
            write(4313,"(A6)")"1.0000"
            do i = 1, 3
                write(4313,"(3F10.6)")(pstruct(tag)%lat(i,j),j=1,3)
            end do
            write(4313,*) (name_element(i)//" ",i=1,pstruct(tag)%ntyp)
            write(4313, *) (pstruct(tag)%nelement(i), i=1,pstruct(tag)%ntyp)
            write(4313,"(A6)")"Direct"
            do i = 1, pstruct(tag)%natom
                write(4313, "(3F10.6)")(pstruct(tag)%pos(j, i),j=1,3)
            end do
        end do
        close(4313)
    end subroutine write_vasp_all
    
    subroutine read_vasp(tag)
        use kinds
        use parameters, only : pstruct, ESflag
        use constants, only: max_atom, max_type
        use Elec_Str_mod, only : ES_gap
        implicit none
        integer(i4b), intent(in) :: tag
        integer(i4b) :: i,j,k,natom,lth, f1
        real(dp) :: inf = 1e5
        real(dp) :: tmp, pos(3, max_atom * max_type), lat_matrix(3,3)
        real(dp) :: energy
        character(len=200) :: strin
        character(len=40) :: nametag, number
        logical :: alive
        call system("cd vasp_"//trim(itoa(tag))//"; cp CONTCAR OUTCAR EIGENVAL OSZICAR ..;"//&
        & "rm -rf ../Matrix_elements; (if [ -s Matrix_elements ]; then cp Matrix_elements ..; fi) ")
        inquire(file = "CONTCAR", exist = alive)
        if(alive) then
            open(4311, file="CONTCAR", status = "old")
        else
            write(1224, *) "no CONTCAR of pop ", tag
            pstruct(tag) % energy = inf
            return
        end if
        read(4311, *, iostat = f1)
        if(f1 /= 0) then
            pstruct(tag) % energy = inf
            close(4311)
            return
        end if
        read(4311, *, iostat = f1) tmp
        if(f1 /= 0) then
            pstruct(tag) % energy = inf
            close(4311)
            return
        end if
        do i = 1, 3
            read(4311, *, iostat = f1)(lat_matrix(i,j),j = 1, 3)
            if(f1 /= 0) then
                pstruct(tag) % energy = inf
                close(4311)
                return
            end if
        end do
        lat_matrix(: ,:) = tmp * lat_matrix(:, :)
        read(4311, *, iostat = f1)
        if(f1 /= 0) then
            pstruct(tag) % energy = inf
            close(4311)
            return
        end if
        read(4311, *, iostat = f1)
        if(f1 /= 0) then
            pstruct(tag) % energy = inf
            close(4311)
            return
        end if
        read(4311, *, iostat = f1)
        if(f1 /= 0) then
            pstruct(tag) % energy = inf
            close(4311)
            return
        end if
        natom = pstruct(tag) % natom
        do i = 1, natom
            read(4311, *, iostat = f1) (pos(j,i), j = 1, 3)
            if(f1 /= 0) then
            pstruct(tag) % energy = inf
                close(4311)
                return
            end if
            pos(1,i) = pos(1,i) - floor(pos(1,i))
            pos(2,i) = pos(2,i) - floor(pos(2,i))
            pos(3,i) = pos(3,i) - floor(pos(3,i))
        end do
        close(4311)
        pstruct(tag) % lat = lat_matrix
        pstruct(tag) % pos = pos
        
        energy = inf
        inquire(file = "OUTCAR", exist = alive)
        if(alive) then
            open(4312, file="OUTCAR", status = "old")
            energy = inf
        else
            write(1224, *) "no OUTCAR of pop ", tag
            pstruct(tag) % energy = inf
            return
        end if
        do while(.true.)
            read(4312, fmt="(A200)", iostat=f1) strin
            if(f1 /= 0) exit
            lth = index(strin, '=')
            if(lth /= 0 .and. lth /= 1) then
                read(strin(:lth-1), "(A40)") nametag
                read(strin(lth+1:), "(A40)") number
                if(trim(nametag) == trim("  free  energy   TOTEN  ")) then
                    read(number, *) energy
                end if
            end if
        end do
        pstruct(tag) % energy = energy
        write(*, *) tag, "energy", energy
        close(4312)
        
        if(ESflag) then
            call ES_gap(tag)
        end if
    end subroutine read_vasp
    
    subroutine bulkmodu(tag)
        use kinds, only : dp, i4b
        use parameters, only : pstruct, rcut, ionicity
        use constants
        implicit none
        integer(i4b), intent(in) :: tag
        integer(i4b) :: i, j, k, natomele(max_type), natom, n, bondnum
        integer(i4b) :: dx(7), dy(7), dz(7)
        real(dp) :: lat_matrix(3,3), pos(3,max_atom * max_type), carpos(3,max_atom * max_type)
        real(dp) :: tmp, totbondlen, bondlen, Nc
        real(dp) :: Hv
        dx = (/1,0,0,1,1,0,1/)
        dy = (/0,1,0,1,0,1,1/)
        dz = (/0,0,1,0,1,1,1/)
        
        natom = pstruct(tag) % natom
        lat_matrix = pstruct(tag) % lat
        pos = pstruct(tag) % pos
        natomele = pstruct(tag) % nelement
        
        n = natom * 8
        bondnum = 0
        totbondlen = 0.0
        do i = 1, 7
            do j = 1, natom
                k = natom * i + j
                pos(1, k) = pos(1, j) + dx(i)
                pos(2, k) = pos(2, j) + dy(i)
                pos(3, k) = pos(3, j) + dz(i)
            end do
        end do
        carpos = 0.0
        do i = 1, n
            do j = 1, 3
                do k = 1, 3
                    carpos(j,i) = carpos(j,i) + pos(k,i) * lat_matrix(k,j)
                end do
            end do
        end do
        do i = 1, natom
            do j = i + 1, n
                tmp = 0
                do k = 1, 3
                    tmp = tmp + (carpos(k,i) - carpos(k,j)) ** 2.0
                end do
                tmp = sqrt(tmp)
                if(tmp < rcut) then
                    bondnum = bondnum + 1
                    totbondlen = totbondlen + tmp
                end if
            end do
        end do
        if (bondnum .eq. 0) then
            pstruct(tag) % hardness = 1000
            return
        end if
        bondlen = totbondlen / real(bondnum) / 10.0
        Nc = real(bondnum) / real(natom) * 2.0
        Hv = Nc / 4.0 * (0.624 - 0.07 * ionicity) / (bondlen ** 3.5)
        pstruct(tag) % hardness = -1.0 * Hv
    end subroutine bulkmodu
    
    subroutine ES_fitness(tag)
        use kinds
        use parameters, only: pstruct, ES_mod, ES_Eg
        implicit none
        integer(i4b), intent(in) :: tag
        real(dp) :: bandgap, direct_gap, weight
        
        bandgap = pstruct(tag) % Eg_id
        direct_gap = pstruct(tag) % Eg_d
        if(ES_mod .EQ. 1) then
            ! 1: maximal bandgap
            pstruct(tag) % hardness = -bandgap
            if(bandgap .GT. 1.0d6) then
                pstruct(tag) % hardness = 1.0d6
            end if
        elseif(ES_mod .EQ. 2) then
            ! 2: minimal bandgap
            pstruct(tag) % hardness = bandgap
        elseif(ES_mod .EQ. 3) then
            ! 3: a target bandgap
            pstruct(tag) % hardness = abs(bandgap - ES_Eg)
        elseif(ES_mod .EQ. 4) then
            ! 4: only direct gap, largest
            weight = 8.0d0
            pstruct(tag) % hardness = -(bandgap - weight * (direct_gap - bandgap))
            if(bandgap .GT. 1.0d6) then
                pstruct(tag) % hardness  = 1.0d6
            end if
        elseif (ES_mod .EQ. 5) then
            ! 5: only direct gap, largest
            IF (bandgap .GT. 1.0d6)THEN
                pstruct(tag) % hardness = 1.0d6 !vasp results are not complete            
            ELSEIF (bandgap .GT. 9000)THEN
                pstruct(tag) % hardness = 9000 !number of electrons is odd          
            ELSEIF (bandgap .GT. 4000)THEN
                !no transition between gap 
                IF (direct_gap - bandgap .GT. 0.03)THEN
                    !indirect gap system
                    weight = 8.0d0
                    pstruct(tag) % hardness = -(bandgap - 5000 - weight * (direct_gap - bandgap))
                    pstruct(tag) % hardness = pstruct(tag) % hardness + 10.0d0
                ELSE
                    pstruct(tag) % hardness = -(bandgap - 5000) + 10.0d0
                END IF
            ELSE
                IF (direct_gap - bandgap .GT. 0.03)THEN
                !indirect gap system                                                                      
                weight = 8.0d0
                pstruct(tag) % hardness = -(bandgap - weight * (direct_gap - bandgap))
            ELSE
                pstruct(tag) % hardness = -bandgap
                END IF
            END IF
        elseif(ES_mod .EQ. 6) then
            ! 6: a target direct gap
            pstruct(tag) % hardness = abs(bandgap - ES_Eg) + abs(direct_gap - bandgap)
        elseif(ES_mod .EQ. 7) then
            if(bandgap .GT. 0.70d0 .AND. abs(bandgap) .LT. 2.0d0) then
                pstruct(tag) % hardness = 4 * abs(direct_gap - bandgap)
            else
                pstruct(tag) % hardness = abs(bandgap - 0.7) * 20
            end if
        end if
    end subroutine ES_fitness

end module de_tools

