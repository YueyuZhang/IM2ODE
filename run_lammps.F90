!     
! File:   run_lammps.F90
! Author: zyy
!
! Created on 2013年4月22日, 下午2:34
!

module run_lammps
    contains
    subroutine lammps()
        use kinds
        use parameters, only : population
        implicit none
        integer(i4b) :: i,j
        do i = 1, population
            call print_ptl(i)
            call system('./lmp_serial < in.put > log.out')
            call read_lmp(i)
        end do
    end subroutine lammps
    
    subroutine print_ptl(tag)
        use kinds
        use parameters, only : pstruct
        implicit none
        integer(i4b), intent(in) :: tag
        integer(i4b) :: i, j, n, atom_typ
        real(dp) :: lattice(3,3)
        n = pstruct(tag) % natom
        atom_typ = pstruct(tag) % ntyp
        lattice = pstruct(tag) % lat
        open(3111, file = "data.input")
        write(3111, "(1X, A)") "lammps data_input for Si bulk"
        write(3111, *)
        write(3111, "(1X, I5, A)") n, " atoms"
        write(3111, "(1X, I5, A)") atom_typ, " atom types"
        write(3111, *)
        write(3111, "(1X, 2F12.6, A)") 0.0, lattice(1,1), " xlo xhi"
        write(3111, "(1X, 2F12.6, A)") 0.0, lattice(2,2), " ylo yhi"
        write(3111, "(1X, 2F12.6, A)") 0.0, lattice(3,3), " zlo zhi"
        write(3111, *)
        write(3111, *) "Atoms"
        write(3111, *)
        do i = 1, n
            write(3111, "(1X, 2I5, 3F12.6)") i, pstruct(tag) % ptype(i), pstruct(tag) % pos(:, i)
        end do
        close(3111)
    end subroutine print_ptl
    
    subroutine read_lmp(tag)
        use kinds
        use parameters, only : pstruct
        implicit none
        integer(i4b), intent(in) :: tag
        integer(i4b) :: i,j,k
        integer(i4b) :: nline1, nline2
        integer(i4b) :: status = 0
        character(len = 200) :: str1, str2
        real(dp) :: a,b,c, energy1, energy2, ss(3)
        nline1 = 0
        nline2 = 0
        do
            call system('sleep 0.03')
            status = 0
            nline1 = 0
            open(3211, file = "log.lammps", status = "old")
            do 
                nline1 = nline1 + 1
                read(3211, fmt = "(A200)", iostat=status) str1
                j = index(str1, 'Loop time of ')
                if(j /= 0) then
                    read(str2, *) j, energy1, energy2, a, b, c
                end if
                str2 = str1
                if(status /= 0) exit
            end do
            if(nline1 .eq. nline2) exit
            nline2 = nline1
            close(3211)
        end do
        status = 0
        !write(*, *) tag, energy1
        pstruct(tag) % energy = energy1
        pstruct(tag) % lat(1,1) = a
        pstruct(tag) % lat(2,2) = b
        pstruct(tag) % lat(3,3) = c
        open(3212, file = "dump.atom", status = "old")
        nline1 = pstruct(tag) % natom
        status = 0
        do
            read(3212, fmt = "(A200)", iostat = status) str1
            if(status /= 0) exit
            read(3212, *)
            do i = 1, nline1
                read(3212, *) j, ss(1), ss(2), ss(3)
                pstruct(tag) % pos(:, i) = ss(:)
                !write(*, *) ss(:)
            end do
        end do
        close(3212)
    end subroutine read_lmp
    
    subroutine print_lmp_res()
        use kinds
        use parameters, only : pool, population
        implicit none
        integer(i4b) :: i,j,k
        character(len=20) :: fname
        call system('mkdir results')
        do i = 1, population
            write(*, *), i, "  energy = ",pool(i) % energy 
            write(fname, "(A7, I3.3)") "results/POSCAR_", i
            open(3113, file = fname)
            do j = 1, 3
                write(3113, "(3F12.6)") pool(i) % lat(j,:)
            end do
            do j = 1, pool(i) % ntyp
                write(3113, *) pool(i) % nelement(j)
            end do
            do j = 1, pool(i) % natom
                do k = 1, 3
                    pool(i) % pos(k, j) = pool(i) % pos(k, j) / pool(i) % lat(k, k)
                end do
                write(3113, "(3F12.6)") pool(i) % pos(:, j)
            end do
            close(3113)
        end do
    end subroutine print_lmp_res
end module run_lammps
    

