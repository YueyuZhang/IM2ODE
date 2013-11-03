!     
! File:   init_mod.F90
! Author: zyy
!
! Created on 2013年4月20日, 下午4:19
!

module init_mod
    use kinds
    use parameters, only : spacegroup_log
    use spacegroup_data, only : init_spg
    use parameters, only : cluster_substrate
    use de_tools, only : read_init_struct
    implicit none
    integer(i4b) :: flag
    contains
    subroutine init_run()
        use kinds
        use parameters, only : pool, population
        implicit none
        integer(i4b) :: i
        real(dp) :: inf = 1e8
        call init_rand()
        call system("rm -rf results")
        call system("mkdir results")
        open(1224,file="results/run_de.out")
        write(1224, *) "start de searching ......"
        call read_input()
        call write_input()
        spacegroup_log = 0
        
        do i = 1, population
            pool(i) % energy = inf
            pool(i) % hardness = inf
            pool(i) % Eg_id = inf
            pool(i) % Eg_d = inf
        end do
        
        call init_spg
        
        if(cluster_substrate) then
            call read_init_struct(flag)
            if(flag == 0) then
                write(1224, *) "Error: wrong format if input structure"
                stop
            end if
        end if
        
    end subroutine init_run
    
    subroutine init_rand()
        use kinds
        integer(i4b) :: i, n, clock
        real(dp) :: x, y, time
        integer(i4b), dimension(:), allocatable :: seed
        call cpu_time(time)
        call random_seed(size = n)
        allocate(seed(n))
        call system_clock(count=clock)
        seed = clock + 37 * (/(i - 1, i = 1, n)/)
        call random_seed(put = seed)
        deallocate(seed)
    end subroutine init_rand
    
    subroutine read_input()
        use parameters, only : sys_name, num_species, num_ele, name_element, atom_dis, volumn
        use parameters, only : population, max_step, de_ratio, symmetry, spg_front, spg_rear
        use parameters, only : mode, hardness, rcut, ionicity
        use parameters, only : ESflag, ES_mod, ES_Eg, Es_opt
        use parameters, only : fix_lat, fix_a, fix_b, fix_c, fix_alpha, fix_beta, fix_gama
        use parameters, only : Q2D, vacuum_layer, Area
        use parameters, only : find_defect, defect_type, center_x, center_y, center_z, defect_radius, length_x, length_y, length_z
        use parameters, only : PRESSURE, PSTRESS
        use parameters, only : cluster, init_radius, cluster_substrate, SubstrateElements
        use parameters, only : cluster_ctr_x, cluster_ctr_y, cluster_ctr_z
        use kinds
        implicit none
        
        integer(i4b) :: n, i, j, k, line, f1, lth
        character(len=40) :: nametag(200), number(200)
        character(len=200) :: strin(200), strtmp
        logical :: flag
        
        inquire(file="de.in", exist = flag)
        if(.not. flag) then
            write(1224, *) "no de.in"
            stop
        end if
        open(unit = 5111, file = "de.in", status = "old")
        line = 0
        f1 = 0
        do while(.true.)
            line = line + 1
            read(5111, fmt="(A200)", iostat=f1) strin(line)
            if(f1 /= 0) exit
        end do
        close(5111)
        
        n = 0
        do i = 1, line
            if(len(trim(strin(i))) == 0 .or. len(trim(strin(i))) == 1) then
                continue
            else
                read(strin(i), *) strtmp
                if(strtmp /= '#') then
                    lth = index(strin(i), '=')
                    if(lth /= 0) then
                        n = n + 1
                        read(strin(i)(:lth-1), "(A40)")nametag(n)
                        read(strin(i)(lth+1:), "(A40)")number(n)
                    end if
                end if
            end if
        end do
        
        call find(nametag, 'SystemName', i)
        if(i == 0) then
            write(1224, *) "Input SystemName"
            sys_name = "DE-search"
        else
            read(number(i), *) sys_name
        end if
        
        call find(nametag, 'NumberOfSpecies', i)
        if(i == 0) then
            write(1224, *) "Input NumberOfSpecies"
            num_species = 1
        else
            read(number(i), *) num_species
        end if
        
        call find(nametag, "NumberOfElements", i)
        if(i == 0) then
            write(1224, *) "Input NumberOfElements"
            num_ele(1) = 2
        else
            read(number(i), *) (num_ele(j), j = 1, num_species)
        end if
        
        call find(nametag, "NameOfElements", i)
        if(i == 0) then
            write(1224, *) "Input NameOfElements"
            name_element(1) = "X"
        else
            read(number(i), *) (name_element(j), j = 1, num_species)
        end if
        
        call find(nametag, "Volumn", i)
        if(i == 0) then
            write(1224, *) "Input Volumn"
            volumn = 10
        else
            read(number(i), *) volumn
        end if
        
        call find(nametag, "DistanceOfAtom", i)
        if(i == 0) then
            write(1224, *) "Input DistanceOfAtom"
        else
            do j = 1, num_species
                read(number(i + j), *) (atom_dis(j, k), k = 1, num_species)
            end do
        end if
        
        call find(nametag, "Population", i)
        if(i == 0) then
            write(1224, *) "Input Population"
            population = 10
        else
            read(number(i), *) population
        end if
        
        call find(nametag, "MaxStep", i)
        if(i == 0) then
            write(1224, *) "Input MaxStep"
            max_step = 10
        else
            read(number(i), *) max_step
        end if
        
        call find(nametag, "De_ratio", i)
        if(i == 0) then
            write(1224, *) "Input De_ratio"
            de_ratio = 0.8
        else
            read(number(i), *) de_ratio
        end if
        
        call find(nametag, "Symmetry", i)
        if(i == 0) then
            symmetry = .false.
        else
            read(number(i), *) symmetry
        end if
        
        call find(nametag, "spg_front", i)
        if(i == 0) then
            spg_front = 1
        else
            read(number(i), *) spg_front
        end if
        
        call find(nametag, "spg_rear", i)
        if(i == 0) then
            spg_rear = 230
        else
            read(number(i), *) spg_rear
        end if
        
        call find(nametag, "Multi-Objective", i)
        if(i == 0) then
            mode = .false.
        else
            read(number(i), *) mode
        end if
        
        call find(nametag, "hardness", i)
        if(i == 0) then
            hardness = .false.
        else
            read(number(i), *) hardness
            
            call find(nametag, "rcut", i)
            if(i == 0) then
                write(1224,*) "Input rcut"
                rcut = 1.5
            else
                read(number(i), *) rcut
            end if
            
            call find(nametag, "ionicity", i)
            if(i == 0) then
                write(1224,*) "Input ionicity"
                ionicity = 0.0
            else
                read(number(i), *) ionicity
            end if
        end if
        
        call find(nametag, "ESflag", i)
        if(i == 0) then
            ESflag = .false.
        else
            read(number(i), *) ESflag
            
            call find(nametag, "ES_mod", i)
            if(i == 0) then
                write(1224, *) "Input ES_mod"
                ES_mod = 6
            else
                read(number(i), *) ES_mod
            end if
            
            call find(nametag, "ES_Eg", i)
            if(i == 0) then
                write(1224, *) "Input ES_Eg"
                ES_Eg = 1.0
            else
                read(number(i), *) ES_Eg
            end if
            
            call find(nametag, "ES_opt", i)
            if(i == 0) then
                write(1224, *) "Inout ES_pot"
                ES_opt = 1.2
            else
                read(number(i), *) ES_opt
            end if
        end if
        
        call find(nametag, "Q2D", i)
        if(i == 0) then
            Q2D = .false.
        else
            read(number(i), *) Q2D
            
            call find(nametag, "vacuum_layer", i)
            if(i == 0) then
                vacuum_layer = 10.0
            else
                read(number(i), *) vacuum_layer
            end if
            
            call find(nametag, "Area", i)
            if(i == 0) then
                Area = 100.0
            else
                read(number(i), *) Area
            end if
        end if
        
        call find(nametag, "cluster", i)
        if(i == 0) then
            cluster = .false.
        else
            read(number(i), *) cluster
            
            call find(nametag, "init_radius", i)
            if(i == 0) then
                write(1224, *) "Input init_radius"
                init_radius = 2.5
            else
                read(number(i), *) init_radius
            end if
            
            call find(nametag, "cluster_ctr_x", i)
            if(i == 0) then
                write(1224, *) "Input center of cluster"
                cluster_ctr_x = 0.5
            else
                read(number(i), *) cluster_ctr_x
            end if
            
            call find(nametag, "cluster_ctr_y", i)
            if(i == 0) then
                write(1224, *) "Input center of cluster"
                cluster_ctr_y = 0.5
            else
                read(number(i), *) cluster_ctr_y
            end if
            
            call find(nametag, "cluster_ctr_z", i)
            if(i == 0) then
                write(1224, *) "Input center of cluster"
                cluster_ctr_z = 0.5
            else
                read(number(i), *) cluster_ctr_z
            end if
            
            call find(nametag, "cluster_substrate", i)
            if(i == 0) then
                cluster_substrate = .false.
            else
                read(number(i), *) cluster_substrate
                
                call find(nametag, "SubstrateElements", i)
                if(i == 0) then
                    write(1224, *) "Input SubstrateElements"
                else
                    read(number(i), *) (SubstrateElements(j), j = 1, num_species)
                end if
            end if
        end if
        
        call find(nametag, "fix_lat", i)
        if(i == 0) then
            fix_lat = .false.
        else
            read(number(i), *) fix_lat
            
            call find(nametag, "fix_a", i)
            if(i == 0) then
                write(1224, *) "Input fix_a"
                fix_a = 1.0
            else
                read(number(i), *) fix_a
            end if
            
            call find(nametag, "fix_b", i)
            if(i == 0) then
                write(1224, *) "Input fix_b"
                fix_b = 1.0
            else
                read(number(i), *) fix_b
            end if
            
            call find(nametag, "fix_c", i)
            if(i == 0) then
                write(1224, *) "Input fix_c"
                fix_c = 10.0
            else
                read(number(i), *) fix_c
            end if
            
            call find(nametag, "fix_alpha", i)
            if(i == 0) then
                write(1224, *) "Input fix_alpha"
                fix_alpha = 90
            else
                read(number(i), *) fix_alpha
            end if
            
            call find(nametag, "fix_beta", i)
            if(i == 0) then
                write(1224, *) "Input fix_beta"
                fix_beta = 90
            else
                read(number(i), *) fix_beta
            end if
            
            call find(nametag, "fix_gama", i)
            if(i == 0) then
                write(1224, *) "Input fix_gama"
                fix_gama = 90
            else
                read(number(i), *) fix_gama
            end if
        end if
        
        call find(nametag, "PSTRESS", i)
        if(i == 0) then
            PRESSURE = .false.
        else
            PRESSURE = .true.
            read(number(i), *) PSTRESS
        end if
        
        call find(nametag, "find_defect", i)
        if(i == 0) then
            find_defect = .false.
        else
            read(number(i), *) find_defect
            
            call find(nametag, "defect_type", i)
            if(i == 0) then
                write(1224, *) "Input defect_type"
                defect_type = 1
            else
                read(number(i), *) defect_type
            end if
            
            call find(nametag, "center_x", i)
            if(i == 0) then
                write(1224, *) "Input center_x"
                center_x = 1.0
            else
                read(number(i), *) center_x
            end if
            
            call find(nametag, "center_y", i)
            if(i == 0) then
                write(1224, *) "Input center_y"
                center_y = 1.0
            else
                read(number(i), *) center_y
            end if
            
            call find(nametag, "center_z", i)
            if(i == 0) then
                write(1224, *) "Input center_z"
                center_z = 1.0
            else
                read(number(i), *) center_z
            end if
            
            if(defect_type == 1) then
                call find(nametag, "length_x", i)
                if(i == 0) then
                    write(1224, *) "Input length_x"
                    length_x = 1.0
                else
                    read(number(i), *) length_x
                end if
                
                call find(nametag, "length_y", i)
                if(i == 0) then
                    write(1224, *) "Input length_y"
                    length_y = 1.0
                else
                    read(number(i), *) length_y
                end if
                
                call find(nametag, "length_z", i)
                if(i == 0) then
                    write(1224, *) "Input length_z"
                    length_z = 1.0
                else
                    read(number(i), *) length_z
                end if
            else
                call find(nametag, "defect_radius", i)
                if(i == 0) then
                    write(1224, *) "Input defect_radius"
                    defect_radius = 1.0
                else
                    read(number(i), *) defect_radius
                end if
            end if
            
        end if
        
    end subroutine read_input
    
    subroutine find(a, b, i)
        use kinds
        character(len=40) :: a(200)
        character(len=*) :: b
        integer(i4b) :: i,j
        i = 0
        do j = 1, 200
            if(trim(a(j)) == b) then
                i = j
                exit
            end if
        end do
    end subroutine find
    
    subroutine write_input()
        use kinds
        use parameters
        implicit none
        integer(i4b) :: i, j
        write(1224, *) "---write input parameters---"
        write(1224, *) "system name: ", sys_name
        write(1224, *) "NumberOfSpecies: ", num_species
        write(1224, *) "NumberOfElements: ", (num_ele(i), i = 1, num_species)
        write(1224, *) "NameOfElements: ", (' '//name_element(i), i = 1, num_species)
        write(1224, *) "Volumn: ", volumn
        write(1224, *) "DistanceOfAtom: "
        do i = 1, num_species
            write(1224, *) (atom_dis(i, j), j = 1, num_species)
        end do
        write(1224, *) "Population: ", population
        write(1224, *) "MaxStep: ", max_step
        write(1224, *) "De_ratio: ", de_ratio
        write(1224, *) "Symmetry: ", symmetry
        write(1224, *) "Q2D: ", Q2D
        if(Q2D) then
            write(1224, *) "vacuum_layer: ", vacuum_layer
            write(1224, *) "Area: ", Area
        end if
        write(1224, *) "Multi-objective: ", mode
        write(1224, *) "cluster: ", cluster
        if(cluster) then
            write(1224, *) "init_radius: ", init_radius
            write(1224, *) "cluster_ctr_x: ", cluster_ctr_x
            write(1224, *) "cluster_ctr_y: ", cluster_ctr_y
            write(1224, *) "cluster_ctr_z: ", cluster_ctr_z
            write(1224, *) "cluster_substrate: ", cluster_substrate
            if(cluster_substrate) then
                write(1224, *) "SubstrateElements: ", (SubstrateElements(i), i = 1, num_species)
            end if
        end if
        write(1224, *) "fix_lat: ", fix_lat
        if(fix_lat) then
            write(1224, *) "fix_a: ", fix_a
            write(1224, *) "fix_b: ", fix_b
            write(1224, *) "fix_c: ", fix_c
            write(1224, *) "fix_alpha: ", fix_alpha
            write(1224, *) "fix_beta: ", fix_beta
            write(1224, *) "fix_gama: ", fix_gama
        end if
        if(PRESSURE) then
            write(1224, *) "PSTRESS: ", PSTRESS
        end if
        write(1224, *) "find_defect: ", find_defect
        if(find_defect) then
            write(1224, *) "defect_type: ", defect_type
            write(1224, *) "center: ", center_x, center_y, center_z
        end if
        write(1224, *) "---end write input---"
    end subroutine write_input
        
end module init_mod

