!     
! File:   init_mod.F90
! Author: zyy
!
! Created on 2013年4月20日, 下午4:19
!

module init_mod
    use kinds
    use parameters, only : spacegroup_log, Pickup
    use spacegroup_data, only : init_spg
    use parameters, only : cluster_substrate, find_defect, model_surface, cluster
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
        open(1224,file="run_de.out")
        call read_input()
        if(.not. Pickup) then
            call system("rm -rf results")
            call system("mkdir results")
            write(1224, *) "start de searching ......"
        else
            write(1224, *) "======CONTINUE======="
        end if
        !call read_input()
        call write_input()
        spacegroup_log = 0
        
        do i = 1, population
            pool(i) % energy = inf
            pool(i) % hardness = inf
            pool(i) % Eg_id = inf
            pool(i) % Eg_d = inf
        end do
        
        call init_spg
        
        if(model_surface) then
            cluster = .true.
            cluster_substrate = .true.
        end if
        
        if(cluster_substrate) then
            call read_init_struct(flag)
            if(flag == 0) then
                write(1224, *) "Error: wrong format of input structure"
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
        use parameters, only : Pickup, Pickup_step, Mystruct
        use parameters, only : mode, hardness, rcut, ionicity
        use parameters, only : ESflag, ES_mod, ES_Eg, Es_opt
        use parameters, only : HSE, HSE_population, LDA_population, energy_cut, gap_cut
        use parameters, only : LDA_ES_Eg, LDA_Es_opt, HSE_ES_Eg, HSE_Es_opt
        use parameters, only : fix_lat, fix_a, fix_b, fix_c, fix_alpha, fix_beta, fix_gama
        use parameters, only : Q2D, vacuum_layer, Area, Layer_hight
        use parameters, only : find_defect, defect_type, center_x, center_y, center_z, defect_radius, length_x, length_y, length_z
        use parameters, only : PRESSURE, PSTRESS
        use parameters, only : cluster, model_ball, model_shell, model_plate
        use parameters, only : init_radius, cluster_substrate, SubstrateElements
        use parameters, only : model_surface, surface_height
        use parameters, only : cluster_ctr_x, cluster_ctr_y, cluster_ctr_z
        use parameters, only : shell_radius_in, shell_radius_out
        use parameters, only : shell_ctr_x, shell_ctr_y, shell_ctr_z
        use parameters, only : plate_radius, plate_height
        use parameters, only : plate_ctr_x, plate_ctr_y, plate_ctr_z
        use parameters, only : selective_dynamics
        use parameters, only : dimer_mode, dimer_dis
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
                cycle
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
        
        call find(nametag, "Pickup", i)
        if(i == 0) then
            write(1224, *) "Input Pickup"
            Pickup = .false.
        else
            read(number(i), *) Pickup
            if(Pickup) then
                call find(nametag, "Pickup_step", i)
                if(i == 0) then
                    write(1224, *) "Input Pickup_step"
                    Pickup_step = 1
                else
                    read(number(i), *) Pickup_step
                end if
            end if
        end if
        
        call find(nametag, "Mystruct", i)
        if(i == 0) then
            write(1224, *) "Input Mystruct"
            Mystruct = 0
        else
            read(number(i), *) Mystruct
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
                write(1224, *) "Input ES_pot"
                ES_opt = 1.2
            else
                read(number(i), *) ES_opt
            end if
            
            call find(nametag, "HSE", i)
            if(i == 0) then
                write(1224, *) "Input HSE"
                HSE = .false.
            else
                read(number(i), *) HSE
                
                call find(nametag, "HSE_population", i)
                if(i == 0) then
                    write(1224, *) "Input HSE_population"
                    HSE_population = 0
                else
                    read(number(i), *) HSE_population
                end if
                
                call find(nametag, "LDA_population", i)
                if(i == 0) then
                    write(1224, *) "Input LDA_population"
                    LDA_population = 6
                else
                    read(number(i), *) LDA_population
                end if
                
                call find(nametag, "energy_cut", i)
                if(i == 0) then
                    write(1224, *) "Input energy_cut"
                    energy_cut = -3
                else
                    read(number(i), *) energy_cut
                end if
                
                call find(nametag, "gap_cut", i)
                if(i == 0) then
                    write(1224, *) "Input gap_cut"
                    gap_cut = 5
                else
                    read(number(i), *) gap_cut
                end if
                
                call find(nametag, "LDA_ES_Eg", i)
                if(i == 0) then
                    write(1224, *) "Input LDA_ES_Eg"
                    LDA_ES_Eg = 0.5
                else
                    read(number(i), *) LDA_ES_Eg
                end if
                
                call find(nametag, "LDA_Es_opt", i)
                if(i == 0) then
                    write(1224, *) "Input LDA_Es_opt"
                    LDA_Es_opt = 2.0
                else
                    read(number(i), *) LDA_Es_opt
                end if
                
                call find(nametag, "HSE_ES_Eg", i)
                if(i == 0) then
                    write(1224, *) "Input HSE_ES_Eg"
                    HSE_ES_Eg = 1.4
                else
                    read(number(i), *) HSE_ES_Eg
                end if
                
                call find(nametag, "HSE_Es_opt", i)
                if(i == 0) then
                    write(1224, *) "Input HSE_Es_opt"
                    HSE_ES_opt = 3.0
                else
                    read(number(i), *) HSE_Es_opt
                end if
                
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

            call find(nametag, "Layer_hight", i)
            if(i == 0) then
                Layer_hight = 2.0
            else
                read(number(i), *) Layer_hight
            end if
        end if
        
        call find(nametag, "cluster", i)
        if(i == 0) then
            cluster = .false.
        else
            read(number(i), *) cluster
            
            call find(nametag, "model_ball", i)
            if(i == 0) then
                write(1224, *) "Input model_ball"
                model_ball = 0.5
            else
                read(number(i), *) model_ball
            end if
            
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
            
            call find(nametag, "model_shell", i)
            if(i == 0) then
                write(1224, *) "Input model_shell"
                model_shell = 0.3
            else
                read(number(i), *) model_shell
            end if
            
            call find(nametag, "shell_radius_out", i)
            if(i == 0) then
                write(1224, *) "Input shell_radius_out"
                shell_radius_out = 1.0
            else
                read(number(i), *) shell_radius_out
            end if
            
            call find(nametag, "shell_radius_in", i)
            if(i == 0) then
                write(1224, *) "Input shell_radius_in"
                shell_radius_in = 0.8
            else
                read(number(i), *) shell_radius_in
            end if
            
            call find(nametag, "shell_ctr_x", i)
            if(i == 0) then
                write(1224, *) "Input shell_ctr_x"
                shell_ctr_x = 0.5
            else
                read(number(i), *) shell_ctr_x
            end if
            
            call find(nametag, "shell_ctr_y", i)
            if(i == 0) then
                write(1224, *) "Input shell_ctr_y"
                shell_ctr_y = 0.5
            else
                read(number(i), *) shell_ctr_y
            end if
            
            call find(nametag, "shell_ctr_z", i)
            if(i == 0) then
                write(1224, *) "Input shell_ctr_z"
                shell_ctr_z = 0.5
            else
                read(number(i), *) shell_ctr_z
            end if
            
            call find(nametag, "model_plate", i)
            if(i == 0) then
                write(1224, *) "Input model_plate"
                model_plate = 0.2
            else
                read(number(i), *) model_plate
            end if
            
            call find(nametag, "plate_radius", i)
            if(i == 0) then
                write(1224, *) "Input plate_radius"
                plate_radius = 1.0
            else
                read(number(i), *) plate_radius
            end if
            
            call find(nametag, "plate_height", i)
            if(i == 0) then
                write(1224, *) "Input plate_height"
                plate_height = 0.2
            else
                read(number(i), *) plate_height
            end if
            
            call find(nametag, "plate_ctr_x", i)
            if(i == 0) then
                write(1224, *) "Input plate_ctr_x"
                plate_ctr_x = 0.5
            else
                read(number(i), *) plate_ctr_x
            end if
            
            call find(nametag, "plate_ctr_y", i)
            if(i == 0) then
                write(1224, *) "Input plate_ctr_y"
                plate_ctr_y = 0.5
            else
                read(number(i), *) plate_ctr_y
            end if
            
            call find(nametag, "plate_ctr_z", i)
            if(i == 0) then
                write(1224, *) "Input plate_ctr_z"
                plate_ctr_z = 0.5
            else
                read(number(i), *) plate_ctr_z
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
        
        call find(nametag, "model_surface", i)
        if(i == 0) then
            model_surface = .false.
        else
            read(number(i), *) model_surface
            
            call find(nametag, "surface_height", i)
            if(i == 0) then
                write(1224, *) "input surface_height"
                surface_height = 0.5
            else
                read(number(i), *) surface_height
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
        
        call find(nametag, "SelectiveDynamics", i)
        if(i == 0) then
            selective_dynamics = .false.
        else
            read(number(i), *) selective_dynamics
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
        
        call find(nametag, "dimer_mode", i)
        if(i == 0) then
            dimer_mode = .false.
        else
            read(number(i), *) dimer_mode
            
            call find(nametag, "dimer_dis", i)
            if(i == 0) then
                dimer_dis = 1.6
            else
                read(number(i), *) dimer_dis
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
        write(1224, *) "Mystruct: ", Mystruct
        write(1224, *) "Pickup: ", Pickup
        if(Pickup) then
            write(1224, *) "Pickup_step: ", Pickup_step
        end if
        write(1224, *) "Symmetry: ", symmetry
        write(1224, *) "spg_front: ", spg_front
        write(1224, *) "spg_rear: ", spg_rear
        write(1224, *) "Q2D: ", Q2D
        if(Q2D) then
            write(1224, *) "vacuum_layer: ", vacuum_layer
            write(1224, *) "Area: ", Area
            write(1224, *) "Layer_hight: ", Layer_hight
        end if
        write(1224, *) "Multi-objective: ", mode
        if(mode) then
            write(1224, *) "hardness: ", hardness
            if(hardness) then
                write(1224, *) "rcut", rcut
                write(1224, *) "ionicity", ionicity
            end if
            write(1224, *) "ESflag: ", ESflag
            if(ESflag) then
                write(1224, *) "ES_mod: ", ES_mod
                write(1224, *) "ES_Eg: ", ES_Eg
                write(1224, *) "ES_opt: ", ES_opt
                if(HSE) then
                    write(1224, *) "HSE_population: ", HSE_population
                    write(1224, *) "LDA_population: ", LDA_population
                    write(1224, *) "energy_cut: ", energy_cut
                    write(1224, *) "gap_cut: ", gap_cut               
                    write(1224, *) "LDA_ES_Eg: ", LDA_ES_Eg
                    write(1224, *) "LDA_Es_opt: ", LDA_Es_opt
                    write(1224, *) "HSE_ES_Eg: ", HSE_ES_Eg
                    write(1224, *) "HSE_Es_opt: ", HSE_Es_opt
                end if
            end if
        end if
        write(1224, *) "cluster: ", cluster
        if(cluster) then
            write(1224, *) "model_ball: ", model_ball
            write(1224, *) "init_radius: ", init_radius
            write(1224, *) "cluster_ctr_x: ", cluster_ctr_x
            write(1224, *) "cluster_ctr_y: ", cluster_ctr_y
            write(1224, *) "cluster_ctr_z: ", cluster_ctr_z
            write(1224, *) "model_shell: ", model_shell
            write(1224, *) "shell_radius_out: ", shell_radius_out
            write(1224, *) "shell_radius_in: ", shell_radius_in
            write(1224, *) "shell_ctr_x: ", shell_ctr_x
            write(1224, *) "shell_ctr_y: ", shell_ctr_y
            write(1224, *) "shell_ctr_z: ", shell_ctr_z
            write(1224, *) "model_plate: ", model_plate
            write(1224, *) "plate_radius: ", plate_radius
            write(1224, *) "plate_height: ", plate_height
            write(1224, *) "plate_ctr_x: ", plate_ctr_x
            write(1224, *) "plate_ctr_y: ", plate_ctr_y
            write(1224, *) "plate_ctr_z: ", plate_ctr_z
            write(1224, *) "cluster_substrate: ", cluster_substrate
            if(cluster_substrate) then
                write(1224, *) "SubstrateElements: ", (SubstrateElements(i), i = 1, num_species)
            end if
        end if
        
        write(1224, *) "model_surface: ", model_surface
        if(model_surface) then
            write(1224, *) "surface_height: ", surface_height
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
        write(1224, *) "selective_dynamics: ", selective_dynamics
        write(1224, *) "find_defect: ", find_defect
        if(find_defect) then
            write(1224, *) "defect_type: ", defect_type
            write(1224, *) "center: ", center_x, center_y, center_z
        end if
        if(dimer_mode) then
            write(1224, *) "dimer_mode: ", dimer_mode
            write(1224, *) "dimer_dis: ", dimer_dis
        end if
        write(1224, *) "---end write input---"
        
    end subroutine write_input
        
end module init_mod
    
