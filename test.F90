!     
! File:   newfortranFreeFormatFile.F90
! Author: zyy
!
! Created on 2013年4月18日, 下午5:08
!

program test
use kinds
use init_mod, only : init_run
!use lattice_mod
use init_struct_all, only : init_struct_tot
use init_struct_spg, only : init_struct_sym
use de_tools, only: read_vasp, write_vasp, write_pwmat, write_input_all, write_vasp_all
use de_tools, only: write_input_HSE_all, write_vasp_HSE_all, read_vasp_Pickup, read_mystruct
use de_tools, only: bulkmodu, ES_fitness
!use run_lammps
use differencial_evolution, only : run_de_vasp, run_mode_vasp
use sort, only : sort_results_mode
use parameters, only : Pickup, Pickup_step, Mystruct
use parameters, only : PWMAT
use parameters, only : population, pstruct, pool, max_step, max_step, mode, cluster
use parameters, only : hardness, ESflag, HSE
use parameters, only : HSE_population, LDA_population
use parameters, only : ES_Eg, ES_opt, LDA_ES_Eg, LDA_Es_opt, HSE_ES_Eg, HSE_Es_opt
implicit none
integer(i4b) :: a, i, j
integer(i4b) :: step
logical :: flag
real(dp) :: tmp

call init_run()
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
step = Pickup_step
if(Pickup) then
    call read_vasp_Pickup(step, flag)
    if(.not. flag) then
        Pickup = .false.
        write(1224, *) "Pickup failed"
    else
        write(1224, *) "Pickup Succeed"
    end if
end if
if(.not. Pickup) then
    do i = 1, population
        !call init_lat_sym(i)
        write(*, *) "init_struct: ", i
        call init_struct_tot(i)
        !if(cluster) then
        !    call init_struct_cluster(i)
        !else
        !    call init_struct_sym(i)
        !end if
    end do
end if
do i = 1, Mystruct
    call read_mystruct(i, flag)
    if(flag) then
        write(1224, *) "Mystruct: ", i
    else
        write(1224, *) "Failed in reading mystruct: ", i
    end if
end do
step = 1
if(Pickup) step = Pickup_step
do j = step, max_step
    write(1224, "(1X, A6, I5)") "Step =", j
    if(.not. Pickup) then
        do i = 1, population
            call write_vasp(i)
            if(PWMAT) then 
                call write_pwmat(i)
            end if
        end do
        call write_input_all(j)
        write(*, *)"run pvasp"
        if(PWMAT) then
            call system("bash run_pwmat.sh")
        else
            call system("bash run_pvasp.sh")
        end if
        do i = 1, population
            call read_vasp(i)
            if(hardness) then
                call bulkmodu(i)
            end if
            if(ESflag) then
                call ES_fitness(i)
            end if
        end do
        call write_vasp_all(j)
        if(HSE) then
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
            call sort_results_mode()
            population = HSE_population
            ES_Eg = HSE_ES_Eg
            ES_opt = HSE_Es_opt
            do i = 1, population
                call write_vasp(i)
            end do
            call write_input_HSE_all(j)
            write(*, *) "run pvasp HSE"
            call system("bash run_pvasp_HSE.sh")
            do i = 1, population
                call read_vasp(i)
                call ES_fitness(i)
            end do
            call write_vasp_HSE_all(j)
        end if
    end if
    Pickup = .false.
    write(*, *) "de operation"
    if(mode) then
        call run_mode_vasp(j)
    else
        call run_de_vasp(j)
    end if
    if(HSE) then
        write(*, *) "in HSE, generate new structure"
        population = LDA_population
        ES_Eg = LDA_ES_Eg
        ES_opt = LDA_Es_opt
        do i = HSE_population, population
            call init_struct_tot(i)
        end do
    end if
end do
close(1224)
!call write_vasp_all(1)
!call write_vasp_all(2)
    
end program test

