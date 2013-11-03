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
use de_tools, only: read_vasp, write_vasp, write_input_all, write_vasp_all
use de_tools, only: bulkmodu, ES_fitness
!use run_lammps
use differencial_evolution, only : run_de_vasp, run_mode_vasp
!use sort, only : sort_results_mode
use parameters, only : population, pstruct, pool, max_step, max_step, mode, cluster
use parameters, only : hardness, ESflag
implicit none
integer(i4b) :: a, i, j
!real(dp) :: max_e
!a = 1
!write(*, *) a
call init_run()
!do i = 1, population
!    call init_lat_nosym(i)
!end do
!call init_struct_nosym()
!call lammps()
!do i = 1, population
!    pool(i) = pstruct(i)
!end do
!do j = 1, max_step
!    call sort_results()
!    call run_de_lammps(i)
!    call lammps()
!    write(*, *) "step = ", j
!    max_e = 0
!    do i = 1, population
!        if(pstruct(i) % energy < pool(i) % energy) then
!            pool(i) = pstruct(i)
!        end if
!        if(pstruct(i) % energy < max_e) then
!            max_e = pstruct(i) % energy
!        end if
!    end do
!    write(*, *) max_e
!end do
!call print_lmp_res()
!close(1224)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

do j = 1, max_step
    write(1224, "(1X, A6, I5)") "Step =", j
    do i = 1, population
        call write_vasp(i)
    end do
    call write_input_all(j)
    call system("bash run_pvasp.sh")
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
!    call sort_results_mode()
    if(mode) then
        call run_mode_vasp(j)
    else
        call run_de_vasp(j)
    end if
end do
close(1224)
!call write_vasp_all(1)
!call write_vasp_all(2)
    
end program test

