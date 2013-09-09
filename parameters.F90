!     
! File:   parameters.F90
! Author: zyy
!
! Created on 2013年4月19日, 上午9:43
!

module parameters
    
    use kinds
    use constants
    
    implicit none
    public
    
    type :: struct_info
        real(dp) :: lat(3,3)
        integer(i4b) :: nelement(max_type), natom, ntyp
        real(dp) :: pos(3, max_atom * max_type)
        integer(i4b) :: ptype(max_atom * max_type)
        integer(i4b) :: id, spg_idx
        integer(i4b) :: frontier
        logical :: flag
        real(dp) :: energy
        real(dp) :: hardness
        real(dp) :: Eg_id, Eg_d
    end type
    type(struct_info), dimension(max_struct) :: pstruct
    type(struct_info), dimension(max_struct) :: pool
    
    ! input de info
    character(len = 20) :: sys_name
    integer(i4b) :: num_species
    integer(i4b) :: num_ele(max_type)
    character(len = 2) :: name_element(max_type)
    real(dp) :: volumn
    integer(i4b) :: max_step
    integer(i4b) :: population
    real(dp) :: atom_dis(max_type, max_type)
    real(dp) :: de_ratio
    
    ! sys_control_sym
    logical :: symmetry
    integer(i4b) :: spg_front, spg_rear !specify which space group you want to use. range: [1, 230], spg_front <= spg_rear
    integer(i4b) :: spacegroup_log(230)
    integer(i4b) :: spg_index
    logical :: mode
    
    !hardness
    logical :: hardness
    real(dp) :: rcut
    real(dp) :: ionicity
    
    !Es_mod(bandgap)
    logical :: ESflag
    integer(i4b) :: ES_mod
    real(dp) :: ES_Eg
    real(dp) :: Es_opt    ! works if and only if ES_mod=6, set the highest boundary for the optical gap
    
    !fix-lattice(control & parameters)
    logical :: fix_lat
    real(dp) :: fix_a, fix_b, fix_c, fix_alpha, fix_beta, fix_gama
    
    !Quansi-2D(control & parameters)
    logical :: Q2D
    real(dp) :: vacuum_layer
    
    
end module parameters
    
    

