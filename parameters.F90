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
        logical :: SelectiveDynamics(3, max_atom * max_type)
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
    logical :: cluster
    logical :: selective_dynamics
    logical :: Pickup
    integer :: Pickup_step
    integer :: Mystruct ! number of structs input by my own
    
    !hardness
    logical :: hardness
    real(dp) :: rcut
    real(dp) :: ionicity
    
    !Es_mod(bandgap)
    logical :: ESflag
    integer(i4b) :: ES_mod
    real(dp) :: ES_Eg
    real(dp) :: Es_opt    ! works if and only if ES_mod=6, set the highest boundary for the optical gap
    !USE HSE to calculate Gap
    logical :: HSE
    integer(i4b) :: HSE_population, LDA_population
    real(dp) :: energy_cut, gap_cut
    real(dp) :: LDA_ES_Eg, LDA_Es_opt, HSE_ES_Eg, HSE_Es_opt
    
    !cluster
    !If cluster mode is used, lattice parameters must be fixed
    logical :: cluster_substrate ! add substrate
    integer(i4b) :: SubstrateElements(max_type)
    type(struct_info) :: substrate
    real(dp) :: model_ball, model_shell, model_plate
    real(dp) :: init_radius
    real(dp) :: cluster_ctr_x, cluster_ctr_y, cluster_ctr_z
    real(dp) :: shell_radius_in, shell_radius_out
    real(dp) :: shell_ctr_x, shell_ctr_y, shell_ctr_z
    real(dp) :: plate_radius, plate_height
    real(dp) :: plate_ctr_x, plate_ctr_y, plate_ctr_z
    
    !fix-lattice(control & parameters)
    logical :: fix_lat
    real(dp) :: fix_a, fix_b, fix_c, fix_alpha, fix_beta, fix_gama
    
    !for high pressure system
    logical :: PRESSURE
    integer(i4b) :: PSTRESS
    
    !Quansi-2D(control & parameters)
    logical :: Q2D
    real(dp) :: vacuum_layer, Area, Layer_hight
    
    !defect structure
    logical :: find_defect
    type(struct_info) :: defect_struct
    integer(i4b) :: defect_type ! geometry type of defect, 1: box, 2: spherecal
    real(dp) :: center_x, center_y, center_z, defect_radius, length_x, length_y, length_z
    
    
end module parameters
    
    

