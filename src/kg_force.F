!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2000  CP2K developers group                                !
!-----------------------------------------------------------------------------!
!!****** cp2k/kg_force [1.0] *
!!
!!   NAME
!!     kg_force
!!
!!   FUNCTION
!!     kg force calculation routine      
!!
!!   AUTHOR
!!     gloria,30.09.2002 
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE kg_force

  USE atomic_kind_types, ONLY: atomic_kind_type, get_atomic_kind,&
                               get_atomic_kind_set
  USE coefficient_types, ONLY: coeff_allocate, coeff_copy, &
       coeff_deallocate, &
       coeff_transform_space,&
       coeff_type,&
       coeff_zero,&
       coeff_sumup
  USE cube_utils,               ONLY: cube_info_type,&
       destroy_cube_info,&
       init_cube_info
  USE dft_types, ONLY : dft_control_type 
  USE external_potential_types, ONLY: epc_potential_type,&
                                      get_potential
  USE greens_fn,                ONLY: pw_green_fn_init
  USE gaussian_gridlevels, ONLY: gridlevel_info_type, &
      init_gaussian_gridlevel, destroy_gaussian_gridlevel
  USE global_types, ONLY : global_environment_type
  USE hartree,                  ONLY: calculate_hartree
  USE l_utils,                  ONLY: destroy_l_info,&
       init_l_info,&
       l_info_type
  USE kg_core, ONLY: calculate_ecore_rspace, calculate_eself
  USE kg_density, ONLY : calculate_density, calculate_epc_density,&
                         calculate_total_density, calculate_v_rspace_forces,&
                         calculate_epc_rspace_forces 
  USE kg_environment_types, ONLY : get_kg_env, kg_environment_type, &
                                   kg_thermodynamic_type
  USE kg_force_types, ONLY : kg_force_type, allocate_kg_forces, &
                             deallocate_kg_forces, init_kg_forces
  USE kg_kxc, ONLY : calculate_kxc_potential
  USE kinds, ONLY : dbl
  USE mathconstants, ONLY : pi, zero
  USE message_passing, ONLY : mp_sum, mp_max, mp_sync,&
                  mp_range, mp_allgather
  USE molecule_types, ONLY : molecule_structure_type, particle_node_type, &
                             linklist_exclusion, linklist_bonds, &
                             linklist_bends, linklist_torsions
  USE particle_types, ONLY : particle_type
  USE pw_grid_types, ONLY : pw_grid_type, HALFSPACE
  USE pw_grids, ONLY : pw_find_cutoff, pw_grid_construct, &   
                       pw_grid_change, pw_grid_setup
  USE pws, ONLY: init_pw_poisson_solver, pw_poisson_solver 
  USE pw_types,  ONLY: COMPLEXDATA1D, COMPLEXDATA3D,&
       REALDATA3D, REALSPACE, RECIPROCALSPACE, pw_type
  USE simulation_cell, ONLY : cell_type, get_hinv, get_cell
  USE termination, ONLY : stop_program, stop_memory
  USE timings, ONLY : timeset, timestop

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: kg_force_control, debug_variables_type, get_grids, &
            kg_force_release

! debug variables
  TYPE debug_variables_type
     REAL ( dbl ) :: pot_bond_corr, pot_hartree, pot_xc, pot_rspace_core
     REAL ( dbl ), DIMENSION ( :, : ), POINTER :: f_hartree, &
           f_rspace_core, f_bc, f_xc
  END TYPE debug_variables_type

! definition of all variables (densities, potentials, grids....)
! needed to perform a force calculation 
  TYPE kg_rho_v_type
     ! definition of all grid types used (working grids and multigrids)
     TYPE(pw_grid_type), DIMENSION(:), POINTER :: pw_grid
     ! number of the grid used as reference grid for the density
     integer :: auxbas_grid
     logical :: initialized     
     REAL(dbl) :: total_rho_gspace, total_rho_core_rspace, &
                  total_rho_rspace  
     TYPE(coeff_type)                             :: rho_gspace, &
                                                      v_gspace, &
                                                      v_rspace, &
                                                      v_xc_rspace
     ! multigrids in qs-style
     TYPE(coeff_type), DIMENSION(:), POINTER        :: mgrid_gspace, &
                                                      mgrid_rspace
     ! rho_elec_rspace is supposed to contain the electron density (for cube output)
     TYPE(coeff_type)                :: rho_elec_rspace
     ! hartree is supposed to contain the hartree potential (for cube ouput)
     TYPE(coeff_type)          :: v_hartree_rspace
     ! the core density (in KG, density of the Zeff nuclear charges)
     TYPE(coeff_type)          :: rho_core
     TYPE(gridlevel_info_type) :: gridlevel_info
     TYPE(l_info_type)         :: l_info
     TYPE(cube_info_type), DIMENSION(:), POINTER :: cube_info
  END TYPE kg_rho_v_type
!!***
  TYPE ( kg_rho_v_type ), POINTER :: kg_rho_v 
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** kg_force/kg_force_control [1.1] *
!!
!!   NAME
!!     kg_force_control
!!
!!   FUNCTION
!!     Calculates the total potential energy, total force, and the
!!     total pressure tensor from the potentials
!!
!!   AUTHOR
!!     gt         
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE kg_force_control ( kg_env, globenv, debug )

  IMPLICIT NONE

! Arguments
  TYPE ( kg_environment_type ), INTENT ( INOUT ):: kg_env
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
  TYPE ( debug_variables_type ), INTENT ( OUT ), OPTIONAL :: debug

! Local parameters
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE kg_force_control (MODULE kg_force)"
! Locals
  TYPE ( atomic_kind_type ), DIMENSION(:), POINTER:: atomic_kind_set
  TYPE ( atomic_kind_type ), POINTER   :: atomic_kind
  TYPE ( cell_type ), POINTER :: box
  TYPE ( cell_type ), POINTER :: box_ref
  TYPE ( coeff_type ), DIMENSION (:), POINTER :: coef_pos
  TYPE ( coeff_type ), DIMENSION (:), POINTER :: coef_vel
  TYPE ( coeff_type ), DIMENSION (:), POINTER :: coef_force
  TYPE ( dft_control_type ), POINTER :: dft_control
  TYPE ( molecule_structure_type ), POINTER :: molecule ( : )
  TYPE ( particle_type ), POINTER :: part ( : )
  TYPE ( particle_node_type ), POINTER :: pnode ( : )
  TYPE ( kg_thermodynamic_type ), POINTER :: thermo
  TYPE ( kg_force_type ), DIMENSION(:), POINTER ::kg_forces
  TYPE(pw_grid_type), DIMENSION(:), POINTER :: pw_multi_grid
  INTEGER :: id, i, ii, natoms, nnodes, handle, isos, n, nkind
  INTEGER :: iatom, ikind 
  REAL ( dbl ) :: pot_nonbond, pot_bond, pot_bend, pot_torsion, vg_coulomb
  REAL ( dbl ) :: ehartree, exc
  REAL ( dbl ) :: total_rho0_rspace
  REAL ( dbl ), DIMENSION ( 3,3 ) :: pv_nonbond, pv_bond, pv_bend, pv_torsion
  REAL ( dbl ), DIMENSION ( :,: ), ALLOCATABLE :: fg_coulomb
  REAL ( dbl ), DIMENSION ( :,: ), ALLOCATABLE :: f_total
  REAL ( dbl ), DIMENSION ( 3,3 ) :: pv_g, ident
  REAL ( dbl ), DIMENSION ( 3,3 ) :: pv_bc
  LOGICAL :: first_time
  INTEGER, DIMENSION(:), POINTER   :: atom_list
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: natom_of_kind
  CHARACTER(LEN=8) :: force_type

! associcating local pointers
  molecule => kg_env % molecule
  part => kg_env % part
  pnode => kg_env % pnode
  box => kg_env % box
  box_ref => kg_env % box_ref

!------------------------------------------------------------------------------

  CALL timeset ( 'FORCE','I',' ',handle )

  CALL get_kg_env( kg_env=kg_env, atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control, thermo=thermo )

  first_time=.FALSE. 
  nnodes = SIZE ( pnode )
  natoms = SIZE ( part )
  nkind = SIZE( atomic_kind_set )
  isos = 0

! allocates forces                    
  IF ( .NOT. ALLOCATED ( natom_of_kind ) ) THEN
     ALLOCATE (natom_of_kind(nkind),STAT=isos)
     IF (isos /= 0) CALL stop_memory(routine,"natom_of_kind",nkind)
     CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                              natom_of_kind=natom_of_kind)
     first_time = .TRUE.
  END IF

! initialization       

  IF ( first_time ) THEN

     ALLOCATE ( kg_rho_v, STAT = isos )
     IF ( isos /= 0 )  CALL stop_memory ( 'kg_force_control', 'kg_rho_v', 0 )

     CALL init_kg_force_calculation ( kg_rho_v, kg_env, globenv )
 
     CALL get_hinv ( box )
 
     CALL calculate_eself ( kg_env, thermo % e_self )

  END IF

! allocate density/pot arrays
 
  CALL allocate_kg_densities (kg_rho_v, globenv )

! allocate forces

  NULLIFY (kg_forces)
  CALL allocate_kg_forces( kg_forces, natom_of_kind )

! initialize forces and energies                    
 
  CALL init_kg_forces(kg_forces)

  DO i = 1, natoms
     part ( i ) % f ( 1 ) = 0.0_dbl
     part ( i ) % f ( 2 ) = 0.0_dbl
     part ( i ) % f ( 3 ) = 0.0_dbl
  END DO
  thermo % pot = 0.0_dbl
  thermo % e_hartree = 0.0_dbl
  thermo % e_rspace_core = 0.0_dbl
  thermo % e_bond_corr = 0.0_dbl
  thermo % e_xc = 0.0_dbl
  thermo % virial % pv_virial = 0.0_dbl

! Calculate the rspace energy of the core charge distribution
! and its force contributions

  CALL calculate_ecore_rspace(kg_env,globenv,thermo%e_rspace_core,kg_forces)

! Calculate the density generated by the effective point core charges qeff

  CALL calculate_epc_density(kg_rho_v%rho_core, kg_rho_v%cube_info(1),&
                         kg_rho_v%l_info,&
                         kg_rho_v%total_rho_core_rspace, kg_env, globenv)

! calculate the frozen electronic pseudo density

  CALL calculate_density ( kg_rho_v%rho_elec_rspace, kg_rho_v%rho_gspace, &
            kg_rho_v%mgrid_rspace, kg_rho_v%mgrid_gspace,&
            kg_rho_v%gridlevel_info, kg_rho_v%cube_info,&
            kg_rho_v%l_info, total_rho0_rspace, kg_env, globenv )

  CALL coeff_zero(kg_rho_v%v_rspace)
! from the frozen rho, calculate the thomas-fermi and xc potentials
  CALL calculate_kxc_potential (kg_rho_v % rho_elec_rspace, kg_rho_v % rho_gspace, &
                                kg_rho_v%v_xc_rspace, thermo % e_xc, dft_control)

! calculate the total density
  kg_rho_v%total_rho_rspace = total_rho0_rspace+&
             kg_rho_v%total_rho_core_rspace

! add the core density to the frozen electronic density

  CALL coeff_sumup(kg_rho_v%rho_core,kg_rho_v%rho_gspace)

  kg_rho_v%total_rho_gspace = calculate_total_density(kg_rho_v%rho_gspace)

! calculate electrostatic potential
  CALL calculate_hartree(kg_rho_v%rho_gspace, thermo % e_hartree,&
           kg_rho_v%v_gspace)

!  calclation of the forces on the ions

  CALL coeff_transform_space(kg_rho_v%v_gspace,kg_rho_v%v_hartree_rspace)

  kg_rho_v%v_hartree_rspace%pw%cr3d(:,:,:) =&
               kg_rho_v%v_hartree_rspace%pw%pw_grid%dvol*&
               kg_rho_v%v_hartree_rspace%pw%cr3d(:,:,:)

  CALL calculate_epc_rspace_forces(kg_rho_v % v_hartree_rspace,&
                   kg_rho_v%cube_info(1),&
                   kg_rho_v%l_info, kg_env, kg_forces)

! if we are doing debugging calculate separately the contribution
! to the forces due by tfxc and hartree 
! else, sumup all the realspace potentials and get the total
! contribution to the forces

  IF ( PRESENT ( debug )  ) THEN
    force_type = 'hartree'
    CALL calculate_v_rspace_forces(kg_rho_v%v_hartree_rspace,kg_rho_v%mgrid_rspace,&
               kg_rho_v%mgrid_gspace, kg_rho_v%auxbas_grid,&
               kg_rho_v%gridlevel_info,&
               kg_rho_v%cube_info, kg_rho_v%l_info,&
               kg_env, kg_forces, globenv, force_type)
    force_type = 'xc'
    CALL calculate_v_rspace_forces(kg_rho_v%v_xc_rspace,kg_rho_v%mgrid_rspace,&
               kg_rho_v%mgrid_gspace, kg_rho_v%auxbas_grid,&
               kg_rho_v%gridlevel_info,&
               kg_rho_v%cube_info, kg_rho_v%l_info,&
               kg_env, kg_forces, globenv, force_type)
  ELSE 
    kg_rho_v%v_rspace%pw%cr3d(:,:,:) =&
               kg_rho_v%v_xc_rspace%pw%cr3d(:,:,:) +&
               kg_rho_v%v_hartree_rspace%pw%cr3d(:,:,:)

    CALL calculate_v_rspace_forces(kg_rho_v%v_rspace,kg_rho_v%mgrid_rspace,&
               kg_rho_v%mgrid_gspace, kg_rho_v%auxbas_grid,&
               kg_rho_v%gridlevel_info,&
               kg_rho_v%cube_info, kg_rho_v%l_info,&
               kg_env, kg_forces, globenv)
  END IF

! add up all the potential energies

  thermo % pot = thermo % e_self + thermo % e_rspace_core + thermo % e_xc + &
                 thermo % e_hartree + thermo % e_bond_corr 
     
! add up all the forces

  DO ikind = 1, nkind         
    kg_forces(ikind)%f_total(:,:) = kg_forces(ikind)%f_rho(:,:) +&
                                kg_forces(ikind)%f_rspace_core(:,:) +&
                                kg_forces(ikind)%f_hartree_core(:,:) 
  END DO


  DO ikind = 1, nkind         
    atomic_kind => atomic_kind_set(ikind)
    CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           atom_list=atom_list) 
    DO iatom=1, SIZE(atom_list)
      i = atom_list(iatom)
! fill in the forces in part
      part(i)%f(1) = kg_forces(ikind)%f_total(1,iatom)
      part(i)%f(2) = kg_forces(ikind)%f_total(2,iatom)
      part(i)%f(3) = kg_forces(ikind)%f_total(3,iatom)
! if we are doing debugging, check if variables are present and assign
      IF ( PRESENT ( debug )  ) THEN
        debug % f_hartree (:,i) = kg_forces(ikind)%f_hartree(:,iatom) + &
                              kg_forces(ikind)%f_hartree_core(:,iatom)
        debug % f_rspace_core (:,i) = kg_forces(ikind)%f_rspace_core(:,iatom)
        debug % f_bc (:,i) = kg_forces(ikind)%f_bond_corr(:,iatom)
        debug % f_xc (:,i) = kg_forces(ikind)%f_xc(:,iatom)
      END IF
    END DO
  END DO

  IF ( PRESENT ( debug )  ) THEN
    debug % pot_hartree = thermo % e_hartree
    debug % pot_rspace_core = thermo % e_rspace_core
    debug % pot_bond_corr = thermo % e_bond_corr
    debug % pot_xc = thermo % e_xc
  END IF

#if defined ( __parallel )
!  CALL mp_sum ( f_total, globenv % group )
#endif

! deallocating all local variables

 CALL deallocate_kg_forces( kg_forces )
 CALL deallocate_kg_densities ( kg_rho_v, globenv )

  IF ( ALLOCATED ( f_total ) ) THEN
    DEALLOCATE ( f_total, STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory ( 'kg_force_control', 'f_total' )
  END IF

  CALL timestop ( zero, handle )

END SUBROUTINE kg_force_control

!!*** **********************************************************************
  SUBROUTINE allocate_kg_densities(kg_rho_v, globenv )

! Arguments
    TYPE(kg_rho_v_type), INTENT(OUT) :: kg_rho_v
    TYPE(global_environment_type), intent(in), target :: globenv

! Locals
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE allocate_kg_densities (MODULE kg_force)"
    TYPE(pw_grid_type), POINTER, DIMENSION (:) :: pw_grid 
    INTEGER :: auxbas, igrid_level, istat, ngrid_level
    INTEGER :: handle

!---------------------------------------------------------------------------

    CALL timeset("init_grid","I","",handle)

    NULLIFY(kg_rho_v%mgrid_gspace, kg_rho_v%mgrid_rspace)    

    auxbas = kg_rho_v % auxbas_grid
    pw_grid => kg_rho_v % pw_grid
    ngrid_level = SIZE ( pw_grid )

! Allocate the rho's on the largest reference grid

    CALL coeff_allocate(kg_rho_v%rho_elec_rspace,pw_grid(auxbas),REALDATA3D)
    kg_rho_v%rho_elec_rspace%pw%in_space = REALSPACE
    CALL coeff_zero (kg_rho_v%rho_elec_rspace)

    CALL coeff_allocate(kg_rho_v%rho_gspace,pw_grid(auxbas),COMPLEXDATA1D)
    kg_rho_v%rho_gspace%pw%in_space = RECIPROCALSPACE
    CALL coeff_zero (kg_rho_v%rho_gspace)

    CALL coeff_allocate(kg_rho_v%v_xc_rspace,pw_grid(auxbas),REALDATA3D)
    kg_rho_v%v_xc_rspace%pw%in_space = REALSPACE
    CALL coeff_zero (kg_rho_v%v_xc_rspace)

    CALL coeff_allocate(kg_rho_v%v_rspace,pw_grid(auxbas),REALDATA3D)
    kg_rho_v%v_rspace%pw%in_space = REALSPACE
    CALL coeff_zero (kg_rho_v%v_rspace)

    CALL coeff_allocate(kg_rho_v%v_gspace,pw_grid(auxbas),COMPLEXDATA1D)
    kg_rho_v%v_gspace%pw%in_space = RECIPROCALSPACE
    CALL coeff_zero (kg_rho_v%v_gspace)

    CALL coeff_allocate(kg_rho_v%v_hartree_rspace,pw_grid(auxbas),REALDATA3D)
    kg_rho_v%v_hartree_rspace%pw%in_space = REALSPACE
    CALL coeff_zero (kg_rho_v%v_hartree_rspace)

    CALL coeff_allocate(kg_rho_v%rho_core,pw_grid(auxbas),COMPLEXDATA1D)
    kg_rho_v%rho_core%pw%in_space = RECIPROCALSPACE
    CALL coeff_zero (kg_rho_v%rho_core)

! allocate the multi-grid rho's
    ALLOCATE (kg_rho_v%mgrid_rspace(ngrid_level),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"mgrid_rspace",0)
    ALLOCATE (kg_rho_v%mgrid_gspace(ngrid_level),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"mgrid_gspace",0)

! setup the multigrids
    DO igrid_level=1,ngrid_level
       CALL coeff_allocate(kg_rho_v%mgrid_rspace(igrid_level),pw_grid(igrid_level),&
            REALDATA3D)
       kg_rho_v%mgrid_rspace(igrid_level)%pw%in_space = REALSPACE
       CALL coeff_zero (kg_rho_v % mgrid_rspace(igrid_level))
       CALL coeff_allocate(kg_rho_v%mgrid_gspace(igrid_level),pw_grid(igrid_level),&
            COMPLEXDATA1D)
       kg_rho_v%mgrid_gspace(igrid_level)%pw%in_space = RECIPROCALSPACE
       CALL coeff_zero (kg_rho_v % mgrid_gspace(igrid_level))
    END DO

    CALL timestop(0.0_dbl,handle)

END SUBROUTINE allocate_kg_densities

!!*** **********************************************************************
  SUBROUTINE deallocate_kg_densities(kg_rho_v, globenv )

! Arguments
    TYPE(kg_rho_v_type), INTENT(OUT) :: kg_rho_v
    TYPE(global_environment_type), intent(in), target :: globenv

! Locals
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE deallocate_kg_densities (MODULE kg_rho_v)"
    TYPE(pw_grid_type), POINTER, DIMENSION (:) :: pw_grid 
    INTEGER :: auxbas, igrid_level, istat, ngrid_level
    INTEGER :: handle

!---------------------------------------------------------------------------

    CALL timeset("init_grid","I","",handle)

    auxbas = kg_rho_v % auxbas_grid
    pw_grid => kg_rho_v % pw_grid
    ngrid_level = SIZE ( pw_grid )

! DeAllocate the rho's on the largest reference grid

    CALL coeff_deallocate(kg_rho_v%rho_elec_rspace)

    CALL coeff_deallocate(kg_rho_v%rho_gspace)

    CALL coeff_deallocate(kg_rho_v%v_xc_rspace)

    CALL coeff_deallocate(kg_rho_v%v_rspace)

    CALL coeff_deallocate(kg_rho_v%v_gspace)

    CALL coeff_deallocate(kg_rho_v%v_hartree_rspace)

    CALL coeff_deallocate(kg_rho_v%rho_core)


! deallocate the multi-grid rho's

    DO igrid_level=1,ngrid_level
       CALL coeff_deallocate(kg_rho_v%mgrid_rspace(igrid_level))
       CALL coeff_deallocate(kg_rho_v%mgrid_gspace(igrid_level))
    END DO

    DEALLOCATE (kg_rho_v%mgrid_rspace, STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"mgrid_rspace",0)
    DEALLOCATE (kg_rho_v%mgrid_gspace,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"mgrid_gspace",0)

    CALL timestop(0.0_dbl,handle)

END SUBROUTINE deallocate_kg_densities

!!*** **********************************************************************
!!**** kg_force/init_kg_force_calculation[1.0] *
!!
!!   NAME
!!     init_kg_forc_calculation
!!
!!   FUNCTION
!!     initializes the environment needed for calculation of kg forces
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - kg_rho_v: the density/potential env to be initialized
!!     - kg_env: the kg environment
!!     - globenv: the global environment
!!
!!   AUTHOR
!!    gloria 
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
  SUBROUTINE init_kg_force_calculation(kg_rho_v, kg_env, globenv )

  IMPLICIT NONE

! Arguments

    TYPE(kg_rho_v_type), INTENT(OUT) :: kg_rho_v
    TYPE(kg_environment_type), INTENT(IN), TARGET :: kg_env
    TYPE(global_environment_type), intent(in), target :: globenv

! Locals

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE init_kg_force_calculation (MODULE kg_force)"

    INTEGER, PARAMETER :: maxgridpoints = 100
    REAL(dbl)          :: rel_cutoff,maxradius
    INTEGER :: ngrid_level
    REAL(dbl), DIMENSION(:), POINTER :: cutoff
    INTEGER  :: cmax,cmaxl,igrid_level, istat,maxl
    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set

!   ---------------------------------------------------------------------------

    NULLIFY( kg_rho_v%cube_info )


    rel_cutoff =  kg_env%dft_control%kg_control%relative_cutoff
    cutoff => kg_env % dft_control%kg_control%e_cutoff
    atomic_kind_set => kg_env % atomic_kind_set
    ngrid_level = size(cutoff)

    CALL init_gaussian_gridlevel(kg_rho_v%gridlevel_info,&
         ngrid_level,cutoff,rel_cutoff)

    CALL init_grid (kg_rho_v,kg_env,globenv=globenv, &
         ngrid_level=ngrid_level)

    ALLOCATE (kg_rho_v%cube_info(ngrid_level),STAT=istat)
    if (istat /= 0) call stop_memory(routine,"cube_info",0)

    cmax = 0
    do igrid_level=1,ngrid_level
       maxradius = maxgridpoints*maxval(kg_rho_v%pw_grid(igrid_level)%dr)
       call init_cube_info(kg_rho_v%cube_info(igrid_level),&
            kg_rho_v%pw_grid(igrid_level)%dr(:),&
            maxradius,cmaxl)
       cmax = max(cmax,cmaxl)
    end do

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,maxl=maxl)

    ! maxl+1 to allow for the forces
    CALL init_l_info(kg_rho_v%l_info,maxl+1,cmax)

    kg_rho_v%total_rho_core_rspace=0.0_dbl
    kg_rho_v%total_rho_rspace=0.0_dbl
    kg_rho_v%total_rho_gspace=0.0_dbl
    kg_rho_v%initialized=.true.

  END SUBROUTINE init_kg_force_calculation
!***************************************************************************

!!****kg_force/release_kg_rho_v [1.0] *
!!
!!   NAME
!!     release_kg_rho_v
!!
!!   FUNCTION
!!     releases the memory used by the kg_rho_v (but not kg_rho_v)
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
  SUBROUTINE release_kg_rho_v(kg_rho_v,kg_env,globenv)

  IMPLICIT NONE

!arguments
    TYPE(kg_rho_v_type), INTENT(INOUT) :: kg_rho_v
    TYPE(kg_environment_type), INTENT(IN),TARGET :: kg_env
    TYPE(global_environment_type), INTENT(IN), TARGET :: globenv

!Locals
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE release_kg_rho_v (MODULE kg_force)"

    INTEGER :: igrid_level,istat

!   ---------------------------------------------------------------------------
 
    DO igrid_level=1,size(kg_rho_v%cube_info)
      CALL destroy_cube_info(kg_rho_v%cube_info(igrid_level))
    END DO

    DEALLOCATE (kg_rho_v%cube_info,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"cube_info")

    CALL destroy_l_info(kg_rho_v%l_info)
    CALL destroy_gaussian_gridlevel(kg_rho_v%gridlevel_info,globenv)

    kg_rho_v%initialized=.false.
  END SUBROUTINE release_kg_rho_v

!***************************************************************************

  SUBROUTINE init_grid ( kg_rho_v, kg_env, globenv, ngrid_level )

  IMPLICIT NONE

!   Arguments

    TYPE(kg_rho_v_type), INTENT(INOUT) :: kg_rho_v
    TYPE(kg_environment_type), INTENT(IN) :: kg_env
    TYPE(global_environment_type), INTENT(IN), TARGET :: globenv
    INTEGER, INTENT(IN) :: ngrid_level

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE init_grid (MODULE kg_force)"

!   *** Local variables ***

    INTEGER :: handle,i,igrid_level,istat,output_unit ! ,grid_info
    INTEGER :: nrhos, auxbas
    LOGICAL :: ionode
    TYPE(cell_type), POINTER          :: cell
    type(pw_grid_type), dimension(:), pointer :: pw_grid
    REAL(dbl), DIMENSION(:), POINTER :: cutoff

!   ---------------------------------------------------------------------------

    CALL timeset("init_grid","I","",handle)

    ionode = globenv % ionode
    output_unit = globenv % scr
    cell => kg_env % box
    cutoff => kg_env % dft_control % kg_control % e_cutoff

    ! the density grid coincides with the multi-grid having the highest cutoff
    ! auxbas_grid is the index of the auxiliary (REFERENCE) grid

    DO i=1,ngrid_level
      IF ( cutoff(i) == kg_env%dft_control%kg_control%cutoff ) kg_rho_v % auxbas_grid = i
    END DO

!   *** Allocate the multi-grid pointer pw_grid  ***

    ALLOCATE (kg_rho_v%pw_grid(ngrid_level),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pw_grid",0)

!   *** Construct the various grids ***

    pw_grid => kg_rho_v%pw_grid
    DO igrid_level=1,ngrid_level
       CALL pw_grid_construct(pw_grid(igrid_level))
       pw_grid(igrid_level)%grid_span = HALFSPACE
       pw_grid(igrid_level)%para%rs_dims ( 1 ) = globenv % num_pe
       pw_grid(igrid_level)%para%rs_dims ( 2 ) = 1

       IF (igrid_level == 1) THEN
         IF (globenv%print%pw_grid_information) THEN
           CALL pw_grid_setup(cell,pw_grid(igrid_level),cutoff(igrid_level),&
                              pe_group=globenv%group,&
                              info=output_unit)
         ELSE
           CALL pw_grid_setup(cell,pw_grid(igrid_level),cutoff(igrid_level),&
                              pe_group=globenv%group)
         END IF
       ELSE
         IF (globenv%print%pw_grid_information) THEN
           CALL pw_grid_setup(cell,pw_grid(igrid_level),cutoff(igrid_level),&
                              pe_group=globenv%group,&
                              info=output_unit,&
                              ref_grid=pw_grid(igrid_level-1))
         ELSE
           CALL pw_grid_setup(cell,pw_grid(igrid_level),cutoff(igrid_level),&
                              pe_group=globenv%group,&
                              ref_grid=pw_grid(igrid_level-1))
         END IF
       END IF

    END DO

    auxbas = kg_rho_v % auxbas_grid

!   poisson solver initialized only for the reference cutoff grid 
!   Total density and potentials are allocated on the reference grid

    CALL pw_green_fn_init(cell%green,pw_grid(auxbas))

    CALL init_pw_poisson_solver(cell%green)

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE init_grid

! *****************************************************************************

 SUBROUTINE get_grids ( pw_grids, gridlevel_info, l_info, cube_info, ref_grid )

  IMPLICIT NONE

  TYPE ( pw_grid_type ), DIMENSION (:), POINTER, OPTIONAL :: pw_grids
  TYPE ( gridlevel_info_type ), POINTER, OPTIONAL :: gridlevel_info
  TYPE ( l_info_type ), POINTER, OPTIONAL         :: l_info
  TYPE ( cube_info_type ), DIMENSION(:), POINTER, OPTIONAL :: cube_info
  INTEGER, OPTIONAL, INTENT(OUT) :: ref_grid


  IF ( PRESENT ( pw_grids ) ) pw_grids => kg_rho_v % pw_grid
  IF ( PRESENT ( gridlevel_info ) ) gridlevel_info => kg_rho_v % gridlevel_info
  IF ( PRESENT ( l_info ) ) l_info => kg_rho_v % l_info
  IF ( PRESENT ( cube_info ) ) cube_info => kg_rho_v % cube_info
  IF ( PRESENT ( ref_grid ) ) ref_grid =  kg_rho_v % auxbas_grid   

END SUBROUTINE get_grids
!!*** **********************************************************************
!!**** kg_force/kg_force_release[1.0] *
!!
!!   NAME
!!     kg_force_release
!!
!!   FUNCTION
!!     before exiting the run, clean all the allocated/associated     
!!     stuff in kg_force.
!!
!!   AUTHOR
!!    gloria
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
  SUBROUTINE kg_force_release(kg_env, globenv )

  IMPLICIT NONE
! arguments
  TYPE ( kg_environment_type ), INTENT ( INOUT ):: kg_env
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
!locals
  INTEGER :: isos
! ---------------------------------------------------------------------------
  CALL release_kg_rho_v(kg_rho_v,kg_env,globenv)

  IF (ASSOCIATED (kg_rho_v))  DEALLOCATE ( kg_rho_v, STAT = isos )
  IF ( isos /= 0 )  CALL stop_memory ( 'kg_force_release', 'kg_rho_v') 

 END SUBROUTINE  kg_force_release

!******************************************************************************

END MODULE kg_force  

!******************************************************************************
