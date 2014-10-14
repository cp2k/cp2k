!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE qs_ks_types
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE atprop_types,                    ONLY: atprop_type
  USE cell_methods,                    ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_release,&
                                             dft_control_retain,&
                                             dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_deallocate_matrix_set,&
                                             cp_dbcsr_p_type
  USE cp_result_types,                 ONLY: cp_result_type
  USE cp_subsys_types,                 ONLY: cp_subsys_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE f77_blas
  USE kpoint_types,                    ONLY: kpoint_create,&
                                             kpoint_release,&
                                             kpoint_type
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types,                  ONLY: particle_type
  USE pw_types,                        ONLY: pw_p_type,&
                                             pw_release,&
                                             pw_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_kind_types,                   ONLY: qs_kind_type
  USE qs_neighbor_list_types,          ONLY: deallocate_neighbor_list_set,&
                                             neighbor_list_set_p_type
  USE qs_rho_types,                    ONLY: qs_rho_create,&
                                             qs_rho_release,&
                                             qs_rho_type
  USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                             qs_subsys_release,&
                                             qs_subsys_retain,&
                                             qs_subsys_type
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ks_types'

  PUBLIC :: qs_ks_env_type, qs_ks_env_create
  PUBLIC :: qs_ks_release, qs_ks_retain, get_ks_env, set_ks_env

! *****************************************************************************
!> \brief calculation environement to calculate the ks matrix,
!>      holds all the needed vars.
!>      assumes that the core hamiltonian and energy are up to date.
!> \param s_mstruct_changed if the sparsity structure of s in qs_env where this
!>        ks_env lives has changed
!> \param rho_changed if the rho stored in the qs_env has changed
!> \param forces_up_to_date if the forces are up to date
!> \param potential_changed if some external potentian changed (such as due to the constraint)
!>
!> \param matrix_h core hamiltonian (matrix_h(1)), if needed by the calculation also
!>        its derivatives wrt. x,y, and z (matrix_h(2:4))
!> \param matrix_ks Khon-Sham matrix (filtred by the structure of S).
!>        with LSD an array with the differents orbitals
!> \param matrix_s the overlap matrix (matrix_s(1)), if needed by the calculation also
!>        its derivatives wrt. x,y, and z (matrix_s(2:4))
!> \param kinetic kinetic part of h
!> \param gamma_matrix gamma (Coulomb) matrix in DFTB, if needed by the calculation also
!>        its derivatives wrt. x,y, and z
!>
!> \param rho the (old) density, in various representations (ao+grid)
!> \param rho_xc the (old) soft density, used for the GAPW_XC method only
!> \param rho_core a grid with the charges of the cores of the atoms in the
!>        reciprocal (g) space
!> \param vppl a realspace grid with the local pseudopotential
!> \param vee a realspace grid with the external electrostatic potential
!>
!> \param neighbor_list_id the current version of the neighbor_list
!> \param sab_orb ,sac_ae, sac_ppl, sap_ppnl:
!>        neighbor lists for the calculation of the core Hamiltonian matrix
!> \param sap_oce: neighbor lists for the calculation of the expansion
!>        coefficients of the local atomic densities rho1_hard and rho1_soft
!> \param sab_se: neighbor lists for the calculation of the 2 centers
!>        hartree term in semi-empirical methods
!> \param sac_tbe: neighbor lists for the calculation of the tight binding
!>        Ewald terms (DFTB)
!> \param sab_core neighbor lists for the calculation of the core interactions
!> \param sab_all neighbor lists for the calculation of the  matrix element of
!>        non-symmetric operators
!> \param sab_vdw: neighbor lists for the calculation of dispersion interactions
!> \param sab_scp: neighbor lists for the calculation of the dispersion interactions
!>        when self-consistent polarization (SCP) is included
!> \param sab_almo: neighbor lists to create ALMO delocalization template
!>
!> \param kpoints information on the kpoints used
!> \param subsys the particles, molecules,... of this environment
!>
!> \par History
!>      04.2002 created [fawzi]
!>      10.2002 introduced pools, uses updated rho as input,
!>              removed most temporary variables, renamed may vars,
!>              began conversion to LSD [fawzi]
!>      10.2014 Moved many components from qs_env here [Ole Schuett]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE qs_ks_env_type
     PRIVATE
     LOGICAL                                               :: s_mstruct_changed = .TRUE.,&
                                                              rho_changed = .TRUE.,&
                                                              potential_changed = .TRUE.,&
                                                              forces_up_to_date = .FALSE.
     INTEGER                                               :: n_evals = 0,&
                                                              id_nr = -1,&
                                                              ref_count = 1

     ! hartree is supposed to contain the hartree potential (for cube ouput)
     ! ugly to keep it always around only for a cube output...
     TYPE(pw_type), POINTER                                :: v_hartree_rspace => Null()
     TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_h => Null(),&
                                                              matrix_ks => Null(),&
                                                              matrix_vxc => Null(),&
                                                              matrix_s => Null(),&
                                                              kinetic => Null(),&
                                                              matrix_w => Null(),&
                                                              matrix_w_mp2 => Null(),&
                                                              matrix_p_mp2 => Null(),&
                                                              gamma_matrix => Null(),&
                                                              matrix_s_aux_fit => Null(),&
                                                              matrix_s_aux_fit_vs_orb => Null(),&
                                                              matrix_ks_aux_fit => Null()


    TYPE(qs_rho_type), POINTER                            :: rho => Null(),&
                                                             rho_buffer => Null(),&
                                                             rho_xc => Null(),&
                                                             rho_aux_fit => Null()
    TYPE(pw_p_type), POINTER                              :: rho_core => Null(),&
                                                             vppl => Null(),&
                                                             rho_nlcc => Null(),&
                                                             rho_nlcc_g => Null(),&
                                                             vee => Null()

    INTEGER                                               :: neighbor_list_id = -1
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb => Null(),&
                                                             sac_ae => Null(),&
                                                             sac_ppl => Null(),&
                                                             sap_ppnl => Null(),&
                                                             sap_oce => Null(),&
                                                             sab_se => Null(),&
                                                             sab_vdw => Null(),&
                                                             sab_tbe => Null(),&
                                                             sab_core => Null(),&
                                                             sab_all => Null(),&
                                                             sab_scp => Null(),&
                                                             sab_aux_fit => Null(),&
                                                             sab_aux_fit_vs_orb => Null(),&
                                                             sab_aux_fit_asymm => Null(),&
                                                             sab_lrc => Null(),&
                                                             sab_almo => Null()

    TYPE(kpoint_type), POINTER                            :: kpoints => Null()
    TYPE(qs_subsys_type), POINTER                         :: subsys => Null()
    TYPE(dft_control_type), POINTER                       :: dft_control => Null()
  END TYPE qs_ks_env_type

! *****************************************************************************
!> \brief type to build arrays of pointers
!> \param ks_env the ks_env pointer
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE qs_ks_env_p_type
     TYPE(qs_ks_env_type), POINTER :: ks_env
  END TYPE qs_ks_env_p_type

CONTAINS

! *****************************************************************************
!> \brief Allocates a new instance of ks_env.
!> \param ks_env ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE qs_ks_env_create(ks_env, error)
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_ks_env_create', &
      routineP = moduleN//':'//routineN

    IF(ASSOCIATED(ks_env)) STOP routineP//"ks_env already associated"
    ALLOCATE(ks_env)
    CALL qs_rho_create(ks_env%rho, error)
    CALL qs_rho_create(ks_env%rho_aux_fit, error)
    CALL qs_rho_create(ks_env%rho_xc, error)
    CALL qs_rho_create(ks_env%rho_buffer, error)
    CALL kpoint_create(ks_env%kpoints,error)
END SUBROUTINE qs_ks_env_create


! *****************************************************************************
!> \brief ...
!> \param ks_env ...
!> \param v_hartree_rspace ...
!> \param s_mstruct_changed ...
!> \param rho_changed ...
!> \param potential_changed ...
!> \param forces_up_to_date ...
!> \param matrix_h ...
!> \param matrix_ks ...
!> \param matrix_vxc ...
!> \param matrix_ks_aux_fit ...
!> \param kinetic ...
!> \param matrix_s ...
!> \param matrix_s_aux_fit ...
!> \param matrix_s_aux_fit_vs_orb ...
!> \param matrix_w ...
!> \param matrix_w_mp2 ...
!> \param matrix_p_mp2 ...
!> \param gamma_matrix ...
!> \param rho ...
!> \param rho_aux_fit ...
!> \param rho_buffer ...
!> \param rho_xc ...
!> \param vppl ...
!> \param rho_core ...
!> \param rho_nlcc ...
!> \param rho_nlcc_g ...
!> \param vee ...
!> \param neighbor_list_id ...
!> \param sab_orb ...
!> \param sab_aux_fit ...
!> \param sab_aux_fit_asymm ...
!> \param sab_aux_fit_vs_orb ...
!> \param sab_all ...
!> \param sac_ae ...
!> \param sac_ppl ...
!> \param sap_ppnl ...
!> \param sap_oce ...
!> \param sab_lrc ...
!> \param sab_se ...
!> \param sab_tbe ...
!> \param sab_core ...
!> \param sab_scp ...
!> \param sab_vdw ...
!> \param sab_almo ...
!> \param kpoints ...
!> \param do_kpoints ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param cell ...
!> \param cell_ref ...
!> \param use_ref_cell ...
!> \param particle_set ...
!> \param energy ...
!> \param force ...
!> \param local_particles ...
!> \param local_molecules ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param subsys ...
!> \param cp_subsys ...
!> \param virial ...
!> \param results ...
!> \param atprop ...
!> \param nkind ...
!> \param natom ...
!> \param dft_control ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE get_ks_env(ks_env,v_hartree_rspace,&
                        s_mstruct_changed,rho_changed,&
                        potential_changed,forces_up_to_date,&
                        matrix_h,matrix_ks,matrix_vxc,matrix_ks_aux_fit,&
                        kinetic,matrix_s,matrix_s_aux_fit,matrix_s_aux_fit_vs_orb,&
                        matrix_w,matrix_w_mp2,matrix_p_mp2,gamma_matrix,&
                        rho, rho_aux_fit, rho_buffer, rho_xc,&
                        vppl, rho_core, rho_nlcc, rho_nlcc_g, vee,&
                        neighbor_list_id,&
                        sab_orb,sab_aux_fit,sab_aux_fit_asymm, &
                        sab_aux_fit_vs_orb,sab_all,sac_ae,sac_ppl,sap_ppnl,sap_oce,sab_lrc,&
                        sab_se,sab_tbe,sab_core,sab_scp,sab_vdw,sab_almo,&
                        kpoints, do_kpoints,&
                        atomic_kind_set,qs_kind_set,cell,cell_ref,use_ref_cell,&
                        particle_set,energy,force,local_particles,local_molecules,&
                        molecule_kind_set,molecule_set,subsys,cp_subsys,virial,results,atprop,&
                        nkind, natom, dft_control,&
                        error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(pw_type), OPTIONAL, POINTER         :: v_hartree_rspace
    LOGICAL, OPTIONAL                        :: s_mstruct_changed, &
                                                rho_changed, &
                                                potential_changed, &
                                                forces_up_to_date
    TYPE(cp_dbcsr_p_type), DIMENSION(:), OPTIONAL, POINTER :: matrix_h, &
      matrix_ks, matrix_vxc, matrix_ks_aux_fit, kinetic, matrix_s, &
      matrix_s_aux_fit, matrix_s_aux_fit_vs_orb, matrix_w, matrix_w_mp2, &
      matrix_p_mp2, gamma_matrix
    TYPE(qs_rho_type), OPTIONAL, POINTER     :: rho, rho_aux_fit, rho_buffer, &
                                                rho_xc
    TYPE(pw_p_type), OPTIONAL, POINTER       :: vppl, rho_core, rho_nlcc, &
                                                rho_nlcc_g, vee
    INTEGER, OPTIONAL                        :: neighbor_list_id
    TYPE(neighbor_list_set_p_type), DIMENSION(:), OPTIONAL, POINTER :: &
      sab_orb, sab_aux_fit, sab_aux_fit_asymm, sab_aux_fit_vs_orb, sab_all, &
      sac_ae, sac_ppl, sap_ppnl, sap_oce, sab_lrc, sab_se, sab_tbe, sab_core, &
      sab_scp, sab_vdw, sab_almo
    TYPE(kpoint_type), OPTIONAL, POINTER     :: kpoints
    LOGICAL, OPTIONAL                        :: do_kpoints
    TYPE(atomic_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: atomic_kind_set
    TYPE(qs_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: qs_kind_set
    TYPE(cell_type), OPTIONAL, POINTER       :: cell, cell_ref
    LOGICAL, OPTIONAL                        :: use_ref_cell
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particle_set
    TYPE(qs_energy_type), OPTIONAL, POINTER  :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: force
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_particles, &
                                                local_molecules
    TYPE(molecule_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: molecule_set
    TYPE(qs_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: cp_subsys
    TYPE(virial_type), OPTIONAL, POINTER     :: virial
    TYPE(cp_result_type), OPTIONAL, POINTER  :: results
    TYPE(atprop_type), OPTIONAL, POINTER     :: atprop
    INTEGER, OPTIONAL                        :: nkind, natom
    TYPE(dft_control_type), OPTIONAL, &
      POINTER                                :: dft_control
    TYPE(cp_error_type), INTENT(INOUT)       :: error

     IF(.NOT. ASSOCIATED(ks_env)) STOP "get_ks_env: not associated"
     IF(ks_env%ref_count<1) STOP "get_ks_env: ks_env%ref_count<1"
     IF(.NOT. ASSOCIATED(ks_env%subsys)) STOP "get_ks_env: subsys not associated"

     IF (PRESENT(s_mstruct_changed)) s_mstruct_changed = ks_env%s_mstruct_changed
     IF (PRESENT(rho_changed)) rho_changed = ks_env%rho_changed
     IF (PRESENT(potential_changed)) potential_changed = ks_env%potential_changed
     IF (PRESENT(forces_up_to_date)) forces_up_to_date = ks_env%forces_up_to_date
     IF (PRESENT(v_hartree_rspace)) v_hartree_rspace => ks_env%v_hartree_rspace
     IF (PRESENT(matrix_h)) matrix_h => ks_env%matrix_h
     IF (PRESENT(matrix_ks)) matrix_ks => ks_env%matrix_ks
     IF (PRESENT(matrix_vxc)) matrix_vxc => ks_env%matrix_vxc
     IF (PRESENT(matrix_ks_aux_fit)) matrix_ks_aux_fit => ks_env%matrix_ks_aux_fit
     IF (PRESENT(kinetic)) kinetic => ks_env%kinetic
     IF (PRESENT(matrix_s)) matrix_s => ks_env%matrix_s
     IF (PRESENT(matrix_s_aux_fit)) matrix_s_aux_fit => ks_env%matrix_s_aux_fit
     IF (PRESENT(matrix_s_aux_fit_vs_orb)) matrix_s_aux_fit_vs_orb => ks_env%matrix_s_aux_fit_vs_orb
     IF (PRESENT(matrix_w)) matrix_w => ks_env%matrix_w
     IF (PRESENT(matrix_w_mp2)) matrix_w_mp2 => ks_env%matrix_w_mp2
     IF (PRESENT(matrix_p_mp2)) matrix_p_mp2 => ks_env%matrix_p_mp2
     IF (PRESENT(gamma_matrix)) gamma_matrix => ks_env%gamma_matrix
     IF (PRESENT(rho)) rho => ks_env%rho
     IF (PRESENT(rho_aux_fit)) rho_aux_fit => ks_env%rho_aux_fit
     IF (PRESENT(rho_buffer)) rho_buffer => ks_env%rho_buffer
     IF (PRESENT(rho_xc)) rho_xc => ks_env%rho_xc
     IF (PRESENT(rho_core)) rho_core => ks_env%rho_core
     IF (PRESENT(rho_nlcc)) rho_nlcc => ks_env%rho_nlcc
     IF (PRESENT(rho_nlcc_g)) rho_nlcc_g => ks_env%rho_nlcc_g
     IF (PRESENT(vppl)) vppl => ks_env%vppl
     IF (PRESENT(vee)) vee => ks_env%vee

     IF (PRESENT(neighbor_list_id)) neighbor_list_id=ks_env%neighbor_list_id
     IF (PRESENT(sab_orb)) sab_orb => ks_env%sab_orb
     IF (PRESENT(sab_aux_fit)) sab_aux_fit => ks_env%sab_aux_fit
     IF (PRESENT(sab_aux_fit_asymm)) sab_aux_fit_asymm => ks_env%sab_aux_fit_asymm
     IF (PRESENT(sab_aux_fit_vs_orb)) sab_aux_fit_vs_orb => ks_env%sab_aux_fit_vs_orb
     IF (PRESENT(sab_all)) sab_all => ks_env%sab_all
     IF (PRESENT(sab_scp)) sab_scp => ks_env%sab_scp
     IF (PRESENT(sab_vdw)) sab_vdw => ks_env%sab_vdw
     IF (PRESENT(sac_ae)) sac_ae => ks_env%sac_ae
     IF (PRESENT(sac_ppl)) sac_ppl => ks_env%sac_ppl
     IF (PRESENT(sap_ppnl)) sap_ppnl => ks_env%sap_ppnl
     IF (PRESENT(sap_oce)) sap_oce =>  ks_env%sap_oce
     IF (PRESENT(sab_se)) sab_se =>  ks_env%sab_se
     IF (PRESENT(sab_lrc)) sab_lrc =>  ks_env%sab_lrc
     IF (PRESENT(sab_tbe)) sab_tbe =>  ks_env%sab_tbe
     IF (PRESENT(sab_core)) sab_core =>  ks_env%sab_core
     IF (PRESENT(sab_almo)) sab_almo =>  ks_env%sab_almo
     IF (PRESENT(dft_control)) dft_control => ks_env%dft_control

     IF (PRESENT(kpoints)) kpoints => ks_env%kpoints
     IF (PRESENT(do_kpoints)) THEN
       do_kpoints = (ks_env%kpoints%nkp > 0)
     END IF

     IF (PRESENT(subsys)) subsys => ks_env%subsys
     CALL qs_subsys_get(ks_env%subsys,&
                        local_molecules_new=local_molecules,&
                        local_particles=local_particles,&
                        atprop=atprop,&
                        virial=virial,&
                        results=results,&
                        cell=cell,&
                        cell_ref=cell_ref,&
                        use_ref_cell=use_ref_cell,&
                        energy=energy,&
                        force=force,&
                        qs_kind_set=qs_kind_set,&
                        cp_subsys=cp_subsys,&
                        atomic_kind_set=atomic_kind_set,&
                        particle_set=particle_set,&
                        molecule_kind_set=molecule_kind_set,&
                        molecule_set=molecule_set,&
                        natom=natom,&
                        nkind=nkind,&
                        error=error)

  END SUBROUTINE get_ks_env


! *****************************************************************************
!> \brief ...
!> \param ks_env ...
!> \param v_hartree_rspace ...
!> \param s_mstruct_changed ...
!> \param rho_changed ...
!> \param potential_changed ...
!> \param forces_up_to_date ...
!> \param matrix_h ...
!> \param matrix_ks ...
!> \param matrix_vxc ...
!> \param matrix_ks_aux_fit ...
!> \param kinetic ...
!> \param matrix_s ...
!> \param matrix_s_aux_fit ...
!> \param matrix_s_aux_fit_vs_orb ...
!> \param matrix_w ...
!> \param matrix_w_mp2 ...
!> \param matrix_p_mp2 ...
!> \param gamma_matrix ...
!> \param vppl ...
!> \param rho_core ...
!> \param rho_nlcc ...
!> \param rho_nlcc_g ...
!> \param vee ...
!> \param neighbor_list_id ...
!> \param sab_orb ...
!> \param sab_aux_fit ...
!> \param sab_aux_fit_asymm ...
!> \param sab_aux_fit_vs_orb ...
!> \param sab_all ...
!> \param sac_ae ...
!> \param sac_ppl ...
!> \param sap_ppnl ...
!> \param sap_oce ...
!> \param sab_lrc ...
!> \param sab_se ...
!> \param sab_tbe ...
!> \param sab_core ...
!> \param sab_scp ...
!> \param sab_vdw ...
!> \param sab_almo ...
!> \param subsys ...
!> \param dft_control ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE set_ks_env(ks_env,v_hartree_rspace,&
                        s_mstruct_changed,rho_changed,&
                        potential_changed,forces_up_to_date,&
                        matrix_h,matrix_ks,matrix_vxc,matrix_ks_aux_fit,&
                        kinetic,matrix_s,matrix_s_aux_fit,matrix_s_aux_fit_vs_orb,&
                        matrix_w,matrix_w_mp2,matrix_p_mp2,gamma_matrix,&
                        vppl, rho_core, rho_nlcc, rho_nlcc_g, vee,&
                        neighbor_list_id,&
                        sab_orb,sab_aux_fit,sab_aux_fit_asymm, &
                        sab_aux_fit_vs_orb,sab_all,sac_ae,sac_ppl,sap_ppnl,sap_oce,sab_lrc,&
                        sab_se,sab_tbe,sab_core,sab_scp,sab_vdw,sab_almo,&
                        subsys,dft_control,&
                        error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(pw_type), OPTIONAL, POINTER         :: v_hartree_rspace
    LOGICAL, OPTIONAL                        :: s_mstruct_changed, &
                                                rho_changed, &
                                                potential_changed, &
                                                forces_up_to_date
    TYPE(cp_dbcsr_p_type), DIMENSION(:), OPTIONAL, POINTER :: matrix_h, &
      matrix_ks, matrix_vxc, matrix_ks_aux_fit, kinetic, matrix_s, &
      matrix_s_aux_fit, matrix_s_aux_fit_vs_orb, matrix_w, matrix_w_mp2, &
      matrix_p_mp2, gamma_matrix
    TYPE(pw_p_type), OPTIONAL, POINTER       :: vppl, rho_core, rho_nlcc, &
                                                rho_nlcc_g, vee
    INTEGER, OPTIONAL                        :: neighbor_list_id
    TYPE(neighbor_list_set_p_type), DIMENSION(:), OPTIONAL, POINTER :: &
      sab_orb, sab_aux_fit, sab_aux_fit_asymm, sab_aux_fit_vs_orb, sab_all, &
      sac_ae, sac_ppl, sap_ppnl, sap_oce, sab_lrc, sab_se, sab_tbe, sab_core, &
      sab_scp, sab_vdw, sab_almo
    TYPE(qs_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(dft_control_type), OPTIONAL, &
      POINTER                                :: dft_control
    TYPE(cp_error_type), INTENT(INOUT)       :: error

     IF(.NOT. ASSOCIATED(ks_env)) STOP "set_ks_env: not associated"
     IF(ks_env%ref_count<1) STOP "set_ks_env: ks_env%ref_count<1"

     IF (PRESENT(s_mstruct_changed)) ks_env%s_mstruct_changed=s_mstruct_changed
     IF (PRESENT(rho_changed)) ks_env%rho_changed=rho_changed
     IF (PRESENT(potential_changed)) ks_env%potential_changed=potential_changed
     IF (PRESENT(forces_up_to_date)) ks_env%forces_up_to_date=forces_up_to_date
     IF (ks_env%s_mstruct_changed .OR. ks_env%potential_changed .OR. ks_env%rho_changed)&
          ks_env%forces_up_to_date=.FALSE.
     IF (PRESENT(v_hartree_rspace)) ks_env%v_hartree_rspace => v_hartree_rspace
     IF (PRESENT(matrix_h)) ks_env%matrix_h => matrix_h
     IF (PRESENT(matrix_ks)) ks_env%matrix_ks => matrix_ks
     IF (PRESENT(matrix_vxc)) ks_env%matrix_vxc => matrix_vxc
     IF (PRESENT(matrix_ks_aux_fit)) ks_env%matrix_ks_aux_fit => matrix_ks_aux_fit
     IF (PRESENT(matrix_s)) ks_env%matrix_s => matrix_s
     IF (PRESENT(matrix_s_aux_fit)) ks_env%matrix_s_aux_fit => matrix_s_aux_fit
     IF (PRESENT(matrix_s_aux_fit_vs_orb)) ks_env%matrix_s_aux_fit_vs_orb => matrix_s_aux_fit_vs_orb
     IF (PRESENT(matrix_w)) ks_env%matrix_w => matrix_w
     IF (PRESENT(matrix_w_mp2)) ks_env%matrix_w_mp2 => matrix_w_mp2
     IF (PRESENT(matrix_p_mp2)) ks_env%matrix_p_mp2 => matrix_p_mp2
     IF (PRESENT(gamma_matrix)) ks_env%gamma_matrix => gamma_matrix
     IF (PRESENT(kinetic)) ks_env%kinetic => kinetic
     IF (PRESENT(rho_core)) ks_env%rho_core => rho_core
     IF (PRESENT(rho_nlcc)) ks_env%rho_nlcc => rho_nlcc
     IF (PRESENT(rho_nlcc_g)) ks_env%rho_nlcc_g => rho_nlcc_g
     IF (PRESENT(vppl)) ks_env%vppl => vppl
     IF (PRESENT(vee)) ks_env%vee => vee

     IF (PRESENT(neighbor_list_id)) ks_env%neighbor_list_id=neighbor_list_id
     IF (PRESENT(sab_orb)) ks_env%sab_orb => sab_orb
     IF (PRESENT(sab_aux_fit)) ks_env%sab_aux_fit => sab_aux_fit
     IF (PRESENT(sab_aux_fit_asymm)) ks_env%sab_aux_fit_asymm => sab_aux_fit_asymm
     IF (PRESENT(sab_aux_fit_vs_orb)) ks_env%sab_aux_fit_vs_orb => sab_aux_fit_vs_orb
     IF (PRESENT(sab_scp)) ks_env%sab_scp => sab_scp
     IF (PRESENT(sab_vdw)) ks_env%sab_vdw => sab_vdw
     IF (PRESENT(sab_all)) ks_env%sab_all => sab_all
     IF (PRESENT(sac_ae)) ks_env%sac_ae => sac_ae
     IF (PRESENT(sac_ppl)) ks_env%sac_ppl => sac_ppl
     IF (PRESENT(sap_ppnl)) ks_env%sap_ppnl => sap_ppnl
     IF (PRESENT(sap_oce)) ks_env%sap_oce => sap_oce
     IF (PRESENT(sab_se)) ks_env%sab_se => sab_se
     IF (PRESENT(sab_lrc)) ks_env%sab_lrc => sab_lrc
     IF (PRESENT(sab_tbe)) ks_env%sab_tbe => sab_tbe
     IF (PRESENT(sab_core)) ks_env%sab_core => sab_core
     IF (PRESENT(sab_almo)) ks_env%sab_almo => sab_almo

     IF (PRESENT(subsys)) THEN
       CALL qs_subsys_release(ks_env%subsys,error)
       CALL qs_subsys_retain(subsys, error=error)
       ks_env%subsys => subsys
    END IF
    IF (PRESENT(dft_control)) THEN
       CALL dft_control_retain(dft_control, error=error)
       CALL dft_control_release(ks_env%dft_control, error=error)
       ks_env%dft_control => dft_control
    END IF
  END SUBROUTINE set_ks_env

! *****************************************************************************
!> \brief releases the ks_env (see doc/ReferenceCounting.html)
!> \param ks_env the ks_env to be released
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE qs_ks_release(ks_env,error)
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_release', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(ks_env)) THEN
       CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)
       ks_env%ref_count=ks_env%ref_count-1

       IF (ks_env%ref_count<1) THEN
          IF (ASSOCIATED(ks_env%v_hartree_rspace)) &
            CALL pw_release(ks_env%v_hartree_rspace,error=error)
          IF (ASSOCIATED(ks_env%matrix_h)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_h,error)
          IF (ASSOCIATED(ks_env%matrix_ks)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_ks,error)
          IF (ASSOCIATED(ks_env%matrix_vxc)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_vxc,error)
          IF (ASSOCIATED(ks_env%matrix_ks_aux_fit)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_ks_aux_fit,error)
          IF (ASSOCIATED(ks_env%matrix_s)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_s,error)
          IF (ASSOCIATED(ks_env%matrix_s_aux_fit)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_s_aux_fit,error)
          IF (ASSOCIATED(ks_env%matrix_s_aux_fit_vs_orb)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_s_aux_fit_vs_orb,error)
          IF (ASSOCIATED(ks_env%matrix_w)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_w,error)
          IF (ASSOCIATED(ks_env%matrix_w_mp2)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_w_mp2,error)
          IF (ASSOCIATED(ks_env%matrix_p_mp2)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_p_mp2,error)
          IF (ASSOCIATED(ks_env%gamma_matrix)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%gamma_matrix,error)
          IF (ASSOCIATED(ks_env%kinetic))&
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%kinetic,error)
          IF(ASSOCIATED(ks_env%rho))&
            CALL qs_rho_release(ks_env%rho, error=error)
          IF(ASSOCIATED(ks_env%rho_buffer))&
            CALL qs_rho_release(ks_env%rho_buffer,error=error)
          IF(ASSOCIATED(ks_env%rho_xc))&
            CALL qs_rho_release(ks_env%rho_xc, error=error)
          IF(ASSOCIATED(ks_env%rho_aux_fit))&
            CALL qs_rho_release(ks_env%rho_aux_fit, error=error)

          IF (ASSOCIATED(ks_env%rho_nlcc_g)) THEN
             CALL pw_release(ks_env%rho_nlcc_g%pw,error=error)
             DEALLOCATE(ks_env%rho_nlcc_g)
          ENDIF
          IF (ASSOCIATED(ks_env%rho_nlcc)) THEN
             CALL pw_release(ks_env%rho_nlcc%pw,error=error)
             DEALLOCATE(ks_env%rho_nlcc)
          ENDIF
          IF (ASSOCIATED(ks_env%rho_core)) THEN
             CALL pw_release(ks_env%rho_core%pw,error=error)
             DEALLOCATE(ks_env%rho_core)
          ENDIF
          IF (ASSOCIATED(ks_env%vppl)) THEN
             CALL pw_release(ks_env%vppl%pw,error=error)
             DEALLOCATE(ks_env%vppl)
          ENDIF
          IF (ASSOCIATED(ks_env%vee)) THEN
             CALL pw_release(ks_env%vee%pw,error=error)
             DEALLOCATE(ks_env%vee)
          ENDIF

          CALL release_sab(ks_env%sab_orb)
          CALL release_sab(ks_env%sac_ae)
          CALL release_sab(ks_env%sac_ppl)
          CALL release_sab(ks_env%sap_ppnl)
          CALL release_sab(ks_env%sap_oce)
          CALL release_sab(ks_env%sab_se)
          CALL release_sab(ks_env%sab_vdw)
          CALL release_sab(ks_env%sab_tbe)
          CALL release_sab(ks_env%sab_core)
          CALL release_sab(ks_env%sab_all)
          CALL release_sab(ks_env%sab_scp)
          CALL release_sab(ks_env%sab_aux_fit)
          CALL release_sab(ks_env%sab_aux_fit_vs_orb)
          CALL release_sab(ks_env%sab_aux_fit_asymm)
          CALL release_sab(ks_env%sab_lrc)
          CALL release_sab(ks_env%sab_almo)
          CALL dft_control_release(ks_env%dft_control, error=error)

          CALL kpoint_release(ks_env%kpoints,error=error)
          CALL qs_subsys_release(ks_env%subsys,error=error)
          DEALLOCATE(ks_env)
       END IF
    END IF
    NULLIFY(ks_env)
  END SUBROUTINE qs_ks_release


! *****************************************************************************
!> \brief releases an array of neighbor_list_sets
!> \param sab ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE release_sab(sab)
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab

    INTEGER                                  :: i

     IF (ASSOCIATED(sab)) THEN
        DO i=1,SIZE(sab)
           CALL deallocate_neighbor_list_set(sab(i)%neighbor_list_set)
        END DO
        DEALLOCATE(sab)
     END IF
  END SUBROUTINE release_sab


! *****************************************************************************
!> \brief retains the given ks_environment
!> \param ks_env the KohnSham environment to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE qs_ks_retain(ks_env, error)
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(ks_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(ks_env%ref_count>0,cp_failure_level,routineP,error,failure)
     ks_env%ref_count=ks_env%ref_count+1
  END IF
END SUBROUTINE qs_ks_retain

END MODULE qs_ks_types
