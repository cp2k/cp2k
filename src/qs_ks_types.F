!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      11.2002 created [fawzi]
!>      10.2014 Moved many components from qs_env to here [Ole Schuett]
!>      11.2014 unified k-point and gamma-point code [Ole Schuett]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE qs_ks_types
  USE admm_dm_types,                   ONLY: admm_dm_release,&
                                             admm_dm_type
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE atprop_types,                    ONLY: atprop_type
  USE cell_types,                      ONLY: cell_type
  USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                             cp_blacs_env_retain,&
                                             cp_blacs_env_type
  USE cp_control_types,                ONLY: dft_control_release,&
                                             dft_control_retain,&
                                             dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_deallocate_matrix_set,&
                                             cp_dbcsr_distribution_release,&
                                             cp_dbcsr_p_type,&
                                             dbcsr_distribution_obj
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_result_types,                 ONLY: cp_result_type
  USE cp_subsys_types,                 ONLY: cp_subsys_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                             distribution_2d_retain,&
                                             distribution_2d_type
  USE kpoint_transitional,             ONLY: get_1d_pointer,&
                                             get_2d_pointer,&
                                             kpoint_transitional_release,&
                                             kpoint_transitional_type,&
                                             set_1d_pointer,&
                                             set_2d_pointer
  USE kpoint_types,                    ONLY: kpoint_create,&
                                             kpoint_release,&
                                             kpoint_type
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_release,&
                                             pw_env_retain,&
                                             pw_env_type
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
  USE task_list_types,                 ONLY: deallocate_task_list,&
                                             task_list_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ks_types'

  PUBLIC :: qs_ks_env_type, qs_ks_env_create, qs_ks_did_change
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
!> \param sab_scp: neighbor lists for the calculation of self-consistent polarization
!> \param sab_almo: neighbor lists to create ALMO delocalization template
!> \param sab_kp: neighbor lists to create kp image cell lists
!>
!> \param kpoints information on the kpoints used
!> \param subsys the particles, molecules,... of this environment
!> \param dft_control object that contains the values of various parameters
!>        relevant to a dft calculation
!> \param distribution_2d: distribution of the atom pairs between the processors
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

    TYPE(kpoint_transitional_type)                        :: matrix_ks,&
                                                             matrix_s,&
                                                             kinetic,&
                                                             matrix_h,&
                                                             matrix_w,&
                                                             matrix_vxc

    TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks_im => Null(),&
                                                             matrix_w_mp2 => Null(),&
                                                             matrix_p_mp2 => Null(),&
                                                             gamma_matrix => Null(),&
                                                             matrix_s_aux_fit => Null(),&
                                                             matrix_s_aux_fit_vs_orb => Null(),&
                                                             matrix_ks_aux_fit => Null(),&
                                                             matrix_ks_aux_fit_im => Null(),&
                                                             matrix_ks_aux_fit_dft => Null(),&
                                                             matrix_ks_aux_fit_hfx => Null()

    TYPE(qs_rho_type), POINTER                            :: rho => Null(),&
                                                             rho_buffer => Null(),&
                                                             rho_xc => Null(),&
                                                             rho_aux_fit => Null(),&
                                                             rho_aux_fit_buffer => Null()

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
                                                             sab_scp => Null(),&
                                                             sab_tbe => Null(),&
                                                             sab_core => Null(),&
                                                             sab_all => Null(),&
                                                             sab_aux_fit => Null(),&
                                                             sab_aux_fit_vs_orb => Null(),&
                                                             sab_aux_fit_asymm => Null(),&
                                                             sab_lrc => Null(),&
                                                             sab_almo => Null(),&
                                                             sab_kp => Null()

    TYPE(task_list_type), POINTER                         :: task_list => Null()
    TYPE(task_list_type), POINTER                         :: task_list_aux_fit => Null()
    TYPE(task_list_type), POINTER                         :: task_list_soft => Null()

    TYPE(kpoint_type), POINTER                            :: kpoints => Null()
    TYPE(qs_subsys_type), POINTER                         :: subsys => Null()
    TYPE(dft_control_type), POINTER                       :: dft_control => Null()
    TYPE(dbcsr_distribution_obj), POINTER                 :: dbcsr_dist => Null()
    TYPE(distribution_2d_type), POINTER                   :: distribution_2d => Null()
    TYPE(pw_env_type), POINTER                            :: pw_env => Null()
    TYPE(cp_para_env_type), POINTER                       :: para_env => Null()
    TYPE(cp_blacs_env_type), POINTER                      :: blacs_env => Null()
    TYPE(admm_dm_type), POINTER                           :: admm_dm => Null()
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
!> \author Ole Schuett
! *****************************************************************************
SUBROUTINE qs_ks_env_create(ks_env)
    TYPE(qs_ks_env_type), POINTER            :: ks_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_ks_env_create', &
      routineP = moduleN//':'//routineN

    IF(ASSOCIATED(ks_env)) STOP routineP//"ks_env already associated"
    ALLOCATE(ks_env)
    CALL qs_rho_create(ks_env%rho)
    CALL qs_rho_create(ks_env%rho_aux_fit)
    CALL qs_rho_create(ks_env%rho_aux_fit_buffer)
    CALL qs_rho_create(ks_env%rho_xc)
    CALL qs_rho_create(ks_env%rho_buffer)
    CALL kpoint_create(ks_env%kpoints)
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
!> \param matrix_ks_im ...
!> \param matrix_vxc ...
!> \param matrix_ks_aux_fit ...
!> \param matrix_ks_aux_fit_im ...
!> \param matrix_ks_aux_fit_dft ...
!> \param matrix_ks_aux_fit_hfx ...
!> \param kinetic ...
!> \param matrix_s ...
!> \param matrix_s_aux_fit ...
!> \param matrix_s_aux_fit_vs_orb ...
!> \param matrix_w ...
!> \param matrix_w_mp2 ...
!> \param matrix_p_mp2 ...
!> \param gamma_matrix ...
!> \param matrix_h_kp ...
!> \param matrix_ks_kp ...
!> \param matrix_vxc_kp ...
!> \param kinetic_kp ...
!> \param matrix_s_kp ...
!> \param matrix_w_kp ...
!> \param rho ...
!> \param rho_aux_fit ...
!> \param rho_aux_fit_buffer ...
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
!> \param sab_vdw ...
!> \param sab_scp ...
!> \param sab_almo ...
!> \param sab_kp ...
!> \param task_list ...
!> \param task_list_aux_fit ...
!> \param task_list_soft ...
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
!> \param dbcsr_dist ...
!> \param distribution_2d ...
!> \param pw_env ...
!> \param para_env ...
!> \param blacs_env ...
!> \param nelectron_total ...
!> \param nelectron_spin ...
!> \param admm_dm ...
! *****************************************************************************
  SUBROUTINE get_ks_env(ks_env,v_hartree_rspace,&
                        s_mstruct_changed,rho_changed,&
                        potential_changed,forces_up_to_date,&
                        matrix_h,matrix_ks,matrix_ks_im,matrix_vxc,matrix_ks_aux_fit,&
                        matrix_ks_aux_fit_im,matrix_ks_aux_fit_dft,matrix_ks_aux_fit_hfx,&
                        kinetic,matrix_s,matrix_s_aux_fit,matrix_s_aux_fit_vs_orb,&
                        matrix_w,matrix_w_mp2,matrix_p_mp2,gamma_matrix,&
                        matrix_h_kp, matrix_ks_kp, matrix_vxc_kp, kinetic_kp, matrix_s_kp, matrix_w_kp,&
                        rho, rho_aux_fit, rho_aux_fit_buffer, &
                        rho_buffer, rho_xc,&
                        vppl, rho_core, rho_nlcc, rho_nlcc_g, vee,&
                        neighbor_list_id,&
                        sab_orb,sab_aux_fit,sab_aux_fit_asymm, &
                        sab_aux_fit_vs_orb,sab_all,sac_ae,sac_ppl,sap_ppnl,sap_oce,sab_lrc,&
                        sab_se,sab_tbe,sab_core,sab_vdw,sab_scp,sab_almo,sab_kp,&
                        task_list, task_list_aux_fit, task_list_soft,&
                        kpoints, do_kpoints,&
                        atomic_kind_set,qs_kind_set,cell,cell_ref,use_ref_cell,&
                        particle_set,energy,force,local_particles,local_molecules,&
                        molecule_kind_set,molecule_set,subsys,cp_subsys,virial,results,atprop,&
                        nkind, natom, dft_control,dbcsr_dist,distribution_2d,pw_env,&
                        para_env,blacs_env, nelectron_total, nelectron_spin,admm_dm)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(pw_type), OPTIONAL, POINTER         :: v_hartree_rspace
    LOGICAL, OPTIONAL                        :: s_mstruct_changed, &
                                                rho_changed, &
                                                potential_changed, &
                                                forces_up_to_date
    TYPE(cp_dbcsr_p_type), DIMENSION(:), OPTIONAL, POINTER :: matrix_h, &
      matrix_ks, matrix_ks_im, matrix_vxc, matrix_ks_aux_fit, &
      matrix_ks_aux_fit_im, matrix_ks_aux_fit_dft, matrix_ks_aux_fit_hfx, &
      kinetic, matrix_s, matrix_s_aux_fit, matrix_s_aux_fit_vs_orb, matrix_w, &
      matrix_w_mp2, matrix_p_mp2, gamma_matrix
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: matrix_h_kp, matrix_ks_kp, &
                                                matrix_vxc_kp, kinetic_kp, &
                                                matrix_s_kp, matrix_w_kp
    TYPE(qs_rho_type), OPTIONAL, POINTER     :: rho, rho_aux_fit, &
                                                rho_aux_fit_buffer, &
                                                rho_buffer, rho_xc
    TYPE(pw_p_type), OPTIONAL, POINTER       :: vppl, rho_core, rho_nlcc, &
                                                rho_nlcc_g, vee
    INTEGER, OPTIONAL                        :: neighbor_list_id
    TYPE(neighbor_list_set_p_type), DIMENSION(:), OPTIONAL, POINTER :: &
      sab_orb, sab_aux_fit, sab_aux_fit_asymm, sab_aux_fit_vs_orb, sab_all, &
      sac_ae, sac_ppl, sap_ppnl, sap_oce, sab_lrc, sab_se, sab_tbe, sab_core, &
      sab_vdw, sab_scp, sab_almo, sab_kp
    TYPE(task_list_type), OPTIONAL, POINTER  :: task_list, task_list_aux_fit, &
                                                task_list_soft
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
    TYPE(dbcsr_distribution_obj), OPTIONAL, &
      POINTER                                :: dbcsr_dist
    TYPE(distribution_2d_type), OPTIONAL, &
      POINTER                                :: distribution_2d
    TYPE(pw_env_type), OPTIONAL, POINTER     :: pw_env
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(cp_blacs_env_type), OPTIONAL, &
      POINTER                                :: blacs_env
    INTEGER, OPTIONAL                        :: nelectron_total
    INTEGER, DIMENSION(2), OPTIONAL          :: nelectron_spin
    TYPE(admm_dm_type), OPTIONAL, POINTER    :: admm_dm

     IF(.NOT. ASSOCIATED(ks_env)) STOP "get_ks_env: not associated"
     IF(ks_env%ref_count<1) STOP "get_ks_env: ks_env%ref_count<1"
     IF(.NOT. ASSOCIATED(ks_env%subsys)) STOP "get_ks_env: subsys not associated"

     IF (PRESENT(s_mstruct_changed)) s_mstruct_changed = ks_env%s_mstruct_changed
     IF (PRESENT(rho_changed)) rho_changed = ks_env%rho_changed
     IF (PRESENT(potential_changed)) potential_changed = ks_env%potential_changed
     IF (PRESENT(forces_up_to_date)) forces_up_to_date = ks_env%forces_up_to_date
     IF (PRESENT(v_hartree_rspace)) v_hartree_rspace => ks_env%v_hartree_rspace

     IF (PRESENT(matrix_ks)) matrix_ks         => get_1d_pointer(ks_env%matrix_ks)
     IF (PRESENT(matrix_s)) matrix_s           => get_1d_pointer(ks_env%matrix_s)
     IF (PRESENT(kinetic)) kinetic             => get_1d_pointer(ks_env%kinetic)
     IF (PRESENT(matrix_h)) matrix_h           => get_1d_pointer(ks_env%matrix_h)
     IF (PRESENT(matrix_w)) matrix_w           => get_1d_pointer(ks_env%matrix_w)
     IF (PRESENT(matrix_vxc)) matrix_vxc       => get_1d_pointer(ks_env%matrix_vxc)

     IF (PRESENT(matrix_ks_kp)) matrix_ks_kp   => get_2d_pointer(ks_env%matrix_ks)
     IF (PRESENT(matrix_s_kp)) matrix_s_kp     => get_2d_pointer(ks_env%matrix_s)
     IF (PRESENT(matrix_w_kp)) matrix_w_kp     => get_2d_pointer(ks_env%matrix_w)
     IF (PRESENT(kinetic_kp)) kinetic_kp       => get_2d_pointer(ks_env%kinetic)
     IF (PRESENT(matrix_h_kp)) matrix_h_kp     => get_2d_pointer(ks_env%matrix_h)
     IF (PRESENT(matrix_vxc_kp)) matrix_vxc_kp => get_2d_pointer(ks_env%matrix_vxc)

     IF (PRESENT(matrix_ks_im)) matrix_ks_im => ks_env%matrix_ks_im
     IF (PRESENT(matrix_ks_aux_fit)) matrix_ks_aux_fit => ks_env%matrix_ks_aux_fit
     IF (PRESENT(matrix_ks_aux_fit_im)) matrix_ks_aux_fit_im => ks_env%matrix_ks_aux_fit_im
     IF (PRESENT(matrix_ks_aux_fit_dft)) matrix_ks_aux_fit_dft => ks_env%matrix_ks_aux_fit_dft
     IF (PRESENT(matrix_ks_aux_fit_hfx)) matrix_ks_aux_fit_hfx => ks_env%matrix_ks_aux_fit_hfx
     IF (PRESENT(matrix_s_aux_fit)) matrix_s_aux_fit => ks_env%matrix_s_aux_fit
     IF (PRESENT(matrix_s_aux_fit_vs_orb)) matrix_s_aux_fit_vs_orb => ks_env%matrix_s_aux_fit_vs_orb
     IF (PRESENT(matrix_w_mp2)) matrix_w_mp2 => ks_env%matrix_w_mp2
     IF (PRESENT(matrix_p_mp2)) matrix_p_mp2 => ks_env%matrix_p_mp2
     IF (PRESENT(gamma_matrix)) gamma_matrix => ks_env%gamma_matrix
     IF (PRESENT(rho)) rho => ks_env%rho
     IF (PRESENT(rho_aux_fit)) rho_aux_fit => ks_env%rho_aux_fit
     IF (PRESENT(rho_aux_fit_buffer)) rho_aux_fit_buffer => ks_env%rho_aux_fit_buffer
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
     IF (PRESENT(sab_vdw)) sab_vdw => ks_env%sab_vdw
     IF (PRESENT(sab_scp)) sab_scp => ks_env%sab_scp
     IF (PRESENT(sac_ae)) sac_ae => ks_env%sac_ae
     IF (PRESENT(sac_ppl)) sac_ppl => ks_env%sac_ppl
     IF (PRESENT(sap_ppnl)) sap_ppnl => ks_env%sap_ppnl
     IF (PRESENT(sap_oce)) sap_oce =>  ks_env%sap_oce
     IF (PRESENT(sab_se)) sab_se =>  ks_env%sab_se
     IF (PRESENT(sab_lrc)) sab_lrc =>  ks_env%sab_lrc
     IF (PRESENT(sab_tbe)) sab_tbe =>  ks_env%sab_tbe
     IF (PRESENT(sab_core)) sab_core =>  ks_env%sab_core
     IF (PRESENT(sab_almo)) sab_almo =>  ks_env%sab_almo
     IF (PRESENT(sab_kp)) sab_kp =>  ks_env%sab_kp
     IF (PRESENT(dft_control)) dft_control => ks_env%dft_control
     IF (PRESENT(dbcsr_dist)) dbcsr_dist => ks_env%dbcsr_dist
     IF (PRESENT(distribution_2d)) distribution_2d => ks_env%distribution_2d
     IF (PRESENT(pw_env)) pw_env => ks_env%pw_env
     IF (PRESENT(para_env)) para_env => ks_env%para_env
     IF (PRESENT(blacs_env)) blacs_env => ks_env%blacs_env
     IF (PRESENT(admm_dm)) admm_dm => ks_env%admm_dm

     IF (PRESENT(task_list)) task_list => ks_env%task_list
     IF (PRESENT(task_list_aux_fit)) task_list_aux_fit => ks_env%task_list_aux_fit
     IF (PRESENT(task_list_soft)) task_list_soft => ks_env%task_list_soft

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
                        nelectron_total=nelectron_total,&
                        nelectron_spin=nelectron_spin)

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
!> \param matrix_ks_im ...
!> \param matrix_vxc ...
!> \param matrix_ks_aux_fit ...
!> \param matrix_ks_aux_fit_im ...
!> \param matrix_ks_aux_fit_dft ...
!> \param matrix_ks_aux_fit_hfx ...
!> \param kinetic ...
!> \param matrix_s ...
!> \param matrix_s_aux_fit ...
!> \param matrix_s_aux_fit_vs_orb ...
!> \param matrix_w ...
!> \param matrix_w_mp2 ...
!> \param matrix_p_mp2 ...
!> \param gamma_matrix ...
!> \param matrix_h_kp ...
!> \param matrix_ks_kp ...
!> \param matrix_vxc_kp ...
!> \param kinetic_kp ...
!> \param matrix_s_kp ...
!> \param matrix_w_kp ...
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
!> \param sab_vdw ...
!> \param sab_scp ...
!> \param sab_almo ...
!> \param sab_kp ...
!> \param task_list ...
!> \param task_list_aux_fit ...
!> \param task_list_soft ...
!> \param subsys ...
!> \param dft_control ...
!> \param dbcsr_dist ...
!> \param distribution_2d ...
!> \param pw_env ...
!> \param para_env ...
!> \param blacs_env ...
!> \param admm_dm ...
! *****************************************************************************
  SUBROUTINE set_ks_env(ks_env,v_hartree_rspace,&
                        s_mstruct_changed,rho_changed,&
                        potential_changed,forces_up_to_date,&
                        matrix_h,matrix_ks,matrix_ks_im,matrix_vxc,matrix_ks_aux_fit,matrix_ks_aux_fit_im,&
                        matrix_ks_aux_fit_dft,matrix_ks_aux_fit_hfx,kinetic,matrix_s,matrix_s_aux_fit,&
                        matrix_s_aux_fit_vs_orb, matrix_w,matrix_w_mp2,matrix_p_mp2,gamma_matrix,&
                        matrix_h_kp, matrix_ks_kp, matrix_vxc_kp, kinetic_kp, matrix_s_kp, matrix_w_kp,&
                        vppl, rho_core, rho_nlcc, rho_nlcc_g, vee,&
                        neighbor_list_id,&
                        sab_orb,sab_aux_fit,sab_aux_fit_asymm, &
                        sab_aux_fit_vs_orb,sab_all,sac_ae,sac_ppl,sap_ppnl,sap_oce,sab_lrc,&
                        sab_se,sab_tbe,sab_core,sab_vdw,sab_scp,sab_almo,sab_kp,&
                        task_list,task_list_aux_fit,task_list_soft,&
                        subsys,dft_control,dbcsr_dist,distribution_2d,pw_env,&
                        para_env,blacs_env,admm_dm)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(pw_type), OPTIONAL, POINTER         :: v_hartree_rspace
    LOGICAL, OPTIONAL                        :: s_mstruct_changed, &
                                                rho_changed, &
                                                potential_changed, &
                                                forces_up_to_date
    TYPE(cp_dbcsr_p_type), DIMENSION(:), OPTIONAL, POINTER :: matrix_h, &
      matrix_ks, matrix_ks_im, matrix_vxc, matrix_ks_aux_fit, &
      matrix_ks_aux_fit_im, matrix_ks_aux_fit_dft, matrix_ks_aux_fit_hfx, &
      kinetic, matrix_s, matrix_s_aux_fit, matrix_s_aux_fit_vs_orb, matrix_w, &
      matrix_w_mp2, matrix_p_mp2, gamma_matrix
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: matrix_h_kp, matrix_ks_kp, &
                                                matrix_vxc_kp, kinetic_kp, &
                                                matrix_s_kp, matrix_w_kp
    TYPE(pw_p_type), OPTIONAL, POINTER       :: vppl, rho_core, rho_nlcc, &
                                                rho_nlcc_g, vee
    INTEGER, OPTIONAL                        :: neighbor_list_id
    TYPE(neighbor_list_set_p_type), DIMENSION(:), OPTIONAL, POINTER :: &
      sab_orb, sab_aux_fit, sab_aux_fit_asymm, sab_aux_fit_vs_orb, sab_all, &
      sac_ae, sac_ppl, sap_ppnl, sap_oce, sab_lrc, sab_se, sab_tbe, sab_core, &
      sab_vdw, sab_scp, sab_almo, sab_kp
    TYPE(task_list_type), OPTIONAL, POINTER  :: task_list, task_list_aux_fit, &
                                                task_list_soft
    TYPE(qs_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(dft_control_type), OPTIONAL, &
      POINTER                                :: dft_control
    TYPE(dbcsr_distribution_obj), OPTIONAL, &
      POINTER                                :: dbcsr_dist
    TYPE(distribution_2d_type), OPTIONAL, &
      POINTER                                :: distribution_2d
    TYPE(pw_env_type), OPTIONAL, POINTER     :: pw_env
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(cp_blacs_env_type), OPTIONAL, &
      POINTER                                :: blacs_env
    TYPE(admm_dm_type), OPTIONAL, POINTER    :: admm_dm

     IF(.NOT. ASSOCIATED(ks_env)) STOP "set_ks_env: not associated"
     IF(ks_env%ref_count<1) STOP "set_ks_env: ks_env%ref_count<1"

     IF (PRESENT(s_mstruct_changed)) ks_env%s_mstruct_changed=s_mstruct_changed
     IF (PRESENT(rho_changed)) ks_env%rho_changed=rho_changed
     IF (PRESENT(potential_changed)) ks_env%potential_changed=potential_changed
     IF (PRESENT(forces_up_to_date)) ks_env%forces_up_to_date=forces_up_to_date
     IF (ks_env%s_mstruct_changed .OR. ks_env%potential_changed .OR. ks_env%rho_changed)&
          ks_env%forces_up_to_date=.FALSE.
     IF (PRESENT(v_hartree_rspace)) ks_env%v_hartree_rspace => v_hartree_rspace

     IF (PRESENT(matrix_h))    CALL set_1d_pointer(ks_env%matrix_h, matrix_h)
     IF (PRESENT(matrix_h_kp)) CALL set_2d_pointer(ks_env%matrix_h, matrix_h_kp)

     IF (PRESENT(matrix_ks))    CALL set_1d_pointer(ks_env%matrix_ks, matrix_ks)
     IF (PRESENT(matrix_ks_kp)) CALL set_2d_pointer(ks_env%matrix_ks, matrix_ks_kp)

     IF (PRESENT(matrix_vxc))    CALL set_1d_pointer(ks_env%matrix_vxc, matrix_vxc)
     IF (PRESENT(matrix_vxc_kp)) CALL set_2d_pointer(ks_env%matrix_vxc, matrix_vxc_kp)

     IF (PRESENT(matrix_s))    CALL set_1d_pointer(ks_env%matrix_s, matrix_s)
     IF (PRESENT(matrix_s_kp)) CALL set_2d_pointer(ks_env%matrix_s, matrix_s_kp)

     IF (PRESENT(kinetic))    CALL set_1d_pointer(ks_env%kinetic, kinetic)
     IF (PRESENT(kinetic_kp)) CALL set_2d_pointer(ks_env%kinetic, kinetic_kp)

     IF (PRESENT(matrix_w))    CALL set_1d_pointer(ks_env%matrix_w, matrix_w)
     IF (PRESENT(matrix_w_kp)) CALL set_2d_pointer(ks_env%matrix_w, matrix_w_kp)

     IF (PRESENT(matrix_ks_im)) ks_env%matrix_ks_im => matrix_ks_im
     IF (PRESENT(matrix_ks_aux_fit)) ks_env%matrix_ks_aux_fit => matrix_ks_aux_fit
     IF (PRESENT(matrix_ks_aux_fit_im)) ks_env%matrix_ks_aux_fit_im => matrix_ks_aux_fit_im
     IF (PRESENT(matrix_ks_aux_fit_dft)) ks_env%matrix_ks_aux_fit_dft => matrix_ks_aux_fit_dft
     IF (PRESENT(matrix_ks_aux_fit_hfx)) ks_env%matrix_ks_aux_fit_hfx => matrix_ks_aux_fit_hfx
     IF (PRESENT(matrix_s_aux_fit)) ks_env%matrix_s_aux_fit => matrix_s_aux_fit
     IF (PRESENT(matrix_s_aux_fit_vs_orb)) ks_env%matrix_s_aux_fit_vs_orb => matrix_s_aux_fit_vs_orb
     IF (PRESENT(matrix_w_mp2)) ks_env%matrix_w_mp2 => matrix_w_mp2
     IF (PRESENT(matrix_p_mp2)) ks_env%matrix_p_mp2 => matrix_p_mp2
     IF (PRESENT(gamma_matrix)) ks_env%gamma_matrix => gamma_matrix
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
     IF (PRESENT(sab_vdw)) ks_env%sab_vdw => sab_vdw
     IF (PRESENT(sab_scp)) ks_env%sab_scp => sab_scp
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
     IF (PRESENT(sab_kp)) ks_env%sab_kp => sab_kp
     IF (PRESENT(admm_dm)) ks_env%admm_dm => admm_dm

     IF (PRESENT(task_list)) ks_env%task_list=>task_list
     IF (PRESENT(task_list_aux_fit)) ks_env%task_list_aux_fit=>task_list_aux_fit
     IF (PRESENT(task_list_soft)) ks_env%task_list_soft=>task_list_soft

     IF (PRESENT(subsys)) THEN
       CALL qs_subsys_release(ks_env%subsys)
       CALL qs_subsys_retain(subsys)
       ks_env%subsys => subsys
    END IF
    IF (PRESENT(dft_control)) THEN
       CALL dft_control_retain(dft_control)
       CALL dft_control_release(ks_env%dft_control)
       ks_env%dft_control => dft_control
    END IF
    IF (PRESENT(dbcsr_dist)) THEN
       IF(ASSOCIATED(ks_env%dbcsr_dist)) THEN
          CALL cp_dbcsr_distribution_release(ks_env%dbcsr_dist)
          DEALLOCATE(ks_env%dbcsr_dist)
       ENDIF
       ks_env%dbcsr_dist => dbcsr_dist
    ENDIF
    IF (PRESENT(distribution_2d)) THEN
       CALL distribution_2d_retain(distribution_2d)
       CALL distribution_2d_release(ks_env%distribution_2d)
       ks_env%distribution_2d => distribution_2d
    ENDIF
    IF (PRESENT(pw_env)) THEN
      CALL pw_env_retain(pw_env)
      CALL pw_env_release(ks_env%pw_env)
      ks_env%pw_env => pw_env
    END IF
    IF (PRESENT(para_env)) THEN
       CALL cp_para_env_retain(para_env)
       CALL cp_para_env_release(ks_env%para_env)
       ks_env%para_env => para_env
    ENDIF
    IF (PRESENT(blacs_env)) THEN
       CALL cp_blacs_env_retain(blacs_env)
       CALL cp_blacs_env_release(ks_env%blacs_env)
       ks_env%blacs_env => blacs_env
    ENDIF
  END SUBROUTINE set_ks_env

! *****************************************************************************
!> \brief releases the ks_env (see doc/ReferenceCounting.html)
!> \param ks_env the ks_env to be released
!> \par History
!>      05.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE qs_ks_release(ks_env)
    TYPE(qs_ks_env_type), POINTER            :: ks_env

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_release', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(ks_env)) THEN
       CPASSERT(ks_env%ref_count>0)
       ks_env%ref_count=ks_env%ref_count-1

       IF (ks_env%ref_count<1) THEN
          IF (ASSOCIATED(ks_env%v_hartree_rspace)) &
            CALL pw_release(ks_env%v_hartree_rspace)
          IF (ASSOCIATED(ks_env%matrix_ks_im)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_ks_im)

          CALL kpoint_transitional_release(ks_env%matrix_ks)
          CALL kpoint_transitional_release(ks_env%matrix_h)
          CALL kpoint_transitional_release(ks_env%matrix_vxc)
          CALL kpoint_transitional_release(ks_env%matrix_s)
          CALL kpoint_transitional_release(ks_env%matrix_w)
          CALL kpoint_transitional_release(ks_env%kinetic)

          IF (ASSOCIATED(ks_env%matrix_ks_aux_fit)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_ks_aux_fit)
          IF (ASSOCIATED(ks_env%matrix_ks_aux_fit_im)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_ks_aux_fit_im)
          IF (ASSOCIATED(ks_env%matrix_ks_aux_fit_dft)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_ks_aux_fit_dft)
          IF (ASSOCIATED(ks_env%matrix_ks_aux_fit_hfx)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_ks_aux_fit_hfx)
          IF (ASSOCIATED(ks_env%matrix_s_aux_fit)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_s_aux_fit)
          IF (ASSOCIATED(ks_env%matrix_s_aux_fit_vs_orb)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_s_aux_fit_vs_orb)
          IF (ASSOCIATED(ks_env%matrix_w_mp2)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_w_mp2)
          IF (ASSOCIATED(ks_env%matrix_p_mp2)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_p_mp2)
          IF (ASSOCIATED(ks_env%gamma_matrix)) &
            CALL cp_dbcsr_deallocate_matrix_set(ks_env%gamma_matrix)
          IF(ASSOCIATED(ks_env%rho))&
            CALL qs_rho_release(ks_env%rho)
          IF(ASSOCIATED(ks_env%rho_buffer))&
            CALL qs_rho_release(ks_env%rho_buffer)
          IF(ASSOCIATED(ks_env%rho_xc))&
            CALL qs_rho_release(ks_env%rho_xc)
          IF(ASSOCIATED(ks_env%rho_aux_fit))&
            CALL qs_rho_release(ks_env%rho_aux_fit)
          IF(ASSOCIATED(ks_env%rho_aux_fit_buffer))&
            CALL qs_rho_release(ks_env%rho_aux_fit_buffer)
          IF (ASSOCIATED(ks_env%distribution_2d)) &
             CALL distribution_2d_release(ks_env%distribution_2d)
          IF(ASSOCIATED(ks_env%task_list)) &
             CALL deallocate_task_list(ks_env%task_list)
          IF(ASSOCIATED(ks_env%task_list_aux_fit)) &
             CALL deallocate_task_list(ks_env%task_list_aux_fit)
          IF(ASSOCIATED(ks_env%task_list_soft)) &
             CALL deallocate_task_list(ks_env%task_list_soft)

          IF (ASSOCIATED(ks_env%rho_nlcc_g)) THEN
             CALL pw_release(ks_env%rho_nlcc_g%pw)
             DEALLOCATE(ks_env%rho_nlcc_g)
          ENDIF
          IF (ASSOCIATED(ks_env%rho_nlcc)) THEN
             CALL pw_release(ks_env%rho_nlcc%pw)
             DEALLOCATE(ks_env%rho_nlcc)
          ENDIF
          IF (ASSOCIATED(ks_env%rho_core)) THEN
             CALL pw_release(ks_env%rho_core%pw)
             DEALLOCATE(ks_env%rho_core)
          ENDIF
          IF (ASSOCIATED(ks_env%vppl)) THEN
             CALL pw_release(ks_env%vppl%pw)
             DEALLOCATE(ks_env%vppl)
          ENDIF
          IF (ASSOCIATED(ks_env%vee)) THEN
             CALL pw_release(ks_env%vee%pw)
             DEALLOCATE(ks_env%vee)
          ENDIF
          IF(ASSOCIATED(ks_env%dbcsr_dist)) THEN
             CALL cp_dbcsr_distribution_release(ks_env%dbcsr_dist)
             DEALLOCATE(ks_env%dbcsr_dist)
          ENDIF

          CALL release_sab(ks_env%sab_orb)
          CALL release_sab(ks_env%sac_ae)
          CALL release_sab(ks_env%sac_ppl)
          CALL release_sab(ks_env%sap_ppnl)
          CALL release_sab(ks_env%sap_oce)
          CALL release_sab(ks_env%sab_se)
          CALL release_sab(ks_env%sab_vdw)
          CALL release_sab(ks_env%sab_scp)
          CALL release_sab(ks_env%sab_tbe)
          CALL release_sab(ks_env%sab_core)
          CALL release_sab(ks_env%sab_all)
          CALL release_sab(ks_env%sab_aux_fit)
          CALL release_sab(ks_env%sab_aux_fit_vs_orb)
          CALL release_sab(ks_env%sab_aux_fit_asymm)
          CALL release_sab(ks_env%sab_lrc)
          CALL release_sab(ks_env%sab_almo)
          CALL release_sab(ks_env%sab_kp)
          CALL dft_control_release(ks_env%dft_control)
          CALL kpoint_release(ks_env%kpoints)
          CALL qs_subsys_release(ks_env%subsys)
          CALL pw_env_release(ks_env%pw_env)
          CALL cp_para_env_release(ks_env%para_env)
          CALL cp_blacs_env_release(ks_env%blacs_env)
          CALL admm_dm_release(ks_env%admm_dm)

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
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE qs_ks_retain(ks_env)
    TYPE(qs_ks_env_type), POINTER            :: ks_env

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPASSERT(ASSOCIATED(ks_env))
  CPASSERT(ks_env%ref_count>0)
  ks_env%ref_count=ks_env%ref_count+1
END SUBROUTINE qs_ks_retain

! *****************************************************************************
!> \brief tells that some of the things relevant to the ks calculation
!>      did change. has to be called when changes happen otherwise
!>      the calculation will give wrong results.
!> \param ks_env the environement that is informed about the changes
!> \param s_mstruct_changed if true it means that the structure of the
!>        overlap matrix has changed
!>        (atoms have moved)
!> \param rho_changed if true it means that the density has changed
!> \param potential_changed ...
!> \param full_reset if true everything has changed
!> \par History
!>      4.2002 created [fawzi]
!>     12.2014 moved from qs_ks_methods, added deallocation of KS-matrices [Ole Schuett]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE qs_ks_did_change(ks_env,s_mstruct_changed,rho_changed,&
                              potential_changed,full_reset)
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    LOGICAL, INTENT(in), OPTIONAL            :: s_mstruct_changed, &
                                                rho_changed, &
                                                potential_changed, full_reset

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_did_change', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure, my_mstruct_chg

    CALL timeset(routineN,handle)
    failure = .FALSE.
    my_mstruct_chg = .FALSE.

    CPASSERT(ASSOCIATED(ks_env))

    IF (PRESENT(rho_changed)) THEN
       IF (rho_changed) ks_env%rho_changed = .TRUE.
    END IF

    IF (PRESENT(potential_changed)) THEN
       IF (potential_changed) ks_env%potential_changed = .TRUE.
    END IF

    IF (PRESENT(s_mstruct_changed)) THEN
       IF (s_mstruct_changed) my_mstruct_chg = .TRUE.
    END IF

    IF (PRESENT(full_reset)) THEN
       IF (full_reset) THEN
          ks_env%potential_changed = .TRUE.
          my_mstruct_chg = .TRUE.
       ENDIF
    END IF

    IF(my_mstruct_chg) THEN
       ks_env%s_mstruct_changed = .TRUE.
       ! *** deallocate matrices that will have the wrong structure ***
       CALL kpoint_transitional_release(ks_env%matrix_ks)
       CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_ks_aux_fit)
       !TODO: deallocate imaginary parts as well
       !CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_ks_im)
       !CALL cp_dbcsr_deallocate_matrix_set(ks_env%matrix_ks_aux_fit_im)
    ENDIF

    CALL  timestop(handle)

  END SUBROUTINE qs_ks_did_change

END MODULE qs_ks_types
