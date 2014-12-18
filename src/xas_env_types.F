!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief define create destroy get and put information
!>      in xas_env to calculate the x-ray absorption spectra
!> \par History
!>      created 05.2005
!> \author MI (05.2005)
! *****************************************************************************
MODULE xas_env_types

  USE basis_set_types,                 ONLY: deallocate_gto_basis_set,&
                                             gto_basis_set_p_type
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_deallocate_matrix_set,&
                                             cp_dbcsr_p_type
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             fm_pool_give_back_fm
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE kinds,                           ONLY: dp
  USE qs_loc_types,                    ONLY: qs_loc_env_new_type,&
                                             qs_loc_env_release,&
                                             qs_loc_env_retain
  USE qs_scf_types,                    ONLY: qs_scf_env_type,&
                                             scf_env_release,&
                                             scf_env_retain
  USE scf_control_types,               ONLY: scf_c_release,&
                                             scf_c_retain,&
                                             scf_control_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_env_types'

! *** Public data types ***

  PUBLIC :: xas_environment_type

! *** Public subroutines ***

  PUBLIC :: get_xas_env, set_xas_env, xas_env_create, xas_env_release, xas_env_retain

! *****************************************************************************
!> \param ref_count counter for pointers to xas_env
!> \param iter_count counter for the step at which xas is calculated
!> \param nao number of atomic orbitals in the basis
!> \param exc_state state that is now excited (this change atom by atom)
!> \param nvirtual number of empy states to take into account for the spectrum
!> \param state_of_atom for each atom the state that have to be excited (global index)
!>                      dimension is the number of atoms to be excited
!> \param atom_of_state atom to which each state is assigned,
!>                     dimension is the number of states occupied that might be excited
!> \param type_of_state character of the state (1s,2s,2p...)
!> \param spectrum for each excitation the energy and the oscillator strength
!> \param centers_wfn for each wfn the center of charge (optimized by localization)
!> \param groundstate_coeff temporary storage for the original mos coefficients
!> \param ostrength_sm sin and cos integrals computed for the contracted GTO functions
!> \param dip_fm_set fm for the sin and cos integrals to define the pos operator
!> \param fm_work temporary fm
!> \param qs_loc_env environment for the localization procedure
!> \par History
!>       created 05-2005
!> \author MI
! *****************************************************************************
  TYPE xas_environment_type
    INTEGER :: ref_count
    INTEGER :: iter_count
    INTEGER :: nao,exc_state, xas_estate
    INTEGER :: nexc_search, nexc_atoms
    INTEGER :: nvirtual,nvirtual2
    INTEGER :: unoccupied_max_iter
    INTEGER, DIMENSION(:), POINTER :: state_of_atom
    INTEGER, DIMENSION(:), POINTER :: atom_of_state
    INTEGER, DIMENSION(:), POINTER :: type_of_state
    INTEGER, DIMENSION(:), POINTER :: mykind_of_atom
    INTEGER, DIMENSION(:), POINTER :: mykind_of_kind
    INTEGER, DIMENSION(:), POINTER :: exc_atoms
    REAL(dp) :: ip_energy, occ_estate, unoccupied_eps, xas_nelectron, homo_occ
    REAL(dp), DIMENSION(:), POINTER :: all_evals
    REAL(dp), DIMENSION(:), POINTER :: unoccupied_evals
    REAL(dp), DIMENSION(:,:), POINTER :: spectrum
    REAL(dp), DIMENSION(:,:), POINTER :: centers_wfn
    TYPE(cp_2d_r_p_type), DIMENSION(:), POINTER :: stogto_overlap
    TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER :: my_gto_basis
    TYPE( cp_fm_p_type ), DIMENSION(:),  POINTER :: groundstate_coeff
    TYPE( cp_fm_p_type ), DIMENSION(:,:),  POINTER :: dip_fm_set
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
       POINTER                                   :: ao_mo_fm_pools
    TYPE( cp_fm_type ), POINTER :: excvec_coeff
    TYPE( cp_fm_type ), POINTER :: excvec_overlap
    TYPE( cp_fm_type ), POINTER :: unoccupied_orbs
    TYPE( cp_fm_type ), POINTER :: all_vectors
    TYPE( cp_fm_type ), POINTER :: fm_work
    TYPE( cp_dbcsr_p_type ), DIMENSION(:), POINTER :: ostrength_sm
    TYPE( qs_loc_env_new_type ), POINTER :: qs_loc_env
    TYPE(qs_scf_env_type), POINTER                        :: scf_env
    TYPE(scf_control_type), POINTER          :: scf_control


  END TYPE xas_environment_type

CONTAINS
! *****************************************************************************
!> \brief ...
!> \param xas_env ...
!> \param iter_count ...
!> \param exc_state ...
!> \param nao ...
!> \param nvirtual ...
!> \param nvirtual2 ...
!> \param centers_wfn ...
!> \param atom_of_state ...
!> \param exc_atoms ...
!> \param type_of_state ...
!> \param mykind_of_atom ...
!> \param mykind_of_kind ...
!> \param state_of_atom ...
!> \param spectrum ...
!> \param groundstate_coeff ...
!> \param ostrength_sm ...
!> \param dip_fm_set ...
!> \param excvec_coeff ...
!> \param excvec_overlap ...
!> \param fm_work ...
!> \param unoccupied_orbs ...
!> \param unoccupied_evals ...
!> \param unoccupied_max_iter ...
!> \param unoccupied_eps ...
!> \param all_vectors ...
!> \param all_evals ...
!> \param my_gto_basis ...
!> \param qs_loc_env ...
!> \param stogto_overlap ...
!> \param occ_estate ...
!> \param xas_nelectron ...
!> \param xas_estate ...
!> \param nexc_atoms ...
!> \param nexc_search ...
!> \param scf_env ...
!> \param scf_control ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE get_xas_env(xas_env,iter_count,exc_state,nao,nvirtual,nvirtual2,&
             centers_wfn,atom_of_state, exc_atoms, type_of_state,mykind_of_atom,mykind_of_kind,&
             state_of_atom,spectrum,groundstate_coeff,ostrength_sm,&
             dip_fm_set,excvec_coeff,excvec_overlap,fm_work,&
             unoccupied_orbs, unoccupied_evals,unoccupied_max_iter,unoccupied_eps,&
             all_vectors,all_evals, my_gto_basis,qs_loc_env,&
             stogto_overlap,occ_estate,xas_nelectron,xas_estate,nexc_atoms,nexc_search,scf_env,scf_control,error)

    TYPE(xas_environment_type), POINTER      :: xas_env
    INTEGER, INTENT(OUT), OPTIONAL           :: iter_count, exc_state, nao, &
                                                nvirtual, nvirtual2
    REAL(dp), DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: centers_wfn
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: atom_of_state, exc_atoms, &
      type_of_state, mykind_of_atom, mykind_of_kind, state_of_atom
    REAL(dp), DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: spectrum
    TYPE(cp_fm_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: groundstate_coeff
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: ostrength_sm
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: dip_fm_set
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: excvec_coeff, excvec_overlap, &
                                                fm_work, unoccupied_orbs
    REAL(dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: unoccupied_evals
    INTEGER, INTENT(OUT), OPTIONAL           :: unoccupied_max_iter
    REAL(dp), OPTIONAL                       :: unoccupied_eps
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: all_vectors
    REAL(dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: all_evals
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: my_gto_basis
    TYPE(qs_loc_env_new_type), OPTIONAL, &
      POINTER                                :: qs_loc_env
    TYPE(cp_2d_r_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: stogto_overlap
    REAL(dp), INTENT(OUT), OPTIONAL          :: occ_estate, xas_nelectron
    INTEGER, INTENT(OUT), OPTIONAL           :: xas_estate, nexc_atoms, &
                                                nexc_search
    TYPE(qs_scf_env_type), OPTIONAL, POINTER :: scf_env
    TYPE(scf_control_type), OPTIONAL, &
      POINTER                                :: scf_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_xas_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

     failure =.FALSE.
     CPPrecondition(ASSOCIATED(xas_env),cp_failure_level,routineP,error,failure)

     IF(.NOT.failure) THEN

       IF(PRESENT(iter_count)) iter_count=xas_env%iter_count
       IF(PRESENT(exc_state)) exc_state=xas_env%exc_state
       IF(PRESENT(nao)) nao=xas_env%nao
       IF(PRESENT(nvirtual)) nvirtual=xas_env%nvirtual
       IF(PRESENT(nvirtual2)) nvirtual2=xas_env%nvirtual2
       IF(PRESENT(xas_nelectron))  xas_nelectron=xas_env%xas_nelectron
       IF(PRESENT(occ_estate))  occ_estate=xas_env%occ_estate
       IF(PRESENT(xas_estate))  xas_estate=xas_env%xas_estate
       IF(PRESENT(nexc_search))  nexc_search=xas_env%nexc_search
       IF(PRESENT(nexc_atoms))  nexc_atoms=xas_env%nexc_atoms
       IF(PRESENT(unoccupied_eps))  unoccupied_eps=xas_env%unoccupied_eps
       IF(PRESENT(unoccupied_max_iter)) unoccupied_max_iter=xas_env%unoccupied_max_iter
       IF(PRESENT(centers_wfn)) centers_wfn=>xas_env%centers_wfn
       IF(PRESENT(atom_of_state)) atom_of_state=>xas_env%atom_of_state
       IF(PRESENT(exc_atoms)) exc_atoms=>xas_env%exc_atoms
       IF(PRESENT(type_of_state)) type_of_state=>xas_env%type_of_state
       IF(PRESENT(state_of_atom)) state_of_atom=>xas_env%state_of_atom
       IF(PRESENT(mykind_of_atom)) mykind_of_atom=>xas_env%mykind_of_atom
       IF(PRESENT(mykind_of_kind)) mykind_of_kind=>xas_env%mykind_of_kind
       IF(PRESENT(unoccupied_evals)) unoccupied_evals=>xas_env%unoccupied_evals
       IF(PRESENT(all_evals)) all_evals=>xas_env%all_evals
       IF(PRESENT(spectrum)) spectrum=>xas_env%spectrum
       IF(PRESENT(groundstate_coeff)) groundstate_coeff=>xas_env%groundstate_coeff
       IF(PRESENT(ostrength_sm)) ostrength_sm=>xas_env%ostrength_sm
       IF(PRESENT(excvec_overlap)) excvec_overlap=>xas_env%excvec_overlap
       IF(PRESENT(unoccupied_orbs)) unoccupied_orbs=>xas_env%unoccupied_orbs
       IF(PRESENT(fm_work)) fm_work=>xas_env%fm_work
       IF(PRESENT(all_vectors)) all_vectors=>xas_env%all_vectors
       IF(PRESENT(dip_fm_set)) dip_fm_set=>xas_env%dip_fm_set
       IF(PRESENT(qs_loc_env)) qs_loc_env=>xas_env%qs_loc_env
       IF(PRESENT(excvec_coeff)) excvec_coeff=>xas_env%excvec_coeff
       IF(PRESENT(my_gto_basis)) my_gto_basis => xas_env%my_gto_basis
       IF(PRESENT(stogto_overlap)) stogto_overlap => xas_env%stogto_overlap
       IF (PRESENT(scf_env)) scf_env => xas_env%scf_env
       IF (PRESENT(scf_control)) scf_control => xas_env%scf_control
     END IF
  END SUBROUTINE get_xas_env

! *****************************************************************************
!> \brief ...
!> \param xas_env ...
!> \param iter_count ...
!> \param nexc_search ...
!> \param nexc_atoms ...
!> \param nvirtual ...
!> \param nvirtual2 ...
!> \param ip_energy ...
!> \param occ_estate ...
!> \param qs_loc_env ...
!> \param xas_estate ...
!> \param xas_nelectron ...
!> \param homo_occ ...
!> \param scf_env ...
!> \param scf_control ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE set_xas_env(xas_env,iter_count,nexc_search,nexc_atoms,&
             nvirtual,nvirtual2, ip_energy, occ_estate, qs_loc_env,&
             xas_estate, xas_nelectron, homo_occ, scf_env, scf_control,error)

    TYPE(xas_environment_type), POINTER      :: xas_env
    INTEGER, INTENT(IN), OPTIONAL            :: iter_count, nexc_search, &
                                                nexc_atoms, nvirtual, &
                                                nvirtual2
    REAL(dp), INTENT(IN), OPTIONAL           :: ip_energy, occ_estate
    TYPE(qs_loc_env_new_type), OPTIONAL, &
      POINTER                                :: qs_loc_env
    INTEGER, INTENT(IN), OPTIONAL            :: xas_estate
    REAL(dp), INTENT(IN), OPTIONAL           :: xas_nelectron, homo_occ
    TYPE(qs_scf_env_type), OPTIONAL, POINTER :: scf_env
    TYPE(scf_control_type), OPTIONAL, &
      POINTER                                :: scf_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_xas_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure =.FALSE.

    CPPrecondition(ASSOCIATED(xas_env),cp_failure_level,routineP,error,failure)

    IF(.NOT. failure) THEN

      IF (PRESENT(iter_count)) xas_env%iter_count = iter_count
      IF (PRESENT(nexc_search)) xas_env%nexc_search = nexc_search
      IF (PRESENT(nexc_atoms)) xas_env%nexc_atoms = nexc_atoms
      IF (PRESENT(nvirtual)) xas_env%nvirtual = nvirtual
      IF (PRESENT(nvirtual2)) xas_env%nvirtual2 = nvirtual2
      IF (PRESENT(occ_estate)) xas_env%occ_estate = occ_estate
      IF (PRESENT(xas_nelectron)) xas_env%xas_nelectron = xas_nelectron
      IF (PRESENT(homo_occ)) xas_env%homo_occ = homo_occ
      IF (PRESENT(xas_estate)) xas_env%xas_estate = xas_estate
      IF (PRESENT(ip_energy)) xas_env%ip_energy = ip_energy
      IF (PRESENT(qs_loc_env)) THEN
        CALL  qs_loc_env_retain(qs_loc_env, error=error)
        IF(ASSOCIATED(xas_env%qs_loc_env)) &
            CALL qs_loc_env_release(xas_env%qs_loc_env,error=error)
        xas_env%qs_loc_env => qs_loc_env
      END IF
      IF (PRESENT(scf_env)) THEN ! accept also null pointers ?
         CALL scf_env_retain(scf_env,error=error)
         CALL scf_env_release(xas_env%scf_env, error=error)
         xas_env%scf_env => scf_env
      END IF
      IF (PRESENT(scf_control)) THEN ! accept also null pointers?
         CALL scf_c_retain(scf_control,error=error)
         CALL scf_c_release(xas_env%scf_control,error=error)
         xas_env%scf_control => scf_control
      END IF

    END IF

  END SUBROUTINE set_xas_env

! *****************************************************************************
!> \brief ...
!> \param xas_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE xas_env_create(xas_env, error)

    TYPE(xas_environment_type), POINTER      :: xas_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xas_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

    failure=.FALSE.
    ALLOCATE(xas_env,STAT=istat)
    CPPostcondition(istat==0,cp_warning_level,routineP,error,failure)

    IF (.NOT.failure) THEN

      xas_env%ref_count=1
      xas_env%iter_count=0
      xas_env%nvirtual=0
      xas_env%nvirtual2=0

      NULLIFY(xas_env%ao_mo_fm_pools)
      NULLIFY(xas_env%my_gto_basis)
      NULLIFY(xas_env%state_of_atom, xas_env%atom_of_state)
      NULLIFY(xas_env%exc_atoms)
      NULLIFY(xas_env%excvec_coeff,xas_env%excvec_overlap)
      NULLIFY(xas_env%type_of_state,xas_env%mykind_of_atom)
      NULLIFY(xas_env%type_of_state,xas_env%mykind_of_kind)
      NULLIFY(xas_env%fm_work,xas_env%groundstate_coeff,xas_env%dip_fm_set)
      NULLIFY(xas_env%ostrength_sm,xas_env%qs_loc_env,xas_env%spectrum)
      NULLIFY(xas_env%all_evals,xas_env%all_vectors)
      NULLIFY(xas_env%unoccupied_evals,xas_env%unoccupied_orbs)
      NULLIFY(xas_env%stogto_overlap)
      NULLIFY(xas_env%scf_env)
      NULLIFY(xas_env%scf_control)

    END IF

   END SUBROUTINE xas_env_create

! *****************************************************************************
!> \brief ...
!> \param xas_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE xas_env_release(xas_env, error)

    TYPE(xas_environment_type), POINTER      :: xas_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xas_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ik, istat, j
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(xas_env)) THEN
       CPPreconditionNoFail(xas_env%ref_count>0,cp_failure_level,routineP,error)
       xas_env%ref_count = xas_env%ref_count -1
       IF (xas_env%ref_count==0) THEN

         DEALLOCATE(xas_env%state_of_atom,xas_env%atom_of_state,STAT=istat)
         CPPostcondition(istat==0,cp_warning_level,routineP,error,failure)
         DEALLOCATE(xas_env%type_of_state,STAT=istat)
         CPPostcondition(istat==0,cp_warning_level,routineP,error,failure)
         DEALLOCATE(xas_env%mykind_of_atom,STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         DEALLOCATE(xas_env%mykind_of_kind,STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         DEALLOCATE(xas_env%exc_atoms,STAT=istat)
         CPPostcondition(istat==0,cp_warning_level,routineP,error,failure)
         DEALLOCATE(xas_env%centers_wfn,STAT=istat)
         CPPostcondition(istat==0,cp_warning_level,routineP,error,failure)
         IF(ASSOCIATED(xas_env%all_evals)) THEN
           DEALLOCATE(xas_env%all_evals,STAT=istat)
           CPPostcondition(istat==0,cp_warning_level,routineP,error,failure)
         END IF
         IF(ASSOCIATED(xas_env%unoccupied_evals)) THEN
           DEALLOCATE(xas_env%unoccupied_evals,STAT=istat)
           CPPostcondition(istat==0,cp_warning_level,routineP,error,failure)
         END IF
         IF(ASSOCIATED(xas_env%groundstate_coeff)) THEN
           DO i = 1,SIZE(xas_env%groundstate_coeff)
             CALL fm_pool_give_back_fm(xas_env%ao_mo_fm_pools(i)%pool,&
                  xas_env%groundstate_coeff(i)%matrix, error=error)
           END DO
           DEALLOCATE(xas_env%groundstate_coeff,STAT=istat)
           CPPostcondition(istat==0,cp_warning_level,routineP,error,failure)
         END IF
         IF(ASSOCIATED(xas_env%dip_fm_set)) THEN
           DO i = 1,SIZE(xas_env%dip_fm_set,2)
             DO j = 1,SIZE(xas_env%dip_fm_set,1)
               CALL cp_fm_release(xas_env%dip_fm_set(j,i)%matrix, error=error)
             END DO
           END DO
           DEALLOCATE(xas_env%dip_fm_set,STAT=istat)
           CPPostcondition(istat==0,cp_warning_level,routineP,error,failure)
         END IF

         IF(ASSOCIATED(xas_env%excvec_coeff)) THEN
           CALL cp_fm_release(xas_env%excvec_coeff,error=error)
         END IF
         IF(ASSOCIATED(xas_env%excvec_overlap)) THEN
           CALL cp_fm_release(xas_env%excvec_overlap,error=error)
         END IF
         IF(ASSOCIATED(xas_env%unoccupied_orbs)) THEN
           CALL cp_fm_release(xas_env%unoccupied_orbs,error=error)
         END IF
         IF(ASSOCIATED(xas_env%fm_work)) THEN
           CALL cp_fm_release(xas_env%fm_work,error=error)
         END IF
         NULLIFY(xas_env%ao_mo_fm_pools)
         IF(ASSOCIATED(xas_env%all_vectors).AND.xas_env%nvirtual .GT. 0) THEN
           CALL cp_fm_release(xas_env%all_vectors,error=error)
         ELSE
           NULLIFY(xas_env%all_vectors)
         END IF

         IF(ASSOCIATED(xas_env%ostrength_sm)) THEN
            CALL cp_dbcsr_deallocate_matrix_set(xas_env%ostrength_sm,error=error)
         END IF
         IF(ASSOCIATED(xas_env%qs_loc_env)) THEN
            CALL qs_loc_env_release(xas_env%qs_loc_env, error=error)
         END IF

         IF(ASSOCIATED(xas_env%my_gto_basis)) THEN
           DO ik = 1,SIZE(xas_env%my_gto_basis,1)
             CALL deallocate_gto_basis_set(xas_env%my_gto_basis(ik)%gto_basis_set,error)
           END DO
           DEALLOCATE(xas_env%my_gto_basis,STAT=istat)
           CPPostcondition(istat==0,cp_warning_level,routineP,error,failure)
         END IF

         IF(ASSOCIATED(xas_env%stogto_overlap)) THEN
           DO ik = 1,SIZE(xas_env%stogto_overlap,1)
             DEALLOCATE(xas_env%stogto_overlap(ik)%array)
           END DO
           DEALLOCATE(xas_env%stogto_overlap,STAT=istat)
           CPPostcondition(istat==0,cp_warning_level,routineP,error,failure)
         END IF

         CALL scf_env_release(xas_env%scf_env, error=error)
         CALL scf_c_release(xas_env%scf_control, error=error)


         DEALLOCATE(xas_env,STAT=istat)
         CPPostcondition(istat==0,cp_warning_level,routineP,error,failure)
         NULLIFY(xas_env)
       END IF
    END IF

  END SUBROUTINE xas_env_release

! *****************************************************************************
!> \brief ...
!> \param xas_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE xas_env_retain(xas_env, error)

    TYPE(xas_environment_type), POINTER      :: xas_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xas_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(xas_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(xas_env%ref_count>0,cp_failure_level,routineP,error,failure)
       xas_env%ref_count=xas_env%ref_count+1
    END IF
  END SUBROUTINE xas_env_retain

END MODULE xas_env_types

