!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of kinetic energy matrix and forces
!> \par History
!>      JGH: from core_hamiltonian
!>      simplify further [7.2014]
!> \author Juerg Hutter
! *****************************************************************************
MODULE qs_kinetic

  USE ai_contraction,                  ONLY: block_add,&
                                             contraction,&
                                             decontraction,&
                                             force_trace
  USE ai_kinetic,                      ONLY: kinetic
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                             gto_basis_set_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_allocate_matrix_set,&
                                             cp_dbcsr_filter,&
                                             cp_dbcsr_finalize,&
                                             cp_dbcsr_get_block_p,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_type,&
                                             dbcsr_distribution_obj
  USE kinds,                           ONLY: dp
  USE orbital_pointers,                ONLY: ncoset
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_integral_utils,               ONLY: basis_set_list_setup,&
                                             get_memory_usage
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             get_neighbor_list_set_p,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE qs_overlap,                      ONLY: create_sab_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_kinetic'

! *** Public subroutines ***

  PUBLIC :: build_kinetic_matrix

CONTAINS

! *****************************************************************************

!> \brief   Calculation of the kinetic energy matrix over Cartesian Gaussian functions.
!> \param   qs_env the QS environment
!> \param   matrix_t The kinetic energy matrix to be calculated
!> \param   matrix_name The name of the matrix (i.e. for output)
!> \param   basis_set_id basis set to be used
!> \param   sab_nl pair list (must be consistent with basis sets!)
!> \param   calculate_forces (optional)
!> \param   matrix_p density matrix for force calculation (optional)
!> \param   eps_filter Filter final matrix (optional)
!> \param   error for error handling
!> \date    11.10.2010
!> \par     History
!>          Ported from qs_overlap, replaces code in build_core_hamiltonian
!>          Refactoring [07.2014] JGH
!>          Simplify options and use new kinetic energy integral routine
!> \author  JGH
!> \version 1.0
! *****************************************************************************
  SUBROUTINE build_kinetic_matrix(qs_env,matrix_t,matrix_name,&
             basis_set_id,sab_nl,calculate_forces,matrix_p,eps_filter,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_t
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: matrix_name
    INTEGER, INTENT(IN)                      :: basis_set_id
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_nl
    LOGICAL, INTENT(IN), OPTIONAL            :: calculate_forces
    TYPE(cp_dbcsr_type), OPTIONAL, POINTER   :: matrix_p
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: eps_filter
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_kinetic_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, iatom, icol, ikind, inode, irow, iset, &
      istat, jatom, jkind, jset, last_jatom, ldsab, natom, ncoa, ncob, nkind, &
      nseta, nsetb, sgfa, sgfb
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: do_forces, do_symmetric, &
                                                failure, found, new_atom_b, &
                                                trans, use_virial
    REAL(KIND=dp)                            :: f0, ff, rab2, tab
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: kab, pab, qab
    REAL(KIND=dp), DIMENSION(3)              :: force_a, rab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: k_block, p_block, rpgfa, &
                                                rpgfb, scon_a, scon_b, zeta, &
                                                zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: dab
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(dbcsr_distribution_obj), POINTER    :: dbcsr_dist
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(virial_type), POINTER               :: virial

     failure = .FALSE.

     CALL timeset(routineN,handle)

     NULLIFY (atomic_kind_set, p_block)
     CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,natom=natom,&
                     dbcsr_dist=dbcsr_dist,error=error)
     nkind = SIZE(atomic_kind_set)

     ALLOCATE (atom_of_kind(natom),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind)

     do_forces = .FALSE.
     IF (PRESENT(calculate_forces)) do_forces = calculate_forces

     ! check for symmetry
     CPPrecondition(SIZE(sab_nl) > 0,cp_failure_level,routineP,error,failure)
     CALL get_neighbor_list_set_p(neighbor_list_sets=sab_nl,symmetric=do_symmetric)

     CALL cp_dbcsr_allocate_matrix_set(matrix_t,1,error=error)
     CALL create_sab_matrix(qs_env, matrix_t, matrix_name, basis_set_id, basis_set_id, &
                            sab_nl, do_symmetric, error)

     IF (do_forces) THEN
        ! if forces -> maybe virial too
        CALL get_qs_env(qs_env=qs_env,force=force,virial=virial,error=error)
        use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)
        ! we need density matrix for forces
        CPPrecondition(PRESENT(matrix_p),cp_failure_level,routineP,error,failure)
     END IF

     ! *** Allocate work storage ***
     ldsab = get_memory_usage(atomic_kind_set,basis_set_id)
     ALLOCATE (kab(ldsab,ldsab),qab(ldsab,ldsab),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     IF (do_forces) THEN
        ALLOCATE(dab(ldsab,ldsab,3),pab(ldsab,ldsab),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     END IF

     ! prepare basis set
     ALLOCATE (basis_set_list(nkind),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     CALL basis_set_list_setup(basis_set_list,basis_set_id,atomic_kind_set)

     ! Iterate of neighbor list
     CALL neighbor_list_iterator_create(nl_iterator,sab_nl)
     DO WHILE (neighbor_list_iterate(nl_iterator)==0)
        CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,inode=inode,&
                               iatom=iatom,jatom=jatom,r=rab)
        atom_a = atom_of_kind(iatom)
        atom_b = atom_of_kind(jatom)
        basis_set_a => basis_set_list(ikind)%gto_basis_set
        IF (.NOT.ASSOCIATED(basis_set_a)) CYCLE
        basis_set_b => basis_set_list(jkind)%gto_basis_set
        IF (.NOT.ASSOCIATED(basis_set_b)) CYCLE
        ! basis ikind
        first_sgfa   =>  basis_set_a%first_sgf
        la_max       =>  basis_set_a%lmax
        la_min       =>  basis_set_a%lmin
        npgfa        =>  basis_set_a%npgf
        nseta        =   basis_set_a%nset
        nsgfa        =>  basis_set_a%nsgf_set
        rpgfa        =>  basis_set_a%pgf_radius
        set_radius_a =>  basis_set_a%set_radius
        scon_a       =>  basis_set_a%scon
        zeta         =>  basis_set_a%zet
        ! basis jkind
        first_sgfb   =>  basis_set_b%first_sgf
        lb_max       =>  basis_set_b%lmax
        lb_min       =>  basis_set_b%lmin
        npgfb        =>  basis_set_b%npgf
        nsetb        =   basis_set_b%nset
        nsgfb        =>  basis_set_b%nsgf_set
        rpgfb        =>  basis_set_b%pgf_radius
        set_radius_b =>  basis_set_b%set_radius
        scon_b       =>  basis_set_b%scon
        zetb         =>  basis_set_b%zet

        IF(inode==1) last_jatom = 0
        IF (jatom /= last_jatom) THEN
          new_atom_b = .TRUE.
          last_jatom = jatom
        ELSE
          new_atom_b = .FALSE.
        END IF

        IF (new_atom_b) THEN
          IF( do_symmetric ) THEN
            IF (iatom <= jatom) THEN
              irow = iatom
              icol = jatom
            ELSE
              irow = jatom
              icol = iatom
            END IF
            f0 = 2.0_dp
            IF (iatom==jatom) f0 = 1.0_dp
            ff = 2.0_dp
          ELSE
            irow = iatom
            icol = jatom
            f0 = 1.0_dp
            ff = 1.0_dp
          END IF
          NULLIFY (k_block)
          CALL cp_dbcsr_get_block_p(matrix=matrix_t(1)%matrix,&
               row=irow,col=icol,BLOCK=k_block,found=found)
          CPPostcondition(found,cp_failure_level,routineP,error,failure)
          IF (do_forces) THEN
             NULLIFY (p_block)
             CALL cp_dbcsr_get_block_p(matrix=matrix_p,row=irow,col=icol,&
                  block=p_block,found=found)
             CPPostcondition(found,cp_failure_level,routineP,error,failure)
          END IF
        END IF

        rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
        tab = SQRT(rab2)
        trans = do_symmetric .AND. (iatom > jatom)

        DO iset=1,nseta

          ncoa = npgfa(iset)*(ncoset(la_max(iset))-ncoset(la_min(iset)-1))
          sgfa = first_sgfa(1,iset)

          DO jset=1,nsetb

            IF (set_radius_a(iset) + set_radius_b(jset) < tab) CYCLE

            ncob = npgfb(jset)*(ncoset(lb_max(jset))-ncoset(lb_min(jset)-1))
            sgfb = first_sgfb(1,jset)

            IF (do_forces .AND. ASSOCIATED(p_block) .AND. ((iatom/=jatom) .OR. use_virial)) THEN
               ! Decontract P matrix block
               kab = 0.0_dp
               CALL block_add("OUT",kab,nsgfa(iset),nsgfb(jset),p_block,sgfa,sgfb,trans=trans,error=error)
               CALL decontraction(kab,pab,scon_a(:,sgfa:),ncoa,nsgfa(iset),scon_b(:,sgfb:),ncob,nsgfb(jset),&
                                  trans=trans,error=error)
               ! calculate integrals and derivatives
               CALL kinetic(la_max(iset),la_min(iset),npgfa(iset),rpgfa(:,iset),zeta(:,iset),&
                            lb_max(jset),lb_min(jset),npgfb(jset),rpgfb(:,jset),zetb(:,jset),&
                            rab,kab,dab,error=error)
               CALL force_trace(force_a,dab,pab,ncoa,ncob,3,error=error)
               force(ikind)%kinetic(:,atom_a)=force(ikind)%kinetic(:,atom_a) + ff*force_a(:)
               force(jkind)%kinetic(:,atom_b)=force(jkind)%kinetic(:,atom_b) - ff*force_a(:)
               IF( use_virial ) THEN
                  CALL virial_pair_force ( virial%pv_virial, f0, force_a, rab, error)
               END IF
            ELSE
               ! calclulate integrals
               CALL kinetic(la_max(iset),la_min(iset),npgfa(iset),rpgfa(:,iset),zeta(:,iset),&
                            lb_max(jset),lb_min(jset),npgfb(jset),rpgfb(:,jset),zetb(:,jset),&
                            rab,kab,error=error)
            END IF
            ! Contraction step 
            CALL contraction(kab,qab,ca=scon_a(:,sgfa:),na=ncoa,ma=nsgfa(iset),&
                                     cb=scon_b(:,sgfb:),nb=ncob,mb=nsgfb(jset),&
                                     trans=trans,error=error)
            CALL block_add("IN",qab,nsgfa(iset),nsgfb(jset),k_block,sgfa,sgfb,trans=trans,error=error)

          END DO
        END DO

     END DO
     CALL neighbor_list_iterator_release(nl_iterator)

     CALL cp_dbcsr_finalize(matrix_t(1)%matrix, error=error)
     IF (PRESENT(eps_filter)) THEN
        CALL cp_dbcsr_filter(matrix_t(1)%matrix, eps_filter, error=error)
     END IF

     ! Release work storage
     DEALLOCATE (atom_of_kind,STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     DEALLOCATE (kab,qab,STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     IF(do_forces) THEN
        DEALLOCATE (pab,dab,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     END IF
     DEALLOCATE (basis_set_list,STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

     CALL timestop(handle)

  END SUBROUTINE build_kinetic_matrix

END MODULE qs_kinetic

