!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_allocate_matrix_set,&
                                             cp_dbcsr_filter,&
                                             cp_dbcsr_finalize,&
                                             cp_dbcsr_get_block_p,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_type,&
                                             dbcsr_distribution_obj
  USE kinds,                           ONLY: dp
  USE kpoint_types,                    ONLY: get_kpoint_info,&
                                             kpoint_type
  USE orbital_pointers,                ONLY: ncoset
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_integral_utils,               ONLY: basis_set_list_setup,&
                                             get_memory_usage
  USE qs_kind_types,                   ONLY: qs_kind_type
  USE qs_ks_types,                     ONLY: get_ks_env,&
                                             qs_ks_env_type
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
!> \param   ks_env the QS environment
!> \param   matrix_t The kinetic energy matrix to be calculated (optional)
!> \param   matrixkp_t The kinetic energy matrices to be calculated (kpoints,optional)
!> \param   matrix_name The name of the matrix (i.e. for output)
!> \param   basis_type basis set to be used
!> \param   sab_nl pair list (must be consistent with basis sets!)
!> \param   calculate_forces (optional)
!> \param   matrix_p density matrix for force calculation (optional)
!> \param   matrixkp_p density matrix for force calculation with kpoints (optional)
!> \param   eps_filter Filter final matrix (optional)
!> \date    11.10.2010
!> \par     History
!>          Ported from qs_overlap, replaces code in build_core_hamiltonian
!>          Refactoring [07.2014] JGH
!>          Simplify options and use new kinetic energy integral routine
!>          kpoints [08.2014] JGH
!> \author  JGH
!> \version 1.0
! *****************************************************************************
  SUBROUTINE build_kinetic_matrix(ks_env,matrix_t,matrixkp_t,matrix_name,&
             basis_type,sab_nl,calculate_forces,matrix_p,matrixkp_p,&
             eps_filter)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: matrix_t
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: matrixkp_t
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: matrix_name
    CHARACTER(LEN=*), INTENT(IN)             :: basis_type
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_nl
    LOGICAL, INTENT(IN), OPTIONAL            :: calculate_forces
    TYPE(cp_dbcsr_type), OPTIONAL, POINTER   :: matrix_p
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: matrixkp_p
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: eps_filter

    CHARACTER(len=*), PARAMETER :: routineN = 'build_kinetic_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, iatom, ic, icol, ikind, inode, irow, &
      iset, istat, jatom, jkind, jset, ldsab, natom, ncoa, ncob, nimg, nkind, &
      nseta, nsetb, sgfa, sgfb
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(3)                    :: cell
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    INTEGER, DIMENSION(:, :, :), POINTER     :: cell_to_index
    LOGICAL                                  :: do_forces, do_symmetric, &
                                                dokp, failure, found, trans, &
                                                use_cell_mapping, use_virial
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
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b
    TYPE(kpoint_type), POINTER               :: kpoints
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(virial_type), POINTER               :: virial

     failure = .FALSE.

     CALL timeset(routineN,handle)

     ! test for matrices (kpoints or standard gamma point)
     IF (PRESENT(matrix_t)) THEN
        dokp = .FALSE.
        use_cell_mapping = .FALSE.
     ELSEIF (PRESENT(matrixkp_t)) THEN
        dokp = .TRUE.
        CALL get_ks_env(ks_env=ks_env,kpoints=kpoints)
        CALL get_kpoint_info(kpoint=kpoints,cell_to_index=cell_to_index)
        use_cell_mapping = (SIZE(cell_to_index) > 1)
     ELSE
        CPAssert(.FALSE.,cp_failure_level,routineP,failure)
     END IF

     NULLIFY (atomic_kind_set, qs_kind_set, p_block,dft_control)
     CALL get_ks_env(ks_env,&
                     dft_control=dft_control,&
                     atomic_kind_set=atomic_kind_set,&
                     natom=natom,&
                     qs_kind_set=qs_kind_set,&
                     dbcsr_dist=dbcsr_dist)

     nimg = dft_control%nimages
     nkind = SIZE(atomic_kind_set)

     ALLOCATE (atom_of_kind(natom),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,failure)
     CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind)

     do_forces = .FALSE.
     IF (PRESENT(calculate_forces)) do_forces = calculate_forces

     ! check for symmetry
     CPPrecondition(SIZE(sab_nl) > 0,cp_failure_level,routineP,failure)
     CALL get_neighbor_list_set_p(neighbor_list_sets=sab_nl,symmetric=do_symmetric)

     ! prepare basis set
     ALLOCATE (basis_set_list(nkind),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,failure)
     CALL basis_set_list_setup(basis_set_list,basis_type,qs_kind_set)

     IF(dokp) THEN
        CALL cp_dbcsr_allocate_matrix_set(matrixkp_t,1,nimg)
        CALL create_sab_matrix(ks_env, matrixkp_t, matrix_name, basis_set_list, basis_set_list, &
                               sab_nl, do_symmetric)
     ELSE
        CALL cp_dbcsr_allocate_matrix_set(matrix_t,1)
        CALL create_sab_matrix(ks_env, matrix_t, matrix_name, basis_set_list, basis_set_list, &
                               sab_nl, do_symmetric)
     END IF

     IF (do_forces) THEN
        ! if forces -> maybe virial too
        CALL get_ks_env(ks_env=ks_env,force=force,virial=virial)
        use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)
        ! we need density matrix for forces
        IF(dokp) THEN
           CPPrecondition(PRESENT(matrixkp_p),cp_failure_level,routineP,failure)
        ELSE
           CPPrecondition(PRESENT(matrix_p),cp_failure_level,routineP,failure)
        END IF
     END IF

     ! *** Allocate work storage ***
     ldsab = get_memory_usage(qs_kind_set,basis_type)
     ALLOCATE (kab(ldsab,ldsab),qab(ldsab,ldsab),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,failure)
     IF (do_forces) THEN
        ALLOCATE(dab(ldsab,ldsab,3),pab(ldsab,ldsab),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,failure)
     END IF

     ! Iterate of neighbor list
     CALL neighbor_list_iterator_create(nl_iterator,sab_nl)
     DO WHILE (neighbor_list_iterate(nl_iterator)==0)
        CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,inode=inode,&
                               iatom=iatom,jatom=jatom,r=rab,cell=cell)
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

        IF(use_cell_mapping) THEN
           ic = cell_to_index(cell(1),cell(2),cell(3))
           CPPostcondition(ic > 0,cp_failure_level,routineP,failure)
        ELSE
           ic = 1
        END IF

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
        IF(dokp) THEN
           CALL cp_dbcsr_get_block_p(matrix=matrixkp_t(1,ic)%matrix,&
                row=irow,col=icol,BLOCK=k_block,found=found)
           CPPostcondition(found,cp_failure_level,routineP,failure)
        ELSE
           CALL cp_dbcsr_get_block_p(matrix=matrix_t(1)%matrix,&
                row=irow,col=icol,BLOCK=k_block,found=found)
           CPPostcondition(found,cp_failure_level,routineP,failure)
        END IF

        IF (do_forces) THEN
           NULLIFY (p_block)
           IF(dokp) THEN
              CALL cp_dbcsr_get_block_p(matrix=matrixkp_p(1,ic)%matrix,&
                   row=irow,col=icol,block=p_block,found=found)
              CPPostcondition(found,cp_failure_level,routineP,failure)
           ELSE
              CALL cp_dbcsr_get_block_p(matrix=matrix_p,row=irow,col=icol,&
                   block=p_block,found=found)
              CPPostcondition(found,cp_failure_level,routineP,failure)
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
               CALL block_add("OUT",kab,nsgfa(iset),nsgfb(jset),p_block,sgfa,sgfb,trans=trans)
               CALL decontraction(kab,pab,scon_a(:,sgfa:),ncoa,nsgfa(iset),scon_b(:,sgfb:),ncob,nsgfb(jset),&
                                  trans=trans)
               ! calculate integrals and derivatives
               CALL kinetic(la_max(iset),la_min(iset),npgfa(iset),rpgfa(:,iset),zeta(:,iset),&
                            lb_max(jset),lb_min(jset),npgfb(jset),rpgfb(:,jset),zetb(:,jset),&
                            rab,kab,dab)
               CALL force_trace(force_a,dab,pab,ncoa,ncob,3)
               force(ikind)%kinetic(:,atom_a)=force(ikind)%kinetic(:,atom_a) + ff*force_a(:)
               force(jkind)%kinetic(:,atom_b)=force(jkind)%kinetic(:,atom_b) - ff*force_a(:)
               IF( use_virial ) THEN
                  CALL virial_pair_force ( virial%pv_virial, f0, force_a, rab)
               END IF
            ELSE
               ! calclulate integrals
               CALL kinetic(la_max(iset),la_min(iset),npgfa(iset),rpgfa(:,iset),zeta(:,iset),&
                            lb_max(jset),lb_min(jset),npgfb(jset),rpgfb(:,jset),zetb(:,jset),&
                            rab,kab)
            END IF
            ! Contraction step 
            CALL contraction(kab,qab,ca=scon_a(:,sgfa:),na=ncoa,ma=nsgfa(iset),&
                                     cb=scon_b(:,sgfb:),nb=ncob,mb=nsgfb(jset),&
                                     trans=trans)
            CALL block_add("IN",qab,nsgfa(iset),nsgfb(jset),k_block,sgfa,sgfb,trans=trans)

          END DO
        END DO

     END DO
     CALL neighbor_list_iterator_release(nl_iterator)

     IF(dokp) THEN
        DO ic=1,nimg
           CALL cp_dbcsr_finalize(matrixkp_t(1,ic)%matrix)
           IF (PRESENT(eps_filter)) THEN
              CALL cp_dbcsr_filter(matrixkp_t(1,ic)%matrix, eps_filter)
           END IF
        END DO
     ELSE
        CALL cp_dbcsr_finalize(matrix_t(1)%matrix)
        IF (PRESENT(eps_filter)) THEN
           CALL cp_dbcsr_filter(matrix_t(1)%matrix, eps_filter)
        END IF
     END IF

     ! Release work storage
     DEALLOCATE (atom_of_kind,STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,failure)
     DEALLOCATE (kab,qab,STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,failure)
     IF(do_forces) THEN
        DEALLOCATE (pab,dab,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,failure)
     END IF
     DEALLOCATE (basis_set_list,STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,failure)

     CALL timestop(handle)

  END SUBROUTINE build_kinetic_matrix

END MODULE qs_kinetic

