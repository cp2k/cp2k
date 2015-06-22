!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of overlap matrix, its derivatives and forces
!> \par History
!>      JGH: removed printing routines
!>      JGH: upgraded to unique routine for overlaps
!>      JGH: Add specific routine for 'forces only'
!>           Major refactoring for new overlap routines
!>      JGH: Kpoints
!> \author Matthias Krack (03.09.2001,25.06.2003)
! *****************************************************************************
MODULE qs_overlap
  USE ai_contraction,                  ONLY: block_add,&
                                             contraction,&
                                             decontraction,&
                                             force_trace
  USE ai_overlap,                      ONLY: overlap_ab
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                             gto_basis_set_type
  USE block_p_types,                   ONLY: block_p_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_allocate_matrix_set, cp_dbcsr_create, cp_dbcsr_filter, &
       cp_dbcsr_finalize, cp_dbcsr_get_block_p, cp_dbcsr_init, &
       cp_dbcsr_p_type, cp_dbcsr_type, dbcsr_distribution_obj, &
       dbcsr_type_antisymmetric, dbcsr_type_no_symmetry, dbcsr_type_symmetric
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE kpoint_types,                    ONLY: get_kpoint_info,&
                                             kpoint_type
  USE orbital_pointers,                ONLY: indco,&
                                             ncoset
  USE orbital_symbols,                 ONLY: cgf_symbol
  USE particle_methods,                ONLY: get_particle_set
  USE particle_types,                  ONLY: particle_type
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
  USE string_utilities,                ONLY: compress,&
                                             uppercase
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_overlap'

  INTEGER, DIMENSION(1:56), PARAMETER :: ndod = (/0,1,1,1,0,0,0,0,0,0,1,1,1,&
      1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,&
      1,1,1,1,1,1,1/)

  INTERFACE create_sab_matrix
     MODULE PROCEDURE create_sab_matrix_1d, create_sab_matrix_2d
  END INTERFACE


! *** Public subroutines ***

  PUBLIC :: build_overlap_matrix, build_overlap_force, create_sab_matrix

CONTAINS

! *****************************************************************************

!> \brief   Calculation of the overlap matrix over Cartesian Gaussian functions.
!> \param   ks_env the QS env
!> \param   matrix_s The overlap matrix to be calculated (optional)
!> \param   matrixkp_s The overlap matrices to be calculated (kpoints, optional)
!> \param   matrix_name The name of the overlap matrix (i.e. for output)
!> \param   nderivative Derivative with respect to basis origin
!> \param   basis_set_id_a basis set to be used for bra in <a|b>
!> \param   basis_set_id_b basis set to be used for ket in <a|b>
!> \param   sab_nl pair list (must be consistent with basis sets!)
!> \param   calculate_forces (optional)
!> \param   matrix_p density matrix for force calculation (optional)
!> \param   matrixkp_p density matrix for force calculation with k_points (optional)
!> \param   error CP2K error reporting 
!> \date    11.03.2002
!> \par     History
!>          Enlarged functionality of this routine. Now overlap matrices based
!>          on different basis sets can be calculated, taking into account also
!>          mixed overlaps NOTE: the pointer to the overlap matrix must now be
!>          put into its corresponding env outside of this routine
!>          [Manuel Guidon]
!>          Generalized for derivatives and force calculations [JHU]
!>          Kpoints, returns overlap matrices in real space index form
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE build_overlap_matrix(ks_env, matrix_s, matrixkp_s, matrix_name,&
     nderivative, basis_set_id_a, basis_set_id_b, sab_nl, calculate_forces,&
     matrix_p, matrixkp_p, error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: matrix_s
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: matrixkp_s
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: matrix_name
    INTEGER, INTENT(IN), OPTIONAL            :: nderivative
    INTEGER, INTENT(IN)                      :: basis_set_id_a, basis_set_id_b
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_nl
    LOGICAL, INTENT(IN), OPTIONAL            :: calculate_forces
    TYPE(cp_dbcsr_type), OPTIONAL, POINTER   :: matrix_p
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: matrixkp_p
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_overlap_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, i, iatom, ic, icol, ikind, inode, &
      irow, iset, istat, jatom, jkind, jset, ldsab, maxder, maxs, n1, n2, &
      natom, ncoa, ncob, nder, nimg, nkind, nseta, nsetb, sgfa, sgfb
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
    REAL(KIND=dp)                            :: dab, f, f0, ff, rab2
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: owork, pmat
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: oint
    REAL(KIND=dp), DIMENSION(3)              :: force_a, rab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_block, rpgfa, rpgfb, &
                                                scon_a, scon_b, zeta, zetb
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(block_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: sint
    TYPE(dbcsr_distribution_obj), POINTER    :: dbcsr_dist
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list_a, &
                                                basis_set_list_b
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b
    TYPE(kpoint_type), POINTER               :: kpoints
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(virial_type), POINTER               :: virial

     NULLIFY(dft_control)
     failure = .FALSE.

     CALL timeset(routineN,handle)

     ! test for matrices (kpoints or standard gamma point)
     IF (PRESENT(matrix_s)) THEN
        dokp = .FALSE.
        use_cell_mapping = .FALSE.
     ELSEIF (PRESENT(matrixkp_s)) THEN
        dokp = .TRUE.
        CALL get_ks_env(ks_env=ks_env,kpoints=kpoints,error=error)
        CALL get_kpoint_info(kpoint=kpoints,cell_to_index=cell_to_index,error=error)
        use_cell_mapping = (SIZE(cell_to_index) > 1)
     ELSE
        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
     END IF

     NULLIFY (atomic_kind_set)
     CALL get_ks_env(ks_env,&
                     atomic_kind_set=atomic_kind_set,&
                     natom=natom,&
                     qs_kind_set=qs_kind_set,&
                     dbcsr_dist=dbcsr_dist,&
                     dft_control=dft_control,&
                     error=error)

     nimg = dft_control%nimages
     nkind = SIZE(qs_kind_set)

     ALLOCATE (atom_of_kind(natom),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     CALL get_atomic_kind_set(atomic_kind_set,atom_of_kind=atom_of_kind)

     IF (PRESENT(calculate_forces)) THEN
       do_forces = calculate_forces
     ELSE
       do_forces = .FALSE.
     END IF

     IF (PRESENT(nderivative)) THEN
       nder = nderivative
     ELSE
       nder = 0
     END IF
     maxder = ncoset(nder)

     ! check for symmetry
     CPPrecondition(SIZE(sab_nl) > 0,cp_failure_level,routineP,error,failure)
     CALL get_neighbor_list_set_p(neighbor_list_sets=sab_nl,symmetric=do_symmetric)
     IF (do_symmetric) THEN
       CPPrecondition(basis_set_id_a == basis_set_id_b,cp_failure_level,routineP,error,failure)
     END IF

     IF(dokp) THEN
        CALL cp_dbcsr_allocate_matrix_set(matrixkp_s,maxder,nimg,error=error)
        CALL create_sab_matrix(ks_env, matrixkp_s, matrix_name, basis_set_id_a, basis_set_id_b, &
                               sab_nl, do_symmetric, error)
     ELSE
        CALL cp_dbcsr_allocate_matrix_set(matrix_s,maxder,error=error)
        CALL create_sab_matrix(ks_env, matrix_s, matrix_name, basis_set_id_a, basis_set_id_b, &
                               sab_nl, do_symmetric, error)
     END IF
     maxs = maxder

     IF (do_forces) THEN
        CALL get_ks_env(ks_env=ks_env,force=force,virial=virial,error=error)
        use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)
     END IF

     ldsab = get_memory_usage(qs_kind_set,basis_set_id_a,basis_set_id_b)
     IF (do_forces) THEN
        ! we need density matrix for forces
        IF(dokp) THEN
           CPPrecondition(PRESENT(matrixkp_p),cp_failure_level,routineP,error,failure)
        ELSE
           CPPrecondition(PRESENT(matrix_p),cp_failure_level,routineP,error,failure)
        END IF
        nder = MAX(nder,1)
        ALLOCATE (pmat(ldsab,ldsab),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     END IF
     maxder = ncoset(nder)
     ALLOCATE (oint(ldsab,ldsab,maxder),owork(ldsab,ldsab),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

     ALLOCATE (sint(maxs),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     DO i=1,maxs
       NULLIFY (sint(i)%block)
     END DO

     ALLOCATE (basis_set_list_a(nkind),basis_set_list_b(nkind),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     CALL basis_set_list_setup(basis_set_list_a,basis_set_id_a,qs_kind_set)
     CALL basis_set_list_setup(basis_set_list_b,basis_set_id_b,qs_kind_set)

     CALL neighbor_list_iterator_create(nl_iterator,sab_nl)
     DO WHILE (neighbor_list_iterate(nl_iterator)==0)
        CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,inode=inode,&
                               iatom=iatom,jatom=jatom,r=rab,cell=cell)
        basis_set_a => basis_set_list_a(ikind)%gto_basis_set
        IF (.NOT.ASSOCIATED(basis_set_a)) CYCLE
        basis_set_b => basis_set_list_b(jkind)%gto_basis_set
        IF (.NOT.ASSOCIATED(basis_set_b)) CYCLE
        atom_a = atom_of_kind(iatom)
        atom_b = atom_of_kind(jatom)
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
           CPPostcondition(ic > 0,cp_failure_level,routineP,error,failure)
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
          ff = 2.0_dp
          IF (iatom==jatom) f0 = 1.0_dp
        ELSE
          irow = iatom
          icol = jatom
          f0 = 1.0_dp
          ff = 1.0_dp
        END IF
        DO i=1,maxs
           NULLIFY (sint(i)%block)
           IF(dokp) THEN
              CALL cp_dbcsr_get_block_p(matrix=matrixkp_s(i,ic)%matrix,&
                   row=irow,col=icol,BLOCK=sint(i)%block,found=found)
              CPPostcondition(found,cp_failure_level,routineP,error,failure)
           ELSE
              CALL cp_dbcsr_get_block_p(matrix=matrix_s(i)%matrix,&
                   row=irow,col=icol,BLOCK=sint(i)%block,found=found)
              CPPostcondition(found,cp_failure_level,routineP,error,failure)
           END IF
        END DO
        IF (do_forces) THEN
           NULLIFY (p_block)
           IF(dokp) THEN
              CALL cp_dbcsr_get_block_p(matrix=matrixkp_p(1,ic)%matrix,&
                   row=irow,col=icol,block=p_block,found=found)
              CPPostcondition(found,cp_failure_level,routineP,error,failure)
           ELSE
              CALL cp_dbcsr_get_block_p(matrix=matrix_p,row=irow,col=icol,&
                   block=p_block,found=found)
              CPPostcondition(found,cp_failure_level,routineP,error,failure)
           END IF
        END IF
        trans = do_symmetric .AND. (iatom > jatom)

        rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
        dab = SQRT(rab2)

        DO iset=1,nseta

          ncoa = npgfa(iset)*ncoset(la_max(iset))
          n1 = npgfa(iset)*(ncoset(la_max(iset))-ncoset(la_min(iset)-1))
          sgfa = first_sgfa(1,iset)

          DO jset=1,nsetb

            IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

            ncob = npgfb(jset)*ncoset(lb_max(jset))
            n2 = npgfb(jset)*(ncoset(lb_max(jset))-ncoset(lb_min(jset)-1))
            sgfb = first_sgfb(1,jset)

            ! calculate integrals and derivatives
            SELECT CASE (nder)
               CASE (0)
                  CALL overlap_ab(la_max(iset),la_min(iset),npgfa(iset),rpgfa(:,iset),zeta(:,iset),&
                       lb_max(jset),lb_min(jset),npgfb(jset),rpgfb(:,jset),zetb(:,jset),&
                       rab,sab=oint(:,:,1),error=error)
               CASE (1)
                  CALL overlap_ab(la_max(iset),la_min(iset),npgfa(iset),rpgfa(:,iset),zeta(:,iset),&
                       lb_max(jset),lb_min(jset),npgfb(jset),rpgfb(:,jset),zetb(:,jset),&
                       rab,sab=oint(:,:,1),dab=oint(:,:,2:4),error=error)
               CASE (2)
                  CALL overlap_ab(la_max(iset),la_min(iset),npgfa(iset),rpgfa(:,iset),zeta(:,iset),&
                       lb_max(jset),lb_min(jset),npgfb(jset),rpgfb(:,jset),zetb(:,jset),&
                       rab,sab=oint(:,:,1),dab=oint(:,:,2:4),ddab=oint(:,:,5:10),error=error)
               CASE DEFAULT
                  CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
            END SELECT
            IF (do_forces .AND. ASSOCIATED(p_block) .AND. ((iatom/=jatom) .OR. use_virial)) THEN
               ! Decontract P matrix block
               owork = 0.0_dp
               CALL block_add("OUT",owork,nsgfa(iset),nsgfb(jset),p_block,sgfa,sgfb,trans=trans,error=error)
               CALL decontraction(owork,pmat,scon_a(:,sgfa:),n1,nsgfa(iset),scon_b(:,sgfb:),n2,nsgfb(jset),&
                                  trans=trans,error=error)
               CALL force_trace(force_a,oint(:,:,2:4),pmat,n1,n2,3,error=error)
               force(ikind)%overlap(:,atom_a)=force(ikind)%overlap(:,atom_a) - ff*force_a(:)
               force(jkind)%overlap(:,atom_b)=force(jkind)%overlap(:,atom_b) + ff*force_a(:)
               IF( use_virial ) THEN
                  CALL virial_pair_force ( virial%pv_virial, -f0, force_a, rab, error)
               END IF
            END IF
            ! Contraction
            DO i=1,maxs
               f = 1.0_dp 
               IF (ndod(i) == 1 .AND. trans) f = -1.0_dp
               CALL contraction(oint(:,:,i),owork,ca=scon_a(:,sgfa:),na=n1,ma=nsgfa(iset),&
                    cb=scon_b(:,sgfb:),nb=n2,mb=nsgfb(jset),fscale=f,trans=trans,error=error)
               CALL block_add("IN",owork,nsgfa(iset),nsgfb(jset),sint(i)%block,&
                    sgfa,sgfb,trans=trans,error=error)
            END DO

          END DO
        END DO

     END DO
     CALL neighbor_list_iterator_release(nl_iterator)

     IF(dokp) THEN
        DO i = 1,maxs
           DO ic=1,nimg
              CALL cp_dbcsr_finalize(matrixkp_s(i,ic)%matrix, error=error)
              CALL cp_dbcsr_filter(matrixkp_s(i,ic)%matrix, &
                   dft_control%qs_control%eps_filter_matrix, error=error)
           ENDDO
        ENDDO
     ELSE
        DO i = 1,maxs
           CALL cp_dbcsr_finalize(matrix_s(i)%matrix, error=error)
           CALL cp_dbcsr_filter(matrix_s(i)%matrix, &
                dft_control%qs_control%eps_filter_matrix, error=error)
        ENDDO
     END IF

     ! *** Release work storage ***
     DEALLOCATE (atom_of_kind,STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     DO i=1,maxs
        NULLIFY (sint(i)%block)
     END DO
     DEALLOCATE (sint,STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     DEALLOCATE (basis_set_list_a,basis_set_list_b,STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     DEALLOCATE (oint,owork,STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     IF(do_forces) THEN
        DEALLOCATE (pmat,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     END IF

     CALL timestop(handle)

  END SUBROUTINE build_overlap_matrix

! *****************************************************************************

!> \brief   Calculation of the force contribution from an overlap matrix 
!>          over Cartesian Gaussian functions.
!> \param   ks_env the QS environment
!> \param   force holds the calcuated force Tr(P dS/dR)
!> \param   basis_set_id_a basis set to be used for bra in <a|b>
!> \param   basis_set_id_b basis set to be used for ket in <a|b>
!> \param   sab_nl pair list (must be consistent with basis sets!)
!> \param   matrix_p density matrix for force calculation 
!> \param   error for error handling
!> \date    11.03.2002
!> \par     History
!>          Enlarged functionality of this routine. Now overlap matrices based
!>          on different basis sets can be calculated, taking into account also
!>          mixed overlaps NOTE: the pointer to the overlap matrix must now be
!>          put into its corresponding env outside of this routine
!>          [Manuel Guidon]
!>          Generalized for derivatives and force calculations [JHU]
!>          This special version is only for forces [07.07.2014, JGH]
!> \author  MK/JGH
!> \version 1.0
! *****************************************************************************
  SUBROUTINE build_overlap_force(ks_env, force, basis_set_id_a, basis_set_id_b, &
                                 sab_nl, matrix_p, error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: force
    INTEGER, INTENT(IN)                      :: basis_set_id_a, basis_set_id_b
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_nl
    TYPE(cp_dbcsr_type), POINTER             :: matrix_p
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_overlap_force', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iatom, icol, ikind, inode, irow, iset, istat, jatom, &
      jkind, jset, last_jatom, ldsab, n1, n2, ncoa, ncob, nder, nkind, nseta, &
      nsetb, sgfa, sgfb
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: do_symmetric, failure, found, &
                                                new_atom_b, trans, use_virial
    REAL(KIND=dp)                            :: dab, f0, ff, rab2
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pab, sab
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: drab
    REAL(KIND=dp), DIMENSION(3)              :: force_a, rab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_block, rpgfa, rpgfb, &
                                                scon_a, scon_b, zeta, zetb
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list_a, &
                                                basis_set_list_b
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(virial_type), POINTER               :: virial

     failure = .FALSE.

     CALL timeset(routineN,handle)

     NULLIFY (qs_kind_set)
     CALL get_ks_env(ks_env=ks_env,qs_kind_set=qs_kind_set,error=error)

     nkind = SIZE(qs_kind_set)
     nder = 1

     ! check for symmetry
     CPPrecondition(SIZE(sab_nl) > 0,cp_failure_level,routineP,error,failure)
     CALL get_neighbor_list_set_p(neighbor_list_sets=sab_nl,symmetric=do_symmetric)

     CALL get_ks_env(ks_env=ks_env,virial=virial,error=error)
     use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)

     ! *** Allocate work storage ***
     ldsab = get_memory_usage(qs_kind_set,basis_set_id_a,basis_set_id_b)
     ALLOCATE (sab(ldsab,ldsab),pab(ldsab,ldsab),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     ALLOCATE (drab(ldsab,ldsab,3),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

     ! set up basis sets
     ALLOCATE (basis_set_list_a(nkind),basis_set_list_b(nkind),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     CALL basis_set_list_setup(basis_set_list_a,basis_set_id_a,qs_kind_set)
     CALL basis_set_list_setup(basis_set_list_b,basis_set_id_b,qs_kind_set)

     ! Loop over neighbor list
     CALL neighbor_list_iterator_create(nl_iterator,sab_nl)
     DO WHILE (neighbor_list_iterate(nl_iterator)==0)
        CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,inode=inode,&
                               iatom=iatom,jatom=jatom,r=rab)
        basis_set_a => basis_set_list_a(ikind)%gto_basis_set
        IF (.NOT.ASSOCIATED(basis_set_a)) CYCLE
        basis_set_b => basis_set_list_b(jkind)%gto_basis_set
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
          NULLIFY (p_block)
          CALL cp_dbcsr_get_block_p(matrix=matrix_p,row=irow,col=icol,&
                  block=p_block,found=found)
        END IF
        trans = do_symmetric .AND. (iatom > jatom)

        rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
        dab = SQRT(rab2)

        DO iset=1,nseta

          ncoa = npgfa(iset)*ncoset(la_max(iset))
          n1 = npgfa(iset)*(ncoset(la_max(iset))-ncoset(la_min(iset)-1))
          sgfa = first_sgfa(1,iset)

          DO jset=1,nsetb

            IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

            ncob = npgfb(jset)*ncoset(lb_max(jset))
            n2 = npgfb(jset)*(ncoset(lb_max(jset))-ncoset(lb_min(jset)-1))
            sgfb = first_sgfb(1,jset)

            IF (ASSOCIATED(p_block) .AND. ((iatom/=jatom) .OR. use_virial)) THEN
               ! Decontract P matrix block
               sab = 0.0_dp
               CALL block_add("OUT",sab,nsgfa(iset),nsgfb(jset),p_block,sgfa,sgfb,trans=trans,error=error)
               CALL decontraction(sab,pab,scon_a(:,sgfa:),n1,nsgfa(iset),scon_b(:,sgfb:),n2,nsgfb(jset),&
                                  trans=trans,error=error)
               ! calculate integrals and derivatives
               CALL overlap_ab(la_max(iset),la_min(iset),npgfa(iset),rpgfa(:,iset),zeta(:,iset),&
                               lb_max(jset),lb_min(jset),npgfb(jset),rpgfb(:,jset),zetb(:,jset),&
                               rab,dab=drab,error=error)
               CALL force_trace(force_a,drab,pab,n1,n2,3,error=error)
               force(1:3,iatom) = force(1:3,iatom) - ff*force_a(1:3)
               force(1:3,jatom) = force(1:3,jatom) + ff*force_a(1:3)
               IF( use_virial ) THEN
                  CALL virial_pair_force ( virial%pv_virial, -f0, force_a, rab, error)
               END IF
            END IF

          END DO
        END DO

     END DO
     CALL neighbor_list_iterator_release(nl_iterator)

     ! *** Release work storage ***
     DEALLOCATE (sab,pab,drab,STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     DEALLOCATE (basis_set_list_a,basis_set_list_b,STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

     CALL timestop(handle)

  END SUBROUTINE build_overlap_force

! *****************************************************************************
!> \brief Setup the structure of a sparse matrix based on the overlap 
!>        neighbor list
!> \param ks_env         The QS environment
!> \param matrix_s       Matrices to be constructed
!> \param matrix_name    Matrix base name
!> \param basis_set_id_a Basis set used for <a|
!> \param basis_set_id_b Basis set used for |b>
!> \param sab_nl         Overlap neighbor list
!> \param symmetric      Is symmetry used in the neighbor list?
!> \param error          CP2K error reporting
! *****************************************************************************
  SUBROUTINE create_sab_matrix_1d(ks_env, matrix_s, matrix_name,&
     basis_set_id_a, basis_set_id_b, sab_nl, symmetric, error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: matrix_name
    INTEGER, INTENT(IN)                      :: basis_set_id_a, basis_set_id_b
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_nl
    LOGICAL, INTENT(IN)                      :: symmetric
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_sab_matrix_1d', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=12)                        :: cgfsym
    CHARACTER(LEN=32)                        :: symmetry_string
    CHARACTER(LEN=default_string_length)     :: mname, name
    INTEGER                                  :: i, istat, maxs, natom
    INTEGER, DIMENSION(:), POINTER           :: col_blk_sizes, row_blk_sizes
    LOGICAL                                  :: failure
    TYPE(dbcsr_distribution_obj), POINTER    :: dbcsr_dist
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

     failure = .FALSE.

     CALL get_ks_env(ks_env=ks_env,particle_set=particle_set,&
                     qs_kind_set=qs_kind_set, dbcsr_dist=dbcsr_dist,error=error)

     natom = SIZE(particle_set)

     IF (PRESENT(matrix_name)) THEN
       mname = matrix_name
     ELSE
       mname = "DUMMY"
     END IF

     maxs = SIZE(matrix_s)

     ALLOCATE (row_blk_sizes(natom),col_blk_sizes(natom), STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

     CALL get_particle_set(particle_set, qs_kind_set, nsgf=row_blk_sizes,&
                           basis_set_id=basis_set_id_a,error=error)
     CALL get_particle_set(particle_set, qs_kind_set, nsgf=col_blk_sizes,&
                           basis_set_id=basis_set_id_b,error=error)

     ! prepare for allocation
     IF (symmetric) THEN
        symmetry_string = dbcsr_type_symmetric
     ELSE
        symmetry_string = dbcsr_type_no_symmetry
     END IF

     DO i=1,maxs
       IF (symmetric) THEN
          IF(ndod(i)==1) THEN
             ! odd derivatives are anti-symmetric
             symmetry_string = dbcsr_type_antisymmetric
          ELSE
             symmetry_string = dbcsr_type_symmetric
          END IF
       ELSE
          symmetry_string = dbcsr_type_no_symmetry
       END IF
       cgfsym = cgf_symbol(1,indco(1:3,i))
       IF(i==1) THEN
          name = mname
       ELSE
          name = TRIM(cgfsym(4:))//" DERIVATIVE OF THE "//TRIM(mname)//&
              " W.R.T. THE NUCLEAR COORDINATES"
       END IF
       CALL compress(name)
       CALL uppercase(name)
       ALLOCATE(matrix_s(i)%matrix, STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       CALL cp_dbcsr_init(matrix_s(i)%matrix,error=error)
       CALL cp_dbcsr_create(matrix=matrix_s(i)%matrix, &
            name=TRIM(name), &
            dist=dbcsr_dist, matrix_type=symmetry_string,&
            row_blk_size=row_blk_sizes, col_blk_size=col_blk_sizes, &
            nblks=0, nze=0, error=error)
       CALL cp_dbcsr_alloc_block_from_nbl(matrix_s(i)%matrix,sab_nl,error)
     END DO
     
     DEALLOCATE(row_blk_sizes,col_blk_sizes)

  END SUBROUTINE create_sab_matrix_1d
! *****************************************************************************

! *****************************************************************************
!> \brief Setup the structure of a sparse matrix based on the overlap 
!>        neighbor list, 2d version
!> \param ks_env         The QS environment
!> \param matrix_s       Matrices to be constructed
!> \param matrix_name    Matrix base name
!> \param basis_set_id_a Basis set used for <a|
!> \param basis_set_id_b Basis set used for |b>
!> \param sab_nl         Overlap neighbor list
!> \param symmetric      Is symmetry used in the neighbor list?
!> \param error          CP2K error reporting
! *****************************************************************************
  SUBROUTINE create_sab_matrix_2d(ks_env, matrix_s, matrix_name,&
     basis_set_id_a, basis_set_id_b, sab_nl, symmetric, error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: matrix_s
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: matrix_name
    INTEGER, INTENT(IN)                      :: basis_set_id_a, basis_set_id_b
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_nl
    LOGICAL, INTENT(IN)                      :: symmetric
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_sab_matrix_2d', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=12)                        :: cgfsym
    CHARACTER(LEN=32)                        :: symmetry_string
    CHARACTER(LEN=default_string_length)     :: mname, name
    INTEGER                                  :: i1, i2, istat, maxs1, maxs2, &
                                                natom
    INTEGER, DIMENSION(:), POINTER           :: col_blk_sizes, row_blk_sizes
    LOGICAL                                  :: failure
    TYPE(dbcsr_distribution_obj), POINTER    :: dbcsr_dist
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    failure = .FALSE.

    CALL get_ks_env(ks_env=ks_env,particle_set=particle_set,&
                    qs_kind_set=qs_kind_set, dbcsr_dist=dbcsr_dist,error=error)

    natom = SIZE(particle_set)

    IF (PRESENT(matrix_name)) THEN
      mname = matrix_name
    ELSE
      mname = "DUMMY"
    END IF

    maxs1 = SIZE(matrix_s,1)
    maxs2 = SIZE(matrix_s,2)

    ALLOCATE (row_blk_sizes(natom),col_blk_sizes(natom), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL get_particle_set(particle_set, qs_kind_set, nsgf=row_blk_sizes,&
                          basis_set_id=basis_set_id_a,error=error)
    CALL get_particle_set(particle_set, qs_kind_set, nsgf=col_blk_sizes,&
                          basis_set_id=basis_set_id_b,error=error)

    ! prepare for allocation
    IF (symmetric) THEN
       symmetry_string = dbcsr_type_symmetric
    ELSE
       symmetry_string = dbcsr_type_no_symmetry
    END IF

    DO i2=1,maxs2
       DO i1=1,maxs1
          IF (symmetric) THEN
             IF(ndod(i1)==1) THEN
                ! odd derivatives are anti-symmetric
                symmetry_string = dbcsr_type_antisymmetric
             ELSE
                symmetry_string = dbcsr_type_symmetric
             END IF
          ELSE
             symmetry_string = dbcsr_type_no_symmetry
          END IF
          cgfsym = cgf_symbol(1,indco(1:3,i1))
          IF(i1==1) THEN
             name = mname
          ELSE
             name = TRIM(cgfsym(4:))//" DERIVATIVE OF THE "//TRIM(mname)//&
                 " W.R.T. THE NUCLEAR COORDINATES"
          END IF
          CALL compress(name)
          CALL uppercase(name)
          ALLOCATE(matrix_s(i1,i2)%matrix, STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          CALL cp_dbcsr_init(matrix_s(i1,i2)%matrix,error=error)
          CALL cp_dbcsr_create(matrix=matrix_s(i1,i2)%matrix, &
               name=TRIM(name), &
               dist=dbcsr_dist, matrix_type=symmetry_string,&
               row_blk_size=row_blk_sizes, col_blk_size=col_blk_sizes, &
               nblks=0, nze=0, error=error)
          CALL cp_dbcsr_alloc_block_from_nbl(matrix_s(i1,i2)%matrix,sab_nl,error)
       END DO
    END DO

    DEALLOCATE(row_blk_sizes,col_blk_sizes)

  END SUBROUTINE create_sab_matrix_2d
! *****************************************************************************

END MODULE qs_overlap


