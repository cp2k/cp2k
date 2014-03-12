!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of overlap matrix, its derivatives and forces
!> \par History
!>      JGH: removed printing routines
!>      JGH: upgraded to unique routine for overlaps
!> \author Matthias Krack (03.09.2001,25.06.2003)
! *****************************************************************************
MODULE qs_overlap

  USE ai_overlap_new,                  ONLY: overlap
  USE array_types,                     ONLY: array_i1d_obj,&
                                             array_new,&
                                             array_nullify,&
                                             array_release
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set_info,&
                                             gto_basis_set_p_type,&
                                             gto_basis_set_type
  USE block_p_types,                   ONLY: block_p_type
  USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_create,&
                                             cp_dbcsr_filter,&
                                             cp_dbcsr_finalize,&
                                             cp_dbcsr_get_block_p,&
                                             cp_dbcsr_init
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_allocate_matrix_set
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dbcsr_types,                     ONLY: dbcsr_distribution_obj,&
                                             dbcsr_type_antisymmetric,&
                                             dbcsr_type_no_symmetry,&
                                             dbcsr_type_symmetric
  USE input_constants,                 ONLY: use_aux_fit_basis_set,&
                                             use_orb_basis_set
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE orbital_pointers,                ONLY: indco,&
                                             init_orbital_pointers,&
                                             ncoset
  USE orbital_symbols,                 ONLY: cgf_symbol
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             get_neighbor_list_set,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type,&
                                             neighbor_list_set_type
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

! *** Public subroutines ***

  PUBLIC :: build_overlap_matrix, build_overlap_matrix_simple 

CONTAINS

! *****************************************************************************
!> \brief   Calculation of the overlap matrix over Cartesian Gaussian functions.
!> \author  MK
!> \par     History
!>          Enlarged functionality of this routine. Now overlap matrices based
!>          on different basis sets can be calculated, taking into account also
!>          mixed overlaps NOTE: the pointer to the overlap matrix must now be
!>          put into its corresponding env outside of this routine
!>          [Manuel Guidon]
!>          Generalized for derivatives and force calculations [JHU]
!> \param   qs_env the QS env
!> \param   para_env the parallel env
!> \param   nderivative Derivative with respect to basis origin
!> \param   matrix_s The overlap matrix to be calculated
!> \param   matrix_name The name of the overlap matrix (i.e. for output)
!> \param   basis_set_id_a basis set to be used for bra in <a|b>
!> \param   basis_set_id_b basis set to be used for ket in <a|b>
!> \param   sab_nl pair list (must be consistent with basis sets!)
!> \param   calculate_forces (optional)
!> \param   admm (optional), use force_admm
!> \param   matrix_p density matrix for force calculation (optional)
!> \param   error for error handling
!> \date    11.03.2002
!> \version 1.0
! *****************************************************************************

  SUBROUTINE build_overlap_matrix(qs_env,nderivative,matrix_s,matrix_name,&
     basis_set_id_a, basis_set_id_b, sab_nl, calculate_forces, admm, matrix_p, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN), OPTIONAL            :: nderivative
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: matrix_s
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: matrix_name
    INTEGER, INTENT(IN), OPTIONAL            :: basis_set_id_a, basis_set_id_b
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_nl
    LOGICAL, INTENT(IN), OPTIONAL            :: calculate_forces, admm
    TYPE(cp_dbcsr_type), OPTIONAL, POINTER   :: matrix_p
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_overlap_matrix', &
      routineP = moduleN//':'//routineN
    INTEGER, DIMENSION(1:56), PARAMETER :: ndod = (/0,1,1,1,0,0,0,0,0,0,1,1,1,&
      1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,&
      1,1,1,1,1,1,1/)

    CHARACTER(LEN=12)                        :: cgfsym
    CHARACTER(LEN=32)                        :: symmetry_string
    CHARACTER(LEN=default_string_length)     :: mname, name
    INTEGER :: atom_a, atom_b, handle, i, iatom, iblock, icol, ikind, inode, &
      irow, iset, istat, jatom, jkind, jset, last_jatom, ldai, ldsab, &
      max_maxco, max_maxlgto, max_maxsgf, maxblock, maxco, maxdco, maxdcoa, &
      maxdcob, maxder, maxlgto, maxsgf, natom, ncoa, ncob, nder, &
      neighbor_list_id, nkind, nseta, nsetb, sgfa, sgfb
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(:), POINTER           :: cbs, la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb, rbs
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL :: allocate_matrix, do_forces, do_symmetric, failure, found, &
      is_admm, new_atom_b, return_s_derivatives, use_virial
    REAL(KIND=dp)                            :: dab, f, f0, rab2
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pab, sab, work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: ai_work
    REAL(KIND=dp), DIMENSION(3)              :: force_a, rab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_block, rpgfa, rpgfb, &
                                                sphi_a, sphi_b, zeta, zetb
    TYPE(array_i1d_obj)                      :: col_blk_sizes, row_blk_sizes
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(block_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: sint
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_distribution_obj), POINTER    :: dbcsr_dist
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list_a, &
                                                basis_set_list_b
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(virial_type), POINTER               :: virial

     failure = .FALSE.

     CALL timeset(routineN,handle)

     CPPrecondition(PRESENT(basis_set_id_a),cp_failure_level,routineP,error,failure)
     CPPrecondition(PRESENT(basis_set_id_b),cp_failure_level,routineP,error,failure)

     NULLIFY (atomic_kind_set, particle_set, para_env)
     CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,particle_set=particle_set,&
                     neighbor_list_id=neighbor_list_id,dbcsr_dist=dbcsr_dist,&
                     para_env=para_env,error=error)

     nkind = SIZE(atomic_kind_set)
     natom = SIZE(particle_set)

     ALLOCATE (atom_of_kind(natom),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind)

     IF (PRESENT(matrix_name)) THEN
       mname = matrix_name
     ELSE
       mname = "DUMMY"
     END IF

     IF (PRESENT(calculate_forces)) THEN
       do_forces = calculate_forces
     ELSE
       do_forces = .FALSE.
     END IF

     IF (PRESENT(admm)) THEN
       is_admm = admm
     ELSE
       is_admm = .FALSE.
     END IF

     IF (PRESENT(nderivative)) THEN
       IF ( nderivative < 0 ) THEN
         ! do not allocate any matrices but calculate forces
         CPPrecondition(do_forces,cp_failure_level,routineP,error,failure)
         nder = 0
         allocate_matrix = .FALSE.
       ELSE
         ! allocate and return maxder matrices
         nder = nderivative
         allocate_matrix = .TRUE.
       END IF
     ELSE
       ! default is to allocate and return the overlap matrix
       nder = 0
       allocate_matrix = .TRUE.
     END IF
     maxder = ncoset(nder)

     ! check for symmetry
     CPPrecondition(SIZE(sab_nl) > 0,cp_failure_level,routineP,error,failure)
     DO i = 1, SIZE(sab_nl)
        neighbor_list_set => sab_nl(i)%neighbor_list_set
        IF (ASSOCIATED(neighbor_list_set)) EXIT
     END DO
     CALL get_neighbor_list_set(neighbor_list_set=neighbor_list_set,symmetric=do_symmetric)
     IF (do_symmetric) THEN
       CPPrecondition(basis_set_id_a == basis_set_id_b,cp_failure_level,routineP,error,failure)
     END IF

     IF (allocate_matrix) THEN
        CPPrecondition(PRESENT(matrix_s),cp_failure_level,routineP,error,failure)
        ALLOCATE (rbs(natom),cbs(natom), STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

        CALL get_particle_set(particle_set=particle_set,nsgf=rbs,&
                              basis_set_id=basis_set_id_a,error=error)
        CALL get_particle_set(particle_set=particle_set,nsgf=cbs,&
                              basis_set_id=basis_set_id_b,error=error)

        CALL cp_dbcsr_allocate_matrix_set(matrix_s,maxder,error=error)
        ! prepare for allocation
        CALL array_nullify (row_blk_sizes)
        CALL array_nullify (col_blk_sizes)
        CALL array_new (row_blk_sizes, rbs, gift=.TRUE.)
        CALL array_new (col_blk_sizes, cbs, gift=.TRUE.)
        IF (do_symmetric) THEN
           symmetry_string = dbcsr_type_symmetric
        ELSE
           symmetry_string = dbcsr_type_no_symmetry
        END IF

        ALLOCATE(matrix_s(1)%matrix, STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        CALL cp_dbcsr_init(matrix_s(1)%matrix,error=error)
        CALL cp_dbcsr_create(matrix=matrix_s(1)%matrix, &
             name=TRIM(mname), &
             dist=dbcsr_dist, matrix_type=symmetry_string,&
             row_blk_size=row_blk_sizes, col_blk_size=col_blk_sizes, &
             nblks=0, nze=0, error=error)
        CALL cp_dbcsr_alloc_block_from_nbl(matrix_s(1)%matrix,sab_nl,error)

        DO i=2,maxder
          IF (do_symmetric) THEN
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
          name = TRIM(cgfsym(4:))//" DERIVATIVE OF THE "//TRIM(mname)//&
                 " W.R.T. THE NUCLEAR COORDINATES"
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

        CALL array_release (row_blk_sizes)
        CALL array_release (col_blk_sizes)
     END IF

     IF (do_forces) THEN
        CALL get_qs_env(qs_env=qs_env,force=force,virial=virial,error=error)
        use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)
     END IF

     ! *** Allocate work storage ***
     max_maxco = 0
     max_maxlgto = 0
     max_maxsgf = 0
     CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                              maxco=maxco,&
                              maxlgto=maxlgto,&
                              maxsgf=maxsgf,&
                              basis_set_id=basis_set_id_a)
     max_maxco = MAX(maxco,max_maxco)
     max_maxlgto = MAX(maxlgto,max_maxlgto)
     max_maxsgf = MAX(maxsgf,max_maxsgf)

     CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                              maxco=maxco,&
                              maxlgto=maxlgto,&
                              maxsgf=maxsgf,&
                              basis_set_id=basis_set_id_b)
     max_maxco = MAX(maxco,max_maxco)
     max_maxlgto = MAX(maxlgto,max_maxlgto)
     max_maxsgf = MAX(maxsgf,max_maxsgf)
     ldsab = MAX(max_maxco,max_maxsgf)

     CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                              maxco=maxdcoa,&
                              maxder=MAX(1,nder),&
                              basis_set_id=basis_set_id_a)
     CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                              maxco=maxdcob,&
                              maxder=MAX(1,nder),&
                              basis_set_id=basis_set_id_b)

     IF (do_forces) THEN
        ! we need density matrix for forces
        CPPrecondition(PRESENT(matrix_p),cp_failure_level,routineP,error,failure)
        nder = MAX(nder,1)
     END IF
     return_s_derivatives = (nder > 0)
     maxder = ncoset(nder)
     maxblock = MAX(maxder,4)
     maxdco = MAX(maxdcoa, maxdcob)
     ldai = MAX(ncoset(max_maxlgto+nder+1),ldsab)
     CALL init_orbital_pointers(ldai)

     ALLOCATE (ai_work(ldai,ldai,maxder),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     ALLOCATE (sab(ldsab,ldsab*maxblock),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     sab=0._dp
     ALLOCATE (work(max_maxco,max_maxsgf),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     work=0._dp
     ALLOCATE (sint(maxder),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     DO i=1,maxder
       NULLIFY (sint(i)%block)
     END DO
     IF (do_forces) THEN
        ALLOCATE(pab(max_maxco,max_maxco),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     END IF
     NULLIFY(p_block)

     ALLOCATE (basis_set_list_a(nkind),basis_set_list_b(nkind),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)

       SELECT CASE (basis_set_id_a)
       CASE (use_orb_basis_set)
         CALL get_atomic_kind(atomic_kind=atomic_kind,orb_basis_set=basis_set_a)
       CASE (use_aux_fit_basis_set)
         CALL get_atomic_kind(atomic_kind=atomic_kind,aux_fit_basis_set=basis_set_a)
       END SELECT
       IF (ASSOCIATED(basis_set_a)) THEN
         basis_set_list_a(ikind)%gto_basis_set => basis_set_a
       ELSE
         NULLIFY(basis_set_list_b(ikind)%gto_basis_set)
       END IF

       SELECT CASE (basis_set_id_b)
       CASE (use_orb_basis_set)
         CALL get_atomic_kind(atomic_kind=atomic_kind,orb_basis_set=basis_set_b)
       CASE (use_aux_fit_basis_set)
         CALL get_atomic_kind(atomic_kind=atomic_kind,aux_fit_basis_set=basis_set_b)
       END SELECT
       IF (ASSOCIATED(basis_set_b)) THEN
         basis_set_list_b(ikind)%gto_basis_set => basis_set_b
       ELSE
         NULLIFY(basis_set_list_b(ikind)%gto_basis_set)
       END IF

     END DO
     CALL neighbor_list_iterator_create(nl_iterator,sab_nl)
     DO WHILE (neighbor_list_iterate(nl_iterator)==0)
        CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,inode=inode,&
                               iatom=iatom,jatom=jatom,r=rab)
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
        sphi_a       =>  basis_set_a%sphi
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
        sphi_b       =>  basis_set_b%sphi
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
            IF (iatom==jatom) THEN
               f0 = 1.0_dp
            ELSE
               f0 = 2.0_dp
            END IF
          ELSE
            irow = iatom
            icol = jatom
            f0 = 1.0_dp
          END IF
          IF (allocate_matrix) THEN
             DO i=1,SIZE(matrix_s)
               NULLIFY (sint(i)%block)
               CALL cp_dbcsr_get_block_p(matrix=matrix_s(i)%matrix,&
                    row=irow,col=icol,BLOCK=sint(i)%block,found=found)
                    CPPostcondition(found,cp_failure_level,routineP,error,failure)
             END DO
          END IF
          IF (do_forces) THEN
             NULLIFY (p_block)
             CALL cp_dbcsr_get_block_p(matrix=matrix_p,row=irow,col=icol,&
                  block=p_block,found=found)
          END IF
        END IF

        rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
        dab = SQRT(rab2)

        DO iset=1,nseta

          ncoa = npgfa(iset)*ncoset(la_max(iset))
          sgfa = first_sgfa(1,iset)

          DO jset=1,nsetb

            IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

            ncob = npgfb(jset)*ncoset(lb_max(jset))
            sgfb = first_sgfb(1,jset)

            ! *** Calculate the primitive overlap integrals ***
            IF (do_forces .AND. ASSOCIATED(p_block) .AND. ((iatom/=jatom) .OR. use_virial)) THEN
               ! *** Decontract P matrix block ***
               IF( do_symmetric ) THEN
                 IF (iatom <= jatom) THEN
                   CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                               p_block(sgfa,sgfb),SIZE(p_block,1),0.0_dp,work(1,1),SIZE(work,1))
                 ELSE
                   CALL dgemm("N","T",ncoa,nsgfb(jset),nsgfa(iset),1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                              p_block(sgfb,sgfa),SIZE(p_block,1),0.0_dp,work(1,1),SIZE(work,1))
                 END IF
               ELSE
                 CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                            p_block(sgfa,sgfb),SIZE(p_block,1),0.0_dp,work(1,1),SIZE(work,1))
               END IF
               CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),1.0_dp,work(1,1),SIZE(work,1),&
                          sphi_b(1,sgfb),SIZE(sphi_b,1),0.0_dp,pab(1,1),SIZE(pab,1))
               ! calculate integrals
               CALL overlap(la_max(iset),la_min(iset),npgfa(iset),&
                            rpgfa(:,iset),zeta(:,iset),&
                            lb_max(jset),lb_min(jset),npgfb(jset),&
                            rpgfb(:,jset),zetb(:,jset),&
                            rab,dab,sab,nder,return_s_derivatives,&
                            ai_work,ldai,pab=pab,force_a=force_a)
               IF(do_symmetric) THEN
                 IF (is_admm) THEN
                    force(ikind)%overlap_admm(:,atom_a)=force(ikind)%overlap_admm(:,atom_a) - 2.0_dp*force_a(:)
                    force(jkind)%overlap_admm(:,atom_b)=force(jkind)%overlap_admm(:,atom_b) + 2.0_dp*force_a(:)
                 ELSE
                    force(ikind)%overlap(:,atom_a)=force(ikind)%overlap(:,atom_a) - 2.0_dp*force_a(:)
                    force(jkind)%overlap(:,atom_b)=force(jkind)%overlap(:,atom_b) + 2.0_dp*force_a(:)
                 ENDIF
                 IF( use_virial ) THEN
                    CALL virial_pair_force ( virial%pv_virial, -f0, force_a, rab, error)
                 END IF
               ELSE
                 IF (is_admm) THEN
                    force(ikind)%overlap_admm(:,atom_a)=force(ikind)%overlap_admm(:,atom_a) - force_a(:)
                    force(jkind)%overlap_admm(:,atom_b)=force(jkind)%overlap_admm(:,atom_b) + force_a(:)
                 ELSE
                    force(ikind)%overlap(:,atom_a)=force(ikind)%overlap(:,atom_a) - force_a(:)
                    force(jkind)%overlap(:,atom_b)=force(jkind)%overlap(:,atom_b) + force_a(:)
                 ENDIF
                 IF( use_virial ) THEN
                    CALL virial_pair_force ( virial%pv_virial, -1.0_dp, force_a, rab, error)
                 END IF
               END IF
            ELSE
               CALL overlap(la_max(iset),la_min(iset),npgfa(iset),&
                            rpgfa(:,iset),zeta(:,iset),&
                            lb_max(jset),lb_min(jset),npgfb(jset),&
                            rpgfb(:,jset),zetb(:,jset),&
                            rab,dab,sab,nder,return_s_derivatives,&
                            ai_work,ldai)
            END IF

            IF(allocate_matrix) THEN
               ! *** Contraction step (overlap matrix and its derivatives) ***
               DO i=1,SIZE(matrix_s)
                  iblock = (i - 1)*SIZE(sab,1) + 1
                  CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,1.0_dp,sab(1,iblock),SIZE(sab,1),&
                       sphi_b(1,sgfb),SIZE(sphi_b,1),0.0_dp,work(1,1),SIZE(work,1))
                  f = 1.0_dp
                  IF ( do_symmetric .AND. (iatom > jatom) ) THEN
                     ! *** The odd derivatives are anti-symmetric ***
                     IF (ndod(i) == 1) f = -1.0_dp
                     CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,f,work(1,1),SIZE(work,1),&
                          sphi_a(1,sgfa),SIZE(sphi_a,1),1.0_dp,sint(i)%block(sgfb,sgfa),SIZE(sint(i)%block,1))
                  ELSE
                     CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,f,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                          work(1,1),SIZE(work,1),1.0_dp,sint(i)%block(sgfa,sgfb),SIZE(sint(i)%block,1))
                   END IF
                END DO
             END IF

          END DO
        END DO

     END DO
     CALL neighbor_list_iterator_release(nl_iterator)

     IF(allocate_matrix) THEN
        DO i = 1,SIZE(matrix_s)
           CALL cp_dbcsr_finalize(matrix_s(i)%matrix, error=error)
           CALL cp_dbcsr_filter(matrix_s(i)%matrix, &
                qs_env%dft_control%qs_control%eps_filter_matrix, error=error)
        ENDDO
     END IF

     ! *** Release work storage ***
     DEALLOCATE (atom_of_kind,STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     DEALLOCATE (ai_work,sab,work,STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     DO i=1,maxder
        NULLIFY (sint(i)%block)
     END DO
     DEALLOCATE (sint,STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     IF(do_forces) THEN
        DEALLOCATE (pab,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     END IF
     DEALLOCATE (basis_set_list_a,basis_set_list_b,STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

     CALL timestop(handle)

  END SUBROUTINE build_overlap_matrix

! *****************************************************************************

  SUBROUTINE build_overlap_matrix_simple(qs_env,basis_set_list_a,basis_set_list_b,&
     nderivative,matrix_s,matrix_name,sab_nl,eps_filter,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list_a, &
                                                basis_set_list_b
    INTEGER, INTENT(IN), OPTIONAL            :: nderivative
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: matrix_s
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: matrix_name
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_nl
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: eps_filter
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_overlap_matrix_simple', &
      routineP = moduleN//':'//routineN
    INTEGER, DIMENSION(1:56), PARAMETER :: ndod = (/0,1,1,1,0,0,0,0,0,0,1,1,1,&
      1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,&
      1,1,1,1,1,1,1/)

    CHARACTER(LEN=32)                        :: symmetry_string
    CHARACTER(LEN=default_string_length)     :: mname, name
    INTEGER :: handle, i, iatom, iblock, icol, ikind, inode, irow, iset, &
      istat, jatom, jkind, jset, last_jatom, ldai, ldsab, max_maxco, &
      max_maxlgto, max_maxsgf, maxblock, maxco, maxdco, maxdcoa, maxdcob, &
      maxder, maxlgto, maxsgf, natom, ncoa, ncob, nder, nseta, nsetb, sgfa, &
      sgfb
    INTEGER, DIMENSION(:), POINTER           :: cbs, la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb, rbs
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: do_symmetric, failure, found, &
                                                new_atom_b, &
                                                return_s_derivatives
    REAL(KIND=dp)                            :: dab, f, f0, rab2
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: sab, work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: ai_work
    REAL(KIND=dp), DIMENSION(3)              :: rab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                zeta, zetb
    TYPE(array_i1d_obj)                      :: col_blk_sizes, row_blk_sizes
    TYPE(block_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: sint
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_distribution_obj), POINTER    :: dbcsr_dist
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

     failure = .FALSE.

     CALL timeset(routineN,handle)

     IF (PRESENT(matrix_name)) THEN
        mname = matrix_name
     ELSE
        mname = "DUMMY"
     END IF

     IF (PRESENT(nderivative)) THEN
        nder = nderivative
     ELSE
        nder = 0
     END IF
     maxder = ncoset(nder)

     NULLIFY (particle_set, para_env)
     CALL get_qs_env(qs_env=qs_env,particle_set=particle_set,&
                     dbcsr_dist=dbcsr_dist,para_env=para_env,error=error)

     ! check for symmetry
     CPPrecondition(SIZE(sab_nl) > 0,cp_failure_level,routineP,error,failure)
     DO i = 1, SIZE(sab_nl)
        neighbor_list_set => sab_nl(i)%neighbor_list_set
        IF (ASSOCIATED(neighbor_list_set)) EXIT
     END DO
     CALL get_neighbor_list_set(neighbor_list_set=neighbor_list_set,symmetric=do_symmetric)

     natom = SIZE(particle_set)
     ALLOCATE (rbs(natom),cbs(natom), STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

     CALL get_particle_set(particle_set=particle_set,nsgf=rbs,basis=basis_set_list_a,error=error)
     CALL get_particle_set(particle_set=particle_set,nsgf=cbs,basis=basis_set_list_b,error=error)

     CALL cp_dbcsr_allocate_matrix_set(matrix_s,maxder,error=error)
     ! prepare for allocation
     CALL array_nullify (row_blk_sizes)
     CALL array_nullify (col_blk_sizes)
     CALL array_new (row_blk_sizes, rbs, gift=.TRUE.)
     CALL array_new (col_blk_sizes, cbs, gift=.TRUE.)

     DO i=1,maxder
       IF (do_symmetric) THEN
          IF(ndod(i)==1) THEN
             ! odd derivatives are anti-symmetric
             symmetry_string = dbcsr_type_antisymmetric
          ELSE
             symmetry_string = dbcsr_type_symmetric
          END IF
       ELSE
          symmetry_string = dbcsr_type_no_symmetry
       END IF
       ALLOCATE(matrix_s(i)%matrix, STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       CALL cp_dbcsr_init(matrix_s(i)%matrix,error=error)
       CALL cp_dbcsr_create(matrix=matrix_s(i)%matrix, &
            name=TRIM(mname), &
            dist=dbcsr_dist, matrix_type=symmetry_string,&
            row_blk_size=row_blk_sizes, col_blk_size=col_blk_sizes, &
            nblks=0, nze=0, error=error)
       CALL cp_dbcsr_alloc_block_from_nbl(matrix_s(i)%matrix,sab_nl,error)
     END DO

     CALL array_release (row_blk_sizes)
     CALL array_release (col_blk_sizes)

     ! *** Allocate work storage ***
     max_maxco = 0
     max_maxlgto = 0
     max_maxsgf = 0
     CALL get_gto_basis_set_info(gto_basis=basis_set_list_a,maxco=maxco,maxl=maxlgto,maxsgf=maxsgf)
     max_maxco = MAX(maxco,max_maxco)
     max_maxlgto = MAX(maxlgto,max_maxlgto)
     max_maxsgf = MAX(maxsgf,max_maxsgf)

     CALL get_gto_basis_set_info(gto_basis=basis_set_list_b,maxco=maxco,maxl=maxlgto,maxsgf=maxsgf)
     max_maxco = MAX(maxco,max_maxco)
     max_maxlgto = MAX(maxlgto,max_maxlgto)
     max_maxsgf = MAX(maxsgf,max_maxsgf)
     ldsab = MAX(max_maxco,max_maxsgf)

     CALL get_gto_basis_set_info(gto_basis=basis_set_list_a,maxco=maxdcoa,maxder=MAX(1,nder))
     CALL get_gto_basis_set_info(gto_basis=basis_set_list_b,maxco=maxdcob,maxder=MAX(1,nder))

     return_s_derivatives = (nder > 0)
     maxblock = MAX(maxder,4)
     maxdco = MAX(maxdcoa, maxdcob)
     ldai = MAX(ncoset(max_maxlgto+nder+1),ldsab)
     CALL init_orbital_pointers(ldai)

     ALLOCATE (ai_work(ldai,ldai,maxder),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     ALLOCATE (sab(ldsab,ldsab*maxblock),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     sab=0._dp
     ALLOCATE (work(max_maxco,max_maxsgf),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     work=0._dp
     ALLOCATE (sint(maxder),STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     DO i=1,maxder
       NULLIFY (sint(i)%block)
     END DO

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
        sphi_a       =>  basis_set_a%sphi
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
        sphi_b       =>  basis_set_b%sphi
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
            IF (iatom==jatom) THEN
               f0 = 1.0_dp
            ELSE
               f0 = 2.0_dp
            END IF
          ELSE
            irow = iatom
            icol = jatom
            f0 = 1.0_dp
          END IF
          DO i=1,SIZE(matrix_s)
            NULLIFY (sint(i)%block)
            CALL cp_dbcsr_get_block_p(matrix=matrix_s(i)%matrix,&
                 row=irow,col=icol,BLOCK=sint(i)%block,found=found)
                 CPPostcondition(found,cp_failure_level,routineP,error,failure)
          END DO
        END IF

        rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
        dab = SQRT(rab2)

        DO iset=1,nseta

          ncoa = npgfa(iset)*ncoset(la_max(iset))
          sgfa = first_sgfa(1,iset)

          DO jset=1,nsetb

            IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

            ncob = npgfb(jset)*ncoset(lb_max(jset))
            sgfb = first_sgfb(1,jset)

            ! *** Calculate the primitive overlap integrals ***
            CALL overlap(la_max(iset),la_min(iset),npgfa(iset),&
                         rpgfa(:,iset),zeta(:,iset),&
                         lb_max(jset),lb_min(jset),npgfb(jset),&
                         rpgfb(:,jset),zetb(:,jset),&
                         rab,dab,sab,nder,return_s_derivatives,&
                         ai_work,ldai)

            ! *** Contraction step (overlap matrix and its derivatives) ***
            DO i=1,SIZE(matrix_s)
               iblock = (i - 1)*SIZE(sab,1) + 1
               CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,1.0_dp,sab(1,iblock),SIZE(sab,1),&
                    sphi_b(1,sgfb),SIZE(sphi_b,1),0.0_dp,work(1,1),SIZE(work,1))
               f = 1.0_dp
               IF ( do_symmetric .AND. (iatom > jatom) ) THEN
                  ! *** The odd derivatives are anti-symmetric ***
                  IF (ndod(i) == 1) f = -1.0_dp
                  CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,f,work(1,1),SIZE(work,1),&
                       sphi_a(1,sgfa),SIZE(sphi_a,1),1.0_dp,sint(i)%block(sgfb,sgfa),SIZE(sint(i)%block,1))
               ELSE
                  CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,f,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                       work(1,1),SIZE(work,1),1.0_dp,sint(i)%block(sgfa,sgfb),SIZE(sint(i)%block,1))
                END IF
             END DO

          END DO
        END DO

     END DO
     CALL neighbor_list_iterator_release(nl_iterator)

     DO i = 1,SIZE(matrix_s)
        CALL cp_dbcsr_finalize(matrix_s(i)%matrix, error=error)
        IF(PRESENT(eps_filter)) THEN
           CALL cp_dbcsr_filter(matrix_s(i)%matrix,eps_filter,error=error)
        END IF
     ENDDO

     ! *** Release work storage ***
     DEALLOCATE (ai_work,sab,work,STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
     DO i=1,maxder
        NULLIFY (sint(i)%block)
     END DO
     DEALLOCATE (sint,STAT=istat)
     CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

     CALL timestop(handle)

  END SUBROUTINE build_overlap_matrix_simple

! *****************************************************************************

END MODULE qs_overlap

