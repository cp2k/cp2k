!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Distribution of the overlap integral matrix. The calculation of the
!>      overlap integral matrix is performed in MODULE core_hamiltonian.
!> \par History
!>      JGH: removed printing routines
!> \author Matthias Krack (03.09.2001,25.06.2003)
! *****************************************************************************
MODULE qs_overlap

  USE ai_overlap_new,                  ONLY: overlap
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE block_p_types,                   ONLY: block_p_type
  USE cp_dbcsr_operations,             ONLY: dbcsr_add_block_node,&
                                             dbcsr_allocate_matrix,&
                                             dbcsr_allocate_matrix_set,&
                                             dbcsr_deallocate_matrix_set,&
                                             dbcsr_obj_p_type,&
                                             dbcsr_replicate_matrix_structure,&
                                             sm_from_dbcsr
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dbcsr_util,                      ONLY: dbcsr_finalize
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE input_constants,                 ONLY: use_aux_fit_basis_set,&
                                             use_orb_basis_set
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE orbital_pointers,                ONLY: indco,&
                                             init_orbital_pointers,&
                                             ncoset
  USE orbital_symbols,                 ONLY: cgf_symbol
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_set_type, &
       neighbor_list_type, neighbor_node_type, next
  USE sparse_matrix_output,            ONLY: write_matrix_distribution,&
                                             write_sparse_matrix
  USE sparse_matrix_types,             ONLY: allocate_matrix_set,&
                                             real_matrix_p_type
  USE string_utilities,                ONLY: compress,&
                                             uppercase
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_overlap'

! *** Public subroutines ***

  PUBLIC :: build_overlap_matrix, &
            atom_overlap_matrix

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
!> \param   qs_env the QS env
!> \param   para_env the parallel env
!> \param   nderivative Derivative with respect to basis origin
!> \param   matrix_s The overlap matrix to be calculated
!> \param   matrix_name The name of the overlap matrix (i.e. for output)
!> \param   basis_set_id_a basis set to be used for bra in <a|b>
!> \param   basis_set_id_b basis set to be used for ket in <a|b>
!> \param   neighbor_list_sab pair list (must be consistent with basis sets!)
!> \param   error for error handling
!> \date    11.03.2002
!> \version 1.0
! *****************************************************************************

  SUBROUTINE build_overlap_matrix(qs_env,para_env,nderivative,matrix_s_sm,matrix_name,&
                                  basis_set_id_a, basis_set_id_b, neighbor_list_sab, &
                                  error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN), OPTIONAL            :: nderivative
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s_sm
    CHARACTER(LEN=*), INTENT(IN)             :: matrix_name
    INTEGER, INTENT(IN)                      :: basis_set_id_a, basis_set_id_b
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: neighbor_list_sab
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_overlap_matrix', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=12)                        :: cgfsym
    CHARACTER(LEN=32)                        :: symmetry_string
    CHARACTER(LEN=80)                        :: name
    INTEGER :: handle, i, iab, iatom, iblock, icol, ikind, ilist, inode, &
      irow, iset, istat, iw, jatom, jkind, jset, last_jatom, ldai, max_maxco, &
      max_maxlgto, max_maxsgf, maxco, maxder, maxlgto, maxsgf, natom, ncoa, &
      ncob, nder, neighbor_list_id, nkind, nlist, nnode, nseta, nsetb, &
      nsgf_a, nsgf_b, sgfa, sgfb
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf_a, first_sgf_b, &
                                                last_sgf_a, last_sgf_b
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: do_symmetric, new_atom_b, &
                                                return_s_derivatives
    REAL(KIND=dp)                            :: dab, f, rab2
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: sab, work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: ai_work
    REAL(KIND=dp), DIMENSION(3)              :: rab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                zeta, zetb
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(block_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: sint
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dbcsr_obj_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b
    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    TYPE(neighbor_list_type), POINTER        :: neighbor_list
    TYPE(neighbor_node_type), POINTER        :: neighbor_node
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    CALL timeset(routineN,handle)
    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    NULLIFY ( atomic_kind_set,particle_set,distribution_2d, matrix_s)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    neighbor_list_id=neighbor_list_id, &
                    distribution_2d=distribution_2d,error=error)

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    IF (PRESENT(nderivative)) THEN
      nder = nderivative
      return_s_derivatives = .TRUE.
    ELSE
      nder = 0
      return_s_derivatives = .FALSE.
    END IF

    maxder = ncoset(nder)
    CALL dbcsr_allocate_matrix_set(matrix_s,maxder,error=error)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,nsgf=nsgf_a,basis_set_id=basis_set_id_a)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,nsgf=nsgf_b,basis_set_id=basis_set_id_b)

    ALLOCATE (first_sgf_a(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "first_sgf_a",natom*int_size)
    ALLOCATE (last_sgf_a(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "last_sgf_a",natom*int_size)
    ALLOCATE (first_sgf_b(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "first_sgf_b",natom*int_size)
    ALLOCATE (last_sgf_b(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "last_sgf_b",natom*int_size)


    CALL get_particle_set(particle_set=particle_set,&
                          first_sgf=first_sgf_a,&
                          last_sgf=last_sgf_a,&
                          basis_set_id=basis_set_id_a,error=error)
    CALL get_particle_set(particle_set=particle_set,&
                          first_sgf=first_sgf_b,&
                          last_sgf=last_sgf_b,&
                          basis_set_id=basis_set_id_b,error=error)

   do_symmetric = .TRUE.
   symmetry_string = "symmetric"
   IF( basis_set_id_a /= basis_set_id_b ) THEN
     do_symmetric = .FALSE.
     symmetry_string = "none"
   END IF
   ! XXXXXXX, sparsity_id, does this result in the same overlap matrix
   ! XXXXXXX, as qs_build_core_hamiltonian, if so it is a pitty to have
   ! XXXXXXX, the code duplication,if not, this sparsity_id is wrong.
   CALL dbcsr_allocate_matrix(matrix=matrix_s(1)%matrix,&
                         nblock_row=natom,&
                         nblock_col=natom,&
                         nrow=nsgf_a,&
                         ncol=nsgf_b,&
                         first_row=first_sgf_a,&
                         last_row=last_sgf_a,&
                         first_col=first_sgf_b,&
                         last_col=last_sgf_b,&
                         matrix_name=TRIM(matrix_name),&
                         matrix_symmetry=symmetry_string,&
                         sparsity_id=neighbor_list_id, &
                         distribution_2d=distribution_2d,error=error)

    DEALLOCATE (first_sgf_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "first_sgf_a")
    DEALLOCATE (last_sgf_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "last_sgf_a")

    DEALLOCATE (first_sgf_b,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "first_sgf_b")
    DEALLOCATE (last_sgf_b,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "last_sgf_b")

    symmetry_string = "antisymmetric"
    IF( basis_set_id_a /= basis_set_id_b) THEN
      symmetry_string = "none"
    END IF
    DO i=2,maxder
      cgfsym = cgf_symbol(1,indco(1:3,i))
      name = TRIM(cgfsym(4:))//" DERIVATIVE OF THE "//TRIM(matrix_name)//&
             " W.R.T. THE NUCLEAR COORDINATES"
      CALL compress(name)
      CALL uppercase(name)
      CALL dbcsr_replicate_matrix_structure(source=matrix_s(1)%matrix,&
                                            TARGET=matrix_s(i)%matrix,&
                                            target_name=TRIM(name),&
                                            target_symmetry=TRIM(symmetry_string),error=error)
    END DO

!   *** Allocate work storage ***
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

    ldai = ncoset(max_maxlgto+nder)
    CALL init_orbital_pointers(ldai)

    ALLOCATE (ai_work(ldai,ldai,ncoset(nder)),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineN,moduleN,__LINE__,&
                       "ai_work",ldai*ldai*ncoset(nder)*dp_size)
    END IF

    ALLOCATE (sab(max_maxco,max_maxco*maxder),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "sab",max_maxco*max_maxco*maxder*dp_size)
    sab(:,:) = 0.0_dp

    ALLOCATE (work(max_maxco,max_maxsgf),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "work",max_maxco*max_maxsgf*dp_size)
    work(:,:) = 0.0_dp

    ALLOCATE (sint(maxder),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "sint",maxder*int_size)
    DO i=1,maxder
      NULLIFY (sint(i)%block)
    END DO

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)
      SELECT CASE (basis_set_id_a)
      CASE (use_orb_basis_set)
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             orb_basis_set=basis_set_a)
      CASE (use_aux_fit_basis_set)
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             aux_fit_basis_set=basis_set_a)
      END SELECT

      IF (.NOT.ASSOCIATED(basis_set_a)) CYCLE

      CALL get_gto_basis_set(gto_basis_set=basis_set_a,&
                             first_sgf=first_sgfa,&
                             lmax=la_max,&
                             lmin=la_min,&
                             npgf=npgfa,&
                             nset=nseta,&
                             nsgf_set=nsgfa,&
                             pgf_radius=rpgfa,&
                             set_radius=set_radius_a,&
                             sphi=sphi_a,&
                             zet=zeta)

      DO jkind=1,nkind

        atomic_kind => atomic_kind_set(jkind)
        SELECT CASE (basis_set_id_b)
        CASE (use_orb_basis_set)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               orb_basis_set=basis_set_b)
        CASE (use_aux_fit_basis_set)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               aux_fit_basis_set=basis_set_b)
        END SELECT


        IF (.NOT.ASSOCIATED(basis_set_b)) CYCLE

        CALL get_gto_basis_set(gto_basis_set=basis_set_b,&
                               first_sgf=first_sgfb,&
                               lmax=lb_max,&
                               lmin=lb_min,&
                               npgf=npgfb,&
                               nset=nsetb,&
                               nsgf_set=nsgfb,&
                               pgf_radius=rpgfb,&
                               set_radius=set_radius_b,&
                               sphi=sphi_b,&
                               zet=zetb)

        iab = ikind + nkind*(jkind - 1)

        IF (.NOT.ASSOCIATED(neighbor_list_sab(iab)%neighbor_list_set)) CYCLE

        neighbor_list_set => neighbor_list_sab(iab)%neighbor_list_set

        CALL get_neighbor_list_set(neighbor_list_set=neighbor_list_set,&
                                   nlist=nlist)

        neighbor_list => first_list(neighbor_list_set)

        DO ilist=1,nlist

          CALL get_neighbor_list(neighbor_list=neighbor_list,&
                                 atom=iatom,&
                                 nnode=nnode)

          last_jatom = 0

          neighbor_node => first_node(neighbor_list)

          DO inode=1,nnode

            CALL get_neighbor_node(neighbor_node=neighbor_node,&
                                   neighbor=jatom,&
                                   r=rab(:))

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
              ELSE
                irow = iatom
                icol = jatom
              END IF
              DO i=1,maxder
                NULLIFY (sint(i)%block)
                CALL dbcsr_add_block_node(matrix=matrix_s(i)%matrix,&
                                    block_row=irow,&
                                    block_col=icol,&
                                    BLOCK=sint(i)%block,error=error)
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

!               *** Calculate the primitive overlap integrals ***

                CALL overlap(la_max(iset),la_min(iset),npgfa(iset),&
                             rpgfa(:,iset),zeta(:,iset),&
                             lb_max(jset),lb_min(jset),npgfb(jset),&
                             rpgfb(:,jset),zetb(:,jset),&
                             rab,dab,sab,nder,return_s_derivatives,&
                             ai_work,ldai)

!               *** Contraction step (overlap matrix and its derivatives) ***

                DO i=1,maxder

                  iblock = (i - 1)*SIZE(sab,1) + 1

                  CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                             1.0_dp,sab(1,iblock),SIZE(sab,1),&
                             sphi_b(1,sgfb),SIZE(sphi_b,1),&
                             0.0_dp,work(1,1),SIZE(work,1))

                  IF ( do_symmetric ) THEN
                    IF (iatom <= jatom) THEN

                      CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                                 1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 work(1,1),SIZE(work,1),&
                                 1.0_dp,sint(i)%block(sgfa,sgfb),&
                                 SIZE(sint(i)%block,1))

                    ELSE

  !                   *** The first derivatives are anti-symmetric ***

                      IF (i > 1) THEN
                        f = -1.0_dp
                      ELSE
                        f = 1.0_dp
                      END IF

                      CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                                 f,work(1,1),SIZE(work,1),&
                                 sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 1.0_dp,sint(i)%block(sgfb,sgfa),&
                                 SIZE(sint(i)%block,1))
                    END IF
                  ELSE
                    CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                               1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                               work(1,1),SIZE(work,1),&
                               1.0_dp,sint(i)%block(sgfa,sgfb),&
                               SIZE(sint(i)%block,1))

!                   *** The first derivatives are anti-symmetric ***
!!!!!!******** This part needs to be fixed in case of a non-symmetric overlap matrix****!!!!!
                    IF (i > 1) THEN
                      f = -1.0_dp
                    ELSE
                      f = 1.0_dp
                    END IF

!                    CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
!                               f,work(1,1),SIZE(work,1),&
!                               sphi_a(1,sgfa),SIZE(sphi_a,1),&
!                               1.0_dp,sint(i)%block(sgfb,sgfa),&
!                               SIZE(sint(i)%block,1))

!!!!!!**********************************************************************************!!!!!!
                  END IF
                END DO

              END DO
            END DO

            neighbor_node => next(neighbor_node)

          END DO

          neighbor_list => next(neighbor_list)

        END DO

      END DO
    END DO

    DO i = 1,SIZE(matrix_s)
       CALL dbcsr_finalize(matrix_s(i)%matrix, error=error)
    ENDDO

!   *** Release work storage ***

    DEALLOCATE (ai_work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "ai_work")

    DEALLOCATE (sab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "sab")

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "work")

    DO i=1,maxder
      NULLIFY (sint(i)%block)
    END DO
    DEALLOCATE (sint,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "sint")

    CALL allocate_matrix_set(matrix_s_sm,SIZE(matrix_s),error=error)
    DO i = 1,SIZE(matrix_s)
       CALL sm_from_dbcsr ( matrix_s_sm(i)%matrix, matrix_s(i)%matrix, distribution_2d, error=error )
    ENDDO
    CALL dbcsr_deallocate_matrix_set ( matrix_s, error=error )


!   *** Print the overlap matrix distribution ***

    CALL write_matrix_distribution(matrix_s_sm(1)%matrix,qs_env%input,para_env,error)

!   *** Print the overlap matrix, if requested ***

    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         qs_env%input,"DFT%PRINT%AO_MATRICES/OVERLAP",error=error),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/OVERLAP",&
            extension=".Log",error=error)
       CALL write_sparse_matrix(matrix_s_sm(1)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
       IF (BTEST(cp_print_key_should_output(logger%iter_info,&
            qs_env%input,"DFT%PRINT%AO_MATRICES/DERIVATIVES",error=error),cp_p_file)) THEN
          DO i=2,maxder
             CALL write_sparse_matrix(matrix_s_sm(i)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
          END DO
       END IF
       CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/OVERLAP", error=error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE build_overlap_matrix

! *****************************************************************************
  SUBROUTINE atom_overlap_matrix(iatom,nseta,set_radius_a,npgfa,nsgfa,la_max,&
                                 la_min,first_sgfa,rpgfa,zeta,sphi_a,&
                                 jatom,nsetb,set_radius_b,npgfb,nsgfb,lb_max,&
                                 lb_min,first_sgfb,rpgfb,zetb,sphi_b,&
                                 rab,nder,ldai,sint,sab,work,ai_work,error)

    INTEGER, INTENT(IN)                      :: iatom, nseta
    REAL(dp), DIMENSION(:), INTENT(IN)       :: set_radius_a
    INTEGER, DIMENSION(:), INTENT(IN)        :: npgfa, nsgfa, la_max, la_min
    INTEGER, DIMENSION(:, :), INTENT(IN)     :: first_sgfa
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: rpgfa, zeta, sphi_a
    INTEGER, INTENT(IN)                      :: jatom, nsetb
    REAL(dp), DIMENSION(:), INTENT(IN)       :: set_radius_b
    INTEGER, DIMENSION(:), INTENT(IN)        :: npgfb, nsgfb, lb_max, lb_min
    INTEGER, DIMENSION(:, :), INTENT(IN)     :: first_sgfb
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: rpgfb, zetb, sphi_b
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rab
    INTEGER, INTENT(IN)                      :: nder, ldai
    TYPE(block_p_type), DIMENSION(:)         :: sint
    REAL(KIND=dp), DIMENSION(:, :)           :: sab, work
    REAL(KIND=dp), DIMENSION(:, :, :)        :: ai_work
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, iblock, iset, jset, &
                                                maxder, ncoa, ncob, sgfa, sgfb
    LOGICAL                                  :: return_s_derivatives
    REAL(dp)                                 :: dab, rab2

!-------------------------------------------------------------------------------

    rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
    dab = SQRT(rab2)
    return_s_derivatives = (nder /= 0)
    maxder = ncoset(nder)

    DO iset=1,nseta

       ncoa = npgfa(iset)*ncoset(la_max(iset))
       sgfa = first_sgfa(1,iset)

       DO jset=1,nsetb

          IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

          ncob = npgfb(jset)*ncoset(lb_max(jset))
          sgfb = first_sgfb(1,jset)

!         *** Calculate the primitive overlap integrals ***

          CALL overlap(la_max(iset),la_min(iset),npgfa(iset),&
                       rpgfa(:,iset),zeta(:,iset),&
                       lb_max(jset),lb_min(jset),npgfb(jset),&
                       rpgfb(:,jset),zetb(:,jset),&
                       rab,dab,sab,nder,return_s_derivatives,&
                       ai_work,ldai)

!         *** Contraction step (overlap matrix and its derivatives) ***

          DO i=1,maxder

            iblock = (i - 1)*SIZE(sab,1) + 1

            CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                       1.0_dp,sab(1,iblock),SIZE(sab,1),&
                       sphi_b(1,sgfb),SIZE(sphi_b,1),&
                       0.0_dp,work(1,1),SIZE(work,1))

            CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                       1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                       work(1,1),SIZE(work,1),&
                       1.0_dp,sint(i)%block(sgfa,sgfb),&
                       SIZE(sint(i)%block,1))

          END DO

       END DO

    END DO

  END SUBROUTINE atom_overlap_matrix

END MODULE qs_overlap

