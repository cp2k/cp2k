!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_overlap [1.0] *
!!
!!   NAME
!!     qs_overlap
!!
!!   FUNCTION
!!     Distribution of the overlap integral matrix. The calculation of the
!!     overlap integral matrix is performed in MODULE core_hamiltonian.
!!
!!   AUTHOR
!!     Matthias Krack (03.09.2001,25.06.2003)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE qs_overlap

  USE ai_overlap_new,                  ONLY: overlap
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_type
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE machine,                         ONLY: m_flush
  USE mathlib,                         ONLY: symmetrize_matrix
  USE message_passing,                 ONLY: mp_sum,&
                                             mp_sync
  USE orbital_pointers,                ONLY: indco,&
                                             init_orbital_pointers,&
                                             ncoset,&
                                             nso
  USE orbital_symbols,                 ONLY: cgf_symbol
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             set_qs_env,&
                                             qs_environment_type
  USE qs_neighbor_list_types,          ONLY: first_list,&
                                             first_node,&
                                             get_neighbor_list,&
                                             get_neighbor_list_set,&
                                             get_neighbor_node,&
                                             neighbor_list_set_p_type,&
                                             neighbor_list_set_type,&
                                             neighbor_list_type,&
                                             neighbor_node_type,&
                                             next
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             allocate_matrix,&
                                             allocate_matrix_set,&
                                             copy_local_sm_to_replicated_fm,&
                                             get_matrix_info,&
                                             real_matrix_p_type,&
                                             real_matrix_type,&
                                             replicate_matrix_structure
  USE string_utilities,                ONLY: compress,&
                                             uppercase
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_overlap"

! *** Public subroutines ***

  PUBLIC :: build_overlap_matrix,&
            write_fm_with_basis_info,&
            write_matrix_distribution,&
            write_sparse_matrix

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE build_overlap_matrix(qs_env,globenv,nderivative)

!   Purpose: Calculation of the overlap matrix over Cartesian Gaussian
!            functions.

!   History: - Creation (11.03.2002, MK)

!   ***************************************************************************

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    INTEGER,INTENT(IN), OPTIONAL             :: nderivative

!   *** Local types ***

    TYPE block_p_type
      REAL(KIND = dp), DIMENSION(:,:), POINTER :: block
    END TYPE block_p_type

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "build_overlap_matrix"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER       :: atomic_kind
    TYPE(gto_basis_set_type), POINTER     :: orb_basis_set
    TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set
    TYPE(neighbor_list_type), POINTER     :: neighbor_list
    TYPE(neighbor_node_type), POINTER     :: neighbor_node

    CHARACTER(LEN=12) :: cgfsym
    CHARACTER(LEN=80) :: name
    REAL(KIND = dp)   :: dab,f,rab2
    INTEGER           :: handle,i,iab,iatom,iblock,icol,ikind,ilist,&
                         inode,irow,iset,istat,jatom,jkind,&
                         jset,last_jatom,ldai,maxco,maxder,maxlgto,&
                         maxsgf,natom,ncoa,ncob,nder,nkind,nlist,&
                         nnode,nseta,nsetb,nsgf,sgfa,sgfb
    LOGICAL           :: new_atom_b,return_s_derivatives

    TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
    TYPE(distribution_2d_type), POINTER                   :: distribution_2d
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb
    TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER       :: matrix_s

    REAL(KIND = dp), DIMENSION(3)                :: rab
    REAL(KIND = dp), DIMENSION(:,:), ALLOCATABLE :: sab,work

    INTEGER, DIMENSION(:), ALLOCATABLE :: first_sgf,last_sgf

    REAL(KIND = dp), DIMENSION(:,:,:), ALLOCATABLE :: ai_work

    REAL(KIND = dp), DIMENSION(:), POINTER   :: set_radius_a,set_radius_b
    INTEGER, DIMENSION(:), POINTER           :: la_max,la_min,lb_max,lb_min,&
                                                npgfa,npgfb,nsgfa,nsgfb
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: rpgfa,rpgfb,sphi_a,sphi_b,&
                                                zeta,zetb
    INTEGER, DIMENSION(:,:), POINTER         :: first_sgfa,first_sgfb

    TYPE(block_p_type), DIMENSION(:), ALLOCATABLE :: sint
    INTEGER                                  :: neighbor_list_id

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)

    NULLIFY ( atomic_kind_set,particle_set,matrix_s,sab_orb,distribution_2d)
    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    matrix_s=matrix_s,&
                    sab_orb=sab_orb,&
                    neighbor_list_id=neighbor_list_id, &
                    distribution_2d=distribution_2d)

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
    CALL allocate_matrix_set(matrix_s,maxder)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,nsgf=nsgf)

    ALLOCATE (first_sgf(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "first_sgf",natom*int_size)
    ALLOCATE (last_sgf(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "last_sgf",natom*int_size)

    CALL get_particle_set(particle_set=particle_set,&
                          first_sgf=first_sgf,&
                          last_sgf=last_sgf)
   ! XXXXXXX, sparsity_id, does this result in the same overlap matrix as qs_build_core_hamiltonian
   ! XXXXXXX, if so it is a pitty to have the code duplication,if not, this sparsity_id is wrong.
    CALL allocate_matrix(matrix=matrix_s(1)%matrix,&
                         nblock_row=natom,&
                         nblock_col=natom,&
                         nrow=nsgf,&
                         ncol=nsgf,&
                         first_row=first_sgf,&
                         last_row=last_sgf,&
                         first_col=first_sgf,&
                         last_col=last_sgf,&
                         matrix_name="OVERLAP MATRIX",&
                         matrix_symmetry="symmetric",&
                         sparsity_id=neighbor_list_id, &  
                         distribution_2d=distribution_2d)

    DEALLOCATE (first_sgf,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "first_sgf")
    DEALLOCATE (last_sgf,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "last_sgf")

    DO i=2,maxder
      cgfsym = cgf_symbol(1,indco(1:3,i))
      name = TRIM(cgfsym(4:))//" DERIVATIVE OF THE OVERLAP MATRIX "//&
             "W.R.T. THE NUCLEAR COORDINATES"
      CALL compress(name)
      CALL uppercase(name)
      CALL replicate_matrix_structure(source=matrix_s(1)%matrix,&
                                      TARGET=matrix_s(i)%matrix,&
                                      target_name=TRIM(name),&
                                      target_symmetry="antisymmetric")
    END DO

!   *** Allocate work storage ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             maxlgto=maxlgto,&
                             maxsgf=maxsgf)

    ldai = ncoset(maxlgto+nder)
    CALL init_orbital_pointers(ldai)

    ALLOCATE (ai_work(ldai,ldai,ncoset(nder)),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,&
                       "ai_work",ldai*ldai*ncoset(nder)*dp_size)
    END IF

    ALLOCATE (sab(maxco,maxco*maxder),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "sab",maxco*maxco*maxder*dp_size)
    sab(:,:) = 0.0_dp

    ALLOCATE (work(maxco,maxsgf),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "work",maxco*maxsgf*dp_size)
    work(:,:) = 0.0_dp

    ALLOCATE (sint(maxder),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "sint",maxder*int_size)
    DO i=1,maxder
      NULLIFY (sint(i)%block)
    END DO

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           orb_basis_set=orb_basis_set)

      IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

      CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
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

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             orb_basis_set=orb_basis_set)

        IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
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

        IF (.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE

        neighbor_list_set => sab_orb(iab)%neighbor_list_set

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
              IF (iatom <= jatom) THEN
                irow = iatom
                icol = jatom
              ELSE
                irow = jatom
                icol = iatom
              END IF
              DO i=1,maxder
                NULLIFY (sint(i)%block)
                CALL add_block_node(matrix=matrix_s(i)%matrix,&
                                    block_row=irow,&
                                    block_col=icol,&
                                    BLOCK=sint(i)%block)
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

                END DO

              END DO
            END DO

            neighbor_node => next(neighbor_node)

          END DO

          neighbor_list => next(neighbor_list)

        END DO

      END DO
    END DO

!   *** Release work storage ***

    DEALLOCATE (ai_work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "ai_work")

    DEALLOCATE (sab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "sab")

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "work")

    DO i=1,maxder
      NULLIFY (sint(i)%block)
    END DO
    DEALLOCATE (sint,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "sint")

    CALL set_qs_env(qs_env=qs_env,&
                    matrix_s=matrix_s)
    CALL timestop(0.0_dp,handle)

!   *** Print the overlap matrix distribution ***

    CALL write_matrix_distribution(matrix_s(1)%matrix,globenv)

!   *** Print the overlap matrix, if requested ***

    IF (globenv%print%overlap_matrix) THEN
      CALL write_sparse_matrix(matrix_s(1)%matrix,4,6,qs_env,globenv)
      IF (globenv%print%derivatives) THEN
        DO i=2,maxder
          CALL write_sparse_matrix(matrix_s(i)%matrix,4,6,qs_env,globenv)
        END DO
      END IF
    END IF

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE build_overlap_matrix

! *****************************************************************************

  SUBROUTINE write_fm_with_basis_info(blacs_matrix,before,after,qs_env,globenv,&
                                      first_row,last_row,first_col,last_col)

!   Purpose: Print a spherical matrix of blacs type.

!   History: - Creation (12.06.2001,MK)
!            - Allow for printing of a sub-matrix (01.07.2003,MK)

!   ***************************************************************************

    TYPE(cp_fm_type), POINTER                :: blacs_matrix
    INTEGER, INTENT(IN)                      :: before, after
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_col, last_col

    CHARACTER(LEN=*), PARAMETER :: routine_name = "write_fm_with_basis_info"

    CHARACTER(LEN=60)                        :: matrix_name
    INTEGER                                  :: col1, col2, group, istat, &
                                                ncol_global, nrow_global, &
                                                nsgf, output_unit, row1, row2
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: matrix
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set

!   ---------------------------------------------------------------------------

    group = globenv%group
    output_unit = globenv%scr
    IF (.NOT.ASSOCIATED(blacs_matrix)) RETURN
    CALL cp_fm_get_info(blacs_matrix,name=matrix_name,nrow_global=nrow_global,ncol_global=ncol_global)

    ALLOCATE(matrix(nrow_global,ncol_global))
    CALL cp_fm_get_submatrix(blacs_matrix,matrix)


!   *** Get the matrix dimension and check the optional arguments ***

    CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,nsgf=nsgf)

    IF (PRESENT(first_row)) THEN
      row1 = MAX(1,first_row)
    ELSE
      row1 = 1
    END IF

    IF (PRESENT(last_row)) THEN
      row2 = MIN(nsgf,last_row)
    ELSE
      row2 = nsgf
    END IF

    IF (PRESENT(first_col)) THEN
      col1 = MAX(1,first_col)
    ELSE
      col1 = 1
    END IF

    IF (PRESENT(last_col)) THEN
      col2 = MIN(nsgf,last_col)
    ELSE
      col2 = nsgf
    END IF

    CALL write_matrix(matrix,matrix_name,before,after,qs_env,globenv,&
                      row1,row2,col1,col2)

!   *** Release work storage ***

    IF (ASSOCIATED(matrix)) THEN
      DEALLOCATE (matrix,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "matrix")
    END IF

  END SUBROUTINE write_fm_with_basis_info

! *****************************************************************************

  SUBROUTINE write_matrix(matrix,matrix_name,before,after,qs_env,globenv,&
                          first_row,last_row,first_col,last_col)

!   Purpose: Write a matrix or a sub-matrix to the output unit.

!   History: - Creation (01.07.2003,MK)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: matrix
    CHARACTER(LEN=*), INTENT(IN)             :: matrix_name
    INTEGER, INTENT(IN)                      :: before, after
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    INTEGER, INTENT(IN)                      :: first_row, last_row, &
                                                first_col, last_col

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "write_matrix"

    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=25)                        :: fmtstr1
    CHARACTER(LEN=35)                        :: fmtstr2
    CHARACTER(LEN=6), DIMENSION(:), POINTER  :: sgf_symbol
    INTEGER :: from, group, iatom, icol, irow, iset, isgf, ishell, iso, jcol, &
      l, left, natom, ncol, ndigits, nset, nsgf, output_unit, right, to, width
    INTEGER, DIMENSION(:), POINTER           :: nshell
    INTEGER, DIMENSION(:, :), POINTER        :: lshell
    LOGICAL                                  :: ionode
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

!   ---------------------------------------------------------------------------

    group = globenv%group
    ionode = globenv%ionode
    output_unit = globenv%scr

    CALL m_flush(output_unit)

    IF (ionode) THEN

      CALL get_qs_env(qs_env=qs_env,&
                      atomic_kind_set=atomic_kind_set,&
                      particle_set=particle_set)

      natom = SIZE(particle_set)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,nsgf=nsgf)

!     *** Definition of the variable formats ***

      fmtstr1 = "(/,T2,23X,  (  X,I5,  X))"
      fmtstr2 = "(T2,2I5,2X,A2,1X,A8,   (1X,F  .  ))"

!     *** Write headline ***

      WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") TRIM(matrix_name)

!     *** Write the variable format strings ***

      ndigits = after + globenv%print%precision

      width = before + ndigits + 3
      ncol = INT(56/width)

      right = MAX((ndigits-2),1)
      left =  width - right - 5

      WRITE (UNIT=fmtstr1(11:12),FMT="(I2)") ncol
      WRITE (UNIT=fmtstr1(14:15),FMT="(I2)") left
      WRITE (UNIT=fmtstr1(21:22),FMT="(I2)") right

      WRITE (UNIT=fmtstr2(22:23),FMT="(I2)") ncol
      WRITE (UNIT=fmtstr2(29:30),FMT="(I2)") width - 1
      WRITE (UNIT=fmtstr2(32:33),FMT="(I2)") ndigits

!     *** Write the matrix in the selected format ***

      DO icol=first_col,last_col,ncol
        from = icol
        to = MIN((from+ncol-1),last_col)
        WRITE (UNIT=output_unit,FMT=fmtstr1) (jcol,jcol=from,to)
        irow = 1
        DO iatom=1,natom
          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                               element_symbol=element_symbol,&
                               orb_basis_set=orb_basis_set)
          CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                 nset=nset,&
                                 nshell=nshell,&
                                 l=lshell,&
                                 sgf_symbol=sgf_symbol)
          isgf = 1
          DO iset=1,nset
            DO ishell=1,nshell(iset)
              l = lshell(ishell,iset)
              DO iso=1,nso(l)
                IF ((irow >= first_row).AND.(irow <= last_row)) THEN
                  WRITE (UNIT=output_unit,FMT=fmtstr2)&
                    irow,iatom,element_symbol,sgf_symbol(isgf),&
                    (matrix(irow,jcol),jcol=from,to)
                END IF
                isgf = isgf + 1
                irow = irow + 1
              END DO
            END DO
          END DO
          IF ((irow >= first_row).AND.(irow <= last_row)) THEN
            WRITE (UNIT=output_unit,FMT="(A)")
          END IF
        END DO
      END DO

      WRITE (UNIT=output_unit,FMT="(/)")

    END IF

    CALL mp_sync(group)
    CALL m_flush(output_unit)

  END SUBROUTINE write_matrix

! *****************************************************************************

  SUBROUTINE write_matrix_distribution(matrix,globenv)

!   Purpose: Print the distribution of a sparse matrix.

!   History: - Creation (25.06.2003,MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER          :: matrix
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "write_matrix_distribution"

    CHARACTER(LEN=80)                        :: matrix_name, matrix_symmetry
    INTEGER                                  :: group, ipe, istat, mype, &
                                                natom, nblock_sum, &
                                                nblock_tot, nelement_sum, &
                                                npe, nrow, output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: nblock, nelement
    LOGICAL                                  :: ionode
    REAL(KIND=dp)                            :: occupation

!   ---------------------------------------------------------------------------
!   *** Quick return, if no printing is requested ***

    IF (.NOT.globenv%print%distribution) RETURN

    group = globenv%group
    ionode = globenv%ionode
    mype = globenv%mepos + 1
    npe = globenv%num_pe
    output_unit = globenv%scr

!   *** Allocate work storage ***

    ALLOCATE (nblock(npe),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nblock",npe*int_size)
    nblock(:) = 0

    ALLOCATE (nelement(npe),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nelement",npe*int_size)
    nelement(:) = 0

    matrix_symmetry = ""

    CALL get_matrix_info(matrix=matrix,&
                         matrix_name=matrix_name,&
                         matrix_symmetry=matrix_symmetry,&
                         nblock_allocated=nblock(mype),&
                         nelement_allocated=nelement(mype),&
                         nblock_row=natom,&
                         nrow=nrow)

    CALL mp_sum(nblock,group)
    CALL mp_sum(nelement,group)

    nblock_sum = SUM(nblock)
    nelement_sum = SUM(nelement)

    IF (TRIM(matrix_symmetry) == "symmetric") THEN
      nblock_tot = natom*(natom + 1)/2
    ELSE
      nblock_tot = natom*natom
    END IF

    occupation = 100.0_dp*REAL(nblock_sum,dp)/REAL(nblock_tot,dp)

    IF (ionode) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
        "DISTRIBUTION OF THE "//TRIM(matrix_name)
      WRITE (UNIT=output_unit,FMT="(/,T3,A,/,/,(I9,T27,I10,T55,I10))")&
        "Process    Number of matrix blocks   Number of matrix elements",&
        (ipe-1,nblock(ipe),nelement(ipe),ipe=1,npe)
      WRITE (UNIT=output_unit,FMT="(/,T7,A3,T27,I10,T55,I10)")&
        "Sum",nblock_sum,nelement_sum
      WRITE (UNIT=output_unit,FMT="(/,T7,A3,T27,I10,A,F5.1,A,T55,I10,A,F5.1,A)")&
        " of",nblock_tot," (",occupation," % occupation)"
    END IF

!   *** Release work storage ***

    DEALLOCATE (nblock,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nblock")

    DEALLOCATE (nelement,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "nelement")

  END SUBROUTINE write_matrix_distribution

! *****************************************************************************

  SUBROUTINE write_sparse_matrix(sparse_matrix,before,after,qs_env,globenv,&
                                 first_row,last_row,first_col,last_col,scale)

!   Purpose: Print a spherical matrix of sparse_matrix_type.

!   History: - Creation (07.06.2000,MK)
!            - Allow for printing of a sub-matrix (01.07.2003,MK)

!   ***************************************************************************

    TYPE(real_matrix_type), POINTER          :: sparse_matrix
    INTEGER, INTENT(IN)                      :: before, after
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    INTEGER, INTENT(IN), OPTIONAL            :: first_row, last_row, &
                                                first_col, last_col
    REAL(dp), INTENT(IN), OPTIONAL           :: scale

    CHARACTER(LEN=*), PARAMETER :: routine_name = "write_sparse_matrix"

    CHARACTER(LEN=40)                        :: matrix_symmetry
    CHARACTER(LEN=80)                        :: matrix_name
    INTEGER                                  :: col1, col2, group, istat, &
                                                nsgf, output_unit, row1, row2
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: matrix
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set

!   ---------------------------------------------------------------------------

    group = globenv%group
    output_unit = globenv%scr

    CALL get_matrix_info(matrix=sparse_matrix,&
                         matrix_name=matrix_name,&
                         matrix_symmetry=matrix_symmetry)

    NULLIFY (matrix)

    CALL copy_local_sm_to_replicated_fm(sparse_matrix,matrix)

    IF (matrix_symmetry == "symmetric") THEN
      CALL symmetrize_matrix(matrix,"upper_to_lower")
    ELSE IF (matrix_symmetry == "antisymmetric") THEN
      CALL symmetrize_matrix(matrix,"anti_upper_to_lower")
    ELSE IF (matrix_symmetry /= "none") THEN
      CALL stop_program(routine_name,module_name,__LINE__,&
                        "Invalid matrix structure",globenv)
    END IF

    CALL mp_sum(matrix,group)

!   *** Get the matrix dimension and check the optional arguments ***

    CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,nsgf=nsgf)

    IF (PRESENT(first_row)) THEN
      row1 = MAX(1,first_row)
    ELSE
      row1 = 1
    END IF

    IF (PRESENT(last_row)) THEN
      row2 = MIN(nsgf,last_row)
    ELSE
      row2 = nsgf
    END IF

    IF (PRESENT(first_col)) THEN
      col1 = MAX(1,first_col)
    ELSE
      col1 = 1
    END IF

    IF (PRESENT(last_col)) THEN
      col2 = MIN(nsgf,last_col)
    ELSE
      col2 = nsgf
    END IF

    IF (PRESENT(scale)) THEN
      matrix=matrix*scale
    END IF

    CALL write_matrix(matrix,matrix_name,before,after,qs_env,globenv,&
                      row1,row2,col1,col2)

    IF (ASSOCIATED(matrix)) THEN
      DEALLOCATE (matrix,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "matrix")
    END IF

  END SUBROUTINE write_sparse_matrix

! *****************************************************************************

END MODULE qs_overlap
