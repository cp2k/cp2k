!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
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
!!     Matthias Krack (03.09.2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE qs_overlap

! *****************************************************************************

  USE kinds, ONLY: int_size,&
                   wp => dp,&
                   wp_size => dp_size

  USE ai_overlap,             ONLY: overlap
  USE atomic_kind_types,      ONLY: atomic_kind_type,&
                                    get_atomic_kind,&
                                    get_atomic_kind_set
  USE basis_set_types,        ONLY: get_gto_basis_set,&
                                    gto_basis_set_type
  USE global_types,           ONLY: global_environment_type
  USE message_passing,        ONLY: mp_sum
  USE mathlib,                ONLY: symmetrize_matrix
  USE orbital_pointers,       ONLY: ncoset,nso
  USE particle_types,         ONLY: particle_type
  USE qs_blacs,               ONLY: blacs_matrix_type,&
                                    copy_blacs_to_full_matrix,&
                                    get_blacs_matrix_info
  USE qs_neighbor_list_types, ONLY: first_list,&
                                    first_node,&
                                    get_neighbor_list,&
                                    get_neighbor_node,&
                                    neighbor_list_set_type,&
                                    neighbor_list_set_p_type,&
                                    neighbor_list_type,&
                                    neighbor_node_type,&
                                    next
  USE qs_environment_types,   ONLY: get_qs_env,&
                                    qs_environment_type,&
                                    set_qs_env
  USE sparse_matrix_types,    ONLY: add_block_node,&
                                    allocate_matrix,&
                                    copy_sparse_to_full_matrix,&
                                    deallocate_matrix,&
                                    get_block_node,&
                                    get_matrix_info,&
                                    real_matrix_p_type,&
                                    real_matrix_type
  USE termination,            ONLY: stop_memory
  USE timings,                ONLY: timeset,timestop

  IMPLICIT NONE

  PRIVATE

! *** Public subroutines ***

  PUBLIC :: build_overlap_matrix,&
            distribute_overlap_matrix,&
            write_blacs_matrix,&
            write_sparse_matrix

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE distribute_overlap_matrix(qs_env,globenv)

!   Purpose: Distribute the overlap integral matrix.

!   History: - Creation (26.06.2000, MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(qs_environment_type), INTENT(INOUT)  :: qs_env

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE distribute_overlap_matrix (MODULE qs_overlap)"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(gto_basis_set_type), POINTER :: orb_basis_set
    TYPE(neighbor_list_type), POINTER :: sab_orb_neighbor_list
    TYPE(neighbor_node_type), POINTER :: sab_orb_neighbor_node
    TYPE(real_matrix_p_type), POINTER :: s

    TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb
    TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set

    INTEGER :: ab,group,handle,iatom,ij,ikind,ipe,istat,jatom,jkind,jpe,mype,&
               n,natom,nkind,npack,npe,nsgf,nsgfa,nsgfb,output_unit,&
               sum_nblock_pe
    LOGICAL :: ionode

    INTEGER, DIMENSION(:), ALLOCATABLE :: nblock,nblock_pe,nelement,nelement_pe

!   ---------------------------------------------------------------------------

    CALL timeset("distribute_overlap_matrix","I","",handle)

    group = globenv%group
    ionode = globenv%ionode
    mype = globenv%mepos
    npe = globenv%num_pe
    output_unit = globenv%scr

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    s=s,&
                    sab_orb=sab_orb)

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)
    npack = natom*(natom + 1)/2

!   *** Allocate work storage ***

    ALLOCATE (nblock(npack),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"nblock",npack*int_size)
    nblock(:) = 0

    ALLOCATE (nblock_pe(0:npe-1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"nblock_pe",npe*int_size)
    nblock_pe(:) = 0

    ALLOCATE (nelement(npack),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"nelement",npack*int_size)
    nelement(:) = 0

    ALLOCATE (nelement_pe(0:npe-1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"nelement_pe",npe*int_size)
    nelement_pe(:) = 0

!   *** Allocate the overlap matrix (we start always from scratch) ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             nsgf=nsgf)

    CALL allocate_matrix(matrix_set=s,&
                         nblock_row=natom,&
                         nblock_col=natom,&
                         nrow=nsgf,&
                         ncol=nsgf,&
                         first_row=particle_set(:)%first_sgf,&
                         last_row=particle_set(:)%last_sgf,&
                         first_col=particle_set(:)%first_sgf,&
                         last_col=particle_set(:)%last_sgf,&
                         matrix_name="OVERLAP MATRIX",&
                         matrix_symmetry="symmetric")

!   *** Loop over all pairs of adjacent atoms and distribute work ***

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           orb_basis_set=orb_basis_set)
      IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE
      CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                             nsgf=nsgfa)

      DO jkind=ikind,nkind

        ab = ikind + jkind*(jkind - 1)/2

        IF (ikind == jkind) THEN
          nsgfb = nsgfa
        ELSE
          atomic_kind => atomic_kind_set(jkind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               orb_basis_set=orb_basis_set)
          IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE
          CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                 nsgf=nsgfb)
        END IF

        n = nsgfa*nsgfb

        sab_orb_neighbor_list => first_list(sab_orb(ab)%neighbor_list_set)

        DO WHILE (ASSOCIATED(sab_orb_neighbor_list))

          CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list,&
                                 atom=iatom)

          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list)

          DO WHILE (ASSOCIATED(sab_orb_neighbor_node))

            CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                   neighbor=jatom)

!           *** Map lower triangle blocks to the upper triangle ***

            IF (iatom <= jatom) THEN
              ij = iatom + jatom*(jatom - 1)/2
            ELSE
              ij = jatom + iatom*(iatom - 1)/2
            END IF

            nblock(ij) = nblock(ij) + 1
            nelement(ij) = nelement(ij) + n

            sab_orb_neighbor_node => next(sab_orb_neighbor_node)

          END DO

          sab_orb_neighbor_list => next(sab_orb_neighbor_list)

        END DO

      END DO

    END DO

!   *** Distribute the atom blocks ***

    DO iatom=1,natom
      DO jatom=iatom,natom

        ij = iatom + jatom*(jatom - 1)/2

        IF (nelement(ij) > 0) THEN

          ipe = 0

          DO jpe=0,npe-1
            IF (nelement_pe(jpe) < nelement_pe(ipe)) ipe = jpe
          END DO

          nblock_pe(ipe) = nblock_pe(ipe) + nblock(ij)
          nelement_pe(ipe) = nelement_pe(ipe) + nelement(ij)

          IF (ipe == mype) CALL add_block_node(matrix=s%matrix,&
                                               block_row=iatom,&
                                               block_col=jatom)

        END IF

      END DO
    END DO

!   *** Put the overlap matrix into the QS environment ***

    CALL set_qs_env(qs_env=qs_env,s=s)

!   *** Print the distribution of the overlap matrix ***

    IF (globenv%print%distribution) THEN

      IF (ionode) THEN
        WRITE (UNIT=output_unit,&
               FMT="(/,/,T2,A,/,/,T3,A,/,/,T5,A,/,/,(I6,8X,I8,8X,I10))")&
          "DISTRIBUTION OF THE OVERLAP MATRIX ELEMENTS",&
          "Image atoms included:",&
          "PE   Matrix blocks   Matrix elements",&
          (ipe,nblock_pe(ipe),nelement_pe(ipe),ipe=0,npe-1)
        WRITE (UNIT=output_unit,FMT="(/,T4,A3,8X,I8,8X,I10)")&
          "Sum",SUM(nblock_pe),SUM(nelement_pe)
      END IF

      nblock_pe(:) = 0
      nelement_pe(:) = 0

      CALL get_matrix_info(matrix=s%matrix,&
                           nblock_allocated=nblock_pe(mype),&
                           nelement_allocated=nelement_pe(mype))

      CALL mp_sum(nblock_pe,group)
      CALL mp_sum(nelement_pe,group)

      sum_nblock_pe = SUM(nblock_pe)

      IF (ionode) THEN
        WRITE (UNIT=output_unit,&
               FMT="(/,/,T3,A,/,/,T5,A,/,/,(I6,8X,I8,8X,I10))")&
          "Image atoms not included:",&
          "PE   Matrix blocks   Matrix elements",&
          (ipe,nblock_pe(ipe),nelement_pe(ipe),ipe=0,npe-1)
        WRITE (UNIT=output_unit,FMT="(/,T4,A3,8X,I8,8X,I10)")&
          "Sum",sum_nblock_pe,SUM(nelement_pe)
        WRITE (UNIT=output_unit,FMT="(/,T4,A3,8X,I8,A,F5.1,A)")&
          " of",npack," blocks in the full matrix (",&
          100.0_wp*REAL(sum_nblock_pe,wp)/REAL(npack,wp)," % occupation)"
      END IF

    END IF

!   *** Release work storage ***

    DEALLOCATE (nblock,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"nblock")

    DEALLOCATE (nblock_pe,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"nblock_pe")

    DEALLOCATE (nelement,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"nelement")

    DEALLOCATE (nelement_pe,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"nelement_pe")

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE distribute_overlap_matrix

! *****************************************************************************

  SUBROUTINE build_overlap_matrix(qs_env,globenv)

!   Purpose: Calculation of the overlap matrix over Cartesian Gaussian
!            functions.

!   History: - Creation (11.03.2002, MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(qs_environment_type), INTENT(INOUT)  :: qs_env

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE build_overlap_matrix (MODULE qs_overlap)"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(gto_basis_set_type), POINTER :: orb_basis_set
    TYPE(neighbor_list_type), POINTER :: sab_orb_neighbor_list
    TYPE(neighbor_node_type), POINTER :: sab_orb_neighbor_node
    TYPE(real_matrix_p_type), POINTER :: s

    TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb

    REAL(wp) :: dab,rab2
    INTEGER  :: ab,first_setb,handle,iatom,ikind,ineighbor,ipgf,iset,ishell,&
                istat,jatom,jkind,jpgf,jset,jshell,katom,maxco,maxsgf,ncoa,&
                ncob,nkind,nneighbor,nseta,nsetb,sgfa,sgfb

    REAL(wp), DIMENSION(3)                :: rab
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: sab,work

    REAL(wp), DIMENSION(:), POINTER   :: set_radius_a,set_radius_b
    INTEGER, DIMENSION(:), POINTER    :: la_max,la_min,lb_max,lb_min,npgfa,&
                                         npgfb,nsgfa,nsgfb
    REAL(wp), DIMENSION(:,:), POINTER :: rpgfa,rpgfb,s_block,sphi_a,sphi_b,&
                                         zeta,zetb
    INTEGER, DIMENSION(:,:), POINTER  :: first_sgfa,first_sgfb

!   ---------------------------------------------------------------------------

    CALL timeset("build_overlap_matrix","I","",handle)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    sab_orb=sab_orb,&
                    s=s)

!   *** Allocate work storage ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco,&
                             maxsgf=maxsgf)

    ALLOCATE (sab(maxco,maxco),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"sab",maxco*maxco*wp_size)

    ALLOCATE (work(maxco,maxsgf),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"work",maxco*maxsgf*wp_size)

    nkind = SIZE(atomic_kind_set)

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

      DO jkind=ikind,nkind

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

        ab = ikind + jkind*(jkind - 1)/2

        sab_orb_neighbor_list => first_list(sab_orb(ab)%neighbor_list_set)

        DO WHILE (ASSOCIATED(sab_orb_neighbor_list))

          CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list,&
                                 atom=iatom)

          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list)

          DO WHILE (ASSOCIATED(sab_orb_neighbor_node))

            CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                   neighbor=jatom,&
                                   r=rab(:))

            IF (iatom <= jatom) THEN
              CALL get_block_node(matrix=s%matrix,&
                                  block_row=iatom,&
                                  block_col=jatom,&
                                  block=s_block)
            ELSE
              CALL get_block_node(matrix=s%matrix,&
                                  block_row=jatom,&
                                  block_col=iatom,&
                                  block=s_block)
            END IF

!           *** Check, if the atomic block has to be ***
!           *** calculated by the current processor  ***

            IF (.NOT.ASSOCIATED(s_block)) THEN
              sab_orb_neighbor_node => next(sab_orb_neighbor_node)
              CYCLE
            END IF

            rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
            dab = SQRT(rab2)

            DO iset=1,nseta

              ncoa = npgfa(iset)*ncoset(la_max(iset))
              sgfa = first_sgfa(1,iset)

              IF (iatom == jatom) THEN
                first_setb = iset
              ELSE
                first_setb = 1
              END IF

              DO jset=first_setb,nsetb

                IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

                ncob = npgfb(jset)*ncoset(lb_max(jset))
                sgfb = first_sgfb(1,jset)

!               *** Calculate the primitive overlap integrals ***

                CALL overlap(la_max(iset),zeta(1:npgfa(iset),iset),&
                             rpgfa(1:npgfa(iset),iset),la_min(iset),&
                             lb_max(jset),zetb(1:npgfb(jset),jset),&
                             rpgfb(1:npgfb(jset),jset),lb_min(jset),&
                             rab,rab2,sab)

!               *** Contraction step (overlap matrix) ***

                CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                           1.0_wp,sab(1,1),SIZE(sab,1),&
                           sphi_b(1,sgfb),SIZE(sphi_b,1),&
                           0.0_wp,work(1,1),SIZE(work,1))

                IF (iatom <= jatom) THEN
                  CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                             1.0_wp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                             work(1,1),SIZE(work,1),&
                             1.0_wp,s_block(sgfa,sgfb),SIZE(s_block,1))
                ELSE
                  CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                             1.0_wp,work(1,1),SIZE(work,1),&
                             sphi_a(1,sgfa),SIZE(sphi_a,1),&
                             1.0_wp,s_block(sgfb,sgfa),SIZE(s_block,1))
                END IF

              END DO
            END DO

!           *** Symmetrize the diagonal blocks ***

            IF (iatom == jatom) THEN
              CALL symmetrize_matrix(s_block,"upper_to_lower")
            END IF

            sab_orb_neighbor_node => next(sab_orb_neighbor_node)

          END DO

          sab_orb_neighbor_list => next(sab_orb_neighbor_list)

        END DO

      END DO
    END DO

!   *** Release work storage ***

    IF (ALLOCATED(sab)) THEN
      DEALLOCATE (sab,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"sab")
    END IF

    IF (ALLOCATED(work)) THEN
      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"work")
    END IF

    CALL timestop(0.0_wp,handle)

!   *** Print the overlap matrix, if requested ***

    IF (globenv%print%overlap_matrix) THEN
      CALL write_sparse_matrix(s%matrix,4,6,qs_env,globenv)
    END IF

  END SUBROUTINE build_overlap_matrix

! *****************************************************************************

  SUBROUTINE write_blacs_matrix(blacs_matrix,before,after,context,qs_env,&
                                globenv)

!   Purpose: Write a spherical matrix of blacs_type.

!   History: - Creation (12.06.2001, MK)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: blacs_matrix
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(qs_environment_type), INTENT(IN)     :: qs_env
    INTEGER, INTENT(IN)                       :: after,before,context

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE write_blacs_matrix (MODULE qs_overlap)"

!   *** Local variables ***

    TYPE(gto_basis_set_type), POINTER :: orb_basis_set

    CHARACTER(LEN=60) :: name
    CHARACTER(LEN=40) :: symmetry
    CHARACTER(LEN=25) :: fmtstr1
    CHARACTER(LEN=35) :: fmtstr2
    CHARACTER(LEN=2)  :: element_symbol
    INTEGER           :: from,iatom,icol,irow,iset,isgf,iso,ishell,istat,jcol,&
                         l,left,natom,ncol,nset,nsgf,output_unit,right,to,width

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set

    CHARACTER(LEN=6), DIMENSION(:), POINTER :: sgf_symbol
    INTEGER, DIMENSION(:), POINTER          :: nshell
    REAL(wp), DIMENSION(:,:), POINTER       :: matrix
    INTEGER, DIMENSION(:,:), POINTER        :: lshell

!   ---------------------------------------------------------------------------

    output_unit = globenv%scr

    CALL get_blacs_matrix_info(blacs_matrix,name=name)

    NULLIFY (matrix)
    CALL copy_blacs_to_full_matrix(blacs_matrix,matrix,context,globenv)

    IF (.NOT.globenv%ionode) RETURN

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set)

    natom = SIZE(particle_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,nsgf=nsgf)

!   *** Definition of the variable formats ***

    fmtstr1 = "(/,T2,23X,  (  X,I5,  X))"
    fmtstr2 = "(T2,2I5,2X,A2,1X,A8,   (1X,F  .  ))"

!   *** Write headline ***

    WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") TRIM(name)

!   *** Write the variable format strings ***

    width = before + after + 3
    ncol = INT(56/width)

    right = MAX((after-2),1)
    left =  width - right - 5

    WRITE (UNIT=fmtstr1(11:12),FMT="(I2)") ncol
    WRITE (UNIT=fmtstr1(14:15),FMT="(I2)") left
    WRITE (UNIT=fmtstr1(21:22),FMT="(I2)") right

    WRITE (UNIT=fmtstr2(22:23),FMT="(I2)") ncol
    WRITE (UNIT=fmtstr2(29:30),FMT="(I2)") width - 1
    WRITE (UNIT=fmtstr2(32:33),FMT="(I2)") after

!   *** Write the matrix in the selected format ***

    DO icol=1,nsgf,ncol
      from = icol
      to = MIN((from+ncol-1),nsgf)
      WRITE (UNIT=output_unit,FMT=fmtstr1) (jcol,jcol=from,to)
      irow = 1
      DO iatom=1,natom
        IF (iatom /= 1) WRITE (UNIT=output_unit,FMT="(A)")
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
              WRITE (UNIT=output_unit,FMT=fmtstr2)&
                irow,iatom,element_symbol,sgf_symbol(isgf),&
                (matrix(irow,jcol),jcol=from,to)
              isgf = isgf + 1
              irow = irow + 1
            END DO
          END DO
        END DO
      END DO
    END DO

    WRITE (UNIT=output_unit,FMT="(/)")

!   *** Release work storage ***

    DEALLOCATE (matrix,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"matrix")

  END SUBROUTINE write_blacs_matrix

! *****************************************************************************

  SUBROUTINE write_sparse_matrix(sparse_matrix,before,after,qs_env,globenv)

!   Purpose: Write a spherical matrix.

!   History: - Creation (07.06.2000, MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(qs_environment_type), INTENT(IN)     :: qs_env
    TYPE(real_matrix_type), POINTER           :: sparse_matrix
    INTEGER, INTENT(IN)                       :: after,before

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE write_sparse_matrix (MODULE qs_overlap)"

!   *** Local variables ***

    TYPE(gto_basis_set_type), POINTER :: orb_basis_set

    CHARACTER(LEN=60) :: name
    CHARACTER(LEN=40) :: symmetry
    CHARACTER(LEN=25) :: fmtstr1
    CHARACTER(LEN=35) :: fmtstr2
    CHARACTER(LEN=2)  :: element_symbol
    INTEGER           :: from,iatom,icol,irow,iset,isgf,iso,ishell,istat,jcol,&
                         l,left,natom,ncol,nset,nsgf,output_unit,right,to,width

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set

    CHARACTER(LEN=6), DIMENSION(:), POINTER :: sgf_symbol
    INTEGER, DIMENSION(:), POINTER          :: nshell
    REAL(wp), DIMENSION(:,:), POINTER       :: matrix
    INTEGER, DIMENSION(:,:), POINTER        :: lshell

!   ---------------------------------------------------------------------------

    output_unit = globenv%scr

    CALL get_matrix_info(matrix=sparse_matrix,&
                         matrix_name=name,&
                         matrix_symmetry=symmetry)

    NULLIFY (matrix)
    CALL copy_sparse_to_full_matrix(sparse_matrix,matrix)

    IF (symmetry == "symmetric") THEN
      CALL symmetrize_matrix(matrix,"upper_to_lower")
    ELSE IF (symmetry == "antisymmetric") THEN
      CALL symmetrize_matrix(matrix,"anti_upper_to_lower")
    END IF

    CALL mp_sum(matrix,globenv%group)

    IF (.NOT.globenv%ionode) RETURN

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set)

    natom = SIZE(particle_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,nsgf=nsgf)

!   *** Definition of the variable formats ***

    fmtstr1 = "(/,T2,23X,  (  X,I5,  X))"
    fmtstr2 = "(T2,2I5,2X,A2,1X,A8,   (1X,F  .  ))"

!   *** Write headline ***

    WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") TRIM(name)

!   *** Write the variable format strings ***

    width = before + after + 3
    ncol = INT(56/width)

    right = MAX((after-2),1)
    left =  width - right - 5

    WRITE (UNIT=fmtstr1(11:12),FMT="(I2)") ncol
    WRITE (UNIT=fmtstr1(14:15),FMT="(I2)") left
    WRITE (UNIT=fmtstr1(21:22),FMT="(I2)") right

    WRITE (UNIT=fmtstr2(22:23),FMT="(I2)") ncol
    WRITE (UNIT=fmtstr2(29:30),FMT="(I2)") width - 1
    WRITE (UNIT=fmtstr2(32:33),FMT="(I2)") after

!   *** Write the matrix in the selected format ***

    DO icol=1,nsgf,ncol
      from = icol
      to = MIN((from+ncol-1),nsgf)
      WRITE (UNIT=output_unit,FMT=fmtstr1) (jcol,jcol=from,to)
      irow = 1
      DO iatom=1,natom
        IF (iatom /= 1) WRITE (UNIT=output_unit,FMT="(A)")
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
              WRITE (UNIT=output_unit,FMT=fmtstr2)&
                irow,iatom,element_symbol,sgf_symbol(isgf),&
                (matrix(irow,jcol),jcol=from,to)
              isgf = isgf + 1
              irow = irow + 1
            END DO
          END DO
        END DO
      END DO
    END DO

    WRITE (UNIT=output_unit,FMT="(/)")

!   *** Release work storage ***

    DEALLOCATE (matrix,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"matrix")

  END SUBROUTINE write_sparse_matrix

! *****************************************************************************

END MODULE qs_overlap
