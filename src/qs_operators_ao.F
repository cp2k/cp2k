!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      created 07.2005
!> \author MI (07.2005)
! *****************************************************************************
MODULE qs_operators_ao
  USE ai_moments,                      ONLY: diff_momop,&
                                             diffop,&
                                             moment
  USE ai_os_rr,                        ONLY: os_rr_ovlp
  USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                             gto_basis_set_type
  USE block_p_types,                   ONLY: block_p_type
  USE cell_methods,                    ONLY: cell_type,&
                                             pbc
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_get_block_p,&
                                             cp_dbcsr_get_matrix_type,&
                                             cp_dbcsr_p_type,&
                                             dbcsr_type_antisymmetric,&
                                             dbcsr_type_no_symmetry
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE orbital_pointers,                ONLY: coset,&
                                             init_orbital_pointers,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             get_qs_kind_set,&
                                             qs_kind_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_operators_ao'

! *** Public subroutines ***

  PUBLIC :: p_xyz_ao, rRc_xyz_ao, rRc_xyz_der_ao
  PUBLIC :: build_lin_mom_matrix, build_ang_mom_matrix

CONTAINS

! *****************************************************************************
!> \brief   Calculation of the linear momemtum matrix over
!>          Cartesian Gaussian functions.
!> \param qs_env ...
!> \param matrix ...
!> \param error ...
!> \date    27.02.2009
!> \author  VW
!> \version 1.0
! *****************************************************************************
  SUBROUTINE build_lin_mom_matrix(qs_env,matrix,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_lin_mom_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iatom, icol, ikind, inode, irow, iset, istat, &
      jatom, jkind, jset, last_jatom, ldai, maxco, maxlgto, maxsgf, natom, &
      ncoa, ncob, neighbor_list_id, nkind, nseta, nsetb, sgfa, sgfb
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure, found, new_atom_b
    REAL(KIND=dp)                            :: dab, rab2
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: intab, rr_work
    REAL(KIND=dp), DIMENSION(3)              :: rab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                zeta, zetb
    TYPE(block_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: integral
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_kind_type), POINTER              :: qs_kind

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY(cell,sab_orb,qs_kind_set,particle_set,para_env)
    NULLIFY(logger)

    logger => cp_error_get_logger(error)

    CALL get_qs_env(qs_env=qs_env,&
                    qs_kind_set=qs_kind_set,&
                    particle_set=particle_set,&
                    neighbor_list_id=neighbor_list_id,&
                    para_env=para_env,&
                    sab_orb=sab_orb,&
                    cell=cell,&
                    error=error)

    nkind = SIZE(qs_kind_set)
    natom = SIZE(particle_set)

!   *** Allocate work storage ***

    CALL get_qs_kind_set(qs_kind_set=qs_kind_set,&
                             maxco=maxco,&
                             maxlgto=maxlgto,&
                             maxsgf=maxsgf)

    ldai = ncoset(maxlgto+1)
    CALL init_orbital_pointers(ldai)

    ALLOCATE(rr_work(ldai,ldai,3),intab(maxco,maxco,3),work(maxco,maxsgf),integral(3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    rr_work(:,:,:) = 0.0_dp
    intab(:,:,:) = 0.0_dp
    work(:,:) = 0.0_dp

    ALLOCATE (basis_set_list(nkind),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ikind=1,nkind
      qs_kind => qs_kind_set(ikind)
      CALL get_qs_kind(qs_kind=qs_kind,orb_basis_set=basis_set_a)
      IF (ASSOCIATED(basis_set_a)) THEN
        basis_set_list(ikind)%gto_basis_set => basis_set_a
      ELSE
        NULLIFY(basis_set_list(ikind)%gto_basis_set)
      END IF
    END DO
    CALL neighbor_list_iterator_create(nl_iterator,sab_orb)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,inode=inode,&
                              iatom=iatom,jatom=jatom,r=rab)
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
          IF (iatom <= jatom) THEN
             irow = iatom
             icol = jatom
          ELSE
             irow = jatom
             icol = iatom
          END IF

          DO i=1,3
             NULLIFY(integral(i)%block)
             CALL cp_dbcsr_get_block_p(matrix=matrix(i)%matrix,&
                  row=irow,col=icol,BLOCK=integral(i)%block,found=found)
             CPPostcondition(ASSOCIATED(INTEGRAL(i)%block),cp_failure_level,routineP,error,failure)
          ENDDO
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

             ! *** Calculate the primitive fermi contact integrals ***

             CALL lin_mom(la_max(iset),la_min(iset),npgfa(iset),&
                          rpgfa(:,iset),zeta(:,iset),&
                          lb_max(jset),lb_min(jset),npgfb(jset),&
                          rpgfb(:,jset),zetb(:,jset),&
                          rab,intab,SIZE(rr_work,1),rr_work,error)

             ! *** Contraction step ***

             DO i=1,3

                CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                           1.0_dp,intab(1,1,i),SIZE(intab,1),&
                           sphi_b(1,sgfb),SIZE(sphi_b,1),&
                           0.0_dp,work(1,1),SIZE(work,1))

                IF (iatom <= jatom) THEN

                   CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                              1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                              work(1,1),SIZE(work,1),&
                              1.0_dp,integral(i)%block(sgfa,sgfb),&
                              SIZE(integral(i)%block,1))

                ELSE

                   CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                             -1.0_dp,work(1,1),SIZE(work,1),&
                              sphi_a(1,sgfa),SIZE(sphi_a,1),&
                              1.0_dp,integral(i)%block(sgfb,sgfa),&
                              SIZE(integral(i)%block,1))

                ENDIF

             ENDDO

          ENDDO

       ENDDO

    END DO
    CALL neighbor_list_iterator_release(nl_iterator)

    ! *** Release work storage ***

    DEALLOCATE(intab,work,integral,basis_set_list,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

!   *** Print the spin orbit matrix, if requested ***

    !IF (BTEST(cp_print_key_should_output(logger%iter_info,&
    !     qs_env%input,"DFT%PRINT%AO_MATRICES/LINEAR_MOMENTUM",error=error),cp_p_file)) THEN
    !   iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/LINEA_MOMENTUM",&
    !        extension=".Log",error=error)
    !   CALL cp_dbcsr_write_sparse_matrix(matrix(1)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
    !   CALL cp_dbcsr_write_sparse_matrix(matrix(2)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
    !   CALL cp_dbcsr_write_sparse_matrix(matrix(3)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
    !   CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
    !        "DFT%PRINT%AO_MATRICES/LINEAR_MOMENTUM", error=error)
    !END IF

    CALL timestop(handle)

  END SUBROUTINE build_lin_mom_matrix


! *****************************************************************************
!> \brief   Calculation of the primitive paramagnetic spin orbit integrals over
!>          Cartesian Gaussian-type functions.
!> \param la_max ...
!> \param la_min ...
!> \param npgfa ...
!> \param rpgfa ...
!> \param zeta ...
!> \param lb_max ...
!> \param lb_min ...
!> \param npgfb ...
!> \param rpgfb ...
!> \param zetb ...
!> \param rab ...
!> \param intab ...
!> \param ldrr ...
!> \param rr ...
!> \param error ...
!> \date    02.03.2009
!> \author  VW
!> \version 1.0
! *****************************************************************************
  SUBROUTINE lin_mom(la_max,la_min,npgfa,rpgfa,zeta,lb_max,lb_min,npgfb,rpgfb,zetb,&
                     rab,intab,ldrr,rr,error)
    INTEGER, INTENT(IN)                      :: la_max, la_min, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfa, zeta
    INTEGER, INTENT(IN)                      :: lb_max, lb_min, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfb, zetb
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: intab
    INTEGER, INTENT(IN)                      :: ldrr
    REAL(dp), &
      DIMENSION(0:ldrr-1, 0:ldrr-1, 3), &
      INTENT(INOUT)                          :: rr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'lin_mom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ax, ay, az, bx, by, bz, coa, &
                                                cob, i, ipgf, j, jpgf, la, &
                                                lb, ma, mb, na, nb
    REAL(dp)                                 :: dab, dumx, dumy, dumz, f0, &
                                                rab2, xhi, zet
    REAL(dp), DIMENSION(3)                   :: rap, rbp

! *** Calculate the distance of the centers a and c ***

    rab2 = rab(1)**2+rab(2)**2+rab(3)**2
    dab = SQRT(rab2)

    ! *** Loop over all pairs of primitive Gaussian-type functions ***

    na = 0

    DO ipgf=1,npgfa

       nb = 0

       DO jpgf=1,npgfb

          ! *** Screening ***

          IF(rpgfa(ipgf) + rpgfb(jpgf) < dab) THEN
             DO j=nb+1,nb+ncoset(lb_max)
                DO i=na+1,na+ncoset(la_max)
                   intab(i,j,1) = 0.0_dp
                   intab(i,j,2) = 0.0_dp
                   intab(i,j,3) = 0.0_dp
                ENDDO
             ENDDO
             nb = nb + ncoset(lb_max)
             CYCLE
          ENDIF

          ! *** Calculate some prefactors ***
          zet =  zeta(ipgf)+zetb(jpgf)
          xhi =  zeta(ipgf)*zetb(jpgf)/zet
          rap =  zetb(jpgf)*rab/zet
          rbp = -zeta(ipgf)*rab/zet

          f0 = (pi/zet)**(1.5_dp)*EXP(-xhi*rab2)

          ! *** Calculate the reccurence relation ***

          CALL os_rr_ovlp(rap,la_max+1,rbp,lb_max,zet,ldrr,rr)

          ! *** Calculate the primitive linear momentum integrals ***
          DO lb = lb_min,lb_max
          DO bx = 0,lb
          DO by = 0,lb-bx
             bz = lb - bx - by
             cob = coset(bx,by,bz)
             mb = nb + cob
             DO la = la_min,la_max
             DO ax = 0,la
             DO ay = 0,la-ax
                az = la - ax - ay
                coa = coset(ax,ay,az)
                ma = na + coa
                !
                !
                ! (a|p_x|b) = 2*a*(a+1x|b) - N_x(a)*(a-1_x|b)
                dumx = 2.0_dp * zeta(ipgf) * rr(ax+1,bx,1)
                IF(ax.GT.0) dumx = dumx - REAL(ax,dp) * rr(ax-1,bx,1)
                intab(ma,mb,1) = f0 * dumx * rr(ay,by,2) * rr(az,bz,3)
                !
                ! (a|p_y|b)
                dumy = 2.0_dp * zeta(ipgf) * rr(ay+1,by,2)
                IF(ay.GT.0) dumy = dumy - REAL(ay,dp) * rr(ay-1,by,2)
                intab(ma,mb,2) = f0 * rr(ax,bx,1) * dumy * rr(az,bz,3)
                !
                ! (a|p_z|b)
                dumz = 2.0_dp * zeta(ipgf) * rr(az+1,bz,3)
                IF(az.GT.0) dumz = dumz - REAL(az,dp) * rr(az-1,bz,3)
                intab(ma,mb,3) = f0 * rr(ax,bx,1) * rr(ay,by,2) * dumz
                !
             ENDDO
             ENDDO
             ENDDO !la

          ENDDO
          ENDDO
          ENDDO !lb

          nb = nb + ncoset(lb_max)

       ENDDO

       na = na + ncoset(la_max)

    ENDDO

  END SUBROUTINE lin_mom

! *****************************************************************************
!> \brief   Calculation of the angular momemtum matrix over
!>          Cartesian Gaussian functions.
!> \param qs_env ...
!> \param matrix ...
!> \param rc ...
!> \param error ...
!> \date    27.02.2009
!> \author  VW
!> \version 1.0
! *****************************************************************************

  SUBROUTINE build_ang_mom_matrix(qs_env,matrix,rc,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_ang_mom_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iatom, icol, ikind, inode, irow, iset, istat, &
      jatom, jkind, jset, last_jatom, ldai, maxco, maxlgto, maxsgf, natom, &
      ncoa, ncob, neighbor_list_id, nkind, nseta, nsetb, sgfa, sgfb
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure, found, new_atom_b
    REAL(KIND=dp)                            :: dab, rab2
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: intab, rr_work
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rac, rbc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                zeta, zetb
    TYPE(block_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: integral
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_all
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_kind_type), POINTER              :: qs_kind

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY(cell,sab_all,qs_kind_set,particle_set,para_env)
    NULLIFY(logger)

    logger => cp_error_get_logger(error)

    CALL get_qs_env(qs_env=qs_env,&
                    qs_kind_set=qs_kind_set,&
                    particle_set=particle_set,&
                    neighbor_list_id=neighbor_list_id,&
                    para_env=para_env,&
                    sab_all=sab_all,&
                    cell=cell,&
                    error=error)

    nkind = SIZE(qs_kind_set)
    natom = SIZE(particle_set)

!   *** Allocate work storage ***

    CALL get_qs_kind_set(qs_kind_set=qs_kind_set,&
                             maxco=maxco,&
                             maxlgto=maxlgto,&
                             maxsgf=maxsgf)

    ldai = ncoset(maxlgto+1)
    CALL init_orbital_pointers(ldai)

    ALLOCATE(rr_work(ldai,ldai,3),intab(maxco,maxco,3),work(maxco,maxsgf),integral(3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    rr_work(:,:,:) = 0.0_dp
    intab(:,:,:) = 0.0_dp
    work(:,:) = 0.0_dp

    ALLOCATE (basis_set_list(nkind),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ikind=1,nkind
      qs_kind => qs_kind_set(ikind)
      CALL get_qs_kind(qs_kind=qs_kind,orb_basis_set=basis_set_a)
      IF (ASSOCIATED(basis_set_a)) THEN
        basis_set_list(ikind)%gto_basis_set => basis_set_a
      ELSE
        NULLIFY(basis_set_list(ikind)%gto_basis_set)
      END IF
    END DO
    CALL neighbor_list_iterator_create(nl_iterator,sab_all)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,inode=inode,&
                              iatom=iatom,jatom=jatom,r=rab)
       basis_set_a => basis_set_list(ikind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_a)) CYCLE
       basis_set_b => basis_set_list(jkind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_b)) CYCLE
       ra = pbc(particle_set(iatom)%r,cell)
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
          !IF (iatom <= jatom) THEN
             irow = iatom
             icol = jatom
          !ELSE
          !   irow = jatom
          !   icol = iatom
          !END IF

          DO i=1,3
             NULLIFY(INTEGRAL(i)%block)
             CALL cp_dbcsr_get_block_p(matrix=matrix(i)%matrix,&
                  row=irow,col=icol,BLOCK=INTEGRAL(i)%block,found=found)
             CPPostcondition(ASSOCIATED(INTEGRAL(i)%block),cp_failure_level,routineP,error,failure)
          ENDDO
       END IF

       rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
       dab = SQRT(rab2)

       DO iset=1,nseta

          ncoa = npgfa(iset)*ncoset(la_max(iset))
          sgfa = first_sgfa(1,iset)

          DO jset=1,nsetb

             IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

             !IF(PRESENT(wancen)) THEN
             !   rc = wancen
             rac = pbc(rc,ra,cell)
             rbc = rac + rab
             !ELSE
             !   rc(1:3) = rb(1:3)
             !   rac(1:3) = -rab(1:3)
             !   rbc(1:3) = 0.0_dp
             !ENDIF

             ncob = npgfb(jset)*ncoset(lb_max(jset))
             sgfb = first_sgfb(1,jset)

             ! *** Calculate the primitive angular momentum integrals ***

             CALL ang_mom(la_max(iset),la_min(iset),npgfa(iset),&
                          rpgfa(:,iset),zeta(:,iset),&
                          lb_max(jset),lb_min(jset),npgfb(jset),&
                          rpgfb(:,jset),zetb(:,jset),&
                          rab,rac,intab,SIZE(rr_work,1),rr_work,error)

             ! *** Contraction step ***

             DO i=1,3

                CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                           1.0_dp,intab(1,1,i),SIZE(intab,1),&
                           sphi_b(1,sgfb),SIZE(sphi_b,1),&
                           0.0_dp,work(1,1),SIZE(work,1))

                !IF (iatom <= jatom) THEN

                   CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                              1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                              work(1,1),SIZE(work,1),&
                              1.0_dp,integral(i)%block(sgfa,sgfb),&
                              SIZE(integral(i)%block,1))


                !ELSE
                !
                !   CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                !              -1.0_dp,work(1,1),SIZE(work,1),&
                !              sphi_a(1,sgfa),SIZE(sphi_a,1),&
                !              1.0_dp,integral(i)%block(sgfb,sgfa),&
                !              SIZE(integral(i)%block,1))
                !
                !ENDIF

             ENDDO

          ENDDO

       ENDDO

    END DO
    CALL neighbor_list_iterator_release(nl_iterator)

    ! *** Release work storage ***

    DEALLOCATE(intab,work,integral,basis_set_list,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

!   *** Print the spin orbit matrix, if requested ***

    !IF (BTEST(cp_print_key_should_output(logger%iter_info,&
    !     qs_env%input,"DFT%PRINT%AO_MATRICES/ANGULAR_MOMENTUM",error=error),cp_p_file)) THEN
    !   iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/ANGULAR_MOMENTUM",&
    !        extension=".Log",error=error)
    !   CALL cp_dbcsr_write_sparse_matrix(matrix(1)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
    !   CALL cp_dbcsr_write_sparse_matrix(matrix(2)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
    !   CALL cp_dbcsr_write_sparse_matrix(matrix(3)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
    !   CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
    !        "DFT%PRINT%AO_MATRICES/ANGULAR_MOMENTUM", error=error)
    !END IF

    CALL timestop(handle)

  END SUBROUTINE build_ang_mom_matrix

! *****************************************************************************
!> \brief   Calculation of the primitive paramagnetic spin orbit integrals over
!>          Cartesian Gaussian-type functions.
!> \param la_max ...
!> \param la_min ...
!> \param npgfa ...
!> \param rpgfa ...
!> \param zeta ...
!> \param lb_max ...
!> \param lb_min ...
!> \param npgfb ...
!> \param rpgfb ...
!> \param zetb ...
!> \param rab ...
!> \param rac ...
!> \param intab ...
!> \param ldrr ...
!> \param rr ...
!> \param error ...
!> \date    02.03.2009
!> \author  VW
!> \version 1.0
! *****************************************************************************
  SUBROUTINE ang_mom(la_max,la_min,npgfa,rpgfa,zeta,lb_max,lb_min,npgfb,rpgfb,zetb,&
                     rab,rac,intab,ldrr,rr,error)
    INTEGER, INTENT(IN)                      :: la_max, la_min, npgfa
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfa, zeta
    INTEGER, INTENT(IN)                      :: lb_max, lb_min, npgfb
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: rpgfb, zetb
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab, rac
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: intab
    INTEGER, INTENT(IN)                      :: ldrr
    REAL(dp), &
      DIMENSION(0:ldrr-1, 0:ldrr-1, 3), &
      INTENT(INOUT)                          :: rr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ang_mom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ax, ay, az, bx, by, bz, coa, &
                                                cob, i, ipgf, j, jpgf, la, &
                                                lb, ma, mb, na, nb
    REAL(dp)                                 :: dab, dumx, dumy, dumz, f0, &
                                                rab2, xhi, zet
    REAL(dp), DIMENSION(3)                   :: rap, rbp

! *** Calculate the distance of the centers a and c ***

    rab2 = rab(1)**2+rab(2)**2+rab(3)**2
    dab = SQRT(rab2)

    ! *** Loop over all pairs of primitive Gaussian-type functions ***

    na = 0

    DO ipgf=1,npgfa

       nb = 0

       DO jpgf=1,npgfb

          ! *** Screening ***

          IF(rpgfa(ipgf) + rpgfb(jpgf) < dab) THEN
             DO j=nb+1,nb+ncoset(lb_max)
                DO i=na+1,na+ncoset(la_max)
                   intab(i,j,1) = 0.0_dp
                   intab(i,j,2) = 0.0_dp
                   intab(i,j,3) = 0.0_dp
                ENDDO
             ENDDO
             nb = nb + ncoset(lb_max)
             CYCLE
          ENDIF

          ! *** Calculate some prefactors ***
          zet =  zeta(ipgf)+zetb(jpgf)
          xhi =  zeta(ipgf)*zetb(jpgf)/zet
          rap =  zetb(jpgf)*rab/zet
          rbp = -zeta(ipgf)*rab/zet

          f0 = (pi/zet)**(1.5_dp)*EXP(-xhi*rab2)

          ! *** Calculate the reccurence relation ***

          CALL os_rr_ovlp(rap,la_max+1,rbp,lb_max,zet,ldrr,rr)

          ! *** Calculate the primitive Fermi contact integrals ***

          DO lb = lb_min,lb_max
          DO bx = 0,lb
          DO by = 0,lb-bx
             bz = lb - bx - by
             cob = coset(bx,by,bz)
             mb = nb + cob
             DO la = la_min,la_max
             DO ax = 0,la
             DO ay = 0,la-ax
                az = la - ax - ay
                coa = coset(ax,ay,az)
                ma = na + coa
                !
                dumx = -2.0_dp * zeta(ipgf) * rr(ax+1,bx,1)
                dumy = -2.0_dp * zeta(ipgf) * rr(ay+1,by,2)
                dumz = -2.0_dp * zeta(ipgf) * rr(az+1,bz,3)
                IF(ax.GT.0) dumx = dumx + REAL(ax,dp) * rr(ax-1,bx,1)
                IF(ay.GT.0) dumy = dumy + REAL(ay,dp) * rr(ay-1,by,2)
                IF(az.GT.0) dumz = dumz + REAL(az,dp) * rr(az-1,bz,3)
                !
                ! (a|l_z|b)
                intab(ma,mb,1) = - f0 * rr(ax,bx,1) * ( &
                     &  ( rr(ay+1,by,2) + rac(2) * rr(ay,by,2) ) * dumz &
                     & -( rr(az+1,bz,3) + rac(3) * rr(az,bz,3) ) * dumy )
                !
                ! (a|l_y|b)
                intab(ma,mb,2) = - f0 * rr(ay,by,2) * ( &
                     &  ( rr(az+1,bz,3) + rac(3) * rr(az,bz,3) ) * dumx &
                     & -( rr(ax+1,bx,1) + rac(1) * rr(ax,bx,1) ) * dumz )
                !
                ! (a|l_z|b)
                intab(ma,mb,3) = - f0 * rr(az,bz,3) * ( &
                     &  ( rr(ax+1,bx,1) + rac(1) * rr(ax,bx,1) ) * dumy &
                     & -( rr(ay+1,by,2) + rac(2) * rr(ay,by,2) ) * dumx )
                !
             ENDDO
             ENDDO
             ENDDO !la

          ENDDO
          ENDDO
          ENDDO !lb

          nb = nb + ncoset(lb_max)

       ENDDO

       na = na + ncoset(la_max)

    ENDDO

  END SUBROUTINE ang_mom


! *****************************************************************************
!> \brief Calculation of the components of the dipole operator in the velocity form
!>      The elements of the  sparse matrices are the integrals in the
!>      basis functions
!> \param op matrix representation of the p operator
!>               calculated in terms of the contracted basis functions
!> \param qs_env enviroment for the lists and the basis sets
!> \param minimum_image take into account only the first neighbors in the lists
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2005 created [MI]
!> \author MI
! *****************************************************************************

  SUBROUTINE p_xyz_ao(op,qs_env,minimum_image,error)

    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: op
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN), OPTIONAL            :: minimum_image
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'p_xyz_ao', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iatom, icol, ikind, inode, irow, iset, istat, &
      jatom, jkind, jset, last_jatom, ldab, ldsa, ldsb, ldwork, maxl, ncoa, &
      ncob, nkind, nseta, nsetb, sgfa, sgfb
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure, found, &
                                                my_minimum_image, new_atom_b
    REAL(KIND=dp)                            :: alpha, dab, Lxo2, Lyo2, Lzo2, &
                                                rab2
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                work, zeta, zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: difab
    TYPE(block_p_type), DIMENSION(:), &
      POINTER                                :: op_dip
    TYPE(cell_type), POINTER                 :: cell
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_kind_type), POINTER              :: qs_kind

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY (qs_kind, qs_kind_set)
    NULLIFY (cell, particle_set)
    NULLIFY (sab_orb)
    NULLIFY (difab, op_dip, work)
    NULLIFY (la_max, la_min, lb_max, lb_min, npgfa, npgfb, nsgfa, nsgfb)
    NULLIFY (set_radius_a,set_radius_b,rpgfa, rpgfb, sphi_a,sphi_b,zeta, zetb )

    CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set,&
         cell=cell,particle_set=particle_set,&
         sab_orb=sab_orb,error=error)

    nkind = SIZE(qs_kind_set)

    CALL get_qs_kind_set(qs_kind_set=qs_kind_set,&
                             maxco=ldwork, maxlgto=maxl )

    my_minimum_image = .FALSE.
    IF(PRESENT(minimum_image)) THEN
       my_minimum_image=minimum_image
       Lxo2 = SQRT ( SUM ( cell % hmat ( :, 1 ) ** 2 ) )/2.0_dp
       Lyo2 = SQRT ( SUM ( cell % hmat ( :, 2 ) ** 2 ) )/2.0_dp
       Lzo2 = SQRT ( SUM ( cell % hmat ( :, 3 ) ** 2 ) )/2.0_dp
    END IF

    ldab = ldwork

    ALLOCATE(difab(ldab,ldab,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    difab(1:ldab,1:ldab,1:3) = 0.0_dp
    ALLOCATE(work(ldwork,ldwork),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    work(1:ldwork,1:ldwork) = 0.0_dp
    ALLOCATE (op_dip(3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO i = 1, 3
       NULLIFY (op_dip(i)%block)
    END DO

    ALLOCATE (basis_set_list(nkind),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ikind=1,nkind
      qs_kind => qs_kind_set(ikind)
      CALL get_qs_kind(qs_kind=qs_kind,orb_basis_set=basis_set_a)
      IF (ASSOCIATED(basis_set_a)) THEN
        basis_set_list(ikind)%gto_basis_set => basis_set_a
      ELSE
        NULLIFY(basis_set_list(ikind)%gto_basis_set)
      END IF
    END DO
    CALL neighbor_list_iterator_create(nl_iterator,sab_orb)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,inode=inode,&
                              iatom=iatom,jatom=jatom,r=rab)
       basis_set_a => basis_set_list(ikind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_a)) CYCLE
       basis_set_b => basis_set_list(jkind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_b)) CYCLE
       ra = pbc(particle_set(iatom)%r,cell)
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

       IF(inode==1) THEN
         last_jatom = 0
         alpha = 1.0_dp
       END IF
       ldsa = SIZE(sphi_a,1)
       ldsb = SIZE(sphi_b,1)


       IF(my_minimum_image) THEN
         IF(ABS(rab(1)) > Lxo2 .OR. ABS(rab(2)) > Lyo2 .OR. ABS(rab(3)) > Lzo2) CYCLE
       END IF

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
           alpha = 1.0_dp
         ELSE
           irow = jatom
           icol = iatom
           IF(cp_dbcsr_get_matrix_type(op(1)%matrix)==dbcsr_type_antisymmetric) THEN
           !IF(op(1)%matrix%symmetry=="antisymmetric") THEN
             alpha = -1.0_dp
           END IF
         END IF

         DO i = 1,3
            NULLIFY(op_dip(i)%block)
            CALL cp_dbcsr_get_block_p(matrix=op(i)%matrix,&
                 row=irow,col=icol,block=op_dip(i)%block,found=found)
            CPPostcondition(ASSOCIATED(op_dip(i)%block),cp_failure_level,routineP,error,failure)
         END DO
       END IF  ! new_atom_b
       rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
       dab = SQRT(rab2)

       DO iset=1,nseta

         ncoa = npgfa(iset)*ncoset(la_max(iset))
         sgfa = first_sgfa(1,iset)

         DO jset=1,nsetb

           ncob = npgfb(jset)*ncoset(lb_max(jset))
           sgfb = first_sgfb(1,jset)

      IF (set_radius_a(iset) + set_radius_b(jset) >= dab) THEN

!            *** Calculate the primitive overlap integrals ***
              CALL diffop(la_max(iset),npgfa(iset),zeta(:,iset),&
                   rpgfa(:,iset),la_min(iset), lb_max(jset),npgfb(jset),&
                   zetb(:,jset),rpgfb(:,jset),lb_min(jset),rab,difab)

!            *** Contraction ***
              CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                   alpha,difab(1,1,1),ldab,sphi_b(1,sgfb),ldsb,&
                   0.0_dp,work(1,1),ldwork)
              IF(iatom<=jatom) THEN
                CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                     1.0_dp,sphi_a(1,sgfa),ldsa,&
                     work(1,1),ldwork,&
                     1.0_dp,op_dip(1)%block(sgfa,sgfb),&
                     SIZE(op_dip(1)%block,1))

              ELSE
                CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                     1.0_dp,work(1,1),ldwork,&
                     sphi_a(1,sgfa),ldsa,&
                     1.0_dp,op_dip(1)%block(sgfb,sgfa),&
                     SIZE(op_dip(1)%block,1))

              END IF

!             *** Contraction ***
              CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                   alpha,difab(1,1,2),ldab,sphi_b(1,sgfb),ldsb,&
                   0.0_dp,work(1,1),ldwork)
              IF(iatom<=jatom) THEN
                CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                     1.0_dp,sphi_a(1,sgfa),ldsa,&
                     work(1,1),ldwork,&
                     1.0_dp,op_dip(2)%block(sgfa,sgfb),&
                     SIZE(op_dip(2)%block,1))
              ELSE
                CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                     1.0_dp,work(1,1),ldwork,&
                     sphi_a(1,sgfa),ldsa,&
                     1.0_dp,op_dip(2)%block(sgfb,sgfa),&
                     SIZE(op_dip(2)%block,1))
              END IF

!            *** Contraction ***
              CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                   alpha,difab(1,1,3),ldab,sphi_b(1,sgfb),ldsb,&
                   0.0_dp,work(1,1),ldwork)
              IF(iatom<=jatom) THEN
                CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                     1.0_dp,sphi_a(1,sgfa),ldsa,&
                     work(1,1),ldwork,&
                     1.0_dp,op_dip(3)%block(sgfa,sgfb),&
                     SIZE(op_dip(3)%block,1))
              ELSE
                CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                     1.0_dp,work(1,1),ldwork,&
                     sphi_a(1,sgfa),ldsa,&
                     1.0_dp,op_dip(3)%block(sgfb,sgfa),&
                    SIZE(op_dip(3)%block,1))
             END IF
           END IF  !  >= dab

         END DO ! jset

       END DO ! iset

    END DO
    CALL neighbor_list_iterator_release(nl_iterator)

    DO i = 1,3
       NULLIFY(op_dip(i)%block)
    END DO
    DEALLOCATE(op_dip, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(difab,work,basis_set_list,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE p_xyz_ao

! *****************************************************************************
!> \brief Calculation of the components of the dipole operator in the length form
!>      by taking the relative position operator r-Rc, with respect a reference point Rc
!>      Probably it does not work for PBC, or maybe yes if the wfn are
!>      sufficiently localized
!>      The elements of the  sparse matrices are the integrals in the
!>      basis functions
!> \param op matrix representation of the p operator
!>               calculated in terms of the contracted basis functions
!> \param qs_env enviroment for the lists and the basis sets
!> \param rc reference vector position
!> \param order maximum order of the momentum, for the doipole order = 1
!> \param minimum_image take into account only the first neighbors in the lists
!> \param soft ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      03.2006 created [MI]
!> \author MI
! *****************************************************************************

  SUBROUTINE rRc_xyz_ao(op,qs_env,rc,order,minimum_image,soft,error)

    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: op
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(dp)                                 :: Rc(3)
    INTEGER, INTENT(IN)                      :: order
    LOGICAL, INTENT(IN), OPTIONAL            :: minimum_image, soft
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rRc_xyz_ao', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iatom, icol, ikind, imom, inode, irow, iset, istat, &
      jatom, jkind, jset, last_jatom, ldab, ldsa, ldsb, ldwork, M_dim, maxl, &
      ncoa, ncob, nkind, nrow, nseta, nsetb, sgfa, sgfb
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                npgfa, npgfb, nsgfa, nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure, found, &
                                                my_minimum_image, my_soft, &
                                                new_atom_b
    REAL(KIND=dp)                            :: dab, Lxo2, Lyo2, Lzo2, rab2
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rac, rb, rbc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                work, zeta, zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: mab
    TYPE(block_p_type), DIMENSION(:), &
      POINTER                                :: op_dip
    TYPE(cell_type), POINTER                 :: cell
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_kind_type), POINTER              :: qs_kind

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY (qs_kind, qs_kind_set)
    NULLIFY (cell, particle_set)
    NULLIFY (sab_orb)
    NULLIFY (mab, op_dip, work)
    NULLIFY (la_max, la_min, lb_max,  npgfa, npgfb, nsgfa, nsgfb)
    NULLIFY (set_radius_a,set_radius_b,rpgfa, rpgfb, sphi_a,sphi_b,zeta, zetb )

    CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set,&
         cell=cell,particle_set=particle_set,sab_orb=sab_orb,&
         error=error)

    nkind = SIZE(qs_kind_set)

    CALL get_qs_kind_set(qs_kind_set=qs_kind_set,&
                             maxco=ldwork, maxlgto=maxl )

    my_minimum_image = .FALSE.
    IF(PRESENT(minimum_image)) THEN
       my_minimum_image=minimum_image
       Lxo2 = SQRT ( SUM ( cell % hmat ( :, 1 ) ** 2 ) )/2.0_dp
       Lyo2 = SQRT ( SUM ( cell % hmat ( :, 2 ) ** 2 ) )/2.0_dp
       Lzo2 = SQRT ( SUM ( cell % hmat ( :, 3 ) ** 2 ) )/2.0_dp
    END IF
    my_soft = .FALSE.
    IF(PRESENT(soft)) THEN
      my_soft = soft
    END IF

    ldab = ldwork

    M_dim = ncoset(order)-1
    CPPostcondition(M_dim<=SIZE(op,1),cp_failure_level,routineP,error,failure)

    ALLOCATE(mab(ldab,ldab,M_dim),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    mab(1:ldab,1:ldab,1:M_dim) = 0.0_dp
    ALLOCATE(work(ldwork,ldwork),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    work(1:ldwork,1:ldwork) = 0.0_dp
    ALLOCATE (op_dip(M_dim),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO imom = 1, M_dim
      NULLIFY (op_dip(imom)%block)
    END DO

    ALLOCATE (basis_set_list(nkind),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ikind=1,nkind
      qs_kind => qs_kind_set(ikind)
      CALL get_qs_kind(qs_kind=qs_kind,softb=my_soft,orb_basis_set=basis_set_a)
      IF (ASSOCIATED(basis_set_a)) THEN
        basis_set_list(ikind)%gto_basis_set => basis_set_a
      ELSE
        NULLIFY(basis_set_list(ikind)%gto_basis_set)
      END IF
    END DO
    CALL neighbor_list_iterator_create(nl_iterator,sab_orb)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,inode=inode,&
                              iatom=iatom,jatom=jatom,r=rab)
       basis_set_a => basis_set_list(ikind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_a)) CYCLE
       basis_set_b => basis_set_list(jkind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_b)) CYCLE
       ra = pbc(particle_set(iatom)%r,cell)
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
       npgfb        =>  basis_set_b%npgf
       nsetb        =   basis_set_b%nset
       nsgfb        =>  basis_set_b%nsgf_set
       rpgfb        =>  basis_set_b%pgf_radius
       set_radius_b =>  basis_set_b%set_radius
       sphi_b       =>  basis_set_b%sphi
       zetb         =>  basis_set_b%zet

       ldsa = SIZE(sphi_a,1)
       ldsb = SIZE(sphi_b,1)
       IF(inode==1) last_jatom = 0

       IF(my_minimum_image) THEN
         IF(ABS(rab(1)) > Lxo2 .OR. ABS(rab(2)) > Lyo2 .OR. ABS(rab(3)) > Lzo2) CYCLE
       END IF

       rb = rab + ra

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

         DO imom = 1, M_dim
            NULLIFY(op_dip(imom)%block)
            CALL cp_dbcsr_get_block_p(matrix=op(imom)%matrix,&
                 row=irow,col=icol,block=op_dip(imom)%block,found=found)
            CPPostcondition(ASSOCIATED(op_dip(imom)%block),cp_failure_level,routineP,error,failure)
         END DO   ! imom
       END IF  ! new_atom_b

       rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
       dab = SQRT(rab2)

       nrow = 0
       DO iset=1,nseta

         ncoa = npgfa(iset)*ncoset(la_max(iset))
         sgfa = first_sgfa(1,iset)

         DO jset=1,nsetb

           ncob = npgfb(jset)*ncoset(lb_max(jset))
           sgfb = first_sgfb(1,jset)

           IF (set_radius_a(iset) + set_radius_b(jset) >= dab) THEN

             rac = pbc(rc,ra,cell)
             rbc = pbc(rc,rb,cell)

!            *** Calculate the primitive overlap integrals ***
             CALL moment(la_max(iset),npgfa(iset),zeta(:,iset),&
               rpgfa(:,iset),la_min(iset),&
               lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),&
               order,rac,rbc,mab)

             DO imom = 1,M_dim
!                 *** Contraction ***
               CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                   1.0_dp,mab(1,1,imom),ldab,sphi_b(1,sgfb),ldsb,&
                   0.0_dp,work(1,1),ldwork)
               IF(iatom<=jatom) THEN
                 CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                     1.0_dp,sphi_a(1,sgfa),ldsa,&
                     work(1,1),ldwork,&
                     1.0_dp,op_dip(imom)%block(sgfa,sgfb),&
                     SIZE(op_dip(imom)%block,1))
               ELSE
                 CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                     1.0_dp,work(1,1),ldwork,&
                     sphi_a(1,sgfa),ldsa,&
                     1.0_dp,op_dip(imom)%block(sgfb,sgfa),&
                     SIZE(op_dip(imom)%block,1))
               END IF

             END DO  ! imom
           END IF  !  >= dab

         END DO ! jset

         nrow = nrow + ncoa

       END DO ! iset

    END DO
    CALL neighbor_list_iterator_release(nl_iterator)

    DO i = 1,3
       NULLIFY(op_dip(i)%block)
    END DO
    DEALLOCATE(op_dip, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(mab,work,basis_set_list,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE   rRc_xyz_ao

! *****************************************************************************
!> \brief Calculation of the  multipole operators integrals
!>      and of its derivatives of the type
!>      [\mu | op | d(\nu)/dR(\nu)]-[d(\mu)/dR(\mu)| op | \nu]
!>      by taking the relative position operator r-Rc, with respect a reference point Rc
!>      The derivative are with respect to the primitive position,
!>      The multipole operator is symmetric and if it does not depend on R(\mu) or R(\nu)
!>      therefore  [\mu | op | d(\nu)/dR(\nu)] = -[d(\mu)/dR(\mu)| op | \nu]
!>        [\mu|op|d(\nu)/dR]-[d(\mu)/dR|op|\nu]=2[\mu|op|d(\nu)/dR]
!>      When it is not the case a correction term is needed
!>
!>     The momentum operator [\mu|M|\nu] is symmetric, the number of components is
!>     detrermined by the order: 3 for order 1 (x,y,x), 9 for order 2(xx,xy,xz,yy,yz,zz)
!>     The derivative of the type [\mu | op | d(\nu)/dR_i(\nu)], where
!>     i indicates the cartesian direction, is antisymmetric only when
!>     the no component M =(r_i) or (r_i r_j) is in the same cartesian
!>     direction of the derivative,  indeed
!>   d([\mu|M|\nu])/dr_i = [d(\mu)/dr_i|M|\nu] + [\mu|M|d(\nu)/dr_i] + [\mu |d(M)/dr_i|\nu]
!>   d([\mu|M|\nu])/dr_i = -[d(\mu)/dR_i(\mu)|M|\nu] -[\mu|M|d(\nu)/dR_i(\nu)] + [\mu |d(M)/dr_i|\nu]
!>     Therefore we cannot use an antisymmetric matrix
!>
!>     The same holds for the derivative with respect to the electronic position r
!>     taking into account that [\mu|op|d(\nu)/dR] = -[\mu|op|d(\nu)/dr]
!> \param op matrix representation of the p operator
!>               calculated in terms of the contracted basis functions
!> \param op_der ...
!> \param qs_env enviroment for the lists and the basis sets
!> \param rc reference vector position
!> \param order maximum order of the momentum, for the dipole order = 1
!> \param minimum_image take into account only the first neighbors in the lists
!> \param soft ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      03.2006 created [MI]
!> \author MI
!> \note
!>      Probably it does not work for PBC, or maybe yes if the wfn are
!>      sufficiently localized
!>      The elements of the  sparse matrices are the integrals in the
!>      basis functions
! *****************************************************************************
  SUBROUTINE rRc_xyz_der_ao(op,op_der,qs_env,rc,order,minimum_image,soft,error)

    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: op
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: op_der
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(dp)                                 :: Rc(3)
    INTEGER, INTENT(IN)                      :: order
    LOGICAL, INTENT(IN), OPTIONAL            :: minimum_image, soft
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rRc_xyz_der_ao', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iatom, icol, idir, ikind, imom, inode, ipgf, irow, &
      iset, istat, j, jatom, jkind, jpgf, jset, last_jatom, lda_min, ldab, &
      ldb_min, ldsa, ldsb, ldwork, M_dim, maxl, na, nb, ncoa, ncob, nda, ndb, &
      nkind, nseta, nsetb, sgfa, sgfb
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure, my_minimum_image, &
                                                my_soft, new_atom_b, &
                                                op_der_found, op_found
    REAL(KIND=dp)                            :: alpha, alpha_der, dab, Lxo2, &
                                                Lyo2, Lzo2, rab2
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rac, rb, rbc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                work, zeta, zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: mab, mab_tmp
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: difmab
    TYPE(block_p_type), DIMENSION(:), &
      POINTER                                :: op_dip
    TYPE(block_p_type), DIMENSION(:, :), &
      POINTER                                :: op_dip_der
    TYPE(cell_type), POINTER                 :: cell
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b, &
                                                orb_basis_set
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_all
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_kind_type), POINTER              :: qs_kind

    CALL timeset(routineN,handle)

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(op),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(op_der),cp_failure_level,routineP,error,failure)
    !IF(.NOT.op_sm_der(1,1)%matrix%symmetry=="none") THEN
    IF(.NOT.cp_dbcsr_get_matrix_type(op_der(1,1)%matrix).EQ.dbcsr_type_no_symmetry) THEN
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END IF

    NULLIFY (qs_kind, qs_kind_set)
    NULLIFY (cell, particle_set, orb_basis_set)
    NULLIFY (sab_all)
    NULLIFY (difmab, mab, mab_tmp)
    NULLIFY (op_dip, op_dip_der, work)
    NULLIFY (la_max, la_min, lb_max, lb_min, npgfa, npgfb, nsgfa, nsgfb)
    NULLIFY (set_radius_a,set_radius_b,rpgfa, rpgfb, sphi_a,sphi_b,zeta, zetb )

    CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set,&
         cell=cell,particle_set=particle_set,&
         sab_all=sab_all,error=error)

    nkind = SIZE(qs_kind_set)

    CALL get_qs_kind_set(qs_kind_set=qs_kind_set,&
                             maxco=ldwork, maxlgto=maxl )

    my_minimum_image = .FALSE.
    IF(PRESENT(minimum_image)) THEN
       my_minimum_image=minimum_image
       Lxo2 = SQRT ( SUM ( cell % hmat ( :, 1 ) ** 2 ) )/2.0_dp
       Lyo2 = SQRT ( SUM ( cell % hmat ( :, 2 ) ** 2 ) )/2.0_dp
       Lzo2 = SQRT ( SUM ( cell % hmat ( :, 3 ) ** 2 ) )/2.0_dp
    END IF
    my_soft = .FALSE.
    IF(PRESENT(soft)) THEN
      my_soft = soft
    END IF

    ldab = ldwork

    M_dim = ncoset(order)-1
    CPPostcondition(M_dim<=SIZE(op,1),cp_failure_level,routineP,error,failure)

    ALLOCATE(mab(ldab,ldab,M_dim),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    mab(1:ldab,1:ldab,1:M_dim) = 0.0_dp
    ALLOCATE(difmab(ldab,ldab,M_dim,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    difmab(1:ldab,1:ldab,1:M_dim,1:3) = 0.0_dp

    ALLOCATE(work(ldwork,ldwork),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    work(1:ldwork,1:ldwork) = 0.0_dp
    ALLOCATE (op_dip(M_dim),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (op_dip_der(M_dim,3),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DO imom = 1, M_dim
       NULLIFY (op_dip(imom)%block)
       DO i =1,3
          NULLIFY (op_dip_der(imom,i)%block)
       END DO
    END DO

    ALLOCATE (basis_set_list(nkind),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ikind=1,nkind
      qs_kind => qs_kind_set(ikind)
      CALL get_qs_kind(qs_kind=qs_kind,softb=my_soft,orb_basis_set=basis_set_a)
      IF (ASSOCIATED(basis_set_a)) THEN
        basis_set_list(ikind)%gto_basis_set => basis_set_a
      ELSE
        NULLIFY(basis_set_list(ikind)%gto_basis_set)
      END IF
    END DO
    CALL neighbor_list_iterator_create(nl_iterator,sab_all)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,inode=inode,&
                              iatom=iatom,jatom=jatom,r=rab)
       basis_set_a => basis_set_list(ikind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_a)) CYCLE
       basis_set_b => basis_set_list(jkind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_b)) CYCLE
       ra = pbc(particle_set(iatom)%r,cell)
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

       ldsa = SIZE(sphi_a,1)
       IF(ldsa.EQ.0) CYCLE
       ldsb = SIZE(sphi_b,1)
       IF(ldsb.EQ.0) CYCLE
       IF(inode==1) last_jatom = 0

       IF(my_minimum_image) THEN
         IF(ABS(rab(1)) > Lxo2 .OR. ABS(rab(2)) > Lyo2 .OR. ABS(rab(3)) > Lzo2) CYCLE
       END IF

       rb = rab + ra

       IF (jatom /= last_jatom) THEN
         new_atom_b = .TRUE.
         last_jatom = jatom
       ELSE
         new_atom_b = .FALSE.
       END IF

       IF (new_atom_b) THEN
           irow = iatom
           icol = jatom
           alpha_der = 2.0_dp

         DO imom = 1, M_dim
            NULLIFY(op_dip(imom)%block)
            CALL cp_dbcsr_get_block_p (matrix=op(imom)%matrix,&
                                    row=irow, col=icol,&
                                    block=op_dip(imom)%block,&
                                    found=op_found)
CPPostcondition(op_found.and.ASSOCIATED(op_dip(imom)%block),cp_failure_level,routineP,error,failure)
            DO idir = 1,3
              NULLIFY(op_dip_der(imom,idir)%block)
              CALL cp_dbcsr_get_block_p (matrix=op_der(imom,idir)%matrix,&
                                      row=irow, col=icol,&
                                      block=op_dip_der(imom,idir)%block,&
                                      found=op_der_found)
CPPostcondition(op_der_found.and.ASSOCIATED(op_dip_der(imom,idir)%block),cp_failure_level,routineP,error,failure)
           END DO  ! idir
         END DO   ! imom
       END IF  ! new_atom_b

       rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
       dab = SQRT(rab2)

       DO iset=1,nseta

         ncoa = npgfa(iset)*ncoset(la_max(iset))
         sgfa = first_sgfa(1,iset)

         DO jset=1,nsetb

           ncob = npgfb(jset)*ncoset(lb_max(jset))
           sgfb = first_sgfb(1,jset)

           IF (set_radius_a(iset) + set_radius_b(jset) >= dab) THEN

              rac = pbc(rc,ra,cell)
              rbc = rac + rab
!                  rac = pbc(rc,ra,cell)
!                  rbc = pbc(rc,rb,cell)

             ALLOCATE(mab_tmp(npgfa(iset)*ncoset(la_max(iset)+1),&
                      npgfb(jset)*ncoset(lb_max(jset)+1),ncoset(order)-1),STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

             lda_min = MAX ( 0, la_min(iset)-1)
             ldb_min = MAX ( 0, lb_min(jset)-1)
!            *** Calculate the primitive overlap integrals ***
             CALL moment(la_max(iset)+1,npgfa(iset),zeta(:,iset),&
                  rpgfa(:,iset),lda_min,&
                  lb_max(jset)+1,npgfb(jset),zetb(:,jset),rpgfb(:,jset),&
                  order,rac,rbc,mab_tmp)

!            *** Calculate the derivatives
             CALL diff_momop(la_max(iset),npgfa(iset),zeta(:,iset),&
                  rpgfa(:,iset),la_min(iset),lb_max(jset),npgfb(jset),&
                  zetb(:,jset),rpgfb(:,jset),lb_min(jset),order,rac,rbc,&
                  difmab,mab_ext=mab_tmp)

! Contract and copy in the sparse matrix
             mab =0.0_dp
             DO imom = 1,M_dim
               na = 0
               nda = 0
               DO ipgf=1,npgfa(iset)
                 nb =0
                 ndb = 0
                 DO jpgf = 1,npgfb(jset)
                   DO j = 1,ncoset(lb_max(jset))
                     DO i = 1,ncoset(la_max(iset))
                       mab(i+na,j+nb,imom) = mab_tmp(i+nda,j+ndb,imom)
                     END DO  ! i
                   END DO  ! j
                   nb = nb + ncoset(lb_max(jset))
                   ndb = ndb + ncoset(lb_max(jset)+1)
                 END DO  ! jpgf
                 na = na + ncoset(la_max(iset))
                 nda = nda + ncoset(la_max(iset)+1)
               END DO  ! ipgf

!                 *** Contraction ***
               CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                   1.0_dp,mab(1,1,imom),ldab,sphi_b(1,sgfb),ldsb,&
                   0.0_dp,work(1,1),ldwork)
               CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                     1.0_dp,sphi_a(1,sgfa),ldsa,&
                     work(1,1),ldwork,&
                     1.0_dp,op_dip(imom)%block(sgfa,sgfb),&
                     SIZE(op_dip(imom)%block,1))

               alpha = -1.0_dp!-alpha_der
               DO idir = 1,3
                 CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                     alpha,difmab(1,1,imom,idir),ldab,sphi_b(1,sgfb),ldsb,&
                     0.0_dp,work(1,1),ldwork)
                 CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                       1.0_dp,sphi_a(1,sgfa),ldsa,&
                       work(1,1),ldwork,&
                       1.0_dp,op_dip_der(imom,idir)%block(sgfa,sgfb),&
                       SIZE(op_dip_der(imom,idir)%block,1))

               END DO  ! idir

             END DO  ! imom

             DEALLOCATE(mab_tmp,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
           END IF  !  >= dab

         END DO ! jset

       END DO ! iset

    END DO
    CALL neighbor_list_iterator_release(nl_iterator)

    DO i = 1,3
       NULLIFY(op_dip(i)%block)
    END DO
    DEALLOCATE(op_dip,op_dip_der, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(mab,difmab,work,basis_set_list,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE   rRc_xyz_der_ao

END MODULE qs_operators_ao

