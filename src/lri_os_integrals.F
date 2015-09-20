!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculates integral matrices for LRIGPW method using the Obara-Saika
!>        (OS) integral scheme
!>        lri : local resolution of the identity
!> \par History
!>      created [06.2015]
!> \author Dorothea Golze  
! *****************************************************************************
MODULE lri_os_integrals
  USE ai_overlap,                      ONLY: overlap
  USE ai_overlap3,                     ONLY: overlap3
  USE ai_overlap_aabb,                 ONLY: overlap_aabb
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE kinds,                           ONLY: dp
  USE lri_debug_integrals,             ONLY: overlap_aabb_test,&
                                             overlap_ab_test,&
                                             overlap_abc_test
  USE orbital_pointers,                ONLY: ncoset
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *****************************************************************************

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lri_os_integrals'

  PUBLIC :: lri_int_ab, lri_int_aba, lri_int_abb, lri_int_aabb

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief calculate overlap integrals (a,b)
!> \param sab integral (a,b)
!> \param dsab derivative of sab with respect to A
!> \param ra ...
!> \param rb ...
!> \param rab ...
!> \param fba basis at center A
!> \param fbb basis at center B
!> \param calculate_forces ...
!> \param debug integrals are debugged by recursive routines if requested
!> \param dmax maximal deviation between integrals when debugging
! *****************************************************************************
  SUBROUTINE lri_int_ab(sab,dsab,ra,rb,rab,fba,fbb,calculate_forces,debug,&
                        dmax)

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: sab
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: dsab
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra, rb, rab
    TYPE(gto_basis_set_type), POINTER        :: fba, fbb
    LOGICAL, INTENT(IN)                      :: calculate_forces, debug
    REAL(KIND=dp), INTENT(INOUT)             :: dmax

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_int_ab', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iset, jset, lds, m1, m2, maxco, maxcoa, maxcob, &
      maxl, maxla, maxlb, ncoa, ncob, nseta, nsetb, sgfa, sgfb, stat
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dab
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: sint
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: devab, swork
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                zeta, zetb

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(la_max, la_min, lb_max, lb_min, npgfa, npgfb, nsgfa, nsgfb,&
            first_sgfa, first_sgfb, set_radius_a, set_radius_b, rpgfa, rpgfb,&
            sphi_a, sphi_b, zeta, zetb)

    ! basis ikind
    first_sgfa   =>  fba%first_sgf
    la_max       =>  fba%lmax
    la_min       =>  fba%lmin
    npgfa        =>  fba%npgf
    nseta        =   fba%nset
    nsgfa        =>  fba%nsgf_set
    rpgfa        =>  fba%pgf_radius
    set_radius_a =>  fba%set_radius
    sphi_a       =>  fba%sphi
    zeta         =>  fba%zet
    ! basis jkind
    first_sgfb   =>  fbb%first_sgf
    lb_max       =>  fbb%lmax
    lb_min       =>  fbb%lmin
    npgfb        =>  fbb%npgf
    nsetb        =   fbb%nset
    nsgfb        =>  fbb%nsgf_set
    rpgfb        =>  fbb%pgf_radius
    set_radius_b =>  fbb%set_radius
    sphi_b       =>  fbb%sphi
    zetb         =>  fbb%zet

    CALL get_gto_basis_set(fba,maxco=maxcoa,maxl=maxla)
    CALL get_gto_basis_set(fbb,maxco=maxcob,maxl=maxlb)
    maxco = MAX(maxcoa,maxcob)
    IF(calculate_forces) THEN
     maxl = MAX(maxla+1,maxlb)
    ELSE
     maxl = MAX(maxla,maxlb)
    ENDIF
    lds = ncoset(maxl)
    ALLOCATE(sint(maxco,maxco),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)
    IF(calculate_forces) THEN
     !derivatives will be stored in devab(:,:,2:4)
     ALLOCATE(swork(lds,lds,4),devab(maxco,maxco,4),STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,failure)
     devab = 0._dp
    ELSE
     ALLOCATE(swork(lds,lds,1),STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,failure)
    ENDIF
    sint  = 0._dp
    swork = 0._dp

    dab = SQRT( SUM(rab**2) )

    DO iset=1,nseta

       ncoa = npgfa(iset)*ncoset(la_max(iset))
       sgfa = first_sgfa(1,iset)

       DO jset=1,nsetb

          IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

          ncob = npgfb(jset)*ncoset(lb_max(jset))
          sgfb = first_sgfb(1,jset)
          m1=sgfa+nsgfa(iset)-1
          m2=sgfb+nsgfb(jset)-1
          sint=0._dp
          swork = 0._dp

          ! calculate integrals
          IF(calculate_forces) THEN
             devab = 0._dp 
             CALL overlap(la_max(iset),la_min(iset),npgfa(iset),rpgfa(:,iset),zeta(:,iset),&
                          lb_max(jset),lb_min(jset),npgfb(jset),rpgfb(:,jset),zetb(:,jset),&
                          rab,dab,sint,0,.FALSE.,swork,lds,sdab=devab)
             DO i=1,3
              !NOTE: devab(:,:,2:4) contains all derivatives for lmin=0 to lmax=lmax 
              !      correct after contraction (multiply with zero for elements l < lmin) 
              CALL ab_contract(dsab(sgfa:m1,sgfb:m2,i),devab(:,:,i+1),sphi_a(:,sgfa:),&
                   sphi_b(:,sgfb:),ncoa,ncob,nsgfa(iset),nsgfb(jset))
             ENDDO

          ELSE
             CALL overlap(la_max(iset),la_min(iset),npgfa(iset),rpgfa(:,iset),zeta(:,iset),&
                          lb_max(jset),lb_min(jset),npgfb(jset),rpgfb(:,jset),zetb(:,jset),&
                          rab,dab,sint,0,.FALSE.,swork,lds)
          ENDIF
          ! debug if requested
          IF(debug) THEN
             CALL overlap_ab_test(la_max(iset),la_min(iset),npgfa(iset),zeta(:,iset),&
                                  lb_max(jset),lb_min(jset),npgfb(jset),zetb(:,jset),&
                                  ra,rb,sint,dmax)
          ENDIF

          CALL ab_contract(sab(sgfa:m1,sgfb:m2),sint,sphi_a(:,sgfa:),sphi_b(:,sgfb:),&
               ncoa,ncob,nsgfa(iset),nsgfb(jset))
       END DO
    END DO

    IF(calculate_forces) THEN
     DEALLOCATE(devab,STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,failure)
    ENDIF

    DEALLOCATE(sint,swork,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)

    CALL timestop(handle)

  END SUBROUTINE lri_int_ab

! *****************************************************************************
!> \brief calculate integrals (a,b,fa)
!> \param abaint integral (a,b,fa)
!> \param dabdaint derivative of abaint with respect to A
!> \param ra ...
!> \param rb ...
!> \param rab ...
!> \param oba orbital basis at center A
!> \param obb orbital basis at center B
!> \param fba auxiliary basis set at center A
!> \param calculate_forces ...
!> \param debug integrals are debugged by recursive routines if requested
!> \param dmax maximal deviation between integrals when debugging
! *****************************************************************************
  SUBROUTINE lri_int_aba(abaint,dabdaint,ra,rb,rab,oba,obb,fba,&
                         calculate_forces,debug,dmax)

    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: abaint
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      OPTIONAL, POINTER                      :: dabdaint
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra, rb, rab
    TYPE(gto_basis_set_type), POINTER        :: oba, obb, fba
    LOGICAL, INTENT(IN)                      :: calculate_forces, debug
    REAL(KIND=dp), INTENT(INOUT)             :: dmax

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_int_aba', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iset, jset, kaset, &
                                                m1, m2, m3, ncoa, ncob, ncoc, &
                                                nseta, nsetb, nsetca, sgfa, &
                                                sgfb, sgfc, stat
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, lca_max, lca_min, &
                                                npgfa, npgfb, npgfca, nsgfa, &
                                                nsgfb, nsgfca
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb, &
                                                first_sgfca
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dab, dac, dbc
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: saba
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: sabda, sdaba, sdabda
    REAL(KIND=dp), DIMENSION(3)              :: rac, rbc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b, &
                                                set_radius_ca
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, rpgfca, sphi_a, &
                                                sphi_b, sphi_ca, zeta, zetb, &
                                                zetca

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(la_max, la_min, lb_max, lb_min, lca_max, lca_min, npgfa, npgfb,&
             npgfca, nsgfa, nsgfb, nsgfca) 
    NULLIFY(first_sgfa, first_sgfb, first_sgfca, set_radius_a, set_radius_b,&
            set_radius_ca, rpgfa, rpgfb, rpgfca, sphi_a, sphi_b, sphi_ca,&
            zeta, zetb, zetca) 

    ! basis ikind
    first_sgfa    =>  oba%first_sgf
    la_max        =>  oba%lmax
    la_min        =>  oba%lmin
    npgfa         =>  oba%npgf
    nseta         =   oba%nset
    nsgfa         =>  oba%nsgf_set
    rpgfa         =>  oba%pgf_radius
    set_radius_a  =>  oba%set_radius
    sphi_a        =>  oba%sphi
    zeta          =>  oba%zet
    ! basis jkind
    first_sgfb    =>  obb%first_sgf
    lb_max        =>  obb%lmax
    lb_min        =>  obb%lmin
    npgfb         =>  obb%npgf
    nsetb         =   obb%nset
    nsgfb         =>  obb%nsgf_set
    rpgfb         =>  obb%pgf_radius
    set_radius_b  =>  obb%set_radius
    sphi_b        =>  obb%sphi
    zetb          =>  obb%zet

    ! basis RI A
    first_sgfca   =>  fba%first_sgf
    lca_max       =>  fba%lmax
    lca_min       =>  fba%lmin
    npgfca        =>  fba%npgf
    nsetca        =   fba%nset
    nsgfca        =>  fba%nsgf_set
    rpgfca        =>  fba%pgf_radius
    set_radius_ca =>  fba%set_radius
    sphi_ca       =>  fba%sphi
    zetca         =>  fba%zet

    dab = SQRT( SUM(rab**2) )

    DO iset=1,nseta

       ncoa = npgfa(iset)*ncoset(la_max(iset))
       sgfa = first_sgfa(1,iset)

       DO jset=1,nsetb

          IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

          ncob = npgfb(jset)*ncoset(lb_max(jset))
          sgfb = first_sgfb(1,jset)
          m1=sgfa+nsgfa(iset)-1
          m2=sgfb+nsgfb(jset)-1

          ! calculate integrals abaint and derivative [d(a,b,a)/dA] dabdaint if requested
          rac = 0._dp
          dac = 0._dp
          rbc = -rab
          dbc = dab

          DO kaset=1,nsetca

             IF (set_radius_b(jset) + set_radius_ca(kaset) < dab) CYCLE

             ncoc = npgfca(kaset)*ncoset(lca_max(kaset))
             sgfc = first_sgfca(1,kaset)
             m3=sgfc+nsgfca(kaset)-1
             IF(ncoa*ncob*ncoc > 0) THEN
                ALLOCATE(saba(ncoa,ncob,ncoc),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,failure)
                saba(:,:,:)    = 0._dp
                ! integrals
                IF(calculate_forces) THEN
                   ALLOCATE(sdaba(ncoa,ncob,ncoc,3),STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,failure)
                   ALLOCATE(sabda(ncoa,ncob,ncoc,3),STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,failure)
                   ALLOCATE(sdabda(ncoa,ncob,ncoc,3),STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,failure)
                   sdaba(:,:,:,:) = 0._dp
                   sabda(:,:,:,:) = 0._dp
                   sdabda(:,:,:,:) = 0._dp
                   CALL overlap3(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                                 lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                                 lca_max(kaset),npgfca(kaset),zetca(:,kaset),rpgfca(:,kaset),lca_min(kaset),&
                                 rab,dab,rac,dac,rbc,dbc,saba,sdaba,sabda)
                   !d(a,b,a)/dA = (da/dA,b,a) + (a,b,da/dA)
                   sdabda(:,:,:,:) = sdaba + sabda

                   DO i=1,3
                    CALL abc_contract(dabdaint(sgfa:m1,sgfb:m2,sgfc:m3,i),sdabda(:,:,:,i),& 
                         sphi_a(:,sgfa:),sphi_b(:,sgfb:),sphi_ca(:,sgfc:),&
                         ncoa,ncob,ncoc,nsgfa(iset),nsgfb(jset),nsgfca(kaset))
                   ENDDO

                   DEALLOCATE(sdaba,sabda,sdabda,STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,failure)
                ELSE
                   CALL overlap3(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                                 lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                                 lca_max(kaset),npgfca(kaset),zetca(:,kaset),rpgfca(:,kaset),lca_min(kaset),&
                                 rab,dab,rac,dac,rbc,dbc,saba)
                ENDIF
                ! debug if requested
                IF(debug) THEN
                   CALL overlap_abc_test(la_max(iset),npgfa(iset),zeta(:,iset),la_min(iset),&
                                         lb_max(jset),npgfb(jset),zetb(:,jset),lb_min(jset),&
                                         lca_max(kaset),npgfca(kaset),zetca(:,kaset),lca_min(kaset),&
                                         ra,rb,ra,saba,dmax)
                ENDIF
                CALL abc_contract(abaint(sgfa:m1,sgfb:m2,sgfc:m3),saba,& 
                     sphi_a(:,sgfa:),sphi_b(:,sgfb:),sphi_ca(:,sgfc:),&
                     ncoa,ncob,ncoc,nsgfa(iset),nsgfb(jset),nsgfca(kaset))
                DEALLOCATE(saba,STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,failure)
             END IF
          END DO
       END DO
    END DO

    CALL timestop(handle)

  END SUBROUTINE lri_int_aba

! *****************************************************************************
!> \brief calculate integrals (a,b,fb)
!> \param abbint integral (a,b,fb)
!> \param dabbint derivative of abbint with respect to A
!> \param ra ...
!> \param rb ...
!> \param rab ...
!> \param oba orbital basis at center A
!> \param obb orbital basis at center B
!> \param fbb auxiliary basis set at center B
!> \param calculate_forces ...
!> \param debug integrals are debugged by recursive routines if requested
!> \param dmax maximal deviation between integrals when debugging
! *****************************************************************************
  SUBROUTINE lri_int_abb(abbint,dabbint,ra,rb,rab,oba,obb,fbb,calculate_forces,&
                         debug,dmax)

    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: abbint
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      OPTIONAL, POINTER                      :: dabbint
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra, rb, rab
    TYPE(gto_basis_set_type), POINTER        :: oba, obb, fbb
    LOGICAL, INTENT(IN)                      :: calculate_forces, debug
    REAL(KIND=dp), INTENT(INOUT)             :: dmax

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_int_abb', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iset, jset, kbset, &
                                                m1, m2, m3, ncoa, ncob, ncoc, &
                                                nseta, nsetb, nsetcb, sgfa, &
                                                sgfb, sgfc, stat
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, lcb_max, lcb_min, &
                                                npgfa, npgfb, npgfcb, nsgfa, &
                                                nsgfb, nsgfcb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb, &
                                                first_sgfcb
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dab, dac, dbc
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: sabb
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: sdabb
    REAL(KIND=dp), DIMENSION(3)              :: rac, rbc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b, &
                                                set_radius_cb
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, rpgfcb, sphi_a, &
                                                sphi_b, sphi_cb, zeta, zetb, &
                                                zetcb

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(la_max, la_min, lb_max, lb_min, lcb_max, lcb_min, npgfa, npgfb,&
            npgfcb, nsgfa, nsgfb, nsgfcb)
    NULLIFY(first_sgfa, first_sgfb, first_sgfcb, set_radius_a, set_radius_b,&
            set_radius_cb, rpgfa, rpgfb, rpgfcb, sphi_a, sphi_b, sphi_cb,&
            zeta, zetb, zetcb) 

    ! basis ikind
    first_sgfa    =>  oba%first_sgf
    la_max        =>  oba%lmax
    la_min        =>  oba%lmin
    npgfa         =>  oba%npgf
    nseta         =   oba%nset
    nsgfa         =>  oba%nsgf_set
    rpgfa         =>  oba%pgf_radius
    set_radius_a  =>  oba%set_radius
    sphi_a        =>  oba%sphi
    zeta          =>  oba%zet
    ! basis jkind
    first_sgfb    =>  obb%first_sgf
    lb_max        =>  obb%lmax
    lb_min        =>  obb%lmin
    npgfb         =>  obb%npgf
    nsetb         =   obb%nset
    nsgfb         =>  obb%nsgf_set
    rpgfb         =>  obb%pgf_radius
    set_radius_b  =>  obb%set_radius
    sphi_b        =>  obb%sphi
    zetb          =>  obb%zet

    ! basis RI B
    first_sgfcb   =>  fbb%first_sgf
    lcb_max       =>  fbb%lmax
    lcb_min       =>  fbb%lmin
    npgfcb        =>  fbb%npgf
    nsetcb        =   fbb%nset
    nsgfcb        =>  fbb%nsgf_set
    rpgfcb        =>  fbb%pgf_radius
    set_radius_cb =>  fbb%set_radius
    sphi_cb       =>  fbb%sphi
    zetcb         =>  fbb%zet

    dab = SQRT( SUM(rab**2) )

    DO iset=1,nseta

       ncoa = npgfa(iset)*ncoset(la_max(iset))
       sgfa = first_sgfa(1,iset)

       DO jset=1,nsetb

          IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

          ncob = npgfb(jset)*ncoset(lb_max(jset))
          sgfb = first_sgfb(1,jset)
          m1=sgfa+nsgfa(iset)-1
          m2=sgfb+nsgfb(jset)-1

          ! calculate integrals abbint and derivative [d(a,b,b)/dA] dabbint if requested
          rac = rab
          dac = dab
          rbc = 0._dp
          dbc = 0._dp

          DO kbset=1,nsetcb
 
             IF (set_radius_a(iset) + set_radius_cb(kbset) < dab) CYCLE

             ncoc = npgfcb(kbset)*ncoset(lcb_max(kbset))
             sgfc = first_sgfcb(1,kbset)
             m3=sgfc+nsgfcb(kbset)-1
             IF(ncoa*ncob*ncoc > 0) THEN
                ALLOCATE(sabb(ncoa,ncob,ncoc),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,failure)
                sabb(:,:,:) = 0._dp
                IF(calculate_forces) THEN
                   ALLOCATE(sdabb(ncoa,ncob,ncoc,3),STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,failure)
                   sdabb(:,:,:,:) = 0._dp
                   CALL overlap3(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                                lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                                lcb_max(kbset),npgfcb(kbset),zetcb(:,kbset),rpgfcb(:,kbset),lcb_min(kbset),&
                                rab,dab,rac,dac,rbc,dbc,sabb,sdabc=sdabb)
                   DO i=1,3
                    CALL abc_contract(dabbint(sgfa:m1,sgfb:m2,sgfc:m3,i),sdabb(:,:,:,i),&
                         sphi_a(:,sgfa:),sphi_b(:,sgfb:),sphi_cb(:,sgfc:),&
                         ncoa,ncob,ncoc,nsgfa(iset),nsgfb(jset),nsgfcb(kbset))
                   ENDDO
                   DEALLOCATE(sdabb,STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,failure)
                ELSE
                   CALL overlap3(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                                lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                                lcb_max(kbset),npgfcb(kbset),zetcb(:,kbset),rpgfcb(:,kbset),lcb_min(kbset),&
                                rab,dab,rac,dac,rbc,dbc,sabb)
                ENDIF
                ! debug if requested
                IF(debug) THEN
                   CALL overlap_abc_test(la_max(iset),npgfa(iset),zeta(:,iset),la_min(iset),&
                                         lb_max(jset),npgfb(jset),zetb(:,jset),lb_min(jset),&
                                         lcb_max(kbset),npgfcb(kbset),zetcb(:,kbset),lcb_min(kbset),&
                                         ra,rb,rb,sabb,dmax)
                ENDIF
                CALL abc_contract(abbint(sgfa:m1,sgfb:m2,sgfc:m3),sabb,&
                     sphi_a(:,sgfa:),sphi_b(:,sgfb:),sphi_cb(:,sgfc:),&
                     ncoa,ncob,ncoc,nsgfa(iset),nsgfb(jset),nsgfcb(kbset))
                DEALLOCATE(sabb,STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,failure)
             ENDIF
          END DO

       END DO
    END DO

    CALL timestop(handle)

  END SUBROUTINE lri_int_abb

! *****************************************************************************
!> \brief calculate overlap integrals (aa,bb)
!> \param saabb integral (aa,bb)
!> \param oba orbital basis at center A
!> \param obb orbital basis at center B
!> \param rab ...
!> \param ra ...
!> \param rb ...
!> \param debug integrals are debugged by recursive routines if requested
!> \param dmax maximal deviation between integrals when debugging
! *****************************************************************************
  SUBROUTINE lri_int_aabb(saabb,oba,obb,rab,ra,rb,debug,dmax)

    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: saabb
    TYPE(gto_basis_set_type), POINTER        :: oba, obb
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab, ra, rb
    LOGICAL, INTENT(IN)                      :: debug
    REAL(KIND=dp), INTENT(INOUT)             :: dmax

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_int_aabb', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iset, isgfa1, jset, jsgfa2, kset, ksgfb1, lds, lset, &
      lsgfb2, m1, m2, m3, m4, maxco, maxcoa, maxcob, maxl, maxla, maxlb, &
      ncoa1, ncoa2, ncob1, ncob2, nseta, nsetb, sgfa1, sgfa2, sgfb1, sgfb2, &
      stat
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: asets_equal, bsets_equal, &
                                                failure
    REAL(KIND=dp)                            :: dab
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: swork
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: sint
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                zeta, zetb

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(la_max, la_min, lb_max, lb_min, npgfa, npgfb, nsgfa, nsgfb,&
            first_sgfa, first_sgfb, set_radius_a, set_radius_b, rpgfa, rpgfb,&
            sphi_a, sphi_b, zeta, zetb)

    ! basis ikind
    first_sgfa   =>  oba%first_sgf
    la_max       =>  oba%lmax
    la_min       =>  oba%lmin
    npgfa        =>  oba%npgf
    nseta        =   oba%nset
    nsgfa        =>  oba%nsgf_set
    rpgfa        =>  oba%pgf_radius
    set_radius_a =>  oba%set_radius
    sphi_a       =>  oba%sphi
    zeta         =>  oba%zet
    ! basis jkind
    first_sgfb   =>  obb%first_sgf
    lb_max       =>  obb%lmax
    lb_min       =>  obb%lmin
    npgfb        =>  obb%npgf
    nsetb        =   obb%nset
    nsgfb        =>  obb%nsgf_set
    rpgfb        =>  obb%pgf_radius
    set_radius_b =>  obb%set_radius
    sphi_b       =>  obb%sphi
    zetb         =>  obb%zet

    CALL get_gto_basis_set(oba,maxco=maxcoa,maxl=maxla)
    CALL get_gto_basis_set(obb,maxco=maxcob,maxl=maxlb)
    maxco = MAX(maxcoa,maxcob)
    maxla = 2*maxla
    maxlb = 2*maxlb
    maxl = MAX(maxla,maxlb)
    lds = ncoset(maxl)
    ALLOCATE(sint(maxco,maxco,maxco,maxco),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)
    ALLOCATE(swork(lds,lds),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)
    sint  = 0._dp
    swork = 0._dp

    dab = SQRT( SUM(rab**2) )

    DO iset=1,nseta

       ncoa1 = npgfa(iset)*ncoset(la_max(iset))
       sgfa1 = first_sgfa(1,iset)
       m1    = sgfa1 + nsgfa(iset)-1

       DO jset=iset,nseta

          ncoa2 = npgfa(jset)*ncoset(la_max(jset))
          sgfa2 = first_sgfa(1,jset)
          m2    = sgfa2 + nsgfa(jset)-1

          DO kset=1,nsetb

             ncob1 = npgfb(kset)*ncoset(lb_max(kset))
             sgfb1 = first_sgfb(1,kset)
             m3    = sgfb1 + nsgfb(kset)-1

             DO lset=kset,nsetb

                ncob2 = npgfb(lset)*ncoset(lb_max(lset))
                sgfb2 = first_sgfb(1,lset)
                m4=sgfb2+nsgfb(lset)-1

                ! check if sets are identical to spare some integral evaluation
                asets_equal = .FALSE.
                IF(iset == jset) asets_equal = .TRUE.
                bsets_equal = .FALSE.
                IF(kset == lset) bsets_equal = .TRUE.
                ! calculate integrals
                CALL overlap_aabb(la_max(iset),la_min(iset),npgfa(iset),rpgfa(:,iset),zeta(:,iset),&
                                  la_max(jset),la_min(jset),npgfa(jset),rpgfa(:,jset),zeta(:,jset),&
                                  lb_max(kset),lb_min(kset),npgfb(kset),rpgfb(:,kset),zetb(:,kset),&
                                  lb_max(lset),lb_min(lset),npgfb(lset),rpgfb(:,lset),zetb(:,lset),&
                                  asets_equal,bsets_equal,rab,dab,sint,swork,lds)
                ! debug if requested
                IF (debug) THEN
                   CALL overlap_aabb_test(la_max(iset),la_min(iset),npgfa(iset),zeta(:,iset),&
                                          la_max(jset),la_min(jset),npgfa(jset),zeta(:,jset),&
                                          lb_max(kset),lb_min(kset),npgfb(kset),zetb(:,kset),&
                                          lb_max(lset),lb_min(lset),npgfb(lset),zetb(:,lset),&
                                          ra,rb,sint,dmax)
                ENDIF

                CALL abcd_contract(saabb(sgfa1:m1,sgfa2:m2,sgfb1:m3,sgfb2:m4),sint,sphi_a(:,sgfa1:),&
                                   sphi_a(:,sgfa2:),sphi_b(:,sgfb1:),sphi_b(:,sgfb2:),ncoa1,ncoa2,&
                                   ncob1,ncob2,nsgfa(iset),nsgfa(jset),nsgfb(kset),nsgfb(lset))

               ! account for the fact that some integrals are alike
                DO isgfa1 = sgfa1,m1
                  DO jsgfa2 = sgfa2,m2
                    DO ksgfb1 = sgfb1,m3
                      DO lsgfb2 = sgfb2,m4
                         saabb(jsgfa2,isgfa1,ksgfb1,lsgfb2) = saabb(isgfa1,jsgfa2,ksgfb1,lsgfb2)
                         saabb(isgfa1,jsgfa2,lsgfb2,ksgfb1) = saabb(isgfa1,jsgfa2,ksgfb1,lsgfb2)
                         saabb(jsgfa2,isgfa1,lsgfb2,ksgfb1) = saabb(isgfa1,jsgfa2,ksgfb1,lsgfb2)
                      END DO
                    END DO
                  END DO
                END DO
               
             END DO
          END DO
       END DO
    END DO

    DEALLOCATE(sint,swork,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)

    CALL timestop(handle)

  END SUBROUTINE lri_int_aabb

! *****************************************************************************
!> \brief contract overlap integrals (a,b) and transfer to spherical Gaussians
!> \param abint ...
!> \param sab ...
!> \param sphi_a ...
!> \param sphi_b ...
!> \param ncoa ...
!> \param ncob ...
!> \param nsgfa ...
!> \param nsgfb ...
! *****************************************************************************
  SUBROUTINE ab_contract(abint,sab,sphi_a,sphi_b,ncoa,ncob,nsgfa,nsgfb)

    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: abint
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: sab, sphi_a, sphi_b
    INTEGER, INTENT(IN)                      :: ncoa, ncob, nsgfa, nsgfb

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ab_contract', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: m1, m2, msphia, msphib, nn, &
                                                stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: cpp

    msphia = SIZE(sphi_a,1)
    msphib = SIZE(sphi_b,1)

    m1 = SIZE(sab,1)
    m2 = SIZE(sab,2)

    nn = SIZE(abint,1)

    ALLOCATE(cpp(nsgfa,m2),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)

    CALL dgemm("T","N",nsgfa,m2,ncoa,1._dp,sphi_a,msphia,sab,m1,0.0_dp,cpp,nsgfa)
    CALL dgemm("N","N",nsgfa,nsgfb,ncob,1._dp,cpp,nsgfa,sphi_b,msphib,0.0_dp,&
               abint,nn)

    DEALLOCATE(cpp,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)

  END SUBROUTINE ab_contract

! *****************************************************************************
!> \brief contract three-center overlap integrals (a,b,c) and transfer
!>        to spherical Gaussians
!> \param abcint ...
!> \param sabc ...
!> \param sphi_a ...
!> \param sphi_b ...
!> \param sphi_c ...
!> \param ncoa ...
!> \param ncob ...
!> \param ncoc ...
!> \param nsgfa ...
!> \param nsgfb ...
!> \param nsgfc ...
! *****************************************************************************
  SUBROUTINE abc_contract(abcint,sabc,sphi_a,sphi_b,sphi_c,ncoa,ncob,ncoc,&
                          nsgfa,nsgfb,nsgfc)

    REAL(KIND=dp), DIMENSION(:, :, :)        :: abcint, sabc
    REAL(KIND=dp), DIMENSION(:, :)           :: sphi_a, sphi_b, sphi_c
    INTEGER, INTENT(IN)                      :: ncoa, ncob, ncoc, nsgfa, &
                                                nsgfb, nsgfc

    CHARACTER(LEN=*), PARAMETER :: routineN = 'abc_contract', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, isgfc, m1, m2, m3, &
                                                msphia, msphib, msphic, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: temp_ccc, work_cpc
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: cpc, cpp

    CALL timeset(routineN,handle)

    msphia = SIZE(sphi_a,1)
    msphib = SIZE(sphi_b,1)
    msphic = SIZE(sphi_c,1)

    m1 = SIZE(sabc,1)
    m2 = SIZE(sabc,2)
    m3 = SIZE(sabc,3)

    ALLOCATE(cpp(nsgfa,m2,m3),cpc(nsgfa,m2,nsgfc),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)
    cpp = 0._dp
    cpc = 0._dp
    ALLOCATE(work_cpc(nsgfa,m2),temp_ccc(nsgfa,nsgfb),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)
    work_cpc(:,:) = 0._dp
    temp_ccc(:,:) = 0._dp

    CALL dgemm("T","N",nsgfa,m2*m3,ncoa,1._dp,sphi_a,msphia,sabc,m1,0.0_dp,cpp,nsgfa)
    CALL dgemm("N","N",nsgfa*m2,nsgfc,ncoc,1._dp,cpp,nsgfa*m2,sphi_c,msphic,0.0_dp,&
               cpc,nsgfa*m2)

    DO isgfc=1,nsgfc
     work_cpc(:,:)=cpc(:,:,isgfc)
     CALL dgemm("N","N",nsgfa,nsgfb,ncob,1._dp,work_cpc,nsgfa,sphi_b,msphib,&
                0.0_dp,temp_ccc,nsgfa)
     abcint(:,:,isgfc)=temp_ccc(:,:)
    END DO

    DEALLOCATE(cpp,cpc,work_cpc,temp_ccc,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)

    CALL timestop(handle)

  END SUBROUTINE abc_contract

! *****************************************************************************
!> \brief contract four-center overlap integrals (a,b,c,d) and transfer
!>        to spherical Gaussians
!> \param abcdint ...
!> \param sabcd ...
!> \param sphi_a ...
!> \param sphi_b ...
!> \param sphi_c ...
!> \param sphi_d ...
!> \param ncoa ...
!> \param ncob ...
!> \param ncoc ...
!> \param ncod ...
!> \param nsgfa ...
!> \param nsgfb ...
!> \param nsgfc ...
!> \param nsgfd ...
! *****************************************************************************
  SUBROUTINE abcd_contract(abcdint,sabcd,sphi_a,sphi_b,sphi_c,sphi_d,ncoa,ncob,&
                           ncoc,ncod,nsgfa,nsgfb,nsgfc,nsgfd)

    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      INTENT(INOUT)                          :: abcdint
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      INTENT(IN)                             :: sabcd
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: sphi_a, sphi_b, sphi_c, sphi_d
    INTEGER, INTENT(IN)                      :: ncoa, ncob, ncoc, ncod, &
                                                nsgfa, nsgfb, nsgfc, nsgfd

    CHARACTER(LEN=*), PARAMETER :: routineN = 'abcd_contract', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, isgfc, isgfd, m1, m2, &
                                                m3, m4, msphia, msphib, &
                                                msphic, msphid, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: temp_cccc, work_cpcc
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: temp_cpcc, work_cppc
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: cpcc, cppc, cppp

    CALL timeset(routineN,handle)

    msphia = SIZE(sphi_a,1)
    msphib = SIZE(sphi_b,1)
    msphic = SIZE(sphi_c,1)
    msphid = SIZE(sphi_d,1)

    m1 = SIZE(sabcd,1)
    m2 = SIZE(sabcd,2)
    m3 = SIZE(sabcd,3)
    m4 = SIZE(sabcd,4)

    ALLOCATE(cppp(nsgfa,m2,m3,m4),cppc(nsgfa,m2,m3,nsgfd),&
             cpcc(nsgfa,m2,nsgfc,nsgfd),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)

    ALLOCATE(work_cppc(nsgfa,m2,m3),temp_cpcc(nsgfa,m2,nsgfc),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)
    work_cppc = 0._dp
    temp_cpcc = 0._dp

    ALLOCATE(work_cpcc(nsgfa,m2),temp_cccc(nsgfa,nsgfb),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)
    work_cpcc = 0._dp
    temp_cccc = 0._dp

    CALL dgemm("T","N",nsgfa,m2*m3*m4,ncoa,1._dp,sphi_a,msphia,sabcd,m1,&
               0.0_dp,cppp,nsgfa)
    CALL dgemm("N","N",nsgfa*m2*m3,nsgfd,ncod,1._dp,cppp,nsgfa*m2*m3,&
               sphi_d,msphid,0.0_dp,cppc,nsgfa*m2*m3)

    DO isgfd=1,nsgfd
      work_cppc(:,:,:) = cppc(:,:,:,isgfd)
      CALL dgemm("N","N",nsgfa*m2,nsgfc,ncoc,1._dp,work_cppc,nsgfa*m2,&
                 sphi_c,msphic,0.0_dp,temp_cpcc,nsgfa*m2)
      cpcc(:,:,:,isgfd) = temp_cpcc(:,:,:)
      DO isgfc=1,nsgfc
        work_cpcc(:,:) = cpcc(:,:,isgfc,isgfd)
        CALL dgemm("N","N",nsgfa,nsgfb,ncob,1._dp,work_cpcc,nsgfa,sphi_b,&
                   msphib,0.0_dp,temp_cccc,nsgfa)
        abcdint(:,:,isgfc,isgfd) = temp_cccc(:,:) 
      END DO
    END DO

    DEALLOCATE(cpcc,cppc,cppp,stat=STAT)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)
    DEALLOCATE(work_cpcc,work_cppc,temp_cpcc,temp_cccc,stat=STAT)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)

    CALL timestop(handle)

  END SUBROUTINE abcd_contract

END MODULE lri_os_integrals
