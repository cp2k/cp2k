!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Several screening methods used in HFX calcualtions 
!> \par History
!>      04.2008 created [Manuel Guidon]
!> \author Manuel Guidon 
! *****************************************************************************
MODULE hfx_screening_methods 
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE f77_blas
  USE hfx_libint_interface,            ONLY: evaluate_eri_screen
  USE hfx_libint_wrapper,              ONLY: lib_int
  USE hfx_types,                       ONLY: hfx_basis_type,&
                                             hfx_potential_type,&
                                             hfx_screen_coeff_type,&
                                             log_zero,&
                                             powell_min_log
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE machine,                         ONLY: default_output_unit
  USE mathconstants
  USE orbital_pointers
  USE powell,                          ONLY: opt_state_type,&
                                             powell_optimize
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_interactions,                 ONLY: exp_radius_very_extended
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
   
  PUBLIC update_pmax_mat, &
         calc_screening_functions,&
         calc_pair_dist_radii
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_screening_methods'

!***
  
  CONTAINS

! *****************************************************************************
!> \brief calculates max values of two-electron integrals in a quartet/shell
!>      w.r.t. different zetas using the library lib_int
!> \param r position
!> \param z zeta
!> \param l angular momentum
!> \param npgf number of primitive cartesian gaussian in actual shell
!> \param max_val schwarz screening value
!> \param potential_parameter contains info for libint
!> \param i_thread thread id of the current task
!> \param tmp_R_1 pgf_based screening coefficients
!> \param rab2 Squared Distance of centers ab
!> \par History
!>     03.2007 created [Manuel Guidon]
!>     02.2009 refactored [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE screen4(lib, ra, rb, zeta, zetb, &
                     la_min, la_max, lb_min, lb_max,&
                     npgfa, npgfb, &
                     max_val, potential_parameter, i_thread, tmp_R_1,&
                     rab2)

    TYPE(lib_int)                            :: lib
    REAL(dp), INTENT(IN)                     :: ra(3), rb(3)
    REAL(dp), DIMENSION(:), INTENT(IN)       :: zeta, zetb
    INTEGER, INTENT(IN)                      :: la_min, la_max, lb_min, &
                                                lb_max, npgfa, npgfb
    REAL(dp), INTENT(INOUT)                  :: max_val
    TYPE(hfx_potential_type)                 :: potential_parameter
    INTEGER, INTENT(IN)                      :: i_thread
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :), POINTER               :: tmp_R_1
    REAL(dp)                                 :: rab2

    CHARACTER(len=*), PARAMETER :: routineN = 'screen4', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ipgf, jpgf, la, lb
    REAL(dp)                                 :: max_val_temp, R1

    max_val_temp = max_val
    DO ipgf = 1,npgfa
      DO jpgf = 1,npgfb
        R1 = MAX(0.0_dp, tmp_R_1(jpgf,ipgf)%x(1)*rab2 + tmp_R_1(jpgf,ipgf)%x(2))
        DO la = la_min,la_max
          DO lb = lb_min, lb_max
            !Build primitives
            CALL evaluate_eri_screen(lib, ra, rb, ra, rb,&
                                     zeta(ipgf), zetb(jpgf), zeta(ipgf), zetb(jpgf),&
                                     la, lb, la ,lb,&
                                     max_val_temp, potential_parameter, R1, R1)

            max_val = MAX(max_val, max_val_temp)
          END DO !lb
        END DO !la
      END DO !jpgf
    END DO !ipgf
  END SUBROUTINE screen4


! *****************************************************************************
!> \brief updates the maximum of the density matrix in compressed form for screening purposes 
!> \param pmax_set buffer to store matrix
!> \note
!>      - updates for each pair of shells the maximum absolute value of p
!> \par History
!>     09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE update_pmax_mat(pmax_set, full_density, full_density_beta, natom, &
                             max_set ,kind_of,basis_parameter, last_sgf_global,&
                             N_pmat, error)

    REAL(dp), DIMENSION(:, :, :, :), POINTER :: pmax_set
    REAL(dp), DIMENSION(:), POINTER          :: full_density, &
                                                full_density_beta
    INTEGER, INTENT(IN)                      :: natom, max_set
    INTEGER                                  :: kind_of(*)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    INTEGER                                  :: last_sgf_global(0:natom)
    INTEGER(int_8), INTENT(IN)               :: N_pmat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_pmax_mat', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, jatom, jkind, jset, katom, kkind, kset, mb, mc, nsetb, &
      nsetc, nsgfb_total, nsgfc_total, pb, pb1, pb2, pc, pc1, pc2, sgfb, sgfc
    INTEGER, DIMENSION(:), POINTER           :: nsgfb, nsgfc
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfb, first_sgfc
    LOGICAL                                  :: failure
    REAL(dp)                                 :: pmax_tmp

    CALL timeset(routineN,handle)

    failure = .FALSE.

    pmax_set = 0.0_dp

    DO jatom=1,natom
      jkind = kind_of(jatom)
      nsgfb_total = basis_parameter(jkind)%nsgf_total
      nsetb = basis_parameter(jkind)%nset
      first_sgfb => basis_parameter(jkind)%first_sgf
      nsgfb => basis_parameter(jkind)%nsgf
      pb1 = last_sgf_global(jatom) - nsgfb_total

      DO katom = 1,natom
        kkind = kind_of(katom)
        nsgfc_total = basis_parameter(kkind)%nsgf_total
        nsetc = basis_parameter(kkind)%nset
        first_sgfc => basis_parameter(kkind)%first_sgf
        nsgfc => basis_parameter(kkind)%nsgf
        pc1 = last_sgf_global(katom) - nsgfc_total
        DO jset = 1,nsetb
          sgfb = first_sgfb(1,jset)
          DO kset = 1,nsetc
            sgfc = first_sgfc(1,kset)
            pmax_tmp = 0.0_dp
            DO mc = 1,nsgfc(kset)
               pc2 = mc+sgfc-1
               pc = pc1 + pc2
               DO mb = 1,nsgfb(jset)
                  pb2 = mb+sgfb-1
                  pb = pb1 + pb2
                  pmax_tmp=MAX(pmax_tmp,ABS(full_density(get_1d_idx(pb,pc,N_pmat))))
                  IF( ASSOCIATED(full_density_beta) ) THEN
                    pmax_tmp=MAX(pmax_tmp,ABS(full_density_beta(get_1d_idx(pb,pc,N_pmat))))
                  END IF
               ENDDO
            ENDDO
            IF(pmax_tmp == 0.0_dp) THEN
              pmax_tmp = log_zero
            ELSE
              pmax_tmp = LOG10(pmax_tmp)
            END IF
            pmax_set(jset,kset,jatom,katom) = pmax_tmp
          END DO
        END DO
      END DO
    END DO

    CALL timestop(handle)

  END SUBROUTINE update_pmax_mat

! *****************************************************************************
!> \brief calculates screening functions for schwarz screening
!> \param qs_env qs_env
!> \param lib structure to libint
!> \param potential_parameter contains infos on potential
!> \param i_thread Thread ID of current task
!> \param coeffs_set set based coefficients
!> \param coeffs_kind kind based coefficients
!> \param coeffs_pgf pgf based coefficients
!> \param radii_pgf coefficients for long-range screening
!> \param max_set Maximum Number of basis set sets in the system
!> \param max_pgf Maximum Number of basis set pgfs in the system
!> \note
!>      This routine calculates (ab|ab) for different distances Rab = |a-b|
!>      and uses the powell optimiztion routines in order to fit the results
!>      in the following form:
!>
!>                 (ab|ab) = (ab|ab)(Rab) = c2*Rab^2 + c0
!>
!>      The missing linear term assures that the functions is monotonically
!>      decaying such that c2 can be used as upper bound when applying the
!>      Minimum Image Convention in the periodic case. Furthermore
!>      it seems to be a good choice to fit the logarithm of the (ab|ab)
!>      The fitting takes place at several levels: kind, set and pgf within
!>      the corresponding ranges of the prodiuct charge distributions
!>      Doing so, we only need arrays of size nkinds^2*2 instead of big
!>      screening matrices
!> \par History
!>     02.2009 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************

  SUBROUTINE calc_screening_functions(qs_env, lib, potential_parameter, i_thread, &
                                      coeffs_set, coeffs_kind, coeffs_pgf, radii_pgf,&
                                      max_set, max_pgf, error) 
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(lib_int)                            :: lib
    TYPE(hfx_potential_type)                 :: potential_parameter
    INTEGER, INTENT(IN)                      :: i_thread
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :, :, :), POINTER         :: coeffs_set
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :), POINTER               :: coeffs_kind
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :, :, :, :, :), POINTER   :: coeffs_pgf, radii_pgf
    INTEGER, INTENT(IN)                      :: max_set, max_pgf
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_screening_functions', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ikind, ipgf, iset, &
                                                jkind, jpgf, jset, la, lb, &
                                                ncoa, ncob, nkind, nseta, &
                                                nsetb, sgfa, sgfb, stat
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure
    REAL(dp) :: DATA(2,0:100), kind_radius_a, kind_radius_b, &
      max_contraction_a, max_contraction_b, max_val, max_val_temp, R1, ra(3), &
      radius, rb(3), x(2)
    REAL(dp), DIMENSION(:), POINTER          :: set_radius_a, set_radius_b
    REAL(dp), DIMENSION(:, :), POINTER       :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                zeta, zetb
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :), POINTER               :: tmp_R_1

    failure = .FALSE.

    CALL timeset(routineN,handle)


    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    error=error)

    nkind = SIZE(atomic_kind_set,1)

    ALLOCATE(coeffs_pgf(max_pgf,max_pgf,max_set,max_set,nkind,nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO ikind = 1,nkind
      DO jkind = 1,nkind
        DO iset=1,max_set
          DO jset=1,max_set
            DO ipgf=1,max_pgf
              DO jpgf=1,max_pgf
                coeffs_pgf(jpgf,ipgf,jset,iset,jkind,ikind)%x(:) = 0.0_dp
              END DO
            END DO   
          END DO
        END DO
      END DO
    END DO
    ra = 0.0_dp
    rb = 0.0_dp
    DATA = 0.0_dp
    DO ikind=1,nkind
      NULLIFY(atom_kind,orb_basis)
      atom_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atom_kind,orb_basis_set=orb_basis)
      NULLIFY(la_max,la_min,npgfa,zeta)
      CALL get_gto_basis_set(gto_basis_set=orb_basis,&
                             lmax=la_max,&
                             lmin=la_min,&
                             npgf=npgfa,&
                             nset=nseta,&
                             zet=zeta,&
                             set_radius=set_radius_a,&
                             first_sgf=first_sgfa,&
                             sphi=sphi_a,&
                             nsgf_set=nsgfa,&
                             pgf_radius=rpgfa)
      DO jkind = 1,nkind
        NULLIFY(atom_kind,orb_basis)
        atom_kind => atomic_kind_set(jkind)
        CALL get_atomic_kind(atomic_kind=atom_kind,orb_basis_set=orb_basis)
        NULLIFY(lb_max,lb_min,npgfb,zetb)
        CALL get_gto_basis_set(gto_basis_set=orb_basis,&
                               lmax=lb_max,&
                               lmin=lb_min,&
                               npgf=npgfb,&
                               nset=nsetb,&
                               zet=zetb,&
                               set_radius=set_radius_b,&
                               first_sgf=first_sgfb,&
                               sphi=sphi_b,&
                               nsgf_set=nsgfb,&
                               pgf_radius=rpgfb)
        DO iset = 1,nseta
          ncoa = npgfa(iset)*ncoset(la_max(iset))
          sgfa = first_sgfa(1,iset)
          max_contraction_a = MAXVAL((/(SUM(ABS(sphi_a(1:ncoa,i))),i=sgfa,sgfa+nsgfa(iset)-1)/))
          DO jset = 1,nsetb
            ncob = npgfb(jset)*ncoset(lb_max(jset))
            sgfb = first_sgfb(1,jset)
            max_contraction_b = MAXVAL((/(SUM(ABS(sphi_b(1:ncob,i))),i=sgfb,sgfb+nsgfb(jset)-1)/))
            radius=set_radius_a(iset)+set_radius_b(jset)
            DO ipgf = 1,npgfa(iset)
              DO jpgf = 1,npgfb(jset)
                DATA = 0.0_dp
                radius=rpgfa(ipgf,iset)+rpgfb(jpgf,jset)
                DO i=0,100 
                  rb(1) = 0.0_dp + REAL(i,dp) * 0.01_dp * radius
                  max_val = 0.0_dp
                  R1 = MAX(0.0_dp, radii_pgf(jpgf,ipgf,jset,iset,jkind,ikind)%x(1)*rb(1)**2 + &
                       radii_pgf(jpgf,ipgf,jset,iset,jkind,ikind)%x(2))
                  DO la = la_min(iset),la_max(iset)
                    DO lb = lb_min(jset), lb_max(jset)
                      !Build primitives
                      max_val_temp = 0.0_dp
                      CALL evaluate_eri_screen(lib, ra, rb, ra, rb,&
                                               zeta(ipgf,iset), zetb(jpgf,jset), zeta(ipgf,iset), zetb(jpgf,jset),&
                                               la, lb, la ,lb,&
                                               max_val_temp, potential_parameter, R1, R1)
                      max_val = MAX(max_val, max_val_temp)
                    END DO !lb
                  END DO !la
                  max_val = SQRT(max_val)
                  max_val = max_val * max_contraction_a * max_contraction_b
                  DATA(1,i) = rb(1)
                  IF(max_val == 0.0_dp ) THEN
                    DATA(2,i) = powell_min_log
                  ELSE
                    DATA(2,i) = LOG10((max_val))
                  END IF
                END DO

                CALL optimize_it(DATA,x,powell_min_log)
                coeffs_pgf(jpgf,ipgf,jset,iset,jkind,ikind)%x = x 

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
  
    ALLOCATE(coeffs_set(max_set,max_set,nkind,nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO ikind = 1,nkind
      DO jkind = 1,nkind
        DO iset=1,max_set
          DO jset=1,max_set
            coeffs_set(jset,iset,jkind,ikind)%x(:) = 0.0_dp
          END DO
        END DO
      END DO
    END DO
    ra = 0.0_dp
    rb = 0.0_dp
    DATA = 0.0_dp
    DO ikind=1,nkind
      NULLIFY(atom_kind,orb_basis)
      atom_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atom_kind,orb_basis_set=orb_basis)
      NULLIFY(la_max,la_min,npgfa,zeta)
      CALL get_gto_basis_set(gto_basis_set=orb_basis,&
                             lmax=la_max,&
                             lmin=la_min,&
                             npgf=npgfa,&
                             nset=nseta,&
                             zet=zeta,&
                             set_radius=set_radius_a,&
                             first_sgf=first_sgfa,&
                             sphi=sphi_a,&
                             nsgf_set=nsgfa)
      DO jkind = 1,nkind
        NULLIFY(atom_kind,orb_basis)
        atom_kind => atomic_kind_set(jkind)
        CALL get_atomic_kind(atomic_kind=atom_kind,orb_basis_set=orb_basis)
        NULLIFY(lb_max,lb_min,npgfb,zetb)
        CALL get_gto_basis_set(gto_basis_set=orb_basis,&
                               lmax=lb_max,&
                               lmin=lb_min,&
                               npgf=npgfb,&
                               nset=nsetb,&
                               zet=zetb,&
                               set_radius=set_radius_b,&
                               first_sgf=first_sgfb,&
                               sphi=sphi_b,&
                               nsgf_set=nsgfb)
        DO iset = 1,nseta
          ncoa = npgfa(iset)*ncoset(la_max(iset))
          sgfa = first_sgfa(1,iset)
          max_contraction_a = MAXVAL((/(SUM(ABS(sphi_a(1:ncoa,i))),i=sgfa,sgfa+nsgfa(iset)-1)/))
          DO jset = 1,nsetb
            ncob = npgfb(jset)*ncoset(lb_max(jset))
            sgfb = first_sgfb(1,jset)
            max_contraction_b = MAXVAL((/(SUM(ABS(sphi_b(1:ncob,i))),i=sgfb,sgfb+nsgfb(jset)-1)/))
            radius=set_radius_a(iset)+set_radius_b(jset)
            DATA = 0.0_dp
            tmp_R_1 => radii_pgf(:,:,jset,iset,jkind,ikind)
            DO i=0,100 
              rb(1) = 0.0_dp + REAL(i,dp) * 0.01_dp * radius
              max_val = 0.0_dp
              CALL screen4(lib, ra, rb, &
                           zeta(:,iset), zetb(:,jset), &
                           la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                           npgfa(iset), npgfb(jset), &
                           max_val, potential_parameter, i_thread, tmp_R_1, rb(1)**2)
              max_val = SQRT(max_val)
              max_val = max_val * max_contraction_a * max_contraction_b
              DATA(1,i) = rb(1)
              IF(max_val == 0.0_dp ) THEN
                DATA(2,i) = powell_min_log
              ELSE
                DATA(2,i) = LOG10((max_val))
              END IF
            END DO

            CALL optimize_it(DATA,x,powell_min_log)
            coeffs_set(jset,iset,jkind,ikind)%x = x 

          END DO
        END DO
       
      END DO
    END DO

    ! ** now kinds
    ALLOCATE(coeffs_kind(nkind,nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO ikind = 1,nkind
      DO jkind = 1,nkind
        coeffs_kind(jkind,ikind)%x(:) = 0.0_dp
      END DO
    END DO
    ra = 0.0_dp
    rb = 0.0_dp
    DO ikind=1,nkind
      NULLIFY(atom_kind,orb_basis)
      atom_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atom_kind,orb_basis_set=orb_basis)
      NULLIFY(la_max,la_min,npgfa,zeta)
      CALL get_gto_basis_set(gto_basis_set=orb_basis,&
                             lmax=la_max,&
                             lmin=la_min,&
                             npgf=npgfa,&
                             nset=nseta,&
                             zet=zeta,&
                             kind_radius=kind_radius_a,&
                             first_sgf=first_sgfa,&
                             sphi=sphi_a,&
                             nsgf_set=nsgfa)
      DO jkind = 1,nkind
        NULLIFY(atom_kind,orb_basis)
        atom_kind => atomic_kind_set(jkind)
        CALL get_atomic_kind(atomic_kind=atom_kind,orb_basis_set=orb_basis)
        NULLIFY(lb_max,lb_min,npgfb,zetb)
        CALL get_gto_basis_set(gto_basis_set=orb_basis,&
                               lmax=lb_max,&
                               lmin=lb_min,&
                               npgf=npgfb,&
                               nset=nsetb,&
                               zet=zetb,&
                               kind_radius=kind_radius_b,&
                               first_sgf=first_sgfb,&
                               sphi=sphi_b,&
                               nsgf_set=nsgfb)

        radius=kind_radius_a + kind_radius_b
        DATA = -100.0_dp
        DO iset=1,nseta
          ncoa = npgfa(iset)*ncoset(la_max(iset))
          sgfa = first_sgfa(1,iset)
          max_contraction_a = MAXVAL((/(SUM(ABS(sphi_a(1:ncoa,i))),i=sgfa,sgfa+nsgfa(iset)-1)/))
          DO jset=1,nsetb
            ncob = npgfb(jset)*ncoset(lb_max(jset))
            sgfb = first_sgfb(1,jset)
            max_contraction_b = MAXVAL((/(SUM(ABS(sphi_b(1:ncob,i))),i=sgfb,sgfb+nsgfb(jset)-1)/))
            DO i=0,100 
              rb(1) = 0.0_dp + REAL(i,dp) * 0.01_dp * radius
              max_val = 0.0_dp
              tmp_R_1 => radii_pgf(:,:,jset,iset,jkind,ikind)
              CALL screen4(lib, ra, rb, &
                           zeta(:,iset), zetb(:,jset), &
                           la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                           npgfa(iset), npgfb(jset), &
                           max_val, potential_parameter, i_thread, tmp_R_1, rb(1)**2)
              DATA(1,i) = rb(1)
              max_val = SQRT(max_val)
              max_val = max_val * max_contraction_a*max_contraction_b
              IF(max_val == 0.0_dp ) THEN
                DATA(2,i) = MAX(powell_min_log,DATA(2,i))
              ELSE
                DATA(2,i) = MAX(LOG10(max_val),DATA(2,i))
              END IF
            END DO
          END DO
        END DO

        CALL optimize_it(DATA,x,powell_min_log)
        coeffs_kind(jkind,ikind)%x = x 

      END DO
    END DO

    CALL timestop(handle)

  END SUBROUTINE calc_screening_functions

! *****************************************************************************
!> \brief calculates radius functions for longrange screening
!> \param qs_env qs_env
!> \param lib structure to libint
!> \param potential_parameter contains infos on potential
!> \param i_thread Thread ID of current task
!> \param radii_pgf pgf based coefficients
!> \param max_set Maximum Number of basis set sets in the system
!> \param max_pgf Maximum Number of basis set pgfs in the system
!> \note
!>      This routine calculates the pair-distribution radius of a product
!>      Gaussian and uses the powell optimiztion routines in order to fit 
!>      the results in the following form:
!>
!>                 (ab| = (ab(Rab) = c2*Rab^2 + c0
!>
!> \par History
!>     02.2009 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************

  SUBROUTINE calc_pair_dist_radii(qs_env, lib, potential_parameter, i_thread, &
                  radii_pgf, max_set, max_pgf, eps_schwarz, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(lib_int)                            :: lib
    TYPE(hfx_potential_type)                 :: potential_parameter
    INTEGER, INTENT(IN)                      :: i_thread
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :, :, :, :, :), POINTER   :: radii_pgf
    INTEGER, INTENT(IN)                      :: max_set, max_pgf
    REAL(dp)                                 :: eps_schwarz
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_pair_dist_radii', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, ikind, ipgf, isave, iset, jkind, jpgf, jset, la, &
      lb, ncoa, ncob, nkind, nseta, nsetb, sgfa, sgfb, stat
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure
    REAL(dp) :: cutoff, DATA(2,0:100), ff, max_contraction_a, &
      max_contraction_b, prefactor, R1, R_max, ra(3), rab(3), rab2, radius, &
      rap(3), rb(3), rp(3), x(2), zetp
    REAL(dp), DIMENSION(:), POINTER          :: set_radius_a, set_radius_b
    REAL(dp), DIMENSION(:, :), POINTER       :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                zeta, zetb
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis

    failure = .FALSE.

    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    error=error)

    nkind = SIZE(atomic_kind_set,1)
    isave=0
    ra = 0.0_dp
    rb = 0.0_dp
    DATA = 0.0_dp
    ALLOCATE(radii_pgf(max_pgf,max_pgf,max_set,max_set,nkind,nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO ikind=1,nkind
      NULLIFY(atom_kind,orb_basis)
      atom_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atom_kind,orb_basis_set=orb_basis)
      NULLIFY(la_max,la_min,npgfa,zeta)
      CALL get_gto_basis_set(gto_basis_set=orb_basis,&
                             lmax=la_max,&
                             lmin=la_min,&
                             npgf=npgfa,&
                             nset=nseta,&
                             zet=zeta,&
                             set_radius=set_radius_a,&
                             first_sgf=first_sgfa,&
                             sphi=sphi_a,&
                             nsgf_set=nsgfa,&
                             pgf_radius=rpgfa)
      DO jkind = 1,nkind
        NULLIFY(atom_kind,orb_basis)
        atom_kind => atomic_kind_set(jkind)
        CALL get_atomic_kind(atomic_kind=atom_kind,orb_basis_set=orb_basis)
        NULLIFY(lb_max,lb_min,npgfb,zetb)
        CALL get_gto_basis_set(gto_basis_set=orb_basis,&
                               lmax=lb_max,&
                               lmin=lb_min,&
                               npgf=npgfb,&
                               nset=nsetb,&
                               zet=zetb,&
                               set_radius=set_radius_b,&
                               first_sgf=first_sgfb,&
                               sphi=sphi_b,&
                               nsgf_set=nsgfb,&
                               pgf_radius=rpgfb)
        DO iset = 1,nseta
          ncoa = npgfa(iset)*ncoset(la_max(iset))
          sgfa = first_sgfa(1,iset)
          max_contraction_a = MAXVAL((/(SUM(ABS(sphi_a(1:ncoa,i))),i=sgfa,sgfa+nsgfa(iset)-1)/))
          DO jset = 1,nsetb
            ncob = npgfb(jset)*ncoset(lb_max(jset))
            sgfb = first_sgfb(1,jset)
            max_contraction_b = MAXVAL((/(SUM(ABS(sphi_b(1:ncob,i))),i=sgfb,sgfb+nsgfb(jset)-1)/))
            DATA = 0.0_dp
            DO ipgf = 1,npgfa(iset)
              DO jpgf = 1,npgfb(jset)
                DATA = 0.0_dp
                radius=rpgfa(ipgf,iset)+rpgfb(jpgf,jset)
                DO i=0,100
                  rb(1) = 0.0_dp + 0.01_dp * radius*i
                  R_max = 0.0_dp
                  DO la = la_min(iset),la_max(iset)
                    DO lb = lb_min(jset), lb_max(jset)
                      zetp = zeta(ipgf,iset) + zetb(jpgf,jset)
                      ff = zetb(jpgf,jset)/zetp
                      rab = 0.0_dp
                      rab(1) = rb(1)
                      rab2 = rb(1)**2
                      prefactor = EXP(-zeta(ipgf,iset)*ff*rab2)
                      rap(:) =ff*rab(:)
                      rp(:) = ra(:) + rap(:) 
                      rb(:) = ra(:) + rab(:)
                      cutoff = 1.0_dp
                      R1=exp_radius_very_extended(la,la,lb,lb,ra=ra,rb=rb,rp=rp,&
                             zetp=zetp,eps=eps_schwarz,prefactor=prefactor,cutoff=cutoff,epsin=1.0E-12_dp)
                      R_max = MAX(R_max,R1)
                    END DO
                  END DO
                  DATA(1,i) = rb(1)
                  DATA(2,i) = R_max
                END DO
                ! the radius can not be negative, we take that into account in the code as well by using a MAX
                ! the functional form we use for fitting does not seem particularly accurate 
                CALL optimize_it(DATA,x,0.0_dp)
                radii_pgf(jpgf,ipgf,jset,iset,jkind,ikind)%x = x 
              END DO !jpgf
            END DO !ipgf
          END DO
        END DO
      END DO
    END DO

    CALL timestop(handle)

  END SUBROUTINE calc_pair_dist_radii

!
!
! little drive routine for the powell minimizer
! data is the data to fit, x is of the form (x(1)*DATA(1)**2+x(2))
! only values of DATA(2) larger than fmin are taken into account
! it constructs an approximate upper bound of the fitted function
!
!
  SUBROUTINE optimize_it(DATA,x,fmin)

    REAL(KIND=dp), INTENT(IN)                :: DATA(2,0:100)
    REAL(KIND=dp), INTENT(OUT)               :: x(2)
    REAL(KIND=dp), INTENT(IN)                :: fmin

    INTEGER                                  :: i, k
    REAL(KIND=dp)                            :: f, large_weight, &
                                                small_weight, weight
    TYPE(opt_state_type)                     :: opt_state

! initial values

     x(1) = 0.0_dp
     x(2) = 0.0_dp

     ! we go in two steps, first we do the symmetric weight to get a good, unique initial guess
     ! we restart afterwards for the assym. 
     ! the assym function appears to have several local minima, depending on the data to fit
     ! the loop over k can make the switch gradual, but there is not much need, seemingly
     DO k=0,4,5

        small_weight=1.0_dp
        large_weight=small_weight*(10.0_dp**k)

        ! init opt run
        opt_state%state = 0
        opt_state%nvar = 2
        opt_state%iprint = 3
        opt_state%unit = default_output_unit
        opt_state%maxfun = 100000
        opt_state%rhobeg = 0.1_dp
        opt_state%rhoend = 0.000001_dp

        DO

          ! compute function
          IF ( opt_state%state == 2 ) THEN
            opt_state%f = 0.0_dp
            DO i=0,100
              f = x(1)*DATA(1,i)**2 +  x(2)
              IF( f > DATA(2,i) ) THEN
                weight = small_weight
              ELSE
                weight = large_weight
              END IF
              IF( DATA(2,i) > fmin ) opt_state%f = opt_state%f + weight * (f-DATA(2,i))**2
            END DO
          END IF

          IF ( opt_state%state == -1 ) EXIT
          CALL powell_optimize (opt_state%nvar, x, opt_state)
        END DO

        ! dealloc mem
        opt_state%state = 8
        CALL powell_optimize (opt_state%nvar, x, opt_state)
 
     ENDDO

  END SUBROUTINE optimize_it

! *****************************************************************************
!> \brief Given a 2d index pair, this function returns a 1d index pair for
!>        a symmetric upper triangle NxN matrix
!>        The compiler should inline this function, therefore it appears in
!>        several modules
!> \param i,j 2d index
!> \param N matrix size
!> \par History
!>      03.2009 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  PURE FUNCTION get_1D_idx(i,j,N)
    INTEGER, INTENT(IN)                      :: i, j
    INTEGER(int_8), INTENT(IN)               :: N
    INTEGER(int_8)                           :: get_1D_idx

    INTEGER(int_8)                           :: min_ij

    min_ij = MIN(i,j)
    get_1D_idx = min_ij*N + MAX(i,j) - (min_ij-1)*min_ij/2 - N

  END FUNCTION get_1D_idx

END MODULE hfx_screening_methods
