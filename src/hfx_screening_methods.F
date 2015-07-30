!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Several screening methods used in HFX calcualtions
!> \par History
!>      04.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
MODULE hfx_screening_methods
  USE ao_util,                         ONLY: exp_radius_very_extended
  USE basis_set_types,                 ONLY: gto_basis_set_type
  USE hfx_libint_interface,            ONLY: evaluate_eri_screen
  USE hfx_libint_wrapper_types,        ONLY: lib_int
  USE hfx_types,                       ONLY: hfx_basis_type,&
                                             hfx_p_kind,&
                                             hfx_potential_type,&
                                             hfx_screen_coeff_type,&
                                             log_zero,&
                                             powell_min_log
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE machine,                         ONLY: default_output_unit
  USE orbital_pointers,                ONLY: ncoset
  USE powell,                          ONLY: opt_state_type,&
                                             powell_optimize
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_pmax_mat, &
         calc_screening_functions,&
         calc_pair_dist_radii

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_screening_methods'

!***

  CONTAINS

! *****************************************************************************
!> \brief calculates max values of two-electron integrals in a quartet/shell
!>      w.r.t. different zetas using the library lib_int
!> \param lib ...
!> \param ra position
!> \param rb position
!> \param zeta zeta
!> \param zetb zeta
!> \param la_min angular momentum
!> \param la_max angular momentum
!> \param lb_min angular momentum
!> \param lb_max angular momentum
!> \param npgfa number of primitive cartesian gaussian in actual shell
!> \param npgfb number of primitive cartesian gaussian in actual shell
!> \param max_val schwarz screening value
!> \param potential_parameter contains info for libint
!> \param tmp_R_1 pgf_based screening coefficients
!> \param rab2 Squared Distance of centers ab
!> \param p_work ...
!> \par History
!>     03.2007 created [Manuel Guidon]
!>     02.2009 refactored [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE screen4(lib, ra, rb, zeta, zetb, &
                     la_min, la_max, lb_min, lb_max,&
                     npgfa, npgfb, &
                     max_val, potential_parameter, tmp_R_1,&
                     rab2, p_work)

    TYPE(lib_int)                            :: lib
    REAL(dp), INTENT(IN)                     :: ra(3), rb(3)
    REAL(dp), DIMENSION(:), INTENT(IN)       :: zeta, zetb
    INTEGER, INTENT(IN)                      :: la_min, la_max, lb_min, &
                                                lb_max, npgfa, npgfb
    REAL(dp), INTENT(INOUT)                  :: max_val
    TYPE(hfx_potential_type)                 :: potential_parameter
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :), POINTER               :: tmp_R_1
    REAL(dp)                                 :: rab2
    REAL(dp), DIMENSION(:), POINTER          :: p_work

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
                                     max_val_temp, potential_parameter, R1, R1, p_work)

            max_val = MAX(max_val, max_val_temp)
          END DO !lb
        END DO !la
      END DO !jpgf
    END DO !ipgf
  END SUBROUTINE screen4


! *****************************************************************************
!> \brief updates the maximum of the density matrix in compressed form for screening purposes
!> \param pmax_set buffer to store matrix
!> \param map_atom_to_kind_atom ...
!> \param set_offset ...
!> \param atomic_block_offset ...
!> \param pmax_atom ...
!> \param full_density_alpha ...
!> \param full_density_beta ...
!> \param image ...
!> \param natom ...
!> \param kind_of ...
!> \param basis_parameter ...
!> \param nkind ...
!> \param is_assoc_atomic_block ...
!> \par History
!>     09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      - updates for each pair of shells the maximum absolute value of p
! *****************************************************************************
  SUBROUTINE update_pmax_mat(pmax_set, map_atom_to_kind_atom, set_offset, atomic_block_offset, &
                             pmax_atom, full_density_alpha, full_density_beta, image, natom, &
                             kind_of,basis_parameter, &
                             nkind, is_assoc_atomic_block)

    TYPE(hfx_p_kind), DIMENSION(:), POINTER  :: pmax_set
    INTEGER, DIMENSION(:), POINTER           :: map_atom_to_kind_atom
    INTEGER, DIMENSION(:, :, :, :), POINTER  :: set_offset
    INTEGER, DIMENSION(:, :), POINTER        :: atomic_block_offset
    REAL(dp), DIMENSION(:, :), POINTER       :: pmax_atom, &
                                                full_density_alpha, &
                                                full_density_beta
    INTEGER, INTENT(IN)                      :: image, natom
    INTEGER                                  :: kind_of(*)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    INTEGER                                  :: nkind
    INTEGER, DIMENSION(:, :)                 :: is_assoc_atomic_block

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_pmax_mat', &
      routineP = moduleN//':'//routineN

    INTEGER :: act_atomic_block_offset, act_set_offset, handle, i, jatom, &
      jkind, jset, katom, kind_kind_idx, kkind, kset, mb, mc, nsetb, nsetc
    INTEGER, DIMENSION(:), POINTER           :: nsgfb, nsgfc
    REAL(dp)                                 :: pmax_tmp

    CALL timeset(routineN,handle)


    DO i = 1,SIZE(pmax_set)
      pmax_set(i)%p_kind = 0.0_dp
    END DO

    pmax_atom = log_zero

    DO jatom=1,natom
      jkind = kind_of(jatom)
      nsetb = basis_parameter(jkind)%nset
      nsgfb => basis_parameter(jkind)%nsgf

      DO katom = 1,natom
        IF( .NOT. is_assoc_atomic_block(jatom,katom) >=1 ) CYCLE
        kkind = kind_of(katom)
        IF( kkind < jkind ) CYCLE
        nsetc = basis_parameter(kkind)%nset
        nsgfc => basis_parameter(kkind)%nsgf
        act_atomic_block_offset = atomic_block_offset(katom,jatom)
        DO jset = 1,nsetb
          DO kset = 1,nsetc
            IF( katom >= jatom ) THEN
              pmax_tmp = 0.0_dp
              act_set_offset = set_offset(kset,jset,kkind,jkind)
              i = act_set_offset + act_atomic_block_offset - 1
              DO mc = 1,nsgfc(kset)
                DO mb = 1,nsgfb(jset)
                  pmax_tmp=MAX(pmax_tmp,ABS(full_density_alpha(i,image)))
                  IF( ASSOCIATED(full_density_beta) ) THEN
                    pmax_tmp=MAX(pmax_tmp,ABS(full_density_beta(i,image)))
                  END IF
                  i=i+1
                ENDDO
              ENDDO
              IF(pmax_tmp == 0.0_dp) THEN
                pmax_tmp = log_zero
              ELSE
                pmax_tmp = LOG10(pmax_tmp)
              END IF
              kind_kind_idx = INT(get_1D_idx(jkind,kkind,INT(nkind,int_8)))
              pmax_set(kind_kind_idx)%p_kind(jset, &
                       kset, &
                       map_atom_to_kind_atom(jatom), &
                       map_atom_to_kind_atom(katom)) = pmax_tmp
            ELSE
              pmax_tmp = 0.0_dp
              act_set_offset = set_offset(jset,kset,jkind,kkind)
              DO mc = 1,nsgfc(kset)
                i = act_set_offset + act_atomic_block_offset -1 + mc -1
                DO mb = 1,nsgfb(jset)
                  pmax_tmp=MAX(pmax_tmp,ABS(full_density_alpha(i,image)))
                  IF( ASSOCIATED(full_density_beta) ) THEN
                    pmax_tmp=MAX(pmax_tmp,ABS(full_density_beta(i,image)))
                  END IF
                  i= i + nsgfc(kset)
                ENDDO
              ENDDO
              IF(pmax_tmp == 0.0_dp) THEN
                pmax_tmp = log_zero
              ELSE
                pmax_tmp = LOG10(pmax_tmp)
              END IF
              kind_kind_idx = INT(get_1D_idx(jkind,kkind,INT(nkind,int_8)))
              pmax_set(kind_kind_idx)%p_kind(jset, &
                       kset, &
                       map_atom_to_kind_atom(jatom), &
                       map_atom_to_kind_atom(katom)) = pmax_tmp
            END IF
          END DO
        END DO
      END DO
    END DO


    DO jatom=1,natom
      jkind = kind_of(jatom)
      nsetb = basis_parameter(jkind)%nset
      DO katom = 1,natom
        IF( .NOT. is_assoc_atomic_block(jatom,katom) >=1 ) CYCLE
        kkind = kind_of(katom)
        IF( kkind < jkind ) CYCLE
        nsetc = basis_parameter(kkind)%nset
        pmax_tmp = log_zero
        DO jset = 1,nsetb
          DO kset = 1,nsetc
            kind_kind_idx = INT(get_1D_idx(jkind,kkind,INT(nkind,int_8)))
            pmax_tmp = MAX(pmax_set(kind_kind_idx)%p_kind(jset, &
                           kset, &
                           map_atom_to_kind_atom(jatom), &
                           map_atom_to_kind_atom(katom)), pmax_tmp)
          END DO
        END DO
        pmax_atom(jatom,katom) = pmax_tmp
        pmax_atom(katom,jatom) = pmax_tmp
      END DO
    END DO

    CALL timestop(handle)

  END SUBROUTINE update_pmax_mat

! *****************************************************************************
!> \brief calculates screening functions for schwarz screening
!> \param qs_env qs_env
!> \param basis_parameter ...
!> \param lib structure to libint
!> \param potential_parameter contains infos on potential
!> \param coeffs_set set based coefficients
!> \param coeffs_kind kind based coefficients
!> \param coeffs_pgf pgf based coefficients
!> \param radii_pgf coefficients for long-range screening
!> \param max_set Maximum Number of basis set sets in the system
!> \param max_pgf Maximum Number of basis set pgfs in the system
!> \param n_threads ...
!> \param i_thread Thread ID of current task
!> \param p_work ...
!> \param error ...
!> \par History
!>     02.2009 created [Manuel Guidon]
!> \author Manuel Guidon
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
! *****************************************************************************

  SUBROUTINE calc_screening_functions(qs_env, basis_parameter, lib, potential_parameter, &
                                      coeffs_set, coeffs_kind, coeffs_pgf, radii_pgf,&
                                      max_set, max_pgf, n_threads, i_thread,&
                                      p_work, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(lib_int)                            :: lib
    TYPE(hfx_potential_type)                 :: potential_parameter
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :, :, :), POINTER         :: coeffs_set
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :), POINTER               :: coeffs_kind
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :, :, :, :, :), POINTER   :: coeffs_pgf, radii_pgf
    INTEGER, INTENT(IN)                      :: max_set, max_pgf, n_threads, &
                                                i_thread
    REAL(dp), DIMENSION(:), POINTER          :: p_work
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
    REAL(dp) :: kind_radius_a, kind_radius_b, max_contraction_a, &
      max_contraction_b, max_val, max_val_temp, R1, ra(3), radius, rb(3), x(2)
    REAL(dp), DIMENSION(:), POINTER          :: set_radius_a, set_radius_b
    REAL(dp), DIMENSION(:, :), POINTER       :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                zeta, zetb
    REAL(dp), SAVE                           :: DATA(2,0:100)
    TYPE(gto_basis_set_type), POINTER        :: orb_basis
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :), POINTER               :: tmp_R_1
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_kind_type), POINTER              :: qs_kind

    failure = .FALSE.

!$OMP MASTER
    CALL timeset(routineN,handle)
!$OMP END MASTER

    CALL get_qs_env(qs_env=qs_env,&
                    qs_kind_set=qs_kind_set,&
                    error=error)

    nkind = SIZE(qs_kind_set,1)

!$OMP MASTER
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

!$OMP END MASTER
!$OMP BARRIER
    ra = 0.0_dp
    rb = 0.0_dp
    DO ikind=1,nkind
      NULLIFY(qs_kind,orb_basis)
      qs_kind => qs_kind_set(ikind)
      CALL get_qs_kind(qs_kind=qs_kind,orb_basis_set=orb_basis)
      NULLIFY(la_max,la_min,npgfa,zeta)

      la_max => basis_parameter(ikind)%lmax
      la_min => basis_parameter(ikind)%lmin
      npgfa => basis_parameter(ikind)%npgf
      nseta = basis_parameter(ikind)%nset
      zeta => basis_parameter(ikind)%zet
      set_radius_a => basis_parameter(ikind)%set_radius
      first_sgfa => basis_parameter(ikind)%first_sgf
      sphi_a => basis_parameter(ikind)%sphi
      nsgfa => basis_parameter(ikind)%nsgf
      rpgfa => basis_parameter(ikind)%pgf_radius

      DO jkind = 1,nkind
        NULLIFY(qs_kind,orb_basis)
        qs_kind => qs_kind_set(jkind)
        CALL get_qs_kind(qs_kind=qs_kind,orb_basis_set=orb_basis)
        NULLIFY(lb_max,lb_min,npgfb,zetb)

        lb_max => basis_parameter(jkind)%lmax
        lb_min => basis_parameter(jkind)%lmin
        npgfb => basis_parameter(jkind)%npgf
        nsetb = basis_parameter(jkind)%nset
        zetb => basis_parameter(jkind)%zet
        set_radius_b => basis_parameter(jkind)%set_radius
        first_sgfb => basis_parameter(jkind)%first_sgf
        sphi_b => basis_parameter(jkind)%sphi
        nsgfb => basis_parameter(jkind)%nsgf
        rpgfb => basis_parameter(jkind)%pgf_radius

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
                radius=rpgfa(ipgf,iset)+rpgfb(jpgf,jset)
                DO i=i_thread,100,n_threads
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
                                               max_val_temp, potential_parameter, R1, R1, p_work)
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
!$OMP BARRIER
!$OMP MASTER
                CALL optimize_it(DATA,x,powell_min_log)
                coeffs_pgf(jpgf,ipgf,jset,iset,jkind,ikind)%x = x
!$OMP END MASTER
!$OMP BARRIER

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

!$OMP MASTER
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
!$OMP END MASTER
!$OMP BARRIER
    ra = 0.0_dp
    rb = 0.0_dp
    DO ikind=1,nkind
      NULLIFY(qs_kind,orb_basis)
      qs_kind => qs_kind_set(ikind)
      CALL get_qs_kind(qs_kind=qs_kind,orb_basis_set=orb_basis)
      NULLIFY(la_max,la_min,npgfa,zeta)
!      CALL get_gto_basis_set(gto_basis_set=orb_basis,&
!                             lmax=la_max,&
!                             lmin=la_min,&
!                             npgf=npgfa,&
!                             nset=nseta,&
!                             zet=zeta,&
!                             set_radius=set_radius_a,&
!                             first_sgf=first_sgfa,&
!                             sphi=sphi_a,&
!                             nsgf_set=nsgfa)
      la_max => basis_parameter(ikind)%lmax
      la_min => basis_parameter(ikind)%lmin
      npgfa => basis_parameter(ikind)%npgf
      nseta = basis_parameter(ikind)%nset
      zeta => basis_parameter(ikind)%zet
      set_radius_a => basis_parameter(ikind)%set_radius
      first_sgfa => basis_parameter(ikind)%first_sgf
      sphi_a => basis_parameter(ikind)%sphi
      nsgfa => basis_parameter(ikind)%nsgf

      DO jkind = 1,nkind
        NULLIFY(qs_kind,orb_basis)
        qs_kind => qs_kind_set(jkind)
        CALL get_qs_kind(qs_kind=qs_kind,orb_basis_set=orb_basis)
        NULLIFY(lb_max,lb_min,npgfb,zetb)

        lb_max => basis_parameter(jkind)%lmax
        lb_min => basis_parameter(jkind)%lmin
        npgfb => basis_parameter(jkind)%npgf
        nsetb = basis_parameter(jkind)%nset
        zetb => basis_parameter(jkind)%zet
        set_radius_b => basis_parameter(jkind)%set_radius
        first_sgfb => basis_parameter(jkind)%first_sgf
        sphi_b => basis_parameter(jkind)%sphi
        nsgfb => basis_parameter(jkind)%nsgf

        DO iset = 1,nseta
          ncoa = npgfa(iset)*ncoset(la_max(iset))
          sgfa = first_sgfa(1,iset)
          max_contraction_a = MAXVAL((/(SUM(ABS(sphi_a(1:ncoa,i))),i=sgfa,sgfa+nsgfa(iset)-1)/))
          DO jset = 1,nsetb
            ncob = npgfb(jset)*ncoset(lb_max(jset))
            sgfb = first_sgfb(1,jset)
            max_contraction_b = MAXVAL((/(SUM(ABS(sphi_b(1:ncob,i))),i=sgfb,sgfb+nsgfb(jset)-1)/))
            radius=set_radius_a(iset)+set_radius_b(jset)
            tmp_R_1 => radii_pgf(:,:,jset,iset,jkind,ikind)
            DO i=i_thread,100,n_threads
              rb(1) = 0.0_dp + REAL(i,dp) * 0.01_dp * radius
              max_val = 0.0_dp
              CALL screen4(lib, ra, rb, &
                           zeta(:,iset), zetb(:,jset), &
                           la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                           npgfa(iset), npgfb(jset), &
                           max_val, potential_parameter, tmp_R_1, rb(1)**2, p_work)
              max_val = SQRT(max_val)
              max_val = max_val * max_contraction_a * max_contraction_b
              DATA(1,i) = rb(1)
              IF(max_val == 0.0_dp ) THEN
                DATA(2,i) = powell_min_log
              ELSE
                DATA(2,i) = LOG10((max_val))
              END IF
            END DO
!$OMP BARRIER
!$OMP MASTER
            CALL optimize_it(DATA,x,powell_min_log)
            coeffs_set(jset,iset,jkind,ikind)%x = x
!$OMP END MASTER
!$OMP BARRIER
          END DO
        END DO

      END DO
    END DO

    ! ** now kinds
!$OMP MASTER
    ALLOCATE(coeffs_kind(nkind,nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO ikind = 1,nkind
      DO jkind = 1,nkind
        coeffs_kind(jkind,ikind)%x(:) = 0.0_dp
      END DO
    END DO
!$OMP END MASTER
    ra = 0.0_dp
    rb = 0.0_dp
    DO ikind=1,nkind
      NULLIFY(qs_kind,orb_basis)
      qs_kind => qs_kind_set(ikind)
      CALL get_qs_kind(qs_kind=qs_kind,orb_basis_set=orb_basis)
      NULLIFY(la_max,la_min,npgfa,zeta)

      la_max => basis_parameter(ikind)%lmax
      la_min => basis_parameter(ikind)%lmin
      npgfa => basis_parameter(ikind)%npgf
      nseta = basis_parameter(ikind)%nset
      zeta => basis_parameter(ikind)%zet
      kind_radius_a = basis_parameter(ikind)%kind_radius
      first_sgfa => basis_parameter(ikind)%first_sgf
      sphi_a => basis_parameter(ikind)%sphi
      nsgfa => basis_parameter(ikind)%nsgf

      DO jkind = 1,nkind
        NULLIFY(qs_kind,orb_basis)
        qs_kind => qs_kind_set(jkind)
        CALL get_qs_kind(qs_kind=qs_kind,orb_basis_set=orb_basis)
        NULLIFY(lb_max,lb_min,npgfb,zetb)

        lb_max => basis_parameter(jkind)%lmax
        lb_min => basis_parameter(jkind)%lmin
        npgfb => basis_parameter(jkind)%npgf
        nsetb = basis_parameter(jkind)%nset
        zetb => basis_parameter(jkind)%zet
        kind_radius_b = basis_parameter(jkind)%kind_radius
        first_sgfb => basis_parameter(jkind)%first_sgf
        sphi_b => basis_parameter(jkind)%sphi
        nsgfb => basis_parameter(jkind)%nsgf


        radius=kind_radius_a + kind_radius_b
        DO iset=1,nseta
          ncoa = npgfa(iset)*ncoset(la_max(iset))
          sgfa = first_sgfa(1,iset)
          max_contraction_a = MAXVAL((/(SUM(ABS(sphi_a(1:ncoa,i))),i=sgfa,sgfa+nsgfa(iset)-1)/))
          DO jset=1,nsetb
            ncob = npgfb(jset)*ncoset(lb_max(jset))
            sgfb = first_sgfb(1,jset)
            max_contraction_b = MAXVAL((/(SUM(ABS(sphi_b(1:ncob,i))),i=sgfb,sgfb+nsgfb(jset)-1)/))
            DO i=i_thread,100,n_threads
              rb(1) = 0.0_dp + REAL(i,dp) * 0.01_dp * radius
              max_val = 0.0_dp
              tmp_R_1 => radii_pgf(:,:,jset,iset,jkind,ikind)
              CALL screen4(lib, ra, rb, &
                           zeta(:,iset), zetb(:,jset), &
                           la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                           npgfa(iset), npgfb(jset), &
                           max_val, potential_parameter, tmp_R_1, rb(1)**2, p_work)
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
!$OMP BARRIER
!$OMP MASTER
        CALL optimize_it(DATA,x,powell_min_log)
        coeffs_kind(jkind,ikind)%x = x
!$OMP END MASTER
!$OMP BARRIER
      END DO
    END DO

!$OMP MASTER
    CALL timestop(handle)
!$OMP END MASTER

  END SUBROUTINE calc_screening_functions

! *****************************************************************************
!> \brief calculates radius functions for longrange screening
!> \param qs_env qs_env
!> \param basis_parameter ...
!> \param radii_pgf pgf based coefficients
!> \param max_set Maximum Number of basis set sets in the system
!> \param max_pgf Maximum Number of basis set pgfs in the system
!> \param eps_schwarz ...
!> \param n_threads ...
!> \param i_thread Thread ID of current task
!> \param error ...
!> \par History
!>     02.2009 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      This routine calculates the pair-distribution radius of a product
!>      Gaussian and uses the powell optimiztion routines in order to fit
!>      the results in the following form:
!>
!>                 (ab| = (ab(Rab) = c2*Rab^2 + c0
!>
! *****************************************************************************

  SUBROUTINE calc_pair_dist_radii(qs_env, basis_parameter, &
                                  radii_pgf, max_set, max_pgf, eps_schwarz, &
                                  n_threads, i_thread, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(hfx_screen_coeff_type), &
      DIMENSION(:, :, :, :, :, :), POINTER   :: radii_pgf
    INTEGER, INTENT(IN)                      :: max_set, max_pgf
    REAL(dp)                                 :: eps_schwarz
    INTEGER, INTENT(IN)                      :: n_threads, i_thread
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_pair_dist_radii', &
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
    REAL(dp) :: cutoff, ff, max_contraction_a, max_contraction_b, prefactor, &
      R1, R_max, ra(3), rab(3), rab2, radius, rap(3), rb(3), rp(3), x(2), zetp
    REAL(dp), DIMENSION(:, :), POINTER       :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                zeta, zetb
    REAL(dp), SAVE                           :: DATA(2,0:100)
    TYPE(gto_basis_set_type), POINTER        :: orb_basis
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_kind_type), POINTER              :: qs_kind

    failure = .FALSE.


!$OMP MASTER
    CALL timeset(routineN,handle)
!$OMP END MASTER
    CALL get_qs_env(qs_env=qs_env,&
                    qs_kind_set=qs_kind_set,&
                    error=error)

    nkind = SIZE(qs_kind_set,1)
    ra = 0.0_dp
    rb = 0.0_dp
!$OMP MASTER
    ALLOCATE(radii_pgf(max_pgf,max_pgf,max_set,max_set,nkind,nkind),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO ikind = 1,nkind
      DO jkind = 1,nkind
        DO iset=1,max_set
          DO jset=1,max_set
            DO ipgf=1,max_pgf
              DO jpgf=1,max_pgf
                radii_pgf(jpgf,ipgf,jset,iset,jkind,ikind)%x(:) = 0.0_dp
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    DATA = 0.0_dp
!$OMP END MASTER
!$OMP BARRIER

    DO ikind=1,nkind
      NULLIFY(qs_kind,orb_basis)
      qs_kind => qs_kind_set(ikind)
      CALL get_qs_kind(qs_kind=qs_kind,orb_basis_set=orb_basis)
      NULLIFY(la_max,la_min,npgfa,zeta)

      la_max => basis_parameter(ikind)%lmax
      la_min => basis_parameter(ikind)%lmin
      npgfa => basis_parameter(ikind)%npgf
      nseta = basis_parameter(ikind)%nset
      zeta => basis_parameter(ikind)%zet
      first_sgfa => basis_parameter(ikind)%first_sgf
      sphi_a => basis_parameter(ikind)%sphi
      nsgfa => basis_parameter(ikind)%nsgf
      rpgfa => basis_parameter(ikind)%pgf_radius


      DO jkind = 1,nkind
        NULLIFY(qs_kind,orb_basis)
        qs_kind => qs_kind_set(jkind)
        CALL get_qs_kind(qs_kind=qs_kind,orb_basis_set=orb_basis)
        NULLIFY(lb_max,lb_min,npgfb,zetb)

        lb_max => basis_parameter(jkind)%lmax
        lb_min => basis_parameter(jkind)%lmin
        npgfb => basis_parameter(jkind)%npgf
        nsetb = basis_parameter(jkind)%nset
        zetb => basis_parameter(jkind)%zet
        first_sgfb => basis_parameter(jkind)%first_sgf
        sphi_b => basis_parameter(jkind)%sphi
        nsgfb => basis_parameter(jkind)%nsgf
        rpgfb => basis_parameter(jkind)%pgf_radius

        DO iset = 1,nseta
          ncoa = npgfa(iset)*ncoset(la_max(iset))
          sgfa = first_sgfa(1,iset)
          max_contraction_a = MAXVAL((/(SUM(ABS(sphi_a(1:ncoa,i))),i=sgfa,sgfa+nsgfa(iset)-1)/))
          DO jset = 1,nsetb
            ncob = npgfb(jset)*ncoset(lb_max(jset))
            sgfb = first_sgfb(1,jset)
            max_contraction_b = MAXVAL((/(SUM(ABS(sphi_b(1:ncob,i))),i=sgfb,sgfb+nsgfb(jset)-1)/))
            DO ipgf = 1,npgfa(iset)
              DO jpgf = 1,npgfb(jset)
                radius=rpgfa(ipgf,iset)+rpgfb(jpgf,jset)
                DO i = i_thread, 100, n_threads
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
!$OMP BARRIER
!$OMP MASTER
                CALL optimize_it(DATA,x,0.0_dp)
                radii_pgf(jpgf,ipgf,jset,iset,jkind,ikind)%x = x
!$OMP END MASTER
!$OMP BARRIER
              END DO !jpgf
            END DO !ipgf
          END DO
        END DO
      END DO
    END DO
!$OMP MASTER
    CALL timestop(handle)
!$OMP END MASTER
END SUBROUTINE calc_pair_dist_radii

!
!
! little driver routine for the powell minimizer
! data is the data to fit, x is of the form (x(1)*DATA(1)**2+x(2))
! only values of DATA(2) larger than fmin are taken into account
! it constructs an approximate upper bound of the fitted function
!
!
! *****************************************************************************
!> \brief ...
!> \param DATA ...
!> \param x ...
!> \param fmin ...
! *****************************************************************************
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
     DO k=0,4,4

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
!> \param i 2d index
!> \param j 2d index
!> \param N matrix size
!> \retval get_1D_idx ...
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
