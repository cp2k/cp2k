!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of charge equilibration method
!> \author JGH
! **************************************************************************************************
MODULE eeq_method
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell,&
                                              pbc,&
                                              plane_distance
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_fm_basic_linalg,              ONLY: cp_fm_invert,&
                                              cp_fm_matvec,&
                                              cp_fm_solve
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_unit_nr
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE ewald_environment_types,         ONLY: ewald_env_create,&
                                              ewald_env_get,&
                                              ewald_env_set,&
                                              ewald_environment_type,&
                                              read_ewald_section_tb
   USE ewald_pw_types,                  ONLY: ewald_pw_create,&
                                              ewald_pw_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: oorootpi
   USE message_passing,                 ONLY: mp_para_env_type
   USE molecule_types,                  ONLY: molecule_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: bohr
   USE pw_poisson_types,                ONLY: do_ewald_spme
   USE qs_dispersion_pairpot,           ONLY: cnumber_init,&
                                              dcnum_type
   USE qs_dispersion_types,             ONLY: qs_dispersion_release,&
                                              qs_dispersion_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_neighbor_lists,               ONLY: atom2d_build,&
                                              atom2d_cleanup,&
                                              build_neighbor_lists,&
                                              local_atoms_type,&
                                              pair_radius_setup
   USE spme,                            ONLY: spme_potential
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'eeq_method'

   ! ==========================================================================
   ! 2019-2020 Sebastian Ehlert : EEQ gfn0-xTB parameters
   INTEGER, PARAMETER                                    :: maxElem = 86
   ! COULOMB DATA
   !> Electronegativity for the EEQ model
   REAL(KIND=dp), PARAMETER :: eeqChi(1:maxElem) = [&
      & 1.2500000_dp, 1.2912463_dp, 0.8540050_dp, 1.1723939_dp, 1.1094487_dp, &
      & 1.3860275_dp, 1.5341534_dp, 1.5378836_dp, 1.5890750_dp, 1.2893646_dp, &
      & 0.7891208_dp, 0.9983021_dp, 0.9620847_dp, 1.0441134_dp, 1.4789559_dp, &
      & 1.3926377_dp, 1.4749100_dp, 1.2250415_dp, 0.8162292_dp, 1.1252036_dp, &
      & 0.9641451_dp, 0.8810155_dp, 0.9741986_dp, 1.1029038_dp, 1.0076949_dp, &
      & 0.7744353_dp, 0.7554040_dp, 1.0182630_dp, 1.0316167_dp, 1.6317474_dp, &
      & 1.1186739_dp, 1.0345958_dp, 1.3090772_dp, 1.4119283_dp, 1.4500674_dp, &
      & 1.1746889_dp, 0.6686200_dp, 1.0744648_dp, 0.9107813_dp, 0.7876056_dp, &
      & 1.0039889_dp, 0.9225265_dp, 0.9035515_dp, 1.0332301_dp, 1.0293975_dp, &
      & 1.0549549_dp, 1.2356867_dp, 1.2793315_dp, 1.1145650_dp, 1.1214927_dp, &
      & 1.2123167_dp, 1.4003158_dp, 1.4255511_dp, 1.1640198_dp, 0.4685133_dp, &
      & 1.0687873_dp, 0.9335398_dp, 1.0573550_dp, 1.0532043_dp, 1.0490537_dp, &
      & 1.0449031_dp, 1.0407524_dp, 1.0366018_dp, 1.0324512_dp, 1.0283005_dp, &
      & 1.0241499_dp, 1.0199992_dp, 1.0158486_dp, 1.0116980_dp, 1.0075473_dp, &
      & 1.0033967_dp, 0.8612827_dp, 1.0422031_dp, 0.7633168_dp, 0.6019707_dp, &
      & 0.7499393_dp, 0.9511744_dp, 0.9357472_dp, 1.3555382_dp, 1.2006726_dp, &
      & 1.2092025_dp, 1.1736669_dp, 1.1936584_dp, 1.3045488_dp, 1.1964604_dp, &
      & 1.2653792_dp]
   !> Chemical hardness for the EEQ model
   REAL(KIND=dp), PARAMETER :: eeqGam(1:maxElem) = [&
      &-0.3023159_dp, 0.7743046_dp, 0.5303164_dp, 0.2176474_dp, 0.1956176_dp, &
      & 0.0308461_dp, 0.0559522_dp, 0.0581228_dp, 0.1574017_dp, 0.6825784_dp, &
      & 0.3922376_dp, 0.5581866_dp, 0.3017510_dp, 0.1039137_dp, 0.2124917_dp, &
      & 0.0580720_dp, 0.2537467_dp, 0.5780354_dp, 0.3920658_dp, -0.0024897_dp, &
      &-0.0061520_dp, 0.1663252_dp, 0.1051751_dp, 0.0009900_dp, 0.0976543_dp, &
      & 0.0612028_dp, 0.0561526_dp, 0.0899774_dp, 0.1313171_dp, 0.5728071_dp, &
      & 0.1741615_dp, 0.2671888_dp, 0.2351989_dp, 0.0718104_dp, 0.3458143_dp, &
      & 0.8203265_dp, 0.4287770_dp, 0.2667067_dp, 0.0873658_dp, 0.0599431_dp, &
      & 0.1581972_dp, 0.1716374_dp, 0.2721649_dp, 0.2817608_dp, 0.1391572_dp, &
      & 0.1175925_dp, 0.2316104_dp, 0.2256303_dp, 0.1230459_dp, 0.0141941_dp, &
      & 0.0188612_dp, 0.0230207_dp, 0.3644113_dp, 0.1668461_dp, 0.5167533_dp, &
      & 0.1979578_dp, 0.0345176_dp, 0.0240233_dp, 0.0246333_dp, 0.0252433_dp, &
      & 0.0258532_dp, 0.0264632_dp, 0.0270732_dp, 0.0276832_dp, 0.0282931_dp, &
      & 0.0289031_dp, 0.0295131_dp, 0.0301230_dp, 0.0307330_dp, 0.0313430_dp, &
      & 0.0319529_dp, 0.0262881_dp, 0.1715396_dp, 0.1803633_dp, 0.3631824_dp, &
      & 0.3010980_dp, 0.1100299_dp, 0.0277514_dp, 0.0554975_dp, 0.7723231_dp, &
      & 0.1287718_dp, 0.1034598_dp, 0.0114935_dp, 0.0160842_dp, 0.3369611_dp, &
      & 0.1844179_dp]

   !> Coordination number dependence of the EN in the EEQ model
   REAL(KIND=dp), PARAMETER :: eeqkCN(1:maxElem) = [&
      & 0.0248762_dp, 0.1342276_dp, 0.0103048_dp, -0.0352374_dp, -0.0980031_dp, &
      & 0.0643920_dp, 0.1053273_dp, 0.1394809_dp, 0.1276675_dp, -0.1081936_dp, &
      &-0.0008132_dp, -0.0279860_dp, -0.0521436_dp, -0.0257206_dp, 0.1651461_dp, &
      & 0.0914418_dp, 0.1213634_dp, -0.0636298_dp, -0.0045838_dp, 0.0007509_dp, &
      &-0.0307730_dp, -0.0286150_dp, -0.0341465_dp, -0.0419655_dp, -0.0088536_dp, &
      &-0.1001069_dp, -0.1190502_dp, -0.0726233_dp, -0.0219233_dp, 0.0641913_dp, &
      &-0.0103130_dp, 0.0262628_dp, 0.0222202_dp, 0.0709954_dp, 0.0422244_dp, &
      &-0.0308245_dp, 0.0086249_dp, -0.0237146_dp, -0.0721798_dp, -0.0848810_dp, &
      &-0.0402828_dp, -0.0372396_dp, -0.0027043_dp, 0.0525839_dp, 0.0051192_dp, &
      & 0.0188401_dp, 0.0103998_dp, 0.0000549_dp, 0.0087717_dp, -0.0237228_dp, &
      & 0.0169656_dp, 0.0924186_dp, 0.0352884_dp, -0.0091444_dp, 0.0192916_dp, &
      &-0.0154483_dp, -0.0736833_dp, -0.0064191_dp, -0.0093012_dp, -0.0121833_dp, &
      &-0.0150654_dp, -0.0179475_dp, -0.0208296_dp, -0.0237117_dp, -0.0265938_dp, &
      &-0.0294759_dp, -0.0323580_dp, -0.0352400_dp, -0.0381221_dp, -0.0410042_dp, &
      &-0.0438863_dp, -0.0894776_dp, -0.0333583_dp, -0.0154963_dp, -0.0121092_dp, &
      &-0.0744239_dp, 0.0050138_dp, -0.0153757_dp, -0.0029221_dp, 0.0239125_dp, &
      & 0.0183012_dp, -0.0238011_dp, -0.0268025_dp, 0.0136505_dp, -0.0132199_dp, &
      &-0.0439890_dp]
   !> Charge width in the EEQ model
   REAL(KIND=dp), PARAMETER :: eeqAlp(1:maxElem) = [&
      & 0.7490227_dp, 0.4196569_dp, 1.4256190_dp, 2.0698743_dp, 1.7358798_dp, &
      & 1.8288757_dp, 1.9346081_dp, 1.6974795_dp, 0.8169179_dp, 0.6138441_dp, &
      & 1.7294046_dp, 1.7925036_dp, 1.2156739_dp, 1.5314457_dp, 1.3730859_dp, &
      & 1.7936326_dp, 2.4255996_dp, 1.5891656_dp, 2.1829647_dp, 1.4177623_dp, &
      & 1.5181399_dp, 1.9919805_dp, 1.7171675_dp, 2.0655063_dp, 1.3318009_dp, &
      & 1.3660068_dp, 1.5694128_dp, 1.2762644_dp, 1.0039549_dp, 0.7338863_dp, &
      & 3.2596250_dp, 1.7530299_dp, 1.5281792_dp, 2.1837813_dp, 2.1642027_dp, &
      & 2.7280594_dp, 0.7838049_dp, 1.4274742_dp, 1.8023947_dp, 1.6093288_dp, &
      & 1.3834349_dp, 1.1740977_dp, 1.5768259_dp, 1.3205263_dp, 1.4259466_dp, &
      & 1.1499748_dp, 0.7013009_dp, 1.2374416_dp, 1.3799991_dp, 1.8528424_dp, &
      & 1.8497568_dp, 2.0159294_dp, 1.2903708_dp, 2.0199161_dp, 0.9530522_dp, &
      & 1.5015025_dp, 2.1917012_dp, 1.9134370_dp, 1.9897910_dp, 2.0661450_dp, &
      & 2.1424991_dp, 2.2188531_dp, 2.2952071_dp, 2.3715611_dp, 2.4479151_dp, &
      & 2.5242691_dp, 2.6006231_dp, 2.6769771_dp, 2.7533312_dp, 2.8296852_dp, &
      & 2.9060392_dp, 1.6423047_dp, 1.3567622_dp, 1.8966648_dp, 0.8253100_dp, &
      & 0.7412219_dp, 1.0350883_dp, 0.9692278_dp, 1.0048087_dp, 2.3138674_dp, &
      & 2.8055966_dp, 3.0968677_dp, 1.6597596_dp, 3.2191575_dp, 1.5388148_dp, &
      & 2.1222013_dp]

   ! covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197)
   ! values for metals decreased by 10 %
   REAL(KIND=dp), PARAMETER :: rcov(1:maxElem) = [&
      & 0.32_dp, 0.46_dp, 1.20_dp, 0.94_dp, 0.77_dp, 0.75_dp, 0.71_dp, 0.63_dp, &
      & 0.64_dp, 0.67_dp, 1.40_dp, 1.25_dp, 1.13_dp, 1.04_dp, 1.10_dp, 1.02_dp, &
      & 0.99_dp, 0.96_dp, 1.76_dp, 1.54_dp, 1.33_dp, 1.22_dp, 1.21_dp, 1.10_dp, &
      & 1.07_dp, 1.04_dp, 1.00_dp, 0.99_dp, 1.01_dp, 1.09_dp, 1.12_dp, 1.09_dp, &
      & 1.15_dp, 1.10_dp, 1.14_dp, 1.17_dp, 1.89_dp, 1.67_dp, 1.47_dp, 1.39_dp, &
      & 1.32_dp, 1.24_dp, 1.15_dp, 1.13_dp, 1.13_dp, 1.08_dp, 1.15_dp, 1.23_dp, &
      & 1.28_dp, 1.26_dp, 1.26_dp, 1.23_dp, 1.32_dp, 1.31_dp, 2.09_dp, 1.76_dp, &
      & 1.62_dp, 1.47_dp, 1.58_dp, 1.57_dp, 1.56_dp, 1.55_dp, 1.51_dp, 1.52_dp, &
      & 1.51_dp, 1.50_dp, 1.49_dp, 1.49_dp, 1.48_dp, 1.53_dp, 1.46_dp, 1.37_dp, &
      & 1.31_dp, 1.23_dp, 1.18_dp, 1.16_dp, 1.11_dp, 1.12_dp, 1.13_dp, 1.32_dp, &
      & 1.30_dp, 1.30_dp, 1.36_dp, 1.31_dp, 1.38_dp, 1.42_dp]

   PUBLIC :: eeq_solver

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param charges ...
! **************************************************************************************************
   SUBROUTINE eeq_charges(qs_env, charges)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: charges

      CHARACTER(len=*), PARAMETER                        :: routineN = 'eeq_charges'

      INTEGER                                            :: handle, iatom, ikind, iunit, jkind, &
                                                            natom, nkind, za, zb
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      INTEGER, DIMENSION(3)                              :: periodic
      LOGICAL                                            :: do_ewald, do_sparse
      REAL(KIND=dp)                                      :: ala, alb, eeq_energy, kappa, lambda, &
                                                            scn, totalcharge, xi
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: chia, cnumbers, gam
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gab
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell, cell_ref
      TYPE(dcnum_type), ALLOCATABLE, DIMENSION(:)        :: dcnum
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: ewald_section, poisson_section, &
                                                            print_section

      CALL timeset(routineN, handle)

      iunit = cp_logger_get_default_unit_nr()

      CALL get_qs_env(qs_env, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set, &
                      particle_set=particle_set, &
                      para_env=para_env, &
                      cell=cell, &
                      dft_control=dft_control)
      CALL get_qs_env(qs_env, nkind=nkind, natom=natom)

      totalcharge = dft_control%charge

      CALL get_cnumbers(qs_env, cnumbers)

      ! gamma[a,b]
      ALLOCATE (gab(nkind, nkind), gam(nkind))
      gab = 0.0_dp
      gam = 0.0_dp
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), zatom=za)
         ala = eeqAlp(za)
         gam(ikind) = eeqGam(za)
         DO jkind = 1, nkind
            CALL get_qs_kind(qs_kind_set(jkind), zatom=zb)
            alb = eeqAlp(zb)
            !
            gab(ikind, jkind) = SQRT(1._dp/(ala*ala + alb*alb))
            !
         END DO
      END DO

      ! Chi[a,a]
      ALLOCATE (chia(natom))
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind_set(ikind), zatom=za)
         xi = eeqChi(za)
         kappa = eeqkCN(za)
         !
         scn = SQRT(cnumbers(iatom)) + 1.0e-14_dp
         chia(iatom) = xi - kappa*scn
         !
      END DO

      CALL cnumber_release(cnumbers, dcnum, .FALSE.)

      CALL get_cell(cell, periodic=periodic)
      do_ewald = .NOT. ALL(periodic == 0)
      do_sparse = .FALSE.
      IF (do_ewald) THEN
         ALLOCATE (ewald_env)
         CALL ewald_env_create(ewald_env, para_env)
         poisson_section => section_vals_get_subs_vals(qs_env%input, "DFT%POISSON")
         CALL ewald_env_set(ewald_env, poisson_section=poisson_section)
         ewald_section => section_vals_get_subs_vals(poisson_section, "EWALD")
         print_section => section_vals_get_subs_vals(qs_env%input, "PRINT%GRID_INFORMATION")
         CALL get_qs_env(qs_env, cell_ref=cell_ref)
         CALL read_ewald_section_tb(ewald_env, ewald_section, cell_ref%hmat)
         ALLOCATE (ewald_pw)
         CALL ewald_pw_create(ewald_pw, ewald_env, cell, cell_ref, print_section=print_section)
         !
         CALL eeq_solver(qs_env, charges, lambda, eeq_energy, &
                         particle_set, kind_of, cell, chia, gam, gab, &
                         totalcharge=totalcharge, ewald=do_ewald, sparse=do_sparse, &
                         ewald_env=ewald_env, ewald_pw=ewald_pw, iounit=iunit)
         !
         CALL ewald_env_release(ewald_env)
         CALL ewald_pw_release(ewald_pw)
      ELSE
         CALL eeq_solver(qs_env, charges, lambda, eeq_energy, &
                         particle_set, kind_of, cell, chia, gam, gab, &
                         totalcharge=totalcharge, iounit=iunit)
      END IF

      DEALLOCATE (gab, gam, chia)

      CALL timestop(handle)

   END SUBROUTINE eeq_charges

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param cnumbers ...
! **************************************************************************************************
   SUBROUTINE get_cnumbers(qs_env, cnumbers)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cnumbers

      INTEGER                                            :: ikind, natom, nkind, za
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: default_present
      REAL(KIND=dp)                                      :: subcells
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: c_radius
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: pair_radius
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dcnum_type), ALLOCATABLE, DIMENSION(:)        :: dcnum
      TYPE(distribution_1d_type), POINTER                :: distribution_1d
      TYPE(distribution_2d_type), POINTER                :: distribution_2d
      TYPE(local_atoms_type), ALLOCATABLE, DIMENSION(:)  :: atom2d
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_cn
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_dispersion_type), POINTER                  :: disp
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL get_qs_env(qs_env, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set, &
                      particle_set=particle_set, &
                      cell=cell, &
                      distribution_2d=distribution_2d, &
                      local_particles=distribution_1d, &
                      molecule_set=molecule_set)
      CALL get_qs_env(qs_env, nkind=nkind, natom=natom)

      ! Check for dispersion_env and sab_cn needed for cnumbers
      ALLOCATE (disp)
      disp%k1 = 16.0_dp
      disp%k2 = 4._dp/3._dp
      disp%eps_cn = 1.E-6_dp
      disp%max_elem = maxElem
      ALLOCATE (disp%rcov(maxElem))
      disp%rcov(1:maxElem) = bohr*disp%k2*rcov(1:maxElem)
      ! Build the neighbor lists for the CN
      NULLIFY (sab_cn)
      ALLOCATE (c_radius(nkind), default_present(nkind), pair_radius(nkind, nkind))
      c_radius(:) = 0.0_dp
      default_present = .TRUE.
      DO ikind = 1, nkind
         CALL get_atomic_kind(atomic_kind_set(ikind), z=za)
         c_radius(ikind) = 4._dp*rcov(za)*bohr
      END DO
      ALLOCATE (atom2d(nkind))
      CALL atom2d_build(atom2d, distribution_1d, distribution_2d, atomic_kind_set, &
                        molecule_set, .FALSE., particle_set=particle_set)
      CALL pair_radius_setup(default_present, default_present, c_radius, c_radius, pair_radius)
      CALL build_neighbor_lists(sab_cn, particle_set, atom2d, cell, pair_radius, &
                                subcells=subcells, operator_type="PP", nlname="sab_cn")
      disp%sab_cn => sab_cn
      DEALLOCATE (c_radius, pair_radius, default_present)
      CALL atom2d_cleanup(atom2d)

      ! Calculate coordination numbers
      CALL cnumber_init(qs_env, cnumbers, dcnum, 2, .FALSE.)

      CALL qs_dispersion_release(disp)

   END SUBROUTINE get_cnumbers

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param charges ...
!> \param lambda ...
!> \param eeq_energy ...
!> \param particle_set ...
!> \param kind_of ...
!> \param cell ...
!> \param chia ...
!> \param gam ...
!> \param gab ...
!> \param totalcharge ...
!> \param ewald ...
!> \param sparse ...
!> \param ewald_env ...
!> \param ewald_pw ...
!> \param iounit ...
! **************************************************************************************************
   SUBROUTINE eeq_solver(qs_env, charges, lambda, eeq_energy, particle_set, kind_of, cell, &
                         chia, gam, gab, totalcharge, ewald, sparse, ewald_env, ewald_pw, iounit)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: charges
      REAL(KIND=dp), INTENT(INOUT)                       :: lambda, eeq_energy
      TYPE(particle_type), DIMENSION(:), INTENT(IN)      :: particle_set
      INTEGER, DIMENSION(:), INTENT(IN)                  :: kind_of
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: chia, gam
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: gab
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: totalcharge
      LOGICAL, INTENT(IN), OPTIONAL                      :: ewald, sparse
      TYPE(ewald_environment_type), OPTIONAL, POINTER    :: ewald_env
      TYPE(ewald_pw_type), OPTIONAL, POINTER             :: ewald_pw
      INTEGER, INTENT(IN), OPTIONAL                      :: iounit

      CHARACTER(len=*), PARAMETER                        :: routineN = 'eeq_solver'
      REAL(KIND=dp), PARAMETER                           :: eps_pcg = 1.E-7_dp

      INTEGER                                            :: ewald_type, handle, i, ia, iac, iar, ic, &
                                                            ikind, ir, iunit, ix, iy, iz, jkind, &
                                                            natom, ncloc, ncvloc, nkind, nrloc, &
                                                            nrvloc, ns
      INTEGER, DIMENSION(3)                              :: cvec, ncell, periodic
      INTEGER, DIMENSION(:), POINTER                     :: cind, cvind, rind, rvind
      LOGICAL                                            :: do_ewald, do_sparse
      REAL(KIND=dp)                                      :: ad, alk, alpha, betak, deth, dr, grc, &
                                                            grc1, grc2, papk, qtot, rcut, rdl, &
                                                            res, rmax, rzk, rzkm, rzkp, xr
      REAL(KIND=dp), DIMENSION(3)                        :: ri, rij, rijl, rj
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      REAL(KIND=dp), DIMENSION(:), POINTER               :: av0, pv0, rv0, xv0, zv0
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: mat_struct, vec_struct
      TYPE(cp_fm_type)                                   :: eeq_mat, mmat, pmat, rhs_vec
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      do_sparse = .FALSE.
      IF (PRESENT(sparse)) do_sparse = sparse
      !
      do_ewald = .FALSE.
      IF (PRESENT(ewald)) do_ewald = ewald
      !
      qtot = 0.0_dp
      IF (PRESENT(totalcharge)) qtot = totalcharge
      !
      iunit = -1
      IF (PRESENT(iounit)) iunit = iounit

      ! sparse option NYA
      IF (do_sparse) THEN
         CALL cp_abort(__LOCATION__, "EEQ: Sparse option not yet available")
      END IF
      !
      natom = SIZE(particle_set)
      nkind = SIZE(gam)
      CALL get_qs_env(qs_env=qs_env, para_env=para_env, blacs_env=blacs_env)
      !
      IF (do_ewald) THEN
         CPASSERT(PRESENT(ewald_env))
         CPASSERT(PRESENT(ewald_pw))
         CALL get_cell(cell=cell, deth=deth)
         CALL ewald_env_get(ewald_env, alpha=alpha, rcut=rcut, ewald_type=ewald_type)
         ad = 2.0_dp*alpha*oorootpi
         IF (ewald_type /= do_ewald_spme) THEN
            CALL cp_abort(__LOCATION__, "Only SPME Ewald method available with EEQ.")
         END IF
         ! check alpha for MIC
         !
         rmax = 2.0_dp*rcut
         ! max cells used
         CALL get_cell(cell, h=hmat, periodic=periodic)
         ncell(1) = CEILING(rmax/plane_distance(1, 0, 0, cell))
         ncell(2) = CEILING(rmax/plane_distance(0, 1, 0, cell))
         ncell(3) = CEILING(rmax/plane_distance(0, 0, 1, cell))
         IF (periodic(1) == 0) ncell(1) = 0
         IF (periodic(2) == 0) ncell(2) = 0
         IF (periodic(3) == 0) ncell(3) = 0
         !
         ns = natom
         CALL cp_fm_struct_create(mat_struct, context=blacs_env, para_env=para_env, &
                                  nrow_global=ns, ncol_global=ns)
         CALL cp_fm_create(eeq_mat, mat_struct)
         CALL cp_fm_set_all(eeq_mat, 0.0_dp, 0.0_dp)
         CALL cp_fm_create(pmat, mat_struct)
         CALL cp_fm_set_all(pmat, 0.0_dp, 0.0_dp)
         CALL cp_fm_create(mmat, mat_struct)
         CALL cp_fm_set_all(mmat, 0.0_dp, 0.0_dp)
         CALL cp_fm_get_info(eeq_mat, nrow_local=nrloc, ncol_local=ncloc, &
                             row_indices=rind, col_indices=cind)
         ! response matrix
         DO ir = 1, nrloc
            iar = rind(ir)
            ikind = kind_of(iar)
            ri(1:3) = particle_set(iar)%r(1:3)
            DO ic = 1, ncloc
               iac = cind(ic)
               jkind = kind_of(iac)
               rj(1:3) = particle_set(iac)%r(1:3)
               rij(1:3) = ri(1:3) - rj(1:3)
               rij = pbc(rij, cell)
               DO ix = -ncell(1), ncell(1)
                  DO iy = -ncell(2), ncell(2)
                     DO iz = -ncell(3), ncell(3)
                        cvec = [ix, iy, iz]
                        rijl = rij + MATMUL(hmat, cvec)
                        dr = SQRT(SUM(rijl**2))
                        IF (dr > rmax) CYCLE
                        IF (iar == iac) THEN
                           grc2 = gam(ikind) + 2.0_dp*gab(ikind, ikind)*oorootpi
                           grc1 = grc2 - ad
                        ELSE
                           grc1 = erf(gab(ikind, jkind)*dr)/dr - erf(alpha*dr)/dr
                           IF (dr < rmax) THEN
                              grc2 = erf(gab(ikind, jkind)*dr)/dr
                           ELSE
                              grc2 = 0.0_dp
                           END IF
                        END IF
                        eeq_mat%local_data(ir, ic) = eeq_mat%local_data(ir, ic) + grc1
                        pmat%local_data(ir, ic) = pmat%local_data(ir, ic) + grc2
                     END DO
                  END DO
               END DO
            END DO
         END DO
         !
         ! preconditioner invers
         CALL cp_fm_invert(pmat, mmat)
         !
         ALLOCATE (av0(natom), xv0(natom), rv0(natom), pv0(natom), zv0(natom))
         ! Initial guess
         CALL cp_fm_matvec(mmat, -chia, xv0)
         grc = (SUM(xv0) - qtot)/REAL(natom, KIND=dp)
         xv0 = xv0 - grc
         CALL apply_potential(eeq_mat, ewald_env, ewald_pw, cell, particle_set, xv0, rv0)
         rv0 = -(chia + rv0)
         grc = SUM(rv0)/REAL(natom, KIND=dp)
         rv0 = rv0 - grc
         res = SQRT(SUM(rv0*rv0))
         IF (iunit > 0) THEN
            WRITE (iunit, '(A)', advance='NO') " EEQ| "
         END IF
         ! PCG
         CALL cp_fm_matvec(mmat, rv0, zv0)
         grc = SUM(zv0)/REAL(natom, KIND=dp)
         zv0 = zv0 - grc
         pv0 = zv0
         DO i = 1, natom
            rzk = DOT_PRODUCT(rv0, zv0)
            CALL apply_potential(eeq_mat, ewald_env, ewald_pw, cell, particle_set, pv0, av0)
            grc = SUM(av0)/REAL(natom, KIND=dp)
            av0 = av0 - grc
            papk = DOT_PRODUCT(pv0, av0)
            alk = rzk/papk
            xv0 = xv0 + alk*pv0
            rv0 = rv0 - alk*av0
            rzkm = DOT_PRODUCT(rv0, zv0)
            res = SQRT(DOT_PRODUCT(rv0, rv0))
            !
            CALL apply_potential(eeq_mat, ewald_env, ewald_pw, cell, particle_set, xv0, zv0)
            IF (iunit > 0) THEN
               rdl = MAX(-LOG10(res), 0.0_dp)
               WRITE (iunit, '(I1)', advance='NO') MIN(NINT(rdl), 9)
               IF (MOD(i, 75) == 0) THEN
                  WRITE (iunit, *)
                  WRITE (iunit, '(A)', advance='NO') " EEQ| "
               END IF
            END IF
            IF (res < eps_pcg) EXIT
            CALL cp_fm_matvec(mmat, rv0, zv0)
            grc = SUM(zv0)/REAL(natom, KIND=dp)
            zv0 = zv0 - grc
            rzkp = DOT_PRODUCT(rv0, zv0)
            betak = rzkp/rzk - rzkm/rzk
            pv0 = zv0 + betak*pv0
         END DO
         !
         charges = xv0
         CALL apply_potential(eeq_mat, ewald_env, ewald_pw, cell, particle_set, xv0, rv0)
         rv0 = -(chia + rv0)
         lambda = SUM(rv0)/REAL(natom, KIND=dp)
         rv0 = rv0 - lambda
         res = SQRT(SUM(rv0*rv0))
         IF (iunit > 0) WRITE (iunit, *)
         IF (res > eps_pcg) THEN
            IF (iunit > 0) THEN
               WRITE (iunit, '(A,T61,E20.6)') " EEQ| WARNING PCG did not converge", res
            END IF
         END IF
         !
         DEALLOCATE (av0, xv0, rv0, pv0, zv0)
         !
      ELSE
         ns = natom + 1
         CALL cp_fm_struct_create(mat_struct, context=blacs_env, para_env=para_env, &
                                  nrow_global=ns, ncol_global=ns)
         CALL cp_fm_create(eeq_mat, mat_struct)
         CALL cp_fm_get_info(eeq_mat, nrow_local=nrloc, ncol_local=ncloc, &
                             row_indices=rind, col_indices=cind)
         CALL cp_fm_struct_create(vec_struct, context=blacs_env, para_env=para_env, &
                                  nrow_global=ns, ncol_global=1)
         CALL cp_fm_create(rhs_vec, vec_struct)
         CALL cp_fm_get_info(rhs_vec, nrow_local=nrvloc, ncol_local=ncvloc, &
                             row_indices=rvind, col_indices=cvind)
         !
         ! set up matrix
         CALL cp_fm_set_all(eeq_mat, 1.0_dp, 0.0_dp)
         CALL cp_fm_set_all(rhs_vec, 0.0_dp)
         DO ir = 1, nrloc
            iar = rind(ir)
            ikind = kind_of(iar)
            IF (iar > natom) CYCLE
            ri(1:3) = particle_set(iar)%r(1:3)
            DO ic = 1, ncloc
               iac = cind(ic)
               jkind = kind_of(iac)
               IF (iac > natom) CYCLE
               rj(1:3) = particle_set(iac)%r(1:3)
               IF (iar == iac) THEN
                  grc = gam(ikind) + 2.0_dp*gab(ikind, ikind)*oorootpi
               ELSE
                  rij(1:3) = ri(1:3) - rj(1:3)
                  rij = pbc(rij, cell)
                  dr = SQRT(SUM(rij**2))
                  grc = erf(gab(ikind, jkind)*dr)/dr
               END IF
               eeq_mat%local_data(ir, ic) = grc
            END DO
         END DO
         ! set up rhs vector
         DO ir = 1, nrvloc
            iar = rvind(ir)
            DO ic = 1, ncvloc
               iac = cvind(ic)
               ia = MAX(iar, iac)
               IF (ia > natom) THEN
                  xr = qtot
               ELSE
                  xr = -chia(ia)
               END IF
               rhs_vec%local_data(ir, ic) = xr
            END DO
         END DO
         !
         CALL cp_fm_solve(eeq_mat, rhs_vec)
         !
         charges = 0.0_dp
         lambda = 0.0_dp
         DO ir = 1, nrvloc
            iar = rvind(ir)
            DO ic = 1, ncvloc
               iac = cvind(ic)
               ia = MAX(iar, iac)
               IF (ia <= natom) THEN
                  xr = rhs_vec%local_data(ir, ic)
                  charges(ia) = xr
               ELSE
                  lambda = rhs_vec%local_data(ir, ic)
               END IF
            END DO
         END DO
         CALL para_env%sum(lambda)
         CALL para_env%sum(charges)
         !
      END IF
      !
      ! energy:   0.5*(q^T.X - lambda*totalcharge)
      eeq_energy = 0.5*SUM(charges(1:natom)*chia(1:natom)) - 0.5_dp*lambda*qtot

      IF (do_ewald) THEN
         CALL cp_fm_struct_release(mat_struct)
         CALL cp_fm_release(eeq_mat)
         CALL cp_fm_release(pmat)
         CALL cp_fm_release(mmat)
      ELSE
         CALL cp_fm_struct_release(mat_struct)
         CALL cp_fm_struct_release(vec_struct)
         CALL cp_fm_release(eeq_mat)
         CALL cp_fm_release(rhs_vec)
      END IF

      CALL timestop(handle)

   END SUBROUTINE eeq_solver

! **************************************************************************************************
!> \brief ...
!> \param fm_mat ...
!> \param ewald_env ...
!> \param ewald_pw ...
!> \param cell ...
!> \param particle_set ...
!> \param charges ...
!> \param potential ...
! **************************************************************************************************
   SUBROUTINE apply_potential(fm_mat, ewald_env, ewald_pw, cell, particle_set, charges, potential)
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), INTENT(IN)      :: particle_set
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: charges
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: potential

      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL ewald_env_get(ewald_env, para_env=para_env)
      potential = 0.0_dp
      CALL spme_potential(ewald_env, ewald_pw, cell, particle_set, charges, particle_set, potential)
      CALL para_env%sum(potential)
      CALL cp_fm_matvec(fm_mat, charges, potential, alpha=1.0_dp, beta=1.0_dp)

   END SUBROUTINE apply_potential

END MODULE eeq_method
