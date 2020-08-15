!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Simplified Tamm Dancoff approach (sTDA).
! **************************************************************************************************
MODULE qs_tddfpt2_stda_types

   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cp_control_types,                ONLY: stda_control_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE kinds,                           ONLY: dp
   USE physcon,                         ONLY: evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   TYPE stda_kind_type
      CHARACTER(LEN=2)                            :: symbol !element symbol
      INTEGER                                     :: z !atomic_number
      INTEGER                                     :: kind_number !kind number
      REAL(KIND=dp)                               :: hardness_param !hardness_parameter eta
      REAL(KIND=dp)                               :: rcut !cutoff radius for short-range Coulomb
   END TYPE

   TYPE stda_kind_p_type
      TYPE(stda_kind_type), POINTER               :: kind_param
   END TYPE

   TYPE stda_env_type
      !fraction of non-local HF exchange !! REDUNDANT also in stda_control
      REAL(KIND=dp)                                :: hfx_fraction
      ! empirical parameters
      REAL(KIND=dp)                                :: alpha_param
      REAL(KIND=dp)                                :: beta_param
      ! filter for TD matrix
      REAL(KIND=dp)                                :: eps_td_filter
      TYPE(stda_kind_p_type), DIMENSION(:), POINTER:: kind_param_set
      !number of atomic orbitals, number of occupied orbitals
      INTEGER, DIMENSION(2)                        :: n_ao
      INTEGER, DIMENSION(2)                        :: nactive
   END TYPE

   !PARAMETERS
!&<
   INTEGER, PARAMETER, PRIVATE :: nelem = 103
   !   H                                                                      He
   !   Li Be                                                 B  C  N  O  F    Ne
   !   Na Mg                                                 Al Si P  S  Cl   Ar
   !   K  Ca Sc                Ti V  Cr Mn Fe Co Ni Cu Zn    Ga Ge As Se Br   Kr
   !   Rb Sr Y                 Zr Nb Mo Tc Ru Rh Pd Ag Cd    In Sn Sb Te I    Xe
   !   Cs Ba La Ce-Lu          Hf Ta W  Re Os Ir Pt Au Hg    Tl Pb Bi Po At   Rn
   !   Fr Ra Ac Th Pa U        Np Pu Am Cm Bk Cf Es Fm Md    No Lr  103

   REAL(KIND=dp), DIMENSION(1:nelem), PARAMETER, PRIVATE:: hardness = &
         (/6.4299544220_dp, 12.544911890_dp, & ! 2 H-He
           2.3745866560_dp, 3.4967633530_dp, 4.6190089720_dp, 5.7409789220_dp, &
           6.8624665290_dp, 7.9854357010_dp, 9.1064753720_dp, 10.23034050_dp, & ! 8 Li-Ne
           2.444141360_dp, 3.0146513830_dp, 3.5849070740_dp, 4.15513090_dp, &
           4.7258039740_dp, 5.2959792410_dp, 5.8661864840_dp, 6.4366187140_dp, & ! 8 Na-Ar
           2.3273178360_dp, 2.7587238140_dp, 2.8581921140_dp, 2.9578300430_dp, &
           3.0573410060_dp, 3.1567254290_dp, 3.2563827230_dp, 3.3559314050_dp, &
           3.4556091170_dp, 3.5550133130_dp, 3.6544183480_dp, 3.7541601450_dp, &
           4.1855197930_dp, 4.6166272460_dp, 5.0662145070_dp, 5.4794960970_dp, &
           5.9110996450_dp, 6.3418467680_dp, & ! 18 K-Kr
           2.1204582570_dp, 2.5373700480_dp, 2.6335468980_dp, 2.7297528930_dp, &
           2.8259738860_dp, 2.9221296040_dp, 3.0183708780_dp, 3.1145981770_dp, &
           3.210756280_dp, 3.3069474480_dp, 3.4031948570_dp, 3.4993761390_dp, &
           3.9163692460_dp, 4.3332332190_dp, 4.7500787860_dp, 5.1669793270_dp, &
           5.5838871020_dp, 6.000897330_dp, & ! 18 Rb-Xe
           0.6829150240_dp, 0.9200946840_dp, 1.1570887860_dp, 1.39427570_dp, &
           1.6314731730_dp, 1.8684389980_dp, 2.1056577930_dp, 2.3426646420_dp, &
           2.5798149820_dp, 2.8170264230_dp, 3.0540365330_dp, 3.2911692310_dp, &
           3.5282971610_dp, 3.7655249290_dp, 4.0025547030_dp, 4.2394783410_dp, &
           4.4765830210_dp, 4.7065224490_dp, 4.9508466940_dp, 5.1879311720_dp, &
           5.4256076210_dp, 5.6619144310_dp, 5.900042920_dp, 6.1367145320_dp, &
           6.3741299770_dp, 6.6102656130_dp, 1.7043485810_dp, 1.9413526120_dp, &
           2.178491510_dp, 2.4158121060_dp, 2.6527780840_dp, 2.8899554570_dp, & ! 32 Cs-Rn
           0.9882529880_dp, 1.2819499970_dp, 1.3497250380_dp, 1.4175257380_dp, &
           1.9368567520_dp, 2.2305576050_dp, 2.5241204960_dp, 3.0436128480_dp, &
           3.4168675260_dp, 3.4049844440_dp, 3.9244199680_dp, 4.2180813280_dp, &
           4.5115926320_dp, 4.8050928950_dp, 5.0989816210_dp, 5.3926054620_dp, &
           5.4606987930_dp/) ! 17 Fr-Lr
!&>

   REAL(KIND=dp), DIMENSION(2), PARAMETER, PRIVATE:: alpha = (/1.420_dp, 0.480_dp/)
   REAL(KIND=dp), DIMENSION(2), PARAMETER, PRIVATE:: beta = (/0.200_dp, 1.830_dp/)

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_stda_types'

   PUBLIC:: stda_env_type, &
            allocate_stda_env, deallocate_stda_env, stda_init_param

CONTAINS

! **************************************************************************************************
!> \brief Get the parameters needed for an sTDA calculation
!> \param qs_env ...
!> \param stda_kernel ...
!> \param stda_control ...
! **************************************************************************************************
   SUBROUTINE stda_init_param(qs_env, stda_kernel, stda_control)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(stda_env_type)                                :: stda_kernel
      TYPE(stda_control_type)                            :: stda_control

      INTEGER                                            :: ikind, iounit, nkind
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(gto_basis_set_type), POINTER                  :: tmp_basis_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(stda_kind_type), POINTER                      :: kind_param

      NULLIFY (logger)
      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      CPASSERT(ASSOCIATED(stda_kernel%kind_param_set))

      NULLIFY (atomic_kind_set, qs_kind_set) !get element symbol and atomic number
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set)
      nkind = SIZE(qs_kind_set)

      DO ikind = 1, nkind
         atomic_kind => atomic_kind_set(ikind)
         kind_param => stda_kernel%kind_param_set(ikind)%kind_param
         ! element symbol, kind_number, atomic number
         CALL get_atomic_kind(atomic_kind, &
                              element_symbol=stda_kernel%kind_param_set(ikind)%kind_param%symbol, &
                              kind_number=stda_kernel%kind_param_set(ikind)%kind_param%kind_number, &
                              z=stda_kernel%kind_param_set(ikind)%kind_param%z)
      ENDDO

      ! HFx_fraction
      stda_kernel%hfx_fraction = stda_control%hfx_fraction
      ! alpha and beta parameter
      stda_kernel%alpha_param = alpha(1) + stda_kernel%hfx_fraction*alpha(2)
      stda_kernel%beta_param = beta(1) + stda_kernel%hfx_fraction*beta(2)

      ! TD Filter
      stda_kernel%eps_td_filter = stda_control%eps_td_filter

      DO ikind = 1, nkind
         ! hardness parameter
         stda_kernel%kind_param_set(ikind)%kind_param%hardness_param = &
            hardness(stda_kernel%kind_param_set(ikind)%kind_param%z)*2.0_dp/evolt
         ! rcut parameter
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind, basis_set=tmp_basis_set)
         CALL get_gto_basis_set(tmp_basis_set, &
                                kind_radius=stda_kernel%kind_param_set(ikind)%kind_param%rcut)
      ENDDO

   END SUBROUTINE stda_init_param
! **************************************************************************************************
!> \brief Allocate the sTDA environment
!> \param qs_env ...
!> \param stda_kernel ...
!> \param n_ao ...
!> \param nactive ...
! **************************************************************************************************
   SUBROUTINE allocate_stda_env(qs_env, stda_kernel, n_ao, nactive)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(stda_env_type)                                :: stda_kernel
      INTEGER, INTENT(IN)                                :: n_ao
      INTEGER, DIMENSION(:), INTENT(IN)                  :: nactive

      INTEGER                                            :: ii, nkind
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(stda_kind_type), POINTER                      :: kind_param

      stda_kernel%hfx_fraction = 0.0d0
      stda_kernel%alpha_param = 0.0d0
      stda_kernel%beta_param = 0.0d0
      stda_kernel%nactive = 0
      stda_kernel%nactive(1:2) = nactive(1:2)
      stda_kernel%n_ao = n_ao
      NULLIFY (stda_kernel%kind_param_set)

      ! initialize stda_kind_parameters
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      nkind = SIZE(qs_kind_set)

      ALLOCATE (stda_kernel%kind_param_set(nkind))
      DO ii = 1, nkind
         NULLIFY (kind_param)
         CALL allocate_stda_kind_param(kind_param)
         stda_kernel%kind_param_set(ii)%kind_param => kind_param
      ENDDO

   END SUBROUTINE allocate_stda_env
! **************************************************************************************************
!> \brief Deallocate the sTDA environment
!> \param stda_kernel ...
! **************************************************************************************************
   SUBROUTINE deallocate_stda_env(stda_kernel)

      TYPE(stda_env_type)                                :: stda_kernel

      INTEGER                                            :: ii
      TYPE(stda_kind_type), POINTER                      :: kind_param

      ! deallocate stda_kind_parameters
      IF (ASSOCIATED(stda_kernel%kind_param_set)) THEN
         DO ii = 1, SIZE(stda_kernel%kind_param_set)
            kind_param => stda_kernel%kind_param_set(ii)%kind_param
            CALL deallocate_stda_kind_param(kind_param)
         ENDDO
         DEALLOCATE (stda_kernel%kind_param_set)
         NULLIFY (stda_kernel%kind_param_set)
      ENDIF

   END SUBROUTINE deallocate_stda_env
! **************************************************************************************************
!> \brief Allocate sTDA kind parameter
!> \param kind_param ...
! **************************************************************************************************
   SUBROUTINE allocate_stda_kind_param(kind_param)

      TYPE(stda_kind_type), POINTER                      :: kind_param

      IF (ASSOCIATED(kind_param)) &
         CALL deallocate_stda_kind_param(kind_param)

      ALLOCATE (kind_param)

      kind_param%symbol = ""
      kind_param%z = 0
      kind_param%kind_number = 0
      kind_param%hardness_param = 0.0_dp
      kind_param%rcut = 0.0_dp

   END SUBROUTINE allocate_stda_kind_param
! **************************************************************************************************
!> \brief Deallocate sTDA kind parameter
!> \param kind_param ...
! **************************************************************************************************
   SUBROUTINE deallocate_stda_kind_param(kind_param)

      TYPE(stda_kind_type), POINTER                      :: kind_param

      CPASSERT(ASSOCIATED(kind_param))
      DEALLOCATE (kind_param)
      NULLIFY (kind_param)

   END SUBROUTINE deallocate_stda_kind_param

END MODULE qs_tddfpt2_stda_types
