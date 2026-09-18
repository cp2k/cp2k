!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Set disperson types for DFT calculations
!> \author JGH (04.2014)
! **************************************************************************************************
MODULE qs_gcp_utils

   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cp_parser_methods,               ONLY: parser_get_next_line,&
                                              parser_get_object
   USE cp_parser_types,                 ONLY: cp_parser_type,&
                                              parser_create,&
                                              parser_release
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE mathconstants,                   ONLY: pi
   USE message_passing,                 ONLY: mp_para_env_type
   USE periodic_table,                  ONLY: get_ptable_info
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_gcp_types,                    ONLY: qs_gcp_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE sto_ng,                          ONLY: get_sto_ng
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_gcp_utils'

   INTEGER, DIMENSION(106) :: nshell = [ &
                              1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, & ! 1-30
                              4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, & ! 31-60
                              6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, & ! 61-90
                              7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7] ! 91-106

   INTEGER, DIMENSION(106) :: nll = [ &
                              1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, & ! 1-30
                              2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, & ! 31-60
                              4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, & ! 61-90
                              4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 0, 0] ! 91-106

   !
   ! Slater exponents for valence states
   ! Aleksander Herman
   ! Empirically adjusted and consistent set of EHT valence orbital parameters for all elements of the periodic table
   ! Modelling and Simulation in Materials Science and Engineering, 12, 21-32 (2004)
   !
   ! Hydrogen uses 1.2000, not the original 1.0000
   !
   REAL(KIND=dp), DIMENSION(4, 106), PARAMETER, PRIVATE :: sexp = RESHAPE([ &
                                 1.2000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 1.6469_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, & !   1   2
                                 0.6534_dp, 0.5305_dp, 0.0000_dp, 0.0000_dp, 1.0365_dp, 0.8994_dp, 0.0000_dp, 0.0000_dp, & !   3   4
                                 1.3990_dp, 1.2685_dp, 0.0000_dp, 0.0000_dp, 1.7210_dp, 1.6105_dp, 0.0000_dp, 0.0000_dp, & !   5   6
                                 2.0348_dp, 1.9398_dp, 0.0000_dp, 0.0000_dp, 2.2399_dp, 2.0477_dp, 0.0000_dp, 0.0000_dp, & !   7   8
                                 2.5644_dp, 2.4022_dp, 0.0000_dp, 0.0000_dp, 2.8812_dp, 2.7421_dp, 0.0000_dp, 0.0000_dp, & !   9  10
                                 0.8675_dp, 0.6148_dp, 0.0000_dp, 0.0000_dp, 1.1935_dp, 0.8809_dp, 0.0000_dp, 0.0000_dp, & !  11  12
                                 1.5143_dp, 1.1660_dp, 0.0000_dp, 0.0000_dp, 1.7580_dp, 1.4337_dp, 0.0000_dp, 0.0000_dp, & !  13  14
                                 1.9860_dp, 1.6755_dp, 0.0000_dp, 0.0000_dp, 2.1362_dp, 1.7721_dp, 0.0000_dp, 0.0000_dp, & !  15  16
                                 2.3617_dp, 2.0176_dp, 0.0000_dp, 0.0000_dp, 2.5796_dp, 2.2501_dp, 0.0000_dp, 0.0000_dp, & !  17  18
                                 0.9362_dp, 0.6914_dp, 0.0000_dp, 0.0000_dp, 1.2112_dp, 0.9329_dp, 0.0000_dp, 0.0000_dp, & !  19  20
                                 1.2870_dp, 0.9828_dp, 2.4341_dp, 0.0000_dp, 1.3416_dp, 1.0104_dp, 2.6439_dp, 0.0000_dp, & !  21  22
                                 1.3570_dp, 0.9947_dp, 2.7809_dp, 0.0000_dp, 1.3804_dp, 0.9784_dp, 2.9775_dp, 0.0000_dp, & !  23  24
                                 1.4761_dp, 1.0641_dp, 3.2208_dp, 0.0000_dp, 1.5465_dp, 1.1114_dp, 3.4537_dp, 0.0000_dp, & !  25  26
                                 1.5650_dp, 1.1001_dp, 3.6023_dp, 0.0000_dp, 1.5532_dp, 1.0594_dp, 3.7017_dp, 0.0000_dp, & !  27  28
                                 1.5791_dp, 1.0527_dp, 3.8962_dp, 0.0000_dp, 1.7778_dp, 1.2448_dp, 0.0000_dp, 0.0000_dp, & !  29  30
                                 2.0675_dp, 1.5073_dp, 0.0000_dp, 0.0000_dp, 2.2702_dp, 1.7680_dp, 0.0000_dp, 0.0000_dp, & !  31  32
                                 2.4546_dp, 1.9819_dp, 0.0000_dp, 0.0000_dp, 2.5680_dp, 2.0548_dp, 0.0000_dp, 0.0000_dp, & !  33  34
                                 2.7523_dp, 2.2652_dp, 0.0000_dp, 0.0000_dp, 2.9299_dp, 2.4617_dp, 0.0000_dp, 0.0000_dp, & !  35  36
                                 1.0963_dp, 0.7990_dp, 0.0000_dp, 0.0000_dp, 1.3664_dp, 1.0415_dp, 0.0000_dp, 0.0000_dp, & !  37  38
                                 1.4613_dp, 1.1100_dp, 2.1576_dp, 0.0000_dp, 1.5393_dp, 1.1647_dp, 2.3831_dp, 0.0000_dp, & !  39  40
                                 1.5926_dp, 1.1738_dp, 2.6256_dp, 0.0000_dp, 1.6579_dp, 1.2186_dp, 2.8241_dp, 0.0000_dp, & !  41  42
                                 1.6930_dp, 1.2490_dp, 2.9340_dp, 0.0000_dp, 1.7347_dp, 1.2514_dp, 3.1524_dp, 0.0000_dp, & !  43  44
                                 1.7671_dp, 1.2623_dp, 3.3113_dp, 0.0000_dp, 1.6261_dp, 1.1221_dp, 3.0858_dp, 0.0000_dp, & !  45  46
                                 1.8184_dp, 1.2719_dp, 3.6171_dp, 0.0000_dp, 1.9900_dp, 1.4596_dp, 0.0000_dp, 0.0000_dp, & !  47  48
                                 2.4649_dp, 1.6848_dp, 0.0000_dp, 0.0000_dp, 2.4041_dp, 1.9128_dp, 0.0000_dp, 0.0000_dp, & !  49  50
                                 2.5492_dp, 2.0781_dp, 0.0000_dp, 0.0000_dp, 2.6576_dp, 2.1718_dp, 0.0000_dp, 0.0000_dp, & !  51  52
                                 2.8080_dp, 2.3390_dp, 0.0000_dp, 0.0000_dp, 2.9595_dp, 2.5074_dp, 0.0000_dp, 0.0000_dp, & !  53  54
                                 1.1993_dp, 0.8918_dp, 0.0000_dp, 0.0000_dp, 1.4519_dp, 1.1397_dp, 0.0000_dp, 0.0000_dp, & !  55  56
                                 1.5331_dp, 1.1979_dp, 2.2743_dp, 4.4161_dp, 1.5379_dp, 1.1930_dp, 2.2912_dp, 4.9478_dp, & !  57  58
                                 1.5162_dp, 1.1834_dp, 2.0558_dp, 4.8982_dp, 1.5322_dp, 1.1923_dp, 2.0718_dp, 5.0744_dp, & !  59  60
                                 1.5486_dp, 1.2018_dp, 2.0863_dp, 5.2466_dp, 1.5653_dp, 1.2118_dp, 2.0999_dp, 5.4145_dp, & !  61  62
                                 1.5762_dp, 1.2152_dp, 2.0980_dp, 5.5679_dp, 1.6703_dp, 1.2874_dp, 2.4862_dp, 5.9888_dp, & !  63  64
                                 1.6186_dp, 1.2460_dp, 2.1383_dp, 5.9040_dp, 1.6358_dp, 1.2570_dp, 2.1472_dp, 6.0598_dp, & !  65  66
                                 1.6536_dp, 1.2687_dp, 2.1566_dp, 6.2155_dp, 1.6723_dp, 1.2813_dp, 2.1668_dp, 6.3703_dp, & !  67  68
                                 1.6898_dp, 1.2928_dp, 2.1731_dp, 6.5208_dp, 1.7063_dp, 1.3030_dp, 2.1754_dp, 6.6686_dp, & !  69  70
                                 1.6647_dp, 1.2167_dp, 2.3795_dp, 0.0000_dp, 1.8411_dp, 1.3822_dp, 2.7702_dp, 0.0000_dp, & !  71  72
                                 1.9554_dp, 1.4857_dp, 3.0193_dp, 0.0000_dp, 2.0190_dp, 1.5296_dp, 3.1936_dp, 0.0000_dp, & !  73  74
                                 2.0447_dp, 1.5276_dp, 3.3237_dp, 0.0000_dp, 2.1361_dp, 1.6102_dp, 3.5241_dp, 0.0000_dp, & !  75  76
                                 2.2167_dp, 1.6814_dp, 3.7077_dp, 0.0000_dp, 2.2646_dp, 1.6759_dp, 3.8996_dp, 0.0000_dp, & !  77  78
                                 2.3185_dp, 1.7126_dp, 4.0525_dp, 0.0000_dp, 2.4306_dp, 1.8672_dp, 0.0000_dp, 0.0000_dp, & !  79  80
                                 2.5779_dp, 1.9899_dp, 0.0000_dp, 0.0000_dp, 2.7241_dp, 2.1837_dp, 0.0000_dp, 0.0000_dp, & !  81  82
                                 2.7869_dp, 2.2146_dp, 0.0000_dp, 0.0000_dp, 2.9312_dp, 2.3830_dp, 0.0000_dp, 0.0000_dp, & !  83  84
                                 3.1160_dp, 2.6200_dp, 0.0000_dp, 0.0000_dp, 3.2053_dp, 2.6866_dp, 0.0000_dp, 0.0000_dp, & !  85  86
                                 1.4160_dp, 1.0598_dp, 0.0000_dp, 0.0000_dp, 1.6336_dp, 1.3011_dp, 0.0000_dp, 0.0000_dp, & !  87  88
                                 1.6540_dp, 1.2890_dp, 2.3740_dp, 3.7960_dp, 1.8381_dp, 1.4726_dp, 2.6584_dp, 4.3613_dp, & !  89  90
                                 1.7770_dp, 1.4120_dp, 2.5710_dp, 4.5540_dp, 1.8246_dp, 1.4588_dp, 2.6496_dp, 4.7702_dp, & !  91  92
                                 1.8451_dp, 1.4739_dp, 2.6940_dp, 4.9412_dp, 1.7983_dp, 1.4366_dp, 2.5123_dp, 4.9882_dp, & !  93  94
                                 1.8011_dp, 1.4317_dp, 2.5170_dp, 5.1301_dp, 1.8408_dp, 1.4418_dp, 2.7349_dp, 5.3476_dp, & !  95  96
                                 1.8464_dp, 1.4697_dp, 2.5922_dp, 5.4596_dp, 1.8647_dp, 1.4838_dp, 2.6205_dp, 5.6140_dp, & !  97  98
                                 1.8890_dp, 1.5050_dp, 2.6590_dp, 5.7740_dp, 1.9070_dp, 1.5190_dp, 2.6850_dp, 5.9220_dp, & !  99 100
                                 1.9240_dp, 1.5320_dp, 2.7090_dp, 6.0690_dp, 1.9400_dp, 1.5440_dp, 2.7300_dp, 6.2130_dp, & ! 101 102
                                 2.1300_dp, 1.7200_dp, 2.9900_dp, 0.0000_dp, 1.9200_dp, 1.4500_dp, 2.9700_dp, 0.0000_dp, & ! 103 104
                                0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp], & ! 105 106
                                                                          [4, 106])

   PUBLIC :: qs_gcp_env_set, qs_gcp_init

! **************************************************************************************************
CONTAINS
! **************************************************************************************************
!> \brief ...
!> \param gcp_env ...
!> \param xc_section ...
! **************************************************************************************************
   SUBROUTINE qs_gcp_env_set(gcp_env, xc_section)
      TYPE(qs_gcp_type), POINTER                         :: gcp_env
      TYPE(section_vals_type), POINTER                   :: xc_section

      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), POINTER                           :: tmpstringlist
      INTEGER                                            :: i_rep, n_rep
      LOGICAL                                            :: explicit
      REAL(dp), POINTER                                  :: params(:)
      TYPE(section_vals_type), POINTER                   :: gcp_section

      CPASSERT(ASSOCIATED(gcp_env))

      gcp_section => section_vals_get_subs_vals(xc_section, "GCP_POTENTIAL")
      CALL section_vals_get(gcp_section, explicit=explicit)
      IF (explicit) THEN
         CALL section_vals_val_get(gcp_section, "VERBOSE", l_val=gcp_env%verbose)
         gcp_env%do_gcp = .TRUE.
         CALL section_vals_val_get(gcp_section, "PARAMETER_FILE_NAME", &
                                   c_val=gcp_env%parameter_file_name)
         CALL section_vals_val_get(gcp_section, "GLOBAL_PARAMETERS", r_vals=params)
         gcp_env%sigma = params(1)
         gcp_env%alpha = params(2)
         gcp_env%beta = params(3)
         gcp_env%eta = params(4)
         ! eamiss definitions
         CALL section_vals_val_get(gcp_section, "DELTA_ENERGY", n_rep_val=n_rep)
         IF (n_rep > 0) THEN
            ALLOCATE (gcp_env%kind_type(n_rep))
            ALLOCATE (gcp_env%ea(n_rep))
            DO i_rep = 1, n_rep
               CALL section_vals_val_get(gcp_section, "DELTA_ENERGY", i_rep_val=i_rep, &
                                         c_vals=tmpstringlist)
               READ (tmpstringlist(1), *) gcp_env%kind_type(i_rep)
               READ (tmpstringlist(2), *) gcp_env%ea(i_rep)
            END DO
         END IF
      ELSE
         gcp_env%do_gcp = .FALSE.
      END IF

   END SUBROUTINE qs_gcp_env_set

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param gcp_env ...
! **************************************************************************************************
   SUBROUTINE qs_gcp_init(qs_env, gcp_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_gcp_type), POINTER                         :: gcp_env

      REAL(KIND=dp), PARAMETER                           :: epsc = 1.e-6_dp

      CHARACTER(LEN=10)                                  :: aname
      CHARACTER(LEN=2)                                   :: element_symbol
      INTEGER                                            :: i, ikind, nbas, nel, nkind, nsto, za
      LOGICAL                                            :: at_end
      REAL(KIND=dp)                                      :: ea
      REAL(KIND=dp), DIMENSION(10)                       :: al, cl
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      IF (gcp_env%do_gcp) THEN
         CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, nkind=nkind)
         ALLOCATE (gcp_env%gcp_kind(nkind))
         DO ikind = 1, nkind
            qs_kind => qs_kind_set(ikind)
            gcp_env%gcp_kind(ikind)%rcsto = 0.0_dp
            CALL get_qs_kind(qs_kind, element_symbol=element_symbol)
            CALL get_ptable_info(element_symbol, number=za)
            gcp_env%gcp_kind(ikind)%za = za
            gcp_env%gcp_kind(ikind)%asto = gcp_env%eta*SUM(sexp(1:4, za))/REAL(nll(za), KIND=dp)
            gcp_env%gcp_kind(ikind)%nq = nshell(za)
            gcp_env%gcp_kind(ikind)%rcsto = ((nshell(za) - 1)*2.5_dp - LOG(epsc))/gcp_env%gcp_kind(ikind)%asto
            ! basis
            NULLIFY (orb_basis)
            CALL get_qs_kind(qs_kind, basis_set=orb_basis, basis_type="ORB")
            CALL get_gto_basis_set(gto_basis_set=orb_basis, nsgf=nbas)
            nel = SUM(qs_kind%elec_conf)
            gcp_env%gcp_kind(ikind)%nbvirt = REAL(nbas, KIND=dp) - 0.5_dp*REAL(nel, KIND=dp)
            ! STO-nG
            nsto = SIZE(gcp_env%gcp_kind(ikind)%al)
            CALL get_sto_ng(gcp_env%gcp_kind(ikind)%asto, nsto, nshell(za), 0, al, cl)
            DO i = 1, nsto
               gcp_env%gcp_kind(ikind)%al(i) = al(i)
               gcp_env%gcp_kind(ikind)%cl(i) = cl(i)*(2._dp*al(i)/pi)**0.75_dp
            END DO
         END DO
         ! eamiss from data file
         IF (gcp_env%parameter_file_name /= "---") THEN
            BLOCK
               TYPE(cp_parser_type) :: parser
               CALL get_qs_env(qs_env, para_env=para_env)
               DO ikind = 1, nkind
                  qs_kind => qs_kind_set(ikind)
                  CALL get_qs_kind(qs_kind, element_symbol=element_symbol)
                  CALL get_ptable_info(element_symbol, number=za)
                  !
                  CALL parser_create(parser, gcp_env%parameter_file_name, para_env=para_env)
                  ea = 0.0_dp
                  DO
                     at_end = .FALSE.
                     CALL parser_get_next_line(parser, 1, at_end)
                     IF (at_end) EXIT
                     CALL parser_get_object(parser, aname)
                     IF (TRIM(aname) == element_symbol) THEN
                        CALL parser_get_object(parser, ea)
                        EXIT
                     END IF
                  END DO
                  CALL parser_release(parser)
                  gcp_env%gcp_kind(ikind)%eamiss = ea
               END DO
            END BLOCK
         END IF
         !
         ! eamiss from input
         IF (ASSOCIATED(gcp_env%kind_type)) THEN
            DO i = 1, SIZE(gcp_env%kind_type)
               IF (TRIM(gcp_env%kind_type(i)) == "XX") CYCLE
               element_symbol = TRIM(gcp_env%kind_type(i))
               CALL get_ptable_info(element_symbol, number=za)
               ea = gcp_env%ea(i)
               DO ikind = 1, nkind
                  IF (za == gcp_env%gcp_kind(ikind)%za) THEN
                     gcp_env%gcp_kind(ikind)%eamiss = ea
                  END IF
               END DO
            END DO
         END IF
      END IF

   END SUBROUTINE qs_gcp_init
! **************************************************************************************************

END MODULE qs_gcp_utils

