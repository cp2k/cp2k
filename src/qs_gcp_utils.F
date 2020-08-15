!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Set disperson types for DFT calculations
!> \author JGH (04.2014)
! **************************************************************************************************
MODULE qs_gcp_utils

   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cp_para_types,                   ONLY: cp_para_env_type
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

   INTEGER, DIMENSION(106) :: nshell = (/ &
                              1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, & ! 1-30
                              4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, & ! 31-60
                              6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, & ! 61-90
                              7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7/) ! 91-106

   INTEGER, DIMENSION(106) :: nll = (/ &
                              1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, & ! 1-30
                              2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, & ! 31-60
                              4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, & ! 61-90
                              4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 0, 0/) ! 91-106

   !
   ! Slater exponents for valence states
   ! Aleksander Herman
   ! Empirically adjusted and consistent set of EHT valence orbital parameters for all elements of the periodic table
   ! Modelling and Simulation in Materials Science and Engineering, 12, 21-32 (2004)
   !
   ! Hydrogen uses 1.2000, not the original 1.0000
   !
   REAL(KIND=dp), DIMENSION(4, 106) :: sexp = RESHAPE((/ &
                           1.2000, 0.0000, 0.0000, 0.0000, 1.6469, 0.0000, 0.0000, 0.0000, 0.6534, 0.5305, 0.0000, 0.0000, & ! 1 2 3
                           1.0365, 0.8994, 0.0000, 0.0000, 1.3990, 1.2685, 0.0000, 0.0000, 1.7210, 1.6105, 0.0000, 0.0000, & ! 4 5 6
                           2.0348, 1.9398, 0.0000, 0.0000, 2.2399, 2.0477, 0.0000, 0.0000, 2.5644, 2.4022, 0.0000, 0.0000, & ! 7 8 9
                        2.8812, 2.7421, 0.0000, 0.0000, 0.8675, 0.6148, 0.0000, 0.0000, 1.1935, 0.8809, 0.0000, 0.0000, & ! 10 11 12
                        1.5143, 1.1660, 0.0000, 0.0000, 1.7580, 1.4337, 0.0000, 0.0000, 1.9860, 1.6755, 0.0000, 0.0000, & ! 13 14 15
                        2.1362, 1.7721, 0.0000, 0.0000, 2.3617, 2.0176, 0.0000, 0.0000, 2.5796, 2.2501, 0.0000, 0.0000, & ! 16 17 18
                        0.9362, 0.6914, 0.0000, 0.0000, 1.2112, 0.9329, 0.0000, 0.0000, 1.2870, 0.9828, 2.4341, 0.0000, & ! 19 20 21
                        1.3416, 1.0104, 2.6439, 0.0000, 1.3570, 0.9947, 2.7809, 0.0000, 1.3804, 0.9784, 2.9775, 0.0000, & ! 22 23 24
                        1.4761, 1.0641, 3.2208, 0.0000, 1.5465, 1.1114, 3.4537, 0.0000, 1.5650, 1.1001, 3.6023, 0.0000, & ! 25 26 27
                        1.5532, 1.0594, 3.7017, 0.0000, 1.5791, 1.0527, 3.8962, 0.0000, 1.7778, 1.2448, 0.0000, 0.0000, & ! 28 29 30
                        2.0675, 1.5073, 0.0000, 0.0000, 2.2702, 1.7680, 0.0000, 0.0000, 2.4546, 1.9819, 0.0000, 0.0000, & ! 31 32 33
                        2.5680, 2.0548, 0.0000, 0.0000, 2.7523, 2.2652, 0.0000, 0.0000, 2.9299, 2.4617, 0.0000, 0.0000, & ! 34 35 36
                        1.0963, 0.7990, 0.0000, 0.0000, 1.3664, 1.0415, 0.0000, 0.0000, 1.4613, 1.1100, 2.1576, 0.0000, & ! 37 38 39
                        1.5393, 1.1647, 2.3831, 0.0000, 1.5926, 1.1738, 2.6256, 0.0000, 1.6579, 1.2186, 2.8241, 0.0000, & ! 40 41 42
                        1.6930, 1.2490, 2.9340, 0.0000, 1.7347, 1.2514, 3.1524, 0.0000, 1.7671, 1.2623, 3.3113, 0.0000, & ! 43 44 45
                        1.6261, 1.1221, 3.0858, 0.0000, 1.8184, 1.2719, 3.6171, 0.0000, 1.9900, 1.4596, 0.0000, 0.0000, & ! 46 47 48
                        2.4649, 1.6848, 0.0000, 0.0000, 2.4041, 1.9128, 0.0000, 0.0000, 2.5492, 2.0781, 0.0000, 0.0000, & ! 49 50 51
                        2.6576, 2.1718, 0.0000, 0.0000, 2.8080, 2.3390, 0.0000, 0.0000, 2.9595, 2.5074, 0.0000, 0.0000, & ! 52 53 54
                        1.1993, 0.8918, 0.0000, 0.0000, 1.4519, 1.1397, 0.0000, 0.0000, 1.5331, 1.1979, 2.2743, 4.4161, & ! 55 56 57
                        1.5379, 1.1930, 2.2912, 4.9478, 1.5162, 1.1834, 2.0558, 4.8982, 1.5322, 1.1923, 2.0718, 5.0744, & ! 58 59 60
                        1.5486, 1.2018, 2.0863, 5.2466, 1.5653, 1.2118, 2.0999, 5.4145, 1.5762, 1.2152, 2.0980, 5.5679, & ! 61 62 63
                        1.6703, 1.2874, 2.4862, 5.9888, 1.6186, 1.2460, 2.1383, 5.9040, 1.6358, 1.2570, 2.1472, 6.0598, & ! 64 65 66
                        1.6536, 1.2687, 2.1566, 6.2155, 1.6723, 1.2813, 2.1668, 6.3703, 1.6898, 1.2928, 2.1731, 6.5208, & ! 67 68 69
                        1.7063, 1.3030, 2.1754, 6.6686, 1.6647, 1.2167, 2.3795, 0.0000, 1.8411, 1.3822, 2.7702, 0.0000, & ! 70 71 72
                        1.9554, 1.4857, 3.0193, 0.0000, 2.0190, 1.5296, 3.1936, 0.0000, 2.0447, 1.5276, 3.3237, 0.0000, & ! 73 74 75
                        2.1361, 1.6102, 3.5241, 0.0000, 2.2167, 1.6814, 3.7077, 0.0000, 2.2646, 1.6759, 3.8996, 0.0000, & ! 76 77 78
                        2.3185, 1.7126, 4.0525, 0.0000, 2.4306, 1.8672, 0.0000, 0.0000, 2.5779, 1.9899, 0.0000, 0.0000, & ! 79 80 81
                        2.7241, 2.1837, 0.0000, 0.0000, 2.7869, 2.2146, 0.0000, 0.0000, 2.9312, 2.3830, 0.0000, 0.0000, & ! 82 83 84
                        3.1160, 2.6200, 0.0000, 0.0000, 3.2053, 2.6866, 0.0000, 0.0000, 1.4160, 1.0598, 0.0000, 0.0000, & ! 85 86 87
                        1.6336, 1.3011, 0.0000, 0.0000, 1.6540, 1.2890, 2.3740, 3.7960, 1.8381, 1.4726, 2.6584, 4.3613, & ! 88 89 90
                        1.7770, 1.4120, 2.5710, 4.5540, 1.8246, 1.4588, 2.6496, 4.7702, 1.8451, 1.4739, 2.6940, 4.9412, & ! 91 92 93
                        1.7983, 1.4366, 2.5123, 4.9882, 1.8011, 1.4317, 2.5170, 5.1301, 1.8408, 1.4418, 2.7349, 5.3476, & ! 94 95 96
                        1.8464, 1.4697, 2.5922, 5.4596, 1.8647, 1.4838, 2.6205, 5.6140, 1.8890, 1.5050, 2.6590, 5.7740, & ! 97 98 99
                     1.9070, 1.5190, 2.6850, 5.9220, 1.9240, 1.5320, 2.7090, 6.0690, 1.9400, 1.5440, 2.7300, 6.2130, & ! 100 101 102
                     2.1300, 1.7200, 2.9900, 0.0000, 1.9200, 1.4500, 2.9700, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, & ! 103 104 105
                                                      0.0000, 0.0000, 0.0000, 0.0000/), & ! 106
                                                      (/4, 106/))

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
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_parser_type), POINTER                      :: parser
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis
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
            CALL get_qs_env(qs_env, para_env=para_env)
            DO ikind = 1, nkind
               qs_kind => qs_kind_set(ikind)
               CALL get_qs_kind(qs_kind, element_symbol=element_symbol)
               CALL get_ptable_info(element_symbol, number=za)
               !
               NULLIFY (parser)
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

