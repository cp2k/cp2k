!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief SOC parameters of Hartwigsten-Goedecker-Hutter pseudopotentials 
!>        (temporary file until storage of SOC parameters in data directory is there
!> \par History
!>      06.2023 created [Jan Wilhelm]
! **************************************************************************************************
MODULE external_potential_soc_parameters
   USE external_potential_types,        ONLY: gth_potential_p_type,&
                                              gth_potential_type,&
                                              sgp_potential_p_type,&
                                              sgp_potential_type, init_potential
  USE kinds,                           ONLY: dp,&
                                              int_8
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'external_potential_soc_parameters'

   PUBLIC :: get_kprj_SOC_parameters

CONTAINS

   SUBROUTINE get_kprj_SOC_parameters(gth_potential, xc_functional,element_symbol)

      TYPE(gth_potential_type), POINTER                  :: gth_potential
      CHARACTER(LEN=*)                                   :: xc_functional
     CHARACTER(LEN=2)           :: element_symbol

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_kprj_SOC_parameters'

      INTEGER :: handle, charge

      CALL timeset(routineN, handle)

      charge = NINT(gth_potential%zeff)

      gth_potential%kprj_ppnl(:,:,:) = 0.0_dp 

      SELECT CASE (xc_functional)
      CASE("LDA")
        CALL get_kprj_SOC_parameters_LDA(gth_potential, element_symbol,charge)
      CASE("PBE")

      CASE DEFAULT
        CPABORT("Unknown xc functional in SOC part of pseudopotential.")
      END SELECT

      ! write SOC parameters from kprj_ppnl into wprj_ppnl
      CALL init_potential(gth_potential)

      CALL timestop(handle)

   END SUBROUTINE get_kprj_SOC_parameters

   SUBROUTINE get_kprj_SOC_parameters_LDA(gth_potential, element_symbol, charge)

     TYPE(gth_potential_type), POINTER                  :: gth_potential

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_kprj_SOC_parameters_LDA'

      INTEGER :: handle, charge
     CHARACTER(LEN=2)           :: element_symbol

      CALL timeset(routineN, handle)

IF(element_symbol == "Ag" .AND. charge == 1) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.01969167_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.00288247_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.00288247_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.00682116_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00945505_dp

ELSE IF(element_symbol == "Ag" .AND. charge == 11) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.20186901_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.13510134_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.13510134_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.31970813_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.04519317_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.02273792_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.02273792_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.05156476_dp

ELSE IF(element_symbol == "Ag" .AND. charge == 19) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.26670728_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.68349949_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.68349949_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.80872750_dp
   gth_potential%kprj_ppnl(1,1,2) = -0.03218230_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.06201267_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.06201267_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.07031576_dp

ELSE IF(element_symbol == "Al" .AND. charge == 3) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.00615404_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.00166789_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.00166789_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.00394694_dp

ELSE IF(element_symbol == "Ar" .AND. charge == 8) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.01639467_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.01232694_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.01232694_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.02917087_dp

ELSE IF(element_symbol == "As" .AND. charge == 5) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.05246616_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.00868894_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.00868894_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.02056179_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00427335_dp

ELSE IF(element_symbol == "At" .AND. charge == 7) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.48077365_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.01477071_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.01477071_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.03495387_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.03754438_dp

ELSE IF(element_symbol == "Au" .AND. charge == 1) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.03934865_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.05618995_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.05618995_dp
   gth_potential%kprj_ppnl(1,3,1) = -0.04561301_dp
   gth_potential%kprj_ppnl(3,1,1) = -0.04561301_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.13296969_dp
   gth_potential%kprj_ppnl(2,3,1) = 0.10794007_dp
   gth_potential%kprj_ppnl(3,2,1) = 0.10794007_dp
   gth_potential%kprj_ppnl(3,3,1) = -0.15342717_dp

ELSE IF(element_symbol == "Au" .AND. charge == 11) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.96005489_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.54790119_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.54790119_dp
   gth_potential%kprj_ppnl(2,2,1) = 1.29657086_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.14848354_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.07473921_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.07473921_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.16949260_dp

ELSE IF(element_symbol == "Au" .AND. charge == 19) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.62807493_dp
   gth_potential%kprj_ppnl(1,2,1) = 1.81668773_dp
   gth_potential%kprj_ppnl(2,1,1) = 1.81668773_dp
   gth_potential%kprj_ppnl(1,3,1) = -1.65141098_dp
   gth_potential%kprj_ppnl(3,1,1) = -1.65141098_dp
   gth_potential%kprj_ppnl(2,2,1) = -4.89888134_dp
   gth_potential%kprj_ppnl(2,3,1) = 3.90795165_dp
   gth_potential%kprj_ppnl(3,2,1) = 3.90795165_dp
   gth_potential%kprj_ppnl(3,3,1) = -2.77740200_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.18880114_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.03175114_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.03175114_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.03600240_dp

ELSE IF(element_symbol == "B " .AND. charge == 3) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.00087786_dp

ELSE IF(element_symbol == "Ba" .AND. charge == 10) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.41508327_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.25344615_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.25344615_dp
   gth_potential%kprj_ppnl(1,3,1) = 0.27561093_dp
   gth_potential%kprj_ppnl(3,1,1) = 0.27561093_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.59976305_dp
   gth_potential%kprj_ppnl(2,3,1) = -0.65221449_dp
   gth_potential%kprj_ppnl(3,2,1) = -0.65221449_dp
   gth_potential%kprj_ppnl(3,3,1) = 0.92706460_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.01766026_dp
   gth_potential%kprj_ppnl(1,1,3) = 0.00615085_dp

ELSE IF(element_symbol == "Ba" .AND. charge == 2) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.04101721_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.03752472_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.03752472_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.08879970_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00532055_dp

ELSE IF(element_symbol == "Be" .AND. charge == 2) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.00012877_dp

ELSE IF(element_symbol == "Bi" .AND. charge == 5) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.30531353_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.00977575_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.00977575_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.02313366_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.02921737_dp

ELSE IF(element_symbol == "Br" .AND. charge == 7) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.07400663_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.02906775_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.02906775_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.06878684_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00714446_dp

ELSE IF(element_symbol == "C " .AND. charge == 4) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.00410365_dp

ELSE IF(element_symbol == "Ca" .AND. charge == 10) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.02170091_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.02394376_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.02394376_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.05666128_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00037106_dp

ELSE IF(element_symbol == "Ca" .AND. charge == 2) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.00336201_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.00540007_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.00540007_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.01277891_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00081449_dp

ELSE IF(element_symbol == "Cd" .AND. charge == 12) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.24140932_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.15923958_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.15923958_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.37682962_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.04932879_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.02594377_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.02594377_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.05883495_dp

ELSE IF(element_symbol == "Cd" .AND. charge == 2) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.03771366_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.00084877_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.00084877_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.00200856_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00310866_dp

ELSE IF(element_symbol == "Ce" .AND. charge == 12) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.55905862_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.30893198_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.30893198_dp
   gth_potential%kprj_ppnl(1,3,1) = 0.22367871_dp
   gth_potential%kprj_ppnl(3,1,1) = 0.22367871_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.73792420_dp
   gth_potential%kprj_ppnl(2,3,1) = -0.52932045_dp
   gth_potential%kprj_ppnl(3,2,1) = -0.52932045_dp
   gth_potential%kprj_ppnl(3,3,1) = 0.37619085_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.01403180_dp
   gth_potential%kprj_ppnl(1,1,3) = 0.00781271_dp

ELSE IF(element_symbol == "Cl" .AND. charge == 7) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.01002014_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.01089588_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.01089588_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.02578435_dp

ELSE IF(element_symbol == "Co" .AND. charge == 17) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.20791470_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.04574353_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.04574353_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.10824895_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.01147529_dp

ELSE IF(element_symbol == "Co" .AND. charge == 9) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.01752129_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.02483306_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.02483306_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.05876575_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00713682_dp

ELSE IF(element_symbol == "Cr" .AND. charge == 14) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.16978070_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.00017378_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.00017378_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.00041123_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00900712_dp

ELSE IF(element_symbol == "Cr" .AND. charge == 6) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.01317649_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.01505450_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.01505450_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.03562545_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00351409_dp

ELSE IF(element_symbol == "Cs" .AND. charge == 1) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.10762747_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.05837136_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.05837136_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.13813185_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00325918_dp

ELSE IF(element_symbol == "Cs" .AND. charge == 9) THEN
   gth_potential%kprj_ppnl(1,1,1) = -2.74852756_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.76700684_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.76700684_dp
   gth_potential%kprj_ppnl(1,3,1) = 0.35818758_dp
   gth_potential%kprj_ppnl(3,1,1) = 0.35818758_dp
   gth_potential%kprj_ppnl(2,2,1) = 1.81506946_dp
   gth_potential%kprj_ppnl(2,3,1) = -0.84762652_dp
   gth_potential%kprj_ppnl(3,2,1) = -0.84762652_dp
   gth_potential%kprj_ppnl(3,3,1) = 1.20482534_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.01084084_dp
   gth_potential%kprj_ppnl(1,1,3) = 0.00475984_dp

ELSE IF(element_symbol == "Cu" .AND. charge == 1) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.01079247_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.00284563_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.00284563_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.00673399_dp
   gth_potential%kprj_ppnl(1,1,2) = -0.00073000_dp

ELSE IF(element_symbol == "Cu" .AND. charge == 11) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.02406730_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.03231927_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.03231927_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.07648136_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.01048914_dp

ELSE IF(element_symbol == "Cu" .AND. charge == 19) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.13410535_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.21399818_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.21399818_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.25320605_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.01450302_dp

ELSE IF(element_symbol == "Dy" .AND. charge == 20) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.12437116_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.56093057_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.56093057_dp
   gth_potential%kprj_ppnl(1,3,1) = -0.32468735_dp
   gth_potential%kprj_ppnl(3,1,1) = -0.32468735_dp
   gth_potential%kprj_ppnl(2,2,1) = 1.32740399_dp
   gth_potential%kprj_ppnl(2,3,1) = 0.76835051_dp
   gth_potential%kprj_ppnl(3,2,1) = 0.76835051_dp
   gth_potential%kprj_ppnl(3,3,1) = -1.09214158_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.05828571_dp
   gth_potential%kprj_ppnl(1,1,3) = 0.01793827_dp

ELSE IF(element_symbol == "Er" .AND. charge == 22) THEN
   gth_potential%kprj_ppnl(1,1,1) = 1.25513474_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.37692191_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.37692191_dp
   gth_potential%kprj_ppnl(1,3,1) = 0.22288650_dp
   gth_potential%kprj_ppnl(3,1,1) = 0.22288650_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.89196004_dp
   gth_potential%kprj_ppnl(2,3,1) = -0.52744572_dp
   gth_potential%kprj_ppnl(3,2,1) = -0.52744572_dp
   gth_potential%kprj_ppnl(3,3,1) = 0.74971696_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.09186175_dp
   gth_potential%kprj_ppnl(1,1,3) = 0.02178337_dp

ELSE IF(element_symbol == "Eu" .AND. charge == 17) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.25225840_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.24692195_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.24692195_dp
   gth_potential%kprj_ppnl(1,3,1) = -0.14287549_dp
   gth_potential%kprj_ppnl(3,1,1) = -0.14287549_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.58432399_dp
   gth_potential%kprj_ppnl(2,3,1) = 0.33810512_dp
   gth_potential%kprj_ppnl(3,2,1) = 0.33810512_dp
   gth_potential%kprj_ppnl(3,3,1) = -0.48058621_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.05102753_dp
   gth_potential%kprj_ppnl(1,1,3) = 0.01326718_dp

ELSE IF(element_symbol == "F " .AND. charge == 7) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.01510589_dp

ELSE IF(element_symbol == "Fe" .AND. charge == 16) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.20145020_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.03204342_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.03204342_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.07582857_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.01032228_dp

ELSE IF(element_symbol == "Fe" .AND. charge == 8) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.01490858_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.02104142_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.02104142_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.04979310_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00572210_dp

ELSE IF(element_symbol == "Ga" .AND. charge == 13) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.07436967_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.04610490_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.04610490_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.05455205_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.01550429_dp

ELSE IF(element_symbol == "Ga" .AND. charge == 3) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.02960710_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.00036899_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.00036899_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.00087320_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00148617_dp

ELSE IF(element_symbol == "Gd" .AND. charge == 18) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.53553977_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.00930685_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.00930685_dp
   gth_potential%kprj_ppnl(1,3,1) = -0.01614340_dp
   gth_potential%kprj_ppnl(3,1,1) = -0.01614340_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.02202404_dp
   gth_potential%kprj_ppnl(2,3,1) = 0.03820226_dp
   gth_potential%kprj_ppnl(3,2,1) = 0.03820226_dp
   gth_potential%kprj_ppnl(3,3,1) = -0.05430110_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.05312761_dp
   gth_potential%kprj_ppnl(1,1,3) = 0.01466612_dp

ELSE IF(element_symbol == "Ge" .AND. charge == 4) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.04398121_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.00414190_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.00414190_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.00980152_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00291785_dp

ELSE IF(element_symbol == "Hf" .AND. charge == 12) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.60750427_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.14014208_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.14014208_dp
   gth_potential%kprj_ppnl(1,3,1) = -0.03597879_dp
   gth_potential%kprj_ppnl(3,1,1) = -0.03597879_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.33163669_dp
   gth_potential%kprj_ppnl(2,3,1) = 0.08514135_dp
   gth_potential%kprj_ppnl(3,2,1) = 0.08514135_dp
   gth_potential%kprj_ppnl(3,3,1) = -0.12102083_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.22211861_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.05348066_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.05348066_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.12128275_dp

ELSE IF(element_symbol == "Hg" .AND. charge == 12) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.08498860_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.03075146_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.03075146_dp
   gth_potential%kprj_ppnl(1,3,1) = 0.19423658_dp
   gth_potential%kprj_ppnl(3,1,1) = 0.19423658_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.07277124_dp
   gth_potential%kprj_ppnl(2,3,1) = -0.45964764_dp
   gth_potential%kprj_ppnl(3,2,1) = -0.45964764_dp
   gth_potential%kprj_ppnl(3,3,1) = 0.65334804_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.15575928_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.05392128_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.05392128_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.12228196_dp

ELSE IF(element_symbol == "Hg" .AND. charge == 2) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.09232980_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.00047238_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.00047238_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.00111786_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.02093075_dp

ELSE IF(element_symbol == "Ho" .AND. charge == 21) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.89594677_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.11627162_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.11627162_dp
   gth_potential%kprj_ppnl(1,3,1) = 0.06116207_dp
   gth_potential%kprj_ppnl(3,1,1) = 0.06116207_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.27514887_dp
   gth_potential%kprj_ppnl(2,3,1) = -0.14473587_dp
   gth_potential%kprj_ppnl(3,2,1) = -0.14473587_dp
   gth_potential%kprj_ppnl(3,3,1) = 0.20572910_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.06948258_dp
   gth_potential%kprj_ppnl(1,1,3) = 0.01983007_dp

ELSE IF(element_symbol == "I " .AND. charge == 7) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.21396743_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.03943200_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.03943200_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.09331314_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.01018039_dp

ELSE IF(element_symbol == "In" .AND. charge == 13) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.30648842_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.20259839_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.20259839_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.47943530_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.05876513_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.03237261_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.03237261_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.07341418_dp

ELSE IF(element_symbol == "In" .AND. charge == 3) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.06636653_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.00235065_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.00235065_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.00556265_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00444800_dp

ELSE IF(element_symbol == "Ir" .AND. charge == 17) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.93012073_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.16052213_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.16052213_dp
   gth_potential%kprj_ppnl(1,3,1) = -0.07810385_dp
   gth_potential%kprj_ppnl(3,1,1) = -0.07810385_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.37986468_dp
   gth_potential%kprj_ppnl(2,3,1) = 0.18482745_dp
   gth_potential%kprj_ppnl(3,2,1) = 0.18482745_dp
   gth_potential%kprj_ppnl(3,3,1) = -0.26271570_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.21951741_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.04945132_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.04945132_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.11214505_dp

ELSE IF(element_symbol == "Ir" .AND. charge == 9) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.56534735_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.36325582_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.36325582_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.85962018_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.12719880_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.06700352_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.06700352_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.15194970_dp

ELSE IF(element_symbol == "K " .AND. charge == 1) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.00434261_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.00433626_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.00433626_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.01026146_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00035395_dp

ELSE IF(element_symbol == "K " .AND. charge == 9) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.01579489_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.02067379_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.02067379_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.04892313_dp

ELSE IF(element_symbol == "Kr" .AND. charge == 8) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.08701098_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.03634480_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.03634480_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.08600750_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00926681_dp

ELSE IF(element_symbol == "La" .AND. charge == 11) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.52462281_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.01305818_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.01305818_dp
   gth_potential%kprj_ppnl(1,3,1) = 0.04223862_dp
   gth_potential%kprj_ppnl(3,1,1) = 0.04223862_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.03090129_dp
   gth_potential%kprj_ppnl(2,3,1) = -0.09995482_dp
   gth_potential%kprj_ppnl(3,2,1) = -0.09995482_dp
   gth_potential%kprj_ppnl(3,3,1) = 0.14207685_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.02089998_dp
   gth_potential%kprj_ppnl(1,1,3) = 0.00719271_dp

ELSE IF(element_symbol == "Li" .AND. charge == 1) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.00001876_dp

ELSE IF(element_symbol == "Lu" .AND. charge == 25) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.15245040_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.58967098_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.58967098_dp
   gth_potential%kprj_ppnl(1,3,1) = -0.36827142_dp
   gth_potential%kprj_ppnl(3,1,1) = -0.36827142_dp
   gth_potential%kprj_ppnl(2,2,1) = 1.39541622_dp
   gth_potential%kprj_ppnl(2,3,1) = 0.87148924_dp
   gth_potential%kprj_ppnl(3,2,1) = 0.87148924_dp
   gth_potential%kprj_ppnl(3,3,1) = -1.23874407_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.07213023_dp
   gth_potential%kprj_ppnl(1,1,3) = 0.02800562_dp

ELSE IF(element_symbol == "Mg" .AND. charge == 10) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.09627735_dp

ELSE IF(element_symbol == "Mg" .AND. charge == 2) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.00515157_dp

ELSE IF(element_symbol == "Mn" .AND. charge == 15) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.16418798_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.01664776_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.01664776_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.03939579_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00958994_dp

ELSE IF(element_symbol == "Mn" .AND. charge == 7) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.01368526_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.01814464_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.01814464_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.04293805_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00453641_dp

ELSE IF(element_symbol == "Mo" .AND. charge == 14) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.37868133_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.05263679_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.05263679_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.12456138_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.01445963_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.00651263_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.00651263_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.01476925_dp

ELSE IF(element_symbol == "Mo" .AND. charge == 6) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.10412389_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.06503729_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.06503729_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.15390632_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.06897223_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.03333236_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.03333236_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.07559069_dp

ELSE IF(element_symbol == "N " .AND. charge == 5) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.00313068_dp

ELSE IF(element_symbol == "Na" .AND. charge == 1) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.00262341_dp

ELSE IF(element_symbol == "Na" .AND. charge == 9) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.03838599_dp

ELSE IF(element_symbol == "Nb" .AND. charge == 13) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.24682052_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.03662021_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.03662021_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.08665924_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.01406431_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.00531597_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.00531597_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.01205549_dp

ELSE IF(element_symbol == "Nb" .AND. charge == 5) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.08081558_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.05292491_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.05292491_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.12524320_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.02565344_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.01390823_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.01390823_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.03154090_dp

ELSE IF(element_symbol == "Nd" .AND. charge == 14) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.21439589_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.46502350_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.46502350_dp
   gth_potential%kprj_ppnl(1,3,1) = -0.24754787_dp
   gth_potential%kprj_ppnl(3,1,1) = -0.24754787_dp
   gth_potential%kprj_ppnl(2,2,1) = 1.10044645_dp
   gth_potential%kprj_ppnl(2,3,1) = 0.58580517_dp
   gth_potential%kprj_ppnl(3,2,1) = 0.58580517_dp
   gth_potential%kprj_ppnl(3,3,1) = -0.83266969_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.61141305_dp
   gth_potential%kprj_ppnl(1,1,3) = 0.00959761_dp

ELSE IF(element_symbol == "Ne" .AND. charge == 8) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.01033644_dp

ELSE IF(element_symbol == "Ni" .AND. charge == 10) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.02038398_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.02906079_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.02906079_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.06877038_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00870832_dp

ELSE IF(element_symbol == "Ni" .AND. charge == 18) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.26957208_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.06061531_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.06061531_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.14344200_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.01353850_dp

ELSE IF(element_symbol == "O " .AND. charge == 6) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.00447578_dp

ELSE IF(element_symbol == "Os" .AND. charge == 16) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.87081735_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.17447632_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.17447632_dp
   gth_potential%kprj_ppnl(1,3,1) = -0.04147424_dp
   gth_potential%kprj_ppnl(3,1,1) = -0.04147424_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.41288633_dp
   gth_potential%kprj_ppnl(2,3,1) = 0.09814596_dp
   gth_potential%kprj_ppnl(3,2,1) = 0.09814596_dp
   gth_potential%kprj_ppnl(3,3,1) = -0.13950571_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.21644008_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.04551632_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.04551632_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.10322131_dp

ELSE IF(element_symbol == "Os" .AND. charge == 8) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.43074635_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.29654416_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.29654416_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.70175157_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.11694065_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.06162898_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.06162898_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.13976139_dp

ELSE IF(element_symbol == "P " .AND. charge == 5) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.00254425_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.00756218_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.00756218_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.01789537_dp

ELSE IF(element_symbol == "Pb" .AND. charge == 4) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.20771072_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.00547158_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.00547158_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.01294813_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.02925557_dp

ELSE IF(element_symbol == "Pd" .AND. charge == 10) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.13690890_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.09330707_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.09330707_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.22080482_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.03379697_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.02107700_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.02107700_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.04779814_dp

ELSE IF(element_symbol == "Pd" .AND. charge == 18) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.40698386_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.04680631_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.04680631_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.11076394_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.02037404_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.01178781_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.01178781_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.02673224_dp

ELSE IF(element_symbol == "Pm" .AND. charge == 15) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.33995543_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.57428953_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.57428953_dp
   gth_potential%kprj_ppnl(1,3,1) = -0.34066424_dp
   gth_potential%kprj_ppnl(3,1,1) = -0.34066424_dp
   gth_potential%kprj_ppnl(2,2,1) = 1.35901707_dp
   gth_potential%kprj_ppnl(2,3,1) = 0.80615873_dp
   gth_potential%kprj_ppnl(3,2,1) = 0.80615873_dp
   gth_potential%kprj_ppnl(3,3,1) = -1.14588259_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.06404352_dp
   gth_potential%kprj_ppnl(1,1,3) = 0.01061852_dp

ELSE IF(element_symbol == "Po" .AND. charge == 6) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.39635435_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.01329517_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.01329517_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.03146213_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.03388641_dp

ELSE IF(element_symbol == "Pr" .AND. charge == 13) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.31427984_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.14831690_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.14831690_dp
   gth_potential%kprj_ppnl(1,3,1) = -0.04912915_dp
   gth_potential%kprj_ppnl(3,1,1) = -0.04912915_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.35098185_dp
   gth_potential%kprj_ppnl(2,3,1) = 0.11626078_dp
   gth_potential%kprj_ppnl(3,2,1) = 0.11626078_dp
   gth_potential%kprj_ppnl(3,3,1) = -0.16525431_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.01090482_dp
   gth_potential%kprj_ppnl(1,1,3) = 0.00854696_dp

ELSE IF(element_symbol == "Pt" .AND. charge == 10) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.76329581_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.45041788_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.45041788_dp
   gth_potential%kprj_ppnl(2,2,1) = 1.06588325_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.14691239_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.07465713_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.07465713_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.16930645_dp

ELSE IF(element_symbol == "Pt" .AND. charge == 18) THEN
   gth_potential%kprj_ppnl(1,1,1) = 1.01705973_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.14714667_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.14714667_dp
   gth_potential%kprj_ppnl(1,3,1) = -0.09867767_dp
   gth_potential%kprj_ppnl(3,1,1) = -0.09867767_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.34821257_dp
   gth_potential%kprj_ppnl(2,3,1) = 0.23351399_dp
   gth_potential%kprj_ppnl(3,2,1) = 0.23351399_dp
   gth_potential%kprj_ppnl(3,3,1) = -0.33191927_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.22647170_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.05042175_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.05042175_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.11434578_dp

ELSE IF(element_symbol == "Rb" .AND. charge == 1) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.04344330_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.02445691_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.02445691_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.05787561_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00370757_dp

ELSE IF(element_symbol == "Rb" .AND. charge == 9) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.86737860_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.52295274_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.52295274_dp
   gth_potential%kprj_ppnl(2,2,1) = 1.23753204_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00868497_dp

ELSE IF(element_symbol == "Re" .AND. charge == 15) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.78871453_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.20705623_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.20705623_dp
   gth_potential%kprj_ppnl(1,3,1) = -0.00505749_dp
   gth_potential%kprj_ppnl(3,1,1) = -0.00505749_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.48998446_dp
   gth_potential%kprj_ppnl(2,3,1) = 0.01196821_dp
   gth_potential%kprj_ppnl(3,2,1) = 0.01196821_dp
   gth_potential%kprj_ppnl(3,3,1) = -0.01701174_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.20973732_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.04535809_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.04535809_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.10286248_dp

ELSE IF(element_symbol == "Re" .AND. charge == 7) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.37058909_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.26063064_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.26063064_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.61676466_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.11155664_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.05802813_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.05802813_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.13159543_dp

ELSE IF(element_symbol == "Rh" .AND. charge == 17) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.38962932_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.05191608_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.05191608_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.12285587_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.01887482_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.01007017_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.01007017_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.02283699_dp

ELSE IF(element_symbol == "Rh" .AND. charge == 9) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.11435695_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.08194750_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.08194750_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.19392318_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.13503645_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.05633249_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.05633249_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.12775008_dp

ELSE IF(element_symbol == "Rn" .AND. charge == 8) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.54989620_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.01004055_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.01004055_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.02376028_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.04548755_dp

ELSE IF(element_symbol == "Ru" .AND. charge == 16) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.33732182_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.04372298_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.04372298_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.10346747_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.01777313_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.00869809_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.00869809_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.01972541_dp

ELSE IF(element_symbol == "Ru" .AND. charge == 8) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.09545431_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.06941124_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.06941124_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.16425697_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.04665238_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.02725483_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.02725483_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.06180814_dp

ELSE IF(element_symbol == "S " .AND. charge == 6) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.00537161_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.00932309_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.00932309_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.02206245_dp

ELSE IF(element_symbol == "Sb" .AND. charge == 5) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.13922153_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.00993605_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.00993605_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.02351298_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00943232_dp

ELSE IF(element_symbol == "Sc" .AND. charge == 11) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.08305262_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.03930672_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.03930672_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.09301668_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00481191_dp

ELSE IF(element_symbol == "Sc" .AND. charge == 3) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.00846271_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.00883730_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.00883730_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.02091287_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00142986_dp

ELSE IF(element_symbol == "Se" .AND. charge == 6) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.06219585_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.02742836_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.02742836_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.06490734_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00578395_dp

ELSE IF(element_symbol == "Si" .AND. charge == 4) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.00037280_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.00610073_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.00610073_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.01443697_dp

ELSE IF(element_symbol == "Sm" .AND. charge == 16) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.24029994_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.50745883_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.50745883_dp
   gth_potential%kprj_ppnl(1,3,1) = -0.31336038_dp
   gth_potential%kprj_ppnl(3,1,1) = -0.31336038_dp
   gth_potential%kprj_ppnl(2,2,1) = 1.20086678_dp
   gth_potential%kprj_ppnl(2,3,1) = 0.74154601_dp
   gth_potential%kprj_ppnl(3,2,1) = 0.74154601_dp
   gth_potential%kprj_ppnl(3,3,1) = -1.05404138_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.06335216_dp
   gth_potential%kprj_ppnl(1,1,3) = 0.01192432_dp

ELSE IF(element_symbol == "Sn" .AND. charge == 4) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.10393057_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.00213712_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.00213712_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.00505735_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00706571_dp

ELSE IF(element_symbol == "Sr" .AND. charge == 10) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.57626460_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.41837761_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.41837761_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.99006212_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00899140_dp

ELSE IF(element_symbol == "Sr" .AND. charge == 2) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.00402188_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.00938423_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.00938423_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.02220714_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00284555_dp

ELSE IF(element_symbol == "Ta" .AND. charge == 13) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.64999172_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.14214268_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.14214268_dp
   gth_potential%kprj_ppnl(1,3,1) = -0.03012256_dp
   gth_potential%kprj_ppnl(3,1,1) = -0.03012256_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.33637097_dp
   gth_potential%kprj_ppnl(2,3,1) = 0.07128298_dp
   gth_potential%kprj_ppnl(3,2,1) = 0.07128298_dp
   gth_potential%kprj_ppnl(3,3,1) = -0.10132239_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.20534356_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.04513357_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.04513357_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.10235332_dp

ELSE IF(element_symbol == "Ta" .AND. charge == 5) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.48563525_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.29375358_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.29375358_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.69514785_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.08671559_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.04173012_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.04173012_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.09463501_dp

ELSE IF(element_symbol == "Tb" .AND. charge == 19) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.00872499_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.48343103_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.48343103_dp
   gth_potential%kprj_ppnl(1,3,1) = -0.25579926_dp
   gth_potential%kprj_ppnl(3,1,1) = -0.25579926_dp
   gth_potential%kprj_ppnl(2,2,1) = 1.14400661_dp
   gth_potential%kprj_ppnl(2,3,1) = 0.60533153_dp
   gth_potential%kprj_ppnl(3,2,1) = 0.60533153_dp
   gth_potential%kprj_ppnl(3,3,1) = -0.86042467_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.05175355_dp
   gth_potential%kprj_ppnl(1,1,3) = 0.01619671_dp

ELSE IF(element_symbol == "Tc" .AND. charge == 15) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.34079130_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.02786799_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.02786799_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.06594771_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.01578982_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.00726905_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.00726905_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.01648465_dp

ELSE IF(element_symbol == "Tc" .AND. charge == 7) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.08075975_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.05944302_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.05944302_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.14066787_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.02602474_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.01842157_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.01842157_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.04177619_dp

ELSE IF(element_symbol == "Te" .AND. charge == 6) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.17299713_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.02139985_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.02139985_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.05064130_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.01080862_dp

ELSE IF(element_symbol == "Ti" .AND. charge == 12) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.12239544_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.02427591_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.02427591_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.05744728_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00582242_dp

ELSE IF(element_symbol == "Ti" .AND. charge == 4) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.00933253_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.01068729_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.01068729_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.02529075_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00200953_dp

ELSE IF(element_symbol == "Tl" .AND. charge == 13) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.01223345_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.01338045_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.01338045_dp
   gth_potential%kprj_ppnl(1,3,1) = 0.30299168_dp
   gth_potential%kprj_ppnl(3,1,1) = 0.30299168_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.03166391_dp
   gth_potential%kprj_ppnl(2,3,1) = -0.71700919_dp
   gth_potential%kprj_ppnl(3,2,1) = -0.71700919_dp
   gth_potential%kprj_ppnl(3,3,1) = 1.01916448_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.18684918_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.07524992_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.07524992_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.17065079_dp

ELSE IF(element_symbol == "Tl" .AND. charge == 3) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.16864149_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.00188436_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.00188436_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.00445921_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.02294091_dp

ELSE IF(element_symbol == "Tm" .AND. charge == 23) THEN
   gth_potential%kprj_ppnl(1,1,1) = 2.04536983_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.99188883_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.99188883_dp
   gth_potential%kprj_ppnl(1,3,1) = 0.56567230_dp
   gth_potential%kprj_ppnl(3,1,1) = 0.56567230_dp
   gth_potential%kprj_ppnl(2,2,1) = -2.34723738_dp
   gth_potential%kprj_ppnl(2,3,1) = -1.33862499_dp
   gth_potential%kprj_ppnl(3,2,1) = -1.33862499_dp
   gth_potential%kprj_ppnl(3,3,1) = 1.90273579_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.11929895_dp
   gth_potential%kprj_ppnl(1,1,3) = 0.02383347_dp

ELSE IF(element_symbol == "V " .AND. charge == 13) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.15640785_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.00339312_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.00339312_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.00802959_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00654774_dp

ELSE IF(element_symbol == "V " .AND. charge == 5) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.01097329_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.01302199_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.01302199_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.03081565_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00271707_dp

ELSE IF(element_symbol == "W " .AND. charge == 14) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.71547973_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.16283543_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.16283543_dp
   gth_potential%kprj_ppnl(1,3,1) = -0.02291919_dp
   gth_potential%kprj_ppnl(3,1,1) = -0.02291919_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.38533896_dp
   gth_potential%kprj_ppnl(2,3,1) = 0.05423670_dp
   gth_potential%kprj_ppnl(3,2,1) = 0.05423670_dp
   gth_potential%kprj_ppnl(3,3,1) = -0.07709262_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.20586005_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.04411598_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.04411598_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.10004563_dp

ELSE IF(element_symbol == "W " .AND. charge == 6) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.50979956_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.31766785_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.31766785_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.75173934_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.07577234_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.03800751_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.03800751_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.08619292_dp

ELSE IF(element_symbol == "Xe" .AND. charge == 8) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.23660318_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.04583840_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.04583840_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.10847345_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.01336223_dp

ELSE IF(element_symbol == "Y " .AND. charge == 11) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.68452714_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.49749106_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.49749106_dp
   gth_potential%kprj_ppnl(2,2,1) = 1.17727871_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.01045632_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.00135148_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.00135148_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.00306487_dp

ELSE IF(element_symbol == "Y " .AND. charge == 3) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.04333632_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.03380154_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.03380154_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.07998905_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00919822_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.00514005_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.00514005_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.01165653_dp

ELSE IF(element_symbol == "Yb" .AND. charge == 24) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.34946961_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.87654953_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.87654953_dp
   gth_potential%kprj_ppnl(1,3,1) = -0.53937986_dp
   gth_potential%kprj_ppnl(3,1,1) = -0.53937986_dp
   gth_potential%kprj_ppnl(2,2,1) = 2.07429479_dp
   gth_potential%kprj_ppnl(2,3,1) = 1.27640570_dp
   gth_potential%kprj_ppnl(3,2,1) = 1.27640570_dp
   gth_potential%kprj_ppnl(3,3,1) = -1.81429663_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.07007567_dp
   gth_potential%kprj_ppnl(1,1,3) = 0.02571848_dp

ELSE IF(element_symbol == "Zn" .AND. charge == 12) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.02575928_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.03841864_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.03841864_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.09091511_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.01276688_dp

ELSE IF(element_symbol == "Zn" .AND. charge == 2) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.01213948_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.00206035_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.00206035_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.00487568_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00022546_dp

ELSE IF(element_symbol == "Zn" .AND. charge == 20) THEN
   gth_potential%kprj_ppnl(1,1,1) = 0.06195779_dp
   gth_potential%kprj_ppnl(1,2,1) = 0.19739272_dp
   gth_potential%kprj_ppnl(2,1,1) = 0.19739272_dp
   gth_potential%kprj_ppnl(2,2,1) = -0.23355822_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.01770736_dp

ELSE IF(element_symbol == "Zr" .AND. charge == 12) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.51261281_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.44782121_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.44782121_dp
   gth_potential%kprj_ppnl(2,2,1) = 1.05973841_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00756230_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.00237820_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.00237820_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.00539325_dp

ELSE IF(element_symbol == "Zr" .AND. charge == 4) THEN
   gth_potential%kprj_ppnl(1,1,1) = -0.05748649_dp
   gth_potential%kprj_ppnl(1,2,1) = -0.04415700_dp
   gth_potential%kprj_ppnl(2,1,1) = -0.04415700_dp
   gth_potential%kprj_ppnl(2,2,1) = 0.10449452_dp
   gth_potential%kprj_ppnl(1,1,2) = 0.00938027_dp
   gth_potential%kprj_ppnl(1,2,2) = 0.00527960_dp
   gth_potential%kprj_ppnl(2,1,2) = 0.00527960_dp
   gth_potential%kprj_ppnl(2,2,2) = -0.01197300_dp
END IF

      CALL timestop(handle)

END SUBROUTINE get_kprj_SOC_parameters_LDA

END MODULE external_potential_soc_parameters
