!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief SOC parameters of Hartwigsten-Goedecker-Hutter pseudopotentials
!>        (temporary file before storage of SOC parameters in data directory)
!> \par History
!>      06.2023 created [Jan Wilhelm]
! **************************************************************************************************
MODULE external_potential_soc_parameters
   USE external_potential_types,        ONLY: gth_potential_type,&
                                              init_potential
   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'external_potential_soc_parameters'

   PUBLIC :: get_kprj_SOC_parameters

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param gth_potential ...
!> \param xc_functional ...
!> \param element_symbol ...
! **************************************************************************************************
   SUBROUTINE get_kprj_SOC_parameters(gth_potential, xc_functional, element_symbol)

      TYPE(gth_potential_type), POINTER                  :: gth_potential
      CHARACTER(LEN=*)                                   :: xc_functional
      CHARACTER(LEN=2)                                   :: element_symbol

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_kprj_SOC_parameters'

      INTEGER                                            :: charge, handle

      CALL timeset(routineN, handle)

      charge = NINT(gth_potential%zeff)

      gth_potential%kprj_ppnl(:, :, :) = 0.0_dp
      gth_potential%wprj_ppnl(:, :) = 0.0_dp

      SELECT CASE (xc_functional)
      CASE ("LDA")
         CALL get_kprj_SOC_parameters_LDA(gth_potential, element_symbol, charge)
      CASE ("PBE")
         CALL get_kprj_SOC_parameters_PBE(gth_potential, element_symbol, charge)
      CASE DEFAULT
         CPABORT("Unknown xc functional in SOC part of pseudopotential.")
      END SELECT

      ! write SOC parameters from kprj_ppnl into wprj_ppnl
      CALL init_potential(gth_potential)

      CALL timestop(handle)

   END SUBROUTINE get_kprj_SOC_parameters

! **************************************************************************************************
!> \brief ...
!> \param gth_potential ...
!> \param element_symbol ...
!> \param charge ...
! **************************************************************************************************
   SUBROUTINE get_kprj_SOC_parameters_LDA(gth_potential, element_symbol, charge)

      TYPE(gth_potential_type), POINTER                  :: gth_potential
      CHARACTER(LEN=2)                                   :: element_symbol
      INTEGER                                            :: charge

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_kprj_SOC_parameters_LDA'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (element_symbol == "Ag" .AND. charge == 1) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.01969167_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.00288247_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.00288247_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.00682116_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00945505_dp

      ELSE IF (element_symbol == "Ag" .AND. charge == 11) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.20186901_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.13510134_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.13510134_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.31970813_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.04519317_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.02273792_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.02273792_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.05156476_dp

      ELSE IF (element_symbol == "Ag" .AND. charge == 19) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.26670728_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.68349949_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.68349949_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.80872750_dp
         gth_potential%kprj_ppnl(1, 1, 2) = -0.03218230_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.06201267_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.06201267_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.07031576_dp

      ELSE IF (element_symbol == "Al" .AND. charge == 3) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.00615404_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.00166789_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.00166789_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.00394694_dp

      ELSE IF (element_symbol == "Ar" .AND. charge == 8) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.01639467_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.01232694_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.01232694_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.02917087_dp

      ELSE IF (element_symbol == "As" .AND. charge == 5) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.05246616_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.00868894_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.00868894_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.02056179_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00427335_dp

      ELSE IF (element_symbol == "At" .AND. charge == 7) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.48077365_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.01477071_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.01477071_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.03495387_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.03754438_dp

      ELSE IF (element_symbol == "Au" .AND. charge == 1) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.03934865_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.05618995_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.05618995_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -0.04561301_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -0.04561301_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.13296969_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 0.10794007_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 0.10794007_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -0.15342717_dp

      ELSE IF (element_symbol == "Au" .AND. charge == 11) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.96005489_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.54790119_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.54790119_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 1.29657086_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.14848354_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.07473921_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.07473921_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.16949260_dp

      ELSE IF (element_symbol == "Au" .AND. charge == 19) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.62807493_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 1.81668773_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 1.81668773_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -1.65141098_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -1.65141098_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -4.89888134_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 3.90795165_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 3.90795165_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -2.77740200_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.18880114_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.03175114_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.03175114_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.03600240_dp

      ELSE IF (element_symbol == "B " .AND. charge == 3) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.00087786_dp

      ELSE IF (element_symbol == "Ba" .AND. charge == 10) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.41508327_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.25344615_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.25344615_dp
         gth_potential%kprj_ppnl(1, 3, 1) = 0.27561093_dp
         gth_potential%kprj_ppnl(3, 1, 1) = 0.27561093_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.59976305_dp
         gth_potential%kprj_ppnl(2, 3, 1) = -0.65221449_dp
         gth_potential%kprj_ppnl(3, 2, 1) = -0.65221449_dp
         gth_potential%kprj_ppnl(3, 3, 1) = 0.92706460_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01766026_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.00615085_dp

      ELSE IF (element_symbol == "Ba" .AND. charge == 2) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.04101721_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.03752472_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.03752472_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.08879970_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00532055_dp

      ELSE IF (element_symbol == "Be" .AND. charge == 2) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.00012877_dp

      ELSE IF (element_symbol == "Bi" .AND. charge == 5) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.30531353_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.00977575_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.00977575_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.02313366_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.02921737_dp

      ELSE IF (element_symbol == "Br" .AND. charge == 7) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.07400663_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.02906775_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.02906775_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.06878684_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00714446_dp

      ELSE IF (element_symbol == "C " .AND. charge == 4) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.00410365_dp

      ELSE IF (element_symbol == "Ca" .AND. charge == 10) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.02170091_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.02394376_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.02394376_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.05666128_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00037106_dp

      ELSE IF (element_symbol == "Ca" .AND. charge == 2) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.00336201_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.00540007_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.00540007_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.01277891_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00081449_dp

      ELSE IF (element_symbol == "Cd" .AND. charge == 12) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.24140932_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.15923958_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.15923958_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.37682962_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.04932879_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.02594377_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.02594377_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.05883495_dp

      ELSE IF (element_symbol == "Cd" .AND. charge == 2) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.03771366_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.00084877_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.00084877_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.00200856_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00310866_dp

      ELSE IF (element_symbol == "Ce" .AND. charge == 12) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.55905862_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.30893198_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.30893198_dp
         gth_potential%kprj_ppnl(1, 3, 1) = 0.22367871_dp
         gth_potential%kprj_ppnl(3, 1, 1) = 0.22367871_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.73792420_dp
         gth_potential%kprj_ppnl(2, 3, 1) = -0.52932045_dp
         gth_potential%kprj_ppnl(3, 2, 1) = -0.52932045_dp
         gth_potential%kprj_ppnl(3, 3, 1) = 0.37619085_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01403180_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.00781271_dp

      ELSE IF (element_symbol == "Cl" .AND. charge == 7) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.01002014_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.01089588_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.01089588_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.02578435_dp

      ELSE IF (element_symbol == "Co" .AND. charge == 17) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.20791470_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.04574353_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.04574353_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.10824895_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01147529_dp

      ELSE IF (element_symbol == "Co" .AND. charge == 9) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.01752129_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.02483306_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.02483306_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.05876575_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00713682_dp

      ELSE IF (element_symbol == "Cr" .AND. charge == 14) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.16978070_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.00017378_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.00017378_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.00041123_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00900712_dp

      ELSE IF (element_symbol == "Cr" .AND. charge == 6) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.01317649_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.01505450_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.01505450_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.03562545_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00351409_dp

      ELSE IF (element_symbol == "Cs" .AND. charge == 1) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.10762747_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.05837136_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.05837136_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.13813185_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00325918_dp

      ELSE IF (element_symbol == "Cs" .AND. charge == 9) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -2.74852756_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.76700684_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.76700684_dp
         gth_potential%kprj_ppnl(1, 3, 1) = 0.35818758_dp
         gth_potential%kprj_ppnl(3, 1, 1) = 0.35818758_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 1.81506946_dp
         gth_potential%kprj_ppnl(2, 3, 1) = -0.84762652_dp
         gth_potential%kprj_ppnl(3, 2, 1) = -0.84762652_dp
         gth_potential%kprj_ppnl(3, 3, 1) = 1.20482534_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01084084_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.00475984_dp

      ELSE IF (element_symbol == "Cu" .AND. charge == 1) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.01079247_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.00284563_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.00284563_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.00673399_dp
         gth_potential%kprj_ppnl(1, 1, 2) = -0.00073000_dp

      ELSE IF (element_symbol == "Cu" .AND. charge == 11) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.02406730_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.03231927_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.03231927_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.07648136_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01048914_dp

      ELSE IF (element_symbol == "Cu" .AND. charge == 19) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.13410535_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.21399818_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.21399818_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.25320605_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01450302_dp

      ELSE IF (element_symbol == "Dy" .AND. charge == 20) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.12437116_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.56093057_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.56093057_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -0.32468735_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -0.32468735_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 1.32740399_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 0.76835051_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 0.76835051_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -1.09214158_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.05828571_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.01793827_dp

      ELSE IF (element_symbol == "Er" .AND. charge == 22) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 1.25513474_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.37692191_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.37692191_dp
         gth_potential%kprj_ppnl(1, 3, 1) = 0.22288650_dp
         gth_potential%kprj_ppnl(3, 1, 1) = 0.22288650_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.89196004_dp
         gth_potential%kprj_ppnl(2, 3, 1) = -0.52744572_dp
         gth_potential%kprj_ppnl(3, 2, 1) = -0.52744572_dp
         gth_potential%kprj_ppnl(3, 3, 1) = 0.74971696_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.09186175_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.02178337_dp

      ELSE IF (element_symbol == "Eu" .AND. charge == 17) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.25225840_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.24692195_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.24692195_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -0.14287549_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -0.14287549_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.58432399_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 0.33810512_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 0.33810512_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -0.48058621_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.05102753_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.01326718_dp

      ELSE IF (element_symbol == "F " .AND. charge == 7) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.01510589_dp

      ELSE IF (element_symbol == "Fe" .AND. charge == 16) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.20145020_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.03204342_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.03204342_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.07582857_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01032228_dp

      ELSE IF (element_symbol == "Fe" .AND. charge == 8) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.01490858_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.02104142_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.02104142_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.04979310_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00572210_dp

      ELSE IF (element_symbol == "Ga" .AND. charge == 13) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.07436967_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.04610490_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.04610490_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.05455205_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01550429_dp

      ELSE IF (element_symbol == "Ga" .AND. charge == 3) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.02960710_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.00036899_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.00036899_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.00087320_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00148617_dp

      ELSE IF (element_symbol == "Gd" .AND. charge == 18) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.53553977_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.00930685_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.00930685_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -0.01614340_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -0.01614340_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.02202404_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 0.03820226_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 0.03820226_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -0.05430110_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.05312761_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.01466612_dp

      ELSE IF (element_symbol == "Ge" .AND. charge == 4) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.04398121_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.00414190_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.00414190_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.00980152_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00291785_dp

      ELSE IF (element_symbol == "Hf" .AND. charge == 12) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.60750427_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.14014208_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.14014208_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -0.03597879_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -0.03597879_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.33163669_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 0.08514135_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 0.08514135_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -0.12102083_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.22211861_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.05348066_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.05348066_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.12128275_dp

      ELSE IF (element_symbol == "Hg" .AND. charge == 12) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.08498860_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.03075146_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.03075146_dp
         gth_potential%kprj_ppnl(1, 3, 1) = 0.19423658_dp
         gth_potential%kprj_ppnl(3, 1, 1) = 0.19423658_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.07277124_dp
         gth_potential%kprj_ppnl(2, 3, 1) = -0.45964764_dp
         gth_potential%kprj_ppnl(3, 2, 1) = -0.45964764_dp
         gth_potential%kprj_ppnl(3, 3, 1) = 0.65334804_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.15575928_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.05392128_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.05392128_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.12228196_dp

      ELSE IF (element_symbol == "Hg" .AND. charge == 2) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.09232980_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.00047238_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.00047238_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.00111786_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.02093075_dp

      ELSE IF (element_symbol == "Ho" .AND. charge == 21) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.89594677_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.11627162_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.11627162_dp
         gth_potential%kprj_ppnl(1, 3, 1) = 0.06116207_dp
         gth_potential%kprj_ppnl(3, 1, 1) = 0.06116207_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.27514887_dp
         gth_potential%kprj_ppnl(2, 3, 1) = -0.14473587_dp
         gth_potential%kprj_ppnl(3, 2, 1) = -0.14473587_dp
         gth_potential%kprj_ppnl(3, 3, 1) = 0.20572910_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.06948258_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.01983007_dp

      ELSE IF (element_symbol == "I " .AND. charge == 7) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.21396743_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.03943200_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.03943200_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.09331314_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01018039_dp

      ELSE IF (element_symbol == "In" .AND. charge == 13) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.30648842_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.20259839_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.20259839_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.47943530_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.05876513_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.03237261_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.03237261_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.07341418_dp

      ELSE IF (element_symbol == "In" .AND. charge == 3) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.06636653_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.00235065_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.00235065_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.00556265_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00444800_dp

      ELSE IF (element_symbol == "Ir" .AND. charge == 17) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.93012073_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.16052213_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.16052213_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -0.07810385_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -0.07810385_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.37986468_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 0.18482745_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 0.18482745_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -0.26271570_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.21951741_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.04945132_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.04945132_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.11214505_dp

      ELSE IF (element_symbol == "Ir" .AND. charge == 9) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.56534735_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.36325582_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.36325582_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.85962018_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.12719880_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.06700352_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.06700352_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.15194970_dp

      ELSE IF (element_symbol == "K " .AND. charge == 1) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.00434261_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.00433626_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.00433626_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.01026146_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00035395_dp

      ELSE IF (element_symbol == "K " .AND. charge == 9) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.01579489_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.02067379_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.02067379_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.04892313_dp

      ELSE IF (element_symbol == "Kr" .AND. charge == 8) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.08701098_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.03634480_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.03634480_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.08600750_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00926681_dp

      ELSE IF (element_symbol == "La" .AND. charge == 11) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.52462281_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.01305818_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.01305818_dp
         gth_potential%kprj_ppnl(1, 3, 1) = 0.04223862_dp
         gth_potential%kprj_ppnl(3, 1, 1) = 0.04223862_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.03090129_dp
         gth_potential%kprj_ppnl(2, 3, 1) = -0.09995482_dp
         gth_potential%kprj_ppnl(3, 2, 1) = -0.09995482_dp
         gth_potential%kprj_ppnl(3, 3, 1) = 0.14207685_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.02089998_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.00719271_dp

      ELSE IF (element_symbol == "Li" .AND. charge == 1) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.00001876_dp

      ELSE IF (element_symbol == "Lu" .AND. charge == 25) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.15245040_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.58967098_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.58967098_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -0.36827142_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -0.36827142_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 1.39541622_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 0.87148924_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 0.87148924_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -1.23874407_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.07213023_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.02800562_dp

      ELSE IF (element_symbol == "Mg" .AND. charge == 10) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.09627735_dp

      ELSE IF (element_symbol == "Mg" .AND. charge == 2) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.00515157_dp

      ELSE IF (element_symbol == "Mn" .AND. charge == 15) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.16418798_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.01664776_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.01664776_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.03939579_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00958994_dp

      ELSE IF (element_symbol == "Mn" .AND. charge == 7) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.01368526_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.01814464_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.01814464_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.04293805_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00453641_dp

      ELSE IF (element_symbol == "Mo" .AND. charge == 14) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.37868133_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.05263679_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.05263679_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.12456138_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01445963_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.00651263_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.00651263_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.01476925_dp

      ELSE IF (element_symbol == "Mo" .AND. charge == 6) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.10412389_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.06503729_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.06503729_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.15390632_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.06897223_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.03333236_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.03333236_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.07559069_dp

      ELSE IF (element_symbol == "N " .AND. charge == 5) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.00313068_dp

      ELSE IF (element_symbol == "Na" .AND. charge == 1) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.00262341_dp

      ELSE IF (element_symbol == "Na" .AND. charge == 9) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.03838599_dp

      ELSE IF (element_symbol == "Nb" .AND. charge == 13) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.24682052_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.03662021_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.03662021_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.08665924_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01406431_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.00531597_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.00531597_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.01205549_dp

      ELSE IF (element_symbol == "Nb" .AND. charge == 5) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.08081558_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.05292491_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.05292491_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.12524320_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.02565344_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.01390823_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.01390823_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.03154090_dp

      ELSE IF (element_symbol == "Nd" .AND. charge == 14) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.21439589_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.46502350_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.46502350_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -0.24754787_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -0.24754787_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 1.10044645_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 0.58580517_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 0.58580517_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -0.83266969_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.61141305_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.00959761_dp

      ELSE IF (element_symbol == "Ne" .AND. charge == 8) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.01033644_dp

      ELSE IF (element_symbol == "Ni" .AND. charge == 10) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.02038398_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.02906079_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.02906079_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.06877038_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00870832_dp

      ELSE IF (element_symbol == "Ni" .AND. charge == 18) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.26957208_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.06061531_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.06061531_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.14344200_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01353850_dp

      ELSE IF (element_symbol == "O " .AND. charge == 6) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.00447578_dp

      ELSE IF (element_symbol == "Os" .AND. charge == 16) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.87081735_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.17447632_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.17447632_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -0.04147424_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -0.04147424_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.41288633_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 0.09814596_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 0.09814596_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -0.13950571_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.21644008_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.04551632_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.04551632_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.10322131_dp

      ELSE IF (element_symbol == "Os" .AND. charge == 8) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.43074635_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.29654416_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.29654416_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.70175157_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.11694065_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.06162898_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.06162898_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.13976139_dp

      ELSE IF (element_symbol == "P " .AND. charge == 5) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.00254425_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.00756218_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.00756218_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.01789537_dp

      ELSE IF (element_symbol == "Pb" .AND. charge == 4) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.20771072_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.00547158_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.00547158_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.01294813_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.02925557_dp

      ELSE IF (element_symbol == "Pd" .AND. charge == 10) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.13690890_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.09330707_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.09330707_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.22080482_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.03379697_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.02107700_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.02107700_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.04779814_dp

      ELSE IF (element_symbol == "Pd" .AND. charge == 18) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.40698386_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.04680631_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.04680631_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.11076394_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.02037404_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.01178781_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.01178781_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.02673224_dp

      ELSE IF (element_symbol == "Pm" .AND. charge == 15) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.33995543_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.57428953_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.57428953_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -0.34066424_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -0.34066424_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 1.35901707_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 0.80615873_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 0.80615873_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -1.14588259_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.06404352_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.01061852_dp

      ELSE IF (element_symbol == "Po" .AND. charge == 6) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.39635435_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.01329517_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.01329517_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.03146213_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.03388641_dp

      ELSE IF (element_symbol == "Pr" .AND. charge == 13) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.31427984_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.14831690_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.14831690_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -0.04912915_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -0.04912915_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.35098185_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 0.11626078_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 0.11626078_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -0.16525431_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01090482_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.00854696_dp

      ELSE IF (element_symbol == "Pt" .AND. charge == 10) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.76329581_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.45041788_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.45041788_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 1.06588325_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.14691239_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.07465713_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.07465713_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.16930645_dp

      ELSE IF (element_symbol == "Pt" .AND. charge == 18) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 1.01705973_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.14714667_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.14714667_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -0.09867767_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -0.09867767_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.34821257_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 0.23351399_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 0.23351399_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -0.33191927_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.22647170_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.05042175_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.05042175_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.11434578_dp

      ELSE IF (element_symbol == "Rb" .AND. charge == 1) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.04344330_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.02445691_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.02445691_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.05787561_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00370757_dp

      ELSE IF (element_symbol == "Rb" .AND. charge == 9) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.86737860_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.52295274_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.52295274_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 1.23753204_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00868497_dp

      ELSE IF (element_symbol == "Re" .AND. charge == 15) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.78871453_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.20705623_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.20705623_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -0.00505749_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -0.00505749_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.48998446_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 0.01196821_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 0.01196821_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -0.01701174_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.20973732_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.04535809_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.04535809_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.10286248_dp

      ELSE IF (element_symbol == "Re" .AND. charge == 7) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.37058909_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.26063064_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.26063064_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.61676466_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.11155664_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.05802813_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.05802813_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.13159543_dp

      ELSE IF (element_symbol == "Rh" .AND. charge == 17) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.38962932_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.05191608_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.05191608_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.12285587_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01887482_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.01007017_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.01007017_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.02283699_dp

      ELSE IF (element_symbol == "Rh" .AND. charge == 9) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.11435695_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.08194750_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.08194750_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.19392318_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.13503645_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.05633249_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.05633249_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.12775008_dp

      ELSE IF (element_symbol == "Rn" .AND. charge == 8) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.54989620_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.01004055_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.01004055_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.02376028_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.04548755_dp

      ELSE IF (element_symbol == "Ru" .AND. charge == 16) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.33732182_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.04372298_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.04372298_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.10346747_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01777313_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.00869809_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.00869809_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.01972541_dp

      ELSE IF (element_symbol == "Ru" .AND. charge == 8) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.09545431_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.06941124_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.06941124_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.16425697_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.04665238_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.02725483_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.02725483_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.06180814_dp

      ELSE IF (element_symbol == "S " .AND. charge == 6) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.00537161_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.00932309_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.00932309_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.02206245_dp

      ELSE IF (element_symbol == "Sb" .AND. charge == 5) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.13922153_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.00993605_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.00993605_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.02351298_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00943232_dp

      ELSE IF (element_symbol == "Sc" .AND. charge == 11) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.08305262_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.03930672_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.03930672_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.09301668_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00481191_dp

      ELSE IF (element_symbol == "Sc" .AND. charge == 3) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.00846271_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.00883730_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.00883730_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.02091287_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00142986_dp

      ELSE IF (element_symbol == "Se" .AND. charge == 6) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.06219585_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.02742836_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.02742836_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.06490734_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00578395_dp

      ELSE IF (element_symbol == "Si" .AND. charge == 4) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.00037280_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.00610073_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.00610073_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.01443697_dp

      ELSE IF (element_symbol == "Sm" .AND. charge == 16) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.24029994_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.50745883_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.50745883_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -0.31336038_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -0.31336038_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 1.20086678_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 0.74154601_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 0.74154601_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -1.05404138_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.06335216_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.01192432_dp

      ELSE IF (element_symbol == "Sn" .AND. charge == 4) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.10393057_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.00213712_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.00213712_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.00505735_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00706571_dp

      ELSE IF (element_symbol == "Sr" .AND. charge == 10) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.57626460_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.41837761_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.41837761_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.99006212_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00899140_dp

      ELSE IF (element_symbol == "Sr" .AND. charge == 2) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.00402188_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.00938423_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.00938423_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.02220714_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00284555_dp

      ELSE IF (element_symbol == "Ta" .AND. charge == 13) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.64999172_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.14214268_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.14214268_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -0.03012256_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -0.03012256_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.33637097_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 0.07128298_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 0.07128298_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -0.10132239_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.20534356_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.04513357_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.04513357_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.10235332_dp

      ELSE IF (element_symbol == "Ta" .AND. charge == 5) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.48563525_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.29375358_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.29375358_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.69514785_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.08671559_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.04173012_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.04173012_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.09463501_dp

      ELSE IF (element_symbol == "Tb" .AND. charge == 19) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.00872499_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.48343103_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.48343103_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -0.25579926_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -0.25579926_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 1.14400661_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 0.60533153_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 0.60533153_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -0.86042467_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.05175355_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.01619671_dp

      ELSE IF (element_symbol == "Tc" .AND. charge == 15) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.34079130_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.02786799_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.02786799_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.06594771_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01578982_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.00726905_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.00726905_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.01648465_dp

      ELSE IF (element_symbol == "Tc" .AND. charge == 7) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.08075975_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.05944302_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.05944302_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.14066787_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.02602474_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.01842157_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.01842157_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.04177619_dp

      ELSE IF (element_symbol == "Te" .AND. charge == 6) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.17299713_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.02139985_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.02139985_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.05064130_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01080862_dp

      ELSE IF (element_symbol == "Ti" .AND. charge == 12) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.12239544_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.02427591_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.02427591_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.05744728_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00582242_dp

      ELSE IF (element_symbol == "Ti" .AND. charge == 4) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.00933253_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.01068729_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.01068729_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.02529075_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00200953_dp

      ELSE IF (element_symbol == "Tl" .AND. charge == 13) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.01223345_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.01338045_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.01338045_dp
         gth_potential%kprj_ppnl(1, 3, 1) = 0.30299168_dp
         gth_potential%kprj_ppnl(3, 1, 1) = 0.30299168_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.03166391_dp
         gth_potential%kprj_ppnl(2, 3, 1) = -0.71700919_dp
         gth_potential%kprj_ppnl(3, 2, 1) = -0.71700919_dp
         gth_potential%kprj_ppnl(3, 3, 1) = 1.01916448_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.18684918_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.07524992_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.07524992_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.17065079_dp

      ELSE IF (element_symbol == "Tl" .AND. charge == 3) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.16864149_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.00188436_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.00188436_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.00445921_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.02294091_dp

      ELSE IF (element_symbol == "Tm" .AND. charge == 23) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 2.04536983_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.99188883_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.99188883_dp
         gth_potential%kprj_ppnl(1, 3, 1) = 0.56567230_dp
         gth_potential%kprj_ppnl(3, 1, 1) = 0.56567230_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -2.34723738_dp
         gth_potential%kprj_ppnl(2, 3, 1) = -1.33862499_dp
         gth_potential%kprj_ppnl(3, 2, 1) = -1.33862499_dp
         gth_potential%kprj_ppnl(3, 3, 1) = 1.90273579_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.11929895_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.02383347_dp

      ELSE IF (element_symbol == "V " .AND. charge == 13) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.15640785_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.00339312_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.00339312_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.00802959_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00654774_dp

      ELSE IF (element_symbol == "V " .AND. charge == 5) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.01097329_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.01302199_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.01302199_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.03081565_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00271707_dp

      ELSE IF (element_symbol == "W " .AND. charge == 14) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.71547973_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.16283543_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.16283543_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -0.02291919_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -0.02291919_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.38533896_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 0.05423670_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 0.05423670_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -0.07709262_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.20586005_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.04411598_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.04411598_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.10004563_dp

      ELSE IF (element_symbol == "W " .AND. charge == 6) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.50979956_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.31766785_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.31766785_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.75173934_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.07577234_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.03800751_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.03800751_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.08619292_dp

      ELSE IF (element_symbol == "Xe" .AND. charge == 8) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.23660318_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.04583840_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.04583840_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.10847345_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01336223_dp

      ELSE IF (element_symbol == "Y " .AND. charge == 11) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.68452714_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.49749106_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.49749106_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 1.17727871_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01045632_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.00135148_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.00135148_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.00306487_dp

      ELSE IF (element_symbol == "Y " .AND. charge == 3) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.04333632_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.03380154_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.03380154_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.07998905_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00919822_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.00514005_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.00514005_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.01165653_dp

      ELSE IF (element_symbol == "Yb" .AND. charge == 24) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.34946961_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.87654953_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.87654953_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -0.53937986_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -0.53937986_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 2.07429479_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 1.27640570_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 1.27640570_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -1.81429663_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.07007567_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.02571848_dp

      ELSE IF (element_symbol == "Zn" .AND. charge == 12) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.02575928_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.03841864_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.03841864_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.09091511_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01276688_dp

      ELSE IF (element_symbol == "Zn" .AND. charge == 2) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.01213948_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.00206035_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.00206035_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.00487568_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00022546_dp

      ELSE IF (element_symbol == "Zn" .AND. charge == 20) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.06195779_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.19739272_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.19739272_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.23355822_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01770736_dp

      ELSE IF (element_symbol == "Zr" .AND. charge == 12) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.51261281_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.44782121_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.44782121_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 1.05973841_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00756230_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.00237820_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.00237820_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.00539325_dp

      ELSE IF (element_symbol == "Zr" .AND. charge == 4) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.05748649_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.04415700_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.04415700_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.10449452_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00938027_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.00527960_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.00527960_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.01197300_dp
      END IF

      CALL timestop(handle)

   END SUBROUTINE get_kprj_SOC_parameters_LDA

! **************************************************************************************************
!> \brief ...
!> \param gth_potential ...
!> \param element_symbol ...
!> \param charge ...
! **************************************************************************************************
   SUBROUTINE get_kprj_SOC_parameters_PBE(gth_potential, element_symbol, charge)

      TYPE(gth_potential_type), POINTER                  :: gth_potential
      CHARACTER(LEN=2)                                   :: element_symbol
      INTEGER                                            :: charge

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_kprj_SOC_parameters_PBE'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (element_symbol == "Ag" .AND. charge == 11) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.25749375_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.41140313_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.41140313_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.48677875_dp
         gth_potential%kprj_ppnl(1, 1, 2) = -0.04679548_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.08682477_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.08682477_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.09845003_dp

      ELSE IF (element_symbol == "Ag" .AND. charge == 19) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.24535051_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.64883901_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.64883901_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.76771666_dp
         gth_potential%kprj_ppnl(1, 1, 2) = -0.03347751_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.06201267_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.06201267_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.07031576_dp

      ELSE IF (element_symbol == "Al" .AND. charge == 3) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.01323750_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.01265070_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.01265070_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.01496851_dp

      ELSE IF (element_symbol == "Ar" .AND. charge == 8) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.05709588_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.03927810_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.03927810_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.04647447_dp

      ELSE IF (element_symbol == "As" .AND. charge == 5) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.07884496_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.03308528_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.03308528_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.03914703_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00335570_dp

      ELSE IF (element_symbol == "At" .AND. charge == 7) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.66100264_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.11918597_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.11918597_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.14102275_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.03023997_dp

      ELSE IF (element_symbol == "Au" .AND. charge == 11) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 1.21911104_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -2.38187511_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -2.38187511_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 2.81827264_dp
         gth_potential%kprj_ppnl(1, 1, 2) = -0.12113041_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.26246301_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.26246301_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.29760507_dp

      ELSE IF (element_symbol == "Au" .AND. charge == 19) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.18370022_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 2.64213942_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 2.64213942_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -2.01606276_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -2.01606276_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -6.48265855_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 4.77087526_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 4.77087526_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -3.39068639_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.15609346_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.06062427_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.06062427_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.06874145_dp

      ELSE IF (element_symbol == "B " .AND. charge == 3) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.00085498_dp

      ELSE IF (element_symbol == "Ba" .AND. charge == 10) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 1.38677712_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -2.22926159_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -2.22926159_dp
         gth_potential%kprj_ppnl(1, 3, 1) = 0.55122185_dp
         gth_potential%kprj_ppnl(3, 1, 1) = 0.55122185_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 3.55539820_dp
         gth_potential%kprj_ppnl(2, 3, 1) = -1.30442898_dp
         gth_potential%kprj_ppnl(3, 2, 1) = -1.30442898_dp
         gth_potential%kprj_ppnl(3, 3, 1) = 0.92706460_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01663250_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.00607912_dp

      ELSE IF (element_symbol == "Bi" .AND. charge == 15) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 1.61410367_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -2.26162040_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -2.26162040_dp
         gth_potential%kprj_ppnl(1, 3, 1) = 1.01570961_dp
         gth_potential%kprj_ppnl(3, 1, 1) = 1.01570961_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 4.36698694_dp
         gth_potential%kprj_ppnl(2, 3, 1) = -2.40360764_dp
         gth_potential%kprj_ppnl(3, 2, 1) = -2.40360764_dp
         gth_potential%kprj_ppnl(3, 3, 1) = 1.70825672_dp
         gth_potential%kprj_ppnl(1, 1, 2) = -0.17397466_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.35675867_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.35675867_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.40452630_dp

      ELSE IF (element_symbol == "Bi" .AND. charge == 5) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.33910609_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.02108594_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.02108594_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.02494923_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.02035678_dp

      ELSE IF (element_symbol == "Br" .AND. charge == 7) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.17300601_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.09956873_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.09956873_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.11781131_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00654644_dp

      ELSE IF (element_symbol == "C " .AND. charge == 4) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.00207319_dp

      ELSE IF (element_symbol == "Ca" .AND. charge == 10) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.15992263_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.19967894_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.19967894_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.23626331_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00046474_dp

      ELSE IF (element_symbol == "Cd" .AND. charge == 12) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.28931079_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.45543886_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.45543886_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.53888253_dp
         gth_potential%kprj_ppnl(1, 1, 2) = -0.06128372_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.11410996_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.11410996_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.12938854_dp

      ELSE IF (element_symbol == "Ce" .AND. charge == 12) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.53180154_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.20234053_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.20234053_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -0.09655667_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -0.09655667_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.07866042_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 0.22849478_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 0.22849478_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -0.16239245_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01651587_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.00785811_dp

      ELSE IF (element_symbol == "Ce" .AND. charge == 30) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 6.58733940_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -7.76347086_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -7.76347086_dp
         gth_potential%kprj_ppnl(1, 3, 1) = 0.01837566_dp
         gth_potential%kprj_ppnl(3, 1, 1) = 0.01837566_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 9.21645527_dp
         gth_potential%kprj_ppnl(2, 3, 1) = -0.04348474_dp
         gth_potential%kprj_ppnl(3, 2, 1) = -0.04348474_dp
         gth_potential%kprj_ppnl(3, 3, 1) = 0.03090484_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.36360640_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.01337882_dp

      ELSE IF (element_symbol == "Cl" .AND. charge == 7) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.04528923_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.03257832_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.03257832_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.03854719_dp

      ELSE IF (element_symbol == "Co" .AND. charge == 17) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.03467390_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.13534705_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.13534705_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.16014479_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01143107_dp

      ELSE IF (element_symbol == "Cr" .AND. charge == 14) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.19042071_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.04824721_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.04824721_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.05708687_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00894423_dp

      ELSE IF (element_symbol == "Cs" .AND. charge == 9) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 1.47901467_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -4.82006821_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -4.82006821_dp
         gth_potential%kprj_ppnl(1, 3, 1) = 0.71637516_dp
         gth_potential%kprj_ppnl(3, 1, 1) = 0.71637516_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 6.89583701_dp
         gth_potential%kprj_ppnl(2, 3, 1) = -1.69525304_dp
         gth_potential%kprj_ppnl(3, 2, 1) = -1.69525304_dp
         gth_potential%kprj_ppnl(3, 3, 1) = 1.20482534_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01060475_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.00322244_dp

      ELSE IF (element_symbol == "Cu" .AND. charge == 11) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.09147262_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.11244872_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.11244872_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.13305111_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01046090_dp

      ELSE IF (element_symbol == "Cu" .AND. charge == 19) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.14035346_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.18322786_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.18322786_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.21679813_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01478075_dp

      ELSE IF (element_symbol == "Dy" .AND. charge == 20) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.44352171_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 3.08608168_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 3.08608168_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -3.44531996_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -3.44531996_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -9.38743338_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 8.15311509_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 8.15311509_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -5.79446221_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.05555208_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.01775081_dp

      ELSE IF (element_symbol == "Er" .AND. charge == 22) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 2.03813991_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -4.27357152_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -4.27357152_dp
         gth_potential%kprj_ppnl(1, 3, 1) = 5.04555798_dp
         gth_potential%kprj_ppnl(3, 1, 1) = 5.04555798_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 13.45664255_dp
         gth_potential%kprj_ppnl(2, 3, 1) = -11.93996941_dp
         gth_potential%kprj_ppnl(3, 2, 1) = -11.93996941_dp
         gth_potential%kprj_ppnl(3, 3, 1) = 8.48579969_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.07930876_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.02173132_dp

      ELSE IF (element_symbol == "Eu" .AND. charge == 17) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 1.65792601_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -3.71639117_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -3.71639117_dp
         gth_potential%kprj_ppnl(1, 3, 1) = 4.37001603_dp
         gth_potential%kprj_ppnl(3, 1, 1) = 4.37001603_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 11.67270356_dp
         gth_potential%kprj_ppnl(2, 3, 1) = -10.34134539_dp
         gth_potential%kprj_ppnl(3, 2, 1) = -10.34134539_dp
         gth_potential%kprj_ppnl(3, 3, 1) = 7.34964910_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.05828324_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.01315523_dp

      ELSE IF (element_symbol == "F " .AND. charge == 7) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.01261401_dp

      ELSE IF (element_symbol == "Fe" .AND. charge == 16) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.09786820_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.08070002_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.08070002_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.09548555_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01036288_dp

      ELSE IF (element_symbol == "Ga" .AND. charge == 13) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.07922219_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.05277747_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.05277747_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.06244715_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01545711_dp

      ELSE IF (element_symbol == "Ga" .AND. charge == 3) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.04378597_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.01162286_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.01162286_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.01375235_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00090534_dp

      ELSE IF (element_symbol == "Gd" .AND. charge == 18) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.54157080_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.12917284_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.12917284_dp
         gth_potential%kprj_ppnl(1, 3, 1) = 0.16099098_dp
         gth_potential%kprj_ppnl(3, 1, 1) = 0.16099098_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.42086479_dp
         gth_potential%kprj_ppnl(2, 3, 1) = -0.38097419_dp
         gth_potential%kprj_ppnl(3, 2, 1) = -0.38097419_dp
         gth_potential%kprj_ppnl(3, 3, 1) = 0.27076038_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.05014620_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.01468660_dp

      ELSE IF (element_symbol == "Ge" .AND. charge == 4) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.07515712_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.03426241_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.03426241_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.04053983_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00209896_dp

      ELSE IF (element_symbol == "Hf" .AND. charge == 12) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.25433226_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 2.89231345_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 2.89231345_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -3.44888647_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -3.44888647_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -9.16410141_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 8.16155501_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 8.16155501_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -5.80046050_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.07470137_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.01580845_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.01580845_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.01792510_dp

      ELSE IF (element_symbol == "Hg" .AND. charge == 12) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.97627823_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -1.29430950_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -1.29430950_dp
         gth_potential%kprj_ppnl(1, 3, 1) = 0.55490433_dp
         gth_potential%kprj_ppnl(3, 1, 1) = 0.55490433_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 2.45527872_dp
         gth_potential%kprj_ppnl(2, 3, 1) = -1.31314330_dp
         gth_potential%kprj_ppnl(3, 2, 1) = -1.31314330_dp
         gth_potential%kprj_ppnl(3, 3, 1) = 0.93325792_dp
         gth_potential%kprj_ppnl(1, 1, 2) = -0.09794385_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.26202008_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.26202008_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.29710284_dp

      ELSE IF (element_symbol == "Ho" .AND. charge == 21) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.21869071_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 1.32948068_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 1.32948068_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -1.61141491_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -1.61141491_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -4.25582283_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 3.81330367_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 3.81330367_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -2.71013518_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.06631442_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.01966354_dp

      ELSE IF (element_symbol == "I " .AND. charge == 7) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.27569304_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.10630683_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.10630683_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.12578393_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01290065_dp

      ELSE IF (element_symbol == "In" .AND. charge == 13) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.34177179_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.54199431_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.54199431_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.64129631_dp
         gth_potential%kprj_ppnl(1, 1, 2) = -0.07755078_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.13327922_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.13327922_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.15112443_dp

      ELSE IF (element_symbol == "In" .AND. charge == 3) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.10576704_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.01722946_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.01722946_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.02038617_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00226126_dp

      ELSE IF (element_symbol == "Ir" .AND. charge == 17) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.31484605_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 3.24759365_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 3.24759365_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -3.24130779_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -3.24130779_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -9.23888776_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 7.67033419_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 7.67033419_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -5.45134725_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.10055806_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.00337529_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.00337529_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.00382722_dp

      ELSE IF (element_symbol == "Ir" .AND. charge == 9) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.67464477_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -1.13818413_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -1.13818413_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 1.34671763_dp
         gth_potential%kprj_ppnl(1, 1, 2) = -0.15231734_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.28618874_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.28618874_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.32450753_dp

      ELSE IF (element_symbol == "K " .AND. charge == 9) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.08936707_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.06986008_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.06986008_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.08265956_dp

      ELSE IF (element_symbol == "Kr" .AND. charge == 8) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.19991623_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.11006240_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.11006240_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.13022759_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00829734_dp

      ELSE IF (element_symbol == "La" .AND. charge == 11) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.64154939_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.31941115_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.31941115_dp
         gth_potential%kprj_ppnl(1, 3, 1) = 0.08447724_dp
         gth_potential%kprj_ppnl(3, 1, 1) = 0.08447724_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.51857409_dp
         gth_potential%kprj_ppnl(2, 3, 1) = -0.19990965_dp
         gth_potential%kprj_ppnl(3, 2, 1) = -0.19990965_dp
         gth_potential%kprj_ppnl(3, 3, 1) = 0.14207685_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.02046237_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.00738047_dp

      ELSE IF (element_symbol == "Lu" .AND. charge == 25) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -1.45526864_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 6.93159024_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 6.93159024_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -7.42062918_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -7.42062918_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -20.55578405_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 17.56041370_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 17.56041370_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -12.48027930_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.06950855_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.02740377_dp

      ELSE IF (element_symbol == "Mg" .AND. charge == 10) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.10308305_dp

      ELSE IF (element_symbol == "Mg" .AND. charge == 2) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.00414275_dp

      ELSE IF (element_symbol == "Mn" .AND. charge == 15) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.18314922_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.04673645_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.04673645_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.05529931_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00951135_dp

      ELSE IF (element_symbol == "Mo" .AND. charge == 14) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.10103602_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.15762775_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.15762775_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.18650767_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00136004_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.02273862_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.02273862_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.02578317_dp

      ELSE IF (element_symbol == "N " .AND. charge == 5) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.00414750_dp

      ELSE IF (element_symbol == "Na" .AND. charge == 9) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.05324287_dp

      ELSE IF (element_symbol == "Nb" .AND. charge == 13) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.18261561_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.04502076_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.04502076_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.05326929_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00621330_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.01529808_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.01529808_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.01734639_dp

      ELSE IF (element_symbol == "Nd" .AND. charge == 14) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.59294947_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.61371631_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.61371631_dp
         gth_potential%kprj_ppnl(1, 3, 1) = 0.67792729_dp
         gth_potential%kprj_ppnl(3, 1, 1) = 0.67792729_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 1.85480449_dp
         gth_potential%kprj_ppnl(2, 3, 1) = -1.60426878_dp
         gth_potential%kprj_ppnl(3, 2, 1) = -1.60426878_dp
         gth_potential%kprj_ppnl(3, 3, 1) = 1.14016235_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.41574225_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.00947171_dp

      ELSE IF (element_symbol == "Ne" .AND. charge == 8) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.01296979_dp

      ELSE IF (element_symbol == "Ni" .AND. charge == 18) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.13310966_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.15168114_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.15168114_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.17947154_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01349066_dp

      ELSE IF (element_symbol == "O " .AND. charge == 6) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.00755927_dp

      ELSE IF (element_symbol == "Os" .AND. charge == 16) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.54785930_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 2.46435737_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 2.46435737_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -2.55981011_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -2.55981011_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -7.17756040_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 6.05761634_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 6.05761634_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -4.30518011_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.08081206_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.01626524_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.01626524_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.01844305_dp

      ELSE IF (element_symbol == "Os" .AND. charge == 8) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.65252379_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -1.16198770_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -1.16198770_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 1.37488239_dp
         gth_potential%kprj_ppnl(1, 1, 2) = -0.12910182_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.26738080_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.26738080_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.30318133_dp

      ELSE IF (element_symbol == "P " .AND. charge == 5) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.02606888_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.02118456_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.02118456_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.02506591_dp

      ELSE IF (element_symbol == "Pb" .AND. charge == 14) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 1.61384844_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -2.17491530_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -2.17491530_dp
         gth_potential%kprj_ppnl(1, 3, 1) = 1.71513508_dp
         gth_potential%kprj_ppnl(3, 1, 1) = 1.71513508_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 5.42883283_dp
         gth_potential%kprj_ppnl(2, 3, 1) = -4.05875039_dp
         gth_potential%kprj_ppnl(3, 2, 1) = -4.05875039_dp
         gth_potential%kprj_ppnl(3, 3, 1) = 2.88457546_dp
         gth_potential%kprj_ppnl(1, 1, 2) = -0.14916233_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.31992302_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.31992302_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.36275861_dp

      ELSE IF (element_symbol == "Pb" .AND. charge == 4) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.27466972_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.02305844_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.02305844_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.02728312_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01870296_dp

      ELSE IF (element_symbol == "Pd" .AND. charge == 10) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.20640805_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.35787163_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.35787163_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.42343942_dp
         gth_potential%kprj_ppnl(1, 1, 2) = -0.06113102_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.10023162_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.10023162_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.11365198_dp

      ELSE IF (element_symbol == "Pd" .AND. charge == 18) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.44546453_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.88057476_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.88057476_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -1.04191011_dp
         gth_potential%kprj_ppnl(1, 1, 2) = -0.03119207_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.06201267_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.06201267_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.07031576_dp

      ELSE IF (element_symbol == "Pm" .AND. charge == 15) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.21858858_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 1.91346180_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 1.91346180_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -2.21408582_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -2.21408582_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -5.95015377_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 5.23948335_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 5.23948335_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -3.72372865_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.05970555_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.01053306_dp

      ELSE IF (element_symbol == "Po" .AND. charge == 6) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.23215645_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.02511538_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.02511538_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.02971692_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.02706723_dp

      ELSE IF (element_symbol == "Pr" .AND. charge == 13) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.11738149_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.91802470_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.91802470_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -1.36977830_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -1.36977830_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -3.36669344_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 3.24148708_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 3.24148708_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -2.30374209_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01209488_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.00863072_dp

      ELSE IF (element_symbol == "Pt" .AND. charge == 10) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.68131890_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -1.13737593_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -1.13737593_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 1.34576135_dp
         gth_potential%kprj_ppnl(1, 1, 2) = -0.13179902_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.26337958_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.26337958_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.29864437_dp

      ELSE IF (element_symbol == "Pt" .AND. charge == 18) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.35156689_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 2.38906445_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 2.38906445_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -2.06064115_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -2.06064115_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -6.25743244_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 4.87636699_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 4.87636699_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -3.46565992_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.18032724_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.02920945_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.02920945_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.03312040_dp

      ELSE IF (element_symbol == "Rb" .AND. charge == 9) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 1.41180234_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -2.61453234_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -2.61453234_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 3.09355639_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00696032_dp

      ELSE IF (element_symbol == "Re" .AND. charge == 15) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.48781828_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 2.51603934_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 2.51603934_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -2.67688069_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -2.67688069_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -7.43361600_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 6.33465589_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 6.33465589_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -4.50207359_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.07649146_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.01299626_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.01299626_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.01473637_dp

      ELSE IF (element_symbol == "Re" .AND. charge == 7) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.69313998_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -1.31294961_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -1.31294961_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 1.55350293_dp
         gth_potential%kprj_ppnl(1, 1, 2) = -0.11238482_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.25543977_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.25543977_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.28964148_dp

      ELSE IF (element_symbol == "Rh" .AND. charge == 17) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.18648086_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.14716083_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.14716083_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.17412304_dp
         gth_potential%kprj_ppnl(1, 1, 2) = -0.02161325_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.05073757_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.05073757_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.05753100_dp

      ELSE IF (element_symbol == "Rh" .AND. charge == 9) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.18034883_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.29258721_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.29258721_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.34619385_dp
         gth_potential%kprj_ppnl(1, 1, 2) = -0.03265865_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.10951165_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.10951165_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.12417454_dp

      ELSE IF (element_symbol == "Rn" .AND. charge == 8) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.64375030_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.05907890_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.05907890_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.06990310_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.03732898_dp

      ELSE IF (element_symbol == "Ru" .AND. charge == 16) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.15600415_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.48648402_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.48648402_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.57561565_dp
         gth_potential%kprj_ppnl(1, 1, 2) = -0.01536873_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.04203938_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.04203938_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.04766818_dp

      ELSE IF (element_symbol == "Ru" .AND. charge == 8) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.24365860_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.44245368_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.44245368_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.52351826_dp
         gth_potential%kprj_ppnl(1, 1, 2) = -0.02436658_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.09258568_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.09258568_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.10498229_dp

      ELSE IF (element_symbol == "S " .AND. charge == 6) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.03492359_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.02660790_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.02660790_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.03148290_dp

      ELSE IF (element_symbol == "Sb" .AND. charge == 5) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.23747176_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.10907432_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.10907432_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.12905847_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00696671_dp

      ELSE IF (element_symbol == "Sc" .AND. charge == 11) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.20196646_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.05799329_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.05799329_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.06861858_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00480094_dp

      ELSE IF (element_symbol == "Se" .AND. charge == 6) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.15386757_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.09508784_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.09508784_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.11250945_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00489231_dp

      ELSE IF (element_symbol == "Si" .AND. charge == 4) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.01877687_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.01639758_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.01639758_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.01940188_dp

      ELSE IF (element_symbol == "Sm" .AND. charge == 16) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.68601615_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 3.27621519_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 3.27621519_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -3.50012699_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -3.50012699_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -9.70364773_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 8.28281221_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 8.28281221_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -5.88663864_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.04985656_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.01177755_dp

      ELSE IF (element_symbol == "Sn" .AND. charge == 4) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.22160101_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.11203353_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.11203353_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.13255986_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00451584_dp

      ELSE IF (element_symbol == "Sr" .AND. charge == 10) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 1.34629459_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -2.37169277_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -2.37169277_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 2.80622473_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00779193_dp

      ELSE IF (element_symbol == "Ta" .AND. charge == 13) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 1.12324816_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.67250512_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.67250512_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -0.97730022_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -0.97730022_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -2.42277460_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 2.31271443_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 2.31271443_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -1.64365843_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.07060129_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.00876170_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.00876170_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.00993484_dp

      ELSE IF (element_symbol == "Ta" .AND. charge == 5) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.67406021_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -1.28529463_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -1.28529463_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 1.52078112_dp
         gth_potential%kprj_ppnl(1, 1, 2) = -0.08463039_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.21860129_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.21860129_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.24787056_dp

      ELSE IF (element_symbol == "Tb" .AND. charge == 19) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.85647669_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 4.27218592_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 4.27218592_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -4.86408800_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -4.86408800_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -13.15288326_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 11.51053306_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 11.51053306_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -8.18059700_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.05530956_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.01609128_dp

      ELSE IF (element_symbol == "Tc" .AND. charge == 15) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.00672609_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.27805622_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.27805622_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.32900056_dp
         gth_potential%kprj_ppnl(1, 1, 2) = -0.00760164_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.03277719_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.03277719_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.03716584_dp

      ELSE IF (element_symbol == "Te" .AND. charge == 6) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.28709255_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.13816104_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.13816104_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.16347435_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00929869_dp

      ELSE IF (element_symbol == "Ti" .AND. charge == 12) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.20121752_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.06088315_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.06088315_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.07203791_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00584456_dp

      ELSE IF (element_symbol == "Tl" .AND. charge == 13) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 1.35297827_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -2.28901237_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -2.28901237_dp
         gth_potential%kprj_ppnl(1, 3, 1) = 0.68067754_dp
         gth_potential%kprj_ppnl(3, 1, 1) = 0.68067754_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 3.84162025_dp
         gth_potential%kprj_ppnl(2, 3, 1) = -1.61077704_dp
         gth_potential%kprj_ppnl(3, 2, 1) = -1.61077704_dp
         gth_potential%kprj_ppnl(3, 3, 1) = 1.14478780_dp
         gth_potential%kprj_ppnl(1, 1, 2) = -0.11559254_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.30671016_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.30671016_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.34777664_dp

      ELSE IF (element_symbol == "Tl" .AND. charge == 3) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.18274699_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.00921547_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.00921547_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.01090389_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01256903_dp

      ELSE IF (element_symbol == "Tm" .AND. charge == 23) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 2.55564771_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -6.24092488_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -6.24092488_dp
         gth_potential%kprj_ppnl(1, 3, 1) = 7.42426380_dp
         gth_potential%kprj_ppnl(3, 1, 1) = 7.42426380_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 19.74462887_dp
         gth_potential%kprj_ppnl(2, 3, 1) = -17.56901479_dp
         gth_potential%kprj_ppnl(3, 2, 1) = -17.56901479_dp
         gth_potential%kprj_ppnl(3, 3, 1) = 12.48639214_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.11677823_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.02372319_dp

      ELSE IF (element_symbol == "V " .AND. charge == 13) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.19032152_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.05946067_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.05946067_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.07035481_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.00648186_dp

      ELSE IF (element_symbol == "W " .AND. charge == 14) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 1.42485618_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.04791036_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.04791036_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -0.23348205_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -0.23348205_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.33202370_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 0.55251937_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 0.55251937_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -0.39267845_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.07015958_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.01318853_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.01318853_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.01495438_dp

      ELSE IF (element_symbol == "W " .AND. charge == 6) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.70941707_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -1.36284306_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -1.36284306_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 1.61253765_dp
         gth_potential%kprj_ppnl(1, 1, 2) = -0.10122950_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.22114149_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.22114149_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.25075088_dp

      ELSE IF (element_symbol == "Xe" .AND. charge == 8) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.27633099_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.09115677_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.09115677_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.10785815_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01674755_dp

      ELSE IF (element_symbol == "Y " .AND. charge == 11) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 1.23578364_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -2.38648117_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -2.38648117_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 2.82372260_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01304100_dp
         gth_potential%kprj_ppnl(1, 2, 2) = -0.00010639_dp
         gth_potential%kprj_ppnl(2, 1, 2) = -0.00010639_dp
         gth_potential%kprj_ppnl(2, 2, 2) = 0.00012063_dp

      ELSE IF (element_symbol == "Yb" .AND. charge == 24) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.03061526_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 2.18978914_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 2.18978914_dp
         gth_potential%kprj_ppnl(1, 3, 1) = -2.38963165_dp
         gth_potential%kprj_ppnl(3, 1, 1) = -2.38963165_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -6.56936572_dp
         gth_potential%kprj_ppnl(2, 3, 1) = 5.65490060_dp
         gth_potential%kprj_ppnl(3, 2, 1) = 5.65490060_dp
         gth_potential%kprj_ppnl(3, 3, 1) = -4.01896790_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.06796310_dp
         gth_potential%kprj_ppnl(1, 1, 3) = 0.02565096_dp

      ELSE IF (element_symbol == "Zn" .AND. charge == 12) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 0.10017431_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -0.10693277_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -0.10693277_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 0.12652456_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01272732_dp

      ELSE IF (element_symbol == "Zn" .AND. charge == 20) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = -0.00982281_dp
         gth_potential%kprj_ppnl(1, 2, 1) = 0.28902245_dp
         gth_potential%kprj_ppnl(2, 1, 1) = 0.28902245_dp
         gth_potential%kprj_ppnl(2, 2, 1) = -0.34197598_dp
         gth_potential%kprj_ppnl(1, 1, 2) = 0.01730011_dp

      ELSE IF (element_symbol == "Zr" .AND. charge == 12) THEN
         gth_potential%kprj_ppnl(1, 1, 1) = 1.00490626_dp
         gth_potential%kprj_ppnl(1, 2, 1) = -1.41611759_dp
         gth_potential%kprj_ppnl(2, 1, 1) = -1.41611759_dp
         gth_potential%kprj_ppnl(2, 2, 1) = 1.67557293_dp
         gth_potential%kprj_ppnl(1, 1, 2) = -0.00255322_dp
         gth_potential%kprj_ppnl(1, 2, 2) = 0.01079629_dp
         gth_potential%kprj_ppnl(2, 1, 2) = 0.01079629_dp
         gth_potential%kprj_ppnl(2, 2, 2) = -0.01224184_dp
      END IF

      CALL timestop(handle)

   END SUBROUTINE get_kprj_SOC_parameters_PBE

END MODULE external_potential_soc_parameters
