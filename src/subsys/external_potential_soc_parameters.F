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
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'external_potential_soc_parameters'

   PUBLIC :: get_kprj_SOC_parameters

CONTAINS

   SUBROUTINE get_kprj_SOC_parameters(gth_potential, xc_functional)

      TYPE(gth_potential_type), POINTER                  :: gth_potential
      CHARACTER(LEN=*)                                   :: xc_functional

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_kprj_SOC_parameters'

      INTEGER :: handle

      CALL timeset(routineN, handle)

      gth_potential%kprj_ppnl(:,:,:) = 0.0_dp 

      SELECT CASE (xc_functional)
      CASE("LDA")

 !        CALL get_kprj_SOC_parameters_LDA(gth_potential)

      CASE("PBE")

      CASE DEFAULT
        CPABORT("Unknown xc functional for retrieving SOC parameter of pseudopotential.")
      END SELECT


   END SUBROUTINE get_kprj_SOC_parameters

END MODULE external_potential_soc_parameters
