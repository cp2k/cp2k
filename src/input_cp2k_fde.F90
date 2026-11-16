!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE input_cp2k_fde
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_create,&
                                              section_type
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_fde'

   PUBLIC :: create_fde_section

CONTAINS

! **************************************************************************************************
!> \brief Create the FDE (frozen-density embedding) input section
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_fde_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="FDE", &
                          description="Parameters for in-DFT frozen-density embedding (FDE). "// &
                          "CP2K drives a freeze-and-thaw optimization and spawns another "// &
                          "program to perform another calculation on the subsystem. The two "// &
                          "codes interact via the embedding potential and densities exchanged "// &
                          "on the real-space grid. The nonadditive kinetic-energy functional "// &
                          "is taken from the KG_METHOD section.", &
                          n_keywords=13, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="_SECTION_PARAMETERS_", &
                          description="Controls the activation of the FDE driver.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FREEZE_AND_THAW", &
                          description="Run the freeze-and-thaw loop, relaxing both the "// &
                          "environment and the subsystem alternately. If .FALSE., perform "// &
                          "a one-shot embedding of the subsystem in the fixed environment "// &
                          "density.", &
                          usage="FREEZE_AND_THAW T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_FT_ITER", &
                          description="Maximum number of freeze-and-thaw iterations.", &
                          usage="MAX_FT_ITER 20", default_i_val=50)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_FT_DENSITY", &
                          description="Convergence threshold on the absolute change of the "// &
                          "density between freeze-and-thaw iterations.", &
                          usage="EPS_FT_DENSITY 1.0E-5", default_r_val=1.0e-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_FT_ENERGY", &
                          description="Convergence threshold on the total energy change "// &
                          "between freeze-and-thaw iterations.", &
                          usage="EPS_FT_ENERGY 1.0E-6", default_r_val=1.0e-6_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DENSITY_MIXING", &
                          description="Linear mixing factor applied to the subsystem density "// &
                          "between macro-iterations. 1.0 = no damping.", &
                          usage="DENSITY_MIXING 0.5", default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COMMAND", &
                          description="Shell command used to run the subsystem solver on the "// &
                          "generated input file. The input file name is passed as the last "// &
                          "argument.", &
                          usage="COMMAND pymolcas", default_c_val="pymolcas")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PROJECT", &
                          description="Project name for the active-subsystem calculation.", &
                          usage="PROJECT sub", default_c_val="sub")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="WORKDIR", &
                          description="Working directory in which the subsystem calculation is "// &
                          "run and the exchange files are written. Defaults to the current "// &
                          "working directory.", &
                          usage="WORKDIR ./sub", default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MOLCAS_INPUT_TEMPLATE", &
                          description="Path to the OpenMolcas input template for the subsystem "// &
                          "(GATEWAY/SEWARD/RASSCF). The driver appends the EMBEDDING block.", &
                          usage="MOLCAS_INPUT_TEMPLATE sub.template.inp", default_c_val="")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="POT_FILE", &
                          description="File name (in WORKDIR) for the embedding potential "// &
                          "written in OpenMolcas format.", &
                          usage="POT_FILE EMBPOT", default_c_val="EMBPOT")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="GRID_FILE", &
                          description="File name (in WORKDIR) for the output grid coordinates "// &
                          "on which the subsystem solver evaluates and returns the density.", &
                          usage="GRID_FILE EMBGRID", default_c_val="EMBGRID")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RHO_FILE", &
                          description="File name (in WORKDIR) with the subsystem density "// &
                          "returned by the subsystem solver on the output grid.", &
                          usage="RHO_FILE EMBRHOA", default_c_val="EMBRHOA")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RESULT_FILE", &
                          description="File name (in WORKDIR) with the scalar results returned "// &
                          "by the subsystem solver (embedded energy, Tr(D*v_emb), active "// &
                          "nuclear repulsion, number of active electrons).", &
                          usage="RESULT_FILE EMBRES", default_c_val="EMBRES")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_fde_section

END MODULE input_cp2k_fde
