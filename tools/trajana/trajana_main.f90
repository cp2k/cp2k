!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM trajana
   USE trajana_command_line,            ONLY: command_name,&
                                              fail,&
                                              has_flag
   USE trajana_dsf_analysis,            ONLY: run_dsf
   USE trajana_geometry_analysis,       ONLY: run_geometry
   USE trajana_hbond_analysis,          ONLY: run_hbond
   USE trajana_rdf_analysis,            ONLY: run_rdf
   USE trajana_text_utils,              ONLY: lower_case
   USE trajana_vacf_analysis,           ONLY: run_vacf

   IMPLICIT NONE

   CHARACTER(LEN=:), ALLOCATABLE :: analysis

   analysis = lower_case(command_name())
   IF (has_flag("--help")) analysis = "help"
   SELECT CASE (analysis)
   CASE ("", "help", "--help", "-h")
      CALL print_help()
   CASE ("geometry")
      CALL run_geometry()
   CASE ("rdf")
      CALL run_rdf()
   CASE ("vacf")
      CALL run_vacf()
   CASE ("hbond")
      CALL run_hbond()
   CASE ("dsf")
      CALL run_dsf()
   CASE DEFAULT
      CALL fail("Unknown trajana analysis: "//analysis)
   END SELECT

CONTAINS

   SUBROUTINE print_help()
      WRITE (*, "(A)") "Usage: trajana.x ANALYSIS [options]"
      WRITE (*, "(A)") ""
      WRITE (*, "(A)") "Analyses:"
      WRITE (*, "(A)") "  geometry  Distances, angles, and torsions from an action file"
      WRITE (*, "(A)") "  rdf       Radial distribution and running coordination number"
      WRITE (*, "(A)") "  vacf      Velocity autocorrelation and optional spectrum"
      WRITE (*, "(A)") "  hbond     Wernet-type hydrogen-bond network statistics"
      WRITE (*, "(A)") "  dsf       Coherent/self scattering and collective current spectra"
      WRITE (*, "(A)") ""
      WRITE (*, "(A)") "Run the tool from tools/trajana and see README.md for complete examples."
   END SUBROUTINE print_help

END PROGRAM trajana
