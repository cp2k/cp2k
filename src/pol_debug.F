!-----------------------------------------------------------------------------
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/pol_debug [1.0] *
!!
!!   NAME
!!     pol_debug
!!
!!   FUNCTION
!!     Logic to debug pol_electrostatics and pol_overlap
!!
!!   AUTHOR
!!     CJM 
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE pol_debug

!-----------------------------------------------------------------------------!
  USE atomic_kinds, ONLY : kind_info_type
  USE coefficient_types, ONLY : coeff_type
  USE dg_types, ONLY : dg_type
  USE empirical_parameters, ONLY : empirical_parameter_type
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE linklist_types, ONLY : linklist_internal_data_type
  USE md, ONLY : thermodynamic_type
  USE molecule_types, ONLY : molecule_structure_type, particle_node_type
  USE pair_potential, ONLY : potentialparm_type
  USE particle_types, ONLY : particle_type
  USE pol_force, ONLY : force_env_types, debug_variables_type, get_pol_grids
  USE pol_force_numer, ONLY : part_electrostatics_numer, coef_electrostatics_numer, &
                              part_overlap_numer, coef_overlap_numer
  USE pw_grid_types, ONLY : pw_grid_type
  USE simulation_cell, ONLY : cell_type
  USE termination, ONLY : stop_memory, stop_program
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: debug_control
!-----------------------------------------------------------------------------!
!
    CONTAINS
!
!-----------------------------------------------------------------------------!
SUBROUTINE debug_control ( globenv, ewald_param, pnode, part, &
           molecule, box, box_ref, drho_basis_info, rho0_basis_info, &
           coef_pos, coef_vel, coef_force, thermo, potparm, empparm,  &
           ll_data )
                        
! Routine to execute numerical tests for debugging 
  IMPLICIT NONE
! Arguments
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
  TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  TYPE ( particle_type ), DIMENSION ( : ), INTENT ( INOUT ) :: part
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( INOUT ) :: pnode
  TYPE ( potentialparm_type ), DIMENSION ( :,: ), INTENT ( IN ) :: potparm
  TYPE ( empirical_parameter_type ), DIMENSION ( : ), INTENT ( IN ) :: empparm
  TYPE ( cell_type ), INTENT ( INOUT ) :: box
  TYPE ( cell_type ), INTENT ( INOUT ) :: box_ref
  TYPE ( coeff_type ), INTENT ( INOUT ) :: coef_pos
  TYPE ( coeff_type ), INTENT ( INOUT ) :: coef_vel
  TYPE ( coeff_type ), INTENT ( INOUT ) :: coef_force
  TYPE ( kind_info_type ), POINTER, DIMENSION ( : ) :: drho_basis_info
  TYPE ( kind_info_type ), POINTER, DIMENSION ( : ) :: rho0_basis_info
  TYPE ( linklist_internal_data_type), DIMENSION ( : ), INTENT ( INOUT ) :: ll_data
  TYPE ( molecule_structure_type ), DIMENSION ( : ), INTENT ( IN ) :: molecule
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo
!
! begin local variable declaration
!
  TYPE ( debug_variables_type ) :: dbg
  TYPE ( dg_type ), DIMENSION ( : ), POINTER :: dg_coef
  TYPE ( dg_type ), DIMENSION ( : ), POINTER :: dg_part
  TYPE ( pw_grid_type ), DIMENSION ( : ), POINTER :: pw_big
  TYPE ( pw_grid_type ), DIMENSION ( : ), POINTER :: pw_small
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: f_part, rel_part, diff_part
  REAL ( dbl ), DIMENSION ( : ), ALLOCATABLE :: rel_coef, diff_coef, f_coef
  INTEGER :: iflag, i, natoms, isos, iatom, ncoeff, icoef, iw, iwf, ir, iref 
  INTEGER, DIMENSION ( 2 ) :: dum
  INTEGER, DIMENSION ( 1 ) :: cdum
  REAL ( dbl ) :: delta, energy_numer
  REAL ( dbl ) :: pv_test ( 3, 3 )
  LOGICAL :: box_change=.FALSE.

  iw = globenv % scr
  iwf = 7
  ir = 5
!
  WRITE ( iw,'(16x,A)' ) &
    '***************************************************'
  WRITE ( iw,'(16x,A)' ) '                   BEGIN POL DEBUG'
  WRITE ( iw,'(16x,A)' ) &
    '***************************************************'
! allocating all arrays involved with debug
  natoms = SIZE ( pnode )
  ncoeff = SIZE ( coef_pos % ao % cr )

  IF ( .NOT. ALLOCATED ( diff_part ) ) ALLOCATE ( diff_part ( 3, natoms ), STAT = isos )
  IF ( .NOT. ALLOCATED ( rel_part ) ) ALLOCATE ( rel_part ( 3, natoms ), STAT = isos )
  IF ( .NOT. ALLOCATED ( f_part ) ) ALLOCATE ( f_part ( 3, natoms ), STAT = isos )
  IF ( .NOT. ALLOCATED ( diff_coef ) ) ALLOCATE ( diff_coef ( ncoeff ), STAT = isos )
  IF ( .NOT. ALLOCATED ( rel_coef ) ) ALLOCATE ( rel_coef ( ncoeff ), STAT = isos )
  IF ( .NOT. ALLOCATED ( f_coef ) ) ALLOCATE ( f_coef ( ncoeff ), STAT = isos )

  ALLOCATE ( dbg % f_elec_part ( 3, natoms ), STAT = isos )
  ALLOCATE ( dbg % f_elec_coef ( ncoeff ), STAT = isos )
  ALLOCATE ( dbg % f_ovlp_part ( 3, natoms ), STAT = isos )
  ALLOCATE ( dbg % f_ovlp_coef ( ncoeff ), STAT = isos )

  IF (isos/=0) CALL stop_program ( 'debug', 'allocation failed' )

  CALL force_env_types ( molecule, pnode, part, box, box_ref, &
     drho_basis_info, rho0_basis_info, coef_pos, coef_vel, coef_force, &
     thermo, potparm, empparm, ewald_param, box_change, globenv, ll_data, dbg)

! get local copies of the polarization grids

  CALL get_pol_grids ( pw_big, pw_small, dg_part, dg_coef, iref ) 
  
!
! Debug electrostatics on ions
!
  WRITE ( iw, '(A)' ) &
    ' DO YOU WANT TO DEBUG ELECTROSTATICS ON IONS (1=yes)?'
  READ ( ir, * ) iflag
  IF ( iflag == 1 ) THEN
! get numerical force for the electrostatic contribution on the ions
    WRITE ( iw, '(A)' ) ' ENTER A DELTA LESS THAN 1'
    READ ( ir, * ) delta
    IF ( delta >= 1._dbl .OR. delta <= 0._dbl ) THEN
      delta = 1.E-5_dbl
      WRITE ( iw, '(A,T71,F10.6)' ) &
        ' DELTA (changed to default) = ', delta
    ELSE
      WRITE ( iw, '(A,T71,F10.6)' ) ' DELTA = ', delta
    END IF
    CALL part_electrostatics_numer ( delta, dg_part, dg_coef, part, &
                                   coef_pos, box, pw_small, pw_big, &
                                   ewald_param, iref, energy_numer, f_part )
    WRITE ( iw, '(A,T61,E20.14)' ) ' ELECTROST NUMER ENER = ', &
            energy_numer
    WRITE ( iw, '(A,T61,E20.14)' ) ' ELECTROST ANAL ENER = ', &
            dbg % pot_electrostatics
! write out all particle numbers and forces
    WRITE ( iw, '( /,A,A )' ) &
          ' DO YOU WANT TO GET THE ELECTROSTATIC FORCES ON ALL PARTICLES ', &
          '(FILE=part_electrostatics, 1=yes)?'
     READ ( ir, * ) iflag
     IF ( iflag == 1 ) THEN
        OPEN (UNIT = iwf, &
              FILE = 'part_electrostatics', STATUS = 'REPLACE', ACTION = 'WRITE')
        DO iatom = 1, natoms
        WRITE ( iwf, '( A,T71,I10 )' ) ' ATOM NUMBER ', iatom
        WRITE ( iwf, '( A,T21,3G20.12 )' ) &
             ' F ANA ELECTROSTATIC =', dbg % f_elec_part ( :, iatom )
        WRITE ( iwf, '( A,T21,3G20.12,/ )' ) &
             ' F NUM ELECTROSTATIC =', f_part ( :, iatom )
        END DO
        CLOSE (iwf)
     END IF
! computing the absolute value of the differences in the forces
    diff_part = ABS ( dbg % f_elec_part - f_part )
    IF ( .NOT. ANY ( dbg % f_elec_part == 0._dbl ) )  &
    rel_part = diff_part / dbg % f_elec_part
! find the maximum difference and the relative and absolute errors.
    WRITE ( iw, '(A,T61,E20.14)' ) ' MAXIMUM ABSOLUTE ERROR = ', &
      MAXVAL ( diff_part )
    IF ( .NOT. ANY ( dbg % f_elec_part == 0._dbl ) )  THEN
      WRITE ( iw,'(A,T61,E20.14)' ) ' MAXIMUM RELATIVE ERROR = ', &
      MAXVAL ( rel_part )
    ENDIF
! write out the particle number and forces of the max absolute and relative error         
    dum = MAXLOC ( diff_part )
    WRITE ( iw, '(A,T71,I10)' ) &
      ' PARTICLE WITH MAX ABSOLUTE ERROR IS ', dum ( 2 )
    WRITE ( iw, '(A,T21,3(2X,G20.14))' ) ' F ANAL ON IONS =', &
      dbg % f_elec_part (:, dum ( 2 ) )
    WRITE ( iw, '(A,T21,3(2X,G20.14))' ) ' F NUMR ON IONS =', &
      f_part (:, dum ( 2 ) )
    IF  ( .NOT. ANY ( dbg % f_elec_part == 0._dbl ) ) THEN
      dum = MAXLOC ( rel_part )
      WRITE ( iw, '(A,T71,I10)' ) &
        ' NON BOND PARTICLE WITH MAX RELATIVE ERROR IS ', dum ( 2 )
      WRITE ( iw, '(A,T21,3G20.14)' ) ' F ANAL ON IONS =', &
        dbg % f_elec_part (:, dum ( 2 ) )
      WRITE ( iw, '(A,T21,3G20.14)' ) ' F NUMR ON IONS =', &
        f_part( :, dum ( 2 ) )
    ENDIF
  END IF
!
! Debug electrostatics on coefficients
! 
  WRITE ( iw, '(A)' ) &
    ' DO YOU WANT TO DEBUG ELECTROSTATICS ON COEFFICENTS (1=yes)?'
  READ ( ir, * ) iflag
  IF ( iflag == 1 ) THEN
! get numerical force for the electrostatic contribution on the Coefficents
    WRITE ( iw, '(A)' ) ' ENTER A DELTA LESS THAN 1'
    READ ( ir, * ) delta
    IF ( delta >= 1 .OR. delta <= 0 ) THEN
      delta = 1.E-5_dbl
      WRITE ( iw, '(A,T71,F10.6)' ) &
        ' DELTA ( changed to default ) = ', delta
    ELSE
      WRITE ( iw, '(A,T71,F10.6)' ) ' DELTA = ', delta
    END IF
    CALL coef_electrostatics_numer ( delta, dg_part, dg_coef, part,  &
                                    coef_pos, box, pw_small, pw_big,  &
                                  ewald_param, iref, energy_numer, f_coef )
! write out all particle numbers and forces
    WRITE ( iw, '( /,A,A )' ) &
          ' DO YOU WANT TO GET THE ELECTROSTATIC FORCES ON A COEFFICIENTS ', &
          '(FILE=coef_electrostatics, 1=yes)?'
     READ ( ir, * ) iflag
     IF ( iflag == 1 ) THEN
        OPEN (UNIT = iwf, &
              FILE = 'coef_electrostatics', STATUS = 'REPLACE', ACTION = 'WRITE')
        DO icoef = 1, ncoeff
        WRITE ( iwf, '( A,T71,I10 )' ) ' COEFF NUMBER ', icoef
        WRITE ( iwf, '( A,T21,3G20.12 )' ) &
             ' F ANA ELECTROSTATIC =', dbg % f_elec_coef ( icoef )
        WRITE ( iwf, '( A,T21,3G20.12,/ )' ) &
             ' F NUM ELECTROSTATIC =', f_coef ( icoef )
        END DO
        CLOSE (iwf)
     END IF
! computing the absolute value of the differences in the forces
    diff_coef = ABS ( dbg % f_elec_coef - f_coef )
    IF ( .NOT. ANY ( dbg % f_elec_coef == 0._dbl ) )  &
    rel_coef = diff_coef / dbg % f_elec_coef
! find the maximum difference and the relative and absolute errors.
    WRITE ( iw, '(A,T61,E20.14)' ) ' MAXIMUM ABSOLUTE ERROR = ', &
      MAXVAL ( diff_coef )
    IF ( .NOT. ANY ( dbg % f_elec_coef == 0._dbl ) ) THEN
      WRITE ( iw,'(A,T61,E20.14)' ) ' MAXIMUM RELATIVE ERROR = ', &
                                    MAXVAL ( rel_coef )
    ENDIF
! write out the particle number and forces of 
! the max absolute and relative error         
    cdum = MAXLOC ( diff_coef )
    WRITE ( iw, '(A,T71,I10)' ) &
      ' COEFFICENT WITH MAX ABSOLUTE ERROR IS ', cdum
    WRITE ( iw, '(A,T21,3G20.14)' ) ' F ANAL ON COEF =', &
      dbg % f_elec_coef ( cdum )
    WRITE ( iw,'(A,T21,3G20.14)' ) ' F NUMR ON COEF =', &
      f_coef ( cdum )
    IF ( .NOT. ANY ( dbg % f_elec_coef == 0._dbl ) ) THEN
      cdum = MAXLOC ( rel_coef )
      WRITE ( iw,'(A,T71,I10)' ) &
        ' COEFFICENT WITH MAX RELATIVE ERROR IS ', cdum   
      WRITE ( iw,'(A,T21,3G20.14)' ) ' F ANAL ON COEF =', &
        dbg % f_elec_coef ( cdum )
      WRITE ( iw,'(A,T21,3G20.14)' ) ' F NUMR ON COEF =', &
        f_coef ( cdum )
    ENDIF
  END IF
!
! Debug overlap forces on ions
!
  WRITE ( iw, '(A)' ) &
    ' DO YOU WANT TO DEBUG OVERLAP FORCES ON IONS (1=yes)?'
  READ ( ir, * ) iflag
  IF ( iflag == 1 ) THEN
! get numerical force for the overlap contribution on the ions
    WRITE ( iw, '(A)' ) ' ENTER A DELTA LESS THAN 1'
    READ ( ir, * ) delta
    IF ( delta >= 1 .OR. delta <= 0 ) THEN
      delta = 1.E-5_dbl
      WRITE ( iw, '(A,T71,F10.6)' ) &
        ' DELTA (changed to default) = ', delta
    ELSE
      WRITE ( iw, '(A,T71,F10.6)' ) ' DELTA = ', delta
    END IF
    CALL part_overlap_numer ( delta, coef_pos, drho_basis_info, part, pnode, &
                              box, empparm, energy_numer, f_part )
    WRITE ( iw, '(A,T61,E20.14)' ) ' OVERLAP NUMER ENER = ', &
            energy_numer
    WRITE ( iw, '(A,T61,E20.14)' ) ' OVERLAP ANAL ENER = ', &
            dbg % pot_overlap
! write out all particle numbers and forces
    WRITE ( iw, '( /,A,A )' ) &
          ' DO YOU WANT TO GET THE OVERLAP FORCES ON ALL PARTICLES ', &
          '(FILE=part_overlap, 1=yes)?'
     READ ( ir, * ) iflag
     IF ( iflag == 1 ) THEN
        OPEN (UNIT = iwf, &
              FILE = 'part_overlap', STATUS = 'REPLACE', ACTION = 'WRITE')
        DO iatom = 1, natoms
        WRITE ( iwf, '( A,T71,I10 )' ) ' ATOM NUMBER ', iatom
        WRITE ( iwf, '( A,T21,3G20.12 )' ) &
             ' F ANA OVERLAP =', dbg % f_ovlp_part ( :, iatom )
        WRITE ( iwf, '( A,T21,3G20.12,/ )' ) &
             ' F NUM OVERLAP =', f_part ( :, iatom )
        END DO
        CLOSE (iwf)
     END IF
! computing the absolute value of the differences in the forces
    diff_part = ABS ( dbg % f_ovlp_part - f_part )
    IF ( .NOT. ANY ( dbg % f_ovlp_part == 0._dbl ) ) &
    rel_part = diff_part / dbg % f_ovlp_part
! find the maximum difference and the relative and absolute errors.
    WRITE ( iw, '(A,T61,E20.14)' ) ' MAXIMUM ABSOLUTE ERROR = ', &
      MAXVAL ( diff_part )
    IF ( .NOT. ANY ( dbg % f_ovlp_part == 0._dbl ) ) THEN
      WRITE ( iw,'(A,T61,E20.14)' ) ' MAXIMUM RELATIVE ERROR = ', &
      MAXVAL ( rel_part )
    ENDIF
! write out the particle number and forces of the max absolute and relative error         
    dum = MAXLOC ( diff_part )
    WRITE ( iw, '(A,T71,I10)' ) &
      ' PARTICLE WITH MAX ABSOLUTE ERROR IS ', dum ( 2 )
    WRITE ( iw, '(A,T21,3(2X,G20.14))' ) ' F ANAL ON IONS =', &
      dbg % f_ovlp_part (:, dum ( 2 ) )
    WRITE ( iw, '(A,T21,3(2X,G20.14))' ) ' F NUMR ON IONS =', &
      f_part (:, dum ( 2 ) )
    IF ( .NOT. ANY ( dbg % f_ovlp_part == 0._dbl ) ) THEN
      dum = MAXLOC ( rel_part )
      WRITE ( iw, '(A,T71,I10)' ) &
        ' NON BOND PARTICLE WITH MAX RELATIVE ERROR IS ', dum ( 2 )
      WRITE ( iw, '(A,T21,3G20.14)' ) ' F ANAL ON IONS =', &
        dbg % f_ovlp_part (:, dum ( 2 ) )
      WRITE ( iw, '(A,T21,3G20.14)' ) ' F NUMR ON IONS =', &
        f_part( :, dum ( 2 ) )
    ENDIF
  END IF
!
! Debug Overlap contribution on coefficients
! 
  WRITE ( iw, '(A)' ) &
    ' DO YOU WANT TO DEBUG OVERLAP ON COEFFICENTS (1=yes)?'
  READ ( ir, * ) iflag
  IF ( iflag == 1 ) THEN
! get numerical force for the overlap contribution on the Coefficents
    WRITE ( iw, '(A)' ) ' ENTER A DELTA LESS THAN 1'
    READ ( ir, * ) delta
    IF ( delta >= 1 .OR. delta <= 0 ) THEN
      delta = 1.E-5_dbl
      WRITE ( iw, '(A,T71,F10.6)' ) &
        ' DELTA ( changed to default ) = ', delta
    ELSE
      WRITE ( iw, '(A,T71,F10.6)' ) ' DELTA = ', delta
    END IF
    CALL coef_overlap_numer ( delta, coef_pos, drho_basis_info, part, pnode, &
                              box, empparm, energy_numer, f_coef )

! write out all particle numbers and forces
    WRITE ( iw, '( /,A,A )' ) &
          ' DO YOU WANT TO GET THE OVERLAP FORCES ON ALL COEFFICIENTS ', &
          '(FILE=coef_overlap, 1=yes)?'
     READ ( ir, * ) iflag
     IF ( iflag == 1 ) THEN
        OPEN (UNIT = iwf, &
              FILE = 'coef_overlap', STATUS = 'REPLACE', ACTION = 'WRITE')
        DO icoef = 1, ncoeff
        WRITE ( iwf, '( A,T71,I10 )' ) ' COEFF NUMBER ', icoef
        WRITE ( iwf, '( A,T21,3G20.12 )' ) &
             ' F ANA OVERLAP =', dbg % f_ovlp_coef ( icoef )
        WRITE ( iwf, '( A,T21,3G20.12,/ )' ) &
             ' F NUM OVERLAP =', f_coef ( icoef )
        END DO
        CLOSE (iwf)
     END IF
! computing the absolute value of the differences in the forces
    diff_coef = ABS ( dbg % f_ovlp_coef - f_coef )
    IF ( .NOT. ANY ( dbg % f_ovlp_coef == 0._dbl ) ) &
    rel_coef = diff_coef / dbg % f_ovlp_coef
! find the maximum difference and the relative and absolute errors.
    WRITE ( iw, '(A,T61,E20.14)' ) ' MAXIMUM ABSOLUTE ERROR = ', &
      MAXVAL ( diff_coef )
    IF ( .NOT. ANY ( dbg % f_ovlp_coef == 0._dbl ) ) THEN
      WRITE ( iw,'(A,T61,E20.14)' ) ' MAXIMUM RELATIVE ERROR = ', &
                                    MAXVAL ( rel_coef )
    ENDIF
! write out the particle number and forces of 
! the max absolute and relative error         
    cdum = MAXLOC ( diff_coef )
    WRITE ( iw, '(A,T71,I10)' ) &
      ' COEFFICENT WITH MAX ABSOLUTE ERROR IS ', cdum
    WRITE ( iw, '(A,T21,3G20.14)' ) ' F ANAL ON COEF =', &
      dbg % f_ovlp_coef ( cdum )
    WRITE ( iw,'(A,T21,3G20.14)' ) ' F NUMR ON COEF =', &
      f_coef ( cdum )
    IF ( .NOT. ANY ( dbg % f_ovlp_coef == 0._dbl ) ) THEN
      cdum = MAXLOC ( rel_coef )
      WRITE ( iw,'(A,T71,I10)' ) &
        ' COEFFICENT WITH MAX RELATIVE ERROR IS ', cdum   
      WRITE ( iw,'(A,T21,3G20.14)' ) ' F ANAL ON COEF =', &
        dbg % f_ovlp_coef ( cdum )
      WRITE ( iw,'(A,T21,3G20.14)' ) ' F NUMR ON COEF =', &
        f_coef ( cdum )
    ENDIF
  END IF
  WRITE ( iw, * ) ' '
  WRITE ( iw, '(16x,A)' ) &
    '***************************************************'
  WRITE ( iw, '(16x,A)' ) '                    END POL DEBUG'
  WRITE ( iw, '(16x,A)' ) &
    '***************************************************'
END SUBROUTINE debug_control
!-----------------------------------------------------------------------------!
END MODULE pol_debug
!-----------------------------------------------------------------------------!

