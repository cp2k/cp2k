!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE initialize_ewalds 
  
  USE dg, ONLY : dg_grid_setup
  USE ewalds, ONLY : ewald_self
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE io_parameters, ONLY : scr
  USE kinds, ONLY : dbl
  USE md, ONLY : thermodynamic_type
  USE molecule_types, ONLY : particle_node_type
  USE particle_types, ONLY : particle_type
  USE pw_grid_types, ONLY : pw_grid_type
  USE pw_grids, ONLY : pw_grid_setup
  USE simulation_cell, ONLY : cell_type
  
  PUBLIC :: initialize_ewald
  PRIVATE
  
CONTAINS

!******************************************************************************

SUBROUTINE initialize_ewald ( part, pnode, ewald_param, box, thermo, &
     ewald_grid, pme_small_grid, pme_big_grid )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( ewald_parameters_type ), INTENT ( OUT ) :: ewald_param
  TYPE ( particle_type ), INTENT ( IN ) :: part ( : )
  TYPE ( particle_node_type ), INTENT ( IN ) :: pnode ( : )
  TYPE ( cell_type ), INTENT ( IN ) :: box
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE ( pw_grid_type ), INTENT ( OUT ), OPTIONAL :: ewald_grid
  TYPE ( pw_grid_type ), INTENT ( OUT ), OPTIONAL :: pme_small_grid
  TYPE ( pw_grid_type ), INTENT ( OUT ), OPTIONAL :: pme_big_grid
  
! Locals
  INTEGER :: natoms, iw, iat, jat, gmax, npts ( 3 )
  REAL ( dbl ) :: qi, qj, i, cutoff
  
!------------------------------------------------------------------------------
  
  iw = scr
  natoms = SIZE ( part )
  
  IF ( ewald_param % ewald_type /= 'NONE' ) THEN
     
     WRITE ( iw, '( A,T71,A )' ) '  Ewald summation is done by:', &
          ewald_param % ewald_type
     WRITE ( iw, '( A,T71,F10.4 )' ) '  Ewald alpha parameter [A]', &
          ewald_param % alpha
     
     SELECT CASE(ewald_param % ewald_TYPE ( 1:3))
     CASE DEFAULT
        WRITE ( iw, '( A,T71,I10 )' ) &
             '  Ewald G-space max. Miller index', ewald_param % gmax
        
     CASE ( 'PME')
        WRITE ( iw, '( A,T71,I10 )' ) &
             ' PME max small-grid points ', ewald_param % ns_max
        WRITE ( iw, '( A,T71,F10.4 )' ) &
             ' PME gaussian tolerance ', ewald_param % epsilon
     END SELECT
     
  ELSE
     
     WRITE( iw, '( A )' ) ' No Ewald summation is performed'
     
  END IF
  
! fire up the reciprocal space and compute self interaction and
! term from the neutralizing background.
  
  IF ( ewald_param % ewald_type /= 'NONE' ) THEN
     
     CALL ewald_self ( ewald_param, thermo, pnode )
     
! set up EWALD_GAUSS
     IF ( PRESENT ( ewald_grid)) THEN
        gmax = ewald_param % gmax
        IF ( gmax == 2 * ( gmax / 2 ) ) THEN
           CALL stop_prg ( "initialize_ewalds", "gmax has to be odd" )
        END IF
        ewald_grid % bounds ( 1, : ) = -gmax / 2
        ewald_grid % bounds ( 2, : ) = +gmax / 2
        npts=( / gmax, gmax, gmax / )
        
        CALL dg_find_cutoff ( npts, box, cutoff )
        
        CALL pw_grid_setup( box, ewald_grid, cutoff)
     END IF
! set up PME_GAUSS
     
!dbg
!      IF ( PRESENT ( pme_small_grid).AND. PRESENT ( pme_big_grid)) 
!        CALL CALL dg_setup(box, npts_s, epsilon, alpha, pme_small_grid, &
!                     pme_big_grid, ewald_param % ewald_type )
!dbg
     
  END IF
  
END SUBROUTINE initialize_ewald

!******************************************************************************

END MODULE initialize_ewalds
