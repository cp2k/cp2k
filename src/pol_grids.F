!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/pol_grids [1.0] *
!!
!!   NAME
!!     pol_grids
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE pol_grids

  USE atomic_kinds, ONLY :  kind_info_type
  USE ao_types, ONLY : ao_type
  USE basis_set_types, ONLY :  gto_basis_set_type
  USE global_types, ONLY :  global_environment_type
  USE dg_types, ONLY : dg_type 
  USE dgs, ONLY : dg_pol_grid_setup, dg_find_cutoff 
  USE dg_rho0s, ONLY :  dg_rho0_type, dg_rho0_setup  
  USE kinds, ONLY : dbl
  USE mathconstants, ONLY : pi, zero
  USE molecule_types, ONLY : particle_node_type
  USE message_passing, ONLY : mp_sum
  USE particle_types, ONLY : particle_type
  USE particle_lists, ONLY : particle_list_type, particle_list_grid
  USE coefficient_lists, ONLY : coef_list_type, coef_list_grid
  USE pw_grid_types, ONLY : pw_grid_type, HALFSPACE
  USE simulation_cell, ONLY : cell_type
  USE structure_factor_types, ONLY : structure_factor_type
  USE structure_factors, ONLY : structure_factor_evaluate, &
                                structure_factor_allocate, &
                                structure_factor_deallocate
  USE termination, ONLY : stop_memory, stop_program
  USE timings, ONLY : timeset, timestop
  USE util, ONLY : matvec_3x3

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: pol_grids_allocate, pol_grids_count, pol_grids_initialize
  PUBLIC :: pol_grids_map, assign_coefs_to_grids, pol_grids_ref_grid                       
  
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** pol_grids/pol_grids_allocate [1.0] *
!!
!!   NAME
!!     pol_grids_allocate
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     gloria
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE pol_grids_allocate ( dg_coef, dg_atom, small_grid, big_grid, & 
                   ntc, nta, ng )
  
  IMPLICIT NONE

! Arguments
  TYPE ( dg_type ), POINTER, DIMENSION (:) :: dg_coef
  TYPE ( dg_type ), POINTER, DIMENSION (:) :: dg_atom
  TYPE ( pw_grid_type ), POINTER, DIMENSION (:) :: small_grid
  TYPE ( pw_grid_type ), POINTER, DIMENSION (:) :: big_grid
  INTEGER, intent(IN) :: ntc
  INTEGER, intent(IN) :: nta
  INTEGER, intent(IN) :: ng
  
! Locals
  INTEGER ::  ios                                         

!---------------------------------------------

  NULLIFY ( dg_coef, dg_atom )
  NULLIFY ( small_grid, big_grid )

  ALLOCATE ( dg_coef ( ntc ), STAT = ios )
  IF ( ios /= 0 ) THEN
    CALL stop_memory ( 'allocate_pol_grids', 'dg_coef', ntc )
  END IF

  ALLOCATE ( dg_atom ( nta ), STAT = ios )
  IF ( ios /= 0 ) THEN
    CALL stop_memory ( 'allocate_pol_grids', 'dg_atom', nta )
  END IF

  ALLOCATE ( small_grid ( ng ), STAT = ios )
  IF ( ios /= 0 ) THEN
    CALL stop_memory ( 'allocate_pol_grids', 'small_grid', ng )
  END IF

  ALLOCATE ( big_grid ( ng ), STAT = ios )
  IF ( ios /= 0 ) THEN
    CALL stop_memory ( 'allocate_pol_grids', 'big_grid', ng )
  END IF

 END SUBROUTINE pol_grids_allocate

!******************************************************************************
!!****** pol_grids/pol_grids_count [1.0] *
!!
!!   NAME
!!     pol_grids_count
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     gloria
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pol_grids_count ( drho_basis_info, rho0_basis_info, ngc, nga, ntc, nta )
  
  IMPLICIT NONE
  
  TYPE ( kind_info_type ), INTENT(IN), DIMENSION ( : ) :: drho_basis_info
  TYPE ( kind_info_type ), INTENT(IN), DIMENSION ( : ) :: rho0_basis_info
  INTEGER, INTENT ( OUT ) :: ngc
  INTEGER, INTENT ( OUT ) :: nga
  INTEGER, INTENT ( OUT ) :: ntc
  INTEGER, INTENT ( OUT ) :: nta

! locals

  INTEGER :: i, iset,  nshell_set, ishell, nkinda, nkindc

  nkindc = size(drho_basis_info)
  nkinda = size(rho0_basis_info)
  nshell_set = 0
  nga = 0
  ngc = 0
  nta = 0
  ntc = 0

! coeffs loop
    DO i=1, nkindc                   
      DO iset = 1, drho_basis_info(i) % orb_basis_set% nset
        nshell_set =  drho_basis_info(i)%orb_basis_set % nshell(iset)
        ngc = ngc + nshell_set
        DO ishell = 1, nshell_set
          ntc = ntc + ( 2 * drho_basis_info ( i ) % orb_basis_set % l ( ishell, iset ) + 1 )
        ENDDO
      ENDDO
    ENDDO

   nshell_set = 0

! atom  loop
    DO i=1, nkinda                   
      DO iset = 1, rho0_basis_info(i) % orb_basis_set%nset
        nshell_set = rho0_basis_info(i)%orb_basis_set % nshell(iset)
        nga = nga + nshell_set
        DO ishell = 1, nshell_set
          nta = nta +(2*rho0_basis_info(i)% orb_basis_set % l(ishell,iset)+1)
        ENDDO
      ENDDO
    ENDDO
    
  END SUBROUTINE pol_grids_count                                       


!******************************************************************************
!!****** pol_grids/pol_grids_map [1.0] *
!!
!!   NAME
!!     pol_grids_map
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     gloria
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pol_grids_map ( dg, ki, igrid )

  IMPLICIT NONE
  
  TYPE ( dg_type ), POINTER, DIMENSION ( : ) :: dg
  TYPE ( kind_info_type ), POINTER, DIMENSION ( : )  :: ki
  INTEGER, INTENT ( INOUT ) :: igrid

! locals
  INTEGER :: igt, ikind, iset, ishell, nshell, ncc, icc, ios, ip
  TYPE ( gto_basis_set_type ), POINTER :: gto_info

  igt = 0
  
  DO ikind=1,SIZE ( ki ) 
  
     gto_info => ki ( ikind ) % orb_basis_set
     
     DO iset = 1, gto_info%nset
     
        nshell = gto_info%nshell ( iset )

        DO ishell = 1, nshell

           igrid = igrid + 1
   
           SELECT CASE ( gto_info % l ( ishell, iset ) )
     
           CASE ( 0 )
   
              igt = igt + 1 

              dg ( igt ) % grid_index = igrid
              dg ( igt ) % dg_rho0 % type  = "S"
              dg ( igt ) % dg_rho0 % grid = igrid 
              dg ( igt ) % dg_rho0 % kind = ikind 
              dg ( igt ) % dg_rho0 % cutoff_radius = gto_info % set_radius ( iset )
              dg ( igt ) % nparts = SIZE ( ki ( ikind ) % atom_list )

              ncc = SIZE ( gto_info % gcc ( :, ishell, iset ) )

              ALLOCATE ( dg ( igt ) % dg_rho0 % gcc ( ncc ), STAT = ios )
              IF (ios /= 0 ) THEN
                 CALL stop_memory ('pol_grids_map', 'gcc', ncc)
              END IF

              ALLOCATE ( dg ( igt ) % dg_rho0 % zet ( ncc ), STAT = ios )
              IF ( ios /= 0 ) THEN
                 CALL stop_memory ( 'pol_grids_map', 'zet', ncc )
              END IF 

              DO icc=1, ncc
                dg ( igt ) % dg_rho0 % zet ( icc ) = gto_info % zet ( icc, iset )
                dg ( igt ) % dg_rho0 % gcc ( icc ) = gto_info % gcc ( icc, ishell, iset )
              ENDDO                                             

           CASE ( 1 )

               DO ip = 1, 3

                 igt = igt + 1
 
                 SELECT CASE ( ip )
                 CASE ( 1 )
                  dg ( igt ) % dg_rho0 % type = "PX"
                 CASE ( 2 )
                  dg ( igt ) % dg_rho0 % type = "PY"
                 CASE ( 3 )
                  dg ( igt ) % dg_rho0 % type = "PZ"
                 END SELECT                       

                 dg ( igt ) % grid_index = igrid
                 dg ( igt ) % dg_rho0 % grid = igrid 
                 dg ( igt ) % dg_rho0 % kind = ikind 
                 dg ( igt ) % dg_rho0 % cutoff_radius = gto_info % set_radius ( iset )
                 dg ( igt ) % nparts = SIZE ( ki ( ikind ) % atom_list )

                 ncc = SIZE ( gto_info % gcc ( :, ishell, iset ) )

                 ALLOCATE ( dg ( igt ) % dg_rho0 % gcc ( ncc ), STAT = ios )
                 IF ( ios /= 0 ) THEN
                    CALL stop_memory ( 'pol_grids_map', 'gcc', ncc )
                 END IF

                 ALLOCATE ( dg ( igt ) % dg_rho0 % zet ( ncc ), STAT = ios )
                 IF ( ios /= 0 ) THEN
                    CALL stop_memory ( 'pol_grids_map', 'zet', ncc )
                 END IF
                  
                 DO icc = 1, ncc
                   dg ( igt ) % dg_rho0 % zet ( icc ) = gto_info % zet ( icc, iset )
                   dg ( igt ) % dg_rho0 % gcc ( icc ) = gto_info % gcc ( icc, ishell, iset )
                 ENDDO                                             
 
                 END DO

           END SELECT
   
        ENDDO

     ENDDO

  ENDDO                   

 END SUBROUTINE pol_grids_map 
!******************************************************************************
!*****
!!****** pol_grids/pol_grids_ref_grid [1.0] *
!!
!!   NAME
!!     pol_grids_ref_grid
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     gloria
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pol_grids_ref_grid ( dg_coef, dg_atom, drho_basis_info, rho0_basis_info, box, &
                             small_grid, big_grid, iref, npts_max, cutoff_radius_max, max_cut )
  
  IMPLICIT NONE

! Arguments
  TYPE ( dg_type ), INTENT ( OUT ), DIMENSION (:) :: dg_coef
  TYPE ( dg_type ), INTENT ( OUT ), DIMENSION (:) :: dg_atom
  TYPE ( kind_info_type ), POINTER, DIMENSION ( : ) :: drho_basis_info
  TYPE ( kind_info_type ), POINTER, DIMENSION ( : ) :: rho0_basis_info
  TYPE ( cell_type ), INTENT ( IN ) :: box
  TYPE ( pw_grid_type ), INTENT ( INOUT ), DIMENSION (:) :: small_grid
  TYPE ( pw_grid_type ), INTENT ( INOUT ), DIMENSION (:) :: big_grid
  INTEGER, INTENT (OUT) :: iref
  INTEGER, INTENT (OUT) :: npts_max ( : )
  REAL ( dbl ), INTENT ( OUT ) :: cutoff_radius_max, max_cut
  
! Locals
  INTEGER ::  npts_s ( 3 )
  INTEGER :: igt, ikind, igrid, ig, igrid_old
  REAL ( dbl ) :: i, cutoff, cutoff_radius, max_cut0
  TYPE (kind_info_type), POINTER :: ki

!------------------------------------------------------------------------------
! initialize igrid_old to negative number as to not
! be activated by accident
  igrid_old = -1

! set up grids

! First: get the grid with max cutoff

  iref = 0
  max_cut0 = 1E-30_dbl


  DO igt = 1, SIZE ( dg_coef )
       
     ikind = dg_coef ( igt ) % dg_rho0 % kind
     ki => drho_basis_info( ikind ) 
     npts_s( : ) = ki % number_of_grid_points
     igrid = dg_coef ( igt ) % dg_rho0 % grid
     cutoff_radius = dg_coef ( igt ) % dg_rho0 % cutoff_radius

     IF ( igrid == igrid_old ) CYCLE

     CALL dg_find_cutoff ( box, npts_s, cutoff_radius, small_grid ( igrid ),  &
                          big_grid ( igrid ), cutoff ) 
     max_cut = MAX ( max_cut0, cutoff )
    
     IF (max_cut /= max_cut0) THEN
       iref = igrid
       max_cut0 = max_cut
       npts_max = npts_s
       cutoff_radius_max = cutoff_radius
     END IF   

     igrid_old = igrid

  ENDDO

  DO igt = 1, SIZE ( dg_atom )
       
     ikind = dg_atom ( igt ) % dg_rho0 % kind
     ki => rho0_basis_info( ikind ) 
     npts_s( : ) = ki % number_of_grid_points
     igrid = dg_atom ( igt ) % dg_rho0 % grid
     cutoff_radius = dg_atom ( igt ) % dg_rho0 % cutoff_radius

     IF ( igrid == igrid_old ) CYCLE

     CALL dg_find_cutoff ( box, npts_s, cutoff_radius, small_grid ( igrid ),  &
                          big_grid ( igrid ), cutoff ) 

     max_cut = MAX ( max_cut0, cutoff )
    
     IF (max_cut /= max_cut0) THEN
       iref = igrid
       max_cut0 = max_cut
       npts_max = npts_s
       cutoff_radius_max = cutoff_radius
     END IF   

     igrid_old = igrid

  ENDDO

END SUBROUTINE pol_grids_ref_grid
!******************************************************************************
!*****
!!****** pol_grids/pol_grids_initialize [1.0] *
!!
!!   NAME
!!     pol_grids_initialize
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     gloria
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE pol_grids_initialize ( dg_coef, dg_atom,  part, pnode, ao, &
     drho_basis_info, rho0_basis_info, pnode_grp, box, small_grid, &
     big_grid, iref, npts_max, cutoff_radius_max, ref_grid )
  
  IMPLICIT NONE

! Arguments
  TYPE ( dg_type ), INTENT ( OUT ), DIMENSION (:) :: dg_coef
  TYPE ( dg_type ), INTENT ( OUT ), DIMENSION (:) :: dg_atom
  TYPE ( particle_type ), DIMENSION ( : ), INTENT ( INOUT ) :: part
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( IN ) :: pnode
  TYPE ( ao_type ), INTENT (inout) :: ao
  TYPE ( kind_info_type ), POINTER, DIMENSION ( : ) :: drho_basis_info
  TYPE ( kind_info_type ), POINTER, DIMENSION ( : ) :: rho0_basis_info
  TYPE ( global_environment_type ), INTENT ( IN ) :: pnode_grp
  TYPE ( cell_type ), INTENT ( IN ) :: box
  TYPE ( pw_grid_type ), INTENT ( INOUT ), DIMENSION (:) :: small_grid
  TYPE ( pw_grid_type ), INTENT ( INOUT ), DIMENSION (:) :: big_grid
  INTEGER, INTENT (IN) :: iref, npts_max ( 3 )
  REAL ( dbl ), INTENT ( IN ) :: cutoff_radius_max
  TYPE ( pw_grid_type ), INTENT ( IN ), OPTIONAL :: ref_grid
  
! Locals
  INTEGER ::  natoms, npts_s ( 3 ), iw
  INTEGER :: igt, ikind, igrid, ig, igrid_old
  REAL ( dbl ) :: i, cutoff, cutoff_radius
  TYPE (kind_info_type), POINTER :: ki

  iw = pnode_grp % scr 
! first: setup the reference grid

  cutoff_radius = cutoff_radius_max
  IF ( PRESENT ( ref_grid ) ) THEN
    CALL dg_pol_grid_setup ( box, npts_max, cutoff_radius, &
         small_grid ( iref ), big_grid ( iref ), pnode_grp, ref_grid )
  ELSE
    CALL dg_pol_grid_setup ( box, npts_max, cutoff_radius, &
         small_grid ( iref ), big_grid ( iref ), pnode_grp )
  ENDIF

! Third: setup rest of the grids w.r.t. the reference grid

! coef loop

  DO igt = 1, SIZE ( dg_coef )
       
     ikind = dg_coef ( igt ) % dg_rho0 % kind
     ki => drho_basis_info( ikind ) 
     npts_s( : ) = ki % number_of_grid_points
     igrid = dg_coef ( igt ) % dg_rho0 % grid
     cutoff_radius = dg_coef ( igt ) % dg_rho0 % cutoff_radius

     CALL coef_list_grid ( ao, dg_coef ( igt ) % clist, part, &
                         dg_coef ( igt ) % ncoefs, igt, igrid )


     IF ( igrid /= igrid_old ) THEN

       IF ( igrid /= iref ) THEN
         IF ( PRESENT ( ref_grid ) ) THEN
           CALL dg_pol_grid_setup ( box, npts_s, cutoff_radius, &
             small_grid ( igrid ), big_grid ( igrid ), &
             pnode_grp, ref_grid  )
         ELSE
           CALL dg_pol_grid_setup ( box, npts_s, cutoff_radius , &
             small_grid ( igrid ), big_grid ( igrid ), &
             pnode_grp, big_grid ( iref )  )
         ENDIF
       ENDIF

     ENDIF

     CALL dg_rho0_setup ( dg_coef ( igt ) % dg_rho0, small_grid ( igrid ) )

     IF ( igrid /= igrid_old ) THEN
       IF ( igrid /= iref ) THEN
         IF ( pnode_grp % ionode ) THEN
           WRITE ( iw, '( A,T63,3I6 )' ) &
            ' POL_GRIDS_INIT| grid number   ', igrid                  
           WRITE ( iw, '( A,T63,3I6 )' ) &
            ' POL_GRIDS_INIT| Small box grid ', small_grid ( igrid ) % npts
           WRITE ( iw, '( A,T63,3I6 )' ) &
            ' POL_GRIDS_INIT| Full box grid ', big_grid ( igrid ) % npts
         END IF
       ELSEIF ( igrid == iref .AND. .NOT.PRESENT ( ref_grid ) ) THEN
         IF ( pnode_grp % ionode ) THEN
           WRITE ( iw, '( A,T63,3I6 )' ) &
            ' POL_GRIDS_INIT| Ref. grid number   ', iref
           WRITE ( iw, '( A,T63,3I6 )' ) &
            ' POL_GRIDS_INIT| Ref. small box grid ', small_grid ( iref ) % npts
           WRITE ( iw, '( A,T63,3I6 )' ) &
            ' POL_GRIDS_INIT| Ref. full box grid ', big_grid ( iref ) % npts
         END IF
       ENDIF
     ENDIF

     igrid_old = igrid

  END DO

! atom Loop

  DO igt = 1, SIZE ( dg_atom )

     ikind = dg_atom ( igt ) % dg_rho0 % kind
     ki => rho0_basis_info ( ikind ) 
     npts_s ( : ) = ki % number_of_grid_points
     igrid = dg_atom ( igt ) % dg_rho0 % grid
     cutoff_radius = dg_atom ( igt ) % dg_rho0 % cutoff_radius

     CALL particle_list_grid ( part, dg_atom ( igt ) % plist, ki % atom_list, &
                 dg_atom ( igt ) % nparts, igrid )
       

     IF ( igrid /= igrid_old ) THEN

       IF ( igrid /= iref ) THEN
         IF ( PRESENT ( ref_grid ) ) THEN
            CALL dg_pol_grid_setup ( box, npts_s, cutoff_radius , &
                small_grid ( igrid ), big_grid ( igrid ), &
                pnode_grp, ref_grid )
         ELSE
            CALL dg_pol_grid_setup ( box, npts_s, cutoff_radius , &
                small_grid ( igrid ), big_grid ( igrid ), &
                pnode_grp, big_grid ( iref ) )
         ENDIF
       ENDIF
     ENDIF

     CALL dg_rho0_setup ( dg_atom ( igt ) % dg_rho0, small_grid ( igrid ) )

     IF ( igrid /= igrid_old ) THEN
       IF ( igrid /= iref ) THEN
         IF ( pnode_grp % ionode ) THEN
           WRITE ( iw, '( A,T63,3I6 )' ) &
            ' POL_GRIDS_INIT| grid number   ', igrid                  
           WRITE ( iw, '( A,T63,3I6 )' ) &
            ' POL_GRIDS_INIT| Small box grid ', small_grid ( igrid ) % npts
           WRITE ( iw, '( A,T63,3I6 )' ) &
            ' POL_GRIDS_INIT| Full box grid ', big_grid ( igrid ) % npts
         END IF
       ELSEIF ( igrid == iref .AND. .NOT.PRESENT ( ref_grid ) ) THEN
         IF ( pnode_grp % ionode ) THEN
           WRITE ( iw, '( A,T63,3I6 )' ) &
            ' POL_GRIDS_INIT| Ref. grid number   ', iref
           WRITE ( iw, '( A,T63,3I6 )' ) &
            ' POL_GRIDS_INIT| Ref. small box grid ', small_grid ( iref ) % npts
           WRITE ( iw, '( A,T63,3I6 )' ) &
            ' POL_GRIDS_INIT| Ref. full box grid ', big_grid ( iref ) % npts
         END IF
       ENDIF
     ENDIF
     igrid_old = igrid

  END DO

END SUBROUTINE pol_grids_initialize

!******************************************************************************
 SUBROUTINE assign_coefs_to_grids( ao, part, ki, ngt )

  IMPLICIT NONE
  
  TYPE (ao_type), INTENT (INOUT) :: ao
  TYPE (particle_type), INTENT(IN), DIMENSION (:) :: part
  TYPE (kind_info_type), DIMENSION (:), INTENT(IN) :: ki
  INTEGER, INTENT(IN) :: ngt

! locals 
  INTEGER :: ngrid_types, ngrid, ngrid_types_kind, ngrid_kind, first_cgf, last_cgf 
  INTEGER :: ii, icoef, icgf, ikind, nshell, ncgf, ipart, iset, ishell
  INTEGER :: grid_type_identifier, grid_identifier        

  ngrid_types = 0
  ngrid = 0
  ngrid_types_kind = 0
  ngrid_kind = 0
  grid_type_identifier = 0
  grid_identifier = 0
  
  DO ikind = 1 , size ( ki )
   
     ngrid_types = ngrid_types + ngrid_types_kind
     ngrid = ngrid + ngrid_kind 
 
     DO ii = 1, ki ( ikind ) % natom

! nshells, total number of shells, gives the number of grids of the kind 
! ncgf, total number of CGF, gives the number of grid types of the kind 
! they are set to 0 here because you are looping on ii, the total
! number of atoms corresponding to each basis kind.
        nshell = 0
        ncgf = 0
        ipart = ki ( ikind ) % atom_list ( ii )

        DO iset = 1, ki ( ikind ) % orb_basis_set% nset
           DO ishell = 1, ki ( ikind ) % orb_basis_set % nshell ( iset )
              nshell = nshell + 1
              grid_identifier = ngrid + nshell
              first_cgf = ki ( ikind ) % orb_basis_set % first_cgf ( ishell, iset )
              last_cgf = ki ( ikind ) % orb_basis_set % last_cgf ( ishell, iset )
              DO icgf = first_cgf, last_cgf
                 ncgf = ncgf + 1
! into this loop, ncgf and icgf are the same, you can use both of them,
! They represent the PARTIAL (i.e. atomic) coefficient index
! icoef is different from icgf, it is the GLOBAL (i.e. whole system) coefficient index
                 grid_type_identifier = ngrid_types + icgf
                 icoef  = part ( ipart ) % coef_list ( icgf )
!!!!!
! here the identifiers are assigned
! gt in pol_electrostatics corresponds to  ao % grid_id ( icoef )
! ig in pol_electrostatics corresponds to  ao % grid_type_id ( icoef )
!!!!!
                 ao % grid_id ( icoef ) = grid_identifier
                 ao % grid_type_id ( icoef ) = grid_type_identifier
              ENDDO
           ENDDO
        ENDDO
    END DO

!here the number of grids corresponding to each basis kind is updated
    ngrid_kind = nshell
!here the number of grid_types corresponding to each basis kind is updated
    ngrid_types_kind = ncgf
 
 END DO                                       

 IF ( grid_type_identifier /= ngt .AND. ngt /= 0 ) THEN
    CALL stop_program ( "assign_coef_to_grid", "inconsistent grid number" )  
 END IF

 END SUBROUTINE assign_coefs_to_grids


!******************************************************************************

END MODULE pol_grids

!******************************************************************************
