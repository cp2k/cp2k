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

  USE ao_types,                        ONLY: ao_type
  USE atomic_kinds,                    ONLY: kind_info_type
  USE basis_set_types,                 ONLY: gto_basis_set_type
  USE coefficient_lists,               ONLY: coef_list_grid
  USE dg_rho0s,                        ONLY: dg_rho0_setup
  USE dg_types,                        ONLY: dg_type
  USE dgs,                             ONLY: dg_find_cutoff,&
                                             dg_pme_grid_setup,&
                                             dg_pol_grid_setup
  USE ewald_parameters_types,          ONLY: ewald_parameters_type
  USE ewalds,                          ONLY: ewald_self
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE mathconstants,                   ONLY: pi
  USE md,                              ONLY: thermodynamic_type
  USE molecule_types,                  ONLY: particle_node_type
  USE particle_lists,                  ONLY: particle_list_grid,&
                                             particle_list_select
  USE particle_types,                  ONLY: particle_type
  USE pw_grid_types,                   ONLY: HALFSPACE,&
                                             pw_grid_type
  USE pw_grids,                        ONLY: pw_find_cutoff,&
                                             pw_grid_construct,&
                                             pw_grid_setup
  USE qs_interactions,                 ONLY: exp_radius
  USE realspace_grid_types,            ONLY: realspace_grid_type,&
                                             rs_grid_setup
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

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
  

    TYPE(dg_type), DIMENSION(:), POINTER     :: dg_coef, dg_atom
    TYPE(pw_grid_type), DIMENSION(:), &
      POINTER                                :: small_grid, big_grid
    INTEGER, intent(IN)                      :: ntc, nta, ng

    INTEGER                                  :: ios

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
  
  
    TYPE(kind_info_type), DIMENSION(:), &
      INTENT(IN)                             :: drho_basis_info, &
                                                rho0_basis_info
    INTEGER, INTENT(OUT)                     :: ngc, nga, ntc, nta

    INTEGER                                  :: i, iset, ishell, nkinda, &
                                                nkindc, nshell_set

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

  
    TYPE(dg_type), DIMENSION(:), POINTER     :: dg
    TYPE(kind_info_type), DIMENSION(:), &
      POINTER                                :: ki
    INTEGER, INTENT(INOUT)                   :: igrid

    INTEGER                                  :: icc, igt, ikind, ios, ip, &
                                                iset, ishell, ncc, nshell
    TYPE(gto_basis_set_type), POINTER        :: gto_info

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

              ncc = SIZE ( gto_info % gcc ,1 )

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

                 ncc = SIZE ( gto_info % gcc ,1 )

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
  

    TYPE(dg_type), DIMENSION(:), INTENT(OUT) :: dg_coef, dg_atom
    TYPE(kind_info_type), DIMENSION(:), &
      POINTER                                :: drho_basis_info, &
                                                rho0_basis_info
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(pw_grid_type), DIMENSION(:), &
      INTENT(INOUT)                          :: small_grid, big_grid
    INTEGER, INTENT(OUT)                     :: iref, npts_max( : )
    REAL(dbl), INTENT(OUT)                   :: cutoff_radius_max, max_cut

    INTEGER                                  :: igrid, igrid_old, igt, ikind, &
                                                npts_s( 3 )
    REAL(dbl)                                :: cutoff, cutoff_radius, &
                                                max_cut0
    TYPE(kind_info_type), POINTER            :: ki

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
SUBROUTINE pol_grids_initialize ( ewald_param, dg_coef, dg_atom, part, &
     pnode, thermo, ao, drho_basis_info, rho0_basis_info, pnode_grp, box,  &
     small_grid, big_grid, iref, npts_max, cutoff_radius_max,  &
     pme_small_grid, pme_big_grid, ewald_grid, dg_ewald, ref_grid )
  

    TYPE(ewald_parameters_type), &
      INTENT(INOUT)                          :: ewald_param
    TYPE(dg_type), DIMENSION(:), INTENT(OUT) :: dg_coef, dg_atom
    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT)                          :: part
    TYPE(particle_node_type), DIMENSION(:), &
      INTENT(IN)                             :: pnode
    TYPE(thermodynamic_type), INTENT(INOUT)  :: thermo
    TYPE(ao_type), INTENT(inout)             :: ao
    TYPE(kind_info_type), DIMENSION(:), &
      POINTER                                :: drho_basis_info, &
                                                rho0_basis_info
    TYPE(global_environment_type), &
      INTENT(IN)                             :: pnode_grp
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(pw_grid_type), DIMENSION(:), &
      INTENT(INOUT)                          :: small_grid, big_grid
    INTEGER, INTENT(IN)                      :: iref, npts_max( 3 )
    REAL(dbl), INTENT(IN)                    :: cutoff_radius_max
    TYPE(pw_grid_type), OPTIONAL, POINTER    :: pme_small_grid, pme_big_grid, &
                                                ewald_grid
    TYPE(dg_type), OPTIONAL, POINTER         :: dg_ewald
    TYPE(pw_grid_type), INTENT(IN), OPTIONAL :: ref_grid

    INTEGER                                  :: gmax, igrid, igrid_old, igt, &
                                                ikind, isos, iw, natoms, &
                                                npts_s( 3 )
    REAL(dbl)                                :: alphasq, cutoff, &
                                                cutoff_radius, ecut, norm
    TYPE(kind_info_type), POINTER            :: ki
    TYPE(realspace_grid_type)                :: rs

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
! Fourth: setup the Ewald (if necessary) with respect to referece grid.

! parallelisation is over atoms (pnodes), so the group of processors
! has to be the same as the group for the pnodes
! This holds only for normal Ewald sums!
! For the grid based methods, G-Space/Real Space are distributed differently
  ewald_param % group = pnode_grp % group

  natoms = SIZE ( part )

  IF ( pnode_grp % ionode ) THEN
    iw = pnode_grp % scr
    IF ( ewald_param % ewald_type /= 'NONE' ) THEN

       WRITE ( iw, '( A,T67,A14 )' ) ' Ewald| Summation is done by:', &
            ADJUSTR(ewald_param % ewald_type)
       WRITE ( iw, '( A,T71,F10.4 )' ) ' Ewald| Alpha parameter [A]', &
            ewald_param % alpha

       SELECT CASE ( ewald_param % ewald_type )
       CASE DEFAULT
          WRITE ( iw, '( A,T71,I10 )' ) &
               ' Ewald| G-space max. Miller index', ewald_param % gmax
       CASE ( 'PME')
          WRITE ( iw, '( A,T71,I10 )' ) &
               ' EWALD| Max small-grid points (input) ', ewald_param % ns_max
          WRITE ( iw, '( A,T71,E10.4 )' ) &
               ' EWALD| Gaussian tolerance (input) ', ewald_param % epsilon
       CASE ( 'SPME' )
          WRITE ( iw, '( A,T71,I10 )' ) &
               ' Ewald| G-space max. Miller index', ewald_param % gmax
          WRITE ( iw, '( A,T71,I10 )' ) &
               ' Ewald| Spline interpolation order ', ewald_param % o_spline
       END SELECT

    ELSE

       WRITE ( iw, '( A, T73, A )' ) ' Ewald| ','not used'

    END IF
  ELSE
    iw = -1
  END IF

! fire up the reciprocal space and compute self interaction and
! term from the neutralizing background.

  IF ( ewald_param % ewald_type /= 'NONE' ) THEN

     CALL ewald_self ( ewald_param, thermo, pnode )
     
     IF ( PRESENT ( ewald_grid ) ) THEN
        ! set up EWALD "Classic"
        NULLIFY ( ewald_grid )
        ALLOCATE ( ewald_grid, STAT = isos )
        IF ( isos /= 0 ) CALL stop_memory ( 'pol_grids', 'ewald_grid' )
        CALL pw_grid_construct ( ewald_grid )
        gmax = ewald_param % gmax
        IF ( gmax == 2 * ( gmax / 2 ) ) THEN
           IF ( pnode_grp % ionode ) &
             CALL stop_program ( "initialize_ewalds", "gmax has to be odd" )
        END IF
        ewald_grid % bounds ( 1, : ) = -gmax/2
        ewald_grid % bounds ( 2, : ) = +gmax/2

        npts_s = (/ gmax, gmax, gmax /)

        ewald_grid % grid_span = HALFSPACE

        CALL pw_find_cutoff ( npts_s, box, cutoff )

        ecut = 0.5_dbl * cutoff * cutoff
        CALL pw_grid_setup( box, ewald_grid, cutoff = ecut, &
                            info = iw, fft_usage = .FALSE. )
     
     ELSEIF ( PRESENT ( pme_small_grid ) .AND. PRESENT ( pme_big_grid ) ) THEN
! set up PME, the double grid method
        NULLIFY ( dg_ewald, pme_small_grid, pme_big_grid )
        ALLOCATE ( dg_ewald, STAT = isos )
        IF ( isos /= 0 ) CALL stop_memory ( 'pol_grids', 'dg_ewald' )
        ALLOCATE ( pme_small_grid, STAT = isos )
        IF ( isos /= 0 ) CALL stop_memory ( 'pol_grids', 'pme_small_grid' )
        ALLOCATE ( pme_big_grid, STAT = isos )
        IF ( isos /= 0 ) CALL stop_memory ( 'pol_grids', 'pme_big_grid' )
        CALL pw_grid_construct ( pme_small_grid )
        CALL pw_grid_construct ( pme_big_grid )

        npts_s ( : ) = ewald_param % ns_max

! compute cut-off radius
        alphasq = ( ewald_param % alpha ) ** 2
        norm = ( 2.0_dbl * alphasq / pi ) ** ( 1.5_dbl )
        cutoff_radius = exp_radius ( 0, 2._dbl * alphasq, ewald_param % epsilon , norm )

        IF ( PRESENT ( ref_grid ) ) THEN
          CALL dg_pme_grid_setup ( box, npts_s, cutoff_radius, &
               pme_small_grid, pme_big_grid, pnode_grp, ref_grid )
        ELSE
          CALL dg_pme_grid_setup ( box, npts_s, cutoff_radius, &
               pme_small_grid, pme_big_grid, pnode_grp, big_grid ( iref ) )
        ENDIF

        CALL rs_grid_setup ( rs, pme_big_grid, -1, pnode_grp % scr )

        ! set up a list with all particles having a charge
        CALL particle_list_select ( pnode, dg_ewald % plist, dg_ewald % nparts )

        IF ( pnode_grp % ionode ) THEN
          WRITE ( iw, '( A,T71,E10.4 )' ) &
             ' EWALD| Gaussian tolerance (effective) ', ewald_param % epsilon
          WRITE ( iw, '( A,T63,3I6 )' ) &
             ' EWALD| Small box grid ', pme_small_grid % npts
          WRITE ( iw, '( A,T63,3I6 )' ) &
             ' EWALD| Full box grid ', pme_big_grid % npts
        END IF
     END IF

  END IF
END SUBROUTINE pol_grids_initialize

!******************************************************************************
 SUBROUTINE assign_coefs_to_grids( ao, part, ki, ngt )

  
    TYPE(ao_type), INTENT(INOUT)             :: ao
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    TYPE(kind_info_type), DIMENSION(:), &
      INTENT(IN)                             :: ki
    INTEGER, INTENT(IN)                      :: ngt

    INTEGER :: first_cgf, grid_identifier, grid_type_identifier, icgf, icoef, &
      ii, ikind, ipart, iset, ishell, last_cgf, ncgf, ngrid, ngrid_kind, &
      ngrid_types, ngrid_types_kind, nshell

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
