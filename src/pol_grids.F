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
  USE md, ONLY : thermodynamic_type
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
  PUBLIC :: pol_grids_map, assign_coefs_to_grids                       
  
!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!*****
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

  NULLIFY (dg_coef, dg_atom)
  NULLIFY (small_grid, big_grid)

  ALLOCATE (dg_coef(ntc), STAT = ios)
  IF (ios /= 0 ) THEN
    CALL stop_memory ('allocate_pol_grids', 'dg_coef', ntc)
  END IF

  ALLOCATE (dg_atom(nta), STAT = ios)
  IF (ios /= 0 ) THEN
    CALL stop_memory ('allocate_pol_grids', 'dg_atom', nta)
  END IF

  ALLOCATE (small_grid(ng), STAT = ios)
  IF (ios /= 0 ) THEN
    CALL stop_memory ('allocate_pol_grids', 'small_grid', ng)
  END IF

  ALLOCATE (big_grid(ng), STAT = ios)
  IF (ios /= 0 ) THEN
    CALL stop_memory ('allocate_pol_grids', 'big_grid', ng)
  END IF

 END SUBROUTINE pol_grids_allocate

!******************************************************************************
!*****
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

SUBROUTINE pol_grids_count (drho_basis_info, rho0_basis_info, ngc, nga, ntc, nta)
  
  IMPLICIT NONE
  
  TYPE ( kind_info_type ), INTENT(IN), DIMENSION ( : ) :: drho_basis_info
  TYPE ( kind_info_type ), INTENT(IN), DIMENSION ( : ) :: rho0_basis_info
  INTEGER, INTENT(out) :: ngc
  INTEGER, INTENT(out) :: nga
  INTEGER, INTENT(out) :: ntc
  INTEGER, INTENT(out) :: nta

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
          ntc = ntc +(2*drho_basis_info(i)% orb_basis_set % l(ishell,iset)+1)
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
!*****
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

SUBROUTINE pol_grids_map ( dg, ki, igrid)

  IMPLICIT NONE
  
  TYPE ( dg_type ), POINTER, DIMENSION (:) :: dg
  TYPE ( kind_info_type ), POINTER, DIMENSION ( : )  :: ki
  INTEGER, INTENT(INOUT) :: igrid

! locals
  INTEGER :: igt, ikind, iset, ishell, nshell, ncc, icc, ios, ip
  TYPE ( gto_basis_set_type ), POINTER :: gto_info

  igt = 0
  
  DO ikind=1,SIZE(ki) 
  
     gto_info => ki(ikind)%orb_basis_set
     
     DO iset = 1, gto_info%nset
     
        nshell = gto_info%nshell(iset)

        DO ishell=1, nshell

           igrid = igrid+1
   
           SELECT CASE(gto_info%l(ishell,iset))
     
           CASE(0)
   
              igt=igt+1 

              dg(igt) % grid_index = igrid
              dg(igt)% dg_rho0 % type  = "S"
              dg(igt)% dg_rho0 % grid = igrid 
              dg(igt)% dg_rho0 % kind = ikind 
              dg(igt)% dg_rho0 % cutoff_radius = gto_info%set_radius(iset)
              dg(igt)% nparts = size (ki(ikind)%atom_list)

              ncc=SIZE(gto_info % gcc(:,ishell,iset))

              ALLOCATE(dg(igt)%dg_rho0%gcc(ncc), STAT=ios)
              IF (ios /= 0 ) THEN
                 CALL stop_memory ('pol_grids_map', 'gcc', ncc)
              END IF

              ALLOCATE(dg(igt)%dg_rho0%zet(ncc), STAT=ios)
              IF (ios /= 0 ) THEN
                 CALL stop_memory ('pol_grids_map', 'zet', ncc)
              END IF 

              DO icc=1, ncc
                dg(igt)% dg_rho0 % zet(icc) = gto_info%zet(icc,iset)
                dg(igt)% dg_rho0 % gcc(icc) = gto_info%gcc(icc,ishell,iset)
              ENDDO                                             

           CASE(1)

               DO ip=1,3

                 igt=igt+1
 
                 SELECT CASE(ip)
                 CASE(1)
                  dg(igt)% dg_rho0% type ="PX"
                 CASE(2)
                  dg(igt)% dg_rho0% type ="PY"
                 CASE(3)
                  dg(igt)% dg_rho0% type ="PZ"
                 END SELECT                       

                 dg(igt) % grid_index = igrid
                 dg(igt)% dg_rho0% grid = igrid 
                 dg(igt)% dg_rho0% kind = ikind 
                 dg(igt)% dg_rho0 % cutoff_radius = gto_info%set_radius(iset)
                 dg(igt)% nparts = size (ki(ikind)%atom_list)

                 ncc=SIZE(gto_info % gcc(:,ishell,iset))

                 ALLOCATE(dg(igt)%dg_rho0%gcc(ncc), STAT=ios)
                 IF (ios /= 0 ) THEN
                    CALL stop_memory ('pol_grids_map', 'gcc', ncc)
                 END IF

                 ALLOCATE(dg(igt)%dg_rho0%zet(ncc), STAT=ios)
                 IF (ios /= 0 ) THEN
                    CALL stop_memory ('pol_grids_map', 'zet', ncc)
                 END IF
                  
                 DO icc=1, ncc
                   dg(igt)% dg_rho0 % zet(icc) = gto_info%zet(icc,iset)
                   dg(igt)% dg_rho0 % gcc(icc) = gto_info%gcc(icc,ishell,iset)
                 ENDDO                                             
 
                 END DO

           END SELECT
   
        ENDDO

     ENDDO

  ENDDO                   

 END SUBROUTINE pol_grids_map 
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
     drho_basis_info, rho0_basis_info, pnode_grp,  box, &
     thermo, small_grid, big_grid, iref )
  
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
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo
  TYPE ( pw_grid_type ), INTENT ( INOUT ), DIMENSION (:) :: small_grid
  TYPE ( pw_grid_type ), INTENT ( INOUT ), DIMENSION (:) :: big_grid
  INTEGER, INTENT (OUT) :: iref
  
! Locals
  INTEGER ::  natoms, npts_s ( 3 ), npts_max ( 3 ), iw
  INTEGER :: igt, ikind, igrid, ig, igrid_old
  REAL ( dbl ) :: qi, qj, i, cutoff, cutoff_radius, cutoff_max,  &
                  max_cut0, max_cut
  TYPE (kind_info_type), POINTER :: ki

!------------------------------------------------------------------------------

! fire up the reciprocal space and compute self interaction and
! term from the neutralizing background.

 iw = pnode_grp % scr  

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
       cutoff_max = cutoff_radius
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
       max_cut0 = max_cut
       npts_max = npts_s
       cutoff_max = cutoff_radius
     END IF   

     igrid_old = igrid

  ENDDO

! Second: setup the reference grid

  CALL dg_pol_grid_setup ( box, npts_max, cutoff_max , &
           small_grid ( iref ), big_grid ( iref ), pnode_grp )

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

       IF ( igrid /= iref ) &
       CALL dg_pol_grid_setup ( box, npts_s, cutoff_radius , &
           small_grid ( igrid ), big_grid ( igrid ), &
           pnode_grp, big_grid ( iref ) )
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
       ELSE
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

       IF ( igrid /= iref ) &
       CALL dg_pol_grid_setup ( box, npts_s, cutoff_radius , &
           small_grid ( igrid ), big_grid ( igrid ), &
           pnode_grp, big_grid ( iref ) )
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
       ELSE
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
  SUBROUTINE assign_coefs_to_grids(ao,part,ki,ngt)

  IMPLICIT NONE
  
  TYPE (ao_type), INTENT (INOUT) :: ao
  TYPE (particle_type), INTENT(IN), DIMENSION (:) :: part
  TYPE (kind_info_type), DIMENSION (:), INTENT(IN) :: ki
  INTEGER, INTENT(IN) :: ngt

! locals 
  INTEGER :: ngrid_types, grid_type_identifier, ikind, ii, ipart
  INTEGER :: ngrid_kind, ncgf, icoef, icgf

  ngrid_types = 0
  ngrid_kind = 0
  
  DO ikind = 1 , size ( ki )
   
    ngrid_types = ngrid_types + ngrid_kind
 
    DO ii = 1, ki ( ikind )%natom

      ipart = ki ( ikind )% atom_list ( ii )
      ncgf = ki ( ikind ) % orb_basis_set % ncgf

      DO icgf = 1, ncgf
         grid_type_identifier = ngrid_types + icgf
         icoef  = part ( ipart ) % coef_list( icgf )
         ao % grid_id ( icoef ) = grid_type_identifier
      ENDDO

    END DO

    ngrid_kind = ncgf
 
 END DO                                       

 IF ( grid_type_identifier /= ngt .AND. ngt /= 0 ) THEN
    CALL stop_program ( "assign_coef_to_grid", "inconsistent grid number" )  
 END IF

 END SUBROUTINE assign_coefs_to_grids


!!*****
!******************************************************************************

END MODULE pol_grids

!******************************************************************************
