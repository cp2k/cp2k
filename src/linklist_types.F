!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/linklist_types [1.0] *
!!
!!   NAME
!!     linklist_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     HAF
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE linklist_types
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: linklist_internal_data_type, set_linklist_internal_data, &
            get_linklist_internal_data, initialize_linklist_data, &
            subcell_data_type

  TYPE subcell_data_type
    INTEGER :: ncell ( 3 )
    INTEGER :: ngrid ( 3 )
    REAL ( dbl ) :: sab_max ( 3 )
    REAL ( dbl ), DIMENSION ( :, :, :, : ), POINTER :: grid_min
    REAL ( dbl ), DIMENSION ( :, :, :, : ), POINTER :: grid_max
    INTEGER, DIMENSION ( :, :, : ), POINTER :: nijk
    INTEGER, DIMENSION ( :, :, :, : ), POINTER :: ijk
  END TYPE subcell_data_type

  TYPE linklist_internal_data_type
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: nonbonded
     TYPE ( neighbor_list_set_p_type ), DIMENSION ( : ), POINTER :: sab_ppnl
     TYPE ( subcell_data_type ), DIMENSION ( :, : ), POINTER :: ll_cell
     TYPE ( global_environment_type ) :: globenv
     INTEGER :: natom_types
     INTEGER :: counter, last_update, num_update
     INTEGER :: print_level
     REAL ( dbl ), DIMENSION ( :, : ), POINTER :: rlist_cut
     REAL ( dbl ), DIMENSION ( :, : ), POINTER :: rlist_cutsq
     REAL ( dbl ) :: subcells, verlet_skin
     REAL ( dbl ) :: lup, aup
     REAL ( dbl ), DIMENSION ( :, : ), POINTER :: r_last_update
  END TYPE linklist_internal_data_type

!!***
!******************************************************************************
  CONTAINS
!******************************************************************************

SUBROUTINE set_linklist_internal_data ( internal_data, globenv, vskin, &
                                        natype, rcut, subcells, count, &
                                        printlevel )


    TYPE(linklist_internal_data_type), &
      INTENT(INOUT)                          :: internal_data
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    REAL(dbl), INTENT(IN), OPTIONAL          :: vskin
    INTEGER, INTENT(IN), OPTIONAL            :: natype
    REAL(dbl), INTENT(IN), OPTIONAL          :: rcut( :, : ), subcells
    INTEGER, INTENT(IN), OPTIONAL            :: count, printlevel

    INTEGER                                  :: iat, isos, jat

!------------------------------------------------------------------------------

  internal_data % globenv = globenv

! Defines the verlet skin
  IF ( PRESENT ( vskin) ) internal_data % verlet_skin = vskin

! Defines the number of atom types
  IF ( PRESENT ( natype ) )  internal_data % natom_types = natype

! Defines the number of atom types
  IF ( PRESENT ( subcells ) )  internal_data % subcells = subcells

! Defines cut-off for particle-particle lists
  IF ( PRESENT ( rcut) ) THEN

! allocating rlist and rlist_cut to upper triangular
     ALLOCATE (internal_data%rlist_cut(internal_data%natom_types, &
               internal_data%natom_types),STAT=isos)
     IF ( isos /= 0 ) CALL stop_memory ( 'set_ll_param', &
          'rlist_cut', internal_data%natom_types ** 2 )

     ALLOCATE ( internal_data%rlist_cutsq ( internal_data%natom_types,  &
                internal_data%natom_types ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'set_ll_param', &
          'rlist_cutsq', internal_data%natom_types ** 2 )
     DO iat = 1, internal_data%natom_types
        DO jat = iat, internal_data%natom_types
           internal_data%rlist_cut(iat,jat) = sqrt(rcut(iat,jat)) + &
              internal_data%verlet_skin
           internal_data%rlist_cutsq(iat,jat) = &
              internal_data%rlist_cut(iat,jat)*internal_data%rlist_cut(iat,jat)
           internal_data%rlist_cut(jat,iat) = internal_data%rlist_cut(iat,jat)
           internal_data%rlist_cutsq(jat,iat) = internal_data%rlist_cutsq(iat,jat)
        END DO
     END DO
  END IF

! Defines the counting variable
  IF ( PRESENT ( count ) ) internal_data%counter = count

! Defines the print level
  IF ( PRESENT ( printlevel ) ) internal_data%print_level = printlevel

END SUBROUTINE set_linklist_internal_data

!******************************************************************************

SUBROUTINE get_linklist_internal_data ( internal_data, last_up, av_up )


    TYPE(linklist_internal_data_type), &
      INTENT(IN)                             :: internal_data
    REAL(dbl), INTENT(OUT)                   :: last_up, av_up

  last_up = internal_data % lup
  av_up = internal_data % aup

END SUBROUTINE get_linklist_internal_data

!******************************************************************************

SUBROUTINE initialize_linklist_data ( internal_data )
    TYPE(linklist_internal_data_type), &
      POINTER                                :: internal_data

    internal_data % natom_types = 0
    internal_data % counter = 0
    internal_data % last_update = 0
    internal_data % num_update = 0
    internal_data % print_level = 1
    NULLIFY ( internal_data % r_last_update )
    NULLIFY ( internal_data % rlist_cut )
    NULLIFY ( internal_data % rlist_cutsq )
    NULLIFY ( internal_data % sab_ppnl )
    NULLIFY ( internal_data % nonbonded )
    NULLIFY ( internal_data % ll_cell )

END SUBROUTINE initialize_linklist_data

!******************************************************************************

END MODULE linklist_types

!******************************************************************************
