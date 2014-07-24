!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Contains routines useful for the application of constraints during MD
!> \par History
!>      none
! *****************************************************************************
MODULE constraint_util
  USE cell_types,                      ONLY: cell_type
  USE colvar_methods,                  ONLY: colvar_eval_mol_f
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE f77_blas
  USE kinds,                           ONLY: dp,&
                                             int_size
  USE mathlib,                         ONLY: matmul_3x3,&
                                             transpose_3d
  USE message_passing,                 ONLY: mp_sum
  USE molecule_kind_types,             ONLY: colvar_constraint_type,&
                                             fixd_constraint_type,&
                                             g3x3_constraint_type,&
                                             g4x6_constraint_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             global_constraint_type,&
                                             local_colvar_constraint_type,&
                                             local_g3x3_constraint_type,&
                                             local_g4x6_constraint_type,&
                                             molecule_type
  USE particle_types,                  ONLY: particle_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: getold,&
            pv_constraint,&
            check_tol,&
            get_roll_matrix,&
            restore_temporary_set,&
            update_temporary_set

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'constraint_util'

CONTAINS

! *****************************************************************************
!> \brief saves all of the old variables
!> \param gci ...
!> \param local_molecules ...
!> \param molecule_set ...
!> \param molecule_kind_set ...
!> \param particle_set ...
!> \param cell ...
!> \param error ...
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE getold( gci, local_molecules, molecule_set, molecule_kind_set,&
       particle_set , cell, error)

    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'getold', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: first_atom, i, ikind, imol, &
                                                n3x3con, n4x6con, nkind, &
                                                nmol_per_kind
    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(g3x3_constraint_type), POINTER      :: g3x3_list( : )
    TYPE(g4x6_constraint_type), POINTER      :: g4x6_list( : )
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(local_g3x3_constraint_type), &
      POINTER                                :: lg3x3( : )
    TYPE(local_g4x6_constraint_type), &
      POINTER                                :: lg4x6( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    NULLIFY(fixd_list)
    nkind = SIZE ( molecule_kind_set )
    ! Intramolecular constraints
    MOL: DO ikind = 1, nkind
       nmol_per_kind = local_molecules%n_el ( ikind )
       DO imol = 1, nmol_per_kind
          i = local_molecules%list ( ikind )%array ( imol )
          molecule => molecule_set ( i )
          molecule_kind => molecule%molecule_kind
          CALL get_molecule_kind ( molecule_kind, ng3x3=n3x3con, ng4x6=n4x6con, &
               colv_list=colv_list, g3x3_list=g3x3_list, g4x6_list=g4x6_list,&
               fixd_list=fixd_list )
          CALL get_molecule ( molecule, first_atom=first_atom, &
               lcolv=lcolv, lg3x3=lg3x3, lg4x6=lg4x6 )
          CALL getold_low(n3x3con, n4x6con, colv_list, g3x3_list, g4x6_list, fixd_list,&
               lcolv, lg3x3, lg4x6, first_atom, particle_set, cell, error)
       END DO
    END DO MOL
    ! Intermolecular constraints
    IF (gci%ntot>0) THEN
       n3x3con= gci%ng3x3
       n4x6con= gci%ng4x6
       colv_list=> gci%colv_list
       g3x3_list=> gci%g3x3_list
       g4x6_list=> gci%g4x6_list
       fixd_list=> gci%fixd_list
       lcolv=> gci%lcolv
       lg3x3=> gci%lg3x3
       lg4x6=> gci%lg4x6
       CALL getold_low(n3x3con, n4x6con, colv_list, g3x3_list, g4x6_list, fixd_list,&
            lcolv, lg3x3, lg4x6, 1, particle_set, cell, error)
    END IF
  END SUBROUTINE getold

! *****************************************************************************
!> \brief saves all of the old variables - Low Level
!> \param n3x3con ...
!> \param n4x6con ...
!> \param colv_list ...
!> \param g3x3_list ...
!> \param g4x6_list ...
!> \param fixd_list ...
!> \param lcolv ...
!> \param lg3x3 ...
!> \param lg4x6 ...
!> \param first_atom ...
!> \param particle_set ...
!> \param cell ...
!> \param error ...
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE getold_low(n3x3con, n4x6con, colv_list, g3x3_list, g4x6_list, fixd_list,&
       lcolv, lg3x3, lg4x6, first_atom, particle_set, cell, error)

    INTEGER, INTENT(IN)                      :: n3x3con, n4x6con
    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(g3x3_constraint_type), POINTER      :: g3x3_list( : )
    TYPE(g4x6_constraint_type), POINTER      :: g4x6_list( : )
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(local_g3x3_constraint_type), &
      POINTER                                :: lg3x3( : )
    TYPE(local_g4x6_constraint_type), &
      POINTER                                :: lg4x6( : )
    INTEGER, INTENT(IN)                      :: first_atom
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'getold_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iconst, index

    IF (ASSOCIATED(colv_list)) THEN
       ! Collective constraints
       DO iconst = 1, SIZE(colv_list)
          CALL colvar_eval_mol_f(lcolv ( iconst )%colvar_old, cell,&
               particles=particle_set, fixd_list=fixd_list, error=error)
       ENDDO
    END IF
    ! 3x3 constraints
    DO iconst = 1, n3x3con
       index = g3x3_list ( iconst )%a + first_atom - 1
       lg3x3 ( iconst )%ra_old = particle_set ( index )%r
       index = g3x3_list ( iconst )%b + first_atom - 1
       lg3x3 ( iconst )%rb_old = particle_set ( index )%r
       index = g3x3_list ( iconst )%c + first_atom - 1
       lg3x3 ( iconst )%rc_old = particle_set ( index )%r
    ENDDO
    ! 4x6 constraints
    DO iconst = 1, n4x6con
       index = g4x6_list ( iconst )%a + first_atom - 1
       lg4x6 ( iconst )%ra_old = particle_set ( index )%r
       index = g4x6_list ( iconst )%b + first_atom - 1
       lg4x6 ( iconst )%rb_old = particle_set ( index )%r
       index = g4x6_list ( iconst )%c + first_atom - 1
       lg4x6 ( iconst )%rc_old = particle_set ( index )%r
       index = g4x6_list ( iconst )%d + first_atom - 1
       lg4x6 ( iconst )%rd_old = particle_set ( index )%r
    ENDDO

  END SUBROUTINE getold_low

! *****************************************************************************
!> \brief ...
!> \param gci ...
!> \param local_molecules ...
!> \param molecule_set ...
!> \param molecule_kind_set ...
!> \param particle_set ...
!> \param virial ...
!> \param group ...
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE pv_constraint ( gci, local_molecules, molecule_set, molecule_kind_set, &
       particle_set, virial, group )

    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(virial_type), INTENT(INOUT)         :: virial
    INTEGER, INTENT(IN)                      :: group

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pv_constraint', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: first_atom, i, ikind, imol, &
                                                ng3x3, ng4x6, nkind, &
                                                nmol_per_kind
    REAL(KIND=dp)                            :: pv( 3, 3 )
    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(g3x3_constraint_type), POINTER      :: g3x3_list( : )
    TYPE(g4x6_constraint_type), POINTER      :: g4x6_list( : )
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(local_g3x3_constraint_type), &
      POINTER                                :: lg3x3( : )
    TYPE(local_g4x6_constraint_type), &
      POINTER                                :: lg4x6( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    pv = 0.0_dp
    nkind = SIZE ( molecule_kind_set )
    ! Intramolecular Constraints
    MOL: DO ikind = 1, nkind
       nmol_per_kind = local_molecules%n_el ( ikind )
       DO imol = 1, nmol_per_kind
          i = local_molecules%list ( ikind )%array ( imol )
          molecule => molecule_set ( i )
          molecule_kind => molecule%molecule_kind
          CALL get_molecule_kind ( molecule_kind, ng3x3 = ng3x3, &
               ng4x6 = ng4x6, g3x3_list=g3x3_list, g4x6_list=g4x6_list, &
               colv_list = colv_list )
          CALL get_molecule ( molecule, first_atom = first_atom, lg3x3=lg3x3, &
               lg4x6 = lg4x6, lcolv = lcolv )
          CALL pv_constraint_low ( ng3x3, ng4x6, g3x3_list, g4x6_list, colv_list,&
               first_atom, lg3x3, lg4x6, lcolv, particle_set, pv )
       END DO
    END DO MOL
    ! Intermolecular constraints
    IF (gci%ntot>0) THEN
       ng3x3= gci%ng3x3
       ng4x6= gci%ng4x6
       colv_list=> gci%colv_list
       g3x3_list=> gci%g3x3_list
       g4x6_list=> gci%g4x6_list
       lcolv=> gci%lcolv
       lg3x3=> gci%lg3x3
       lg4x6=> gci%lg4x6
       CALL pv_constraint_low ( ng3x3, ng4x6, g3x3_list, g4x6_list, colv_list,&
            1, lg3x3, lg4x6, lcolv, particle_set, pv )
    END IF
    CALL mp_sum ( pv, group )
    ! Symmetrize PV
    pv(1,2) = 0.5_dp*(pv(1,2)+pv(2,1))
    pv(2,1) = pv(1,2)
    pv(1,3) = 0.5_dp*(pv(1,3)+pv(3,1))
    pv(3,1) = pv(1,3)
    pv(3,2) = 0.5_dp*(pv(3,2)+pv(2,3))
    pv(2,3) = pv(3,2)
    ! Store in virial type
    virial%pv_constraint = pv

  END SUBROUTINE pv_constraint

! *****************************************************************************
!> \brief ...
!> \param ng3x3 ...
!> \param ng4x6 ...
!> \param g3x3_list ...
!> \param g4x6_list ...
!> \param colv_list ...
!> \param first_atom ...
!> \param lg3x3 ...
!> \param lg4x6 ...
!> \param lcolv ...
!> \param particle_set ...
!> \param pv ...
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE pv_constraint_low ( ng3x3, ng4x6, g3x3_list, g4x6_list, colv_list,&
       first_atom, lg3x3, lg4x6, lcolv, particle_set, pv )

    INTEGER, INTENT(IN)                      :: ng3x3, ng4x6
    TYPE(g3x3_constraint_type), POINTER      :: g3x3_list( : )
    TYPE(g4x6_constraint_type), POINTER      :: g4x6_list( : )
    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    INTEGER, INTENT(IN)                      :: first_atom
    TYPE(local_g3x3_constraint_type), &
      POINTER                                :: lg3x3( : )
    TYPE(local_g4x6_constraint_type), &
      POINTER                                :: lg4x6( : )
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(INOUT)                          :: pv

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pv_constraint_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iconst, index_a, index_b, &
                                                index_c, index_d
    REAL(KIND=dp)                            :: fc1( 3 ), fc2( 3 ), fc3( 3 ), &
                                                fc4( 3 ), lambda_3x3( 3 ), &
                                                lambda_4x6( 6 )

    IF (ASSOCIATED(colv_list)) THEN
       ! Colvar Constraints
       DO iconst = 1, SIZE(colv_list)
          CALL pv_colv_eval(pv, lcolv ( iconst ), particle_set)
       END DO
    END IF
    ! 3x3
    DO iconst = 1, ng3x3
       !  pv gets updated with FULL multiplier
       lambda_3x3 = lg3x3 ( iconst )%lambda

       fc1= lambda_3x3(1)*lg3x3 ( iconst ) %fa+ &
            lambda_3x3(2)*lg3x3 ( iconst ) %fb
       fc2=-lambda_3x3(1)*lg3x3 ( iconst ) %fa+ &
            lambda_3x3(3)*lg3x3 ( iconst ) %fc
       fc3=-lambda_3x3(2)*lg3x3 ( iconst ) %fb- &
            lambda_3x3(3)*lg3x3 ( iconst ) %fc
       index_a = g3x3_list ( iconst )%a + first_atom -1
       index_b = g3x3_list ( iconst )%b + first_atom -1
       index_c = g3x3_list ( iconst )%c + first_atom -1

       !pv(1,1)
       pv(1,1) = pv(1,1) + fc1(1)*particle_set ( index_a )%r ( 1 )
       pv(1,1) = pv(1,1) + fc2(1)*particle_set ( index_b )%r ( 1 )
       pv(1,1) = pv(1,1) + fc3(1)*particle_set ( index_c )%r ( 1 )
       !pv(1,2)
       pv(1,2) = pv(1,2) + fc1(1)*particle_set ( index_a )%r ( 2 )
       pv(1,2) = pv(1,2) + fc2(1)*particle_set ( index_b )%r ( 2 )
       pv(1,2) = pv(1,2) + fc3(1)*particle_set ( index_c )%r ( 2 )
       !pv(1,3)
       pv(1,3) = pv(1,3) + fc1(1)*particle_set ( index_a )%r ( 3 )
       pv(1,3) = pv(1,3) + fc2(1)*particle_set ( index_b )%r ( 3 )
       pv(1,3) = pv(1,3) + fc3(1)*particle_set ( index_c )%r ( 3 )
       !pv(2,1)
       pv(2,1) = pv(2,1) + fc1(2)*particle_set ( index_a )%r ( 1 )
       pv(2,1) = pv(2,1) + fc2(2)*particle_set ( index_b )%r ( 1 )
       pv(2,1) = pv(2,1) + fc3(2)*particle_set ( index_c )%r ( 1 )
       !pv(2,2)
       pv(2,2) = pv(2,2) + fc1(2)*particle_set ( index_a )%r ( 2 )
       pv(2,2) = pv(2,2) + fc2(2)*particle_set ( index_b )%r ( 2 )
       pv(2,2) = pv(2,2) + fc3(2)*particle_set ( index_c )%r ( 2 )
       !pv(2,3)
       pv(2,3) = pv(2,3) + fc1(2)*particle_set ( index_a )%r ( 3 )
       pv(2,3) = pv(2,3) + fc2(2)*particle_set ( index_b )%r ( 3 )
       pv(2,3) = pv(2,3) + fc3(2)*particle_set ( index_c )%r ( 3 )
       !pv(3,1)
       pv(3,1) = pv(3,1) + fc1(3)*particle_set ( index_a )%r ( 1 )
       pv(3,1) = pv(3,1) + fc2(3)*particle_set ( index_b )%r ( 1 )
       pv(3,1) = pv(3,1) + fc3(3)*particle_set ( index_c )%r ( 1 )
       !pv(3,2)
       pv(3,2) = pv(3,2) + fc1(3)*particle_set ( index_a )%r ( 2 )
       pv(3,2) = pv(3,2) + fc2(3)*particle_set ( index_b )%r ( 2 )
       pv(3,2) = pv(3,2) + fc3(3)*particle_set ( index_c )%r ( 2 )
       !pv(3,3)
       pv(3,3) = pv(3,3) + fc1(3)*particle_set ( index_a )%r ( 3 )
       pv(3,3) = pv(3,3) + fc2(3)*particle_set ( index_b )%r ( 3 )
       pv(3,3) = pv(3,3) + fc3(3)*particle_set ( index_c )%r ( 3 )
    END DO

    ! 4x6
    DO iconst = 1, ng4x6
       !  pv gets updated with FULL multiplier
       lambda_4x6 = lg4x6 ( iconst )%lambda

       fc1= lambda_4x6(1)*lg4x6 ( iconst ) %fa+ &
            lambda_4x6(2)*lg4x6 ( iconst ) %fb+ &
            lambda_4x6(3)*lg4x6 ( iconst ) %fc
       fc2=-lambda_4x6(1)*lg4x6 ( iconst ) %fa+ &
            lambda_4x6(4)*lg4x6 ( iconst ) %fd+ &
            lambda_4x6(5)*lg4x6 ( iconst ) %fe
       fc3=-lambda_4x6(2)*lg4x6 ( iconst ) %fb- &
            lambda_4x6(4)*lg4x6 ( iconst ) %fd+ &
            lambda_4x6(6)*lg4x6 ( iconst ) %ff
       fc4=-lambda_4x6(3)*lg4x6 ( iconst ) %fc- &
            lambda_4x6(5)*lg4x6 ( iconst ) %fe- &
            lambda_4x6(6)*lg4x6 ( iconst ) %ff
       index_a = g4x6_list ( iconst )%a + first_atom -1
       index_b = g4x6_list ( iconst )%b + first_atom -1
       index_c = g4x6_list ( iconst )%c + first_atom -1
       index_d = g4x6_list ( iconst )%d + first_atom -1

       !pv(1,1)
       pv(1,1) = pv(1,1) + fc1(1)*particle_set(index_a)%r(1)
       pv(1,1) = pv(1,1) + fc2(1)*particle_set(index_b)%r(1)
       pv(1,1) = pv(1,1) + fc3(1)*particle_set(index_c)%r(1)
       pv(1,1) = pv(1,1) + fc4(1)*particle_set(index_d)%r(1)
       !pv(1,2)
       pv(1,2) = pv(1,2) + fc1(1)*particle_set(index_a)%r(2)
       pv(1,2) = pv(1,2) + fc2(1)*particle_set(index_b)%r(2)
       pv(1,2) = pv(1,2) + fc3(1)*particle_set(index_c)%r(2)
       pv(1,2) = pv(1,2) + fc4(1)*particle_set(index_d)%r(2)
       !pv(1,3)
       pv(1,3) = pv(1,3) + fc1(1)*particle_set(index_a)%r(3)
       pv(1,3) = pv(1,3) + fc2(1)*particle_set(index_b)%r(3)
       pv(1,3) = pv(1,3) + fc3(1)*particle_set(index_c)%r(3)
       pv(1,3) = pv(1,3) + fc4(1)*particle_set(index_d)%r(3)
       !pv(2,1)
       pv(2,1) = pv(2,1) + fc1(2)*particle_set(index_a)%r(1)
       pv(2,1) = pv(2,1) + fc2(2)*particle_set(index_b)%r(1)
       pv(2,1) = pv(2,1) + fc3(2)*particle_set(index_c)%r(1)
       pv(2,1) = pv(2,1) + fc4(2)*particle_set(index_d)%r(1)
       !pv(2,2)
       pv(2,2) = pv(2,2) + fc1(2)*particle_set(index_a)%r(2)
       pv(2,2) = pv(2,2) + fc2(2)*particle_set(index_b)%r(2)
       pv(2,2) = pv(2,2) + fc3(2)*particle_set(index_c)%r(2)
       pv(2,2) = pv(2,2) + fc4(2)*particle_set(index_d)%r(2)
       !pv(2,3)
       pv(2,3) = pv(2,3) + fc1(2)*particle_set(index_a)%r(3)
       pv(2,3) = pv(2,3) + fc2(2)*particle_set(index_b)%r(3)
       pv(2,3) = pv(2,3) + fc3(2)*particle_set(index_c)%r(3)
       pv(2,3) = pv(2,3) + fc4(2)*particle_set(index_d)%r(3)
       !pv(3,1)
       pv(3,1) = pv(3,1) + fc1(3)*particle_set(index_a)%r(1)
       pv(3,1) = pv(3,1) + fc2(3)*particle_set(index_b)%r(1)
       pv(3,1) = pv(3,1) + fc3(3)*particle_set(index_c)%r(1)
       pv(3,1) = pv(3,1) + fc4(3)*particle_set(index_d)%r(1)
       !pv(3,2)
       pv(3,2) = pv(3,2) + fc1(3)*particle_set(index_a)%r(2)
       pv(3,2) = pv(3,2) + fc2(3)*particle_set(index_b)%r(2)
       pv(3,2) = pv(3,2) + fc3(3)*particle_set(index_c)%r(2)
       pv(3,2) = pv(3,2) + fc4(3)*particle_set(index_d)%r(2)
       !pv(3,3)
       pv(3,3) = pv(3,3) + fc1(3)*particle_set(index_a)%r(3)
       pv(3,3) = pv(3,3) + fc2(3)*particle_set(index_b)%r(3)
       pv(3,3) = pv(3,3) + fc3(3)*particle_set(index_c)%r(3)
       pv(3,3) = pv(3,3) + fc4(3)*particle_set(index_d)%r(3)
    END DO

  END SUBROUTINE pv_constraint_low

! *****************************************************************************
!> \brief ...
!> \param pv ...
!> \param lcolv ...
!> \param particle_set ...
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE pv_colv_eval(pv, lcolv, particle_set)
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(INOUT)                          :: pv
    TYPE(local_colvar_constraint_type), &
      INTENT(IN)                             :: lcolv
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    INTEGER                                  :: i, iatm, ind, j
    REAL(KIND=dp)                            :: lambda, tmp
    REAL(KIND=dp), DIMENSION(3)              :: f

    DO iatm = 1, SIZE( lcolv%colvar_old%i_atom )
       ind =   lcolv%colvar_old%i_atom (iatm)
       f   = - lcolv%colvar_old%dsdr(:,iatm)
       !  pv gets updated with FULL multiplier
       lambda = lcolv%lambda
       DO i = 1,3
          tmp = lambda * particle_set(ind)%r(i)
          DO j = 1,3
             pv(j,i) = pv(j,i) + f(j) * tmp
          END DO
       END DO
    END DO

  END SUBROUTINE pv_colv_eval

! *****************************************************************************
!> \brief ...
!> \param roll_tol ...
!> \param iroll ...
!> \param char ...
!> \param matrix ...
!> \param veps ...
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE check_tol(roll_tol,iroll,char,matrix,veps)

    REAL(KIND=dp), INTENT(OUT)               :: roll_tol
    INTEGER, INTENT(INOUT)                   :: iroll
    CHARACTER(LEN=*), INTENT(IN)             :: char
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: matrix, veps

    CHARACTER(LEN=*), PARAMETER :: routineN = 'check_tol', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: local_tol
    REAL(KIND=dp), DIMENSION(3, 3)           :: diff_rattle, diff_shake
    REAL(KIND=dp), DIMENSION(3, 3), SAVE     :: matrix_old, veps_old

    SELECT CASE ( char )
    CASE ( 'SHAKE' )
       IF (iroll==1) THEN
          matrix_old = matrix
          roll_tol = -1.E10_dp
       ELSE
          roll_tol = 0.0_dp
          diff_shake = ABS ( matrix_old - matrix )
          local_tol = MAXVAL ( diff_shake )
          roll_tol = MAX ( roll_tol, local_tol )
          matrix_old = matrix
       END IF
       iroll = iroll + 1
    CASE ( 'RATTLE' )
       IF (iroll==1) THEN
          veps_old = veps
          roll_tol = -1.E+10_dp
       ELSE
          roll_tol = 0.0_dp
          ! compute tolerance on veps
          diff_rattle = ABS ( veps - veps_old )
          local_tol = MAXVAL ( diff_rattle )
          roll_tol = MAX ( roll_tol, local_tol )
          veps_old = veps
       END IF
       iroll = iroll + 1
    END SELECT

  END SUBROUTINE check_tol

! *****************************************************************************
!> \brief ...
!> \param char ...
!> \param r_shake ...
!> \param v_shake ...
!> \param vector_r ...
!> \param vector_v ...
!> \param u ...
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE get_roll_matrix(char,r_shake,v_shake,vector_r,vector_v,u)

    CHARACTER(len=*), INTENT(IN)             :: char
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT), OPTIONAL                  :: r_shake, v_shake
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: vector_r, vector_v
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN), OPTIONAL                   :: u

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_roll_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(KIND=dp), DIMENSION(3, 3)           :: diag

    IF ( PRESENT ( r_shake ) ) r_shake = 0.0_dp
    IF ( PRESENT ( v_shake ) ) v_shake = 0.0_dp
    diag = 0.0_dp

    SELECT CASE ( char )
    CASE ( 'SHAKE' )
       IF ( PRESENT ( u ) .AND. PRESENT ( vector_v ) .AND.  &
            PRESENT ( vector_r ) ) THEN
          diag ( 1, 1 ) = vector_r ( 1 )
          diag ( 2, 2 ) = vector_r ( 2 )
          diag ( 3, 3 ) = vector_r ( 3 )
          r_shake = MATMUL_3X3 ( MATMUL_3X3 ( u, diag ),TRANSPOSE_3D ( u ) )
          diag ( 1, 1 ) = vector_v ( 1 )
          diag ( 2, 2 ) = vector_v ( 2 )
          diag ( 3, 3 ) = vector_v ( 3 )
          v_shake = MATMUL_3X3 ( MATMUL_3X3 ( u, diag ), TRANSPOSE_3D ( u ) )
          diag = MATMUL_3x3 ( r_shake, v_shake )
          r_shake = diag
       ELSEIF ( .NOT. PRESENT ( u ) .AND. PRESENT ( vector_v ) .AND.  &
            PRESENT ( vector_r ) ) THEN
          DO i = 1, 3
             r_shake ( i, i ) = vector_r ( i ) * vector_v ( i )
             v_shake ( i, i ) = vector_v ( i )
          END DO
       ELSE
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Not sufficient parameters")
       END IF
    CASE ( 'RATTLE' )
       IF ( PRESENT ( u ) .AND. PRESENT ( vector_v ) ) THEN
          diag ( 1, 1 ) = vector_v ( 1 )
          diag ( 2, 2 ) = vector_v ( 2 )
          diag ( 3, 3 ) = vector_v ( 3 )
          v_shake = MATMUL_3x3 ( MATMUL_3X3 ( u, diag ), TRANSPOSE_3D ( u ) )
       ELSEIF ( .NOT. PRESENT ( u ) .AND. PRESENT ( vector_v ) ) THEN
          DO i = 1, 3
             v_shake ( i, i ) = vector_v ( i )
          END DO
       ELSE
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Not sufficient parameters")
       END IF
    END SELECT

  END SUBROUTINE get_roll_matrix

! *****************************************************************************
!> \brief ...
!> \param particle_set ...
!> \param local_particles ...
!> \param pos ...
!> \param vel ...
!> \par History
!>      Teodoro Laino [tlaino] 2007
! *****************************************************************************
  SUBROUTINE restore_temporary_set(particle_set, local_particles, pos, vel)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(distribution_1d_type), POINTER      :: local_particles
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: pos, vel

    CHARACTER(LEN=*), PARAMETER :: routineN = 'restore_temporary_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iparticle, iparticle_kind, &
                                                iparticle_local, &
                                                nparticle_local, stat
    LOGICAL                                  :: failure
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: wrk

    failure = .FALSE.
    ALLOCATE (wrk(SIZE(particle_set)),stat=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "wrk",int_size*SIZE(particle_set))
    wrk = .TRUE.
    DO iparticle_kind=1,SIZE(local_particles%n_el)
       nparticle_local = local_particles%n_el(iparticle_kind)
       DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          wrk(iparticle) = .FALSE.
       END DO
    END DO
    IF (PRESENT(vel)) THEN
       DO iparticle=1,SIZE(particle_set)
          IF (wrk(iparticle)) THEN
             vel (:,iparticle) =  0.0_dp
          END IF
       END DO
    END IF
    IF (PRESENT(pos)) THEN
       DO iparticle=1,SIZE(particle_set)
          IF (wrk(iparticle)) THEN
             pos (:,iparticle) =  0.0_dp
          END IF
       END DO
    END IF
    DEALLOCATE (wrk,stat=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"wrk")
  END SUBROUTINE restore_temporary_set

! *****************************************************************************
!> \brief ...
!> \param group ...
!> \param pos ...
!> \param vel ...
!> \par History
!>      Teodoro Laino [tlaino] 2007
! *****************************************************************************
  SUBROUTINE update_temporary_set(group, pos, vel)
    INTEGER, INTENT(IN)                      :: group
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: pos, vel

    IF (PRESENT(pos)) THEN
       CALL mp_sum(pos, group)
    END IF
    IF (PRESENT(vel)) THEN
       CALL mp_sum(vel, group)
    END IF
  END SUBROUTINE update_temporary_set

END MODULE constraint_util
