!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief evaluations of colvar for internal coordinates schemes
!> \par History
!>      05-2007 created [tlaino]
!> \author Teodoro Laino - Zurich University (2007) [tlaino]
! *****************************************************************************
MODULE colvar_utils
  USE cell_types,                      ONLY: cell_type
  USE colvar_methods,                  ONLY: colvar_eval_mol_f
  USE colvar_types,                    ONLY: colvar_counters
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: invert_matrix
  USE message_passing,                 ONLY: mp_sum
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE molecule_kind_types,             ONLY: colvar_constraint_type,&
                                             fixd_constraint_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             global_constraint_type,&
                                             local_colvar_constraint_type,&
                                             molecule_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: number_of_colvar,&
            eval_colvar,&
            set_colvars_target,&
            get_clv_force

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'colvar_utils'

CONTAINS

! *****************************************************************************
!> \brief Gives back the number of colvar defined for a force_eval
!> \author Teodoro Laino 05.2007 [tlaino] - Zurich University
! *****************************************************************************
  FUNCTION number_of_colvar(force_env, only_intra_colvar, unique, error) RESULT(ntot)
    TYPE(force_env_type), POINTER            :: force_env
    LOGICAL, INTENT(IN), OPTIONAL            :: only_intra_colvar, unique
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: ntot

    CHARACTER(LEN=*), PARAMETER :: routineN = 'number_of_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ikind, imol
    LOGICAL                                  :: failure, my_unique, &
                                                skip_inter_colvar
    TYPE(colvar_counters)                    :: ncolv
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), POINTER        :: molecule_kind, &
                                                molecule_kind_set(:)
    TYPE(molecule_type), POINTER             :: molecule, molecule_set( : )

    failure = .FALSE.
    NULLIFY(subsys, molecules, molecule_kind, molecule, molecule_set, gci)
    CALL timeset(routineN,handle)
    skip_inter_colvar = .FALSE.
    my_unique         = .FALSE.
    IF (PRESENT(only_intra_colvar)) skip_inter_colvar = only_intra_colvar
    IF (PRESENT(unique)) my_unique = unique
    IF (.NOT. failure) THEN
       ntot = 0
       CALL force_env_get(force_env=force_env, subsys=subsys, error=error)
       CALL cp_subsys_get(subsys=subsys, molecules_new=molecules,  gci=gci,&
            molecule_kinds_new=molecule_kinds, error=error)

       molecule_set      => molecules % els
       ! Intramolecular Colvar
       IF (my_unique) THEN
          molecule_kind_set => molecule_kinds % els
          DO ikind = 1, molecule_kinds % n_els
             molecule_kind => molecule_kind_set(ikind)
             CALL get_molecule_kind ( molecule_kind, ncolv=ncolv)
             ntot = ntot + ncolv%ntot
          END DO
       ELSE
          MOL:  DO imol = 1, SIZE(molecule_set)
             molecule => molecule_set ( imol )
             molecule_kind => molecule % molecule_kind

             CALL get_molecule_kind ( molecule_kind,&
                                      ncolv = ncolv )
             ntot = ntot + ncolv%ntot
          END DO MOL
       END IF
       ! Intermolecular Colvar
       IF (.NOT.skip_inter_colvar) THEN
          IF (ASSOCIATED(gci)) THEN
             ntot = ntot + gci%ncolv%ntot
          END IF
       END IF
    END IF
    CALL timestop(handle)

  END FUNCTION number_of_colvar

! *****************************************************************************
!> \brief Set the value of target for constraints/restraints
!> \author Teodoro Laino 05.2007 [tlaino] - Zurich University
! *****************************************************************************
  SUBROUTINE set_colvars_target(targets, force_env, error)
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: targets
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_colvars_target', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ikind, ind, nkind
    LOGICAL                                  :: failure
    TYPE(cell_type), POINTER                 :: cell
    TYPE(colvar_constraint_type), &
      DIMENSION(:), POINTER                  :: colv_list
    TYPE(colvar_counters)                    :: ncolv
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    failure = .FALSE.
    NULLIFY(cell, subsys, molecule_kinds, molecule_kind, gci, colv_list)
    CALL timeset(routineN,handle)
    IF (.NOT. failure) THEN
       CALL force_env_get(force_env=force_env, subsys=subsys, cell=cell, error=error)
       CALL cp_subsys_get(subsys=subsys, gci=gci, molecule_kinds_new=molecule_kinds,&
            error=error)

       nkind  =  molecule_kinds % n_els
       ! Set Target for Intramolecular Colvars
       MOL:  DO ikind = 1, nkind
          molecule_kind => molecule_kinds % els (ikind)
          CALL get_molecule_kind ( molecule_kind,&
                                   colv_list = colv_list,&
                                   ncolv = ncolv )
          IF (ncolv%ntot /= 0) THEN
             DO i = 1, SIZE(colv_list)
                ind = colv_list(i)%inp_seq_num
                colv_list(i)%expected_value = targets(ind)
             END DO
          END IF
       END DO MOL
       ! Set Target for Intermolecular Colvars
       IF (ASSOCIATED(gci)) THEN
          IF ( gci%ncolv%ntot /= 0 ) THEN
             colv_list => gci%colv_list
             DO i = 1, SIZE(colv_list)
                ind = colv_list(i)%inp_seq_num
                colv_list(i)%expected_value = targets(ind)
             END DO
          ENDIF
       END IF
    END IF
    CALL timestop(handle)

  END SUBROUTINE set_colvars_target

! *****************************************************************************
!> \brief Computes the values of colvars and the Wilson matrix B and its invers A
!> \author Teodoro Laino 05.2007 [tlaino] - Zurich University
! *****************************************************************************
  SUBROUTINE eval_colvar ( force_env, coords, cvalues, Bmatrix, MassI, Amatrix, error )

    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: coords
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: cvalues
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: Bmatrix
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: MassI
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: Amatrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'eval_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ikind, imol, &
                                                n_tot, natom, nkind, &
                                                nmol_per_kind, offset, stat
    INTEGER, DIMENSION(:), POINTER           :: map, wrk
    LOGICAL                                  :: check, failure
    REAL(KIND=dp)                            :: inv_error
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: bwrk, Gmatrix, Gmatrix_i
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rwrk
    TYPE(cell_type), POINTER                 :: cell
    TYPE(colvar_counters)                    :: ncolv
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule, molecule_set( : )
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), POINTER             :: particle_set( : )

    failure = .FALSE.
    NULLIFY(cell, subsys, local_molecules, molecule_kinds, &
         molecules, molecule_kind, molecule, &
         molecule_set, particles, particle_set, gci)
    IF (PRESENT(Bmatrix)) THEN
       check = ASSOCIATED(Bmatrix)
       CPPrecondition(check, cp_failure_level, routineP, error, failure)
       Bmatrix = 0.0_dp
    END IF
    CALL timeset(routineN,handle)
    IF (.NOT. failure) THEN
       ALLOCATE(map(SIZE(cvalues)), stat=stat)
       CPPrecondition(stat==0, cp_failure_level, routineP, error, failure)
       CALL force_env_get(force_env=force_env, subsys=subsys, cell=cell, error=error)
       n_tot=0
       cvalues = 0.0_dp
       CALL cp_subsys_get(subsys=subsys,&
            particles=particles,&
            molecules_new=molecules,&
            local_molecules_new=local_molecules,&
            gci=gci,&
            molecule_kinds_new=molecule_kinds,&
            error=error)

       nkind             =  molecule_kinds % n_els
       particle_set      => particles%els
       molecule_set      => molecules % els
       ! Intramolecular Colvars
       IF (number_of_colvar(force_env, only_intra_colvar=.TRUE., error=error)/=0) THEN
          MOL:  DO ikind = 1, nkind
             nmol_per_kind = local_molecules % n_el ( ikind )
             DO imol = 1, nmol_per_kind
                i = local_molecules % list ( ikind ) % array ( imol )
                molecule => molecule_set ( i )
                molecule_kind => molecule % molecule_kind

                CALL get_molecule_kind ( molecule_kind,&
                                         ncolv = ncolv )
                offset = get_colvar_offset(i, molecule_set)
                ! Collective variables
                IF ( ncolv%ntot /= 0 ) THEN
                   CALL eval_colv_int( molecule, particle_set, coords, cell, cvalues,&
                                       Bmatrix, offset, n_tot, map, error )
                ENDIF
             END DO
          END DO MOL
          CALL mp_sum(n_tot,force_env%para_env%group)
          CALL mp_sum(cvalues,force_env%para_env%group)
          IF (PRESENT(Bmatrix)) CALL mp_sum(Bmatrix,force_env%para_env%group)
       END IF
       offset = n_tot
       ! Intermolecular Colvars
       IF (ASSOCIATED(gci)) THEN
          IF ( gci%ncolv%ntot /= 0 ) THEN
             CALL eval_colv_ext( gci, particle_set, coords, cell, cvalues,&
                                 Bmatrix, offset, n_tot, map, error )
          ENDIF
       END IF
       CPPrecondition(n_tot==SIZE(cvalues), cp_failure_level, routineP, error, failure)
       ! Sort values of Collective Variables according the order of the input
       ! sections
       ALLOCATE(wrk(SIZE(cvalues)), stat=stat)
       CPPrecondition(stat==0, cp_failure_level, routineP, error, failure)
       ALLOCATE(rwrk(SIZE(cvalues)), stat=stat)
       CPPrecondition(stat==0, cp_failure_level, routineP, error, failure)
       CALL sort(map, SIZE(map), wrk)
       rwrk = cvalues
       DO i = 1, SIZE(wrk)
          cvalues(i) = rwrk(wrk(i))
       END DO
       ! check and sort on Bmatrix
       IF (PRESENT(Bmatrix)) THEN
          check = n_tot==SIZE(Bmatrix,2)
          CPPrecondition(check, cp_failure_level, routineP, error, failure)
          ALLOCATE(bwrk(SIZE(Bmatrix,1),SIZE(Bmatrix,2)), stat=stat)
          CPPrecondition(stat==0, cp_failure_level, routineP, error, failure)
          bwrk(:,:) = Bmatrix
          DO i = 1, SIZE(wrk)
             Bmatrix(:,i) = bwrk(:,wrk(i))
          END DO
          DEALLOCATE(bwrk, stat=stat)
          CPPrecondition(stat==0, cp_failure_level, routineP, error, failure)
       END IF
       DEALLOCATE(rwrk, stat=stat)
       CPPrecondition(stat==0, cp_failure_level, routineP, error, failure)
       DEALLOCATE(wrk, stat=stat)
       CPPrecondition(stat==0, cp_failure_level, routineP, error, failure)
       DEALLOCATE(map, stat=stat)
       CPPrecondition(stat==0, cp_failure_level, routineP, error, failure)
       ! Construction of the Amatrix
       IF (PRESENT(Bmatrix).AND.PRESENT(Amatrix)) THEN
          CPPrecondition(ASSOCIATED(Amatrix), cp_failure_level, routineP, error, failure)
          check = SIZE(Bmatrix,1)==SIZE(Amatrix,2)
          CPPrecondition(check, cp_failure_level, routineP, error, failure)
          check = SIZE(Bmatrix,2)==SIZE(Amatrix,1)
          CPPrecondition(check, cp_failure_level, routineP, error, failure)
          ALLOCATE(Gmatrix(n_tot,n_tot), stat=stat)
          CPPrecondition(stat==0, cp_failure_level, routineP, error, failure)
          ALLOCATE(Gmatrix_i(n_tot,n_tot), stat=stat)
          CPPrecondition(stat==0, cp_failure_level, routineP, error, failure)
          Gmatrix(:,:) = MATMUL(TRANSPOSE(Bmatrix),Bmatrix)
          CALL invert_matrix(Gmatrix, Gmatrix_i, inv_error, error=error)
          CALL cp_assert(ABS(inv_error)<=1.0E-8_dp,cp_warning_level, cp_assertion_failed, routineP,&
               "Error in inverting the Gmatrix larger than 1.0E-8!"//&
 CPSourceFileRef,&
               error)
          Amatrix = MATMUL(Gmatrix_i,TRANSPOSE(Bmatrix))
          DEALLOCATE(Gmatrix_i, stat=stat)
          CPPrecondition(stat==0, cp_failure_level, routineP, error, failure)
          DEALLOCATE(Gmatrix, stat=stat)
          CPPrecondition(stat==0, cp_failure_level, routineP, error, failure)
       END IF
       IF (PRESENT(MassI)) THEN
          natom = SIZE(particle_set)
          CPPrecondition(ASSOCIATED(MassI), cp_failure_level, routineP, error, failure)
          CPPrecondition(SIZE(MassI)==natom*3, cp_failure_level, routineP, error, failure)
          DO i = 1, natom
             MassI((i-1)*3+1) = 1.0_dp/particle_set(i)%atomic_kind%mass
             MassI((i-1)*3+2) = 1.0_dp/particle_set(i)%atomic_kind%mass
             MassI((i-1)*3+3) = 1.0_dp/particle_set(i)%atomic_kind%mass
          END DO
       END IF
    END IF
    CALL timestop(handle)

  END SUBROUTINE eval_colvar

! *****************************************************************************
!> \brief Computes the offset of the colvar for the specific molecule
!> \author Teodoro Laino 05.2007 [tlaino] - Zurich University
! *****************************************************************************
  FUNCTION get_colvar_offset(i, molecule_set) RESULT(offset)
    INTEGER, INTENT(IN)                      :: i
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    INTEGER                                  :: offset

    INTEGER                                  :: j
    TYPE(colvar_counters)                    :: ncolv
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    offset = 0
    DO j = 1, i-1
       molecule => molecule_set ( j )
       molecule_kind => molecule % molecule_kind
       CALL get_molecule_kind ( molecule_kind,&
                                ncolv = ncolv )
       offset = offset + ncolv%ntot
    END DO

  END FUNCTION get_colvar_offset

! *****************************************************************************
!> \brief Computes Intramolecular colvar
!> \author Teodoro Laino 05.2007 [tlaino] - Zurich University
! *****************************************************************************
  SUBROUTINE eval_colv_int( molecule, particle_set, coords, cell, cvalues,&
       Bmatrix, offset, n_tot, map, error )

    TYPE(molecule_type), POINTER             :: molecule
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: coords
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: cvalues
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: Bmatrix
    INTEGER, INTENT(IN)                      :: offset
    INTEGER, INTENT(INOUT)                   :: n_tot
    INTEGER, DIMENSION(:), POINTER           :: map
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'eval_colv_int', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

    NULLIFY(fixd_list)

    failure = .FALSE.
    molecule_kind => molecule % molecule_kind
    CALL get_molecule_kind ( molecule_kind, colv_list = colv_list, fixd_list=fixd_list )
    CALL get_molecule ( molecule, lcolv=lcolv )
    CALL eval_colv_low( colv_list, fixd_list, lcolv, particle_set,&
         coords, cell, cvalues, Bmatrix, offset, n_tot, map, error)

  END SUBROUTINE eval_colv_int

! *****************************************************************************
!> \brief Computes Intermolecular colvar
!> \author Teodoro Laino 05.2007 [tlaino] - Zurich University
! *****************************************************************************
  SUBROUTINE eval_colv_ext( gci, particle_set, coords, cell, cvalues,&
       Bmatrix, offset, n_tot, map, error )
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: coords
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: cvalues
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: Bmatrix
    INTEGER, INTENT(IN)                      :: offset
    INTEGER, INTENT(INOUT)                   :: n_tot
    INTEGER, DIMENSION(:), POINTER           :: map
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'eval_colv_ext', &
      routineP = moduleN//':'//routineN

    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )

    colv_list => gci%colv_list
    fixd_list => gci%fixd_list
    lcolv => gci%lcolv
    CALL eval_colv_low( colv_list, fixd_list, lcolv, particle_set,&
          coords, cell, cvalues, Bmatrix, offset, n_tot, map, error)

  END SUBROUTINE eval_colv_ext

! *****************************************************************************
!> \brief Real evaluation of colvar and of the Wilson-Eliashevich Matrix
!>      B_ik : i: internal  coordinates
!>             k: cartesian coordinates
!> \author Teodoro Laino 05.2007 [tlaino] - Zurich University
! *****************************************************************************
  SUBROUTINE eval_colv_low( colv_list, fixd_list, lcolv, particle_set, coords,&
       cell, cvalues, Bmatrix, offset, n_tot, map, error)

    TYPE(colvar_constraint_type), POINTER    :: colv_list( : )
    TYPE(fixd_constraint_type), &
      DIMENSION(:), POINTER                  :: fixd_list
    TYPE(local_colvar_constraint_type), &
      POINTER                                :: lcolv( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: coords
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: cvalues
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: Bmatrix
    INTEGER, INTENT(IN)                      :: offset
    INTEGER, INTENT(INOUT)                   :: n_tot
    INTEGER, DIMENSION(:), POINTER           :: map
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'eval_colv_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatm, iconst, ind, ival
    LOGICAL                                  :: failure

    failure = .FALSE.
    ival = offset
    DO iconst = 1, SIZE(colv_list)
       n_tot = n_tot  + 1
       ival  = ival   + 1
       ! Update colvar
       IF (PRESENT(coords)) THEN
          CALL colvar_eval_mol_f(lcolv(iconst)%colvar, cell, particles=particle_set,&
               pos=RESHAPE(coords,(/3,SIZE(particle_set)/)), fixd_list=fixd_list, &
               error=error)
       ELSE
          CALL colvar_eval_mol_f(lcolv(iconst)%colvar, cell, particles=particle_set,&
               fixd_list=fixd_list, error=error)
       END IF
       cvalues(ival) = lcolv(iconst)%colvar%ss
       map(ival)     = colv_list(iconst)%inp_seq_num
       ! Build the Wilson-Eliashevich Matrix
       IF (PRESENT(Bmatrix)) THEN
          DO iatm = 1, SIZE( lcolv (iconst) % colvar % i_atom )
             ind = (lcolv(iconst)%colvar%i_atom(iatm)-1)*3
             Bmatrix(ind+1,ival) = lcolv(iconst)%colvar%dsdr(1,iatm)
             Bmatrix(ind+2,ival) = lcolv(iconst)%colvar%dsdr(2,iatm)
             Bmatrix(ind+3,ival) = lcolv(iconst)%colvar%dsdr(3,iatm)
          END DO
       END IF
    END DO

  END SUBROUTINE eval_colv_low

! *****************************************************************************
!> \brief Computes the forces in the frame of collective variables, and additional
!>        also the local metric tensor
!> \author Teodoro Laino 05.2007
! *****************************************************************************
  SUBROUTINE get_clv_force(force_env, forces, coords, nsize_xyz, nsize_int, cvalues,&
       Mmatrix, error)
    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT), OPTIONAL                :: forces, coords
    INTEGER, INTENT(IN)                      :: nsize_xyz, nsize_int
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: cvalues, Mmatrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_clv_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, k, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: tmp
    REAL(KIND=dp), DIMENSION(:), POINTER     :: MassI, wrk
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: Amatrix, Bmatrix

    failure = .FALSE.
    ALLOCATE(Bmatrix(nsize_xyz,nsize_int),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(MassI(nsize_xyz),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ! Transform gradients if requested
    IF (PRESENT(forces)) THEN
       ALLOCATE(wrk(nsize_int),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(Amatrix(nsize_int,nsize_xyz),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Compute the transformation matrices and the invese mass diagonal Matrix
       CALL eval_colvar(force_env, coords, cvalues, Bmatrix, MassI, Amatrix, error)
       wrk = MATMUL(Amatrix,forces)
       forces              = 0.0_dp
       forces(1:nsize_int) = wrk
       DEALLOCATE(Amatrix, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(wrk, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
       ! Compute the transformation matrices and the invese mass diagonal Matrix
       CALL eval_colvar(force_env, coords, cvalues, Bmatrix, MassI, error=error)
    END IF
    ! Compute the Metric Tensor
    DO i = 1, nsize_int
       DO j = 1, i
          tmp = 0.0_dp
          DO k = 1, nsize_xyz
             tmp = tmp + Bmatrix(k,j)*MassI(k)*Bmatrix(k,i)
          END DO
          Mmatrix((i-1)*nsize_int+j) = tmp
          Mmatrix((j-1)*nsize_int+i) = tmp
       END DO
    END DO
    DEALLOCATE(MassI, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(Bmatrix, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE get_clv_force

END MODULE colvar_utils
