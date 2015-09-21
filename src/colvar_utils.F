!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
  USE colvar_types,                    ONLY: &
       colvar_counters, colvar_setup, colvar_type, coord_colvar_id, &
       distance_from_path_colvar_id, gyration_colvar_id, mindist_colvar_id, &
       population_colvar_id, reaction_path_colvar_id, rmsd_colvar_id
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE input_constants,                 ONLY: rmsd_all,&
                                             rmsd_list
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: invert_matrix
  USE memory_utilities,                ONLY: reallocate
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
  USE rmsd,                            ONLY: rmsd3
  USE string_utilities,                ONLY: uppercase
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: number_of_colvar,&
            eval_colvar,&
            set_colvars_target,&
            get_clv_force,&
            post_process_colvar

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'colvar_utils'

CONTAINS

! *****************************************************************************
!> \brief Gives back the number of colvar defined for a force_eval
!> \param force_env ...
!> \param only_intra_colvar ...
!> \param unique ...
!> \retval ntot ...
!> \author Teodoro Laino 05.2007 [tlaino] - Zurich University
! *****************************************************************************
  FUNCTION number_of_colvar(force_env, only_intra_colvar, unique) RESULT(ntot)
    TYPE(force_env_type), POINTER            :: force_env
    LOGICAL, INTENT(IN), OPTIONAL            :: only_intra_colvar, unique
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
    ntot = 0
    CALL force_env_get(force_env=force_env, subsys=subsys)
    CALL cp_subsys_get(subsys=subsys, molecules_new=molecules,  gci=gci,&
         molecule_kinds_new=molecule_kinds)

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
    CALL timestop(handle)

  END FUNCTION number_of_colvar

! *****************************************************************************
!> \brief Set the value of target for constraints/restraints
!> \param targets ...
!> \param force_env ...
!> \author Teodoro Laino 05.2007 [tlaino] - Zurich University
! *****************************************************************************
  SUBROUTINE set_colvars_target(targets, force_env)
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: targets
    TYPE(force_env_type), POINTER            :: force_env

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
    CALL force_env_get(force_env=force_env, subsys=subsys, cell=cell)
    CALL cp_subsys_get(subsys=subsys, gci=gci, molecule_kinds_new=molecule_kinds)

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
    CALL timestop(handle)

  END SUBROUTINE set_colvars_target

! *****************************************************************************
!> \brief Computes the values of colvars and the Wilson matrix B and its invers A
!> \param force_env ...
!> \param coords ...
!> \param cvalues ...
!> \param Bmatrix ...
!> \param MassI ...
!> \param Amatrix ...
!> \author Teodoro Laino 05.2007 [tlaino] - Zurich University
! *****************************************************************************
  SUBROUTINE eval_colvar ( force_env, coords, cvalues, Bmatrix, MassI, Amatrix)

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

    CHARACTER(LEN=*), PARAMETER :: routineN = 'eval_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ikind, imol, &
                                                n_tot, natom, nkind, &
                                                nmol_per_kind, offset
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
       CPASSERT(check)
       Bmatrix = 0.0_dp
    END IF
    CALL timeset(routineN,handle)
    ALLOCATE(map(SIZE(cvalues)))
    map=HUGE(0) ! init all, since used in a sort, but not all set in parallel.
    CALL force_env_get(force_env=force_env, subsys=subsys, cell=cell)
    n_tot=0
    cvalues = 0.0_dp
    CALL cp_subsys_get(subsys=subsys,&
         particles=particles,&
         molecules_new=molecules,&
         local_molecules_new=local_molecules,&
         gci=gci,&
         molecule_kinds_new=molecule_kinds)

    nkind             =  molecule_kinds % n_els
    particle_set      => particles%els
    molecule_set      => molecules % els
    ! Intramolecular Colvars
    IF (number_of_colvar(force_env, only_intra_colvar=.TRUE.)/=0) THEN
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
                                    Bmatrix, offset, n_tot, map)
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
                              Bmatrix, offset, n_tot, map)
       ENDIF
    END IF
    CPASSERT(n_tot==SIZE(cvalues))
    ! Sort values of Collective Variables according the order of the input
    ! sections
    ALLOCATE(wrk(SIZE(cvalues)))
    ALLOCATE(rwrk(SIZE(cvalues)))
    CALL sort(map, SIZE(map), wrk)
    rwrk = cvalues
    DO i = 1, SIZE(wrk)
       cvalues(i) = rwrk(wrk(i))
    END DO
    ! check and sort on Bmatrix
    IF (PRESENT(Bmatrix)) THEN
       check = n_tot==SIZE(Bmatrix,2)
       CPASSERT(check)
       ALLOCATE(bwrk(SIZE(Bmatrix,1),SIZE(Bmatrix,2)))
       bwrk(:,:) = Bmatrix
       DO i = 1, SIZE(wrk)
          Bmatrix(:,i) = bwrk(:,wrk(i))
       END DO
       DEALLOCATE(bwrk)
    END IF
    DEALLOCATE(rwrk)
    DEALLOCATE(wrk)
    DEALLOCATE(map)
    ! Construction of the Amatrix
    IF (PRESENT(Bmatrix).AND.PRESENT(Amatrix)) THEN
       CPASSERT(ASSOCIATED(Amatrix))
       check = SIZE(Bmatrix,1)==SIZE(Amatrix,2)
       CPASSERT(check)
       check = SIZE(Bmatrix,2)==SIZE(Amatrix,1)
       CPASSERT(check)
       ALLOCATE(Gmatrix(n_tot,n_tot))
       ALLOCATE(Gmatrix_i(n_tot,n_tot))
       Gmatrix(:,:) = MATMUL(TRANSPOSE(Bmatrix),Bmatrix)
       CALL invert_matrix(Gmatrix, Gmatrix_i, inv_error)
       CALL cp_assert(ABS(inv_error)<=1.0E-8_dp,cp_warning_level, cp_assertion_failed, routineP,&
            "Error in inverting the Gmatrix larger than 1.0E-8!"//&
CPSourceFileRef)
       Amatrix = MATMUL(Gmatrix_i,TRANSPOSE(Bmatrix))
       DEALLOCATE(Gmatrix_i)
       DEALLOCATE(Gmatrix)
    END IF
    IF (PRESENT(MassI)) THEN
       natom = SIZE(particle_set)
       CPASSERT(ASSOCIATED(MassI))
       CPASSERT(SIZE(MassI)==natom*3)
       DO i = 1, natom
          MassI((i-1)*3+1) = 1.0_dp/particle_set(i)%atomic_kind%mass
          MassI((i-1)*3+2) = 1.0_dp/particle_set(i)%atomic_kind%mass
          MassI((i-1)*3+3) = 1.0_dp/particle_set(i)%atomic_kind%mass
       END DO
    END IF
    CALL timestop(handle)

  END SUBROUTINE eval_colvar

! *****************************************************************************
!> \brief Computes the offset of the colvar for the specific molecule
!> \param i ...
!> \param molecule_set ...
!> \retval offset ...
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
!> \param molecule ...
!> \param particle_set ...
!> \param coords ...
!> \param cell ...
!> \param cvalues ...
!> \param Bmatrix ...
!> \param offset ...
!> \param n_tot ...
!> \param map ...
!> \author Teodoro Laino 05.2007 [tlaino] - Zurich University
! *****************************************************************************
  SUBROUTINE eval_colv_int( molecule, particle_set, coords, cell, cvalues,&
       Bmatrix, offset, n_tot, map)

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
         coords, cell, cvalues, Bmatrix, offset, n_tot, map)

  END SUBROUTINE eval_colv_int

! *****************************************************************************
!> \brief Computes Intermolecular colvar
!> \param gci ...
!> \param particle_set ...
!> \param coords ...
!> \param cell ...
!> \param cvalues ...
!> \param Bmatrix ...
!> \param offset ...
!> \param n_tot ...
!> \param map ...
!> \author Teodoro Laino 05.2007 [tlaino] - Zurich University
! *****************************************************************************
  SUBROUTINE eval_colv_ext( gci, particle_set, coords, cell, cvalues,&
       Bmatrix, offset, n_tot, map)
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
          coords, cell, cvalues, Bmatrix, offset, n_tot, map)

  END SUBROUTINE eval_colv_ext

! *****************************************************************************
!> \brief Real evaluation of colvar and of the Wilson-Eliashevich Matrix
!>      B_ik : i: internal  coordinates
!>             k: cartesian coordinates
!> \param colv_list ...
!> \param fixd_list ...
!> \param lcolv ...
!> \param particle_set ...
!> \param coords ...
!> \param cell ...
!> \param cvalues ...
!> \param Bmatrix ...
!> \param offset ...
!> \param n_tot ...
!> \param map ...
!> \author Teodoro Laino 05.2007 [tlaino] - Zurich University
! *****************************************************************************
  SUBROUTINE eval_colv_low( colv_list, fixd_list, lcolv, particle_set, coords,&
       cell, cvalues, Bmatrix, offset, n_tot, map)

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
               pos=RESHAPE(coords,(/3,SIZE(particle_set)/)), fixd_list=fixd_list)
       ELSE
          CALL colvar_eval_mol_f(lcolv(iconst)%colvar, cell, particles=particle_set,&
               fixd_list=fixd_list)
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
!> \param force_env ...
!> \param forces ...
!> \param coords ...
!> \param nsize_xyz ...
!> \param nsize_int ...
!> \param cvalues ...
!> \param Mmatrix ...
!> \author Teodoro Laino 05.2007
! *****************************************************************************
  SUBROUTINE get_clv_force(force_env, forces, coords, nsize_xyz, nsize_int, cvalues,&
       Mmatrix)
    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT), OPTIONAL                :: forces, coords
    INTEGER, INTENT(IN)                      :: nsize_xyz, nsize_int
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: cvalues, Mmatrix

    CHARACTER(len=*), PARAMETER :: routineN = 'get_clv_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, k
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: tmp
    REAL(KIND=dp), DIMENSION(:), POINTER     :: MassI, wrk
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: Amatrix, Bmatrix

    failure = .FALSE.
    ALLOCATE(Bmatrix(nsize_xyz,nsize_int))
    ALLOCATE(MassI(nsize_xyz))
    ! Transform gradients if requested
    IF (PRESENT(forces)) THEN
       ALLOCATE(wrk(nsize_int))
       ALLOCATE(Amatrix(nsize_int,nsize_xyz))
       ! Compute the transformation matrices and the invese mass diagonal Matrix
       CALL eval_colvar(force_env, coords, cvalues, Bmatrix, MassI, Amatrix)
       wrk = MATMUL(Amatrix,forces)
       forces              = 0.0_dp
       forces(1:nsize_int) = wrk
       DEALLOCATE(Amatrix)
       DEALLOCATE(wrk)
    ELSE
       ! Compute the transformation matrices and the invese mass diagonal Matrix
       CALL eval_colvar(force_env, coords, cvalues, Bmatrix, MassI)
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
    DEALLOCATE(MassI)
    DEALLOCATE(Bmatrix)
  END SUBROUTINE get_clv_force

! *****************************************************************************
!> \brief Complete the description of the COORDINATION colvar when
!>      defined using KINDS
!> \param colvar ...
!> \param particles ...
!> \par History
!>      1.2009 Fabio Sterpone : Added a part for population
!>     10.2014 Moved out of colvar_types.F [Ole Schuett]
!> \author Teodoro Laino - 07.2007
! *****************************************************************************
  SUBROUTINE post_process_colvar(colvar, particles)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particles

    CHARACTER(len=*), PARAMETER :: routineN = 'post_process_colvar', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=3)                         :: name_kind
    INTEGER                                  :: i, ii, j, natoms, nkinds, &
                                                nr_frame, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    natoms = SIZE(particles)
    IF (colvar%type_id == coord_colvar_id) THEN
       IF (colvar%coord_param%use_kinds_from.OR.colvar%coord_param%use_kinds_to) THEN
          ! Atoms from
          IF (colvar%coord_param%use_kinds_from) THEN
             colvar%coord_param%use_kinds_from = .FALSE.
             nkinds = SIZE(colvar%coord_param%c_kinds_from)
             DO i = 1, natoms
                DO j = 1 , nkinds
                  name_kind = TRIM(particles(i)%atomic_kind%name)
                  CALL uppercase(name_kind)
                   IF (TRIM(colvar%coord_param%c_kinds_from(j))==name_kind) THEN
                      CALL reallocate(colvar%coord_param%i_at_from, 1, colvar%coord_param%n_atoms_from+1)
                      colvar%coord_param%n_atoms_from = colvar%coord_param%n_atoms_from + 1
                      colvar%coord_param%i_at_from(colvar%coord_param%n_atoms_from) = i
                   END IF
                END DO
             END DO
             stat=colvar%coord_param%n_atoms_from
             CPASSERT(stat/=0)
          END IF
          ! Atoms to
          IF (colvar%coord_param%use_kinds_to) THEN
             colvar%coord_param%use_kinds_to = .FALSE.
             nkinds = SIZE(colvar%coord_param%c_kinds_to)
             DO i = 1, natoms
                DO j = 1 , nkinds
                   name_kind = TRIM(particles(i)%atomic_kind%name)
                   CALL uppercase(name_kind)
                   IF (TRIM(colvar%coord_param%c_kinds_to(j))==name_kind) THEN
                      CALL reallocate(colvar%coord_param%i_at_to, 1, colvar%coord_param%n_atoms_to+1)
                      colvar%coord_param%n_atoms_to = colvar%coord_param%n_atoms_to + 1
                      colvar%coord_param%i_at_to(colvar%coord_param%n_atoms_to) = i
                   END IF
                END DO
             END DO
             stat=colvar%coord_param%n_atoms_to
             CPASSERT(stat/=0)
          END IF
          ! Atoms to b
          IF (colvar%coord_param%use_kinds_to_b) THEN
             colvar%coord_param%use_kinds_to_b = .FALSE.
             nkinds = SIZE(colvar%coord_param%c_kinds_to_b)
             DO i = 1, natoms
                DO j = 1 , nkinds
                  name_kind = TRIM(particles(i)%atomic_kind%name)
                  CALL uppercase(name_kind)
                   IF (TRIM(colvar%coord_param%c_kinds_to_b(j))==name_kind) THEN
                      CALL reallocate(colvar%coord_param%i_at_to_b, 1, colvar%coord_param%n_atoms_to_b+1)
                      colvar%coord_param%n_atoms_to_b = colvar%coord_param%n_atoms_to_b + 1
                      colvar%coord_param%i_at_to_b(colvar%coord_param%n_atoms_to_b) = i
                   END IF
                END DO
             END DO
             stat=colvar%coord_param%n_atoms_to_b
             CPASSERT(stat/=0)
          END IF

          ! Setup the colvar
          CALL colvar_setup(colvar)
       END IF
    END IF

    IF (colvar%type_id == mindist_colvar_id) THEN
       IF (colvar%mindist_param%use_kinds_from.OR.colvar%mindist_param%use_kinds_to) THEN
          ! Atoms from
          IF (colvar%mindist_param%use_kinds_from) THEN
             colvar%mindist_param%use_kinds_from = .FALSE.
             nkinds = SIZE(colvar%mindist_param%k_coord_from)
             DO i = 1, natoms
                DO j = 1 , nkinds
                  name_kind = TRIM(particles(i)%atomic_kind%name)
                  CALL uppercase(name_kind)
                   IF (TRIM(colvar%mindist_param%k_coord_from(j))==name_kind) THEN
                      CALL reallocate(colvar%mindist_param%i_coord_from, 1, colvar%mindist_param%n_coord_from+1)
                      colvar%mindist_param%n_coord_from = colvar%mindist_param%n_coord_from + 1
                      colvar%mindist_param%i_coord_from(colvar%mindist_param%n_coord_from) = i
                   END IF
                END DO
             END DO
             stat=colvar%mindist_param%n_coord_from
             CPASSERT(stat/=0)
          END IF
          ! Atoms to
          IF (colvar%mindist_param%use_kinds_to) THEN
             colvar%mindist_param%use_kinds_to = .FALSE.
             nkinds = SIZE(colvar%mindist_param%k_coord_to)
             DO i = 1, natoms
                DO j = 1 , nkinds
                   name_kind = TRIM(particles(i)%atomic_kind%name)
                   CALL uppercase(name_kind)
                   IF (TRIM(colvar%mindist_param%k_coord_to(j))==name_kind) THEN
                      CALL reallocate(colvar%mindist_param%i_coord_to, 1, colvar%mindist_param%n_coord_to+1)
                      colvar%mindist_param%n_coord_to = colvar%mindist_param%n_coord_to + 1
                      colvar%mindist_param%i_coord_to(colvar%mindist_param%n_coord_to) = i
                   END IF
                END DO
             END DO
             stat=colvar%mindist_param%n_coord_to
             CPASSERT(stat/=0)
          END IF
          ! Setup the colvar
          CALL colvar_setup(colvar)
       END IF
    END IF

    IF (colvar%type_id == population_colvar_id) THEN

       IF (colvar%population_param%use_kinds_from.OR.colvar%population_param%use_kinds_to) THEN
          ! Atoms from
          IF (colvar%population_param%use_kinds_from) THEN
             colvar%population_param%use_kinds_from = .FALSE.
             nkinds = SIZE(colvar%population_param%c_kinds_from)
             DO i = 1, natoms
                DO j = 1 , nkinds
                  name_kind = TRIM(particles(i)%atomic_kind%name)
                  CALL uppercase(name_kind)
                   IF (TRIM(colvar%population_param%c_kinds_from(j))==name_kind) THEN
                      CALL reallocate(colvar%population_param%i_at_from, 1, colvar%population_param%n_atoms_from+1)
                      colvar%population_param%n_atoms_from = colvar%population_param%n_atoms_from + 1
                      colvar%population_param%i_at_from(colvar%population_param%n_atoms_from) = i
                   END IF
                END DO
             END DO
             stat=colvar%population_param%n_atoms_from
             CPASSERT(stat/=0)
          END IF
          ! Atoms to
          IF (colvar%population_param%use_kinds_to) THEN
             colvar%population_param%use_kinds_to = .FALSE.
             nkinds = SIZE(colvar%population_param%c_kinds_to)
             DO i = 1, natoms
                DO j = 1 , nkinds
                  name_kind = TRIM(particles(i)%atomic_kind%name)
                  CALL uppercase(name_kind)
                   IF (TRIM(colvar%population_param%c_kinds_to(j))==name_kind) THEN
                      CALL reallocate(colvar%population_param%i_at_to, 1, colvar%population_param%n_atoms_to+1)
                      colvar%population_param%n_atoms_to = colvar%population_param%n_atoms_to + 1
                      colvar%population_param%i_at_to(colvar%population_param%n_atoms_to) = i
                   END IF
                END DO
             END DO
             stat=colvar%population_param%n_atoms_to
             CPASSERT(stat/=0)
          END IF
          ! Setup the colvar
          CALL colvar_setup(colvar)
       END IF

    END IF

    IF (colvar%type_id == gyration_colvar_id) THEN

       IF (colvar%gyration_param%use_kinds) THEN
          ! Atoms from
          IF (colvar%gyration_param%use_kinds) THEN
             colvar%gyration_param%use_kinds = .FALSE.
             nkinds = SIZE(colvar%gyration_param%c_kinds)
             DO i = 1, natoms
                DO j = 1 , nkinds
                  name_kind = TRIM(particles(i)%atomic_kind%name)
                  CALL uppercase(name_kind)
                   IF (TRIM(colvar%gyration_param%c_kinds(j))==name_kind) THEN
                      CALL reallocate(colvar%gyration_param%i_at, 1, colvar%gyration_param%n_atoms+1)
                      colvar%gyration_param%n_atoms = colvar%gyration_param%n_atoms + 1
                      colvar%gyration_param%i_at(colvar%gyration_param%n_atoms) = i
                   END IF
                END DO
             END DO
             stat=colvar%gyration_param%n_atoms
             CPASSERT(stat/=0)
          END IF
          ! Setup the colvar
          CALL colvar_setup(colvar)
       END IF
    END IF

    IF (colvar%type_id == rmsd_colvar_id) THEN
      IF(colvar%rmsd_param%subset==rmsd_all .OR. colvar%rmsd_param%subset==rmsd_list) THEN
         ! weights are masses
         DO i=1,SIZE(colvar%rmsd_param%i_rmsd)
           ii = colvar%rmsd_param%i_rmsd(i)
           colvar%rmsd_param%weights(ii) = particles(ii)%atomic_kind%mass
         END DO
      END IF

      IF(colvar%rmsd_param%align_frames) THEN
        nr_frame = SIZE(colvar%rmsd_param%r_ref,2)
        DO i = 2,nr_frame
          CALL rmsd3(particles, colvar%rmsd_param%r_ref(:,i), colvar%rmsd_param%r_ref(:,1), -1,&
                     rotate=.TRUE.)
        END DO
      END IF

    END IF

    IF (colvar%type_id == distance_from_path_colvar_id .OR. colvar%type_id == reaction_path_colvar_id) THEN
      IF(colvar%reaction_path_param%dist_rmsd .OR. colvar%reaction_path_param%rmsd) THEN
         IF(colvar%reaction_path_param%align_frames) THEN
           nr_frame = colvar%reaction_path_param%nr_frames
           DO i = 2,nr_frame
             CALL rmsd3(particles, colvar%reaction_path_param%r_ref(:,i), colvar%reaction_path_param%r_ref(:,1), -1,&
                     rotate=.TRUE.)
            END DO
         END IF
      END IF
    END IF

  END SUBROUTINE post_process_colvar

END MODULE colvar_utils
