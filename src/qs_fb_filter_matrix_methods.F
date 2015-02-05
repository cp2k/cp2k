MODULE qs_fb_filter_matrix_methods

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp_dbcsr_interface,              ONLY: &
       array_data, array_i1d_obj, array_new, array_nullify, cp_dbcsr_create, &
       cp_dbcsr_finalize, cp_dbcsr_get_info, cp_dbcsr_get_stored_coordinates, &
       cp_dbcsr_init, cp_dbcsr_put_block, cp_dbcsr_row_block_sizes, &
       cp_dbcsr_type, dbcsr_distribution_obj, dbcsr_type_no_symmetry
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE fermi_utils,                     ONLY: Fermi,&
                                             FermiFixed
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_8
  USE message_passing,                 ONLY: mp_alltoall
  USE particle_types,                  ONLY: particle_type
  USE qs_fb_atomic_halo_types,         ONLY: &
       fb_atomic_halo_atom_global2halo, fb_atomic_halo_create, &
       fb_atomic_halo_get, fb_atomic_halo_list_get, fb_atomic_halo_list_obj, &
       fb_atomic_halo_nullify, fb_atomic_halo_obj, fb_atomic_halo_release, &
       fb_atomic_halo_set
  USE qs_fb_atomic_matrix_methods,     ONLY: fb_atmatrix_calc_size,&
                                             fb_atmatrix_construct
  USE qs_fb_com_tasks_types,           ONLY: &
       TASK_COST, TASK_DEST, TASK_N_RECORDS, TASK_PAIR, TASK_SRC, &
       fb_com_atom_pairs_calc_buffer_sizes, fb_com_atom_pairs_create, &
       fb_com_atom_pairs_decode, fb_com_atom_pairs_get, &
       fb_com_atom_pairs_has_data, fb_com_atom_pairs_init, &
       fb_com_atom_pairs_nullify, fb_com_atom_pairs_obj, &
       fb_com_atom_pairs_release, fb_com_tasks_build_atom_pairs, &
       fb_com_tasks_create, fb_com_tasks_encode_pair, fb_com_tasks_nullify, &
       fb_com_tasks_obj, fb_com_tasks_release, fb_com_tasks_set, &
       fb_com_tasks_transpose_dest_src
  USE qs_fb_trial_fns_types,           ONLY: fb_trial_fns_get,&
                                             fb_trial_fns_obj
  USE string_utilities,                ONLY: compress,&
                                             uppercase
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fb_filter_matrix_methods'

  PUBLIC ::  fb_fltrmat_build, &
             fb_fltrmat_add_blkcol, &
             fb_fltrmat_generate_com_pairs, &
             fb_fltrmat_build_atomic_fltrmat, &
             fb_fltrmat_fermi_dirac_mu, &
             fb_fltrmat_fermi_dirac_ne

CONTAINS


  ! **********************************************************************
  !> \brief Build the filter matrix
  !> \param H_mat : DBCSR system KS matrix
  !> \param S_mat : DBCSR system overlap matrix
  !> \param atomic_halos : list of all local atomic halos, each halo gives
  !>                       one atomic matrix and contributes to one blk
  !>                       col to the filter matrix
  !> \param trial_funs : the trial functions to be used to shrink the
  !>                     size of the new "filtered" basis
  !> \param para_env : cp2k parallel environment
  !> \param particle_set : set of all particles in the system
  !> \param fermi_level : the fermi level used for defining the filter
  !>                      function, which is a Fermi-Dirac distribution
  !>                      function
  !> \param filter_temp : the filter temperature used for defining the
  !>                      filter function
  !> \param name        : name given to the filter matrix
  !> \param filter_mat  : DBCSR format filter matrix
  !> \prarm error       : cp2k error container
  !> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
  ! **********************************************************************
! *****************************************************************************
!> \brief ...
!> \param H_mat ...
!> \param S_mat ...
!> \param atomic_halos ...
!> \param trial_fns ...
!> \param para_env ...
!> \param particle_set ...
!> \param fermi_level ...
!> \param filter_temp ...
!> \param name ...
!> \param filter_mat ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE fb_fltrmat_build(H_mat, &
                              S_mat, &
                              atomic_halos, &
                              trial_fns, &
                              para_env, &
                              particle_set, &
                              fermi_level, &
                              filter_temp, &
                              name, &
                              filter_mat, &
                              error)
    TYPE(cp_dbcsr_type), POINTER             :: H_mat, S_mat
    TYPE(fb_atomic_halo_list_obj), &
      INTENT(in)                             :: atomic_halos
    TYPE(fb_trial_fns_obj), INTENT(in)       :: trial_fns
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), INTENT(in)                :: fermi_level, filter_temp
    CHARACTER(LEN=*), INTENT(in)             :: name
    TYPE(cp_dbcsr_type), POINTER             :: filter_mat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fb_fltrmat_build', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=32)                        :: symmetry_string
    CHARACTER(LEN=default_string_length)     :: name_string
    INTEGER                                  :: handle, iblkcol, ihalo, &
                                                ikind, LT_UNIT, max_nhalos, &
                                                nblkcols_total, nhalos, stat
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size_data, &
                                                dummy_halo_atoms, halo_atoms, &
                                                ntfns
    LOGICAL                                  :: failure, send_data_only
    TYPE(array_i1d_obj)                      :: col_blk_size, row_blk_size
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(dbcsr_distribution_obj)             :: dbcsr_dist
    TYPE(fb_atomic_halo_obj)                 :: dummy_atomic_halo
    TYPE(fb_atomic_halo_obj), DIMENSION(:), &
      POINTER                                :: halos

!BEG:DEBUG:LT:2015/01/20
!END:DEBUG:LT:2015/01/20
!BEG:DEBUG:LT:2015/01/29
!END:DEBUG:LT:2015/01/29
!BEG:DEBUG:LT:2015/01/30

    LT_UNIT = 100 + para_env%mepos
!END:DEBUG:LT:2015/01/30

    CALL timeset(routineN, handle)

    failure = .FALSE.

    NULLIFY(halos, atomic_kind, col_blk_size_data, ntfns, &
            dummy_halo_atoms)
    CALL array_nullify(row_blk_size)
    CALL array_nullify(col_blk_size)
    CALL fb_atomic_halo_nullify(dummy_atomic_halo)

    ! filter_mat must be of a dissassociated status (i.e. brand new)
    CPPrecondition(.NOT.ASSOCIATED(filter_mat), cp_failure_level, routineP, error, failure)

    ! get trial function information
    CALL fb_trial_fns_get(trial_fns=trial_fns, &
                          nfunctions=ntfns, &
                          error=error)

    ! calculate the row_blk_size and col_blk_size arrays for
    ! constructing the filter matrix in DBCSR format
    ! row_blk_size for the filter matrix is the same as H or S
    CALL cp_dbcsr_get_info(H_mat, &
                           nblkcols_total=nblkcols_total, &
                           row_blk_size=row_blk_size, &
                           distribution=dbcsr_dist)
    ALLOCATE(col_blk_size_data(nblkcols_total), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    col_blk_size_data = 0
    DO iblkcol = 1, nblkcols_total
       atomic_kind => particle_set(iblkcol)%atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind, &
                            kind_number=ikind)
       col_blk_size_data(iblkcol) = ntfns(ikind)
    END DO
    ! DO NOT deallocate cbs if gift=.TRUE. as col_blk_sizes will only point to cbs
    CALL array_new(col_blk_size, col_blk_size_data, gift=.TRUE.)
    name_string = name
    CALL compress(name_string)
    CALL uppercase(name_string)
    ! the filter matrix is non-square and is always non-symmetric
    symmetry_string = dbcsr_type_no_symmetry
    ! create empty filter matrix
    ALLOCATE(filter_mat, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    CALL cp_dbcsr_init(filter_mat, error=error)
    CALL cp_dbcsr_create(matrix=filter_mat, &
                         name=name_string, &
                         dist=dbcsr_dist, &
                         matrix_type=symmetry_string, &
                         row_blk_size=row_blk_size, &
                         col_blk_size=col_blk_size, &
                         nblks=0, &
                         nze=0, &
                         error=error)

!BEG:DEBUG:LT:2015/01/30
    ! write (LT_UNIT,*) "|LT| fb_fltrmat_build: atomic_halos associated = ", &
    !       fb_atomic_halo_list_has_data(atomic_halos)
!END:DEBUG:LT:2015/01/30

    CALL fb_atomic_halo_list_get(atomic_halos=atomic_halos, &
                                 nhalos=nhalos, &
                                 max_nhalos=max_nhalos, &
                                 halos=halos, &
                                 error=error)

    ! create dummy empty atomic halo
    CALL fb_atomic_halo_create(dummy_atomic_halo, error)
    ALLOCATE(dummy_halo_atoms(0), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    CALL fb_atomic_halo_set(atomic_halo=dummy_atomic_halo, &
                            owner_atom=0, &
                            owner_id_in_halo=0, &
                            natoms=0, &
                            halo_atoms=dummy_halo_atoms, &
                            nelectrons=0, &
                            sorted=.TRUE., &
                            error=error)

    send_data_only = .FALSE.

    DO ihalo = 1, max_nhalos

!BEG:DEBUG:LT:2015/01/29
       ! write (LT_UNIT, *) "|LT| fb_fltrmat_build: nhalos = ", nhalos
       ! write (LT_UNIT, *) "|LT| fb_fltrmat_build: max_nhalos = ", max_nhalos
!END:DEBUG:LT:2015/01/29

       IF (ihalo > nhalos) THEN
          send_data_only = .TRUE.
       END IF

       ! construct the filter matrix block by block
       IF (send_data_only) THEN

!BEG:DEBUG:LT:2015/01/29
          ! write (LT_UNIT, *) "|LT| fb_fltrmat_build: doing dummy calculation"
!END:DEBUG:LT:2015/01/29

          CALL fb_fltrmat_add_blkcol(H_mat, &
                                     S_mat, &
                                     dummy_atomic_halo, &
                                     trial_fns, &
                                     para_env, &
                                     particle_set, &
                                     fermi_level, &
                                     filter_temp, &
                                     filter_mat, &
                                     error)
       ELSE

!BEG:DEBUG:LT:2015/01/20
          ! CALL fb_atomic_halo_get(halos(ihalo), halo_atoms=halo_atoms, error=error)
          ! write (LT_UNIT, *) "|LT| fb_fltrmat_build: halo ", ihalo, " = ", halo_atoms
!END:DEBUG:LT:2015/01/20

          CALL fb_fltrmat_add_blkcol(H_mat, &
                                     S_mat, &
                                     halos(ihalo), &
                                     trial_fns, &
                                     para_env, &
                                     particle_set, &
                                     fermi_level, &
                                     filter_temp, &
                                     filter_mat, &
                                     error)
       END IF ! send_data_only
    END DO

    ! finalise the filter matrix
    CALL cp_dbcsr_finalize(filter_mat, error=error)

    ! cleanup
    CALL fb_atomic_halo_release(dummy_atomic_halo, error)

    CALL timestop(handle)

  END SUBROUTINE fb_fltrmat_build


  ! **********************************************************************
  !> \brief Add a computed blocks in one column to the filter matrix.
  !>        It does not finalise the filter matrix
  !> \param H_mat : DBCSR system KS matrix
  !> \param S_mat : DBCSR system overlap matrix
  !> \param atomic_halo :  the halo that contributes to the blk
  !>                       col of the filter matrix
  !> \param trial_funs : the trial functions to be used to shrink the
  !>                     size of the new "filtered" basis
  !> \param para_env : cp2k parallel environment
  !> \param particle_set : set of all particles in the system
  !> \param fermi_level : the fermi level used for defining the filter
  !>                      function, which is a Fermi-Dirac distribution
  !>                      function
  !> \param filter_temp : the filter temperature used for defining the
  !>                      filter function
  !> \param filter_mat  : DBCSR format filter matrix
  !> \prarm error       : cp2k error container
  !> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
  ! **********************************************************************
! *****************************************************************************
!> \brief ...
!> \param H_mat ...
!> \param S_mat ...
!> \param atomic_halo ...
!> \param trial_fns ...
!> \param para_env ...
!> \param particle_set ...
!> \param fermi_level ...
!> \param filter_temp ...
!> \param filter_mat ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE fb_fltrmat_add_blkcol(H_mat, &
                                   S_mat, &
                                   atomic_halo, &
                                   trial_fns, &
                                   para_env, &
                                   particle_set, &
                                   fermi_level, &
                                   filter_temp, &
                                   filter_mat, &
                                   error)
    TYPE(cp_dbcsr_type), POINTER             :: H_mat, S_mat
    TYPE(fb_atomic_halo_obj), INTENT(in)     :: atomic_halo
    TYPE(fb_trial_fns_obj), INTENT(in)       :: trial_fns
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), INTENT(in)                :: fermi_level, filter_temp
    TYPE(cp_dbcsr_type), POINTER             :: filter_mat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fb_fltrmat_add_blkcol', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iatom_global, iatom_in_halo, ind, ipair, ipe, itrial, &
      jatom_global, jatom_in_halo, jkind, LT_UNIT, natoms_global, &
      natoms_in_halo, ncols_atmatrix, ncols_blk, nrows_atmatrix, nrows_blk, &
      numprocs, pe, recv_encode, send_encode, stat
    INTEGER(KIND=int_8), DIMENSION(:), &
      POINTER                                :: pairs_recv, pairs_send
    INTEGER, ALLOCATABLE, DIMENSION(:) :: atomic_H_blk_col_start, &
      atomic_H_blk_row_start, atomic_S_blk_col_start, atomic_S_blk_row_start, &
      col_block_size_data, recv_disps, recv_pair_count, recv_pair_disps, &
      recv_sizes, send_disps, send_pair_count, send_pair_disps, send_sizes
    INTEGER, DIMENSION(:), POINTER           :: ntfns, row_block_size_data
    INTEGER, DIMENSION(:, :), POINTER        :: tfns
    LOGICAL                                  :: failure, found
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: recv_buf, send_buf
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: atomic_filter_mat, atomic_H, &
                                                atomic_S
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(fb_com_atom_pairs_obj)              :: com_pairs_recv, com_pairs_send

!BEG:ORIG:LT:2015/01/26
!END:ORIG:LT:2015/01/26
!BEG:REPL:LT:2015/01/26
! REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE, TARGET :: atomic_H, atomic_S
! REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE, TARGET :: atomic_filter_mat
!END:REPL:LT:2015/01/26
!BEG:DEBUG:LT:2015/01/20
! INTEGER :: ii, jj
! CHARACTER(LEN=100) :: fmt
! REAL(KIND=dp), DIMENSION(:), POINTER :: part_mat_pointer_1
! REAL(KIND=dp), DIMENSION(:,:), POINTER :: part_mat_pointer_2
!END:DEBUG:LT:2015/01/20
!BEG:DEBUG:LT:2015/01/29
!END:DEBUG:LT:2015/01/29
!BEG:DEBUG:LT:2015/01/30

    LT_UNIT = 100 + para_env%mepos
!END:DEBUG:LT:2015/01/30

    CALL timeset(routineN, handle)

    failure = .FALSE.

    NULLIFY(pairs_send, pairs_recv, ntfns, tfns, atomic_kind, &
            row_block_size_data)
    CALL fb_com_atom_pairs_nullify(com_pairs_send)
    CALL fb_com_atom_pairs_nullify(com_pairs_recv)


    ! ----------------------------------------------------------------------
    ! Get communication buffers ready
    ! ----------------------------------------------------------------------

    ! generate send and recv atom pairs
    CALL fb_com_atom_pairs_create(com_pairs_send, error)
    CALL fb_com_atom_pairs_create(com_pairs_recv, error)
    CALL fb_fltrmat_generate_com_pairs(filter_mat, &
                                       atomic_halo, &
                                       para_env, &
                                       com_pairs_send, &
                                       com_pairs_recv, &
                                       error)
    CALL fb_com_atom_pairs_get(atom_pairs=com_pairs_send, &
                               natoms_encode=send_encode, &
                               pairs=pairs_send, &
                               error=error)
    CALL fb_com_atom_pairs_get(atom_pairs=com_pairs_recv, &
                               natoms_encode=recv_encode, &
                               pairs=pairs_recv, &
                               error=error)

    ! get para_env info
    numprocs = para_env%num_pe
    ! me = para_env%mepos + 1   ! my process id, starting counting from 1

    ! obtain trail function information
    CALL fb_trial_fns_get(trial_fns=trial_fns, &
                          nfunctions=ntfns, &
                          functions=tfns, &
                          error=error)

!BEG:DEBUG:LT:2015/01/20
    ! write (LT_UNIT, *) "|LT| fb_fltrmat_add_blkcol: ntfns = ", ntfns
    ! DO ipe = 1, SIZE(ntfns)
    !    write (LT_UNIT, *) "|LT| fb_fltrmat_add_blkcol: tfns(kind=ipe) = ", tfns(:,ipe)
    ! END DO
!END:DEBUG:LT:2015/01/20

    ! obtain row and col block size data for filter matrix
    row_block_size_data => array_data(cp_dbcsr_row_block_sizes(H_mat))
    natoms_global = SIZE(particle_set)
    ALLOCATE(col_block_size_data(natoms_global), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DO jatom_global = 1, natoms_global
       atomic_kind => particle_set(jatom_global)%atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind, kind_number=jkind)
       col_block_size_data(jatom_global) = ntfns(jkind)
    END DO

    ! allocate temporary arrays for send
    ALLOCATE(send_sizes(numprocs), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(send_disps(numprocs), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(send_pair_count(numprocs), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(send_pair_disps(numprocs), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ! setup send buffer sizes
    CALL fb_com_atom_pairs_calc_buffer_sizes(com_pairs_send, &
                                             filter_mat, &
                                             numprocs, &
                                             row_block_size_data, &
                                             col_block_size_data, &
                                             send_sizes, &
                                             send_disps, &
                                             send_pair_count, &
                                             send_pair_disps, &
                                             error)
    ! allocate send buffer
    ALLOCATE(send_buf(SUM(send_sizes)), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! allocate temporary array for recv
    ALLOCATE(recv_sizes(numprocs), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(recv_disps(numprocs), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(recv_pair_count(numprocs), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(recv_pair_disps(numprocs), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ! setup recv buffer sizes
    CALL fb_com_atom_pairs_calc_buffer_sizes(com_pairs_recv, &
                                             filter_mat, &
                                             numprocs, &
                                             row_block_size_data, &
                                             col_block_size_data, &
                                             recv_sizes, &
                                             recv_disps, &
                                             recv_pair_count, &
                                             recv_pair_disps, &
                                             error)
    ! allocate recv buffer
    ALLOCATE(recv_buf(SUM(recv_sizes)), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! ----------------------------------------------------------------------
    ! Construct atomic filter matrix for this atomic_halo
    ! ----------------------------------------------------------------------

    CALL fb_atomic_halo_get(atomic_halo=atomic_halo, &
                            natoms=natoms_in_halo, &
                            error=error)

    ! construct atomic matrix for H for atomic_halo
    ALLOCATE(atomic_H_blk_row_start(natoms_in_halo + 1), &
             atomic_H_blk_col_start(natoms_in_halo + 1), &
             STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    CALL fb_atmatrix_calc_size(H_mat, &
                               atomic_halo, &
                               nrows_atmatrix, &
                               ncols_atmatrix, &
                               atomic_H_blk_row_start, &
                               atomic_H_blk_col_start, &
                               error)

!BEG:DEBUG:LT:2015/01/20
    ! write (LT_UNIT, *) "atomic_H_blk_row_start = ", atomic_H_blk_row_start
    ! write (LT_UNIT, *) "atomic_H_blk_col_start = ", atomic_H_blk_col_start
!END:DEBUG:LT:2015/01/20

    ALLOCATE(atomic_H(nrows_atmatrix,ncols_atmatrix), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    CALL fb_atmatrix_construct(H_mat, &
                               atomic_halo, &
                               para_env, &
                               atomic_H, &
                               atomic_H_blk_row_start, &
                               atomic_H_blk_col_start, &
                               error)

 !BEG:DEBUG:LT:2015/01/20
    ! write (LT_UNIT,*) "|LT| fb_fltrmat_add_blkcol: H_mat ="
    ! CALL cp_dbcsr_print(H_mat, unit_nr=LT_UNIT, error=error)
    ! write (LT_UNIT,*) "|LT| fb_fltrmat_add_blkcol: atomic_H ="
    ! write (fmt,*) "(", SIZE(atomic_H,2), "f9.4",")"
    ! DO ipair = 1, SIZE(atomic_H,1)
    !    write (LT_UNIT,fmt=fmt) atomic_H(ipair,:)
    ! END DO
 !END:DEBUG:LT:2015/01/20

    ! construct atomic matrix for S for atomic_halo
    ALLOCATE(atomic_S_blk_row_start(natoms_in_halo + 1), &
             atomic_S_blk_col_start(natoms_in_halo + 1), &
             STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    CALL fb_atmatrix_calc_size(S_mat, &
                               atomic_halo, &
                               nrows_atmatrix, &
                               ncols_atmatrix, &
                               atomic_S_blk_row_start, &
                               atomic_S_blk_col_start, &
                               error)
    ALLOCATE(atomic_S(nrows_atmatrix,ncols_atmatrix), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    CALL fb_atmatrix_construct(S_mat, &
                               atomic_halo, &
                               para_env, &
                               atomic_S, &
                               atomic_S_blk_row_start, &
                               atomic_S_blk_col_start, &
                               error)

 !BEG:DEBUG:LT:2015/01/20
    ! write (LT_UNIT,*) "|LT| fb_fltrmat_add_blkcol: S_mat ="
    ! CALL cp_dbcsr_print(S_mat, unit_nr=LT_UNIT, error=error)
    ! write (LT_UNIT,*) "|LT| fb_fltrmat_add_blkcol: atomic_S ="
    ! write (fmt,*) "(", SIZE(atomic_S,2), "f9.4",")"
    ! DO ipair = 1, SIZE(atomic_S,1)
    !    write (LT_UNIT,fmt=fmt) atomic_S(ipair,:)
    ! END DO
 !END:DEBUG:LT:2015/01/20

    ! construct the atomic filter matrix
    ALLOCATE(atomic_filter_mat(nrows_atmatrix,ncols_atmatrix), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ! calculate atomic filter matrix only if it is non-zero sized
    IF (nrows_atmatrix > 0 .AND. ncols_atmatrix > 0) THEN
       CALL fb_fltrmat_build_atomic_fltrmat(atomic_H, &
                                            atomic_S, &
                                            fermi_level, &
                                            filter_temp, &
                                            atomic_filter_mat, &
                                            error)
!BEG:DEBUG:LT:2015/01/20
       ! write (LT_UNIT, *) "|LT| fb_fltrmat_add_blkcol: after call to fb_fltrmat_build_atomic_filrmat"
       ! write (LT_UNIT, *) "|LT| fb_fltrmat_add_blkcol: atomic_S ="
       ! write (fmt,*) "(", SIZE(atomic_S,2), "f9.4",")"
       ! DO ii = 1, SIZE(atomic_S,1)
       !    write (LT_UNIT,fmt=fmt) atomic_S(ii,:)
       ! END DO
!END:DEBUG:LT:2015/01/20

!BEG:DEBUG:LT:2015/01/20
       ! write (LT_UNIT, *) "|LT| fb_fltrmat_add_blkcol: atomic_filter_mat ="
       ! write (fmt,*) "(", SIZE(atomic_filter_mat,2), "f9.4",")"
       ! DO ipair = 1, SIZE(atomic_filter_mat,1)
       !    write (LT_UNIT,fmt=fmt) atomic_filter_mat(ipair,:)
       ! END DO
!END:DEBUG:LT:2015/01/20

    END IF

    ! ----------------------------------------------------------------------
    ! Construct filter matrix blocks and add to the correct locations
    ! in send_buffer
    ! ----------------------------------------------------------------------

    ! initialise send buffer
    IF (SIZE(send_buf) > 0) send_buf = 0.0_dp
    ! assign values
    DO ipe = 1, numprocs
       send_sizes(ipe) = 0
       DO ipair = 1, send_pair_count(ipe)

!BEG:DEBUG:LT:2015/01/29
          ! write (LT_UNIT, *) "|LT| fb_fltrmat_add_blkcol: ipair = ", ipair
          ! write (LT_UNIT, *) "|LT| fb_fltrmat_add_blkcol: SIZE(pairs_send) = ", SIZE(pairs_send)
          ! write (LT_UNIT, *) "|LT| fb_fltrmat_add_blkcol: ipe = ", ipe
          ! write (LT_UNIT, *) &
          !       "|LT| fb_fltrmat_add_blkcol: pairs_send(send_pair_disps(ipe) + ipair) = ", &
          !       pairs_send(send_pair_disps(ipe) + ipair)
!END:DEBUG:LT:2015/01/29

          CALL fb_com_atom_pairs_decode(pairs_send(send_pair_disps(ipe) + ipair), &
                                        pe, iatom_global, jatom_global, &
                                        send_encode)
          CALL fb_atomic_halo_atom_global2halo(atomic_halo, &
                                               iatom_global, iatom_in_halo, &
                                               found, error)
          CPPostcondition(found, cp_failure_level, routineP, error, failure)
          CALL fb_atomic_halo_atom_global2halo(atomic_halo, &
                                               jatom_global, jatom_in_halo, &
                                               found, error)
          CPPostcondition(found, cp_failure_level, routineP, error, failure)
          atomic_kind => particle_set(jatom_global)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atomic_kind, &
                               kind_number=jkind)
          nrows_blk = row_block_size_data(iatom_global)
          ncols_blk = ntfns(jkind)

          ! do it column-wise one trial function at a time

!BEG:DEBUG:LT:2015/01/20
          ! write (LT_UNIT, *) "iatom_global, jatom_global", iatom_global, jatom_global
          ! write (LT_UNIT, *) "iatom_in_halo, jatom_in_halo", iatom_in_halo, jatom_in_halo
          ! write (LT_UNIT, *) &
          !       "|LT| fb_fltrmat_add_blkcol: part_atomic_filter (iatom_in_halo, row_start)", &
          !       iatom_in_halo, atomic_H_blk_row_start(iatom_in_halo)
          ! part_mat_pointer_2 => atomic_filter_mat(atomic_H_blk_row_start(iatom_in_halo) : &
          !                                         atomic_H_blk_row_start(iatom_in_halo+1)-1, &
          !                                         1:ncols_atmatrix)
          ! write (fmt,*) "(", SIZE(part_mat_pointer_2, 2), "f9.4",")"
          ! DO ii = 1, SIZE(part_mat_pointer_2, 1)
          !    write (LT_UNIT,fmt=fmt) part_mat_pointer_2(ii,:)
          ! END DO
!END:DEBUG:LT:2015/01/20

!BEG:DEBUG:LT:2015/01/20
          ! write (LT_UNIT, *) "|LT| fb_fltrmat_add_blkcol: atomic_S ="
          ! write (fmt,*) "(", SIZE(atomic_S,2), "f9.4",")"
          ! DO ii = 1, SIZE(atomic_S,1)
          !    write (LT_UNIT,fmt=fmt) atomic_S(ii,:)
          ! END DO
!END:DEBUG:LT:2015/01/20

          DO itrial = 1, ntfns(jkind)

!BEG:DEBUG:LT:2015/01/20
             ! write (LT_UNIT, *) "atomic_S_blk_col_start = ", atomic_S_blk_col_start
             ! write (LT_UNIT, *) "part_atomic_s = (jatom_global, itrial)", &
             !                    jatom_global, itrial
             ! part_mat_pointer_1 => atomic_S(1 : nrows_atmatrix, &
             !                                atomic_S_blk_col_start(jatom_in_halo) + &
             !                                tfns(itrial,jkind) - 1)
             ! DO ii = 1, SIZE(part_mat_pointer_1,1)
             !    write (LT_UNIT,fmt="(1f10.4)") part_mat_pointer_1(ii)
             ! END DO
!END:DEBUG:LT:2015/01/20

             ind = send_disps(ipe) + send_sizes(ipe) + (itrial-1) * nrows_blk
!BEG:ORIG:LT:2015/02/03
             CALL dgemv("N",                                           &
                         nrows_blk,                                    &
                         ncols_atmatrix,                               &
                         1.0_dp,                                       &
                         atomic_filter_mat(                            &
                           atomic_H_blk_row_start(iatom_in_halo) :     &
                           atomic_H_blk_row_start(iatom_in_halo+1)-1,  &
                           1 : ncols_atmatrix                          &
                         ),                                            &
                         nrows_blk,                                    &
                         atomic_S(                                     &
                           1 : nrows_atmatrix,                         &
                           atomic_S_blk_col_start(jatom_in_halo) +     &
                           tfns(itrial,jkind) - 1                      &
                         ),                                            &
                         1,                                            &
                         0.0_dp,                                       &
                         send_buf(ind + 1 : ind + nrows_blk),          &
                         1)
!END:ORIG:LT:2015/02/03
!BEG:REPL:LT:2015/02/04
             ! CALL BL_GEMV(TRANS = "N",                                           &
             !              M      = nrows_blk,                                    &
             !              N      = ncols_atmatrix,                               &
             !              ALPHA  = 1.0_dp,                                       &
             !              A      = part_mat_pointer_2,                           &
             !              LDA    = nrows_blk,                                    &
             !              X      = part_mat_pointer_1,                           &
             !              INCX   = 1,                                            &
             !              BETA   = 0.0_dp,                                       &
             !              Y      = send_buf(ind + 1 : ind + nrows_blk),          &
             !              INCY   = 1)
!END:REPL:LT:2015/02/04
          END DO ! itrial

!BEG:DEBUG:LT:2015/02/03
          ! write (LT_UNIT, *) "|LT| fb_fltrmat_add_blkcol: send_sizes = ", send_sizes
          ! write (LT_UNIT, *) &
          !       "|LT| fb_fltrmat_add_blkcol: send_buf (",iatom_global,",",jatom_global,"):"
          ! DO ii = 1, nrows_blk
          !    DO jj = 1, ncols_blk
          !       write (UNIT=LT_UNIT,FMT="(f9.4)",ADVANCE="no") &
          !             send_buf(send_disps(ipe) + send_sizes(ipe) + ii + nrows_blk*(jj-1))
          !    END DO ! jj
          !    write (LT_UNIT,*) ""
          ! END DO ! ii
!END:DEBUG:LT:2015/02/03

          send_sizes(ipe) = send_sizes(ipe) + nrows_blk * ncols_blk
       END DO ! ipair
    END DO  ! ipe

!BEG:DEBUG:LT:2015/02/03
    ! write (LT_UNIT, *) "|LT| fb_fltrmat_add_blkcol: send_disps = ", send_disps
!END:DEBUG:LT:2015/02/03

    DEALLOCATE(atomic_H, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(atomic_H_blk_row_start, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(atomic_S, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(atomic_S_blk_row_start, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(atomic_filter_mat, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! ----------------------------------------------------------------------
    ! Do communication
    ! ----------------------------------------------------------------------

    CALL mp_alltoall(send_buf, send_sizes, send_disps, &
                     recv_buf, recv_sizes, recv_disps, &
                     para_env%group)

    DEALLOCATE(send_buf, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(send_sizes, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(send_disps, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(send_pair_count, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(send_pair_disps, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! ----------------------------------------------------------------------
    ! Unpack the recv buffer and add the blocks to correct parts of
    ! the DBCSR filter matrix
    ! ----------------------------------------------------------------------

    DO ipe = 1, numprocs
       recv_sizes(ipe) = 0
       DO ipair = 1, recv_pair_count(ipe)
          CALL fb_com_atom_pairs_decode(pairs_recv(recv_pair_disps(ipe) + ipair), &
                                        pe, iatom_global, jatom_global, &
                                        recv_encode)
          nrows_blk = row_block_size_data(iatom_global)
          ncols_blk = col_block_size_data(jatom_global)
          ind = recv_disps(ipe) + recv_sizes(ipe)
          CALL cp_dbcsr_put_block(filter_mat, &
                                  iatom_global, jatom_global, &
                                  recv_buf((ind+1) : (ind+nrows_blk*ncols_blk)))

!BEG:DEBUG:LT:2015/02/03
          ! write (LT_UNIT, *) &
          !       "|LT| fb_fltrmat_add_blkcol: recv_buf (",iatom_global,",",jatom_global,"):"
          ! write (LT_UNIT, *) recv_buf((ind+1) : (ind+nrows_blk*ncols_blk))
          ! DO ii = 1, nrows_blk
          !    DO jj = 1, ncols_blk
          !       write (UNIT=LT_UNIT,FMT="(f9.4)",ADVANCE="no") &
          !             recv_buf(ind + ii + nrows_blk*(jj-1))
          !    END DO ! jj
          !    write (LT_UNIT,*) ""
          ! END DO ! ii
!END:DEBUG:LT:2015/02/03

       recv_sizes(ipe) = recv_sizes(ipe) + nrows_blk * ncols_blk
       END DO ! ipair
    END DO ! ipe

    ! cleanup rest of the temporary arrays
    DEALLOCATE(recv_buf, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(recv_sizes, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(recv_pair_count, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(recv_pair_disps, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    CALL fb_com_atom_pairs_release(com_pairs_send, error)
    CALL fb_com_atom_pairs_release(com_pairs_recv, error)

    ! cannot finalise the matrix until all blocks has been added

    CALL timestop(handle)

  END SUBROUTINE fb_fltrmat_add_blkcol


  ! **********************************************************************
  !> \brief generate the list of blocks (atom pairs) to be sent and received
  !>        in order to construct the filter matrix for each atomic halo
  !> \param filter_mat : DBCSR formated filter matrix
  !> \param atomic_halo :  the halo that contributes to a blk
  !>                       col of the filter matrix
  !> \param para_env : cp2k parallel environment
  !> \param atom_pairs_send : list of blocks to be sent
  !> \param atom_pairs_recv : list of blocks to be received
  !> \prarm error       : cp2k error container
  !> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
  ! **********************************************************************
! *****************************************************************************
!> \brief ...
!> \param filter_mat ...
!> \param atomic_halo ...
!> \param para_env ...
!> \param atom_pairs_send ...
!> \param atom_pairs_recv ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE fb_fltrmat_generate_com_pairs(filter_mat, &
                                           atomic_halo, &
                                           para_env, &
                                           atom_pairs_send, &
                                           atom_pairs_recv, &
                                           error)
    TYPE(cp_dbcsr_type), POINTER             :: filter_mat
    TYPE(fb_atomic_halo_obj), INTENT(in)     :: atomic_halo
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(fb_com_atom_pairs_obj), INTENT(out) :: atom_pairs_send, &
                                                atom_pairs_recv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'fb_fltrmat_generate_com_pairs', &
      routineP = moduleN//':'//routineN

    INTEGER :: dest, handle, iatom_global, iatom_in_halo, iatom_stored, &
      itask, jatom_global, jatom_stored, natoms_in_halo, nblkrows_total, &
      ntasks_send, src, stat
    INTEGER(KIND=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks_send
    INTEGER, DIMENSION(:), POINTER           :: halo_atoms
    LOGICAL                                  :: failure, transpose
    TYPE(fb_com_tasks_obj)                   :: com_tasks_recv, com_tasks_send

    CALL timeset(routineN, handle)

    failure = .FALSE.

    NULLIFY(tasks_send)
    CALL fb_com_tasks_nullify(com_tasks_send)
    CALL fb_com_tasks_nullify(com_tasks_recv)

    ! initialise atom_pairs_send and atom_pairs_recv
    IF (fb_com_atom_pairs_has_data(atom_pairs_send)) THEN
       CALL fb_com_atom_pairs_init(atom_pairs_send, error)
    ELSE
       CALL fb_com_atom_pairs_create(atom_pairs_send, error)
    END IF
    IF (fb_com_atom_pairs_has_data(atom_pairs_recv)) THEN
       CALL fb_com_atom_pairs_init(atom_pairs_recv, error)
    ELSE
       CALL fb_com_atom_pairs_create(atom_pairs_recv, error)
    END IF

    ! source is always the local processor
    src = para_env%mepos

    ! The total number of filter matrix blocks each processor is going
    ! to construct equals to the total number of halo atoms in all of
    ! the atomic halos local to the processor. The number of send
    ! tasks will not exceed this. We do one halo (col) at a time, and
    ! each call of this subroutine will only work on one filter matrix
    ! col corresponding to atomic_halo.

    ! The col atom block index for each filter matrix block are the
    ! owner atom of each halo. The row atom block index for each
    ! filter matrix block corresponding to each col are the halo atoms
    ! of the corresponding halos. Filter matrix is non-symmetric: it
    ! is non-square, because the blocks themselves are non-sqaure

    CALL fb_atomic_halo_get(atomic_halo=atomic_halo, &
                            owner_atom=jatom_global, &
                            natoms=natoms_in_halo, &
                            halo_atoms=halo_atoms, &
                            error=error)
    ntasks_send = natoms_in_halo

    ! allocate send tasks
    ALLOCATE(tasks_send(TASK_N_RECORDS,ntasks_send), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! Get the total number of atoms, this can be obtained from the
    ! total number of block rows in the DBCSR filter matrix.  We
    ! assumes that before calling this subroutine, the filter_mat has
    ! already been created and initialised: i.e. using
    ! cp_dbcsr_create_new. Even if the matrix is at the moment empty,
    ! the attribute nblkrows_total is already assigned from the dbcsr
    ! distribution data
    CALL cp_dbcsr_get_info(filter_mat, &
                           nblkrows_total=nblkrows_total)

    ! construct send tasks
    itask = 1
    DO iatom_in_halo = 1, natoms_in_halo
       iatom_global = halo_atoms(iatom_in_halo)
       iatom_stored = iatom_global
       jatom_stored = jatom_global
       transpose = .FALSE.
       ! find where the constructed block of filter matrix belongs to
       CALL cp_dbcsr_get_stored_coordinates(filter_mat, &
                                            iatom_stored, &
                                            jatom_stored, &
                                            transpose=transpose, &
                                            processor=dest)
       ! create the send tasks
       tasks_send(TASK_DEST,itask) = dest
       tasks_send(TASK_SRC,itask) = src
       CALL fb_com_tasks_encode_pair(tasks_send(TASK_PAIR,itask), &
                                     iatom_global, jatom_global, &
                                     nblkrows_total)
       ! calculation of cost not implemented at the moment
       tasks_send(TASK_COST,itask) = 0
       itask = itask + 1
    END DO ! iatom_in_halo

    CALL fb_com_tasks_create(com_tasks_recv, error)
    CALL fb_com_tasks_create(com_tasks_send, error)

    CALL fb_com_tasks_set(com_tasks=com_tasks_send, &
                          task_dim=TASK_N_RECORDS, &
                          ntasks=ntasks_send, &
                          nencode=nblkrows_total, &
                          tasks=tasks_send, &
                          error=error)

    ! genearte the recv task list (tasks_recv) from the send task list
    CALL fb_com_tasks_transpose_dest_src(com_tasks_recv, "<", com_tasks_send, &
                                         para_env, error)

    ! task lists are now complete, now construct the atom_pairs_send
    ! and atom_pairs_recv from the tasks lists
    CALL fb_com_tasks_build_atom_pairs(com_tasks=com_tasks_send, &
                                       atom_pairs=atom_pairs_send, &
                                       natoms_encode=nblkrows_total, &
                                       send_or_recv="send", &
                                       symmetric=.FALSE., &
                                       error=error)
    CALL fb_com_tasks_build_atom_pairs(com_tasks=com_tasks_recv, &
                                       atom_pairs=atom_pairs_recv, &
                                       natoms_encode=nblkrows_total, &
                                       send_or_recv="recv", &
                                       symmetric=.FALSE., &
                                       error=error)

    ! cleanup
    CALL fb_com_tasks_release(com_tasks_recv, error)
    CALL fb_com_tasks_release(com_tasks_send, error)

    CALL timestop(handle)

  END SUBROUTINE fb_fltrmat_generate_com_pairs


  ! **********************************************************************
  !> \brief Build the atomic filter matrix for each atomic halo
  !> \param atomic_H : atomic KS matrix
  !> \param atomic_S : atomic overlap matrix
  !> \param fermi_level : fermi level used to construct the Fermi-Dirac
  !>                      filter function
  !> \param filter_temp : temperature used to construct the Fermi-Dirac
  !>                      filter function
  !> \param atomic_filter_mat : the atomic filter matrix
  !> \prarm error       : cp2k error container
  !> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
  ! **********************************************************************
! *****************************************************************************
!> \brief ...
!> \param atomic_H ...
!> \param atomic_S ...
!> \param fermi_level ...
!> \param filter_temp ...
!> \param atomic_filter_mat ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE fb_fltrmat_build_atomic_fltrmat(atomic_H, &
                                             atomic_S, &
                                             fermi_level, &
                                             filter_temp, &
                                             atomic_filter_mat, &
                                             error)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(in)                             :: atomic_H, atomic_S
    REAL(KIND=dp), INTENT(in)                :: fermi_level, filter_temp
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(out)                            :: atomic_filter_mat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'fb_fltrmat_build_atomic_fltrmat', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ii, info, jj, &
                                                LT_UNIT, mat_dim, stat, &
                                                work_array_size
    LOGICAL                                  :: check_ok, failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eigenvalues, filter_function, &
                                                work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: atomic_S_copy, eigenvectors, &
                                                filtered_eigenvectors

!BEG:DEBUG:LT:2015/01/26
! CHARACTER(LEN=100) :: fmt
! INTEGER :: kk
! REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: inner_products, temp_mat
!END:DEBUG:LT:2015/01/26
!BEG:DEBUG:LT:2015/01/29
!END:DEBUG:LT:2015/01/29
!BEG:DEBUG:LT:2015/01/30

    LT_UNIT = 200
!END:DEBUG:LT:2015/01/30

    CALL timeset(routineN, handle)

    failure = .FALSE.

    ! This subroutine assumes atomic_filter_mat is not zero size, in
    ! other words, it really has to be constructed, instead of just
    ! being a dummy

    check_ok = SIZE(atomic_filter_mat, 1) > 0 .AND. &
               SIZE(atomic_filter_mat, 2) > 0
    CPPostcondition(check_ok, cp_failure_level, routineP, error, failure)

    ! initialise
    atomic_filter_mat = 0.0_dp
    mat_dim = SIZE(atomic_H, 1)

    ! diagonalise using LAPACK
    ALLOCATE(eigenvalues(mat_dim), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ! get optimal work array size
    ALLOCATE(work(1), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ! dsygv will overwrite part of atomic_H and atomic_S, thus need to copy them
    ALLOCATE(atomic_S_copy(SIZE(atomic_S,1),SIZE(atomic_S,2)), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    atomic_S_copy = atomic_S
    ALLOCATE(eigenvectors(SIZE(atomic_H,1),SIZE(atomic_H,2)), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    eigenvectors = atomic_H
    CALL dsygv(1,             &
               'V',           &
               'U',           &
               mat_dim,       &
               eigenvectors,  &
               mat_dim,       &
               atomic_S_copy, &
               mat_dim,       &
               eigenvalues,   &
               work,          &
               -1,            &
               info)
    work_array_size = work(1)
    ! now allocate work array
    DEALLOCATE(work, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(work(work_array_size), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    work = 0.0_dp
    ! do calculation
    atomic_S_copy = atomic_S
    eigenvectors = atomic_H
    CALL dsygv(1,               &
               'V',             &
               'U',             &
               mat_dim,         &
               eigenvectors,    &
               mat_dim,         &
               atomic_S_copy,   &
               mat_dim,         &
               eigenvalues,     &
               work,            &
               work_array_size, &
               info)
    ! check if diagonalisation is successful
    IF (info .NE. 0) THEN
       WRITE (*,*) "DSYGV ERROR MESSAGE: ", info
       CPErrorMessage(cp_failure_level, routineP, "DSYGV failed", error)
       CPPrecondition(.FALSE., cp_failure_level, routineP, error, failure)
    END IF

    DEALLOCATE(work, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(atomic_S_copy, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! first get the filter function
    ALLOCATE(filter_function(mat_dim), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    filter_function = 0.0_dp
    CALL fb_fltrmat_fermi_dirac_mu(filter_function, &
                                   eigenvalues, &
                                   filter_temp, &
                                   fermi_level)

!BEG:DEBUG:LT:2015/01/26
    ! write (LT_UNIT,*) &
    !       "|LT| fb_fltrmat_build_atomic_fltrmat: filter_function = ", filter_function
    ! write (LT_UNIT,*) &
    !       "|LT| fb_fltrmat_build_atomic_fltrmat: eigenvalues = ", eigenvalues
!END:DEBUG:LT:2015/01/26

    DEALLOCATE(eigenvalues, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! atomic_H has the eigenvectors, construct the version of it
    ! filtered through the filter function
    ALLOCATE(filtered_eigenvectors(mat_dim,mat_dim), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DO jj = 1, mat_dim
       DO ii = 1, mat_dim
          filtered_eigenvectors(ii,jj) = &
               filter_function(jj) * eigenvectors(ii,jj)
       END DO ! ii
    END DO ! jj

!BEG:DEBUG:LT:2015/01/26
    ! write (LT_UNIT, *) "|LT| fb_fltrmat_build_atomic_fltrmat: eigenvectors = "
    ! write (fmt,*) "(", SIZE(eigenvectors,2), "f9.4", ")"
    ! DO ii = 1, SIZE(eigenvectors,1)
    !    write (LT_UNIT,fmt=fmt) eigenvectors(ii,:)
    ! END DO
!END:DEBUG:LT:2015/01/26

!BEG:DEBUG:LT:2015/01/26
    ! write (LT_UNIT, *) "|LT| fb_fltrmat_build_atomic_fltrmat: check orthonormality of eigenvectors"
    ! ALLOCATE(inner_products(SIZE(eigenvectors,1),SIZE(eigenvectors,2)), STAT=stat)
    ! inner_products = 0.0_dp
    ! ALLOCATE(temp_mat(SIZE(eigenvectors,1),SIZE(eigenvectors,2)), STAT=stat)
    ! temp_mat = 0.0_dp
    ! DO jj = 1, SIZE(eigenvectors,2)
    !    DO ii = 1, SIZE(eigenvectors,1)
    !       DO kk = 1, SIZE(eigenvectors,2)
    !          temp_mat(ii,jj) = temp_mat(ii,jj) + &
    !                            eigenvectors(kk,ii) * atomic_S(kk,jj)
    !       END DO
    !    END DO
    ! END DO
    ! DO jj = 1, SIZE(eigenvectors,2)
    !    DO ii = 1, SIZE(eigenvectors,1)
    !       DO kk = 1, SIZE(eigenvectors,2)
    !          inner_products(ii,jj) = inner_products(ii,jj) + &
    !                                  temp_mat(ii,kk) * eigenvectors(kk,jj)
    !       END DO
    !    END DO
    ! END DO
    ! write (fmt,*) "(", SIZE(inner_products,2), "f9.4", ")"
    ! DO ii = 1, SIZE(inner_products,1)
    !    write (LT_UNIT,fmt=fmt) inner_products(ii,:)
    ! END DO
    ! DEALLOCATE(inner_products)
    ! DEALLOCATE(temp_mat)
!END:DEBUG:LT:2015/01/26

!BEG:DEBUG:LT:2015/01/26
    ! write (LT_UNIT, *) "|LT| fb_fltrmat_build_atomic_fltrmat: filtered_eigenvalues = "
    ! write (fmt,*) "(", SIZE(filtered_eigenvectors,2), "f9.4", ")"
    ! DO ii = 1, SIZE(filtered_eigenvectors,1)
    !    write (LT_UNIT,fmt=fmt) filtered_eigenvectors(ii,:)
    ! END DO
!END:DEBUG:LT:2015/01/26

    DEALLOCATE(filter_function, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    ! construct atomic filter matrix
    CALL dgemm("N",                   &
               "T",                   &
               mat_dim,               &
               mat_dim,               &
               mat_dim,               &
               1.0_dp,                &
               filtered_eigenvectors, &
               mat_dim,               &
               eigenvectors,          &
               mat_dim,               &
               0.0_dp,                &
               atomic_filter_mat,     &
               mat_dim)

!BEG:DEBUG:LT:2015/01/26
    ! write (LT_UNIT, *) "|LT| fb_fltrmat_build_atomic_fltrmat: atomic_filter_mat = "
    ! write (fmt,*) "(", SIZE(atomic_filter_mat,2), "f9.4", ")"
    ! DO ii = 1, SIZE(atomic_filter_mat,1)
    !    write (LT_UNIT,fmt=fmt) atomic_filter_mat(ii,:)
    ! END DO
!END:DEBUG:LT:2015/01/26

    ! remove small negative terms due to numerical error, the filter
    ! matrix must not be negative definite
    DO jj = 1, SIZE(atomic_filter_mat,2)
       DO ii = 1, SIZE(atomic_filter_mat,1)
          IF (ABS(atomic_filter_mat(ii,jj)) < EPSILON(0.0_dp)) THEN
             atomic_filter_mat(ii,jj) = 0.0_dp
          END IF
       END DO
    END DO

    DEALLOCATE(filtered_eigenvectors, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(eigenvectors, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

    CALL timestop(handle)

  END SUBROUTINE fb_fltrmat_build_atomic_fltrmat


  ! **********************************************************************
  !> \brief get values of Fermi-Dirac distribution based on a given fermi
  !>        level at a given set of energy eigenvalues
  !> \param f : the Fermi-Dirac distribution function values
  !> \param eigenvals : set of energy eigenvalues
  !> \param T : temperature
  !> \param mu : the fermi level
  !> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
  ! **********************************************************************
! *****************************************************************************
!> \brief ...
!> \param f ...
!> \param eigenvals ...
!> \param T ...
!> \param mu ...
! *****************************************************************************
  SUBROUTINE fb_fltrmat_fermi_dirac_mu(f, eigenvals, T, mu)
    REAL(KIND=dp), DIMENSION(:), INTENT(out) :: f
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: eigenvals
    REAL(KIND=dp), INTENT(in)                :: T, mu

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_fltrmat_fermi_dirac_mu', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: kTS, ne

! we want fermi function to normalise to 1, so maxocc = 1 here

    CALL Fermi(f, ne, kTS, eigenvals, mu, T, 1.0_dp)
  END SUBROUTINE fb_fltrmat_fermi_dirac_mu


  ! **********************************************************************
  !> \brief get values of Fermi-Dirac distribution based on a given electron
  !>        number at a given set of energy eigenvales
  !> \param f : the Fermi-Dirac distribution function values
  !> \param eigenvals : set of energy eigenvalues
  !> \param T : temperature
  !> \param ne : number of electrons
  !> \param maxocc : maximum occupancy per orbital
  !> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
  ! **********************************************************************
! *****************************************************************************
!> \brief ...
!> \param f ...
!> \param eigenvals ...
!> \param T ...
!> \param ne ...
!> \param maxocc ...
! *****************************************************************************
  SUBROUTINE fb_fltrmat_fermi_dirac_ne(f, eigenvals, T, ne, maxocc)
    REAL(KIND=dp), DIMENSION(:), INTENT(out) :: f
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: eigenvals
    REAL(KIND=dp), INTENT(in)                :: T, ne, maxocc

    CHARACTER(len=*), PARAMETER :: routineN = 'fb_fltrmat_fermi_dirac_ne', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: kTS, mu

! mu is the calculated fermi level
! kTS is the calculated entropic contribution to the energy i.e. -TS
! kTS= kT*[f ln f + (1-f) ln (1-f)]

    failure = .FALSE.
    CALL FermiFixed(f, mu, kTS, eigenvals, ne, T, maxocc)
  END SUBROUTINE fb_fltrmat_fermi_dirac_ne


END MODULE qs_fb_filter_matrix_methods
