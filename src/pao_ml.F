!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Main module for PAO Machine Learning
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_ml
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_create,&
                                              cell_type
   USE cp_dbcsr_interface,              ONLY: cp_dbcsr_iterator,&
                                              cp_dbcsr_iterator_blocks_left,&
                                              cp_dbcsr_iterator_next_block,&
                                              cp_dbcsr_iterator_start,&
                                              cp_dbcsr_iterator_stop,&
                                              cp_dbcsr_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE machine,                         ONLY: m_flush
   USE message_passing,                 ONLY: mp_bcast,&
                                              mp_sum
   USE pao_input,                       ONLY: id2str,&
                                              pao_ml_gp,&
                                              pao_ml_lazy,&
                                              pao_ml_nn,&
                                              pao_ml_prior_mean,&
                                              pao_ml_prior_zero,&
                                              pao_rotinv_param
   USE pao_io,                          ONLY: pao_ioblock_type,&
                                              pao_iokind_type,&
                                              pao_kinds_ensure_equal,&
                                              pao_read_raw
   USE pao_ml_descriptor,               ONLY: pao_ml_calc_descriptor
   USE pao_ml_gaussprocess,             ONLY: pao_ml_gp_gradient,&
                                              pao_ml_gp_predict,&
                                              pao_ml_gp_train
   USE pao_ml_neuralnet,                ONLY: pao_ml_nn_gradient,&
                                              pao_ml_nn_predict,&
                                              pao_ml_nn_train
   USE pao_types,                       ONLY: pao_env_type,&
                                              training_matrix_type
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_ml'

   PUBLIC :: pao_ml_init, pao_ml_predict, pao_ml_forces

   ! linked list used to group training points by kind
   TYPE training_point_type
      TYPE(training_point_type), POINTER       :: next => Null()
      REAL(dp), DIMENSION(:), ALLOCATABLE      :: input
      REAL(dp), DIMENSION(:), ALLOCATABLE      :: output
   END TYPE training_point_type

   TYPE training_list_type
      CHARACTER(LEN=default_string_length)     :: kindname = ""
      TYPE(training_point_type), POINTER       :: head => Null()
      INTEGER                                  :: npoints = 0
   END TYPE training_list_type

CONTAINS

! **************************************************************************************************
!> \brief Initializes the learning machinery
!> \param pao ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE pao_ml_init(pao, qs_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_ml_init', routineP = moduleN//':'//routineN

      INTEGER                                            :: i
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(training_list_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: training_lists

      IF (SIZE(pao%ml_training_set) == 0) RETURN

      IF (pao%iw > 0) WRITE (pao%iw, *) 'PAO|ML| Initializing maschine learning...'

      IF (pao%parameterization /= pao_rotinv_param) &
         CPABORT("PAO maschine learning requires ROTINV parametrization")

      CALL get_qs_env(qs_env, para_env=para_env, atomic_kind_set=atomic_kind_set)

      ! create training-set data-structure
      ALLOCATE (training_lists(SIZE(atomic_kind_set)))
      DO i = 1, SIZE(training_lists)
         CALL get_atomic_kind(atomic_kind_set(i), name=training_lists(i)%kindname)
      ENDDO

      ! parses training files, calculates descriptors and stores all training-points as linked lists
      DO i = 1, SIZE(pao%ml_training_set)
         CALL add_to_training_list(pao, qs_env, training_lists, filename=pao%ml_training_set(i)%fn)
      ENDDO

      ! ensure there there are training points for all kinds that use pao
      CALL sanity_check(qs_env, training_lists)

      ! turns linked lists into matrices and syncs them across ranks
      CALL training_list2matrix(training_lists, pao%ml_training_matrices, para_env)

      ! calculate and substract prior
      CALL pao_ml_substract_prior(pao%ml_prior, pao%ml_training_matrices)

      ! print some statistics about the training set and dump it upon request
      CALL pao_ml_print(pao, pao%ml_training_matrices)

      ! use training-set to train model
      CALL pao_ml_train(pao)

   END SUBROUTINE pao_ml_init

! **************************************************************************************************
!> \brief Reads the given file and adds it's training points to linked lists.
!> \param pao ...
!> \param qs_env ...
!> \param training_lists ...
!> \param filename ...
! **************************************************************************************************
   SUBROUTINE add_to_training_list(pao, qs_env, training_lists, filename)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(training_list_type), DIMENSION(:)             :: training_lists
      CHARACTER(LEN=default_path_length)                 :: filename

      CHARACTER(len=*), PARAMETER :: routineN = 'add_to_training_list', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: param
      INTEGER                                            :: iatom, ikind, natoms, nkinds, nparams
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom2kind, kindsmap
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: hmat, positions
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(pao_ioblock_type), ALLOCATABLE, DIMENSION(:)  :: xblocks
      TYPE(pao_iokind_type), ALLOCATABLE, DIMENSION(:)   :: kinds
      TYPE(particle_type), DIMENSION(:), POINTER         :: my_particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(training_point_type), POINTER                 :: new_point

      NULLIFY (new_point, cell)

      IF (pao%iw > 0) WRITE (pao%iw, '(A,A)') " PAO|ML| Reading training frame from file: ", TRIM(filename)

      CALL get_qs_env(qs_env, para_env=para_env)

      ! parse training data on first rank
      IF (para_env%mepos == para_env%source) THEN
         CALL pao_read_raw(filename, param, hmat, kinds, atom2kind, positions, xblocks)

         ! check parametrization
         IF (TRIM(param) .NE. TRIM(ADJUSTL(id2str(pao%parameterization)))) &
            CPABORT("Restart PAO parametrization does not match")

         ! map read-in kinds onto kinds of this run
         CALL match_kinds(pao, qs_env, kinds, kindsmap)
         nkinds = SIZE(kindsmap)
         natoms = SIZE(positions, 1)
      ENDIF

      ! broadcast parsed raw training data
      CALL mp_bcast(nkinds, para_env%source, para_env%group)
      CALL mp_bcast(natoms, para_env%source, para_env%group)
      IF (para_env%mepos /= para_env%source) THEN
         ALLOCATE (hmat(3, 3))
         ALLOCATE (kindsmap(nkinds))
         ALLOCATE (positions(natoms, 3))
         ALLOCATE (atom2kind(natoms))
      ENDIF
      CALL mp_bcast(hmat, para_env%source, para_env%group)
      CALL mp_bcast(kindsmap, para_env%source, para_env%group)
      CALL mp_bcast(atom2kind, para_env%source, para_env%group)
      CALL mp_bcast(positions, para_env%source, para_env%group)

      ! create cell from read-in h-matrix
      CALL cell_create(cell, hmat)

      ! create a particle_set based on read-in positions and refere to kinds of this run
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set)
      ALLOCATE (my_particle_set(natoms))
      DO iatom = 1, natoms
         ikind = kindsmap(atom2kind(iatom))
         my_particle_set(iatom)%atomic_kind => atomic_kind_set(ikind)
         my_particle_set(iatom)%r = positions(iatom, :)
      ENDDO

      ! fill linked list with training points
      ! Afterwards all ranks will have lists with the same number of entries,
      ! however the input and output arrays will only be allocated on one rank per entry.
      ! We farm out the expensive calculation of the descriptor across ranks.
      DO iatom = 1, natoms
         ALLOCATE (new_point)

         ! training-point input, calculate descriptor only on one rank
         IF (MOD(iatom-1, para_env%num_pe) == para_env%mepos) THEN
            CALL pao_ml_calc_descriptor(pao, &
                                        my_particle_set, &
                                        qs_kind_set, &
                                        cell, &
                                        iatom=iatom, &
                                        descriptor=new_point%input)
         ENDIF

         ! copy training-point output on first rank
         IF (para_env%mepos == para_env%source) THEN
            nparams = SIZE(xblocks(iatom)%p, 1)
            ALLOCATE (new_point%output(nparams))
            new_point%output(:) = xblocks(iatom)%p(:, 1)
         ENDIF

         ! add to linked list
         ikind = kindsmap(atom2kind(iatom))
         training_lists(ikind)%npoints = training_lists(ikind)%npoints+1
         new_point%next => training_lists(ikind)%head
         training_lists(ikind)%head => new_point
      ENDDO

      DEALLOCATE (cell, my_particle_set, hmat, kindsmap, positions, atom2kind)

   END SUBROUTINE add_to_training_list

! **************************************************************************************************
!> \brief Make read-in kinds on to atomic-kinds of this run
!> \param pao ...
!> \param qs_env ...
!> \param kinds ...
!> \param kindsmap ...
! **************************************************************************************************
   SUBROUTINE match_kinds(pao, qs_env, kinds, kindsmap)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pao_iokind_type), DIMENSION(:)                :: kinds
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kindsmap

      CHARACTER(LEN=default_string_length)               :: name
      INTEGER                                            :: ikind, jkind
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set

      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set)

      CPASSERT(.NOT. ALLOCATED(kindsmap))
      ALLOCATE (kindsmap(SIZE(kinds)))
      kindsmap(:) = -1

      DO ikind = 1, SIZE(kinds)
         DO jkind = 1, SIZE(atomic_kind_set)
            CALL get_atomic_kind(atomic_kind_set(jkind), name=name)
            ! match kinds via their name
            IF (TRIM(kinds(ikind)%name) .EQ. TRIM(name)) THEN
               CALL pao_kinds_ensure_equal(pao, qs_env, jkind, kinds(ikind))
               kindsmap(ikind) = jkind
               EXIT
            ENDIF
         ENDDO
      ENDDO

      IF (ANY(kindsmap < 1)) &
         CPABORT("PAO: Could not match all kinds from training set")
   END SUBROUTINE match_kinds

! **************************************************************************************************
!> \brief Checks that there is at least one training point per pao-enabled kind
!> \param qs_env ...
!> \param training_lists ...
! **************************************************************************************************
   SUBROUTINE sanity_check(qs_env, training_lists)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(training_list_type), DIMENSION(:), TARGET     :: training_lists

      INTEGER                                            :: ikind, pao_basis_size
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(training_list_type), POINTER                  :: training_list

      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)

      DO ikind = 1, SIZE(training_lists)
         training_list => training_lists(ikind)
         IF (training_list%npoints > 0) CYCLE ! it's ok
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set, pao_basis_size=pao_basis_size)
         IF (pao_basis_size /= basis_set%nsgf) & ! if this kind has pao enabled...
            CPABORT("Found no training-points for kind: "//TRIM(training_list%kindname))
      END DO

   END SUBROUTINE sanity_check

! **************************************************************************************************
!> \brief Turns the linked lists of training points into matrices
!> \param training_lists ...
!> \param training_matrices ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE training_list2matrix(training_lists, training_matrices, para_env)
      TYPE(training_list_type), ALLOCATABLE, &
         DIMENSION(:), TARGET                            :: training_lists
      TYPE(training_matrix_type), ALLOCATABLE, &
         DIMENSION(:), TARGET                            :: training_matrices
      TYPE(cp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: i, ikind, inp_size, ninputs, noutputs, &
                                                            npoints, out_size
      TYPE(training_list_type), POINTER                  :: training_list
      TYPE(training_matrix_type), POINTER                :: training_matrix
      TYPE(training_point_type), POINTER                 :: cur_point, prev_point

      CPASSERT(ALLOCATED(training_lists) .AND. .NOT. ALLOCATED(training_matrices))

      ALLOCATE (training_matrices(SIZE(training_lists)))

      DO ikind = 1, SIZE(training_lists)
         training_list => training_lists(ikind)
         training_matrix => training_matrices(ikind)
         training_matrix%kindname = training_list%kindname ! copy kindname
         npoints = training_list%npoints ! number of points
         IF (npoints == 0) THEN
            ALLOCATE (training_matrix%inputs(0, 0))
            ALLOCATE (training_matrix%outputs(0, 0))
            CYCLE
         ENDIF

         ! figure out size of input and output
         inp_size = 0; out_size = 0
         IF (ALLOCATED(training_list%head%input)) &
            inp_size = SIZE(training_list%head%input)
         IF (ALLOCATED(training_list%head%output)) &
            out_size = SIZE(training_list%head%output)
         CALL mp_sum(inp_size, para_env%group)
         CALL mp_sum(out_size, para_env%group)

         ! allocate matices to hold all training points
         ALLOCATE (training_matrix%inputs(inp_size, npoints))
         ALLOCATE (training_matrix%outputs(out_size, npoints))
         training_matrix%inputs(:, :) = 0.0_dp
         training_matrix%outputs(:, :) = 0.0_dp

         ! loop over all training points, consume linked-list in the process
         ninputs = 0; noutputs = 0
         cur_point => training_list%head
         NULLIFY (training_list%head)
         DO i = 1, npoints
            IF (ALLOCATED(cur_point%input)) THEN
               training_matrix%inputs(:, i) = cur_point%input(:)
               ninputs = ninputs+1
            ENDIF
            IF (ALLOCATED(cur_point%output)) THEN
               training_matrix%outputs(:, i) = cur_point%output(:)
               noutputs = noutputs+1
            ENDIF
            ! advance to next entry and deallocate the current one
            prev_point => cur_point
            cur_point => cur_point%next
            DEALLOCATE (prev_point)
         ENDDO
         training_list%npoints = 0 ! list is now empty

         ! sync training_matrix across ranks
         CALL mp_sum(training_matrix%inputs, para_env%group)
         CALL mp_sum(training_matrix%outputs, para_env%group)

         ! sanity check
         CALL mp_sum(noutputs, para_env%group)
         CALL mp_sum(ninputs, para_env%group)
         CPASSERT(noutputs == npoints .AND. ninputs == npoints)
      ENDDO

   END SUBROUTINE training_list2matrix

! **************************************************************************************************
!> \brief TODO
!> \param ml_prior ...
!> \param training_matrices ...
! **************************************************************************************************
   SUBROUTINE pao_ml_substract_prior(ml_prior, training_matrices)
      INTEGER, INTENT(IN)                                :: ml_prior
      TYPE(training_matrix_type), DIMENSION(:), TARGET   :: training_matrices

      INTEGER                                            :: i, ikind, npoints, out_size
      TYPE(training_matrix_type), POINTER                :: training_matrix

      DO ikind = 1, SIZE(training_matrices)
         training_matrix => training_matrices(ikind)
         out_size = SIZE(training_matrix%outputs, 1)
         npoints = SIZE(training_matrix%outputs, 2)
         IF (npoints == 0) CYCLE
         ALLOCATE (training_matrix%prior(out_size))

         ! calculate prior
         SELECT CASE (ml_prior)
         CASE (pao_ml_prior_zero)
            training_matrix%prior(:) = 0.0_dp
         CASE (pao_ml_prior_mean)
            training_matrix%prior(:) = SUM(training_matrix%outputs, 2)/REAL(npoints, dp)
         CASE DEFAULT
            CPABORT("PAO: unknown prior")
         END SELECT

         ! substract prior from all training points
         DO i = 1, npoints
            training_matrix%outputs(:, i) = training_matrix%outputs(:, i)-training_matrix%prior
         ENDDO
      ENDDO

   END SUBROUTINE pao_ml_substract_prior

! **************************************************************************************************
!> \brief Print some statistics about the training set and dump it upon request
!> \param pao ...
!> \param training_matrices ...
! **************************************************************************************************
   SUBROUTINE pao_ml_print(pao, training_matrices)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(training_matrix_type), DIMENSION(:), TARGET   :: training_matrices

      INTEGER                                            :: i, ikind, N, npoints
      TYPE(training_matrix_type), POINTER                :: training_matrix

      ! dump training data
      IF (pao%iw_mldata > 0) THEN
         DO ikind = 1, SIZE(training_matrices)
            training_matrix => training_matrices(ikind)
            npoints = SIZE(training_matrix%outputs, 2)
            DO i = 1, npoints
               WRITE (pao%iw_mldata, *) "PAO|ML| training-point kind: ", TRIM(training_matrix%kindname), &
                  " point:", i, " in:", training_matrix%inputs(:, i), &
                  " out:", training_matrix%outputs(:, i)
            ENDDO
         ENDDO
         CALL m_flush(pao%iw_mldata)
      ENDIF

      ! print stats
      IF (pao%iw > 0) THEN
         DO ikind = 1, SIZE(training_matrices)
            training_matrix => training_matrices(ikind)
            N = SIZE(training_matrix%inputs)
            IF (N == 0) CYCLE
            WRITE (pao%iw, "(A,I3,A,E10.1,1X,E10.1,1X,E10.1)") " PAO|ML| Descriptor for kind: "// &
               TRIM(training_matrix%kindname)//" size: ", &
               SIZE(training_matrix%inputs, 1), " min/mean/max: ", &
               MINVAL(training_matrix%inputs), &
               SUM(training_matrix%inputs)/REAL(N, dp), &
               MAXVAL(training_matrix%inputs)
         ENDDO
      ENDIF

   END SUBROUTINE pao_ml_print

! **************************************************************************************************
!> \brief Calls the actuall learning algorthim to traing on the given matrices
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_ml_train(pao)
      TYPE(pao_env_type), POINTER                        :: pao

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_ml_train'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      SELECT CASE (pao%ml_method)
      CASE (pao_ml_gp)
         CALL pao_ml_gp_train(pao)
      CASE (pao_ml_nn)
         CALL pao_ml_nn_train(pao)
      CASE (pao_ml_lazy)
         ! nothing to do
      CASE DEFAULT
         CPABORT("PAO: unkown machine learning scheme")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE pao_ml_train

! **************************************************************************************************
!> \brief Fills pao%matrix_X based on machine learning predictions
!> \param pao ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE pao_ml_predict(pao, qs_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_ml_predict'

      INTEGER                                            :: acol, arow, handle, iatom, ikind, natoms
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: descriptor, variances
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_X
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_dbcsr_iterator)                            :: iter
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      para_env=para_env, &
                      cell=cell, &
                      particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      natom=natoms)

      ! fill matrix_X
      ALLOCATE (variances(natoms))
      variances(:) = 0.0_dp
!$OMP PARALLEL DEFAULT(NONE) SHARED(pao,qs_env,particle_set,qs_kind_set,cell,variances) &
!$OMP PRIVATE(iter,arow,acol,iatom,ikind,descriptor,block_X)
      CALL cp_dbcsr_iterator_start(iter, pao%matrix_X)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
         CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_X)
         iatom = arow; CPASSERT(arow == acol)
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         IF (SIZE(block_X) == 0) CYCLE ! pao disabled for iatom

         ! calculate decriptor
         CALL pao_ml_calc_descriptor(pao, &
                                     particle_set, &
                                     qs_kind_set, &
                                     cell, &
                                     iatom, &
                                     descriptor)

         ! call actual machine learning for prediction
         CALL pao_ml_predict_low(pao, ikind=ikind, &
                                 descriptor=descriptor, &
                                 output=block_X(:, 1), &
                                 variance=variances(iatom))

         DEALLOCATE (descriptor)

         !add prior
         block_X(:, 1) = block_X(:, 1)+pao%ml_training_matrices(ikind)%prior
      ENDDO
      CALL cp_dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      ! print variances
      CALL mp_sum(variances, para_env%group)
      IF (pao%iw_mlvar > 0) THEN
         DO iatom = 1, natoms
            WRITE (pao%iw_mlvar, *) "PAO|ML| atom:", iatom, " prediction variance:", variances(iatom)
         ENDDO
         CALL m_flush(pao%iw_mlvar)
      ENDIF

      ! one-line summary
      IF (pao%iw > 0) WRITE (pao%iw, "(A,E20.10,A,T71,I10)") " PAO|ML| max prediction variance:", &
         MAXVAL(variances), " for atom:", MAXLOC(variances)

      IF (MAXVAL(variances) > pao%ml_tolerance) &
         CPABORT("Variance of prediction above ML_TOLERANCE.")

      DEALLOCATE (variances)

      CALL timestop(handle)

   END SUBROUTINE pao_ml_predict

! **************************************************************************************************
!> \brief Queries the actuall learning algorthim to make a prediction
!> \param pao ...
!> \param ikind ...
!> \param descriptor ...
!> \param output ...
!> \param variance ...
! **************************************************************************************************
   SUBROUTINE pao_ml_predict_low(pao, ikind, descriptor, output, variance)
      TYPE(pao_env_type), POINTER                        :: pao
      INTEGER, INTENT(IN)                                :: ikind
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: descriptor
      REAL(dp), DIMENSION(:), INTENT(OUT)                :: output
      REAL(dp), INTENT(OUT)                              :: variance

      SELECT CASE (pao%ml_method)
      CASE (pao_ml_gp)
         CALL pao_ml_gp_predict(pao, ikind, descriptor, output, variance)
      CASE (pao_ml_nn)
         CALL pao_ml_nn_predict(pao, ikind, descriptor, output, variance)
      CASE (pao_ml_lazy)
         output = 0.0_dp ! let's be really lazy and just rely on the prior
         variance = 0
      CASE DEFAULT
         CPABORT("PAO: unkown machine learning scheme")
      END SELECT

   END SUBROUTINE pao_ml_predict_low

! **************************************************************************************************
!> \brief Calculate forces contributed by machine learning
!> \param pao ...
!> \param qs_env ...
!> \param matrix_G ...
!> \param forces ...
! **************************************************************************************************
   SUBROUTINE pao_ml_forces(pao, qs_env, matrix_G, forces)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_dbcsr_type)                                :: matrix_G
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: forces

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_ml_forces'

      INTEGER                                            :: acol, arow, handle, iatom, ikind
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: descr_grad, descriptor
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_G
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_dbcsr_iterator)                            :: iter
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      para_env=para_env, &
                      cell=cell, &
                      particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set)

!$OMP PARALLEL DEFAULT(NONE) SHARED(pao,matrix_G,particle_set,qs_kind_set,cell) &
!$OMP REDUCTION(+:forces) &
!$OMP PRIVATE(iter,arow,acol,iatom,ikind,block_G,descriptor,descr_grad)
      CALL cp_dbcsr_iterator_start(iter, matrix_G)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
         CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_G)
         iatom = arow; CPASSERT(arow == acol)
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         IF (SIZE(block_G) == 0) CYCLE ! pao disabled for iatom

         ! calculate decriptor
         CALL pao_ml_calc_descriptor(pao, &
                                     particle_set, &
                                     qs_kind_set, &
                                     cell, &
                                     iatom=iatom, &
                                     descriptor=descriptor)

         ! calcaulte derivate of machine learning prediction
         CALL pao_ml_gradient_low(pao, ikind=ikind, &
                                  descriptor=descriptor, &
                                  outer_deriv=block_G(:, 1), &
                                  gradient=descr_grad)

         ! calculate force contributions from descriptor
         CALL pao_ml_calc_descriptor(pao, &
                                     particle_set, &
                                     qs_kind_set, &
                                     cell, &
                                     iatom=iatom, &
                                     descr_grad=descr_grad, &
                                     forces=forces)

         DEALLOCATE (descriptor, descr_grad)
      ENDDO
      CALL cp_dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      CALL timestop(handle)

   END SUBROUTINE pao_ml_forces

! **************************************************************************************************
!> \brief Calculate gradient of machine learning algorithm
!> \param pao ...
!> \param ikind ...
!> \param descriptor ...
!> \param outer_deriv ...
!> \param gradient ...
! **************************************************************************************************
   SUBROUTINE pao_ml_gradient_low(pao, ikind, descriptor, outer_deriv, gradient)
      TYPE(pao_env_type), POINTER                        :: pao
      INTEGER, INTENT(IN)                                :: ikind
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: descriptor, outer_deriv
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: gradient

      ALLOCATE (gradient(SIZE(descriptor)))

      SELECT CASE (pao%ml_method)
      CASE (pao_ml_gp)
         CALL pao_ml_gp_gradient(pao, ikind, descriptor, outer_deriv, gradient)
      CASE (pao_ml_nn)
         CALL pao_ml_nn_gradient(pao, ikind, descriptor, outer_deriv, gradient)
      CASE (pao_ml_lazy)
         gradient = 0.0_dp
      CASE DEFAULT
         CPABORT("PAO: unkown machine learning scheme")
      END SELECT

   END SUBROUTINE pao_ml_gradient_low

END MODULE pao_ml
