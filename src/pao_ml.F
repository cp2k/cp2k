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
   USE cell_types,                      ONLY: cell_create,&
                                              cell_type
   USE cp_dbcsr_interface,              ONLY: cp_dbcsr_get_block_p,&
                                              cp_dbcsr_iterator,&
                                              cp_dbcsr_iterator_blocks_left,&
                                              cp_dbcsr_iterator_next_block,&
                                              cp_dbcsr_iterator_start,&
                                              cp_dbcsr_iterator_stop,&
                                              cp_dbcsr_p_type
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_bcast,&
                                              mp_sum
   USE pao_input,                       ONLY: id2str,&
                                              pao_exp_param,&
                                              pao_ml_gp,&
                                              pao_ml_nn
   USE pao_io,                          ONLY: pao_ioblock_type,&
                                              pao_iokind_type,&
                                              pao_kinds_ensure_equal,&
                                              pao_read_raw
   USE pao_ml_codec,                    ONLY: pao_ml_decode_xblock,&
                                              pao_ml_encode_xblock
   USE pao_ml_descriptor,               ONLY: pao_ml_calc_descriptor
   USE pao_ml_gaussprocess,             ONLY: pao_ml_gp_predict,&
                                              pao_ml_gp_train
   USE pao_ml_neuralnet,                ONLY: pao_ml_nn_predict,&
                                              pao_ml_nn_train
   USE pao_types,                       ONLY: pao_env_type,&
                                              training_matrix_type
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_ml'

   PUBLIC :: pao_machine_learning_init

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
   SUBROUTINE pao_machine_learning_init(pao, qs_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_machine_learning_init', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: line
      INTEGER                                            :: ikind, istat, unit_nr
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(training_list_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: training_lists
      TYPE(training_matrix_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: training_matrices

      CPASSERT(LEN_TRIM(pao%ml_training_file) > 0)

      IF (pao%iw > 0) WRITE (pao%iw, *) 'PAO|ML| Initializing maschine learning...'

      IF (pao%parameterization /= pao_exp_param) &
         CPABORT("PAO maschine learning requires EXP parametrization")

      CALL get_qs_env(qs_env, para_env=para_env, atomic_kind_set=atomic_kind_set)

      ! create training-set data-structure
      ALLOCATE (training_lists(SIZE(atomic_kind_set)))
      DO ikind = 1, SIZE(training_lists)
         CALL get_atomic_kind(atomic_kind_set(ikind), name=training_lists(ikind)%kindname)
      ENDDO

      IF (para_env%mepos == para_env%source) THEN
         CALL open_file(file_name=pao%ml_training_file, file_status="OLD", file_form="FORMATTED", &
                        file_action="READ", unit_number=unit_nr)
         DO WHILE (.TRUE.)
            READ (unit_nr, IOSTAT=istat, FMT='(A)') line
            line = ADJUSTL(line)
            IF (istat /= 0) EXIT ! end of file reached
            IF (LEN_TRIM(line) == 0) CYCLE ! skip empty line
            IF (line(1:1) == "!" .OR. line(1:1) == "#") CYCLE ! skipt comment

            CALL add_to_training_list(pao, qs_env, training_lists, filename=line)
         ENDDO

         CALL close_file(unit_number=unit_nr)

         CALL training_list2matrix(training_lists, training_matrices)
         DEALLOCATE (training_lists)
      ENDIF

      CALL bcast_training(para_env, training_matrices)
      CALL pao_ml_train(pao, training_matrices)
      CALL pao_ml_predict(pao, qs_env, training_matrices)

      DEALLOCATE (training_matrices)

   END SUBROUTINE pao_machine_learning_init

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
      INTEGER                                            :: iatom, ikind, natoms
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom2kind, kindsmap
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: hmat, positions
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(pao_ioblock_type), ALLOCATABLE, DIMENSION(:)  :: u0blocks, xblocks
      TYPE(pao_iokind_type), ALLOCATABLE, DIMENSION(:)   :: kinds
      TYPE(particle_type), DIMENSION(:), POINTER         :: my_particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(training_point_type), POINTER                 :: new_point

      NULLIFY (new_point, cell)

      IF (pao%iw > 0) WRITE (pao%iw, '(A,A)') " PAO|ML| Reading training frame from file: ", TRIM(filename)

      ! parse training data
      CALL pao_read_raw(filename, param, hmat, kinds, atom2kind, positions, xblocks, u0blocks)

      ! check parametrization
      IF (TRIM(param) .NE. TRIM(ADJUSTL(id2str(pao%parameterization)))) &
         CPABORT("Restart PAO parametrization does not match")

      ! create cell from read-in h-matrix
      CALL cell_create(cell, hmat)

      ! map read-in kinds onto kinds of this run
      CALL match_kinds(pao, qs_env, kinds, kindsmap)

      ! create a particle_set based on read-in positions and refere to kinds of this run
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set)
      natoms = SIZE(positions, 1)
      ALLOCATE (my_particle_set(natoms))
      DO iatom = 1, natoms
         ikind = kindsmap(atom2kind(iatom))
         my_particle_set(iatom)%atomic_kind => atomic_kind_set(ikind)
         my_particle_set(iatom)%r = positions(iatom, :)
      ENDDO

      IF (.NOT. ALLOCATED(u0blocks)) &
         CPABORT("PAO-ML requires U0_blocks e.g. from EXP parametization")

      ! build training points
      DO ikind = 1, SIZE(kindsmap)
         DO iatom = 1, SIZE(xblocks)
            IF (kindsmap(atom2kind(iatom)) == ikind) THEN
               ALLOCATE (new_point)

               CALL pao_ml_calc_descriptor(my_particle_set, &
                                           qs_kind_set, &
                                           cell, &
                                           iatom=iatom, &
                                           descriptor=new_point%input)

               CALL pao_ml_encode_xblock(my_particle_set, &
                                         qs_kind_set, &
                                         cell, &
                                         iatom=iatom, &
                                         u0block=u0blocks(iatom)%p, &
                                         xblock=xblocks(iatom)%p, &
                                         output=new_point%output)

               ! add to linked list
               training_lists(ikind)%npoints = training_lists(ikind)%npoints+1
               new_point%next => training_lists(ikind)%head
               training_lists(ikind)%head => new_point
            ENDIF
         ENDDO
      ENDDO

      DEALLOCATE (cell, my_particle_set)
      ! ALLOCATABLEs deallocate themselfs

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

   END SUBROUTINE match_kinds

! **************************************************************************************************
!> \brief Turns the linked lists of training points into matrices
!> \param training_lists ...
!> \param training_matrices ...
! **************************************************************************************************
   SUBROUTINE training_list2matrix(training_lists, training_matrices)
      TYPE(training_list_type), ALLOCATABLE, &
         DIMENSION(:), TARGET                            :: training_lists
      TYPE(training_matrix_type), ALLOCATABLE, &
         DIMENSION(:), TARGET                            :: training_matrices

      INTEGER                                            :: i, ikind, inp_size, npoints, out_size
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
         CPASSERT(npoints > 0)

         inp_size = SIZE(training_list%head%input)
         out_size = SIZE(training_list%head%output)

         ! copy from link-list into arrays
         ALLOCATE (training_matrix%inputs(npoints, inp_size))
         ALLOCATE (training_matrix%outputs(npoints, out_size))

         cur_point => training_list%head
         NULLIFY (training_list%head)
         DO i = 1, npoints
            training_matrix%inputs(i, :) = cur_point%input(:)
            training_matrix%outputs(i, :) = cur_point%output(:)
            prev_point => cur_point
            cur_point => cur_point%next
            DEALLOCATE (prev_point) ! deallocate linked-list
         ENDDO
         training_list%npoints = 0 ! list is now empty
      ENDDO

   END SUBROUTINE training_list2matrix

! **************************************************************************************************
!> \brief Broadcasts the matrices of trainging points across all MPI ranks.
!> \param para_env ...
!> \param training_matrices ...
! **************************************************************************************************
   SUBROUTINE bcast_training(para_env, training_matrices)
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(training_matrix_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: training_matrices

      INTEGER                                            :: ikind, nkinds, s1, s2

      IF (para_env%mepos == para_env%source) &
         nkinds = SIZE(training_matrices)
      CALL mp_bcast(nkinds, para_env%source, para_env%group)
      IF (para_env%mepos /= para_env%source) &
         ALLOCATE (training_matrices(nkinds))

      DO ikind = 1, nkinds
         CALL mp_bcast(training_matrices(ikind)%kindname, para_env%source, para_env%group)

         IF (para_env%mepos == para_env%source) THEN
            s1 = SIZE(training_matrices(ikind)%outputs, 1)
            s2 = SIZE(training_matrices(ikind)%outputs, 2)
         ENDIF
         CALL mp_bcast(s1, para_env%source, para_env%group)
         CALL mp_bcast(s2, para_env%source, para_env%group)
         IF (para_env%mepos /= para_env%source) &
            ALLOCATE (training_matrices(ikind)%outputs(s1, s2))
         CALL mp_bcast(training_matrices(ikind)%outputs, para_env%source, para_env%group)

         IF (para_env%mepos == para_env%source) THEN
            s1 = SIZE(training_matrices(ikind)%inputs, 1)
            s2 = SIZE(training_matrices(ikind)%inputs, 2)
         ENDIF
         CALL mp_bcast(s1, para_env%source, para_env%group)
         CALL mp_bcast(s2, para_env%source, para_env%group)
         IF (para_env%mepos /= para_env%source) &
            ALLOCATE (training_matrices(ikind)%inputs(s1, s2))
         CALL mp_bcast(training_matrices(ikind)%inputs, para_env%source, para_env%group)
      ENDDO

   END SUBROUTINE bcast_training

! **************************************************************************************************
!> \brief Calls the actuall learning algorthim to traing on the given matrices
!> \param pao ...
!> \param training_matrices ...
! **************************************************************************************************
   SUBROUTINE pao_ml_train(pao, training_matrices)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(training_matrix_type), DIMENSION(:)           :: training_matrices

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_ml_train', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      SELECT CASE (pao%ml_method)
      CASE (pao_ml_gp)
         CALL pao_ml_gp_train(pao, training_matrices)
      CASE (pao_ml_nn)
         CALL pao_ml_nn_train(pao, training_matrices)
      CASE DEFAULT
         CPABORT("PAO: unkown machine learning scheme")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE pao_ml_train

! **************************************************************************************************
!> \brief Fills pao%matrix_X based on machine learning predictions
!> \param pao ...
!> \param qs_env ...
!> \param training_matrices ...
! **************************************************************************************************
   SUBROUTINE pao_ml_predict(pao, qs_env, training_matrices)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(training_matrix_type), DIMENSION(:)           :: training_matrices

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_ml_predict'

      INTEGER                                            :: acol, arow, handle, iatom, ikind, natoms
      LOGICAL                                            :: found
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: descriptor, output, variances
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_U0, block_X
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_dbcsr_iterator)                            :: iter
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      para_env=para_env, &
                      cell=cell, &
                      matrix_s=matrix_s, &
                      particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      natom=natoms)

      ! fill matrix_X
      ALLOCATE (variances(natoms))
      variances = 0.0_dp
      CALL cp_dbcsr_iterator_start(iter, pao%matrix_X)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
         CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_X)
         iatom = arow; CPASSERT(arow == acol)
         CALL cp_dbcsr_get_block_p(matrix=pao%matrix_U0, row=iatom, col=iatom, block=block_U0, found=found)
         CPASSERT(ASSOCIATED(block_U0))
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)

         ! calculate decriptor
         CALL pao_ml_calc_descriptor(particle_set, &
                                     qs_kind_set, &
                                     cell, &
                                     iatom, &
                                     descriptor)

         ! call actual machine learning for prediction
         CALL pao_ml_predict_low(pao, training_matrices(ikind), descriptor, output, variances(iatom))

         ! decode output into block_X
         CALL pao_ml_decode_xblock(particle_set, &
                                   qs_kind_set, &
                                   cell, &
                                   iatom=iatom, &
                                   u0block=block_U0, &
                                   output=output, &
                                   xblock=block_X)

         DEALLOCATE (output, descriptor)
      ENDDO
      CALL cp_dbcsr_iterator_stop(iter)

      ! print variances
      CALL mp_sum(variances, para_env%group)
      IF (pao%iw > 0) THEN
         DO iatom = 1, natoms
            WRITE (pao%iw, *) "PAO|ML| predicting atom:", iatom, "variance:", variances(iatom)
         ENDDO
         WRITE (pao%iw, *) "PAO|ML| max prediction variance:", MAXVAL(variances)
      ENDIF

      IF (MAXVAL(variances) > pao%ml_tolerance) &
         CPABORT("Variance of prediction above ML_TOLERANCE.")

      DEALLOCATE (variances)

      CALL timestop(handle)

   END SUBROUTINE pao_ml_predict

! **************************************************************************************************
!> \brief Queries the actuall learning algorthim to make a prediction
!> \param pao ...
!> \param training_matrix ...
!> \param descriptor ...
!> \param output ...
!> \param variance ...
! **************************************************************************************************
   SUBROUTINE pao_ml_predict_low(pao, training_matrix, descriptor, output, variance)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(training_matrix_type), INTENT(IN)             :: training_matrix
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: descriptor
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: output
      REAL(dp), INTENT(OUT)                              :: variance

      INTEGER                                            :: npoints, width_out

      CPASSERT(.NOT. ALLOCATED(output))

      npoints = SIZE(training_matrix%outputs, 1)
      CPASSERT(npoints > 0)

      width_out = SIZE(training_matrix%outputs, 2)
      ALLOCATE (output(width_out))

      SELECT CASE (pao%ml_method)
      CASE (pao_ml_gp)
         CALL pao_ml_gp_predict(pao, training_matrix, descriptor, output, variance)
      CASE (pao_ml_nn)
         CALL pao_ml_nn_predict(training_matrix, descriptor, output, variance)
      CASE DEFAULT
         CPABORT("PAO: unkown machine learning scheme")
      END SELECT

   END SUBROUTINE pao_ml_predict_low

END MODULE pao_ml
