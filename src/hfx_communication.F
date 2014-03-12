!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for data exchange between MPI processes
!> \par History
!>      04.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
MODULE hfx_communication
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_get_block_p,&
                                             cp_dbcsr_iterator_blocks_left,&
                                             cp_dbcsr_iterator_next_block,&
                                             cp_dbcsr_iterator_start,&
                                             cp_dbcsr_iterator_stop
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_iterator,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE hfx_types,                       ONLY: hfx_2D_map,&
                                             hfx_basis_type,&
                                             hfx_type
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE mathconstants
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_isendrecv,&
                                             mp_max,&
                                             mp_sync,&
                                             mp_waitall
  USE orbital_pointers
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC get_full_density,&
         distribute_ks_matrix,&
         scale_and_add_fock_to_ks_matrix,&
         get_atomic_block_maps
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_communication'

!***

  CONTAINS

! *****************************************************************************
!> \brief - Collects full density matrix from all CPUs
!> \param para_en v
!> \param full_density The full Density matrix
!> \param rho Distributed density
!> \param number_of_p_entries Maximal buffer size
!> \param is_assoc_atomic_block tells us which node owns which atomic pairs
!> \param particle_se t
!> \param nato m
!> \param last_sgf_global helper array for mapping local to global indices
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      - Communication with left/right node only
!>        added a mp_sync before and after the ring of isendrecv. This *speed up* the
!>        communication, and might protect against idle neighbors flooding a busy node
!>        with messages [Joost]
!> \par History
!>      11.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE get_full_density(para_env, full_density, rho, number_of_p_entries, &
                              block_offset, natom,  kind_of, basis_parameter,&
                              get_max_vals_spin, rho_beta, antisymmetric, error)

    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(dp), DIMENSION(:), POINTER          :: full_density
    TYPE(cp_dbcsr_type), POINTER             :: rho
    INTEGER, INTENT(IN)                      :: number_of_p_entries
    INTEGER, DIMENSION(:), POINTER           :: block_offset
    INTEGER, INTENT(IN)                      :: natom
    INTEGER                                  :: kind_of(*)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    LOGICAL, INTENT(IN)                      :: get_max_vals_spin
    TYPE(cp_dbcsr_type), OPTIONAL, POINTER   :: rho_beta
    LOGICAL, INTENT(IN)                      :: antisymmetric
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_full_density', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk, block_size, data_from, dest, i, iatom, icpu, ikind, iset, &
      istat, jatom, jkind, jset, mepos, ncpu, nseta, nsetb, pa, pa1, pb, pb1, &
      req(2), source, source_cpu
    INTEGER, DIMENSION(:), POINTER           :: nsgfa, nsgfb
    LOGICAL                                  :: failure, found
    REAL(dp)                                 :: symmfac
    REAL(dp), DIMENSION(:), POINTER          :: recbuffer, sendbuffer, &
                                                swapbuffer
    REAL(dp), DIMENSION(:, :), POINTER       :: sparse_block, &
                                                sparse_block_beta
    TYPE(cp_dbcsr_iterator)                  :: iter

    failure = .FALSE.

    full_density = 0.0_dp
    ALLOCATE(sendbuffer(number_of_p_entries),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(recbuffer(number_of_p_entries),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    i=1
    CALL cp_dbcsr_iterator_start(iter, rho, shared=.FALSE.)
    DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
      CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, sparse_block, blk)
      ! the resulting vector will eb only the upper triangle.
      ! in case of antisymmetry take care to change signs if a lower block gets copied
      symmfac=1.0_dp
      IF(antisymmetric.AND.iatom>jatom)symmfac=-1.0_dp
      ikind = kind_of(iatom)
      nseta = basis_parameter(ikind)%nset
      nsgfa => basis_parameter(ikind)%nsgf
        jkind = kind_of(jatom)
        nsetb = basis_parameter(jkind)%nset
        nsgfb => basis_parameter(jkind)%nsgf
        IF( get_max_vals_spin ) THEN
          CALL cp_dbcsr_get_block_p(rho_beta,&
               row=iatom,col=jatom,BLOCK=sparse_block_beta,found=found)
          pa = 0
          DO iset = 1,nseta
            pb = 0
            DO jset = 1,nsetb
              DO pa1 = pa+1, pa + nsgfa(iset)
                DO pb1 = pb+1, pb + nsgfb(jset)
                  sendbuffer(i) = MAX(ABS(sparse_block(pa1,pb1)),ABS(sparse_block_beta(pa1,pb1)))
                  i = i + 1
                END DO
              END DO
              pb = pb + nsgfb(jset)
            END DO
            pa = pa + nsgfa(iset)
          END DO
        ELSE
          pa = 0
          DO iset = 1,nseta
            pb = 0
            DO jset = 1,nsetb
              DO pa1 = pa+1, pa + nsgfa(iset)
                DO pb1 = pb+1, pb + nsgfb(jset)
                  sendbuffer(i) = sparse_block(pa1,pb1)*symmfac
                  i = i + 1
                END DO
              END DO
              pb = pb + nsgfb(jset)
            END DO
            pa = pa + nsgfa(iset)
          END DO
        END IF
    END DO
    CALL cp_dbcsr_iterator_stop(iter)

    ! sync before/after a ring of isendrecv
    CALL mp_sync(para_env%group)
    ncpu = para_env%num_pe
    mepos = para_env%mepos
    dest  =MODULO(mepos+1,ncpu)
    source=MODULO(mepos-1,ncpu)
    DO icpu = 0, ncpu-1
      IF (icpu .NE. ncpu-1) THEN
        CALL mp_isendrecv(sendbuffer,dest,recbuffer,source, &
                          para_env%group,req(1),req(2),13)
      ENDIF
      data_from = MODULO(mepos - icpu,ncpu)
      source_cpu = MODULO(data_from,ncpu) + 1
      block_size = block_offset(source_cpu+1) - block_offset(source_cpu)
      full_density(block_offset(source_cpu):block_offset(source_cpu)+block_size-1) = sendbuffer(1:block_size)

      IF (icpu .NE. ncpu -1) THEN
        CALL mp_waitall(req)
      ENDIF
      swapbuffer=>sendbuffer
      sendbuffer=>recbuffer
      recbuffer=>swapbuffer
    ENDDO
    DEALLOCATE(sendbuffer, recbuffer, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ! sync before/after a ring of isendrecv
    CALL mp_sync(para_env%group)

  END SUBROUTINE get_full_density

! *****************************************************************************
!> \brief - Distributes the local full Kohn-Sham matrix to all CPUS
!> \param para_en v
!> \param full_ks The full Kohn-Sham matrix
!> \param ks_matrix Distributed Kohn-Sham matrix
!> \param number_of_p_entries Maximal buffer size
!> \param is_assoc_atomic_block tells us which node owns which atomic pairs
!> \param particle_se t
!> \param nato m
!> \param last_sgf_global helper array for mapping local to global indices
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      - Communication with left/right node only
!> \par History
!>      11.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE distribute_ks_matrix(para_env, full_ks, ks_matrix, number_of_p_entries, &
                                  block_offset, natom, kind_of, basis_parameter, &
                                  off_diag_fac, diag_fac, error)

    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(dp), DIMENSION(:)                   :: full_ks
    TYPE(cp_dbcsr_type), POINTER             :: ks_matrix
    INTEGER, INTENT(IN)                      :: number_of_p_entries
    INTEGER, DIMENSION(:), POINTER           :: block_offset
    INTEGER, INTENT(IN)                      :: natom
    INTEGER                                  :: kind_of(*)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    REAL(dp), INTENT(IN), OPTIONAL           :: off_diag_fac, diag_fac
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'distribute_ks_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk, block_size, data_to, dest, dest_cpu, i, iatom, icpu, &
      ikind, iset, istat, jatom, jkind, jset, mepos, ncpu, nseta, nsetb, pa, &
      pa1, pb, pb1, req(2), source
    INTEGER, DIMENSION(:), POINTER           :: nsgfa, nsgfb
    LOGICAL                                  :: failure
    REAL(dp)                                 :: my_fac, myd_fac
    REAL(dp), DIMENSION(:), POINTER          :: recbuffer, sendbuffer, &
                                                swapbuffer
    REAL(dp), DIMENSION(:, :), POINTER       :: sparse_block
    TYPE(cp_dbcsr_iterator)                  :: iter

    my_fac = 1.0_dp; myd_fac=1.0_dp
    IF( PRESENT( off_diag_fac ) ) my_fac = off_diag_fac
    IF( PRESENT( diag_fac ) ) myd_fac = diag_fac

    ALLOCATE(sendbuffer(number_of_p_entries), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    sendbuffer =0.0_dp
    ALLOCATE(recbuffer(number_of_p_entries), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    recbuffer =0.0_dp

    ncpu = para_env%num_pe
    mepos = para_env%mepos
    dest  =MODULO(mepos+1,ncpu)
    source=MODULO(mepos-1,ncpu)

    ! sync before/after a ring of isendrecv
    CALL mp_sync(para_env%group)
    DO icpu = 1,ncpu
      i = 1
      data_to = mepos-icpu
      dest_cpu = MODULO(data_to, ncpu) + 1
      block_size = block_offset(dest_cpu+1) - block_offset(dest_cpu)
      sendbuffer(1:block_size) = sendbuffer(1:block_size) + full_ks(block_offset(dest_cpu):block_offset(dest_cpu)+block_size-1)
      IF ( icpu .EQ. ncpu ) EXIT
      CALL mp_isendrecv(sendbuffer,dest,recbuffer,source, &
                        para_env%group,req(1),req(2),13)

      CALL mp_waitall(req)
      swapbuffer=>sendbuffer
      sendbuffer=>recbuffer
      recbuffer=>swapbuffer
    ENDDO
    ! sync before/after a ring of isendrecv
    CALL mp_sync(para_env%group)

    i = 1
    CALL cp_dbcsr_iterator_start(iter, ks_matrix, shared=.FALSE.)
    DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
       CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, sparse_block, blk)

      ikind = kind_of(iatom)
      nseta = basis_parameter(ikind)%nset
      nsgfa => basis_parameter(ikind)%nsgf
        jkind = kind_of(jatom)
        nsetb = basis_parameter(jkind)%nset
        nsgfb => basis_parameter(jkind)%nsgf
        pa = 0
        DO iset = 1,nseta
          pb = 0
          DO jset = 1,nsetb
            DO pa1 = pa+1,pa + nsgfa(iset)
              DO pb1 = pb+1,pb + nsgfb(jset)
                IF(iatom == jatom .AND. pa1==pb1) THEN
                  sparse_block(pa1,pb1) = sendbuffer(i)*myd_fac + sparse_block(pa1,pb1)
                ELSE
                  sparse_block(pa1,pb1) = sendbuffer(i)*my_fac + sparse_block(pa1,pb1)
                END IF
                i = i + 1
              END DO
            END DO
            pb = pb + nsgfb(jset)
          END DO
          pa = pa + nsgfa(iset)
        END DO
    END DO
    CALL cp_dbcsr_iterator_stop(iter)

    DEALLOCATE(sendbuffer, recbuffer, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE distribute_ks_matrix

! *****************************************************************************
!> \brief - Distributes the local full Kohn-Sham matrix to all CPUS. Is called in
!>        case of adiabatic rescaling. This is just a refactored version of
!>        distribute_ks_matrix
!> \param para_en v
!> \param full_ks The full Kohn-Sham matrix
!> \param ks_matrix Distributed Kohn-Sham matrix
!> \param number_of_p_entries Maximal buffer size
!> \param is_assoc_atomic_block tells us which node owns which atomic pairs
!> \param particle_se t
!> \param nato m
!> \param last_sgf_global helper array for mapping local to global indices
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      - Communication with left/right node only
!> \par History
!>      11.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE scale_and_add_fock_to_ks_matrix(para_env, qs_env, ks_matrix, irep ,&
                                             scaling_factor, error)

    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    INTEGER, INTENT(IN)                      :: irep
    REAL(dp), INTENT(IN)                     :: scaling_factor
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'scale_and_add_fock_to_ks_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, ikind, istat, natom, &
                                                nspins, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: kind_of, last_sgf_global
    LOGICAL                                  :: failure
    REAL(dp), DIMENSION(:), POINTER          :: full_ks
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(hfx_type), POINTER                  :: actual_x_data
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

!! All shared data is saved in i_thread = 1!

    actual_x_data => qs_env%x_data(irep, 1)
    basis_parameter => actual_x_data%basis_parameter

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    error=error)

    nspins = qs_env%dft_control%nspins
    natom = SIZE(particle_set,1)
    ALLOCATE(kind_of(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             kind_of=kind_of)

    ALLOCATE(last_sgf_global(0:natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    last_sgf_global(0)=0
    DO iatom=1,natom
      ikind = kind_of(iatom)
      last_sgf_global(iatom) = last_sgf_global(iatom-1)+ basis_parameter(ikind)%nsgf_total
    END DO
    full_ks => actual_x_data%full_ks_alpha
    IF(scaling_factor /= 1.0_dp ) THEN
       full_ks = full_ks * scaling_factor
    END IF
    CALL distribute_ks_matrix(para_env, full_ks, ks_matrix(1)%matrix, actual_x_data%number_of_p_entries, &
                              actual_x_data%block_offset, natom, kind_of, basis_parameter, &
                              off_diag_fac=0.5_dp,error=error)
    DEALLOCATE(actual_x_data%full_ks_alpha, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF(nspins == 2) THEN
      full_ks => actual_x_data%full_ks_beta
      IF(scaling_factor /= 1.0_dp ) THEN
         full_ks = full_ks * scaling_factor
      END IF
      CALL distribute_ks_matrix(para_env, full_ks, ks_matrix(2)%matrix, actual_x_data%number_of_p_entries, &
                              actual_x_data%block_offset, natom, kind_of, basis_parameter, &
                              off_diag_fac=0.5_dp,error=error)

      DEALLOCATE(actual_x_data%full_ks_beta, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    DEALLOCATE(kind_of, last_sgf_global, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE scale_and_add_fock_to_ks_matrix

! *****************************************************************************
!> \brief Given a 2d index pair, this function returns a 1d index pair for
!>        a symmetric upper triangle NxN matrix
!>        The compiler should inline this function, therefore it appears in
!>        several modules
!> \param i,j 2d index
!> \param N matrix size
!> \par History
!>      03.2009 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  PURE FUNCTION get_1D_idx(i,j,N)
    INTEGER, INTENT(IN)                      :: i, j
    INTEGER(int_8), INTENT(IN)               :: N
    INTEGER(int_8)                           :: get_1D_idx

    INTEGER(int_8)                           :: min_ij

    min_ij = MIN(i,j)
    get_1D_idx = min_ij*N + MAX(i,j) - (min_ij-1)*min_ij/2 - N

  END FUNCTION get_1D_idx

! *****************************************************************************
!> \brief create a several maps array that reflects the ks matrix sparsity
!> \param para_en v
!> \par History
!>      11.2007 refactored [Joost VandeVondele]
!>      07.2009 add new maps
!> \notes
!>      is_assoc_atomic_block returns the mpi rank + 1 for associated blocks,
!>      zero for unassiated blocks
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE get_atomic_block_maps(matrix,basis_parameter,kind_of,&
                                   is_assoc_atomic_block,number_of_p_entries,&
                                   para_env,atomic_block_offset, set_offset, &
                                   block_offset, map_atoms_to_cpus, nkind, error)

    TYPE(cp_dbcsr_type), POINTER             :: matrix
    TYPE(hfx_basis_type), DIMENSION(:)       :: basis_parameter
    INTEGER, DIMENSION(:)                    :: kind_of
    INTEGER, DIMENSION(:, :), INTENT(OUT)    :: is_assoc_atomic_block
    INTEGER, INTENT(OUT)                     :: number_of_p_entries
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, DIMENSION(:, :), POINTER        :: atomic_block_offset
    INTEGER, DIMENSION(:, :, :, :), POINTER  :: set_offset
    INTEGER, DIMENSION(:), POINTER           :: block_offset
    TYPE(hfx_2D_map), DIMENSION(:), POINTER  :: map_atoms_to_cpus
    INTEGER                                  :: nkind
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_atomic_block_maps', &
      routineP = moduleN//':'//routineN

    INTEGER :: blk, handle, iatom, ibuf, icpu, ikind, ilist, iset, itask, &
      jatom, jkind, jset, natom, ncpu, nseta, nsetb, number_of_p_blocks, &
      offset, stat, tmp(2)
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: buffer_in, buffer_out, &
                                                counter, rcount, rdispl
    INTEGER, DIMENSION(:), POINTER           :: iatom_list, jatom_list, &
                                                nsgfa, nsgfb
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: sparse_block
    TYPE(cp_dbcsr_iterator)                  :: iter

    failure=.FALSE.

    CALL timeset(routineN,handle)

    is_assoc_atomic_block = 0
    number_of_p_entries = 0
    number_of_p_blocks = 0

    !
    ! count number_of_p_blocks and number_of_p_entries
    !
    CALL cp_dbcsr_iterator_start(iter, matrix, shared=.FALSE.)
    DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
       CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, sparse_block, blk)
       ikind = kind_of(iatom)
        jkind = kind_of(jatom)
        number_of_p_blocks = number_of_p_blocks + 1
        number_of_p_entries =  number_of_p_entries + &
              basis_parameter(ikind)%nsgf_total*basis_parameter(jkind)%nsgf_total
    END DO
    CALL cp_dbcsr_iterator_stop(iter)

    tmp=(/number_of_p_entries,number_of_p_blocks/)
    CALL mp_max(tmp, para_env%group)
    number_of_p_entries=tmp(1)
    number_of_p_blocks=tmp(2)
    !
    ! send this info around, so we can construct is_assoc_atomic_block
    ! pack all data in buffers and use allgatherv
    !
    ALLOCATE(buffer_in(3*number_of_p_blocks),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(buffer_out(3*number_of_p_blocks*para_env%num_pe),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(rcount(para_env%num_pe),rdispl(para_env%num_pe),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    buffer_in=0
    ibuf=0

    CALL cp_dbcsr_iterator_start(iter, matrix, shared=.FALSE.)
    DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
       CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, sparse_block, blk)

        buffer_in(ibuf+1)=iatom
        buffer_in(ibuf+2)=jatom
        buffer_in(ibuf+3)=para_env%mepos + 1
        ibuf=ibuf+3
    END DO
    CALL cp_dbcsr_iterator_stop(iter)

    rcount=SIZE(buffer_in)
    rdispl(1)=0
    DO icpu=2,para_env%num_pe
       rdispl(icpu)=rdispl(icpu-1)+rcount(icpu-1)
    ENDDO
    CALL mp_allgather(buffer_in,buffer_out,rcount,rdispl,para_env%group)

    DO ibuf=0,number_of_p_blocks*para_env%num_pe*3-3,3
       itask=buffer_out(ibuf+3)
       ! buffer_out can be 0 if buffer_in contained less elements than the max number of atom pairs
       ! is_assoc_atomic_block is a map for atom pairs to a processor (assumes symmetry, i,j on the ame as j,i)
       IF (itask.NE.0) THEN
           iatom=buffer_out(ibuf+1)
           jatom=buffer_out(ibuf+2)
           is_assoc_atomic_block(iatom,jatom)=itask
           is_assoc_atomic_block(jatom,iatom)=itask
       ENDIF
    ENDDO

    IF( ASSOCIATED(map_atoms_to_cpus) ) THEN
      DO icpu=1,para_env%num_pe
        DEALLOCATE(map_atoms_to_cpus(icpu)%iatom_list)
        DEALLOCATE(map_atoms_to_cpus(icpu)%jatom_list)
      END DO
      DEALLOCATE(map_atoms_to_cpus)
    END IF

    natom = SIZE(is_assoc_atomic_block,1)
    ALLOCATE(map_atoms_to_cpus(para_env%num_pe), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(counter(para_env%num_pe), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    counter = 0

    DO iatom=1,natom
      DO jatom=iatom,natom
        icpu = is_assoc_atomic_block(jatom,iatom)
        IF( icpu>0 ) counter(icpu) = counter(icpu) + 1
      END DO
    END DO
    DO icpu = 1,para_env%num_pe
      ALLOCATE(map_atoms_to_cpus(icpu)%iatom_list(counter(icpu)),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(map_atoms_to_cpus(icpu)%jatom_list(counter(icpu)),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END DO
    counter = 0
    DO iatom=1,natom
      DO jatom=iatom,natom
        icpu = is_assoc_atomic_block(jatom,iatom)
        IF( icpu >0 ) THEN
          counter(icpu) = counter(icpu) + 1
          map_atoms_to_cpus(icpu)%jatom_list(counter(icpu))=jatom
          map_atoms_to_cpus(icpu)%iatom_list(counter(icpu))=iatom
        END IF
      END DO
    END DO

    DEALLOCATE(counter, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ncpu = para_env%num_pe
    offset = 1
    atomic_block_offset = 0
    block_offset = 0
    DO icpu = 1, ncpu
      iatom_list => map_atoms_to_cpus(icpu)%iatom_list
      jatom_list => map_atoms_to_cpus(icpu)%jatom_list
      block_offset(icpu) = offset
      DO ilist = 1,SIZE(iatom_list)
        iatom = iatom_list(ilist)
        ikind = kind_of(iatom)
        jatom = jatom_list(ilist)
        jkind = kind_of(jatom)
        atomic_block_offset(iatom,jatom) = offset
        atomic_block_offset(jatom,iatom) = offset
        offset = offset + basis_parameter(ikind)%nsgf_total * basis_parameter(jkind)%nsgf_total
      END DO
    END DO
    block_offset(ncpu+1) = offset
    set_offset = 0
    DO ikind=1,nkind
      nseta = basis_parameter(ikind)%nset
      nsgfa => basis_parameter(ikind)%nsgf
      DO jkind=1,nkind
        nsetb = basis_parameter(jkind)%nset
        nsgfb => basis_parameter(jkind)%nsgf
        offset = 1
        DO iset = 1,nseta
          DO jset = 1,nsetb
            set_offset(jset,iset,jkind,ikind) = offset
            offset = offset + nsgfa(iset)*nsgfb(jset)
          END DO
        END DO
      END DO
    END DO

    CALL timestop(handle)

  END SUBROUTINE get_atomic_block_maps

END MODULE hfx_communication
