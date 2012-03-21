!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Rountines to calculate MP2 energy
!> \par History
!>      05.2011 created [Mauro Del Ben]
!> \author Mauro Del Ben
! *****************************************************************************
MODULE mp2
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type
  USE cp_fm_types,                     ONLY: cp_fm_get_submatrix,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_env,                     ONLY: cp_para_env_create,&
                                             cp_para_env_release
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE hfx_types,                       ONLY: hfx_basis_type,&
                                             hfx_release,&
                                             hfx_type
  USE input_constants,                 ONLY: do_mp2_potential_TShPSC
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE machine,                         ONLY: m_flush,&
                                             m_memory,&
                                             m_walltime
  USE message_passing,                 ONLY: mp_comm_split_direct,&
                                             mp_max,&
                                             mp_sum,&
                                             mp_sync
  USE mp2_direct_method,               ONLY: mp2_canonical_direct_single_batch
  USE mp2_gpw,                         ONLY: mp2_gpw_main
  USE mp2_types,                       ONLY: mp2_biel_type,&
                                             mp2_method_direct,&
                                             mp2_method_gpw,&
                                             mp2_method_laplace,&
                                             mp2_type
  USE particle_types,                  ONLY: particle_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2'

  PUBLIC :: mp2_main

  CONTAINS

! *****************************************************************************
!> \brief the main entry point for MP2 calculations
!> \author Mauro Del Ben
! *****************************************************************************
  SUBROUTINE mp2_main(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mp2_main', &
      routineP = moduleN//':'//routineN

    INTEGER :: dimen, handle, i, iatom, ikind, max_nset, nao, natom, &
      nelec_alpha, nelec_beta, nelectron, nkind, nmo, nspins, stat, unit_nr
    INTEGER(KIND=int_8)                      :: mem
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: kind_of
    LOGICAL                                  :: failure
    REAL(KIND=dp) :: Emp2, Emp2_AA, Emp2_AA_Cou, Emp2_AA_ex, Emp2_AB, &
      Emp2_AB_Cou, Emp2_AB_ex, Emp2_BB, Emp2_BB_Cou, Emp2_BB_ex, Emp2_Cou, &
      Emp2_ex, mem_real, t1, t2
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Auto, Auto_alpha, Auto_beta
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: C, C_alpha, C_beta
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(hfx_type), POINTER                  :: actual_x_data
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(mp2_biel_type)                      :: mp2_biel
    TYPE(mp2_type), POINTER                  :: mp2_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(section_vals_type), POINTER         :: hfx_sections, input

    CALL timeset(routineN,handle)
    failure=.FALSE.
    logger => cp_error_get_logger(error)

    CALL get_qs_env(qs_env=qs_env, input=input, atomic_kind_set=atomic_kind_set, &
                  particle_set=particle_set, para_env=para_env,   &
                  energy=energy, rho=rho, mos=mos, scf_env=scf_env,&
                  matrix_ks=matrix_ks, matrix_s=matrix_s, mp2_env=mp2_env, error=error)



    unit_nr = cp_print_key_unit_nr(logger,input,"DFT%XC%MP2%MP2_INFO",&
                                   extension=".mp2Log",error=error)

    IF (unit_nr>0) THEN
       WRITE(unit_nr,*)
       WRITE(unit_nr,*)
       WRITE(unit_nr,'(T2,A)') 'MP2 section'
       WRITE(unit_nr,'(T2,A)') '-----------'
       WRITE(unit_nr,*)
    ENDIF


    IF(mp2_env%mp2_num_proc<=0 .OR. mp2_env%mp2_num_proc>para_env%num_pe .OR. MOD(para_env%num_pe,mp2_env%mp2_num_proc).NE.0) THEN
       IF (unit_nr>0 .AND. mp2_env%mp2_num_proc.NE.-1) &
           WRITE(unit_nr,'(T3,A,T76,I5)') 'Requested number of processes per group:', mp2_env%mp2_num_proc
       mp2_env%mp2_num_proc=para_env%num_pe
    ENDIF
    IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T76,I5)')     'Used number of processes per group:', mp2_env%mp2_num_proc
    IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T69,F9.2,A3)') 'Maximum allowed memory usage per MPI processes during MP2:',&
                                                    mp2_env%mp2_memory, ' MB'

    IF(mp2_env%method.NE.mp2_method_gpw) THEN
      mem=m_memory()
      mem_real=(mem+1024*1024-1)/(1024*1024)
      CALL mp_max(mem_real,para_env%group)
      mp2_env%mp2_memory=mp2_env%mp2_memory-mem_real
      IF(mp2_env%mp2_memory<0.0_dp) mp2_env%mp2_memory=1.0_dp

      IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T69,F9.2,A3)') 'Available memory per MPI processes for MP2:',&
                                                      mp2_env%mp2_memory, ' MB'
    END IF

    IF (unit_nr>0) CALL m_flush(unit_nr)

    nspins=qs_env%dft_control%nspins

    natom = SIZE(particle_set,1)

    ALLOCATE(kind_of(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               kind_of=kind_of)
    nkind = SIZE(atomic_kind_set,1)

    actual_x_data => qs_env%x_data(1,1)
    basis_parameter   => actual_x_data%basis_parameter

    dimen=0
    max_nset=0
    DO iatom=1, natom
      ikind=kind_of(iatom)
      dimen=dimen + SUM(basis_parameter(ikind)%nsgf)
      max_nset=MAX(max_nset,basis_parameter(ikind)%nset)
    END DO

    CALL get_mo_set( mo_set=mos(1)%mo_set,nao=nao)
    CPPostcondition(dimen==nao,cp_failure_level,routineP,error,failure)

    hfx_sections => section_vals_get_subs_vals(input,"DFT%XC%HF",error=error)

    !   build the table of index
    t1=m_walltime()
    ALLOCATE(mp2_biel%index_table(natom,max_nset),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL build_index_table(natom,max_nset,mp2_biel%index_table,basis_parameter,kind_of)

    Emp2=0.D+00
    Emp2_Cou=0.D+00
    Emp2_ex=0.D+00
    t1=m_walltime()
    SELECT CASE(mp2_env%method)
     CASE(mp2_method_laplace)
       CALL cp_unimplemented_error(fromWhere=routineP, &
              message="laplace not implemented",&
              error=error, error_level=cp_failure_level)
     CASE (mp2_method_direct)
       DO i=1,SIZE(mos)
          CALL get_mo_set( mo_set=mos(i)%mo_set,&
                      nmo=nmo, nao=nao, mo_coeff=mo_coeff)
          IF (nmo.NE.nao) THEN
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                 " Direct MP2 needs the full set of virtual MOs, use ADDED_MOS in the input"//&
                 CPSourceFileRef,&
                 only_ionode=.TRUE.)
          ENDIF
       ENDDO

       IF (unit_nr>0) WRITE(unit_nr,*)

       IF(nspins==2) THEN
         IF (unit_nr>0) WRITE(unit_nr,'(T3,A)') 'Unrestricted Canonical Direct Methods:'
    ! for now, require the mos to be always present

         ! get the alpha coeff and eigenvalues
         CALL get_mo_set( mo_set=mos(1)%mo_set,&
                      nelectron=nelec_alpha,&
                      eigenvalues=mo_eigenvalues,&
                      mo_coeff=mo_coeff)
         ALLOCATE(C_alpha(dimen,dimen),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

         ALLOCATE(Auto_alpha(dimen),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

         CALL cp_fm_get_submatrix(mo_coeff,C_alpha, 1, 1, dimen, dimen, .FALSE., error)
         Auto_alpha(:)=mo_eigenvalues(:)

         ! get the beta coeff and eigenvalues
         CALL get_mo_set( mo_set=mos(2)%mo_set,&
                      nelectron=nelec_beta,&
                      eigenvalues=mo_eigenvalues,&
                      mo_coeff=mo_coeff)
         ALLOCATE(C_beta(dimen,dimen),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

         ALLOCATE(Auto_beta(dimen),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

         CALL cp_fm_get_submatrix(mo_coeff,C_beta, 1, 1, dimen, dimen, .FALSE., error)
         Auto_beta(:)=mo_eigenvalues(:)

         ! calculate the alpha-alpha MP2
         Emp2_AA=0.0_dp
         Emp2_AA_Cou=0.0_dp
         Emp2_AA_ex=0.0_dp
         CALL mp2_direct_energy(dimen,nelec_alpha,nelec_alpha,mp2_biel,mp2_env,C_alpha,Auto_alpha,Emp2_AA,Emp2_AA_Cou,Emp2_AA_ex,&
                                kind_of,basis_parameter,&
                                qs_env,matrix_ks,rho,hfx_sections,para_env, &
                                unit_nr,error=error)
         IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T56,F25.14)') 'MP2 Energy Alpha-Alpha = ', Emp2_AA
         IF (unit_nr>0) WRITE(unit_nr,*)

         Emp2_BB=0.0_dp
         Emp2_BB_Cou=0.0_dp
         Emp2_BB_ex=0.0_dp
         CALL mp2_direct_energy(dimen,nelec_beta,nelec_beta,mp2_biel,mp2_env,C_beta,Auto_beta,Emp2_BB,Emp2_BB_Cou,Emp2_BB_ex,&
                                kind_of,basis_parameter,&
                                qs_env,matrix_ks,rho,hfx_sections,para_env, &
                                unit_nr,error=error)
         IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T56,F25.14)') 'MP2 Energy Beta-Beta= ', Emp2_BB
         IF (unit_nr>0) WRITE(unit_nr,*)

         Emp2_AB=0.0_dp
         Emp2_AB_Cou=0.0_dp
         Emp2_AB_ex=0.0_dp
         CALL mp2_direct_energy(dimen,nelec_alpha,nelec_beta,mp2_biel,mp2_env,C_alpha,&
                                Auto_alpha,Emp2_AB,Emp2_AB_Cou,Emp2_AB_ex,&
                                kind_of,basis_parameter,&
                                qs_env,matrix_ks,rho,hfx_sections,para_env, &
                                unit_nr,C_beta,Auto_beta,error=error)
         IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T56,F25.14)') 'MP2 Energy Alpha-Beta= ', Emp2_AB
         IF (unit_nr>0) WRITE(unit_nr,*)

         Emp2=Emp2_AA+Emp2_BB+Emp2_AB*2.0_dp !+Emp2_BA
         Emp2_Cou=Emp2_AA_Cou+Emp2_BB_Cou+Emp2_AB_Cou*2.0_dp !+Emp2_BA
         Emp2_ex=Emp2_AA_ex+Emp2_BB_ex+Emp2_AB_ex*2.0_dp !+Emp2_BA

       ELSE

         IF (unit_nr>0) WRITE(unit_nr,'(T3,A)') 'Canonical Direct Methods:'
         CALL get_mo_set( mo_set=mos(1)%mo_set,&
                      nelectron=nelectron,&
                      eigenvalues=mo_eigenvalues,&
                      mo_coeff=mo_coeff)
         ALLOCATE(C(dimen,dimen),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

         ALLOCATE(Auto(dimen),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)


         CALL cp_fm_get_submatrix(mo_coeff,C, 1, 1, dimen, dimen, .FALSE., error)
         Auto(:)=mo_eigenvalues(:)

         CALL mp2_direct_energy(dimen,nelectron/2,nelectron/2,mp2_biel,mp2_env,C,Auto,Emp2,Emp2_Cou,Emp2_ex,&
                                kind_of,basis_parameter,&
                                qs_env,matrix_ks,rho,hfx_sections,para_env, &
                                unit_nr,error=error)


       END IF

     CASE (mp2_method_gpw)
       ! free Hartree Fock stuff
       CALL hfx_release(qs_env%x_data, error=error)
       mem=m_memory()
       mem_real=(mem+1024*1024-1)/(1024*1024)
       CALL mp_max(mem_real,para_env%group)
       mp2_env%mp2_memory=mp2_env%mp2_memory-mem_real
       IF(mp2_env%mp2_memory<0.0_dp) mp2_env%mp2_memory=1.0_dp

       IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T69,F9.2,A3)') 'Available memory per MPI processes for MP2:',&
                                                       mp2_env%mp2_memory, ' MB'
       ! go with mp2_gpw
       CALL  mp2_gpw_main(qs_env,mp2_env,Emp2,Emp2_Cou,Emp2_EX,para_env,unit_nr,error)

     CASE DEFAULT
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
   END SELECT
   t2=m_walltime()
   IF (unit_nr>0) WRITE(unit_nr,*)
   IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T56,F25.6)')  'Time:',t2-t1
   IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T56,F25.14)') 'Direct MP2 Coulomb Energy = ', Emp2_Cou/2.0_dp
   IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T56,F25.14)') 'Direct MP2 Exchange Energy = ', Emp2_ex
   IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T56,F25.14)') 'Direct MP2 Canonical Energy = ', Emp2
   IF (unit_nr>0) WRITE(unit_nr,*)

   ! we have it !!!!
   energy%mp2=Emp2
   energy%total=energy%total+Emp2

   CALL cp_print_key_finished_output(unit_nr,logger,input,&
                                     "DFT%XC%MP2%MP2_INFO", error=error)

   CALL timestop(handle)

  END SUBROUTINE mp2_main

  SUBROUTINE build_index_table(natom,max_nset,index_table,basis_parameter,kind_of)
    INTEGER                                  :: natom, max_nset
    INTEGER, DIMENSION(natom, max_nset)      :: index_table
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    INTEGER, DIMENSION(natom)                :: kind_of

    INTEGER                                  :: counter, iatom, ikind, iset, &
                                                nset

  index_table=-HUGE(0)
  counter=0
  DO iatom=1, natom
    ikind=kind_of(iatom)
    nset = basis_parameter(ikind)%nset
    DO iset=1, nset
      index_table(iatom,iset)=counter+1
      counter=counter+basis_parameter(ikind)%nsgf(iset)
    END DO
  END DO


  END SUBROUTINE build_index_table

  SUBROUTINE mp2_direct_energy(dimen,occ_i,occ_j,mp2_biel,mp2_env,C_i,Auto_i,Emp2,Emp2_Cou,Emp2_ex,&
                               kind_of,basis_parameter,&
                               qs_env,matrix_ks,rho,hfx_sections,para_env, &
                               unit_nr,C_j,Auto_j,error)
    INTEGER                                  :: dimen, occ_i, occ_j
    TYPE(mp2_biel_type)                      :: mp2_biel
    TYPE(mp2_type), POINTER                  :: mp2_env
    REAL(KIND=dp), DIMENSION(dimen, dimen)   :: C_i
    REAL(KIND=dp), DIMENSION(dimen)          :: Auto_i
    REAL(KIND=dp)                            :: Emp2, Emp2_Cou, Emp2_ex
    INTEGER, DIMENSION(:)                    :: kind_of
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: hfx_sections
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: unit_nr
    REAL(KIND=dp), DIMENSION(dimen, dimen), &
      OPTIONAL                               :: C_j
    REAL(KIND=dp), DIMENSION(dimen), &
      OPTIONAL                               :: Auto_j
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_direct_energy', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: due = 2.D+00, zero = 0.D+00

    INTEGER :: a, batch_number, color_sub, comm_sub, counter, &
      elements_ij_proc, group_counter, handle, i, i_batch, i_batch_start, &
      i_group_counter, j, j_batch_start, j_group_counter, last_batch, &
      max_batch_number, max_batch_size, max_set, minimum_memory_needed, &
      my_batch_size, my_I_batch_size, my_I_occupied_end, my_I_occupied_start, &
      my_J_batch_size, my_J_occupied_end, my_J_occupied_start, natom, &
      Ni_occupied, Nj_occupied, number_groups, number_i_subset, &
      number_j_subset, one, sqrt_number_groups, stat, &
      total_I_size_batch_group, total_J_size_batch_group, virt_i, virt_j
    INTEGER, ALLOCATABLE, DIMENSION(:) :: batch_sizes, batch_sizes_tmp, &
      vector_batch_I_size_group, vector_batch_J_size_group
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: ij_list_proc, &
                                                ij_list_proc_temp, ij_matrix
    LOGICAL                                  :: alpha_bata_case = .FALSE., &
                                                failure
    TYPE(cp_para_env_type), POINTER          :: para_env_sub

    CALL timeset(routineN,handle)
    failure=.FALSE.

    IF(PRESENT(C_j).AND.PRESENT(Auto_j)) alpha_bata_case=.TRUE.

    IF (unit_nr>0.AND.mp2_env%potential_parameter%potential_type==do_mp2_potential_TShPSC) THEN
      IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T64,F12.6,A5)') 'Truncated MP2 method, Rt=',&
                                                     mp2_env%potential_parameter%truncation_radius,' Bohr'
    END IF

    ! create the local para env
    ! each para_env_sub corresponds to a group that is going to compute
    ! all the integrals. To each group a batch I is assigned and the
    ! communication takes place only inside the group
    number_groups=para_env%num_pe/mp2_env%mp2_num_proc
    IF(number_groups*mp2_env%mp2_num_proc/=para_env%num_pe) THEN
      CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
          " The number of processors needs to be a multiple of the processors per group. "//&
          CPSourceFileRef,&
          only_ionode=.TRUE.)
    END IF
    IF(number_groups>occ_i*occ_j) THEN
      IF (unit_nr>0) WRITE(unit_nr,'(T3,A)') 'Number of groups greater then the number of IJ pairs!'
      IF (unit_nr>0) WRITE(unit_nr,'(T3,A)') 'Consider using more processors per group for improved efficiency'
    END IF

    color_sub=para_env%mepos/mp2_env%mp2_num_proc
    CALL mp_comm_split_direct(para_env%group,comm_sub,color_sub)
    NULLIFY(para_env_sub)
    CALL cp_para_env_create(para_env_sub,comm_sub,error=error)

    ! calculate the maximal size of the batch, according to the maximum RS size
    max_set=SIZE(mp2_biel%index_table,2)
    minimum_memory_needed=(8*(max_set**4))/1024**2
    IF(minimum_memory_needed>mp2_env%mp2_memory) THEN
      IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T67,F12.6,A2)') 'Memory required below the minimum, new memory:',&
                                                          minimum_memory_needed,'MB'
      mp2_env%mp2_memory=minimum_memory_needed
    END IF

    ! Distribute the batches over the groups in
    ! a rectangular fashion, bigger size for J index
    ! the sizes of the I batches should be as small as possible
    sqrt_number_groups=INT(SQRT(REAL(number_groups,KIND=dp)))
    DO i=1, number_groups
      IF(MOD(number_groups,i)==0) THEN
        IF(sqrt_number_groups/i<=1) THEN
          number_j_subset=i
          EXIT
        END IF
      END IF
    END DO
    number_i_subset=number_groups/number_j_subset

    IF(number_i_subset<number_j_subset) THEN
      number_i_subset=number_j_subset
      number_j_subset=number_groups/number_i_subset
    END IF

    ! Distribute the I index and the J index over groups
    total_I_size_batch_group=occ_i/number_i_subset
    IF(total_I_size_batch_group<1) total_I_size_batch_group=1
    ALLOCATE(vector_batch_I_size_group(0:number_i_subset-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    vector_batch_I_size_group=0
    DO i=0, number_i_subset-1
      vector_batch_I_size_group(i)=total_I_size_batch_group
    END DO
    IF(SUM(vector_batch_I_size_group)/=occ_i) THEN
      one=1
      IF(SUM(vector_batch_I_size_group)>occ_i) one=-1
        i=-1
        DO
          i=i+1
          vector_batch_I_size_group(i)=vector_batch_I_size_group(i)+one
          IF(SUM(vector_batch_I_size_group)==occ_i) EXIT
          IF(i==number_i_subset-1) i=-1
        END DO
    END IF

    total_J_size_batch_group=occ_j/number_j_subset
    IF(total_J_size_batch_group<1) total_J_size_batch_group=1
    ALLOCATE(vector_batch_J_size_group(0:number_j_subset-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    vector_batch_J_size_group=0
    DO i=0, number_J_subset-1
      vector_batch_J_size_group(i)=total_J_size_batch_group
    END DO
    IF(SUM(vector_batch_J_size_group)/=occ_j) THEN
      one=1
      IF(SUM(vector_batch_J_size_group)>occ_j) one=-1
        i=-1
        DO
          i=i+1
          vector_batch_J_size_group(i)=vector_batch_J_size_group(i)+one
          IF(SUM(vector_batch_J_size_group)==occ_j) EXIT
          IF(i==number_J_subset-1) i=-1
        END DO
    END IF

    ! now the starting and ending I and J occupied orbitals are assigned to each group
    group_counter=0
    i_group_counter=0
    my_I_occupied_start=1
    DO i=0, number_i_subset-1
      my_J_occupied_start=1
      j_group_counter=0
      DO j=0, number_j_subset-1
        group_counter=group_counter+1
        IF(color_sub==group_counter-1) EXIT
        my_J_occupied_start=my_J_occupied_start+vector_batch_J_size_group(j)
        j_group_counter=j_group_counter+1
      END DO
      IF(color_sub==group_counter-1) EXIT
      my_I_occupied_start=my_I_occupied_start+vector_batch_I_size_group(i)
      i_group_counter=i_group_counter+1
    END DO
    my_I_occupied_end=my_I_occupied_start+vector_batch_I_size_group(i_group_counter)-1
    my_I_batch_size=vector_batch_I_size_group(i_group_counter)
    my_J_occupied_end=my_J_occupied_start+vector_batch_J_size_group(j_group_counter)-1
    my_J_batch_size=vector_batch_J_size_group(j_group_counter)

    DEALLOCATE(vector_batch_I_size_group)
    DEALLOCATE(vector_batch_J_size_group)


    max_batch_size=MIN(  &
                   MAX(1,&
                   INT(mp2_env%mp2_memory*INT(1024,KIND=int_8)**2/&
                       (8*(2*dimen-occ_i)*INT(dimen,KIND=int_8)*my_J_batch_size/para_env_sub%num_pe))) &
                   ,my_I_batch_size)
    IF(max_batch_size<1) THEN
       max_batch_size=(8*(occ_i+1)*INT(dimen,KIND=int_8)**2/para_env%num_pe)/1024**2
       IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T73,I6,A2)') 'More memory required, at least:',max_batch_size,'MB'
       max_batch_size=1
    END IF

    ! create the size of the batches inside the group
    my_batch_size=my_I_batch_size
    ALLOCATE(batch_sizes(my_batch_size),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    batch_sizes=-HUGE(0)
    batch_number=0
    DO i=1, my_batch_size
      IF(i*max_batch_size>my_batch_size) EXIT
      batch_number=batch_number+1
      batch_sizes(i)=max_batch_size
    END DO
    last_batch=my_batch_size-max_batch_size*batch_number
    IF(last_batch>0) THEN
      batch_number=batch_number+1
      batch_sizes(batch_number)=last_batch
    END IF

    ALLOCATE(batch_sizes_tmp(batch_number),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    batch_sizes_tmp(1:batch_number)=batch_sizes(1:batch_number)
    DEALLOCATE(batch_sizes)
    ALLOCATE(batch_sizes(batch_number),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    batch_sizes=batch_sizes_tmp
    DEALLOCATE(batch_sizes_tmp)

    max_batch_size=MAXVAL(batch_sizes)
    CALL mp_max(max_batch_size,para_env%group)
    max_batch_number=batch_number
    CALL mp_max(max_batch_number,para_env%group)
    IF (unit_nr>0) THEN
      WRITE(unit_nr,'(T3,A,T76,I5)') 'Maximum used batch size: ',max_batch_size
      WRITE(unit_nr,'(T3,A,T76,I5)') 'Number of integral recomputations: ',max_batch_number
      CALL m_flush(unit_nr)
    END IF

    ! Batches sizes exceed the occupied orbitals allocated for group
    CPPostcondition(SUM(batch_sizes)<=my_batch_size,cp_failure_level,routineP,error,failure)

    virt_i=dimen-occ_i
    virt_j=dimen-occ_j
    natom=SIZE(mp2_biel%index_table,1)

    CALL mp_sync(para_env%group)
    Emp2=zero
    Emp2_Cou=zero
    Emp2_ex=zero
    i_batch_start=my_I_occupied_start-1
    j_batch_start=my_J_occupied_start-1
    Nj_occupied=my_J_batch_size
    DO i_batch=1, batch_number

         Ni_occupied=batch_sizes(i_batch)

         counter=-1
         ALLOCATE(ij_matrix(Ni_occupied,Nj_occupied),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

         ij_matrix=0
         DO i=1, Ni_occupied
           DO j=1, Nj_occupied
             counter=counter+1
             IF (MOD(counter,para_env_sub%num_pe)==para_env_sub%mepos) THEN
               ij_matrix(i,j)=ij_matrix(i,j)+1
             END IF
           END DO
         END DO

         ALLOCATE(ij_list_proc_temp(Ni_occupied*occ_j,2),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

         elements_ij_proc=0
         DO i=1, Ni_occupied
           DO j=1, Nj_occupied
             IF(ij_matrix(i,j)==0) CYCLE
             elements_ij_proc=elements_ij_proc+1
             ij_list_proc_temp(elements_ij_proc,1)=i
             ij_list_proc_temp(elements_ij_proc,2)=j
           END DO
         END DO
         DEALLOCATE(ij_matrix)

         ALLOCATE(ij_list_proc(elements_ij_proc,2))
         DO i=1, elements_ij_proc
           ij_list_proc(i,1)=ij_list_proc_temp(i,1)
           ij_list_proc(i,2)=ij_list_proc_temp(i,2)
         END DO
         DEALLOCATE(ij_list_proc_temp)

         IF(.NOT.alpha_bata_case) THEN
           CALL mp2_canonical_direct_single_batch(Emp2,Emp2_Cou,Emp2_ex,mp2_env,qs_env,rho,hfx_sections,para_env_sub,&
                                        mp2_biel,dimen,C_i,Auto_i,i_batch_start,Ni_occupied,occ_i,&
                                        elements_ij_proc, ij_list_proc,Nj_occupied,j_batch_start,&
                                        error=error)
         ELSE
           CALL mp2_canonical_direct_single_batch(Emp2,Emp2_Cou,Emp2_ex,mp2_env,qs_env,rho,hfx_sections,para_env_sub,&
                                        mp2_biel,dimen,C_i,Auto_i,i_batch_start,Ni_occupied,occ_i,&
                                        elements_ij_proc, ij_list_proc,Nj_occupied,j_batch_start,&
                                        occ_j,C_j,Auto_j,error=error)
         END IF

         i_batch_start=i_batch_start+Ni_occupied

         DEALLOCATE(ij_list_proc)

    END DO

    CALL mp_sum(Emp2_Cou,para_env%group)
    CALL mp_sum(Emp2_Ex,para_env%group)
    CALL mp_sum(Emp2,para_env%group)

    CALL cp_para_env_release(para_env_sub,error)

    CALL timestop(handle)

  END SUBROUTINE  mp2_direct_energy

END MODULE mp2


