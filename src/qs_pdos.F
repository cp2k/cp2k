!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation and writing of projected density of  states
!>         The DOS is computed per angular momentum and per kind
!> \par History
!>      -
!> \author Marcella (29.02.2008,MK)
! *****************************************************************************
MODULE qs_pdos
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
  USE cp_fm_diag,                      ONLY: cp_fm_power
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_init_random,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: nso,&
                                             nsoset
  USE orbital_symbols,                 ONLY: l_sym,&
                                             sgf_symbol
  USE particle_types,                  ONLY: particle_type
  USE preconditioner_types,            ONLY: preconditioner_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_p_type,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE qs_collocate_density,            ONLY: calculate_wavefunction
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             get_qs_kind_set,&
                                             qs_kind_type
  USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_type
  USE qs_ot_eigensolver,               ONLY: ot_eigensolver
  USE scf_control_types,               ONLY: scf_control_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_pdos'

! *****************************************************************************
 ! *** Public subroutines ***

  PUBLIC :: calculate_projected_dos


  TYPE ldos_type
     INTEGER :: maxl, nlist
     LOGICAL :: separate_components
     INTEGER, DIMENSION(:), POINTER           :: list_index
     REAL(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: pdos_array
  END TYPE ldos_type

  TYPE r_ldos_type
     INTEGER :: nlist, npoints
     INTEGER, DIMENSION(:,:), POINTER :: index_grid_local
     INTEGER, DIMENSION(:), POINTER           :: list_index
     REAL(KIND=dp), DIMENSION(:), POINTER     :: x_range,y_range,z_range
     REAL(KIND=dp), DIMENSION(:), POINTER     :: eval_range
     REAL(KIND=dp), DIMENSION(:), &
      POINTER                                 :: pdos_array
  END TYPE r_ldos_type

  TYPE ldos_p_type
     TYPE(ldos_type), POINTER :: ldos
  END TYPE ldos_p_type

  TYPE r_ldos_p_type
     TYPE(r_ldos_type), POINTER :: ldos
  END TYPE r_ldos_p_type
CONTAINS

! *****************************************************************************
!> \brief   Compute and write projected density of states
!> \param mo_set ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param qs_env ...
!> \param dft_section ...
!> \param ispin ...
!> \param error ...
!> \date    26.02.2008
!> \par History:
!>       -
!> \par Variables
!>       -
!>       -
!> \author  MI
!> \version 1.0
! *****************************************************************************
  SUBROUTINE calculate_projected_dos(mo_set,atomic_kind_set,qs_kind_set,particle_set,qs_env,&
                                 dft_section,ispin,error)

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: dft_section
    INTEGER, INTENT(IN), OPTIONAL            :: ispin
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_projected_dos', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=16)                        :: fmtstr2
    CHARACTER(LEN=27)                        :: fmtstr1
    CHARACTER(LEN=6), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: tmp_str
    CHARACTER(LEN=default_string_length)     :: kind_name, my_act, my_mittle, &
                                                my_pos, spin(2)
    CHARACTER(LEN=default_string_length), &
      ALLOCATABLE, DIMENSION(:)              :: ldos_index, r_ldos_index
    INTEGER :: handle, homo, i, iatom, ikind, il, ildos, im, imo, in_x, in_y, &
      in_z, ir, irow, iset, isgf, ishell, iso, iterstep, iw, j, jx, jy, jz, &
      k, lcomponent, lshell, maxl, maxlgto, my_spin, n_dependent, n_r_ldos, &
      n_rep, nao, natom, ncol_global, nkind, nldos, nlumo, nmo, np_tot, &
      npoints, nrow_global, nset, nsgf, nvirt, out_each, output_unit, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: firstrow
    INTEGER, DIMENSION(:), POINTER           :: list, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: bo, l
    LOGICAL                                  :: append, do_ldos, do_r_ldos, &
                                                do_virt, failure, ionode, &
                                                separate_components, &
                                                should_output
    LOGICAL, DIMENSION(:, :), POINTER        :: read_r
    REAL(KIND=dp)                            :: dh(3,3), dvol, e_fermi, r(3), &
                                                r_vec(3), ratom(3)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eigenvalues, evals_virt, &
                                                occupation_numbers
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vecbuffer
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: pdos_array
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: s_matrix
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: matrix_shalf, matrix_shalfc, &
                                                matrix_work, mo_coeff, mo_virt
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(ldos_p_type), DIMENSION(:), POINTER :: ldos_p
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: wf_g, wf_r
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(r_ldos_p_type), DIMENSION(:), &
      POINTER                                :: r_ldos_p
    TYPE(section_vals_type), POINTER         :: ldos_section

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    failure = .FALSE.
    should_output = BTEST(cp_print_key_should_output(logger%iter_info,dft_section,&
                          "PRINT%PDOS",error=error),cp_p_file)
    output_unit= cp_logger_get_default_io_unit(logger)

    spin(1)="ALPHA"
    spin(2)="BETA"
    IF ((.NOT.should_output)) RETURN

    NULLIFY (context, s_matrix, orb_basis_set, para_env, pdos_array)
    NULLIFY (eigenvalues,fm_struct_tmp,matrix_work,matrix_shalf,mo_coeff,vecbuffer)
    NULLIFY (ldos_section,list,cell,pw_env,auxbas_pw_pool,mo_virt,evals_virt)
    NULLIFY (occupation_numbers, ldos_p, r_ldos_p,dft_control,occupation_numbers )

    CALL timeset(routineN,handle)
    iterstep = logger%iter_info%iteration(logger%iter_info%n_rlevel)

    IF (output_unit > 0) WRITE(UNIT=output_unit,FMT='(/,(T3,A,T61,I10))')&
      " Calculate PDOS at iteration step ", iterstep
    CALL get_qs_env(qs_env=qs_env,&
         matrix_s=s_matrix,error=error)

    CALL get_atomic_kind_set(atomic_kind_set, natom=natom)
    CALL get_qs_kind_set(qs_kind_set, nsgf=nsgf, maxlgto=maxlgto, error=error)
    nkind = SIZE(atomic_kind_set)

    CALL get_mo_set(mo_set=mo_set, mo_coeff=mo_coeff, homo=homo, nao=nao, nmo=nmo,&
                    mu=e_fermi)
    CALL cp_fm_get_info(mo_coeff,&
         context=context, para_env=para_env,&
         nrow_global=nrow_global,&
         ncol_global=ncol_global,error=error)

    CALL section_vals_val_get(dft_section,"PRINT%PDOS%OUT_EACH_MO",i_val=out_each,error=error)
    IF(out_each==-1) out_each=nao+1
    CALL section_vals_val_get(dft_section,"PRINT%PDOS%NLUMO",i_val=nlumo,error=error)
    IF (nlumo == -1) nlumo = nao - homo
    do_virt = (nlumo>(nmo-homo))
    nvirt = nlumo - (nmo-homo)
    ! Generate virtual orbitals
    IF(do_virt) THEN
       IF(PRESENT(ispin)) THEN
         my_spin = ispin
       ELSE
         my_spin = 1
       END IF

       CALL generate_virtual_mo(qs_env,mo_set, evals_virt, mo_virt, nvirt, ispin=my_spin, error=error )
    ELSE
       NULLIFY (evals_virt,mo_virt)
       nvirt = 0
    END IF

    ! Create S^1/2 : from sparse to full matrix
    CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=context, &
                  nrow_global=nrow_global,ncol_global=nrow_global,error=error)
    CALL cp_fm_create(matrix_shalf, fm_struct_tmp,name="matrix_shalf",error=error)
    CALL cp_fm_create(matrix_work, fm_struct_tmp,name="matrix_work",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    CALL copy_dbcsr_to_fm(s_matrix(1)%matrix,matrix_shalf,error=error)
    CALL cp_fm_power(matrix_shalf,matrix_work,0.5_dp,EPSILON(0.0_dp),n_dependent,error=error)
    CALL cp_fm_release(matrix_work,error=error)

    ! Multiply S^(1/2) time the mOS coefficients to get orthonormalized MOS
    CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=context, &
         nrow_global=nrow_global,ncol_global=ncol_global,error=error)
    CALL cp_fm_create(matrix_shalfc,fm_struct_tmp,name="matrix_shalfc",error=error)
    CALL cp_gemm("N","N",nrow_global,ncol_global,nrow_global, &
                   1.0_dp,matrix_shalf,mo_coeff,0.0_dp,matrix_shalfc,error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)

    IF (do_virt) THEN
      IF (output_unit > 0) WRITE(UNIT=output_unit,FMT='(/,(T3,A,T14,I10,T27,A))')&
      " Compute ", nvirt, " additional unoccupied KS orbitals"
      CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=context, &
            nrow_global=nrow_global,ncol_global=nvirt,error=error)
      CALL cp_fm_create(matrix_work,fm_struct_tmp,name="matrix_shalfc",error=error)
      CALL cp_gemm("N","N",nrow_global,nvirt,nrow_global, &
                     1.0_dp,matrix_shalf,mo_virt,0.0_dp,matrix_work,error=error)
      CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    END IF

    CALL cp_fm_release(matrix_shalf,error=error)
    ! Array to store the PDOS per kind and angular momentum
    do_ldos = .FALSE.
    ldos_section => section_vals_get_subs_vals(dft_section,"PRINT%PDOS%LDOS",error=error)

    CALL section_vals_get(ldos_section,n_repetition=nldos,error=error)
    IF(nldos>0) THEN
      IF (output_unit > 0) WRITE(UNIT=output_unit,FMT='(/,(T3,A,T61,I10))')&
      " Prepare the list of atoms for LDOS.   Number of lists: ", nldos
      do_ldos = .TRUE.
      ALLOCATE(ldos_p(nldos),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(ldos_index(nldos),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DO ildos=1,nldos
        WRITE(ldos_index(ildos),'(I0)') ildos
        ALLOCATE(ldos_p(ildos)%ldos,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        NULLIFY(ldos_p(ildos)%ldos%pdos_array)
        NULLIFY(ldos_p(ildos)%ldos%list_index)

        CALL section_vals_val_get(ldos_section,"LIST",i_rep_section=ildos, n_rep_val=n_rep,&
             error=error)
        IF(n_rep>0) THEN
          ldos_p(ildos)%ldos%nlist = 0
          DO ir = 1,n_rep
            NULLIFY (list)
            CALL section_vals_val_get(ldos_section,"LIST",i_rep_section=ildos,i_rep_val=ir,&
                 i_vals=list, error=error)
            IF(ASSOCIATED(list)) THEN
               CALL reallocate(ldos_p(ildos)%ldos%list_index,1,ldos_p(ildos)%ldos%nlist+SIZE(list))
               DO i = 1,SIZE(list)
                  ldos_p(ildos)%ldos%list_index(i+ldos_p(ildos)%ldos%nlist) = list(i)
               END DO
               ldos_p(ildos)%ldos%nlist = ldos_p(ildos)%ldos%nlist + SIZE(list)
            END IF
          END DO
        ELSE
          ! stop, LDOS without list of atoms is not implemented
        END IF
        
        IF (output_unit > 0) WRITE(UNIT=output_unit,FMT='((T10,A,T18,I6,T25,A,T36,I10,A))')&
         " List ",ildos , " contains ",  ldos_p(ildos)%ldos%nlist , " atoms"
        CALL section_vals_val_get(ldos_section,"COMPONENTS",i_rep_section=ildos,&
             l_val=ldos_p(ildos)%ldos%separate_components,error=error)
        IF(ldos_p(ildos)%ldos%separate_components) THEN
          ALLOCATE(ldos_p(ildos)%ldos%pdos_array(nsoset(maxlgto),nmo+nvirt),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ELSE
          ALLOCATE(ldos_p(ildos)%ldos%pdos_array(0:maxlgto,nmo+nvirt),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END IF
        ldos_p(ildos)%ldos%pdos_array=0.0_dp
        ldos_p(ildos)%ldos%maxl = -1

      END DO
    END IF

    do_r_ldos = .FALSE.
    ldos_section => section_vals_get_subs_vals(dft_section,"PRINT%PDOS%R_LDOS",error=error)
    CALL section_vals_get(ldos_section,n_repetition=n_r_ldos,error=error)
    IF(n_r_ldos>0) THEN
      do_r_ldos = .TRUE.
      IF (output_unit > 0) WRITE(UNIT=output_unit,FMT='(/,(T3,A,T61,I10))')&
      " Prepare the list of points for R_LDOS.   Number of lists: ", n_r_ldos
      ALLOCATE(r_ldos_p(n_r_ldos),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(r_ldos_index(n_r_ldos),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      CALL get_qs_env(qs_env=qs_env,&
           cell=cell,&
           dft_control=dft_control,&
           pw_env=pw_env,error=error)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
                      pw_pools=pw_pools,error=error)

      CALL pw_pool_create_pw(auxbas_pw_pool,wf_r%pw,&
                              use_data = REALDATA3D,&
                              in_space = REALSPACE, error=error)
      CALL pw_pool_create_pw(auxbas_pw_pool,wf_g%pw,&
                              use_data = COMPLEXDATA1D,&
                              in_space = RECIPROCALSPACE, error=error)
      ALLOCATE(read_r(4,n_r_ldos),STAT=stat)
      DO ildos=1,n_r_ldos
        WRITE(r_ldos_index(ildos),'(I0)') ildos
        ALLOCATE(r_ldos_p(ildos)%ldos,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        NULLIFY(r_ldos_p(ildos)%ldos%pdos_array)
        NULLIFY(r_ldos_p(ildos)%ldos%list_index)

        CALL section_vals_val_get(ldos_section,"LIST",i_rep_section=ildos, n_rep_val=n_rep,&
             error=error)
        IF(n_rep>0) THEN
          r_ldos_p(ildos)%ldos%nlist = 0
          DO ir = 1,n_rep
            NULLIFY (list)
            CALL section_vals_val_get(ldos_section,"LIST",i_rep_section=ildos,i_rep_val=ir,&
                 i_vals=list, error=error)
            IF(ASSOCIATED(list)) THEN
               CALL reallocate(r_ldos_p(ildos)%ldos%list_index,1,r_ldos_p(ildos)%ldos%nlist+SIZE(list))
               DO i = 1,SIZE(list)
                  r_ldos_p(ildos)%ldos%list_index(i+r_ldos_p(ildos)%ldos%nlist) = list(i)
               END DO
               r_ldos_p(ildos)%ldos%nlist = r_ldos_p(ildos)%ldos%nlist + SIZE(list)
            END IF
          END DO
        ELSE
          ! stop, LDOS without list of atoms is not implemented
        END IF

        ALLOCATE(r_ldos_p(ildos)%ldos%pdos_array(nmo+nvirt),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        r_ldos_p(ildos)%ldos%pdos_array=0.0_dp
        read_r(1:3, ildos) = .FALSE.
        CALL section_vals_val_get(ldos_section,"XRANGE",i_rep_section=ildos,explicit=read_r(1,ildos),error=error)
        IF (read_r(1,ildos)) THEN
          CALL section_vals_val_get(ldos_section,"XRANGE",i_rep_section=ildos,r_vals=&
             r_ldos_p(ildos)%ldos%x_range, error=error)
        ELSE
          ALLOCATE(r_ldos_p(ildos)%ldos%x_range(2))
           r_ldos_p(ildos)%ldos%x_range=0.0_dp
        END IF
        CALL section_vals_val_get(ldos_section,"YRANGE",i_rep_section=ildos,explicit=read_r(2,ildos),error=error)
        IF (read_r(2,ildos)) THEN
          CALL section_vals_val_get(ldos_section,"YRANGE",i_rep_section=ildos,r_vals=&
             r_ldos_p(ildos)%ldos%y_range, error=error)
        ELSE
          ALLOCATE(r_ldos_p(ildos)%ldos%y_range(2))
           r_ldos_p(ildos)%ldos%y_range=0.0_dp
        END IF
        CALL section_vals_val_get(ldos_section,"ZRANGE",i_rep_section=ildos,explicit=read_r(3,ildos),error=error)
        IF (read_r(3,ildos)) THEN
          CALL section_vals_val_get(ldos_section,"ZRANGE",i_rep_section=ildos,r_vals=&
             r_ldos_p(ildos)%ldos%z_range, error=error)
        ELSE
          ALLOCATE(r_ldos_p(ildos)%ldos%z_range(2))
           r_ldos_p(ildos)%ldos%z_range=0.0_dp
        END IF

        CALL section_vals_val_get(ldos_section,"ERANGE",i_rep_section=ildos,explicit=read_r(4,ildos),error=error)
        IF (read_r(4,ildos)) THEN
          CALL section_vals_val_get(ldos_section,"ERANGE",i_rep_section=ildos,&
               r_vals=r_ldos_p(ildos)%ldos%eval_range, error=error)
        ELSE
          ALLOCATE(r_ldos_p(ildos)%ldos%eval_range(2))
          r_ldos_p(ildos)%ldos%eval_range(1)=-HUGE(0.0_dp)
          r_ldos_p(ildos)%ldos%eval_range(2)=+HUGE(0.0_dp)
        END IF

        bo => wf_r%pw%pw_grid%bounds_local
        dh = wf_r%pw%pw_grid%dh
        dvol = wf_r%pw%pw_grid%dvol
        np_tot = wf_r%pw%pw_grid%npts(1)*wf_r%pw%pw_grid%npts(2)*wf_r%pw%pw_grid%npts(3)
        ALLOCATE(r_ldos_p(ildos)%ldos%index_grid_local(3,np_tot), STAT =stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

        r_ldos_p(ildos)%ldos%npoints = 0
        DO jz=bo(1,3),bo(2,3)
        DO jy=bo(1,2),bo(2,2)
        DO jx=bo(1,1),bo(2,1)
          !compute the position of the grid point
          i = jx -wf_r%pw%pw_grid%bounds(1,1)
          j = jy -wf_r%pw%pw_grid%bounds(1,2)
          k = jz -wf_r%pw%pw_grid%bounds(1,3)
          r(3)=k*dh(3,3)+j*dh(3,2)+i*dh(3,1)
          r(2)=k*dh(2,3)+j*dh(2,2)+i*dh(2,1)
          r(1)=k*dh(1,3)+j*dh(1,2)+i*dh(1,1)
           
          DO il = 1, r_ldos_p(ildos)%ldos%nlist
             iatom = r_ldos_p(ildos)%ldos%list_index(il)
             ratom = particle_set(iatom)%r
             r_vec = pbc(ratom,r,cell)
             IF (cell%orthorhombic) THEN
               IF(cell%perd(1)==0) r_vec(1)=MODULO(r_vec(1),cell%hmat(1,1))
               IF(cell%perd(2)==0) r_vec(2)=MODULO(r_vec(2),cell%hmat(2,2))
               IF(cell%perd(3)==0) r_vec(3)=MODULO(r_vec(3),cell%hmat(3,3))
             END IF

             in_x =0
             in_y =0
             in_z =0
             IF( r_ldos_p(ildos)%ldos%x_range(1)/=0.0_dp) THEN
                IF(r_vec(1)>r_ldos_p(ildos)%ldos%x_range(1) .AND. &
                   r_vec(1)<r_ldos_p(ildos)%ldos%x_range(2)) THEN
                   in_x = 1
                END IF
             ELSE
               in_x = 1
             END IF
             IF( r_ldos_p(ildos)%ldos%y_range(1)/=0.0_dp) THEN
                IF(r_vec(2)>r_ldos_p(ildos)%ldos%y_range(1) .AND. &
                   r_vec(2)<r_ldos_p(ildos)%ldos%y_range(2)) THEN
                   in_y = 1
                END IF
             ELSE
               in_y = 1
             END IF
             IF( r_ldos_p(ildos)%ldos%z_range(1)/=0.0_dp) THEN
                IF(r_vec(3)>r_ldos_p(ildos)%ldos%z_range(1) .AND. &
                   r_vec(3)<r_ldos_p(ildos)%ldos%z_range(2)) THEN
                   in_z = 1
                END IF
             ELSE
               in_z = 1
             END IF
             IF(in_x*in_y*in_z>0) THEN
               r_ldos_p(ildos)%ldos% npoints = r_ldos_p(ildos)%ldos%npoints+1
               r_ldos_p(ildos)%ldos%index_grid_local(1,r_ldos_p(ildos)%ldos%npoints) = jx
               r_ldos_p(ildos)%ldos%index_grid_local(2,r_ldos_p(ildos)%ldos%npoints) = jy
               r_ldos_p(ildos)%ldos%index_grid_local(3,r_ldos_p(ildos)%ldos%npoints) = jz
               EXIT
             END IF
          END DO
        END DO
        END DO
        END DO
        CALL  reallocate(r_ldos_p(ildos)%ldos%index_grid_local,1,3,1,r_ldos_p(ildos)%ldos%npoints)
        npoints =  r_ldos_p(ildos)%ldos%npoints
        CALL mp_sum(npoints,para_env%group)
        IF (output_unit > 0) WRITE(UNIT=output_unit,FMT='((T10,A,T18,I6,T25,A,T36,I10,A))')&
         " List ",ildos , " contains ", npoints , " grid points"
      END DO
    END IF


    CALL section_vals_val_get(dft_section,"PRINT%PDOS%COMPONENTS",l_val=separate_components,&
         error=error)
    IF (separate_components) THEN
      ALLOCATE(pdos_array(nsoset(maxlgto),nkind,nmo+nvirt),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
      ALLOCATE(pdos_array(0:maxlgto,nkind,nmo+nvirt),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(do_virt) THEN
       ALLOCATE(eigenvalues(nmo+nvirt),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       eigenvalues(1:nmo) = mo_set%eigenvalues(1:nmo)
       eigenvalues(nmo+1:nmo+nvirt) = evals_virt(1:nvirt)
       ALLOCATE(occupation_numbers(nmo+nvirt),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       occupation_numbers(:) = 0.0_dp
       occupation_numbers(1:nmo) = mo_set%occupation_numbers(1:nmo)
    ELSE
       eigenvalues => mo_set%eigenvalues
       occupation_numbers => mo_set%occupation_numbers
    END IF

    pdos_array = 0.0_dp
    nao = mo_set%nao
    ALLOCATE(vecbuffer(1,nao),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    vecbuffer = 0.0_dp
    ALLOCATE(firstrow(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    firstrow=0

    !Adjust energy range for r_ldos
    DO ildos=1,n_r_ldos
       IF(eigenvalues(1)>r_ldos_p(ildos)%ldos%eval_range(1)) & 
            r_ldos_p(ildos)%ldos%eval_range(1) = eigenvalues(1)
       IF(eigenvalues(nmo+nvirt)<r_ldos_p(ildos)%ldos%eval_range(2)) &
            r_ldos_p(ildos)%ldos%eval_range(2) = eigenvalues(nmo+nvirt)
    END DO
          
    IF (output_unit > 0) WRITE(UNIT=output_unit,FMT='(/,(T15,A))')&
             "---- PDOS: start iteration on the KS states --- "

    DO imo =  1,nmo+nvirt

      IF (output_unit > 0 .AND. MOD(imo,out_each)==0) WRITE(UNIT=output_unit,FMT='((T20,A,I10))')&
               " KS state index : ", imo
      ! Extract the eigenvector from the distributed full matrix
      IF(imo>nmo) THEN
        CALL cp_fm_get_submatrix(matrix_work,vecbuffer,1,imo-nmo,&
                       nao,1,transpose=.TRUE.,error=error)
      ELSE
        CALL cp_fm_get_submatrix(matrix_shalfc,vecbuffer,1,imo,&
                       nao,1,transpose=.TRUE.,error=error)
      END IF

      ! Calculate the pdos for all the kinds
      irow = 1
      DO iatom = 1,natom
          firstrow(iatom)=irow
          NULLIFY(orb_basis_set)
          CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
          CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set, error=error)
          CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                          nset=nset,&
                          nshell=nshell,&
                          l=l, maxl=maxl)
          IF (separate_components) THEN
            isgf = 1
            DO iset=1,nset
              DO ishell=1,nshell(iset)
                 lshell = l(ishell,iset)
                 DO iso=1,nso(lshell)
                   lcomponent = nsoset(lshell-1)+iso
                   pdos_array(lcomponent,ikind,imo) = &
                        pdos_array(lcomponent,ikind,imo) +&
                        vecbuffer(1,irow) * vecbuffer(1,irow)
                   irow = irow + 1
                 END DO  ! iso
              END DO  ! ishell
            END DO  ! iset
          ELSE
            isgf = 1
            DO iset=1,nset
              DO ishell=1,nshell(iset)
                 lshell = l(ishell,iset)
                 DO iso=1,nso(lshell)
                   pdos_array(lshell,ikind,imo) = &
                        pdos_array(lshell,ikind,imo) +&
                        vecbuffer(1,irow) * vecbuffer(1,irow)
                   irow = irow + 1
                 END DO  ! iso
              END DO  ! ishell
            END DO  ! iset
          END IF
      END DO  ! iatom

      ! Calculate the pdos for all the lists
      DO ildos=1,nldos
        DO il=1,ldos_p(ildos)%ldos%nlist
          iatom = ldos_p(ildos)%ldos%list_index(il)

          irow=firstrow(iatom)
          NULLIFY(orb_basis_set)
          CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
          CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set, error=error)

          CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                          nset=nset,&
                          nshell=nshell,&
                          l=l, maxl=maxl)
          ldos_p(ildos)%ldos%maxl = MAX(ldos_p(ildos)%ldos%maxl,maxl)
          IF(ldos_p(ildos)%ldos%separate_components) THEN
            isgf = 1
            DO iset=1,nset
              DO ishell=1,nshell(iset)
                 lshell = l(ishell,iset)
                 DO iso=1,nso(lshell)
                   lcomponent = nsoset(lshell-1)+iso
                   ldos_p(ildos)%ldos%pdos_array(lcomponent,imo) =&
                      ldos_p(ildos)%ldos%pdos_array(lcomponent,imo) + &
                      vecbuffer(1,irow) * vecbuffer(1,irow)
                   irow = irow + 1
                 END DO  ! iso
              END DO  ! ishell
            END DO  ! iset
          ELSE
            isgf = 1
            DO iset=1,nset
              DO ishell=1,nshell(iset)
                 lshell = l(ishell,iset)
                 DO iso=1,nso(lshell)
                   ldos_p(ildos)%ldos%pdos_array(lshell,imo) = &
                      ldos_p(ildos)%ldos%pdos_array(lshell,imo) + &
                       vecbuffer(1,irow) * vecbuffer(1,irow)
                  irow = irow + 1
                 END DO  ! iso
              END DO  ! ishell
            END DO  ! iset
          END IF
        END DO !il
      END DO !ildos


      ! Calculate the DOS projected in a given volume in real space
      DO ildos = 1,n_r_ldos
        IF(r_ldos_p(ildos)%ldos%eval_range(1)<=eigenvalues(imo) .AND. &
           r_ldos_p(ildos)%ldos%eval_range(2)>=eigenvalues(imo)) THEN

           IF(imo>nmo) THEN
             CALL calculate_wavefunction(mo_virt,imo-nmo,&
                    wf_r, wf_g, atomic_kind_set, qs_kind_set,cell,dft_control,particle_set, &
                  pw_env, error = error )
           ELSE
             CALL calculate_wavefunction(mo_coeff,imo,&
                   wf_r, wf_g, atomic_kind_set,qs_kind_set,cell,dft_control,particle_set, &
                  pw_env, error = error )
           END IF
           r_ldos_p(ildos)%ldos%pdos_array(imo)  = 0.0_dp
           DO il = 1,r_ldos_p(ildos)%ldos%npoints
              j=j+1
              jx=r_ldos_p(ildos)%ldos%index_grid_local(1,il)
              jy=r_ldos_p(ildos)%ldos%index_grid_local(2,il)
              jz=r_ldos_p(ildos)%ldos%index_grid_local(3,il)
              r_ldos_p(ildos)%ldos%pdos_array(imo) = r_ldos_p(ildos)%ldos%pdos_array(imo)+ &
                   wf_r%pw%cr3d(jx,jy,jz)* wf_r%pw%cr3d(jx,jy,jz)
           END DO
           r_ldos_p(ildos)%ldos%pdos_array(imo) = r_ldos_p(ildos)%ldos%pdos_array(imo)*dvol
        END IF
      END DO
    END DO  ! imo

    CALL cp_fm_release(matrix_shalfc,error=error)
    DEALLOCATE(vecbuffer,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)


    CALL section_vals_val_get(dft_section,"PRINT%PDOS%APPEND",l_val=append,error=error)
    IF(append .AND. iterstep > 1 ) THEN
      my_pos = "APPEND"
    ELSE
      my_pos = "REWIND"
    END IF
    my_act = "WRITE"
    DO ikind = 1,nkind

      NULLIFY(orb_basis_set)
      CALL get_atomic_kind(atomic_kind_set(ikind), name=kind_name)
      CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set, error=error)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_set, maxl=maxl)

      ! basis none has no associated maxl, and no pdos
      IF (maxl<0) CYCLE

      IF(PRESENT(ispin)) THEN
        my_mittle =TRIM(spin(ispin))//"_k"//TRIM(ADJUSTL(cp_to_string(ikind)))
        my_spin = ispin
      ELSE
        my_mittle ="k"//TRIM(ADJUSTL(cp_to_string(ikind)))
        my_spin = 1
      END IF

      iw = cp_print_key_unit_nr(logger,dft_section,"PRINT%PDOS",&
           extension=".pdos", file_position=my_pos, file_action=my_act,&
           file_form="FORMATTED", middle_name=TRIM(my_mittle), &
           error=error)
      IF(iw>0) THEN

        fmtstr1 = "(I8,2X,2F16.6,  (2X,F16.8))"
        fmtstr2 = "(A42,  (10X,A8))"
        IF(separate_components) THEN
          WRITE (UNIT=fmtstr1(15:16),FMT="(I2)") nsoset(maxl)
          WRITE (UNIT=fmtstr2(6:7),FMT="(I2)") nsoset(maxl)
        ELSE
          WRITE (UNIT=fmtstr1(15:16),FMT="(I2)") maxl+1
          WRITE (UNIT=fmtstr2(6:7),FMT="(I2)") maxl+1
        END IF

        WRITE (UNIT=iw,FMT="(A,I0,A,F12.6,A)")&
          "# Projected DOS for atomic kind "//TRIM(kind_name)//" at iteration step i = ",&
          iterstep,", E(Fermi) = ",e_fermi," a.u."
        IF (separate_components) THEN
          ALLOCATE(tmp_str(0:0,0:maxl,-maxl:maxl),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          tmp_str = ""
          DO j = 0, maxl
             DO i = -j,j
                tmp_str(0,j,i) = sgf_symbol(0,j,i)
             END DO
          END DO

          WRITE (UNIT=iw,FMT=fmtstr2)&
            "#     MO Eigenvalue [a.u.]      Occupation",&
            ((TRIM(tmp_str(0,il,im)),im=-il,il),il=0,maxl)
          DO imo = 1,nmo+nvirt
            WRITE (UNIT=iw,FMT=fmtstr1) imo,eigenvalues(imo),occupation_numbers(imo),&
              (pdos_array(lshell,ikind,imo),lshell=1,nsoset(maxl))
          END DO
          DEALLOCATE(tmp_str,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ELSE
          WRITE (UNIT=iw,FMT=fmtstr2)&
            "#     MO Eigenvalue [a.u.]      Occupation",&
            (TRIM(l_sym(il)),il=0,maxl)
          DO imo = 1,nmo+nvirt
            WRITE (UNIT=iw,FMT=fmtstr1) imo,eigenvalues(imo),occupation_numbers(imo),&
              (pdos_array(lshell,ikind,imo),lshell=0,maxl)
          END DO
        END IF
      END IF
      CALL cp_print_key_finished_output(iw,logger,dft_section,&
                "PRINT%PDOS", error=error)

    END DO  ! ikind


    ! write the pdos for the lists, each ona different file,
    ! the filenames are indexed with the list number
    DO ildos=1,nldos
     ! basis none has no associated maxl, and no pdos
      IF (ldos_p(ildos)%ldos%maxl>0) THEN

        IF(PRESENT(ispin)) THEN
          my_mittle =TRIM(spin(ispin))//"_list"//TRIM(ldos_index(ildos))
          my_spin = ispin
        ELSE
          my_mittle ="list"//TRIM(ldos_index(ildos))
          my_spin = 1
        END IF

        iw = cp_print_key_unit_nr(logger,dft_section,"PRINT%PDOS",&
             extension=".pdos", file_position=my_pos, file_action=my_act,&
             file_form="FORMATTED", middle_name=TRIM(my_mittle), &
             error=error)
        IF(iw>0) THEN

          fmtstr1 = "(I8,2X,2F16.6,  (2X,F16.8))"
          fmtstr2 = "(A42,  (10X,A8))"
          IF(ldos_p(ildos)%ldos%separate_components) THEN
            WRITE (UNIT=fmtstr1(15:16),FMT="(I2)") nsoset(ldos_p(ildos)%ldos%maxl)
            WRITE (UNIT=fmtstr2(6:7),FMT="(I2)") nsoset(ldos_p(ildos)%ldos%maxl)
          ELSE
            WRITE (UNIT=fmtstr1(15:16),FMT="(I2)") ldos_p(ildos)%ldos%maxl+1
            WRITE (UNIT=fmtstr2(6:7),FMT="(I2)") ldos_p(ildos)%ldos%maxl+1
          END IF

          WRITE (UNIT=iw,FMT="(A,I0,A,I0,A,I0,A,F12.6,A)")&
            "# Projected DOS for list ", ildos, " of ",ldos_p(ildos)%ldos%nlist," atoms, at iteration step i = ",iterstep,&
             ", E(Fermi) = ",e_fermi," a.u."
          IF(ldos_p(ildos)%ldos%separate_components) THEN
            ALLOCATE(tmp_str(0:0,0:ldos_p(ildos)%ldos%maxl,-ldos_p(ildos)%ldos%maxl:ldos_p(ildos)%ldos%maxl),stat=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
            tmp_str = ""
            DO j = 0, ldos_p(ildos)%ldos%maxl
               DO i = -j,j
                  tmp_str(0,j,i) = sgf_symbol(0,j,i)
               END DO
            END DO

            WRITE (UNIT=iw,FMT=fmtstr2)&
              "#     MO Eigenvalue [a.u.]      Occupation",&
              ((TRIM(tmp_str(0,il,im)),im=-il,il),il=0,ldos_p(ildos)%ldos%maxl)
            DO imo = 1,nmo+nvirt
              WRITE (UNIT=iw,FMT=fmtstr1) imo,eigenvalues(imo),occupation_numbers(imo),&
                (ldos_p(ildos)%ldos%pdos_array(lshell,imo),lshell=1,nsoset(ldos_p(ildos)%ldos%maxl))
            END DO
            DEALLOCATE(tmp_str,stat=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ELSE
            WRITE (UNIT=iw,FMT=fmtstr2)&
              "#     MO Eigenvalue [a.u.]      Occupation",&
              (TRIM(l_sym(il)),il=0,ldos_p(ildos)%ldos%maxl)
            DO imo = 1,nmo+nvirt
              WRITE (UNIT=iw,FMT=fmtstr1) imo,eigenvalues(imo),occupation_numbers(imo),&
                (ldos_p(ildos)%ldos%pdos_array(lshell,imo),lshell=0,ldos_p(ildos)%ldos%maxl)
            END DO
          END IF
        END IF
        CALL cp_print_key_finished_output(iw,logger,dft_section,&
                  "PRINT%PDOS", error=error)
      END IF  ! maxl>0
    END DO  ! ildos

    ! write the pdos for the lists, each ona different file,
    ! the filenames are indexed with the list number
    DO ildos=1,n_r_ldos

        npoints =  r_ldos_p(ildos)%ldos%npoints
        CALL mp_sum(npoints,para_env%group)
        CALL mp_sum(np_tot,para_env%group)
        CALL mp_sum(r_ldos_p(ildos)%ldos%pdos_array,para_env%group)
        IF(PRESENT(ispin)) THEN
          my_mittle =TRIM(spin(ispin))//"_r_list"//TRIM(r_ldos_index(ildos))
          my_spin = ispin
        ELSE
          my_mittle ="r_list"//TRIM(r_ldos_index(ildos))
          my_spin = 1
        END IF

        iw = cp_print_key_unit_nr(logger,dft_section,"PRINT%PDOS",&
             extension=".pdos", file_position=my_pos, file_action=my_act,&
             file_form="FORMATTED", middle_name=TRIM(my_mittle), &
             error=error)
        IF(iw>0) THEN
          fmtstr1 = "(I8,2X,2F16.6,  (2X,F16.8))"
          fmtstr2 = "(A42,  (10X,A8))"

          WRITE (UNIT=iw,FMT="(A,I0,A,F12.6,F12.6,A,F12.6,A)")&
            "# Projected DOS in real space, using ",npoints,&
           " points of the grid, and eval in the range",r_ldos_p(ildos)%ldos%eval_range(1:2),&
             " Hartree, E(Fermi) = ",e_fermi," a.u."
            WRITE (UNIT=iw,FMT="(A)")&
              "#     MO Eigenvalue [a.u.]      Occupation      LDOS"
            DO imo = 1,nmo+nvirt
               IF(r_ldos_p(ildos)%ldos%eval_range(1)<=eigenvalues(imo) .AND. &
                 r_ldos_p(ildos)%ldos%eval_range(2)>=eigenvalues(imo)) THEN
                 WRITE (UNIT=iw,FMT="(I8,2X,2F16.6,E20.10,E20.10)") imo,eigenvalues(imo),occupation_numbers(imo),&
                         r_ldos_p(ildos)%ldos%pdos_array(imo), r_ldos_p(ildos)%ldos%pdos_array(imo)*np_tot
               END IF
            END DO

        END IF
        CALL cp_print_key_finished_output(iw,logger,dft_section,&
                  "PRINT%PDOS", error=error)
    END DO


    ! deallocate local variables
    DEALLOCATE(pdos_array, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(firstrow, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(do_ldos) THEN
       DO ildos = 1,nldos
         DEALLOCATE(ldos_p(ildos)%ldos%pdos_array, STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         DEALLOCATE(ldos_p(ildos)%ldos%list_index, STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         DEALLOCATE(ldos_p(ildos)%ldos, STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END DO
       DEALLOCATE(ldos_p, STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(ldos_index, STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(do_r_ldos) THEN
       DO ildos = 1,n_r_ldos
         DEALLOCATE(r_ldos_p(ildos)%ldos%index_grid_local,STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         DEALLOCATE(r_ldos_p(ildos)%ldos%pdos_array,STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         DEALLOCATE(r_ldos_p(ildos)%ldos%list_index,STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         IF(.NOT. read_r(1,ildos)) THEN
           DEALLOCATE(r_ldos_p(ildos)%ldos%x_range)
         END IF
         IF(.NOT. read_r(2,ildos)) THEN
           DEALLOCATE(r_ldos_p(ildos)%ldos%y_range)
         END IF
         IF(.NOT. read_r(3,ildos)) THEN
           DEALLOCATE(r_ldos_p(ildos)%ldos%z_range)
         END IF
         IF(.NOT. read_r(4,ildos)) THEN
           DEALLOCATE(r_ldos_p(ildos)%ldos%eval_range)
         END IF
         DEALLOCATE(r_ldos_p(ildos)%ldos,STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END DO
       DEALLOCATE(read_r,STAT=stat)
       DEALLOCATE(r_ldos_p, STAT=stat)
       DEALLOCATE(r_ldos_index, STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL pw_pool_give_back_pw(auxbas_pw_pool,wf_r%pw, error=error)
       CALL pw_pool_give_back_pw(auxbas_pw_pool,wf_g%pw, error=error)
    END IF
    IF(do_virt) THEN
      DEALLOCATE(evals_virt,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      CALL cp_fm_release(mo_virt,error=error)
      CALL cp_fm_release(matrix_work,error=error)
      DEALLOCATE(eigenvalues,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(occupation_numbers,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    CALL timestop(handle)

  END SUBROUTINE  calculate_projected_dos

! *****************************************************************************
!> \brief   Compute additional virtual states  starting from the available MOS
!> \param qs_env ...
!> \param mo_set ...
!> \param evals_virt ...
!> \param mo_virt ...
!> \param nvirt ...
!> \param ispin ...
!> \param error ...
!> \date    08.03.2008
!> \par History:
!>       -
!> \par Variables
!>       -
!>       -
!> \author  MI
!> \version 1.0
! *****************************************************************************

  SUBROUTINE generate_virtual_mo (qs_env, mo_set, evals_virt, mo_virt, &
             nvirt, ispin, error )

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(mo_set_type), POINTER               :: mo_set
    REAL(KIND=dp), DIMENSION(:), POINTER     :: evals_virt
    TYPE(cp_fm_type), POINTER                :: mo_virt
    INTEGER, INTENT(IN)                      :: nvirt, ispin
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'generate_virtual_mo', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nmo, nrow_global, stat
    LOGICAL                                  :: failure
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix, s_matrix
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(preconditioner_type), POINTER       :: local_preconditioner
    TYPE(scf_control_type), POINTER          :: scf_control

    failure = .FALSE.
    NULLIFY(evals_virt,mo_virt)
    ALLOCATE(evals_virt(nvirt),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env,matrix_ks=ks_matrix, matrix_s=s_matrix, &
         scf_control=scf_control, error=error)
    CALL get_mo_set(mo_set=mo_set, mo_coeff=mo_coeff, nmo=nmo )
    CALL cp_fm_get_info(mo_coeff, context=context, para_env=para_env,&
         nrow_global=nrow_global, error=error)

    CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=context, &
         nrow_global=nrow_global,ncol_global=nvirt,error=error)
    CALL cp_fm_create(mo_virt, fm_struct_tmp,name="virtual",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)
    CALL cp_fm_init_random(mo_virt,nvirt,error=error)

    NULLIFY (local_preconditioner)

    CALL ot_eigensolver(matrix_h=ks_matrix(ispin)%matrix,matrix_s=s_matrix(1)%matrix, &
         matrix_c_fm=mo_virt, matrix_orthogonal_space_fm=mo_coeff, &
         eps_gradient=scf_control%eps_lumos, &
         preconditioner=local_preconditioner, &
         iter_max=scf_control%max_iter_lumos,&
         size_ortho_space=nmo,error=error)

    CALL calculate_subspace_eigenvalues(mo_virt,ks_matrix(ispin)%matrix,&
         evals_virt, error=error)


  END SUBROUTINE generate_virtual_mo

END MODULE qs_pdos

