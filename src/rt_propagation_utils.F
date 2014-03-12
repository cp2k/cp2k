!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines needed for EMD
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

MODULE  rt_propagation_utils



  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             cp_dbcsr_deallocate_matrix_set
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_gemm,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_transpose,&
                                             cp_fm_upper_to_full
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE input_constants,                 ONLY: use_orb_basis_set,&
                                             use_restart_wfn,&
                                             use_rt_restart
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants
  USE orbital_pointers,                ONLY: ncoset
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE qs_dftb_matrices,                ONLY: build_dftb_overlap
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_mo_methods,                   ONLY: calculate_density_matrix
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type,&
                                             read_mo_set,&
                                             read_rt_mos_from_restart
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_overlap,                      ONLY: build_overlap_matrix
  USE rt_propagation_types,            ONLY: get_rtp,&
                                             rt_prop_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC:: calc_S_derivs,&
          get_restart_wfn,&
          calculate_dp_imag,&
          calc_c_mat_force

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_propagation_utils'

!***

  CONTAINS



! *****************************************************************************
!> \brief Calculates dS/dR respectily the velocity weighted derivatves
!>        only needed for ehrenfest MD.
!>
!> \param qs_env the qs environment
!> \param error
!>
!> \par History
!>      02.2009 created [Manuel Guidon]
!> \author Florian Schiffmann
! *****************************************************************************
  SUBROUTINE calc_S_derivs(qs_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_S_derivs', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: name
    INTEGER :: handle, iatom, icol_atom, icol_global, icol_local, irow_atom, &
      irow_global, irow_local, j, m, maxder, n, natom, ncol_local, nder, &
      neighbor_list_id, nkind, nrow_local, nsgf, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, first_sgf, &
                                                last_sgf
    INTEGER, DIMENSION(6, 2)                 :: c_map_mat
    INTEGER, DIMENSION(:), POINTER           :: col_atom_index, col_indices, &
                                                row_atom_index, row_indices
    LOGICAL                                  :: failure, return_s_derivatives
    REAL(KIND=dp)                            :: scale
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_data
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: s_derivs
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: C_mat, exp_H_new
    TYPE(cp_fm_type), POINTER                :: B_mat, tmp_mat, tmp_mat2
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(rt_prop_type), POINTER              :: rtp

    failure = .FALSE.
    CALL timeset(routineN,handle)

    return_s_derivatives=.TRUE.

    NULLIFY (atomic_kind_set)
    NULLIFY (force)
    NULLIFY (particle_set)
    NULLIFY(rtp)
    NULLIFY(s_derivs)

    CALL get_qs_env(qs_env=qs_env,&
                    rtp=rtp,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    neighbor_list_id=neighbor_list_id,&
                    sab_orb=sab_orb,&
                    force=force,&
                    error=error)

    CALL get_rtp(rtp=rtp,B_mat=B_mat,C_mat=C_mat,exp_H_new=exp_H_new,error=error)
    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    ALLOCATE (atom_of_kind(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind)

    nder = 2
    maxder = ncoset(nder)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,nsgf=nsgf)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set)

    ALLOCATE (first_sgf(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (last_sgf(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL get_particle_set(particle_set=particle_set,first_sgf=first_sgf,last_sgf=last_sgf,&
         error=error)

    ! calculate the overlap derivative matrices
    IF(qs_env%dft_control%qs_control%dftb)THEN
      CALL build_dftb_overlap(qs_env,nder,s_derivs,error)
    ELSE
    CALL build_overlap_matrix(qs_env,nderivative=nder,matrix_s=s_derivs,&
         basis_set_id_a=use_orb_basis_set,basis_set_id_b=use_orb_basis_set,sab_nl=sab_orb,error=error)
    END IF

    ! here we have the S derivativs needed, lets calculate the B and C matrix
    CALL cp_fm_create(tmp_mat,&
                 matrix_struct=exp_H_new(1)%matrix%matrix_struct,&
                 name="TMP",&
                 error=error)
    CALL cp_fm_create(tmp_mat2,&
                 matrix_struct=Exp_H_new(1)%matrix%matrix_struct,&
                 name="TMP2",&
                 error=error)

    CALL cp_fm_get_info(B_mat, nrow_local=nrow_local, ncol_local=ncol_local, &
          row_indices=row_indices, col_indices=col_indices,&
          nrow_global=irow_global,ncol_global=icol_global,&
          local_data=local_data,error=error)
    local_data=0._dp
    ALLOCATE(row_atom_index(SIZE(row_indices)))
    DO irow_local=1,nrow_local
       irow_global=row_indices(irow_local)
       DO iatom=1,natom
          IF (first_sgf(iatom)<=irow_global .AND. irow_global <= last_sgf(iatom)) EXIT
       ENDDO
       row_atom_index(irow_local)=iatom
    ENDDO

    ALLOCATE(col_atom_index(SIZE(col_indices)))
    DO icol_local=1,ncol_local
       icol_global=col_indices(icol_local)
       DO iatom=1,natom
          IF (first_sgf(iatom)<=icol_global .AND. icol_global <= last_sgf(iatom)) EXIT
       ENDDO
       col_atom_index(icol_local)=iatom
    ENDDO

    !calculate scalar product v(Rb)*<alpha|d/dRb beta> (B_mat), and store the first derivatives

    CALL cp_fm_set_all(B_mat,0.0_dp,0.0_dp,error)

    DO m=1,3
       CALL cp_fm_set_all(rtp%S_der(m)%matrix,0.0_dp,0.0_dp,error)
       CALL copy_dbcsr_to_fm(s_derivs(m+1)%matrix,tmp_mat2,error=error)
       CALL cp_fm_upper_to_full(tmp_mat2,tmp_mat,error)

       DO icol_local=1,ncol_local
          DO irow_local=1,nrow_local
             irow_global=row_indices(irow_local)
             irow_atom=row_atom_index(irow_local)
             icol_global=col_indices(icol_local)
             icol_atom=col_atom_index(icol_local)
             !signs are right, since s_derivs is <d/dR alpha|beta>
             IF(irow_atom== icol_atom) CYCLE
             IF (irow_global>icol_global) THEN
                scale=1.0_dp
             ELSE
                scale=-1.0_dp
             ENDIF
             local_data(irow_local,icol_local)=local_data(irow_local,icol_local)+&
                  scale*(particle_set(icol_atom)%v(m))*tmp_mat2%local_data(irow_local,icol_local)
             rtp%S_der(m)%matrix%local_data(irow_local,icol_local)=scale*tmp_mat2%local_data(irow_local,icol_local)
          ENDDO
       ENDDO
    END DO



    !calculate C matrix: v(Rb)*<d/dRa alpha| d/dRb beta>
    c_map_mat=0
    n=0
    DO j=1,3
       DO m=j,3
          n=n+1
          c_map_mat(n,1)=j
          IF(m==j)CYCLE
          c_map_mat(n,2)=m
       END DO
    END DO

    DO m=1,3
       CALL cp_fm_set_all(C_mat(m)%matrix,0.0_dp,0.0_dp,error)
    END DO

    DO m=1,6
       CALL copy_dbcsr_to_fm(s_derivs(m+4)%matrix,tmp_mat2,error=error)
       CALL cp_fm_upper_to_full(tmp_mat2,tmp_mat,error)
       DO j=1,2
          IF(c_map_mat(m,j)==0)CYCLE
          CALL cp_fm_scale_and_add(1.0_dp,C_mat(c_map_mat(m,j))%matrix,-1.0_dp,tmp_mat2,error)
       END DO
    END DO
    DO m=1,3
       NULLIFY(local_data)
       CALL cp_fm_get_info(C_mat(m)%matrix,local_data=local_data,error=error)
       DO icol_local=1,ncol_local
          DO irow_local=1,nrow_local
             irow_global=row_indices(irow_local)
             irow_atom=row_atom_index(irow_local)
             icol_global=col_indices(icol_local)
             icol_atom=col_atom_index(icol_local)

             local_data(irow_local,icol_local)=local_data(irow_local,icol_local)*&
                  particle_set(irow_atom)%v(m)
          ENDDO
       ENDDO
    END DO

    CALL cp_fm_release(tmp_mat,error)
    CALL cp_fm_release(tmp_mat2,error)
    CALL cp_dbcsr_deallocate_matrix_set(s_derivs,error)
    DEALLOCATE(row_atom_index,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(col_atom_index,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (first_sgf,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (last_sgf,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE calc_S_derivs


! *****************************************************************************
!> \brief reads the restart file. At the moment only SCF (means only real)
!> \author Florian Schiffmann (02.09)
! *****************************************************************************


  SUBROUTINE get_restart_wfn(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_restart_wfn', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, id_nr, ispin, nspin, &
                                                qs_env_id
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: p_rmpv
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_new, mos_old
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(section_vals_type), POINTER         :: dft_section, input

    NULLIFY(atomic_kind_set,mo_array,particle_set)

    CALL get_qs_env(qs_env,atomic_kind_set=atomic_kind_set,&
         particle_set=particle_set,mos=mo_array, &
         input=input,rtp=rtp,&
         id_nr=qs_env_id, dft_control=dft_control,&
         error=error)

    id_nr=0
    nspin=SIZE(mo_array)
    p_rmpv => qs_env%rho%rho_ao
    dft_section =>  section_vals_get_subs_vals(input,"DFT",error=error)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set)
    CALL get_rtp(rtp=rtp,mos_old=mos_old,mos_new=mos_new,error=error)

    SELECT CASE(dft_control%rtp_control%initial_wfn)
    CASE(use_restart_wfn)
       CALL read_mo_set(mo_array,atomic_kind_set,particle_set,qs_env%para_env,&
            id_nr=id_nr,multiplicity=dft_control%multiplicity,dft_section=dft_section,&
            error=error)
       DO i=1,SIZE(qs_env%mos)
          CALL cp_fm_to_fm(qs_env%mos(i)%mo_set%mo_coeff,mos_old(2*i-1)%matrix,error)
          CALL cp_fm_set_all(mos_new(2*i)%matrix,zero,zero,error)
       END DO
    CASE(use_rt_restart)
       CALL read_rt_mos_from_restart(mo_array,mos_old,atomic_kind_set,particle_set,qs_env%para_env,&
            id_nr,dft_control%multiplicity,dft_section, error)
    END SELECT

! Can't quite see how things end up on mo array in case of rt_restart???
    DO ispin=1,nspin
       CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
            p_rmpv(ispin)%matrix,error=error)
    ENDDO


  END SUBROUTINE get_restart_wfn

! *****************************************************************************
!> \brief calculates the imaginary part of the density matrix
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE calculate_dp_imag(dp_imag,mos_new,mo_set,error)
    TYPE(cp_fm_type), POINTER                :: dp_imag
    TYPE(cp_fm_p_type), DIMENSION(2)         :: mos_new
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_dp_imag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nao, nmo
    REAL(KIND=dp)                            :: alpha
    TYPE(cp_fm_type), POINTER                :: fm_tmp

    CALL timeset(routineN,handle)

    IF ( .NOT. mo_set%uniform_occupation ) THEN ! not all orbitals 1..homo are equally occupied
       NULLIFY(fm_tmp)
       CALL cp_fm_create(fm_tmp,mo_set%mo_coeff%matrix_struct,error=error)
       CALL cp_fm_to_fm(mos_new(2)%matrix,fm_tmp,error=error)
       CALL cp_fm_column_scale(fm_tmp,mo_set%occupation_numbers(1:mo_set%homo))
       alpha=1.0_dp

    ELSE
       alpha=mo_set%maxocc
       CALL cp_fm_create(fm_tmp,mo_set%mo_coeff%matrix_struct,error=error)
       CALL cp_fm_to_fm(mos_new(2)%matrix,fm_tmp,error=error)
    ENDIF
    CALL get_mo_set(mo_set, nao=nao,nmo=nmo)

    CALL cp_fm_gemm("N","T",nao,nao,nmo,alpha,mos_new(1)%matrix,fm_tmp,0.0_dp,dp_imag,error)

    CALL cp_fm_gemm("N","T",nao,nao,nmo,-alpha,fm_tmp,mos_new(1)%matrix,1.0_dp,dp_imag,error)

    CALL cp_fm_release(fm_tmp,error=error)

    CALL timestop(handle)

  END SUBROUTINE calculate_dp_imag

! *****************************************************************************
!> \brief calculates the three additional force contributions needed in EMD
!>        P_imag*C , P_imag*B*S^-1*S_der , P*S^-1*H*S_der
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE calc_c_mat_force(qs_env,dp_imag,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_type), POINTER                :: dp_imag
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_c_mat_force', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, iatom, icol_atom, icol_global, icol_local, ikind, &
      im, irow_atom, irow_global, irow_local, ispin, kind_atom, nao, natom, &
      ncol_local, nmo, nrow_local, re, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, col_atom_index, &
                                                first_sgf, kind_of, last_sgf, &
                                                row_atom_index
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha, cforce(2,3)
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_data
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: C_mat, mos_new, S_der, SinvB, &
                                                SinvH
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: B_mat, S_inv, tmp_mat, &
                                                tmp_mat1, tmp_mat2, tmp_mat3, &
                                                tmp_mo1, tmp_mo2
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(rt_prop_type), POINTER              :: rtp

    failure=.FALSE.

    CALL timeset(routineN,handle)

    NULLIFY(rtp,particle_set,atomic_kind_set,mos)
    CALL get_qs_env(qs_env=qs_env,rtp=rtp,particle_set=particle_set,&
                 atomic_kind_set=atomic_kind_set,mos=mos,force=force,error=error)

    CALL get_rtp(rtp=rtp,C_mat=C_mat,mos_new=mos_new,error=error)

    natom = SIZE(particle_set)
    ALLOCATE (first_sgf(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (last_sgf(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL get_particle_set(particle_set=particle_set,first_sgf=first_sgf,last_sgf=last_sgf,&
         error=error)

    CALL cp_fm_get_info(C_mat(1)%matrix, nrow_local=nrow_local, ncol_local=ncol_local, &
          row_indices=row_indices, col_indices=col_indices, &
          error=error)

    CALL cp_fm_create(tmp_mat,C_mat(1)%matrix%matrix_struct,error=error)
    CALL cp_fm_create(tmp_mat1,C_mat(1)%matrix%matrix_struct,error=error)
    CALL cp_fm_create(tmp_mat2,C_mat(1)%matrix%matrix_struct,error=error)
    CALL cp_fm_create(tmp_mat3,C_mat(1)%matrix%matrix_struct,error=error)

    ALLOCATE (atom_of_kind(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (kind_of(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(row_atom_index(SIZE(row_indices)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO irow_local=1,nrow_local
       irow_global=row_indices(irow_local)
       DO iatom=1,natom
          IF (first_sgf(iatom)<=irow_global .AND. irow_global <= last_sgf(iatom)) EXIT
       ENDDO
       row_atom_index(irow_local)=iatom
    ENDDO

    ALLOCATE(col_atom_index(SIZE(col_indices)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO icol_local=1,ncol_local
       icol_global=col_indices(icol_local)
       DO iatom=1,natom
          IF (first_sgf(iatom)<=icol_global .AND. icol_global <= last_sgf(iatom)) EXIT
       ENDDO
       col_atom_index(icol_local)=iatom
    ENDDO

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind,kind_of=kind_of)


    CALL cp_fm_create(dp_imag,&
         matrix_struct=qs_env%rtp%B_mat%matrix_struct,&
         name="dp_imag",&
         error=error)


    DO ispin=1,SIZE(mos)
       re=2*ispin-1
       im=2*ispin
       CALL calculate_dp_imag(dp_imag,mos_new(re:im),mos(ispin)%mo_set,error)

       alpha=mos(ispin)%mo_set%maxocc

       ! calculate the W matrix analogue for EMD
       ! calculate B_mat force

       CALL get_rtp(rtp=rtp,B_mat=B_mat,S_der=S_der,S_inv=S_inv,SinvH=SinvH,SinvB=SinvB,mos_new=mos_new,error=error)
       CALL cp_fm_get_info(mos_new(re)%matrix,nrow_global=nao,ncol_global=nmo,error=error)

       CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nmo,ncol_global=nao,&
            para_env=mos_new(re)%matrix%matrix_struct%para_env, &
            context=mos_new(re)%matrix%matrix_struct%context,error=error)

       CALL cp_fm_create(tmp_mo1,tmp_fm_struct,error=error)
       CALL cp_fm_create(tmp_mo2,tmp_fm_struct,error=error)
       CALL cp_fm_struct_release(tmp_fm_struct,error)
       CALL cp_fm_gemm("T","T",nmo,nao,nao,alpha,mos_new(im)%matrix,SinvB(ispin)%matrix,0.0_dp,tmp_mo1,error=error)
       CALL cp_fm_gemm("T","T",nmo,nao,nao,alpha,mos_new(re)%matrix,SinvH(ispin)%matrix,1.0_dp,tmp_mo1,error=error)
       CALL cp_fm_gemm("T","T",nmo,nao,nao,-alpha,mos_new(re)%matrix,SinvB(ispin)%matrix,0.0_dp,tmp_mo2,error=error)
       CALL cp_fm_gemm("T","T",nmo,nao,nao,alpha,mos_new(im)%matrix,SinvH(ispin)%matrix,1.0_dp,tmp_mo2,error=error)
       CALL cp_fm_gemm("N","N",nao,nao,nmo,1.0_dp,mos_new(re)%matrix,tmp_mo1,0.0_dp,tmp_mat1,error=error)
       CALL cp_fm_gemm("N","N",nao,nao,nmo,1.0_dp,mos_new(im)%matrix,tmp_mo2,1.0_dp,tmp_mat1,error=error)


       CALL cp_fm_transpose(tmp_mat1,tmp_mat,error)

       !C_mat forces

       CALL cp_fm_transpose(dp_imag,tmp_mat2,error)
       CALL cp_fm_get_info(tmp_mat2,local_data=local_data,error=error)

       !This loop evaluates all forces, and puts them in the right place
       DO i=1,3
          CALL cp_fm_transpose(C_mat(i)%matrix,tmp_mat3,error)

          DO icol_local=1,ncol_local
             DO irow_local=1,nrow_local

                icol_atom=col_atom_index(icol_local)
                irow_atom=row_atom_index(irow_local)
                ikind=kind_of(icol_atom)
                kind_atom=atom_of_kind(icol_atom)

                force(ikind)%ehrenfest(i,kind_atom)=force(ikind)%ehrenfest(i,kind_atom)+&
                     tmp_mat%local_data(irow_local,icol_local)* S_der(i)%matrix%local_data(irow_local,icol_local)-&
                     local_data(irow_local,icol_local)* C_mat(i)%matrix%local_data(irow_local,icol_local)
                ikind=kind_of(irow_atom)
                kind_atom=atom_of_kind(irow_atom)

                force(ikind)%ehrenfest(i,kind_atom)=force(ikind)%ehrenfest(i,kind_atom)-&
                     tmp_mat1%local_data(irow_local,icol_local)* S_der(i)%matrix%local_data(irow_local,icol_local)+ &
                     local_data(irow_local,icol_local)* tmp_mat3%local_data(irow_local,icol_local)
                cforce(kind_atom,i)=cforce(kind_atom,i)+local_data(irow_local,icol_local)* &
                                     tmp_mat3%local_data(irow_local,icol_local)
             END DO
          END DO
       END DO

       ! recall QS forces, at this point have the other sign.
       CALL cp_fm_release(tmp_mo1,error)
       CALL cp_fm_release(tmp_mo2,error)

    END DO
    ! recall QS forces, at this point have the other sign.
    DO i=1,SIZE(force)
       force(i)%ehrenfest(:,:)=- force(i)%ehrenfest(:,:)
    END DO

    CALL cp_fm_release(tmp_mat,error)
    CALL cp_fm_release(tmp_mat1,error)
    CALL cp_fm_release(tmp_mat2,error)
    CALL cp_fm_release(tmp_mat3,error)
    CALL cp_fm_release(dp_imag,error)

    DEALLOCATE(row_atom_index,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(col_atom_index,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (first_sgf,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (last_sgf,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (atom_of_kind,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (kind_of,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle)
  END SUBROUTINE calc_c_mat_force

END MODULE rt_propagation_utils
