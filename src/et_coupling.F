!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief calculates the electron transfer coupling elements
!>      Wu, Van Voorhis, JCP 125, 164105 (2006)
!> \author fschiff (01.2007)
! *****************************************************************************
MODULE et_coupling

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_deallocate_matrix_set,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_invert,&
                                             cp_fm_transpose
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             fm_pool_get_el_struct
  USE cp_fm_struct,                    ONLY: cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_subsys_types,                 ONLY: cp_subsys_get
  USE f77_blas
  USE input_constants,                 ONLY: do_spin_density
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE mathlib,                         ONLY: diamat_all
  USE message_passing,                 ONLY: mp_sum
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE pw_types,                        ONLY: pw_p_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_matrix_pools,                 ONLY: mpools_get
  USE qs_mo_types,                     ONLY: get_mo_set
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'et_coupling'

! *** Public subroutines ***

  PUBLIC :: calc_et_coupling,becke_restraint,becke_force

CONTAINS
! *****************************************************************************
  SUBROUTINE calc_et_coupling(qs_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calc_et_coupling', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iw, j, k, my_id, &
                                                nao, ncol_local, nmo, &
                                                nrow_local
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    REAL(KIND=dp)                            :: Sda, strength, Waa, Wbb, Wda
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: a, b, S_det
    REAL(KIND=dp), DIMENSION(2)              :: eigenv
    REAL(KIND=dp), DIMENSION(2, 2)           :: S_mat, tmp_mat, U, W_mat
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: rest_MO
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: mo_mo_fm_pools
    TYPE(cp_fm_struct_type), POINTER         :: mo_mo_fmstruct
    TYPE(cp_fm_type), POINTER                :: inverse_mat, SMO, Tinverse, &
                                                tmp2
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: et_coupling_section

    NULLIFY(tmp2,SMO,rest_MO,mo_mo_fmstruct, rest_MO)

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    et_coupling_section => section_vals_get_subs_vals(qs_env%input,&
                            "PROPERTIES%ET_COUPLING",error=error)
    my_id=qs_env%dft_control%qs_control%becke_control%density_type

    iw=cp_print_key_unit_nr(logger,et_coupling_section,"PROGRAM_RUN_INFO",&
                                         extension=".log",error=error)

    ALLOCATE(rest_MO(2))
    ALLOCATE(a(qs_env%dft_control%nspins))
    ALLOCATE(b(qs_env%dft_control%nspins))
    ALLOCATE(S_det(qs_env%dft_control%nspins))

    CALL mpools_get(qs_env%mpools, mo_mo_fm_pools=mo_mo_fm_pools,&
         error=error)
    mo_mo_fmstruct => fm_pool_get_el_struct(mo_mo_fm_pools(1)%pool,&
         error=error)
    DO i=1 , qs_env%dft_control%nspins
    mo_mo_fmstruct => fm_pool_get_el_struct(mo_mo_fm_pools(i)%pool,&
         error=error)

    CALL get_mo_set(mo_set=qs_env%mos(i)%mo_set,&
                    nao=nao,&
                    nmo=nmo)

    CALL cp_fm_create(matrix=tmp2,&
         matrix_struct=qs_env%mos(i)%mo_set%mo_coeff%matrix_struct,&
         name="ET_TMP"//TRIM(ADJUSTL(cp_to_string(2)))//"MATRIX",&
         error=error)
    CALL cp_fm_create(matrix=inverse_mat,&
         matrix_struct=mo_mo_fmstruct,&
         name="INVERSE"//TRIM(ADJUSTL(cp_to_string(2)))//"MATRIX",&
         error=error)
    CALL cp_fm_create(matrix=Tinverse,&
         matrix_struct=mo_mo_fmstruct,&
         name="T_INVERSE"//TRIM(ADJUSTL(cp_to_string(2)))//"MATRIX",&
         error=error)
    CALL cp_fm_create(matrix=SMO,&
         matrix_struct=mo_mo_fmstruct,&
         name="ET_SMO"//TRIM(ADJUSTL(cp_to_string(1)))//"MATRIX",&
         error=error)
    DO j=1,2
       NULLIFY(rest_MO(j)%matrix)
       CALL cp_fm_create(matrix=rest_MO(j)%matrix,&
            matrix_struct=mo_mo_fmstruct,&
            name="ET_rest_MO"//TRIM(ADJUSTL(cp_to_string(j)))//"MATRIX",&
            error=error)
    END DO

!   calculate MO-overlap

    CALL cp_dbcsr_sm_fm_multiply(qs_env%matrix_s(1)%matrix,qs_env%et_coupling%et_mo_coeff(i)%matrix,&
         tmp2,nmo,1.0_dp,0.0_dp,error=error)
    CALL cp_fm_gemm('T','N',nmo,nmo,nao,1.0_dp,&
         qs_env%mos(i)%mo_set%mo_coeff,&
         tmp2,0.0_dp,SMO,error=error)

!    calculate the MO-representation of the restraint matrix A

    CALL cp_dbcsr_sm_fm_multiply(qs_env%et_coupling%rest_mat(1)%matrix,&
         qs_env%et_coupling%et_mo_coeff(i)%matrix,&
         tmp2,nmo,1.0_dp,0.0_dp,error=error)

    CALL cp_fm_gemm('T','N',nmo,nmo,nao,1.0_dp,&
         qs_env%mos(i)%mo_set%mo_coeff,&
         tmp2,0.0_dp,rest_MO(1)%matrix,error=error)

!    calculate the MO-representation of the restraint matrix D

    CALL cp_dbcsr_sm_fm_multiply(qs_env%et_coupling%rest_mat(2)%matrix,&
         qs_env%mos(i)%mo_set%mo_coeff,&
         tmp2,nmo,1.0_dp,0.0_dp,error=error)

    CALL cp_fm_gemm('T','N',nmo,nmo,nao,1.0_dp,&
         qs_env%et_coupling%et_mo_coeff(i)%matrix,&
         tmp2,0.0_dp,rest_MO(2)%matrix,error=error)

   CALL cp_fm_invert(SMO,inverse_mat,S_det(i),error=error)

    CALL cp_fm_get_info(inverse_mat, nrow_local=nrow_local, ncol_local=ncol_local,&
                       row_indices=row_indices,col_indices=col_indices,error=error)
    b(i)=0.0_dp

    DO j=1,ncol_local
       DO k=1,nrow_local
          b(i)=b(i)+rest_MO(2)%matrix%local_data(k,j)*inverse_mat%local_data(k,j)
       END DO
    END DO

    CALL cp_fm_transpose(inverse_mat,Tinverse,error)
    a(i)=0.0_dp
    DO j=1,ncol_local
       DO k=1,nrow_local
          a(i)=a(i)+rest_MO(1)%matrix%local_data(k,j)*Tinverse%local_data(k,j)
       END DO
    END DO
    IF((my_id==do_spin_density).AND.i==2)THEN
       a(i)=-a(i)
       b(i)=-b(i)
    END IF
    CALL mp_sum(a(i),qs_env%para_env%group)

    CALL mp_sum(b(i),qs_env%para_env%group)

       CALL cp_fm_release(tmp2,error)
       CALL cp_fm_release(rest_MO(1)%matrix,error)
       CALL cp_fm_release(rest_MO(2)%matrix,error)
       CALL cp_fm_release(SMO,error)
       CALL cp_fm_release(Tinverse,error)
       CALL cp_fm_release(inverse_mat,error)
    END DO
    DEALLOCATE(rest_MO)

!    solve eigenstates for the projector matrix

    IF( qs_env%dft_control%nspins ==2)THEN
       Sda=S_det(1)*S_det(2)
       Wda=((a(1)+a(2))+(b(1)+b(2)))*0.5_dp*Sda
    ELSE
       Sda=S_det(1)**2
       Wda=(a(1)+b(1))*Sda
    END IF

    IF(qs_env%dft_control%qs_control%ddapc_restraint)THEN
    Waa=qs_env%et_coupling%order_p
    Wbb=qs_env%dft_control%qs_control%ddapc_restraint_control(1)%ddapc_restraint_control%ddapc_order_p
    strength=qs_env%dft_control%qs_control%ddapc_restraint_control(1)%ddapc_restraint_control%strength
    END IF
    IF(qs_env%dft_control%qs_control%becke_restraint)THEN
    Waa=qs_env%et_coupling%order_p
    Wbb=qs_env%dft_control%qs_control%becke_control%becke_order_p
    strength=qs_env%dft_control%qs_control%becke_control%strength
    END IF

!!   construct S and W   !!!
    S_mat(1,1)=1.0_dp
    S_mat(2,2)=1.0_dp
    S_mat(2,1)=Sda
    S_mat(1,2)=Sda

    W_mat(1,1)=Wbb
    W_mat(2,2)=Waa
    W_mat(2,1)=Wda
    W_mat(1,2)=Wda

!!  solve WC=SCN
    CALL diamat_all (S_mat,eigenv,.TRUE.,error=error)
    ! U = S**(-1/2)
    U=0.0_dp
    U(1,1)=1.0_dp/SQRT(eigenv(1))
    U(2,2)=1.0_dp/SQRT(eigenv(2))
    tmp_mat=MATMUL(U,TRANSPOSE(S_mat))
    U=MATMUL(S_mat,tmp_mat)
    tmp_mat=MATMUL(W_mat,U)
    W_mat=MATMUL(U,tmp_mat)
    CALL diamat_all (W_mat,eigenv,.TRUE.,error=error)
    tmp_mat=MATMUL(U,W_mat)

    W_mat(1,1)=qs_env%energy%total
    W_mat(2,2)=qs_env%et_coupling%energy
    a(1)=(qs_env%energy%total+strength*Wbb)*Sda-strength*Wda
    a(2)=(qs_env%et_coupling%energy+qs_env%et_coupling%e1*Waa)*Sda-qs_env%et_coupling%e1*Wda
    W_mat(1,2)=(a(1)+a(2))*0.5_dp
    W_mat(2,1)=W_mat(1,2)

    S_mat=MATMUL(W_mat,(tmp_mat))
    W_mat=MATMUL(TRANSPOSE(tmp_mat),S_mat)

    IF (iw>0) THEN
       WRITE(iw,*)
       WRITE(iw,'(T3,A,T60,(3X,F12.6))') 'Strength of constraint A          :',qs_env%et_coupling%e1
       WRITE(iw,'(T3,A,T60,(3X,F12.6))') 'Strength of constraint B          :',strength
       WRITE(iw,'(T3,A,T60,(3X,F12.6))') 'Final target value of constraint A:',Waa
       WRITE(iw,'(T3,A,T60,(3X,F12.6))') 'Final target value of constraint B:',Wbb
       WRITE(iw,*)
       WRITE(iw,'(T3,A,T60,(3X,F12.6))')&
            'Diabatic electronic coupling matrix element(mHartree):',ABS(W_mat(1,2)*1000.0_dp)

    ENDIF

    CALL cp_dbcsr_deallocate_matrix_set(qs_env%et_coupling%rest_mat,error=error)

    CALL cp_print_key_finished_output(iw,logger,et_coupling_section,&
               "PROGRAM_RUN_INFO", error=error)
    CALL timestop(handle)
  END SUBROUTINE calc_et_coupling

! *****************************************************************************
!> \brief calculates a becke contraint
!> \param calc_pot logical if potential has to be calculated or only_energy
!> \author fschiff (01.2007)
! *****************************************************************************
  SUBROUTINE becke_restraint(qs_env,becke_const ,calc_pot,calculate_forces,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(pw_p_type)                          :: becke_const
    LOGICAL                                  :: calc_pot, calculate_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'becke_restraint', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iatom, ip, istat, &
                                                j, jatom, k, np(3)
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: catom
    INTEGER, DIMENSION(2, 3)                 :: bo
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: dE, dist1, dist2, dvol, my1, &
                                                myexp, strength, &
                                                sum_cell_f_all, &
                                                sum_cell_f_constr, target_val
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: cell_functions
    REAL(kind=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: R12
    REAL(kind=dp), DIMENSION(3)              :: cell_v, dist_vec, dr, grid_p, &
                                                r, r1, shift
    REAL(kind=dp), DIMENSION(:), POINTER     :: coeff
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_list_type), POINTER        :: particles

    failure=.FALSE.
      CALL timeset(routineN,handle)

    IF (qs_env%dft_control%qs_control%becke_restraint) THEN

          target_val=qs_env%dft_control%qs_control%becke_control%target
          strength=qs_env%dft_control%qs_control%becke_control%strength
          dvol = becke_const%pw%pw_grid%dvol
          IF(calc_pot)THEN

          NULLIFY(cell,particles,coeff)
          CALL cp_subsys_get(qs_env%subsys,particles=particles,error=error)
          becke_const%pw%cr3d=0.0_dp

          CALL get_qs_env(qs_env,cell=cell,error=error)

          ALLOCATE (catom(qs_env%dft_control%qs_control%becke_control%natoms),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                           "catom",int_size*&
                                           qs_env%dft_control%qs_control%becke_control%natoms)
          DO i=1,qs_env%dft_control%qs_control%becke_control%natoms
             catom(i)=qs_env%dft_control%qs_control%becke_control%atoms(i)
          ENDDO

          bo=becke_const%pw%pw_grid%bounds_local
          np=becke_const%pw%pw_grid%npts
          dvol = becke_const%pw%pw_grid%dvol
          dr=becke_const%pw%pw_grid%dr
          np= becke_const%pw%pw_grid%npts
          shift=-REAL(MODULO(np,2),dp)*dr/2.0_dp
          coeff => qs_env%dft_control%qs_control%becke_control%coeff
          ALLOCATE (cell_functions(particles%n_els),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                           "cell_functions",dp_size*particles%n_els)
          !calculate distances from target atom (only implemented for a diatomic system)

          DO i=1,3
             cell_v(i)=cell%hmat(i,i)
          END DO

          ALLOCATE (R12(particles%n_els,particles%n_els),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                           "R12",dp_size*particles%n_els**2)

          DO iatom=1,particles%n_els-1
             DO jatom=iatom+1,particles%n_els
                r=particles%els(iatom)%r
                r1=particles%els(jatom)%r
                DO i=1,3
                   r(i)=MODULO(r(i),cell%hmat(i,i))-cell%hmat(i,i)/2._dp
                   r1(i)=MODULO(r1(i),cell%hmat(i,i))-cell%hmat(i,i)/2._dp
                END DO
                dist_vec=(r-r1)-ANINT((r-r1)/cell_v)*cell_v
                R12(iatom,jatom)=SQRT(DOT_PRODUCT(dist_vec,dist_vec))
                R12(jatom,iatom)=R12(iatom,jatom)
             END DO
          END DO

          DO k=bo(1,1),bo(2,1)
             DO j=bo(1,2),bo(2,2)
                DO i=bo(1,3),bo(2,3)
                   grid_p(1)=k*dr(1)+shift(1)
                   grid_p(2)=j*dr(2)+shift(2)
                   grid_p(3)=i*dr(3)+shift(3)

                   cell_functions=1.0_dp
                   DO iatom=1,particles%n_els
                      r=particles%els(iatom)%r
                      DO ip=1,3
                         r(ip)=MODULO(r(ip),cell%hmat(ip,ip))-cell%hmat(ip,ip)/2._dp
                      END DO
                      dist_vec=(r-grid_p)-ANINT((r-grid_p)/cell_v)*cell_v
                      dist1=SQRT(DOT_PRODUCT(dist_vec,dist_vec))
                      IF(dist1.LE.6.0_dp)THEN
                      DO jatom=1,particles%n_els
                         IF(jatom.NE.iatom)THEN
                         r1=particles%els(jatom)%r
                         DO ip=1,3
                            r1(ip)=MODULO(r1(ip),cell%hmat(ip,ip))-cell%hmat(ip,ip)/2._dp
                         END DO
                         dist_vec=(r1-grid_p)-ANINT((r1-grid_p)/cell_v)*cell_v
                         dist2=SQRT(DOT_PRODUCT(dist_vec,dist_vec))
                         my1= (dist1-dist2)/R12(iatom,jatom)
                         myexp=1.5_dp*my1-0.5_dp*my1**3
                         myexp=1.5_dp*myexp-0.5_dp*myexp**3
                         myexp=1.5_dp*myexp-0.5_dp*myexp**3
                         cell_functions(iatom)=cell_functions(iatom)*0.5_dp*(1-myexp)
                      END IF
                      END DO
                      ELSE
                         cell_functions(iatom)=0.0_dp
                      END IF
                   END DO

                   sum_cell_f_constr=0.0_dp
                   DO ip=1,qs_env%dft_control%qs_control%becke_control%natoms
                      sum_cell_f_constr=sum_cell_f_constr+cell_functions(catom(ip))*coeff(ip)
                   END DO

                   sum_cell_f_all=0.0_dp
                   DO ip=1,particles%n_els
                      sum_cell_f_all=sum_cell_f_all+cell_functions(ip)
                   END DO
                   IF(ABS(sum_cell_f_all).LE.0.000001)THEN
                      becke_const%pw%cr3d(k,j,i)=0.0_dp
                   ELSE
                      becke_const%pw%cr3d(k,j,i)=sum_cell_f_constr/sum_cell_f_all
                   END IF

                END DO
             END DO
          END DO

          DEALLOCATE(catom,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"catom")
          DEALLOCATE(R12,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"R12")
          DEALLOCATE(cell_functions,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"cell_functions")

       END IF

       dE=0.0_dp
       DO i=1 , qs_env%dft_control%nspins
             dE=dE+accurate_sum(becke_const%pw%cr3d*qs_env%rho%rho_r(i)%pw%cr3d)*dvol
       END DO
       CALL mp_sum(dE,qs_env%para_env%group)
       qs_env%dft_control%qs_control%becke_control%becke_order_p=dE
       qs_env%energy%becke=(dE-target_val)*strength

       IF(calculate_forces) CALL becke_force(qs_env,becke_const,error)

    END IF
    CALL timestop(handle)
  END SUBROUTINE becke_restraint

! *****************************************************************************
!> \brief calculates a becke contraint forces
!> \author fschiff (01.2007)
! *****************************************************************************
  SUBROUTINE becke_force(qs_env,becke_const,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(pw_p_type)                          :: becke_const
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'becke_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iatom, ikind, ip, &
                                                ispin, istat, j, jatom, k, &
                                                natom, np(3)
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, catom, kind_of
    INTEGER, DIMENSION(2, 3)                 :: bo
    LOGICAL                                  :: failure
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: is_constraint
    REAL(kind=dp) :: dist1, dist2, dpolyn, dvol, f_3, my1, polyn, strength, &
      sum_cell_f_constr, sum_Pm, target_val, th, tmp_const
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: coefficients, P_i
    REAL(kind=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: d_sum_const_dR, d_sum_Pm_dR, &
                                                dP_i_dRi, dw_dR, integrated, &
                                                R12, s_my_ij
    REAL(kind=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: dP_i_dRj, ds_dR
    REAL(kind=dp), DIMENSION(3) :: cell_v, dist_vec_i, dist_vec_j, dmy_dR_i, &
      dmy_dR_j, dr, dr1_r2, dr_i_dR, dr_ij_dR, dr_j_dR, grid_p, r, r1, shift
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_rho_type), POINTER               :: rho

    failure=.FALSE.
    CALL timeset(routineN,handle)
    th=1.0e-8_dp
    IF(qs_env%dft_control%qs_control%becke_restraint)THEN

       target_val=qs_env%dft_control%qs_control%becke_control%TARGET
       strength=qs_env%dft_control%qs_control%becke_control%strength

       NULLIFY(cell,particles)
       CALL cp_subsys_get(qs_env%subsys,particles=particles,error=error)

       CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set,&
                       particle_set=particle_set,&
                       cell=cell,rho=rho,force=force,error=error)

       ALLOCATE (catom(qs_env%dft_control%qs_control%becke_control%natoms),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                        "catom",dp_size*qs_env%dft_control%qs_control%becke_control%natoms)
       ALLOCATE (is_constraint(particles%n_els),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                        "is_constraint",dp_size*particles%n_els)
       ALLOCATE (coefficients(particles%n_els),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                        "coefficients",dp_size*particles%n_els)

       is_constraint=.FALSE.
       DO i=1,qs_env%dft_control%qs_control%becke_control%natoms
          catom(i)=qs_env%dft_control%qs_control%becke_control%atoms(i)
          is_constraint(qs_env%dft_control%qs_control%becke_control%atoms(i))=.TRUE.
          coefficients(catom(i))= qs_env%dft_control%qs_control%becke_control%coeff(i)
       ENDDO

       natom = SIZE(particle_set)
       ALLOCATE (atom_of_kind(natom),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                        "atom_of_kind",natom*int_size)

       ALLOCATE (kind_of(natom),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                        "kind_of",natom*int_size)

       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                                atom_of_kind=atom_of_kind,&
                                kind_of=kind_of)

       bo=becke_const%pw%pw_grid%bounds_local
       np=becke_const%pw%pw_grid%npts
       dvol = becke_const%pw%pw_grid%dvol
       dr=becke_const%pw%pw_grid%dr
       np= becke_const%pw%pw_grid%npts
       shift=-REAL(MODULO(np,2),dp)*dr/2.0_dp
       !calculate distances from target atom (only implemented for a diatomic system)

       DO i=1,3
          cell_v(i)=cell%hmat(i,i)
       END DO

       ALLOCATE (R12(particles%n_els,particles%n_els),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                        "R12",dp_size*particles%n_els**2)
       ALLOCATE (s_my_ij(particles%n_els,particles%n_els),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                        "s_my_ij",dp_size*particles%n_els**2)
       ALLOCATE (ds_dR(3,particles%n_els,particles%n_els),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                        "ds_dR",dp_size*3*particles%n_els**2)
       ALLOCATE (P_i(particles%n_els),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                        "P_i",dp_size*particles%n_els)
       ALLOCATE (dw_dR(3,particles%n_els),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                        "dw_dR",dp_size*3*particles%n_els)
       ALLOCATE (integrated(3,particles%n_els),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                        "integrated",dp_size*3*particles%n_els)
       ALLOCATE (d_sum_Pm_dR(3,particles%n_els),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                        "d_sum_Pm_dR",dp_size*3*particles%n_els)
       ALLOCATE (d_sum_const_dR(3,particles%n_els),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                        "d_sum_const_dR",dp_size*3*particles%n_els)
       ALLOCATE (dP_i_dRj(3,particles%n_els,particles%n_els),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                        "dP_i_dRj",dp_size*3*particles%n_els**2)
       ALLOCATE (dP_i_dRi(3,particles%n_els),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                        "dP_i_dRi",dp_size*3*particles%n_els)

       DO iatom=1,particles%n_els-1
          DO jatom=iatom+1,particles%n_els
             r=particles%els(iatom)%r
             r1=particles%els(jatom)%r
             DO i=1,3
                r(i)=MODULO(r(i),cell%hmat(i,i))-cell%hmat(i,i)/2._dp
                r1(i)=MODULO(r1(i),cell%hmat(i,i))-cell%hmat(i,i)/2._dp
             END DO
             dist_vec_i=(r-r1)-ANINT((r-r1)/cell_v)*cell_v
             R12(iatom,jatom)=SQRT(DOT_PRODUCT(dist_vec_i,dist_vec_i))
             R12(jatom,iatom)=R12(iatom,jatom)
          END DO
       END DO

       integrated=0.0_dp

       DO k=bo(1,1),bo(2,1)
          DO j=bo(1,2),bo(2,2)
             DO i=bo(1,3),bo(2,3)
                grid_p(1)=k*dr(1)+shift(1)
                grid_p(2)=j*dr(2)+shift(2)
                grid_p(3)=i*dr(3)+shift(3)

                d_sum_Pm_dR=0.0_dp
                d_sum_const_dR=0.0_dp
                P_i=1.0_dp
                dP_i_dRi=0.0_dp
                DO iatom=1,particles%n_els
                   r=particles%els(iatom)%r
                   DO ip=1,3
                      r(ip)=MODULO(r(ip),cell%hmat(ip,ip))-cell%hmat(ip,ip)/2._dp
                   END DO
                   dist_vec_i=(r-grid_p)-ANINT((r-grid_p)/cell_v)*cell_v
                   dist1=SQRT(DOT_PRODUCT(dist_vec_i,dist_vec_i))
                   IF (dist1.LE.5.0_dp)THEN
                   IF(dist1.le.th)dist1=th
                   dr_i_dR(:)= dist_vec_i(:) / dist1
!                   IF(dist1.le.0.00001)dr_i_dR=0.0_dp
                   DO jatom=1,particles%n_els
                      IF (jatom.NE.iatom)THEN
                         r1=particles%els(jatom)%r
                         DO ip=1,3
                            r1(ip)=MODULO(r1(ip),cell%hmat(ip,ip))-cell%hmat(ip,ip)/2._dp
                         END DO
                         dist_vec_j= (r1-grid_p) -ANINT((r1-grid_p)/cell_v) *cell_v
                         dr1_r2=(r-r1) -ANINT((r-r1)/cell_v) *cell_v
                         dist2=SQRT(DOT_PRODUCT(dist_vec_j,dist_vec_j))
                         IF(dist2.LE.th)dist2=th
                         tmp_const = (R12(iatom,jatom)**3)
                         dr_ij_dR(:)= dr1_r2(:) / tmp_const

                         !derivativ w.r.t. Rj
                         dr_j_dR=dist_vec_j(:) / dist2

                         dmy_dR_j(:)= -(dr_j_dR(:) / R12(iatom,jatom) - (dist1-dist2) * dr_ij_dR(:))

                         !derivativ w.r.t. Ri
                         dmy_dR_i(:)= dr_i_dR(:) / R12(iatom,jatom) - (dist1-dist2) * dr_ij_dR(:)

                         my1= (dist1-dist2)/R12(iatom,jatom)
                         polyn= 1.5_dp * my1 - 0.5_dp * my1**3

                         dpolyn= 1.5_dp - 1.5_dp * my1**2

                         tmp_const = (1.5_dp**2)*(dpolyn-dpolyn*(polyn**2))-&
                                      (1.5_dp)*((1.5_dp*polyn-0.5*(polyn**3))**2)*&
                                      (1.5_dp*dpolyn-1.5_dp*dpolyn*(polyn**2))

                         ds_dR(:,iatom,jatom)=- 0.5_dp * tmp_const * dmy_dR_i(:)
                         ds_dR(:,jatom,iatom)=- 0.5_dp *tmp_const  * dmy_dR_j(:)
                         f_3 =  1.5_dp * polyn - 0.5_dp * polyn**3
                         f_3 =   1.5_dp * f_3 - 0.5_dp * f_3**3
                         s_my_ij(iatom,jatom)= 0.5_dp * (1-f_3)
                         P_i(iatom)= P_i(iatom) * s_my_ij(iatom,jatom)

                      END IF
                   END DO

                   DO jatom = 1,particles%n_els
                      IF(iatom.NE.jatom)THEN
                         IF(ABS( s_my_ij(iatom,jatom) ).LE.th)s_my_ij(iatom,jatom)=s_my_ij(iatom,jatom)+th

                         dP_i_dRi(:,iatom)= dP_i_dRi(:,iatom)+&
                              P_i(iatom) / s_my_ij(iatom,jatom) * ds_dR(:,iatom,jatom)
                         dP_i_dRj(:,iatom,jatom)= P_i(iatom) / s_my_ij(iatom,jatom) * ds_dR(:,jatom,iatom)

                      END IF
                   END DO
                   d_sum_Pm_dR(:,iatom)= d_sum_Pm_dR(:,iatom)+ dP_i_dRi(:,iatom)

                   IF(is_constraint(iatom)) d_sum_const_dR(:,iatom)= d_sum_const_dR(:,iatom)+ dP_i_dRi(:,iatom)*&
                                                                     coefficients(iatom)
                   DO jatom = 1,particles%n_els
                      IF(iatom.NE.jatom)THEN
                         d_sum_Pm_dR(:,jatom)=d_sum_Pm_dR(:,jatom)+ dP_i_dRj(:,iatom,jatom)

                         IF(is_constraint(iatom))d_sum_const_dR(:,jatom)= d_sum_const_dR(:,jatom)+ dP_i_dRj(:,iatom,jatom)*&
                                                                     coefficients(iatom)
                      END IF
                   END DO
                ELSE
                   P_i(iatom)=0.0_dp

                END IF

             END DO

             sum_Pm=0.0_dp
             DO ip=1,particles%n_els
                sum_Pm=sum_Pm+P_i(ip)
             END DO

             sum_cell_f_constr=0.0_dp
             DO ip=1,qs_env%dft_control%qs_control%becke_control%natoms
                sum_cell_f_constr=sum_cell_f_constr+P_i(catom(ip))* qs_env%dft_control%qs_control%becke_control%coeff(ip)
             END DO

             DO iatom=1,particles%n_els
                IF(ABS(sum_Pm).EQ.0.0_dp)THEN
                   dw_dR(:,iatom)=0.0_dp
                ELSE
                   dw_dR(:,iatom)= d_sum_const_dR(:,iatom)/sum_Pm - &
                        sum_cell_f_constr * d_sum_Pm_dR(:,iatom) / (sum_Pm**2)
                END IF
             END DO

                DO iatom=1,particles%n_els
                   DO ispin=1 , qs_env%dft_control%nspins
                        integrated(:,iatom)=integrated(:,iatom) + dw_dR(:,iatom) * rho%rho_r(ispin)%pw%cr3d(k,j,i) * dvol
                   END DO
                END DO
             END DO
          END DO
       END DO

       DEALLOCATE (catom,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"catom")
       DEALLOCATE (R12,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"R12")
       DEALLOCATE (s_my_ij,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"s_my_ij")
       DEALLOCATE (ds_dR,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"ds_dR")
       DEALLOCATE (P_i,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"P_i")
       DEALLOCATE (dw_dR,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"dw_dR")
       DEALLOCATE (is_constraint,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"is_constraint")
       DEALLOCATE (d_sum_Pm_dR,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"d_sum")
       DEALLOCATE (d_sum_const_dR,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"d_sum_const")
       DEALLOCATE (dP_i_dRi,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"dP_i_dRj")
       DEALLOCATE (dP_i_dRj,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"dP_i_dRj")
       DEALLOCATE (coefficients,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"coefficients")
       CALL mp_sum(integrated,qs_env%para_env%group)

       DO iatom=1,natom
          ikind = kind_of(iatom)
          i = atom_of_kind(iatom)
          force(ikind)%rho_elec(:,i) =&
               force(ikind)%rho_elec(:,i) + integrated(:,iatom)*strength
       END DO

       DEALLOCATE (integrated,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"integrated")
       DEALLOCATE (atom_of_kind,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"atom_of_kind")
       DEALLOCATE (kind_of,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"kind_of")

     END IF

     CALL timestop(handle)

  END SUBROUTINE becke_force

END MODULE et_coupling

