!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/et_coupling [1.0] *
!!
!!   NAME
!!     et_coupling
!!
!!   FUNCTION
!!     calculates the electron transfer coupling elements
!!     Wu, Van Voorhis, JCP 125, 164105 (2006)
!!
!!   AUTHOR
!!     fschiff (01.2007)
!!
!!
!!***
!******************************************************************************

MODULE et_coupling

  USE cell_types,                      ONLY: cell_type
  USE coefficient_types,               ONLY: coeff_type
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
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE cp_sm_pool_types,                ONLY: cp_sm_pool_p_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_list_get,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: diamat_all
  USE message_passing,                 ONLY: mp_sum
  USE particle_list_types,             ONLY: particle_list_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_matrix_pools,                 ONLY: mpools_get
  USE qs_mo_types,                     ONLY: get_mo_set
  USE sparse_matrix_types,             ONLY: deallocate_matrix_set
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'et_coupling'

! *** Public subroutines ***

  PUBLIC :: calc_et_coupling,becke_constraint




CONTAINS
  SUBROUTINE calc_et_coupling(qs_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, j, k, nao, ncol_local, &
                                                nmo, nrow_local, iw
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: sign
    REAL(KIND=dp)                            :: Sda, strength, tttt, Waa, &
                                                Wbb, Wda
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
    TYPE(cp_sm_pool_p_type), DIMENSION(:), &
      POINTER                                :: S_sm_pools
    TYPE(section_vals_type), POINTER         :: et_coupling_section

    NULLIFY(tmp2,SMO,rest_MO,mo_mo_fmstruct, S_sm_pools,rest_MO)
    logger => cp_error_get_logger(error)
    et_coupling_section => section_vals_get_subs_vals(qs_env%input,&
                            "PROPERTIES%ET_COUPLING",error=error)
    iw=cp_print_key_unit_nr(logger,et_coupling_section,"PROGRAM_RUN_INFO",&
                                         extension=".log",error=error)

    ALLOCATE(rest_MO(2))
    ALLOCATE(a(qs_env%dft_control%nspins))
    ALLOCATE(b(qs_env%dft_control%nspins))
    ALLOCATE(S_det(qs_env%dft_control%nspins))


    CALL mpools_get(qs_env%mpools, mo_mo_fm_pools=mo_mo_fm_pools,&
         S_sm_pools=S_sm_pools,&
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

    CALL cp_sm_fm_multiply(qs_env%matrix_s(1)%matrix,qs_env%et_coupling%et_mo_coeff(i)%matrix,&
         tmp2,nmo,1.0_dp,0.0_dp,error=error)
    CALL cp_fm_gemm('T','N',nmo,nmo,nao,1.0_dp,&
         qs_env%mos(i)%mo_set%mo_coeff,&
         tmp2,0.0_dp,SMO,error=error)

!    calculate the MO-representation of the restraint matrix A
         
    CALL cp_sm_fm_multiply(qs_env%et_coupling%rest_mat(1)%matrix,&
         qs_env%et_coupling%et_mo_coeff(i)%matrix,&
         tmp2,nmo,1.0_dp,0.0_dp,error=error)   
    
    CALL cp_fm_gemm('T','N',nmo,nmo,nao,1.0_dp,&
         qs_env%mos(i)%mo_set%mo_coeff,&
         tmp2,0.0_dp,rest_MO(1)%matrix,error=error)

!    calculate the MO-representation of the restraint matrix D
             
    CALL cp_sm_fm_multiply(qs_env%et_coupling%rest_mat(2)%matrix,&
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
          tttt=(-1.0_dp)**MOD(row_indices(k)+col_indices(j),2)
          b(i)=b(i)+rest_MO(2)%matrix%local_data(k,j)*inverse_mat%local_data(k,j)*tttt
       END DO
    END DO

    CALL cp_fm_transpose(inverse_mat,Tinverse,error)
    a(i)=0.0_dp
    DO j=1,ncol_local
       DO k=1,nrow_local
          tttt=(-1.0_dp)**MOD(row_indices(k)+col_indices(j),2)
          a(i)=a(i)+rest_MO(1)%matrix%local_data(k,j)*Tinverse%local_data(k,j)*tttt
       END DO
    END DO
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

    Sda=S_det(1)*S_det(2)
    Wda=((a(1)+a(2))+(b(1)+b(2)))*0.5_dp*Sda
    IF(qs_env%dft_control%qs_control%ddapc_restraint)THEN
    Waa=qs_env%et_coupling%order_p
    Wbb=qs_env%dft_control%qs_control%ddapc_restraint_control%ddapc_order_p
    strength=qs_env%dft_control%qs_control%ddapc_restraint_control%strength    
    END IF
    IF(qs_env%dft_control%qs_control%becke_constraint)THEN
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
    CALL diamat_all (S_mat,eigenv,.TRUE.)
    ! U = S**(-1/2)
    U=0.0_dp
    U(1,1)=1.0_dp/SQRT(eigenv(1))
    U(2,2)=1.0_dp/SQRT(eigenv(2))
    tmp_mat=MATMUL(U,TRANSPOSE(S_mat))
    U=MATMUL(S_mat,tmp_mat)
    tmp_mat=MATMUL(W_mat,U)
    W_mat=MATMUL(U,tmp_mat)
    CALL diamat_all (W_mat,eigenv,.TRUE.)
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
            'Diabatic electronic coupling matrix element(mHartree):',W_mat(1,2)*1000.0_dp       
       
    ENDIF

    CALL deallocate_matrix_set(qs_env%et_coupling%rest_mat,error=error)

    CALL cp_print_key_finished_output(iw,logger,et_coupling_section,&
               "PROGRAM_RUN_INFO", error=error)

  END SUBROUTINE calc_et_coupling




!!****** cp2k/becke_constraint [1.0] *
!!
!!   NAME
!!     becke_constraint
!!
!!   FUNCTION
!!     calculates a becke contraint
!!
!!   AUTHOR
!!     fschiff (01.2007)
!!
!!
!!***
!******************************************************************************
  SUBROUTINE becke_constraint(qs_env,becke_const ,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(coeff_type)                         :: becke_const
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'becke_constraint', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iatom, ip, j, jatom, k, &
                                                np(3)
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
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_list_type), POINTER        :: particles

  failure=.FALSE.


  IF(qs_env%dft_control%qs_control%becke_constraint)THEN
     NULLIFY(cell,particles)
     CALL cp_subsys_get(qs_env%subsys,particles=particles,error=error)
     becke_const%pw%cr3d=0.0_dp 
     
     CALL get_qs_env(qs_env,cell=cell,error=error)

     ALLOCATE(catom(qs_env%dft_control%qs_control%becke_control%natoms))
     DO i=1,qs_env%dft_control%qs_control%becke_control%natoms
        catom(i)=qs_env%dft_control%qs_control%becke_control%atoms(i)
     ENDDO

     target_val=qs_env%dft_control%qs_control%becke_control%TARGET
     strength=qs_env%dft_control%qs_control%becke_control%strength
     bo=becke_const%pw%pw_grid%bounds_local
     np=becke_const%pw%pw_grid%npts
     dvol = becke_const%pw%pw_grid%dvol
     dr=becke_const%pw%pw_grid%dr
     np= becke_const%pw%pw_grid%npts
     shift=-REAL(MODULO(np,2),dp)*dr/2.0_dp
     ALLOCATE(cell_functions(particles%n_els))
!calculate distances from target atom (only implemented for a diatomic system)

     DO i=1,3
        cell_v(i)=cell%hmat(i,i)
     END DO
     
     ALLOCATE(R12(particles%n_els,particles%n_els))
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
              DO iatom=1,particles%n_els-1
                 r=particles%els(iatom)%r
                 DO ip=1,3
                    r(ip)=MODULO(r(ip),cell%hmat(ip,ip))-cell%hmat(ip,ip)/2._dp
                 END DO
                 dist_vec=(r-grid_p)-ANINT((r-grid_p)/cell_v)*cell_v
                 dist1=SQRT(DOT_PRODUCT(dist_vec,dist_vec))

                 DO jatom=iatom+1,particles%n_els
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
                    cell_functions(jatom)=cell_functions(jatom)*0.5_dp*(1+myexp)
                 END DO
              END DO

              sum_cell_f_constr=0.0_dp
              DO ip=1,qs_env%dft_control%qs_control%becke_control%natoms
                 sum_cell_f_constr=sum_cell_f_constr+cell_functions(catom(ip))
              END DO

              sum_cell_f_all=0.0_dp
              DO ip=1,particles%n_els
                 sum_cell_f_all=sum_cell_f_all+cell_functions(ip)
              END DO
              
              becke_const%pw%cr3d(k,j,i)=sum_cell_f_constr/sum_cell_f_all
                       
                 
           END DO
        END DO
     END DO

     dE=0.0_dp
     DO i=1 , qs_env%dft_control%nspins
           dE=dE+accurate_sum(becke_const%pw%cr3d*qs_env%rho%rho_r(i)%pw%cr3d)*dvol
     END DO
     CALL mp_sum(dE,qs_env%para_env%group)
     qs_env%dft_control%qs_control%becke_control%becke_order_p=dE
     qs_env%energy%becke=(dE-target_val)*strength
     DEALLOCATE(catom)
     DEALLOCATE(R12)
     DEALLOCATE(cell_functions)
     
     
  END IF
END SUBROUTINE becke_constraint

    
END MODULE et_coupling

