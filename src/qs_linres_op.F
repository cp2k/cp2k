!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculate the operators p rxp and D needed in the optimization
!>      of the different contribution of the firs order response orbitals
!>      in a epr calculation
!> \note
!>      The interactions are considered only within the minimum image convention
!> \par History
!>       created 07-2005 [MI]
!> \author MI
! *****************************************************************************
MODULE qs_linres_op
  USE array_types,                     ONLY: array_i1d_obj,&
                                             array_new,&
                                             array_nullify,&
                                             array_release
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_array_i_utils,                ONLY: cp_2d_i_p_type
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_add_dist2d,&
                                             dbcsr_alloc_block_from_nbl,&
                                             dbcsr_allocate_matrix_set,&
                                             dbcsr_deallocate_matrix_set,&
                                             dbcsr_obj_p_type,&
                                             dbcsr_sm_fm_multiply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: &
       cp_fm_create, cp_fm_get_info, cp_fm_get_submatrix, cp_fm_p_type, &
       cp_fm_release, cp_fm_set_all, cp_fm_set_submatrix, cp_fm_to_fm, &
       cp_fm_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dbcsr_methods,                   ONLY: dbcsr_distribution_release,&
                                             dbcsr_init
  USE dbcsr_operations,                ONLY: dbcsr_copy,&
                                             dbcsr_set
  USE dbcsr_types,                     ONLY: dbcsr_distribution_obj,&
                                             dbcsr_type_antisymmetric,&
                                             dbcsr_type_normal
  USE dbcsr_util,                      ONLY: convert_offsets_to_sizes,&
                                             dbcsr_checksum,&
                                             dbcsr_create,&
                                             dbcsr_dist2d_to_dist
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE orbital_pointers,                ONLY: coset
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE qs_elec_field,                   ONLY: build_efg_matrix
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_fermi_contact,                ONLY: build_fermi_contact_matrix
  USE qs_linres_types,                 ONLY: current_env_type,&
                                             get_current_env,&
                                             get_issc_env,&
                                             issc_env_type,&
                                             linres_control_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_operators_ao,                 ONLY: build_ang_mom_matrix,&
                                             build_lin_mom_matrix,&
                                             rRc_xyz_ao
  USE qs_spin_orbit,                   ONLY: build_pso_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: current_operators,issc_operators,fac_vecp,ind_m2,set_vecp,set_vecp_rev,fm_scale_by_pbc_AC

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_op'

CONTAINS

! *****************************************************************************
!> \brief Calculate the first order hamiltonian applied to the ao
!>      and then apply them to the ground state orbitals,
!>      the h1_psi1 full matrices are then ready to solve the
!>      non-homogeneous linear equations that give the psi1
!>      linear response orbitals.
!> \note
!>      For the operators rxp and D the h1 depends on the psi0 to which
!>      is applied, or better the center of charge of the psi0 is
!>      used to define the position operator
!>      The centers of the orbitals result form the orbital localization procedure
!>      that typicaly uses the berry phase operator to define the Wannier centers.
!> \par History
!>      07.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE current_operators(current_env,qs_env,error)

    TYPE(current_env_type)                   :: current_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'current_operators', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iao, icenter, idir, ii, iii, ispin, istat, istate, j, &
      nao, natom, nbr_center(2), nmo, nsgf, nspins, nstates(2), output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf, last_sgf
    INTEGER, DIMENSION(:), POINTER           :: rbs
    LOGICAL                                  :: failure
    REAL(dp)                                 :: chk(3), ck(3), ckdk(3), dk(3)
    REAL(dp), DIMENSION(:, :), POINTER       :: basisfun_center, vecbuf_c0
    TYPE(array_i1d_obj)                      :: row_blk_sizes
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_2d_i_p_type), DIMENSION(:), &
      POINTER                                :: center_list
    TYPE(cp_2d_r_p_type), DIMENSION(3)       :: vecbuf_RmdC0
    TYPE(cp_2d_r_p_type), DIMENSION(:), &
      POINTER                                :: centers_set
    TYPE(cp_fm_p_type), DIMENSION(3)         :: fm_Rmd_mos
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: psi0_order
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: p_psi0, rxp_psi0
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: fm_work1, mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_distribution_obj)             :: dbcsr_dist
    TYPE(dbcsr_obj_p_type), DIMENSION(:), &
      POINTER                                :: op_ao
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_all
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: lr_section

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY(atomic_kind_set,cell,dft_control,linres_control,&
         &  logger,particle_set,lr_section,fm_work1,&
         &  basisfun_center,centers_set,center_list,p_psi0,&
         &  rxp_psi0,vecbuf_c0,psi0_order,&
         &  mo_coeff,op_ao,sab_all,distribution_2d)

    logger => cp_error_get_logger(error)
    lr_section => section_vals_get_subs_vals(qs_env%input,&
                  "PROPERTIES%LINRES",error=error)

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
                  extension=".linresLog",error=error)
    IF(output_unit>0) THEN
       WRITE(output_unit,FMT="(T2,A,/)")&
            "CURRENT| Calculation of the p and (r-d)xp operators applied to psi0"
    ENDIF

    CALL get_qs_env(qs_env=qs_env, &
                    atomic_kind_set=atomic_kind_set, &
                    cell=cell, &
                    dft_control=dft_control, &
                    linres_control=linres_control, &
                    para_env=para_env, &
                    particle_set=particle_set, &
                    sab_all=sab_all,&
                    distribution_2d=distribution_2d,&
                    error=error)

    nspins = dft_control%nspins

    CALL get_current_env(current_env=current_env, nao=nao, centers_set=centers_set,&
                         center_list=center_list,basisfun_center=basisfun_center, &
                         nbr_center=nbr_center,p_psi0=p_psi0, rxp_psi0=rxp_psi0,&
                         psi0_order=psi0_order,&
                         nstates=nstates,error=error)

    ALLOCATE(vecbuf_c0(1,nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO idir = 1,3
       NULLIFY(vecbuf_Rmdc0(idir)%array)
       ALLOCATE(vecbuf_Rmdc0(idir)%array(1,nao),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDDO

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,nsgf=nsgf)

    natom = SIZE(particle_set,1)
    ALLOCATE (first_sgf(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (last_sgf(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL get_particle_set(particle_set=particle_set,&
                          first_sgf=first_sgf,&
                          last_sgf=last_sgf,error=error)


    ! Calculate the (r - dk)xp operator applied to psi0k
    ! We use the distributive property of the vector product and calculatr
    ! (r-c)xp + (c-d)xp
    ! where c depends on the contracted functions and not on the states
    ! d is the center of a specific state and a loop over states is needed
    ! the second term can be added in a second moment as a correction
    ! notice: (r-c) and p are operators, whereas (c-d) is a multiplicative factor

    !    !First term: operator matrix elements
    !    CALL rmc_x_p_xyz_ao(op_rmd_ao,qs_env,minimum_image=.FALSE.,error=error)
    !************************************************************
    !


    !
    ! prepare for allocation
    CALL dbcsr_dist2d_to_dist (distribution_2d, dbcsr_dist, error)
    ALLOCATE (rbs(natom), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL convert_offsets_to_sizes (first_sgf, rbs, last_sgf)
    CALL array_nullify (row_blk_sizes)
    CALL array_new (row_blk_sizes, rbs, gift=.TRUE.)
    !
    !
    CALL dbcsr_allocate_matrix_set(op_ao,3,error=error)
    ALLOCATE(op_ao(1)%matrix,op_ao(2)%matrix,op_ao(3)%matrix)
    CALL dbcsr_init (op_ao(1)%matrix)
    CALL dbcsr_init (op_ao(2)%matrix)
    CALL dbcsr_init (op_ao(3)%matrix)

    CALL dbcsr_create(matrix=op_ao(1)%matrix, &
         name="op_ao", &
         dist=dbcsr_dist, matrix_type=dbcsr_type_normal,&
         row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
         nblks=0, nze=0, mutable_work=.TRUE., &
         error=error)
    CALL cp_dbcsr_add_dist2d (op_ao(1)%matrix%m, distribution_2d, error)
    CALL dbcsr_alloc_block_from_nbl(op_ao(1)%matrix,sab_all,error=error)
    CALL dbcsr_set(op_ao(1)%matrix,0.0_dp,error=error)

    DO idir=2,3
       CALL dbcsr_copy(op_ao(idir)%matrix,op_ao(1)%matrix,&
            "op_ao"//"-"//TRIM(ADJUSTL(cp_to_string(idir))),error=error)
       CALL dbcsr_set(op_ao(idir)%matrix,0.0_dp,error=error)
    ENDDO


    chk(:) = 0.0_dp
    DO ispin = 1,nspins
       mo_coeff => psi0_order(ispin)%matrix
       nmo = nstates(ispin)
       CALL cp_fm_set_all(p_psi0(ispin,1)%matrix,0.0_dp,error=error)
       CALL cp_fm_set_all(p_psi0(ispin,2)%matrix,0.0_dp,error=error)
       CALL cp_fm_set_all(p_psi0(ispin,3)%matrix,0.0_dp,error=error)
       DO icenter = 1,nbr_center(ispin)
          CALL dbcsr_set(op_ao(1)%matrix,0.0_dp,error=error)
          CALL dbcsr_set(op_ao(2)%matrix,0.0_dp,error=error)
          CALL dbcsr_set(op_ao(3)%matrix,0.0_dp,error=error)
          !CALL rmc_x_p_xyz_ao(op_ao,qs_env,minimum_image=.FALSE.,&
          !     &              wancen=centers_set(ispin)%array(1:3,icenter),&
          !     &              error=error)
          CALL build_ang_mom_matrix(qs_env,op_ao,centers_set(ispin)%array(1:3,icenter),error)
          !
          ! accumulate checksums
          chk(1) = chk(1)+dbcsr_checksum(op_ao(1)%matrix,error)
          chk(2) = chk(2)+dbcsr_checksum(op_ao(2)%matrix,error)
          chk(3) = chk(3)+dbcsr_checksum(op_ao(3)%matrix,error)
          DO idir = 1,3
             CALL cp_fm_set_all(rxp_psi0(ispin,idir)%matrix,0.0_dp,error=error)
             CALL dbcsr_sm_fm_multiply(op_ao(idir)%matrix,mo_coeff,&
                  &                    rxp_psi0(ispin,idir)%matrix,ncol=nmo,&
                  &                    alpha=-1.0_dp,error=error)
             DO j = center_list(ispin)%array(1,icenter),center_list(ispin)%array(1,icenter+1)-1
                istate = center_list(ispin)%array(2,j)
                CALL cp_fm_to_fm(rxp_psi0(ispin,idir)%matrix,&
                     &           p_psi0(ispin,idir)%matrix,1,istate,istate)
             ENDDO
          END DO
       ENDDO
       CALL cp_fm_to_fm(p_psi0(ispin,1)%matrix,rxp_psi0(ispin,1)%matrix,error)
       CALL cp_fm_to_fm(p_psi0(ispin,2)%matrix,rxp_psi0(ispin,2)%matrix,error)
       CALL cp_fm_to_fm(p_psi0(ispin,3)%matrix,rxp_psi0(ispin,3)%matrix,error)
    ENDDO
    !
    CALL dbcsr_deallocate_matrix_set(op_ao,error=error)
    !
    ! print checksums
    IF(output_unit>0) THEN
       WRITE(output_unit,'(T2,A,E23.16)') 'CURRENT| current_operators: CheckSum L_x =',chk(1)
       WRITE(output_unit,'(T2,A,E23.16)') 'CURRENT| current_operators: CheckSum L_y =',chk(2)
       WRITE(output_unit,'(T2,A,E23.16)') 'CURRENT| current_operators: CheckSum L_z =',chk(3)
    ENDIF
    !
    ! Calculate the px py pz operators
    CALL dbcsr_allocate_matrix_set(op_ao,3,error=error)
    ALLOCATE(op_ao(1)%matrix,op_ao(2)%matrix,op_ao(3)%matrix)
    CALL dbcsr_init (op_ao(1)%matrix)
    CALL dbcsr_init (op_ao(2)%matrix)
    CALL dbcsr_init (op_ao(3)%matrix)

    CALL dbcsr_create(matrix=op_ao(1)%matrix, &
         name="op_ao", &
         dist=dbcsr_dist, matrix_type=dbcsr_type_antisymmetric,&
         row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
         nblks=0, nze=0, mutable_work=.TRUE., &
         error=error)
    CALL cp_dbcsr_add_dist2d (op_ao(1)%matrix%m, distribution_2d, error)
    CALL dbcsr_alloc_block_from_nbl(op_ao(1)%matrix,sab_all,error=error)
    CALL dbcsr_set(op_ao(1)%matrix,0.0_dp,error=error)

    DO idir=2,3
       CALL dbcsr_copy(op_ao(idir)%matrix,op_ao(1)%matrix,&
            "op_ao"//"-"//TRIM(ADJUSTL(cp_to_string(idir))),error=error)
       CALL dbcsr_set(op_ao(idir)%matrix,0.0_dp,error=error)
    ENDDO
    !
    CALL build_lin_mom_matrix(qs_env,op_ao,error)
    !CALL p_xyz_ao(op_ao,qs_env,minimum_image=.FALSE.,error=error)
    !
    ! print checksums
    chk(1) = dbcsr_checksum(op_ao(1)%matrix,error)
    chk(2) = dbcsr_checksum(op_ao(2)%matrix,error)
    chk(3) = dbcsr_checksum(op_ao(3)%matrix,error)
    IF(output_unit>0) THEN
       WRITE(output_unit,'(T2,A,E23.16)') 'CURRENT| current_operators: CheckSum P_x =',chk(1)
       WRITE(output_unit,'(T2,A,E23.16)') 'CURRENT| current_operators: CheckSum P_y =',chk(2)
       WRITE(output_unit,'(T2,A,E23.16)') 'CURRENT| current_operators: CheckSum P_z =',chk(3)
    ENDIF
    ! Apply the p operator to the psi0
    DO idir = 1,3
      DO ispin = 1,nspins
         mo_coeff => psi0_order(ispin)%matrix
         nmo = nstates(ispin)
         CALL cp_fm_set_all(p_psi0(ispin,idir)%matrix,0.0_dp,error=error)
         CALL dbcsr_sm_fm_multiply(op_ao(idir)%matrix,mo_coeff,&
              &                    p_psi0(ispin,idir)%matrix,ncol=nmo,&
              &                    alpha=-1.0_dp,error=error)
      END DO
    END DO
    !
    CALL dbcsr_deallocate_matrix_set(op_ao,error=error)
    !
    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         "PRINT%PROGRAM_RUN_INFO",error=error)
    !
    ! Apply the (r-c)xp operator to the psi0
    !DO ispin = 1,nspins
    !  CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff, nmo=nmo, homo=homo)
    !  DO idir = 1,3
    !     CALL cp_fm_set_all(rxp_psi0(ispin,idir)%matrix,0.0_dp,error=error)
    !     CALL cp_sm_fm_multiply(op_rmd_ao(idir)%matrix,mo_coeff,&
    !            rxp_psi0(ispin,idir)%matrix,ncol=nmo,alpha=-1.0_dp,error=error)
    !  END DO
    !END DO

    !Calculate the second term of the operator state by state
    !!!! what follows is a way to avoid calculating the L matrix for each centers.
    !!!! not tested
    IF(.FALSE.) THEN
       DO ispin = 1, nspins
          !   Allocate full matrices as working storage in the calculation
          !   of the rxp operator matrix. 3 matrices for the 3 Cartesian direction
          !   plus one to apply the momentum oprator to the modified mos fm
          mo_coeff => psi0_order(ispin)%matrix
          nmo = nstates(ispin)
          NULLIFY(tmp_fm_struct)
          CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
               ncol_global=nmo,para_env=para_env,&
               context=mo_coeff%matrix_struct%context,error=error)
          DO idir = 1,3
             NULLIFY(fm_Rmd_mos(idir)%matrix)
             CALL cp_fm_create(fm_Rmd_mos(idir)%matrix,tmp_fm_struct,error=error)
          END DO
          CALL cp_fm_create (fm_work1,tmp_fm_struct,error=error)
          CALL cp_fm_struct_release ( tmp_fm_struct ,error=error)
          
          ! This part should be done better, using the full matrix distribution
          DO istate = 1,nmo
             CALL cp_fm_get_submatrix(psi0_order(ispin)%matrix,vecbuf_c0,1,istate,nao,1,&
                  &                   transpose=.TRUE.,error=error)
             !center of the localized psi0 state istate
             dk(1:3) = centers_set(ispin)%array(1:3,istate)
             DO idir = 1,3
                !  This loop should be distributed over the processors
                DO iao = 1,nao
                   ck(1:3) = basisfun_center(1:3,iao)
                   ckdk = pbc(dk,ck,cell)
                   vecbuf_Rmdc0(idir)%array(1,iao) = vecbuf_c0(1,iao)*ckdk(idir)
                END DO ! iao
                CALL cp_fm_set_submatrix(fm_Rmd_mos(idir)%matrix,vecbuf_Rmdc0(idir)%array,&
                     1,istate,nao,1,transpose=.TRUE.,error=error)
             END DO  ! idir
          END DO  ! istate
       
          DO idir = 1,3
             CALL set_vecp(idir,ii,iii)

             !Add the second term to the idir component
             CALL cp_fm_set_all(fm_work1,0.0_dp,error=error)
             CALL dbcsr_sm_fm_multiply(op_ao(iii)%matrix,fm_Rmd_mos(ii)%matrix,&
                  fm_work1,ncol=nmo,alpha=-1.0_dp,error=error)
             CALL cp_fm_scale_and_add(1.0_dp,rxp_psi0(ispin,idir)%matrix,&
                  1.0_dp, fm_work1,error=error)
          
             CALL cp_fm_set_all(fm_work1,0.0_dp,error=error)
             CALL dbcsr_sm_fm_multiply(op_ao(ii)%matrix,fm_Rmd_mos(iii)%matrix,&
                  fm_work1,ncol=nmo,alpha=-1.0_dp,error=error)
             CALL cp_fm_scale_and_add(1.0_dp,rxp_psi0(ispin,idir)%matrix,&
                  -1.0_dp,fm_work1,error=error)

          END DO ! idir

          DO idir = 1,3
             CALL cp_fm_release (fm_Rmd_mos(idir)%matrix,error=error)
          END DO
          CALL cp_fm_release (fm_work1,error=error)

       END DO  ! ispin
    ENDIF

    CALL dbcsr_distribution_release (dbcsr_dist)
    CALL array_release (row_blk_sizes)

    DEALLOCATE (first_sgf,last_sgf,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(vecbuf_c0,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO idir = 1,3
       DEALLOCATE(vecbuf_Rmdc0(idir)%array,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END DO

    CALL timestop(handle)

  END SUBROUTINE current_operators

  SUBROUTINE issc_operators(issc_env,qs_env,iatom,error)

    TYPE(issc_env_type)                      :: issc_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: iatom
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'issc_operators', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idir, ispin, nmo, &
                                                nspins, output_unit
    LOGICAL                                  :: do_dso, do_fc, do_pso, do_sd, &
                                                failure
    REAL(dp)                                 :: chk(20), r_i(3)
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: fc_psi0
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: dso_psi0, efg_psi0, pso_psi0
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_obj_p_type), DIMENSION(:), &
      POINTER                                :: matrix_dso, matrix_efg, &
                                                matrix_fc, matrix_pso
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: lr_section

    failure = .FALSE.

    CALL timeset(routineN,handle)

    NULLIFY(matrix_fc,matrix_pso,matrix_efg)
    NULLIFY(efg_psi0,pso_psi0,fc_psi0)

    logger => cp_error_get_logger(error)
    lr_section => section_vals_get_subs_vals(qs_env%input,&
                  "PROPERTIES%LINRES",error=error)

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
                  extension=".linresLog",error=error)

    CALL get_qs_env(qs_env=qs_env, &
                    atomic_kind_set=atomic_kind_set, &
                    cell=cell, &
                    dft_control=dft_control, &
                    linres_control=linres_control, &
                    para_env=para_env, &
                    mos=mos,&
                    particle_set=particle_set, &
                    error=error)

    nspins = dft_control%nspins

    CALL get_issc_env(issc_env=issc_env, &
                      matrix_efg=matrix_efg, & !this is used only here alloc/dealloc here???
                      matrix_pso=matrix_pso, & !this is used only here alloc/dealloc here???
                      matrix_fc=matrix_fc, & !this is used only here alloc/dealloc here???
                      matrix_dso=matrix_dso, & !this is used only here alloc/dealloc here???
                      efg_psi0=efg_psi0, &
                      pso_psi0=pso_psi0, &
                      dso_psi0=dso_psi0, &
                      fc_psi0=fc_psi0, &
                      do_fc=do_fc,&
                      do_sd=do_sd,&
                      do_pso=do_pso,&
                      do_dso=do_dso,&
                      error=error)
    !
    !
    r_i = particle_set(iatom)%r !pbc(particle_set(iatom)%r,cell)
    !write(*,*) 'issc_operators iatom=',iatom,' r_i=',r_i
    chk = 0.0_dp
    !
    !
    !
    ! Fermi contact integral
    !IF(do_fc) THEN
    IF(.TRUE.) THEN ! for the moment we build it (regs)
       CALL dbcsr_set(matrix_fc(1)%matrix,0.0_dp,error=error)
       CALL build_fermi_contact_matrix(qs_env,matrix_fc,r_i,error)

       chk(1) = dbcsr_checksum(matrix_fc(1)%matrix,error=error)

       IF(output_unit>0) THEN
          WRITE(output_unit,'(T2,A,E23.16)') 'ISSC| fermi_contact: CheckSum =',chk(1)
       ENDIF
    ENDIF
    !
    ! spin-orbit integral
    !IF(do_pso) THEN
    IF(.TRUE.) THEN ! for the moment we build it (regs)
       CALL dbcsr_set(matrix_pso(1)%matrix,0.0_dp,error=error)
       CALL dbcsr_set(matrix_pso(2)%matrix,0.0_dp,error=error)
       CALL dbcsr_set(matrix_pso(3)%matrix,0.0_dp,error=error)
       CALL build_pso_matrix(qs_env,matrix_pso,r_i,error)

       chk(2) = dbcsr_checksum(matrix_pso(1)%matrix,error=error)
       chk(3) = dbcsr_checksum(matrix_pso(2)%matrix,error=error)
       chk(4) = dbcsr_checksum(matrix_pso(3)%matrix,error=error)

       IF(output_unit>0) THEN
          WRITE(output_unit,'(T2,A,E23.16)') 'ISSC| pso_x: CheckSum =',chk(2)
          WRITE(output_unit,'(T2,A,E23.16)') 'ISSC| pso_y: CheckSum =',chk(3)
          WRITE(output_unit,'(T2,A,E23.16)') 'ISSC| pso_z: CheckSum =',chk(4)
       ENDIF
    ENDIF
    !
    ! electric field integral
    !IF(do_sd) THEN
    IF(.TRUE.) THEN ! for the moment we build it (regs)
       CALL dbcsr_set(matrix_efg(1)%matrix,0.0_dp,error=error)
       CALL dbcsr_set(matrix_efg(2)%matrix,0.0_dp,error=error)
       CALL dbcsr_set(matrix_efg(3)%matrix,0.0_dp,error=error)
       CALL dbcsr_set(matrix_efg(4)%matrix,0.0_dp,error=error)
       CALL dbcsr_set(matrix_efg(5)%matrix,0.0_dp,error=error)
       CALL dbcsr_set(matrix_efg(6)%matrix,0.0_dp,error=error)
       CALL build_efg_matrix(qs_env,matrix_efg,r_i,error)

       chk(5) = dbcsr_checksum(matrix_efg(1)%matrix,error=error)
       chk(6) = dbcsr_checksum(matrix_efg(2)%matrix,error=error)
       chk(7) = dbcsr_checksum(matrix_efg(3)%matrix,error=error)
       chk(8) = dbcsr_checksum(matrix_efg(4)%matrix,error=error)
       chk(9) = dbcsr_checksum(matrix_efg(5)%matrix,error=error)
       chk(10)= dbcsr_checksum(matrix_efg(6)%matrix,error=error)

       IF(output_unit>0) THEN
          WRITE(output_unit,'(T2,A,E23.16)') 'ISSC| efg (3xx-rr)/3: CheckSum =',chk(5)
          WRITE(output_unit,'(T2,A,E23.16)') 'ISSC| efg (3yy-rr)/3: CheckSum =',chk(6)
          WRITE(output_unit,'(T2,A,E23.16)') 'ISSC| efg (3zz-rr)/3: CheckSum =',chk(7)
          WRITE(output_unit,'(T2,A,E23.16)') 'ISSC| efg xy: CheckSum =',chk(8)
          WRITE(output_unit,'(T2,A,E23.16)') 'ISSC| efg xz: CheckSum =',chk(9)
          WRITE(output_unit,'(T2,A,E23.16)') 'ISSC| efg yz: CheckSum =',chk(10)
       ENDIF
    ENDIF
    !
    !
    IF(output_unit>0) THEN
       WRITE(output_unit,'(T2,A,E23.16)') 'ISSC| all operator: CheckSum =',SUM(chk(1:10))
    ENDIF
    !
    !>>> debugging only  here we build the dipole matrix... debugging the kernel...
    IF(do_dso) THEN
       CALL dbcsr_set(matrix_dso(1)%matrix,0.0_dp,error=error)
       CALL dbcsr_set(matrix_dso(2)%matrix,0.0_dp,error=error)
       CALL dbcsr_set(matrix_dso(3)%matrix,0.0_dp,error=error)
       CALL rRc_xyz_ao(matrix_dso,qs_env,(/0.0_dp,0.0_dp,0.0_dp/),1,error=error)
    ENDIF
    !
    ! multiply by the mos
    DO ispin = 1,nspins
       !
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff)
       CALL cp_fm_get_info(mo_coeff,ncol_global=nmo,error=error)
       !
       ! EFG
       IF(do_sd) THEN
          DO idir = 1,6
             CALL dbcsr_sm_fm_multiply(matrix_efg(idir)%matrix,mo_coeff,&
                  &                    efg_psi0(ispin,idir)%matrix,ncol=nmo,&
                  &                    alpha=1.0_dp,error=error)
          ENDDO
       ENDIF
       !
       ! PSO
       IF(do_pso) THEN
          DO idir = 1,3
             CALL dbcsr_sm_fm_multiply(matrix_pso(idir)%matrix,mo_coeff,&
                  &                   pso_psi0(ispin,idir)%matrix,ncol=nmo,&
                  &                   alpha=-1.0_dp,error=error)
          ENDDO
       ENDIF
       !
       ! FC
       IF(do_fc) THEN
          CALL dbcsr_sm_fm_multiply(matrix_fc(1)%matrix,mo_coeff,&
               &                    fc_psi0(ispin)%matrix,ncol=nmo,&
               &                    alpha=1.0_dp,error=error)
       ENDIF
       !
       !>>> for debugging only
       IF(do_dso) THEN
          DO idir = 1,3
             CALL dbcsr_sm_fm_multiply(matrix_dso(idir)%matrix,mo_coeff,&
                  &                    dso_psi0(ispin,idir)%matrix,ncol=nmo,&
                  &                    alpha=-1.0_dp,error=error)
          ENDDO
       ENDIF
       !<<< for debugging only
    ENDDO

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         "PRINT%PROGRAM_RUN_INFO",error=error)

    CALL timestop(handle)

  END SUBROUTINE issc_operators

! *****************************************************************************
  FUNCTION fac_vecp(a,b,c) RESULT(factor)

    INTEGER                                  :: a, b, c
    REAL(dp)                                 :: factor

     factor = 0.0_dp

     IF ((b.eq.a+1 .OR. b.eq.a-2) .AND. (c.eq.b+1 .OR. c.eq.b-2)) THEN
       factor = 1.0_dp
     ELSEIF ((b.eq.a-1 .OR. b.eq.a+2) .AND. (c.eq.b-1 .OR. c.eq.b+2)) THEN
       factor = -1.0_dp
     END IF

  END FUNCTION fac_vecp

! *****************************************************************************
  FUNCTION ind_m2(ii,iii) RESULT(i)

    INTEGER                                  :: ii, iii, i

    INTEGER                                  :: l(3)

     i = 0
     l(1:3) = 0
     IF(ii==0) THEN
       l(iii) = 1
     ELSEIF(iii==0) THEN
       l(ii) = 1
     ELSEIF(ii==iii) THEN
       l(ii) = 2
       i = coset(l(1),l(2),l(3)) - 1
     ELSE
       l(ii) = 1
       l(iii) = 1
     ENDIF
     i = coset(l(1),l(2),l(3)) - 1
  END FUNCTION ind_m2

! *****************************************************************************
  SUBROUTINE set_vecp(i1,i2,i3)

    INTEGER, INTENT(IN)                      :: i1
    INTEGER, INTENT(OUT)                     :: i2, i3

    IF(i1==1) THEN
      i2=2
      i3=3
    ELSEIF(i1==2) THEN
      i2=3
      i3=1
    ELSEIF(i1==3) THEN
      i2=1
      i3=2
    ELSE
    END IF

  END SUBROUTINE set_vecp
! *****************************************************************************
  SUBROUTINE set_vecp_rev(i1,i2,i3)

    INTEGER, INTENT(IN)                      :: i1, i2
    INTEGER, INTENT(OUT)                     :: i3

    IF((i1+i2)==3) THEN
      i3=3
    ELSEIF((i1+i2)==4) THEN
      i3=2
    ELSEIF((i1+i2)==5) THEN
      i3=1
    ELSE
    END IF

  END SUBROUTINE set_vecp_rev


! *****************************************************************************
!> \brief scale a matrix as a_ij = a_ij * pbc(rc(:,j),ra(:,i))(ixyz)
!> \author vw
! *****************************************************************************
  SUBROUTINE fm_scale_by_pbc_AC(matrix,ra,rc,cell,ixyz,error)
    TYPE(cp_fm_type), POINTER                :: matrix
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(in)                             :: ra, rc
    TYPE(cell_type), POINTER                 :: cell
    INTEGER, INTENT(IN)                      :: ixyz
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fm_scale_by_pbc_AC', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, icol_global, icol_local, irow_global, irow_local, m, &
      mypcol, myprow, n, ncol_block, ncol_global, ncol_local, npcol, nprow, &
      nrow_block, nrow_global, nrow_local
    REAL(KIND=dp)                            :: dist(3), rra(3), rrc(3)
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a
#if defined(__SCALAPACK)
    INTEGER, EXTERNAL                          :: indxl2g
#endif

    CALL timeset(routineN,handle)

    myprow=matrix%matrix_struct%context%mepos(1)
    mypcol=matrix%matrix_struct%context%mepos(2)
    nprow=matrix%matrix_struct%context%num_pe(1)
    npcol=matrix%matrix_struct%context%num_pe(2)

    nrow_block = matrix%matrix_struct%nrow_block
    ncol_block = matrix%matrix_struct%ncol_block
    nrow_global =matrix%matrix_struct%nrow_global
    ncol_global =matrix%matrix_struct%ncol_global
    nrow_local = matrix%matrix_struct%nrow_locals(myprow)
    ncol_local = matrix%matrix_struct%ncol_locals(mypcol)

    n=SIZE(rc,2)
    m=SIZE(ra,2)

    a => matrix%local_data
#if defined(__SCALAPACK)
    DO icol_local=1,ncol_local
       icol_global = indxl2g(icol_local,ncol_block,mypcol,&
            matrix%matrix_struct%first_p_pos(2),npcol)
       IF(icol_global.GT.n) CYCLE
       rrc = rc(:,icol_global)
       DO irow_local=1,nrow_local
          irow_global = indxl2g(irow_local,nrow_block,myprow,&
               matrix%matrix_struct%first_p_pos(1),nprow)
          IF(irow_global.GT.m) CYCLE
          rra = ra(:,irow_global)
          dist = pbc(rrc,rra,cell)
          a(irow_local,icol_local) = a(irow_local,icol_local) * dist(ixyz)
       ENDDO
    ENDDO
#else
    DO icol_local=1,MIN(n,ncol_global)
       rrc = rc(:,icol_local)
       DO irow_local=1,MIN(m,nrow_global)
          rra = ra(:,irow_local)
          dist = pbc(rrc,rra,cell)
          a(irow_local,icol_local) = a(irow_local,icol_local) * dist(ixyz)
       ENDDO
    ENDDO
#endif

    CALL timestop(handle)

  END SUBROUTINE fm_scale_by_pbc_AC

END MODULE qs_linres_op

