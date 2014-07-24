!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cell_methods,                    ONLY: cell_type,&
                                             pbc
  USE cp_array_i_utils,                ONLY: cp_2d_i_p_type
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_cfm_basic_linalg,             ONLY: cp_cfm_solve
  USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                             cp_cfm_get_info,&
                                             cp_cfm_p_type,&
                                             cp_cfm_release,&
                                             cp_cfm_set_all
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
  USE cp_dbcsr_interface,              ONLY: &
       array_i1d_obj, array_new, array_nullify, array_release, &
       convert_offsets_to_sizes, cp_dbcsr_allocate_matrix_set, &
       cp_dbcsr_checksum, cp_dbcsr_copy, cp_dbcsr_create, &
       cp_dbcsr_deallocate_matrix, cp_dbcsr_deallocate_matrix_set, &
       cp_dbcsr_init, cp_dbcsr_p_type, cp_dbcsr_set, cp_dbcsr_type, &
       dbcsr_distribution_obj, dbcsr_type_antisymmetric, &
       dbcsr_type_no_symmetry
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: &
       cp_fm_create, cp_fm_get_info, cp_fm_get_submatrix, cp_fm_p_type, &
       cp_fm_release, cp_fm_set_all, cp_fm_set_submatrix, cp_fm_to_fm, &
       cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: twopi
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
                                             get_polar_env,&
                                             issc_env_type,&
                                             linres_control_type,&
                                             polar_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_moments,                      ONLY: build_berry_moment_matrix
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_operators_ao,                 ONLY: build_ang_mom_matrix,&
                                             build_lin_mom_matrix,&
                                             rRc_xyz_ao
  USE qs_spin_orbit,                   ONLY: build_pso_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: current_operators,issc_operators,fac_vecp,ind_m2,set_vecp,set_vecp_rev,&
            fm_scale_by_pbc_AC, polar_operators

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_op'

CONTAINS

! *****************************************************************************
!> \brief Calculate the first order hamiltonian applied to the ao
!>      and then apply them to the ground state orbitals,
!>      the h1_psi1 full matrices are then ready to solve the
!>      non-homogeneous linear equations that give the psi1
!>      linear response orbitals.
!> \param current_env ...
!> \param qs_env ...
!> \param error ...
!> \par History
!>      07.2005 created [MI]
!> \author MI
!> \note
!>      For the operators rxp and D the h1 depends on the psi0 to which
!>      is applied, or better the center of charge of the psi0 is
!>      used to define the position operator
!>      The centers of the orbitals result form the orbital localization procedure
!>      that typicaly uses the berry phase operator to define the Wannier centers.
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
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: op_ao
    TYPE(cp_fm_p_type), DIMENSION(3)         :: fm_Rmd_mos
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: psi0_order
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: p_psi0, rxp_psi0
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: fm_work1, mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_distribution_obj), POINTER    :: dbcsr_dist
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_all, sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: lr_section

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY(atomic_kind_set,cell,dft_control,linres_control,&
            logger,particle_set,lr_section,fm_work1,&
            basisfun_center,centers_set,center_list,p_psi0,&
            rxp_psi0,vecbuf_c0,psi0_order,&
            mo_coeff,op_ao,sab_all)

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
                    sab_orb=sab_orb,&
                    dbcsr_dist=dbcsr_dist,&
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
    ! One possible way to go is to use the distributive property of the vector product and calculatr
    ! (r-c)xp + (c-d)xp
    ! where c depends on the contracted functions and not on the states
    ! d is the center of a specific state and a loop over states is needed
    ! the second term can be added in a second moment as a correction
    ! notice: (r-c) and p are operators, whereas (c-d) is a multiplicative factor

    !    !First term: operator matrix elements
    !    CALL rmc_x_p_xyz_ao(op_rmd_ao,qs_env,minimum_image=.FALSE.,error=error)
    !************************************************************
    !
    ! Since many psi0 vector can have the same center, depending on how the center is selected,
    ! the (r - dk)xp operator matrix is computed Ncenter times,
    ! where Ncenter is the total number of different centers
    ! and each time it is multiplied by all the psi0 with center dk to get the rxp_psi0 matrix


    !
    ! prepare for allocation
    ALLOCATE (rbs(natom), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL convert_offsets_to_sizes (first_sgf, rbs, last_sgf)
    CALL array_nullify (row_blk_sizes)
    CALL array_new (row_blk_sizes, rbs, gift=.TRUE.)
    !
    !
    CALL cp_dbcsr_allocate_matrix_set(op_ao,3,error=error)
    ALLOCATE(op_ao(1)%matrix,op_ao(2)%matrix,op_ao(3)%matrix)
    CALL cp_dbcsr_init (op_ao(1)%matrix, error=error)
    CALL cp_dbcsr_init (op_ao(2)%matrix, error=error)
    CALL cp_dbcsr_init (op_ao(3)%matrix, error=error)

    CALL cp_dbcsr_create(matrix=op_ao(1)%matrix, &
         name="op_ao", &
         dist=dbcsr_dist, matrix_type=dbcsr_type_no_symmetry,&
         row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
         nblks=0, nze=0, mutable_work=.TRUE., &
         error=error)
    CALL cp_dbcsr_alloc_block_from_nbl(op_ao(1)%matrix,sab_all,error=error)
    CALL cp_dbcsr_set(op_ao(1)%matrix,0.0_dp,error=error)

    DO idir=2,3
       CALL cp_dbcsr_copy(op_ao(idir)%matrix,op_ao(1)%matrix,&
            "op_ao"//"-"//TRIM(ADJUSTL(cp_to_string(idir))),error=error)
       CALL cp_dbcsr_set(op_ao(idir)%matrix,0.0_dp,error=error)
    ENDDO


    chk(:) = 0.0_dp
    DO ispin = 1,nspins
       mo_coeff => psi0_order(ispin)%matrix
       nmo = nstates(ispin)
       CALL cp_fm_set_all(p_psi0(ispin,1)%matrix,0.0_dp,error=error)
       CALL cp_fm_set_all(p_psi0(ispin,2)%matrix,0.0_dp,error=error)
       CALL cp_fm_set_all(p_psi0(ispin,3)%matrix,0.0_dp,error=error)
       DO icenter = 1,nbr_center(ispin)
          CALL cp_dbcsr_set(op_ao(1)%matrix,0.0_dp,error=error)
          CALL cp_dbcsr_set(op_ao(2)%matrix,0.0_dp,error=error)
          CALL cp_dbcsr_set(op_ao(3)%matrix,0.0_dp,error=error)
          !CALL rmc_x_p_xyz_ao(op_ao,qs_env,minimum_image=.FALSE.,&
          !     &              wancen=centers_set(ispin)%array(1:3,icenter),&
          !     &              error=error)
          CALL build_ang_mom_matrix(qs_env,op_ao,centers_set(ispin)%array(1:3,icenter),error)
          !
          ! accumulate checksums
          chk(1) = chk(1)+cp_dbcsr_checksum(op_ao(1)%matrix,error=error)
          chk(2) = chk(2)+cp_dbcsr_checksum(op_ao(2)%matrix,error=error)
          chk(3) = chk(3)+cp_dbcsr_checksum(op_ao(3)%matrix,error=error)
          DO idir = 1,3
             CALL cp_fm_set_all(rxp_psi0(ispin,idir)%matrix,0.0_dp,error=error)
             CALL cp_dbcsr_sm_fm_multiply(op_ao(idir)%matrix,mo_coeff,&
                                       rxp_psi0(ispin,idir)%matrix,ncol=nmo,&
                                       alpha=-1.0_dp,error=error)
             DO j = center_list(ispin)%array(1,icenter),center_list(ispin)%array(1,icenter+1)-1
                istate = center_list(ispin)%array(2,j)
           ! the p_psi0 fm is used as temporary matrix to store the results for the psi0 centered in dk
                CALL cp_fm_to_fm(rxp_psi0(ispin,idir)%matrix,&
                                 p_psi0(ispin,idir)%matrix,1,istate,istate)
             ENDDO
          END DO
       ENDDO
       CALL cp_fm_to_fm(p_psi0(ispin,1)%matrix,rxp_psi0(ispin,1)%matrix,error)
       CALL cp_fm_to_fm(p_psi0(ispin,2)%matrix,rxp_psi0(ispin,2)%matrix,error)
       CALL cp_fm_to_fm(p_psi0(ispin,3)%matrix,rxp_psi0(ispin,3)%matrix,error)
    ENDDO
    !
    CALL cp_dbcsr_deallocate_matrix_set(op_ao,error=error)
    !
    ! print checksums
    IF(output_unit>0) THEN
       WRITE(output_unit,'(T2,A,E23.16)') 'CURRENT| current_operators: CheckSum L_x =',chk(1)
       WRITE(output_unit,'(T2,A,E23.16)') 'CURRENT| current_operators: CheckSum L_y =',chk(2)
       WRITE(output_unit,'(T2,A,E23.16)') 'CURRENT| current_operators: CheckSum L_z =',chk(3)
    ENDIF
    !
    ! Calculate the px py pz operators
    CALL cp_dbcsr_allocate_matrix_set(op_ao,3,error=error)
    ALLOCATE(op_ao(1)%matrix,op_ao(2)%matrix,op_ao(3)%matrix)
    CALL cp_dbcsr_init (op_ao(1)%matrix, error=error)
    CALL cp_dbcsr_init (op_ao(2)%matrix, error=error)
    CALL cp_dbcsr_init (op_ao(3)%matrix, error=error)

    CALL cp_dbcsr_create(matrix=op_ao(1)%matrix, &
         name="op_ao", &
         dist=dbcsr_dist, matrix_type=dbcsr_type_antisymmetric,&
         row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
         nblks=0, nze=0, mutable_work=.TRUE., &
         error=error)
    CALL cp_dbcsr_alloc_block_from_nbl(op_ao(1)%matrix,sab_orb,error=error)
    CALL cp_dbcsr_set(op_ao(1)%matrix,0.0_dp,error=error)

    DO idir=2,3
       CALL cp_dbcsr_copy(op_ao(idir)%matrix,op_ao(1)%matrix,&
            "op_ao"//"-"//TRIM(ADJUSTL(cp_to_string(idir))),error=error)
       CALL cp_dbcsr_set(op_ao(idir)%matrix,0.0_dp,error=error)
    ENDDO
    !
    CALL build_lin_mom_matrix(qs_env,op_ao,error)
    !CALL p_xyz_ao(op_ao,qs_env,minimum_image=.FALSE.,error=error)
    !
    ! print checksums
    chk(1) = cp_dbcsr_checksum(op_ao(1)%matrix,error=error)
    chk(2) = cp_dbcsr_checksum(op_ao(2)%matrix,error=error)
    chk(3) = cp_dbcsr_checksum(op_ao(3)%matrix,error=error)
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
         CALL cp_dbcsr_sm_fm_multiply(op_ao(idir)%matrix,mo_coeff,&
                                   p_psi0(ispin,idir)%matrix,ncol=nmo,&
                                   alpha=-1.0_dp,error=error)
      END DO
    END DO
    !
    CALL cp_dbcsr_deallocate_matrix_set(op_ao,error=error)
    !
    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         "PRINT%PROGRAM_RUN_INFO",error=error)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  This part is not necessary with the present implementation
    !  the angular momentum operator is computed directly for each dk independently
    !  and multiplied by the proper psi0 (i.e. those centered in dk)
    !  If Wannier centers are used, and no grouping of states with close centers is applied
    !  the (r-dk)xp operator is computed Nstate times and each time applied to only one vector psi0
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
                                      transpose=.TRUE.,error=error)
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
             CALL cp_dbcsr_sm_fm_multiply(op_ao(iii)%matrix,fm_Rmd_mos(ii)%matrix,&
                  fm_work1,ncol=nmo,alpha=-1.0_dp,error=error)
             CALL cp_fm_scale_and_add(1.0_dp,rxp_psi0(ispin,idir)%matrix,&
                  1.0_dp, fm_work1,error=error)

             CALL cp_fm_set_all(fm_work1,0.0_dp,error=error)
             CALL cp_dbcsr_sm_fm_multiply(op_ao(ii)%matrix,fm_Rmd_mos(iii)%matrix,&
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

! *****************************************************************************
!> \brief ...
!> \param issc_env ...
!> \param qs_env ...
!> \param iatom ...
!> \param error ...
! *****************************************************************************
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
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_dso, matrix_efg, &
                                                matrix_fc, matrix_pso
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: fc_psi0
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: dso_psi0, efg_psi0, pso_psi0
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
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
       CALL cp_dbcsr_set(matrix_fc(1)%matrix,0.0_dp,error=error)
       CALL build_fermi_contact_matrix(qs_env,matrix_fc,r_i,error)

       chk(1) = cp_dbcsr_checksum(matrix_fc(1)%matrix,error=error)

       IF(output_unit>0) THEN
          WRITE(output_unit,'(T2,A,E23.16)') 'ISSC| fermi_contact: CheckSum =',chk(1)
       ENDIF
    ENDIF
    !
    ! spin-orbit integral
    !IF(do_pso) THEN
    IF(.TRUE.) THEN ! for the moment we build it (regs)
       CALL cp_dbcsr_set(matrix_pso(1)%matrix,0.0_dp,error=error)
       CALL cp_dbcsr_set(matrix_pso(2)%matrix,0.0_dp,error=error)
       CALL cp_dbcsr_set(matrix_pso(3)%matrix,0.0_dp,error=error)
       CALL build_pso_matrix(qs_env,matrix_pso,r_i,error)

       chk(2) = cp_dbcsr_checksum(matrix_pso(1)%matrix,error=error)
       chk(3) = cp_dbcsr_checksum(matrix_pso(2)%matrix,error=error)
       chk(4) = cp_dbcsr_checksum(matrix_pso(3)%matrix,error=error)

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
       CALL cp_dbcsr_set(matrix_efg(1)%matrix,0.0_dp,error=error)
       CALL cp_dbcsr_set(matrix_efg(2)%matrix,0.0_dp,error=error)
       CALL cp_dbcsr_set(matrix_efg(3)%matrix,0.0_dp,error=error)
       CALL cp_dbcsr_set(matrix_efg(4)%matrix,0.0_dp,error=error)
       CALL cp_dbcsr_set(matrix_efg(5)%matrix,0.0_dp,error=error)
       CALL cp_dbcsr_set(matrix_efg(6)%matrix,0.0_dp,error=error)
       CALL build_efg_matrix(qs_env,matrix_efg,r_i,error)

       chk(5) = cp_dbcsr_checksum(matrix_efg(1)%matrix,error=error)
       chk(6) = cp_dbcsr_checksum(matrix_efg(2)%matrix,error=error)
       chk(7) = cp_dbcsr_checksum(matrix_efg(3)%matrix,error=error)
       chk(8) = cp_dbcsr_checksum(matrix_efg(4)%matrix,error=error)
       chk(9) = cp_dbcsr_checksum(matrix_efg(5)%matrix,error=error)
       chk(10)= cp_dbcsr_checksum(matrix_efg(6)%matrix,error=error)

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
       CALL cp_dbcsr_set(matrix_dso(1)%matrix,0.0_dp,error=error)
       CALL cp_dbcsr_set(matrix_dso(2)%matrix,0.0_dp,error=error)
       CALL cp_dbcsr_set(matrix_dso(3)%matrix,0.0_dp,error=error)
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
             CALL cp_dbcsr_sm_fm_multiply(matrix_efg(idir)%matrix,mo_coeff,&
                                       efg_psi0(ispin,idir)%matrix,ncol=nmo,&
                                       alpha=1.0_dp,error=error)
          ENDDO
       ENDIF
       !
       ! PSO
       IF(do_pso) THEN
          DO idir = 1,3
             CALL cp_dbcsr_sm_fm_multiply(matrix_pso(idir)%matrix,mo_coeff,&
                                      pso_psi0(ispin,idir)%matrix,ncol=nmo,&
                                      alpha=-1.0_dp,error=error)
          ENDDO
       ENDIF
       !
       ! FC
       IF(do_fc) THEN
          CALL cp_dbcsr_sm_fm_multiply(matrix_fc(1)%matrix,mo_coeff,&
                                    fc_psi0(ispin)%matrix,ncol=nmo,&
                                    alpha=1.0_dp,error=error)
       ENDIF
       !
       !>>> for debugging only
       IF(do_dso) THEN
          DO idir = 1,3
             CALL cp_dbcsr_sm_fm_multiply(matrix_dso(idir)%matrix,mo_coeff,&
                                       dso_psi0(ispin,idir)%matrix,ncol=nmo,&
                                       alpha=-1.0_dp,error=error)
          ENDDO
       ENDIF
       !<<< for debugging only
    ENDDO

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         "PRINT%PROGRAM_RUN_INFO",error=error)

    CALL timestop(handle)

  END SUBROUTINE issc_operators


! *****************************************************************************
!> \brief Calculate the Berry phase operator in the AO basis and 
!>         then the derivative of the Berry phase operator with respect to 
!>         the ground state wave function (see paper Putrino et al., JCP, 13, 7102) for the AOs;  
!>        afterwards multiply with the ground state MO coefficients 
!> 
!> \param polar_env ...
!> \param qs_env ...
!> \param error ...
!> \par History
!>      01.2013 created [SL]
!> \author SL
! *****************************************************************************

  SUBROUTINE polar_operators(polar_env,qs_env,error)

    TYPE(polar_env_type)                     :: polar_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'polar_operators', &
      routineP = moduleN//':'//routineN
    COMPLEX(KIND=dp), PARAMETER              :: one = (1.0_dp,0.0_dp) , &
                                                zero = (0.0_dp,0.0_dp)

    COMPLEX(DP)                              :: zdet, zdeta
    INTEGER                                  :: handle, i, idim, ispin, &
                                                istat, nao, nmo, nspins, &
                                                output_unit, tmp_dim, z
    LOGICAL                                  :: do_raman, failure
    REAL(dp)                                 :: kvec(3), maxocc
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_cfm_p_type), DIMENSION(:), &
      POINTER                                :: eigrmat
    TYPE(cp_cfm_p_type), DIMENSION(:, :), &
      POINTER                                :: inv_mat
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_dbcsr_type), POINTER             :: cosmat, sinmat
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: dBerry_psi0, mo_derivs, &
                                                op_fm_set, opvec
    TYPE(cp_fm_p_type), DIMENSION(:, :, :), &
      POINTER                                :: inv_work
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(section_vals_type), POINTER         :: lr_section

!  INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf, last_sgf
!  INTEGER, DIMENSION(:), POINTER           :: list_cubes, rbs

    failure = .FALSE.
    CALL timeset(routineN,handle)

    NULLIFY(dBerry_psi0,inv_work,mo_derivs,&
          op_fm_set,opvec,sinmat,cosmat)

    logger => cp_error_get_logger(error)
    lr_section => section_vals_get_subs_vals(qs_env%input,&
                  "PROPERTIES%LINRES",error=error)

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
                  extension=".linresLog",error=error)
    NULLIFY ( cell, dft_control,  mos,  matrix_s)
    CALL get_qs_env(qs_env=qs_env, &
                    cell=cell, &
                    dft_control=dft_control, &
                    para_env=para_env, &
                    mos=mos,&
                    matrix_s=matrix_s, &
                    error=error)

    nspins = dft_control%nspins


    CALL get_polar_env(polar_env=polar_env, &
                      do_raman=do_raman, &
                      dBerry_psi0=dBerry_psi0, &
                      mo_derivs=mo_derivs, &
                      error=error)
!SL calculate dipole berry phase
    IF(do_raman) THEN

      DO i= 1,3
          DO ispin = 1,nspins
             CALL cp_fm_set_all(dBerry_psi0(i,ispin)%matrix,0.0_dp,error=error)
             CALL cp_fm_set_all(mo_derivs(i,ispin)%matrix,0.0_dp,error=error)
          ENDDO
      END DO

    !
    !

    ! initialize all work matrices needed
    ALLOCATE ( opvec(2, dft_control%nspins ), STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE ( op_fm_set( 2, dft_control%nspins ), STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE ( eigrmat( dft_control%nspins ), STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE ( inv_mat(3, dft_control%nspins ), STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE ( inv_work(2,3, dft_control%nspins ), STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ! A bit to allocate for the wavefunction
    DO ispin = 1, dft_control%nspins
       NULLIFY(tmp_fm_struct,mo_coeff)
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff,nao=nao,nmo=nmo)
       CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nmo,&
            ncol_global=nmo,para_env=para_env,context=mo_coeff%matrix_struct%context,&
            error=error)
       DO i = 1, SIZE ( op_fm_set, 1 )
         CALL cp_fm_create (opvec(i,ispin)%matrix , mo_coeff%matrix_struct ,error=error)
          NULLIFY(op_fm_set(i,ispin)%matrix)
          CALL cp_fm_create (op_fm_set(i,ispin)%matrix , tmp_fm_struct ,error=error)
       END DO
       CALL cp_cfm_create ( eigrmat(ispin)%matrix, op_fm_set(1,ispin)%matrix%matrix_struct ,&
            error=error)
       CALL cp_fm_struct_release(tmp_fm_struct,error=error)
       DO i=1,3
         CALL cp_cfm_create ( inv_mat(i,ispin)%matrix, op_fm_set(1,ispin)%matrix%matrix_struct ,&
            error=error)
         CALL cp_fm_create ( inv_work(2,i,ispin)%matrix, op_fm_set(2,ispin)%matrix%matrix_struct ,&
             error=error)
         CALL cp_fm_create ( inv_work(1,i,ispin)%matrix, op_fm_set(1,ispin)%matrix%matrix_struct ,&
             error=error)

       END DO
    END DO

    !
    NULLIFY (cosmat, sinmat)
    ALLOCATE(cosmat, sinmat)
    CALL cp_dbcsr_init(cosmat, error=error)
    CALL cp_dbcsr_init(sinmat, error=error)
    CALL cp_dbcsr_copy(cosmat,matrix_s(1)%matrix,'COS MOM',error=error)
    CALL cp_dbcsr_copy(sinmat,matrix_s(1)%matrix,'SIN MOM',error=error)
    CALL cp_dbcsr_set(cosmat,0.0_dp,error=error)
    CALL cp_dbcsr_set(sinmat,0.0_dp,error=error)

    DO i=1,3
      kvec(:) = twopi*cell%h_inv(i,:)
      CALL cp_dbcsr_set(cosmat,0.0_dp,error=error)
      CALL cp_dbcsr_set(sinmat,0.0_dp,error=error)
      CALL build_berry_moment_matrix(qs_env,cosmat,sinmat,kvec,error)

      DO ispin=1, dft_control%nspins   ! spin
        CALL get_mo_set(mo_set=mos(ispin)%mo_set,nao=nao,mo_coeff=mo_coeff,nmo=nmo)

        CALL cp_dbcsr_sm_fm_multiply(cosmat, mo_coeff, opvec(1,ispin)%matrix, ncol=nmo, error=error)
        CALL cp_gemm("T","N",nmo,nmo,nao,1.0_dp,mo_coeff,opvec(1,ispin)%matrix,0.0_dp,&
           op_fm_set(1,ispin)%matrix,error=error)
        CALL cp_dbcsr_sm_fm_multiply(sinmat, mo_coeff, opvec(2,ispin)%matrix, ncol=nmo, error=error)
        CALL cp_gemm("T","N",nmo,nmo,nao,1.0_dp,mo_coeff,opvec(2,ispin)%matrix,0.0_dp,&
               op_fm_set(2,ispin)%matrix,error=error)

      ENDDO

        !second step invert C^T S_berry C
        zdet=one
        DO ispin = 1, dft_control%nspins
          CALL cp_cfm_get_info(eigrmat(ispin)%matrix,ncol_local=tmp_dim,error=error)
          DO idim=1,tmp_dim
             eigrmat(ispin)%matrix%local_data(:,idim) = &
                  CMPLX (op_fm_set(1,ispin)%matrix%local_data(:,idim), &
                  -op_fm_set(2,ispin)%matrix%local_data(:,idim),dp)
          END DO
          CALL cp_cfm_set_all (inv_mat(i,ispin)%matrix,zero,one,error)
          CALL cp_cfm_solve ( eigrmat(ispin)%matrix,inv_mat(i,ispin)%matrix, zdeta,error )
        END DO
 !
 !       
 !compute the derivative and add the result to mo_derivatives
        DO ispin=1,dft_control%nspins
          CALL cp_cfm_get_info(eigrmat(ispin)%matrix,ncol_local=tmp_dim,error=error)
          CALL get_mo_set(mo_set=mos(ispin)%mo_set,nao=nao,nmo=nmo,maxocc=maxocc)
          DO z=1,tmp_dim
             inv_work(1,i,ispin)%matrix%local_data(:,z)=REAL(inv_mat(i,ispin)%matrix%local_data(:,z),dp)
             inv_work(2,i,ispin)%matrix%local_data(:,z)=AIMAG(inv_mat(i,ispin)%matrix%local_data(:,z))
          END DO
          CALL cp_gemm("N","N",nao,nmo,nmo,1.0_dp,opvec(1,ispin)%matrix,inv_work(2,i,ispin)%matrix,&
              0.0_dp,mo_derivs(i,ispin)%matrix,error)
          CALL cp_gemm("N","N",nao,nmo,nmo,-1.0_dp,opvec(2,ispin)%matrix,inv_work(1,i,ispin)%matrix,&
              1.0_dp,mo_derivs(i,ispin)%matrix,error)

          CALL cp_gemm("N","N",nao,nmo,nmo,-1.0_dp,opvec(1,ispin)%matrix,inv_work(2,i,ispin)%matrix,&
               0.0_dp,dBerry_psi0(i,ispin)%matrix,error)
          CALL cp_gemm("N","N",nao,nmo,nmo,1.0_dp,opvec(2,ispin)%matrix,inv_work(1,i,ispin)%matrix,&
               1.0_dp,dBerry_psi0(i,ispin)%matrix,error)


        END DO
    END DO !x/y/z-direction
!SL we omit here the multiplication with hmat (this scaling back done at the end of the response calc)

     DO ispin=1,dft_control%nspins
          CALL cp_cfm_release(eigrmat(ispin)%matrix,error=error)
        DO i=1,3
          CALL cp_cfm_release(inv_mat(i,ispin)%matrix,error=error)
          CALL cp_fm_release(inv_work(1,i,ispin)%matrix,error=error)
          CALL cp_fm_release(inv_work(2,i,ispin)%matrix,error=error)
        END DO
          CALL cp_fm_release(opvec(1,ispin)%matrix,error=error)
          CALL cp_fm_release(opvec(2,ispin)%matrix,error=error)
          CALL cp_fm_release(op_fm_set(2,ispin)%matrix,error=error)
          CALL cp_fm_release(op_fm_set(1,ispin)%matrix,error=error)
    END DO

    CALL cp_dbcsr_deallocate_matrix ( cosmat, error=error )
    CALL cp_dbcsr_deallocate_matrix ( sinmat, error=error )




    DEALLOCATE ( inv_work, STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE ( opvec, STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE ( op_fm_set, STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE ( eigrmat, STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE ( inv_mat, STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)



    ENDIF  ! raman
    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         "PRINT%PROGRAM_RUN_INFO",error=error)

    CALL timestop(handle)

  END SUBROUTINE polar_operators

! *****************************************************************************




! *****************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
!> \param c ...
!> \retval factor ...
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
!> \brief ...
!> \param ii ...
!> \param iii ...
!> \retval i ...
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
!> \brief ...
!> \param i1 ...
!> \param i2 ...
!> \param i3 ...
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
!> \brief ...
!> \param i1 ...
!> \param i2 ...
!> \param i3 ...
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
!> \param matrix ...
!> \param ra ...
!> \param rc ...
!> \param cell ...
!> \param ixyz ...
!> \param error ...
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

    INTEGER :: handle, icol_local, irow_local, m, mypcol, myprow, n, &
      ncol_block, ncol_global, ncol_local, npcol, nprow, nrow_block, &
      nrow_global, nrow_local
    REAL(KIND=dp)                            :: dist(3), rra(3), rrc(3)
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: a
#if defined(__SCALAPACK)
    INTEGER                                  :: icol_global, irow_global
    INTEGER, EXTERNAL                        :: indxl2g
#endif

    CALL timeset(routineN,handle)

    myprow = matrix%matrix_struct%context%mepos(1)
    mypcol = matrix%matrix_struct%context%mepos(2)
    nprow = matrix%matrix_struct%context%num_pe(1)
    npcol = matrix%matrix_struct%context%num_pe(2)

    nrow_block  = matrix%matrix_struct%nrow_block
    ncol_block  = matrix%matrix_struct%ncol_block
    nrow_global = matrix%matrix_struct%nrow_global
    ncol_global = matrix%matrix_struct%ncol_global
    nrow_local  = matrix%matrix_struct%nrow_locals(myprow)
    ncol_local  = matrix%matrix_struct%ncol_locals(mypcol)

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

