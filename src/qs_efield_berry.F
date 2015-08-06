!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculates the energy contribution and the mo_derivative of
!>        a static periodic electric field
!> \par History
!>      none
!> \author fschiff (06.2010)
! *****************************************************************************
MODULE qs_efield_berry
  USE ai_moments,                      ONLY: cossin
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                             gto_basis_set_type
  USE block_p_types,                   ONLY: block_p_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_cfm_basic_linalg,             ONLY: cp_cfm_add_fm,&
                                             cp_cfm_solve
  USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                             cp_cfm_p_type,&
                                             cp_cfm_release,&
                                             cp_cfm_set_all
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_copy,&
                                             cp_dbcsr_deallocate_matrix_set,&
                                             cp_dbcsr_get_block_p,&
                                             cp_dbcsr_init,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_set,&
                                             cp_dbcsr_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             copy_fm_to_dbcsr,&
                                             cp_dbcsr_plus_fm_fm_t,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi,&
                                             twopi
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             get_qs_kind_set,&
                                             qs_kind_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_moments,                      ONLY: build_berry_moment_matrix
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE qs_period_efield_types,          ONLY: efield_berry_type,&
                                             init_efield_matrices,&
                                             set_efield_matrices
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_pair_force
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_efield_berry'

  ! *** Public subroutines ***

  PUBLIC :: qs_efield_berry_phase

! *****************************************************************************

  CONTAINS

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param just_energy ...
!> \param calculate_forces ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE qs_efield_berry_phase(qs_env,just_energy,calculate_forces,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: just_energy, calculate_forces
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_efield_berry_phase', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: s_mstruct_changed
    TYPE(dft_control_type), POINTER          :: dft_control

    CALL timeset(routineN,handle)

    NULLIFY(dft_control)
    CALL get_qs_env(qs_env,s_mstruct_changed=s_mstruct_changed,&
                    dft_control=dft_control,error=error)

    IF(dft_control%apply_period_efield)THEN
       IF(s_mstruct_changed) CALL qs_efield_integrals(qs_env,error)
       IF(dft_control%period_efield%displacement_field) THEN
          CALL qs_dispfield_derivatives(qs_env,just_energy,calculate_forces,error)
       ELSE
          CALL qs_efield_derivatives(qs_env,just_energy,calculate_forces,error)
       END IF
    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_efield_berry_phase

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE qs_efield_integrals (qs_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_efield_integrals', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, stat
    LOGICAL                                  :: failure
    REAL(dp), DIMENSION(3)                   :: kvec
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: cosmat, matrix_s, sinmat
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(efield_berry_type), POINTER         :: efield

    CALL timeset(routineN,handle)
    failure = .FALSE.
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,error=error)
    NULLIFY (matrix_s)
    CALL get_qs_env(qs_env=qs_env,efield=efield,cell=cell,matrix_s=matrix_s,error=error)
    CALL init_efield_matrices(efield,error)
    ALLOCATE(cosmat(3),sinmat(3),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO i=1,3
       ALLOCATE(cosmat(i)%matrix,sinmat(i)%matrix,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL cp_dbcsr_init(cosmat(i)%matrix,error=error)
       CALL cp_dbcsr_init(sinmat(i)%matrix, error=error)

       CALL cp_dbcsr_copy(cosmat(i)%matrix,matrix_s(1)%matrix,'COS MAT',error=error)
       CALL cp_dbcsr_copy(sinmat(i)%matrix,matrix_s(1)%matrix,'SIN MAT',error=error)
       CALL cp_dbcsr_set(cosmat(i)%matrix,0.0_dp,error=error)
       CALL cp_dbcsr_set(sinmat(i)%matrix,0.0_dp,error=error)

       kvec(:) = twopi*cell%h_inv(i,:)
       CALL build_berry_moment_matrix(qs_env,cosmat(i)%matrix,sinmat(i)%matrix,kvec,error)
    END DO
    CALL set_efield_matrices(efield=efield,cosmat=cosmat,sinmat=sinmat)
    CALL set_qs_env(qs_env=qs_env,efield=efield,error=error)
    CALL timestop(handle)

  END SUBROUTINE qs_efield_integrals

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param just_energy ...
!> \param calculate_forces ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE qs_efield_derivatives(qs_env,just_energy,calculate_forces,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: just_energy, calculate_forces
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_efield_derivatives', &
      routineP = moduleN//':'//routineN
    COMPLEX(KIND=dp), PARAMETER              :: one = (1.0_dp,0.0_dp), &
                                                zero = (0.0_dp,0.0_dp), &
                                                zone = (0.0_dp,1.0_dp)

    COMPLEX(dp)                              :: zdet, zdeta, zi(3)
    INTEGER :: atom_a, atom_b, handle, i, ia, iatom, icol, idir, ikind, &
      inode, irow, iset, ispin, istat, j, jatom, jkind, jset, ldab, ldsa, &
      ldsb, lsab, n1, n2, nao, natom, ncoa, ncob, nkind, nmo, nseta, nsetb, &
      sgfa, sgfb, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure, found, uniform, &
                                                use_virial
    REAL(dp) :: charge, ci(3), dab, dd, ener_field, f0, fab, fieldpol(3), &
      focc, fpolvec(3), hmat(3,3), occ, qi(3), ti(3)
    REAL(dp), DIMENSION(3)                   :: forcea, forceb, kvec, ra, &
                                                rab, rb, ria
    REAL(dp), DIMENSION(:, :), POINTER       :: cosab, iblock, rblock, sinab, &
                                                work
    REAL(dp), DIMENSION(:, :, :), POINTER    :: dcosab, dsinab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                zeta, zetb
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(block_p_type), DIMENSION(3, 2)      :: dcost, dsint
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_cfm_p_type), DIMENSION(:), &
      POINTER                                :: eigrmat, inv_mat
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s, mo_derivs
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: tempmat
    TYPE(cp_dbcsr_type), POINTER             :: cosmat, mo_coeff_b, sinmat
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_coeff_tmp, mo_derivs_tmp
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: inv_work, op_fm_set, opvec
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(efield_berry_type), POINTER         :: efield
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_kind_type), POINTER              :: qs_kind
    TYPE(virial_type), POINTER               :: virial

    CALL timeset(routineN,handle)

    NULLIFY(dft_control,cell,particle_set)
    CALL get_qs_env(qs_env,dft_control=dft_control,cell=cell,&
         particle_set=particle_set,virial=virial,error=error)
    NULLIFY(qs_kind_set,efield,para_env,sab_orb)
    CALL get_qs_env(qs_env=qs_env,qs_kind_set=qs_kind_set,&
         efield=efield,energy=energy,para_env=para_env,sab_orb=sab_orb,error=error)

    ! calculate stress only if forces requested also
    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)
    use_virial = use_virial .AND. calculate_forces
    ! disable stress calculation
    IF (use_virial) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="Stress tensor for periodic E-field not implemented", &
            error=error, error_level=cp_failure_level)
    END IF

    fieldpol=dft_control%period_efield%polarisation
    fieldpol=fieldpol/SQRT(DOT_PRODUCT(fieldpol,fieldpol))
    fieldpol=-fieldpol*dft_control%period_efield%strength
    hmat = cell%hmat(:,:)/twopi
    DO idir=1,3
       fpolvec(idir)=fieldpol(1)*hmat(1,idir)+fieldpol(2)*hmat(2,idir)+fieldpol(3)*hmat(3,idir)
    END DO

    ! nuclear contribution
    natom = SIZE(particle_set)
    IF (calculate_forces) THEN
       CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,force=force,error=error)
       ALLOCATE (atom_of_kind(natom),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL get_atomic_kind_set(atomic_kind_set,atom_of_kind=atom_of_kind)
    END IF
    zi(:) = CMPLX ( 1._dp, 0._dp, dp )
    DO ia = 1,natom
       CALL get_atomic_kind(particle_set(ia)%atomic_kind,kind_number=ikind)
       CALL get_qs_kind(qs_kind_set(ikind),core_charge=charge,error=error)
       ria = particle_set(ia)%r
       ria = pbc(ria,cell)
       DO idir = 1, 3
          kvec(:) = twopi*cell%h_inv(idir,:)
          dd = SUM ( kvec(:) * ria(:) )
          zdeta = CMPLX(COS(dd),SIN(dd),KIND=dp)**charge
          zi(idir) = zi(idir) * zdeta
       END DO
       IF(calculate_forces)THEN
          IF(para_env%mepos==0) THEN
             iatom = atom_of_kind(ia)
             forcea(:) = fieldpol(:)*charge
             force(ikind)%efield(:,iatom)=force(ikind)%efield(:,iatom)+forcea(:)
          END IF
       END IF
       IF (use_virial) THEN
          IF(para_env%mepos==0) &
             CALL virial_pair_force ( virial%pv_virial, 1.0_dp, forcea, ria, error)
       END IF
    END DO
    qi = AIMAG(LOG(zi))

    ! check uniform occupation
    NULLIFY(mos)
    CALL get_qs_env(qs_env=qs_env,mos=mos,error=error)
    DO ispin = 1,dft_control%nspins
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,maxocc=occ,uniform_occupation=uniform)
       IF (.NOT.uniform) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Berry phase moments for non uniform MOs' occupation numbers not implemented", &
               error=error, error_level=cp_failure_level)
       END IF
    END DO

    NULLIFY (mo_derivs)
    CALL get_qs_env(qs_env=qs_env,mo_derivs=mo_derivs,error=error)
    ! initialize all work matrices needed
    ALLOCATE ( op_fm_set( 2, dft_control%nspins ), STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE ( opvec( 2, dft_control%nspins ), STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE ( eigrmat( dft_control%nspins ), STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE ( inv_mat( dft_control%nspins ), STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE ( inv_work(2, dft_control%nspins ), STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE ( mo_derivs_tmp(SIZE(mo_derivs)), STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE ( mo_coeff_tmp(SIZE(mo_derivs)), STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ! Allocate temp matrices for the wavefunction derivatives
    DO ispin = 1, dft_control%nspins
       NULLIFY(tmp_fm_struct,mo_coeff)
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff,nao=nao,nmo=nmo)
       CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nmo,&
            ncol_global=nmo,para_env=para_env,context=mo_coeff%matrix_struct%context,error=error)
       CALL cp_fm_create (mo_derivs_tmp(ispin)%matrix,mo_coeff%matrix_struct,error=error)
       CALL cp_fm_create (mo_coeff_tmp(ispin)%matrix,mo_coeff%matrix_struct,error=error)
       CALL copy_dbcsr_to_fm(mo_derivs(ispin)%matrix,mo_derivs_tmp(ispin)%matrix,error=error)
       DO i = 1, SIZE(op_fm_set,1)
          CALL cp_fm_create (opvec(i,ispin)%matrix,mo_coeff%matrix_struct,error=error)
          NULLIFY(op_fm_set(i,ispin)%matrix)
          CALL cp_fm_create (op_fm_set(i,ispin)%matrix,tmp_fm_struct,error=error)
          CALL cp_fm_create (inv_work(i,ispin)%matrix,op_fm_set(i,ispin)%matrix%matrix_struct,error=error)
       END DO
       CALL cp_cfm_create (eigrmat(ispin)%matrix,op_fm_set(1,ispin)%matrix%matrix_struct,error=error)
       CALL cp_cfm_create (inv_mat(ispin)%matrix,op_fm_set(1,ispin)%matrix%matrix_struct,error=error)
       CALL cp_fm_struct_release(tmp_fm_struct,error=error)
    END DO
    ! temp matrices for force calculation
    IF(calculate_forces)THEN
       NULLIFY (matrix_s)
       CALL get_qs_env(qs_env=qs_env,matrix_s=matrix_s,error=error)
       ALLOCATE(tempmat(2,dft_control%nspins),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO ispin=1,dft_control%nspins
          ALLOCATE(tempmat(1,ispin)%matrix,tempmat(2,ispin)%matrix,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          CALL cp_dbcsr_init(tempmat(1,ispin)%matrix,error=error)
          CALL cp_dbcsr_init(tempmat(2,ispin)%matrix, error=error)
          CALL cp_dbcsr_copy(tempmat(1,ispin)%matrix,matrix_s(1)%matrix,'TEMPMAT',error=error)
          CALL cp_dbcsr_copy(tempmat(2,ispin)%matrix,matrix_s(1)%matrix,'TEMPMAT',error=error)
          CALL cp_dbcsr_set(tempmat(1,ispin)%matrix,0.0_dp,error=error)
          CALL cp_dbcsr_set(tempmat(2,ispin)%matrix,0.0_dp,error=error)
       END DO
       ! integration
       CALL get_qs_kind_set(qs_kind_set,maxco=ldab,maxsgf=lsab,error=error)
       ALLOCATE(cosab(ldab,ldab),sinab(ldab,ldab),work(ldab,ldab),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(dcosab(ldab,ldab,3),dsinab(ldab,ldab,3),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       lsab = MAX(ldab,lsab)
       DO i=1,3
          ALLOCATE(dcost(i,1)%block(lsab,lsab),dsint(i,1)%block(lsab,lsab),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(dcost(i,2)%block(lsab,lsab),dsint(i,2)%block(lsab,lsab),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END DO
    END IF

    !Start the MO derivative calculation
    !loop over all cell vectors
    DO idir=1,3
       ci(idir) = 0.0_dp
       zi(idir) = zero
       IF(ABS(fpolvec(idir)).GT.1.0E-12_dp)THEN
          cosmat => efield%cosmat(idir)%matrix
          sinmat => efield%sinmat(idir)%matrix
          !evaluate the expression needed for the derivative (S_berry * C  and [C^T S_berry C]^-1)
          !first step S_berry * C  and C^T S_berry C
          DO ispin=1, dft_control%nspins   ! spin
             IF(mos(ispin)%mo_set%use_mo_coeff_b)THEN
                CALL get_mo_set(mo_set=mos(ispin)%mo_set,nao=nao,mo_coeff_b=mo_coeff_b,nmo=nmo)
                CALL copy_dbcsr_to_fm(mo_coeff_b,mo_coeff_tmp(ispin)%matrix,error=error)
             ELSE
                CALL get_mo_set(mo_set=mos(ispin)%mo_set,nao=nao,mo_coeff=mo_coeff_tmp(ispin)%matrix,nmo=nmo)
             END IF
             CALL cp_dbcsr_sm_fm_multiply(cosmat,mo_coeff_tmp(ispin)%matrix,opvec(1,ispin)%matrix,ncol=nmo,error=error)
             CALL cp_gemm("T","N",nmo,nmo,nao,1.0_dp,mo_coeff_tmp(ispin)%matrix,opvec(1,ispin)%matrix,0.0_dp,&
                  op_fm_set(1,ispin)%matrix,error=error)
             CALL cp_dbcsr_sm_fm_multiply(sinmat,mo_coeff_tmp(ispin)%matrix,opvec(2,ispin)%matrix,ncol=nmo,error=error)
             CALL cp_gemm("T","N",nmo,nmo,nao,1.0_dp,mo_coeff_tmp(ispin)%matrix,opvec(2,ispin)%matrix,0.0_dp,&
                  op_fm_set(2,ispin)%matrix,error=error)
          ENDDO
          !second step invert C^T S_berry C
          zdet = one
          DO ispin = 1, dft_control%nspins
             CALL cp_cfm_add_fm(zero,eigrmat(ispin)%matrix,one,op_fm_set(1,ispin)%matrix,error)
             CALL cp_cfm_add_fm(one,eigrmat(ispin)%matrix,-zone,op_fm_set(2,ispin)%matrix,error)
             CALL cp_cfm_set_all(inv_mat(ispin)%matrix,zero,one,error)
             CALL cp_cfm_solve(eigrmat(ispin)%matrix,inv_mat(ispin)%matrix,zdeta,error)
             zdet = zdet*zdeta
          END DO
          zi(idir) = zdet**occ
          ci(idir) = AIMAG(LOG(zdet**occ))

          IF(.NOT.just_energy) THEN
             !compute the orbital derivative
             focc = fpolvec(idir)
             DO ispin=1,dft_control%nspins
                inv_work(1,ispin)%matrix%local_data(:,:)=REAL(inv_mat(ispin)%matrix%local_data(:,:),dp)
                inv_work(2,ispin)%matrix%local_data(:,:)=AIMAG(inv_mat(ispin)%matrix%local_data(:,:))
                CALL get_mo_set(mo_set=mos(ispin)%mo_set,nao=nao,nmo=nmo)
                CALL cp_gemm("N","N",nao,nmo,nmo,focc,opvec(1,ispin)%matrix,inv_work(2,ispin)%matrix,&
                     1.0_dp,mo_derivs_tmp(ispin)%matrix,error)
                CALL cp_gemm("N","N",nao,nmo,nmo,-focc,opvec(2,ispin)%matrix,inv_work(1,ispin)%matrix,&
                     1.0_dp,mo_derivs_tmp(ispin)%matrix,error)
             END DO
          END IF

          !compute nuclear forces
          IF(calculate_forces)THEN
             nkind = SIZE(qs_kind_set)
             natom = SIZE(particle_set)
             kvec(:) = twopi*cell%h_inv(idir,:)

             ! calculate: C [C^T S_berry C]^(-1) C^T
             ! Store this matrix in DBCSR form (only S overlap blocks)
             DO ispin=1,dft_control%nspins
                CALL cp_dbcsr_set(tempmat(1,ispin)%matrix,0.0_dp,error=error)
                CALL cp_dbcsr_set(tempmat(2,ispin)%matrix,0.0_dp,error=error)
                CALL get_mo_set(mo_set=mos(ispin)%mo_set,nao=nao,nmo=nmo)
                CALL cp_gemm("N","N",nao,nmo,nmo,1.0_dp,mo_coeff_tmp(ispin)%matrix,inv_work(1,ispin)%matrix,0.0_dp,&
                     opvec(1,ispin)%matrix,error=error)
                CALL cp_gemm("N","N",nao,nmo,nmo,1.0_dp,mo_coeff_tmp(ispin)%matrix,inv_work(2,ispin)%matrix,0.0_dp,&
                     opvec(2,ispin)%matrix,error=error)
                CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=tempmat(1,ispin)%matrix,&
                     matrix_v=opvec(1,ispin)%matrix,matrix_g=mo_coeff_tmp(ispin)%matrix,ncol=nmo,error=error)
                CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=tempmat(2,ispin)%matrix,&
                     matrix_v=opvec(2,ispin)%matrix,matrix_g=mo_coeff_tmp(ispin)%matrix,ncol=nmo,error=error)
             END DO

             ! Calculation of derivative integrals (da|eikr|b) and (a|eikr|db)
             ALLOCATE (basis_set_list(nkind),STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
             DO ikind=1,nkind
                qs_kind => qs_kind_set(ikind)
                CALL get_qs_kind(qs_kind=qs_kind,basis_set=basis_set_a,error=error)
                IF (ASSOCIATED(basis_set_a)) THEN
                   basis_set_list(ikind)%gto_basis_set => basis_set_a
                ELSE
                   NULLIFY(basis_set_list(ikind)%gto_basis_set)
                END IF
             END DO
             !
             CALL neighbor_list_iterator_create(nl_iterator,sab_orb)
             DO WHILE (neighbor_list_iterate(nl_iterator)==0)
                CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,inode=inode,&
                                       iatom=iatom,jatom=jatom,r=rab)
                basis_set_a => basis_set_list(ikind)%gto_basis_set
                IF (.NOT.ASSOCIATED(basis_set_a)) CYCLE
                basis_set_b => basis_set_list(jkind)%gto_basis_set
                IF (.NOT.ASSOCIATED(basis_set_b)) CYCLE
                ! basis ikind
                first_sgfa   =>  basis_set_a%first_sgf
                la_max       =>  basis_set_a%lmax
                la_min       =>  basis_set_a%lmin
                npgfa        =>  basis_set_a%npgf
                nseta        =   basis_set_a%nset
                nsgfa        =>  basis_set_a%nsgf_set
                rpgfa        =>  basis_set_a%pgf_radius
                set_radius_a =>  basis_set_a%set_radius
                sphi_a       =>  basis_set_a%sphi
                zeta         =>  basis_set_a%zet
                ! basis jkind
                first_sgfb   =>  basis_set_b%first_sgf
                lb_max       =>  basis_set_b%lmax
                lb_min       =>  basis_set_b%lmin
                npgfb        =>  basis_set_b%npgf
                nsetb        =   basis_set_b%nset
                nsgfb        =>  basis_set_b%nsgf_set
                rpgfb        =>  basis_set_b%pgf_radius
                set_radius_b =>  basis_set_b%set_radius
                sphi_b       =>  basis_set_b%sphi
                zetb         =>  basis_set_b%zet

                atom_a = atom_of_kind(iatom)
                atom_b = atom_of_kind(jatom)

                ldsa = SIZE(sphi_a,1)
                ldsb = SIZE(sphi_b,1)
                ra(:) = pbc(particle_set(iatom)%r(:),cell)
                rb(:) = ra+rab
                dab = SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))

                IF (iatom <= jatom) THEN
                   irow = iatom
                   icol = jatom
                ELSE
                   irow = jatom
                   icol = iatom
                END IF

                IF (iatom == jatom .AND. dab < 1.e-10_dp) THEN
                   fab = 1.0_dp * occ
                ELSE
                   fab = 2.0_dp * occ
                END IF

                DO i=1,3
                   dcost(i,1)%block = 0.0_dp
                   dsint(i,1)%block = 0.0_dp
                   dcost(i,2)%block = 0.0_dp
                   dsint(i,2)%block = 0.0_dp
                END DO

                DO iset=1,nseta
                   ncoa = npgfa(iset)*ncoset(la_max(iset))
                   sgfa = first_sgfa(1,iset)
                   DO jset=1,nsetb
                      IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE
                      ncob = npgfb(jset)*ncoset(lb_max(jset))
                      sgfb = first_sgfb(1,jset)
                      ! Calculate the primitive integrals (da|b)
                      CALL cossin(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                           lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                           ra,rb,kvec,cosab,sinab,dcosab,dsinab)
                      DO i=1,3
                         CALL contract_all(dcost(i,1)%block,dsint(i,1)%block,&
                              ncoa,nsgfa(iset),sgfa,sphi_a,ldsa,&
                              ncob,nsgfb(jset),sgfb,sphi_b,ldsb,&
                              dcosab(:,:,i),dsinab(:,:,i),ldab,work,ldab)
                      END DO
                      ! Calculate the primitive integrals (a|db)
                      CALL cossin(lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                           la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                           rb,ra,kvec,cosab,sinab,dcosab,dsinab)
                      DO i=1,3
                         dcosab(1:ncoa,1:ncob,i) = TRANSPOSE(dcosab(1:ncob,1:ncoa,i))
                         dsinab(1:ncoa,1:ncob,i) = TRANSPOSE(dsinab(1:ncob,1:ncoa,i))
                         CALL contract_all(dcost(i,2)%block,dsint(i,2)%block,&
                              ncoa,nsgfa(iset),sgfa,sphi_a,ldsa,&
                              ncob,nsgfb(jset),sgfb,sphi_b,ldsb,&
                              dcosab(:,:,i),dsinab(:,:,i),ldab,work,ldab)
                      END DO
                   END DO
                END DO
                forcea = 0.0_dp
                forceb = 0.0_dp
                DO ispin=1,dft_control%nspins
                   NULLIFY (rblock,iblock)
                   CALL cp_dbcsr_get_block_p(matrix=tempmat(1,ispin)%matrix,&
                        row=irow,col=icol,BLOCK=rblock,found=found)
                   CPPostcondition(found,cp_failure_level,routineP,error,failure)
                   CALL cp_dbcsr_get_block_p(matrix=tempmat(2,ispin)%matrix,&
                        row=irow,col=icol,BLOCK=iblock,found=found)
                   CPPostcondition(found,cp_failure_level,routineP,error,failure)
                   n1 = SIZE(rblock,1)
                   n2 = SIZE(rblock,2)
                   CPPostcondition(SIZE(iblock,1)==n1,cp_failure_level,routineP,error,failure)
                   CPPostcondition(SIZE(iblock,2)==n2,cp_failure_level,routineP,error,failure)
                   CPPostcondition(lsab>=n1,cp_failure_level,routineP,error,failure)
                   CPPostcondition(lsab>=n2,cp_failure_level,routineP,error,failure)
                   IF (iatom <= jatom) THEN
                      DO i=1,3
                         forcea(i) = forcea(i) + SUM(rblock(1:n1,1:n2)*dsint(i,1)%block(1:n1,1:n2)) &
                                               - SUM(iblock(1:n1,1:n2)*dcost(i,1)%block(1:n1,1:n2))
                         forceb(i) = forceb(i) + SUM(rblock(1:n1,1:n2)*dsint(i,2)%block(1:n1,1:n2)) &
                                               - SUM(iblock(1:n1,1:n2)*dcost(i,2)%block(1:n1,1:n2))
                      END DO
                   ELSE
                      DO i=1,3
                         forcea(i) = forcea(i) + SUM(TRANSPOSE(rblock(1:n1,1:n2))*dsint(i,1)%block(1:n2,1:n1)) &
                                               - SUM(TRANSPOSE(iblock(1:n1,1:n2))*dcost(i,1)%block(1:n2,1:n1))
                         forceb(i) = forceb(i) + SUM(TRANSPOSE(rblock(1:n1,1:n2))*dsint(i,2)%block(1:n2,1:n1)) &
                                               - SUM(TRANSPOSE(iblock(1:n1,1:n2))*dcost(i,2)%block(1:n2,1:n1))
                      END DO
                   END IF
                END DO
                force(ikind)%efield(1:3,atom_a)=force(ikind)%efield(1:3,atom_a)-fab*fpolvec(idir)*forcea(1:3)
                force(jkind)%efield(1:3,atom_b)=force(jkind)%efield(1:3,atom_b)-fab*fpolvec(idir)*forceb(1:3)
                IF (use_virial) THEN
                   f0 = -fab*fpolvec(idir)
                   CALL virial_pair_force ( virial%pv_virial, f0, forcea, ra, error)
                   CALL virial_pair_force ( virial%pv_virial, f0, forceb, rb, error)
                END IF

             END DO
             CALL neighbor_list_iterator_release(nl_iterator)
             DEALLOCATE (basis_set_list,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

          END IF
       END IF
    END DO

    ! Energy
    ener_field=0.0_dp
    DO idir=1,3
       ! make sure the total normalized polarization is within [-1:1]
       IF( qi(idir)+ci(idir) > pi ) ci(idir) = ci(idir) - twopi
       IF( qi(idir)+ci(idir) < -pi ) ci(idir) = ci(idir) + twopi
       ener_field = ener_field + fpolvec(idir) * (qi(idir)+ci(idir))
    END DO
    energy%efield=ener_field

    ! check for smoothness of energy surface
    IF(efield%field_energy == 0.0_dp) THEN
       efield%field_energy = ener_field
    ELSEIF( ABS(efield%field_energy-ener_field) > pi*ABS(SUM(fpolvec))) THEN
       CALL cp_error_message(cp_note_level, routineP, &
        "Large change of e-field energy. Non-smooth energy surface", error)
    ELSE
       efield%field_energy = ener_field
    END IF

    IF(.NOT.just_energy) THEN
       ! Add the result to mo_derivativs
       DO ispin=1, dft_control%nspins
          CALL copy_fm_to_dbcsr(mo_derivs_tmp(ispin)%matrix,mo_derivs(ispin)%matrix,error=error)
       END DO
       IF (use_virial) THEN
          ti = 0.0_dp
          DO i=1,3
             DO j=1,3
                ti(j) = ti(j) + hmat(j,i)*(qi(i)+ci(i))
             END DO
          END DO
          DO i=1,3
             DO j=1,3
                virial%pv_virial(i,j) = virial%pv_virial(i,j) - fieldpol(i)*ti(j)
             END DO
          END DO
       END IF
    END IF

    DO ispin = 1, dft_control%nspins
       CALL cp_cfm_release(eigrmat(ispin)%matrix,error=error)
       CALL cp_cfm_release(inv_mat(ispin)%matrix,error=error)
       CALL cp_fm_release(mo_derivs_tmp(ispin)%matrix,error=error)
       CALL cp_fm_release(mo_coeff_tmp(ispin)%matrix,error=error)
       DO i = 1, SIZE ( op_fm_set, 1 )
          CALL cp_fm_release(opvec(i,ispin)%matrix,error=error)
          CALL cp_fm_release(op_fm_set(i,ispin)%matrix,error=error)
          CALL cp_fm_release(inv_work(i,ispin)%matrix,error=error)
       END DO
    END DO
    DEALLOCATE(inv_mat,inv_work,op_fm_set,opvec,eigrmat,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(mo_coeff_tmp,mo_derivs_tmp,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    IF(calculate_forces)THEN
       DO ikind=1,SIZE(atomic_kind_set)
          CALL mp_sum(force(ikind)%efield,para_env%group)
       END DO
       DEALLOCATE(atom_of_kind,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(cosab,sinab,work,dcosab,dsinab,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DO i=1,3
          DEALLOCATE(dcost(i,1)%block,dsint(i,1)%block,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(dcost(i,2)%block,dsint(i,2)%block,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END DO
       CALL cp_dbcsr_deallocate_matrix_set(tempmat,error)
    END IF
    CALL timestop(handle)

  END SUBROUTINE qs_efield_derivatives

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param just_energy ...
!> \param calculate_forces ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE qs_dispfield_derivatives(qs_env,just_energy,calculate_forces,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: just_energy, calculate_forces
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_dispfield_derivatives', &
      routineP = moduleN//':'//routineN
    COMPLEX(KIND=dp), PARAMETER              :: one = (1.0_dp,0.0_dp), &
                                                zero = (0.0_dp,0.0_dp), &
                                                zone = (0.0_dp,1.0_dp)

    COMPLEX(dp)                              :: zdet, zdeta, zi(3)
    INTEGER :: handle, i, ia, iatom, icol, idir, ikind, inode, irow, iset, &
      ispin, istat, jatom, jkind, jset, ldab, ldsa, ldsb, lsab, n1, n2, nao, &
      natom, ncoa, ncob, nkind, nmo, nseta, nsetb, sgfa, sgfb, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure, found, uniform, &
                                                use_virial
    REAL(dp) :: charge, ci(3), dab, dd, di(3), ener_field, fab, fieldpol(3), &
      focc, hmat(3,3), occ, omega, qi(3), rlog(3), zlog(3)
    REAL(dp), DIMENSION(3)                   :: forcea, forceb, kvec, ra, &
                                                rab, rb, ria
    REAL(dp), DIMENSION(:, :), POINTER       :: cosab, iblock, rblock, sinab, &
                                                work
    REAL(dp), DIMENSION(:, :, :), POINTER    :: dcosab, dsinab, force_tmp
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                zeta, zetb
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(block_p_type), DIMENSION(3, 2)      :: dcost, dsint
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_cfm_p_type), DIMENSION(:), &
      POINTER                                :: eigrmat, inv_mat
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s, mo_derivs
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: tempmat
    TYPE(cp_dbcsr_type), POINTER             :: cosmat, mo_coeff_b, sinmat
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_coeff_tmp
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: inv_work, mo_derivs_tmp, &
                                                op_fm_set, opvec
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(efield_berry_type), POINTER         :: efield
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_kind_type), POINTER              :: qs_kind
    TYPE(virial_type), POINTER               :: virial

    CALL timeset(routineN,handle)

    NULLIFY(dft_control,cell,particle_set)
    CALL get_qs_env(qs_env,dft_control=dft_control,cell=cell,&
         particle_set=particle_set,virial=virial,error=error)
    NULLIFY(qs_kind_set,efield,para_env,sab_orb)
    CALL get_qs_env(qs_env=qs_env,qs_kind_set=qs_kind_set,&
         efield=efield,energy=energy,para_env=para_env,sab_orb=sab_orb,error=error)

    ! calculate stress only if forces requested also
    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)
    use_virial = use_virial .AND. calculate_forces
    ! disable stress calculation
    IF (use_virial) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="Stress tensor for periodic D-field not implemented", &
            error=error, error_level=cp_failure_level)
    END IF

    fieldpol=dft_control%period_efield%polarisation
    fieldpol=fieldpol/SQRT(DOT_PRODUCT(fieldpol,fieldpol))
    fieldpol=fieldpol*dft_control%period_efield%strength

    omega = cell%deth
    hmat = cell%hmat(:,:)/(twopi*omega)

    ! nuclear contribution to polarization
    natom = SIZE(particle_set)
    IF (calculate_forces) THEN
       CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,force=force,error=error)
       ALLOCATE (atom_of_kind(natom),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL get_atomic_kind_set(atomic_kind_set,atom_of_kind=atom_of_kind)
       ALLOCATE (force_tmp(natom,3,3),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       force_tmp = 0.0_dp
    END IF
    zi(:) = CMPLX ( 1._dp, 0._dp, dp )
    DO ia = 1,natom
       CALL get_atomic_kind(particle_set(ia)%atomic_kind,kind_number=ikind)
       CALL get_qs_kind(qs_kind_set(ikind),core_charge=charge,error=error)
       ria = particle_set(ia)%r
       ria = pbc(ria,cell)
       DO idir = 1, 3
          kvec(:) = twopi*cell%h_inv(idir,:)
          dd = SUM ( kvec(:) * ria(:) )
          zdeta = CMPLX(COS(dd),SIN(dd),KIND=dp)**charge
          zi(idir) = zi(idir) * zdeta
       END DO
       IF(calculate_forces)THEN
          IF(para_env%mepos==0) THEN
             DO i=1,3
                force_tmp(ia,i,i)=force_tmp(ia,i,i) + charge/omega
             END DO
          END IF
       END IF
    END DO
    rlog = AIMAG(LOG(zi))

    ! check uniform occupation
    NULLIFY(mos)
    CALL get_qs_env(qs_env=qs_env,mos=mos,error=error)
    DO ispin = 1,dft_control%nspins
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,maxocc=occ,uniform_occupation=uniform)
       IF (.NOT.uniform) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Berry phase moments for non uniform MOs' occupation numbers not implemented", &
               error=error, error_level=cp_failure_level)
       END IF
    END DO

    ! initialize all work matrices needed
    NULLIFY (mo_derivs)
    CALL get_qs_env(qs_env=qs_env,mo_derivs=mo_derivs,error=error)
    ALLOCATE ( op_fm_set( 2, dft_control%nspins ), STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE ( opvec( 2, dft_control%nspins ), STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE ( eigrmat( dft_control%nspins ), STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE ( inv_mat( dft_control%nspins ), STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE ( inv_work(2, dft_control%nspins ), STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE ( mo_derivs_tmp(3,SIZE(mo_derivs)), STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE ( mo_coeff_tmp(SIZE(mo_derivs)), STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ! Allocate temp matrices for the wavefunction derivatives
    DO ispin = 1, dft_control%nspins
       NULLIFY(tmp_fm_struct,mo_coeff)
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff,nao=nao,nmo=nmo)
       CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nmo,&
            ncol_global=nmo,para_env=para_env,context=mo_coeff%matrix_struct%context,error=error)
       CALL cp_fm_create (mo_coeff_tmp(ispin)%matrix,mo_coeff%matrix_struct,error=error)
       DO i=1,3
          CALL cp_fm_create (mo_derivs_tmp(i,ispin)%matrix,mo_coeff%matrix_struct,error=error)
          CALL cp_fm_set_all (matrix=mo_derivs_tmp(i,ispin)%matrix,alpha=0.0_dp,error=error)
       END DO
       DO i = 1, SIZE(op_fm_set,1)
          CALL cp_fm_create (opvec(i,ispin)%matrix,mo_coeff%matrix_struct,error=error)
          NULLIFY(op_fm_set(i,ispin)%matrix)
          CALL cp_fm_create (op_fm_set(i,ispin)%matrix,tmp_fm_struct,error=error)
          CALL cp_fm_create (inv_work(i,ispin)%matrix,op_fm_set(i,ispin)%matrix%matrix_struct,error=error)
       END DO
       CALL cp_cfm_create (eigrmat(ispin)%matrix,op_fm_set(1,ispin)%matrix%matrix_struct,error=error)
       CALL cp_cfm_create (inv_mat(ispin)%matrix,op_fm_set(1,ispin)%matrix%matrix_struct,error=error)
       CALL cp_fm_struct_release(tmp_fm_struct,error=error)
    END DO
    ! temp matrices for force calculation
    IF(calculate_forces)THEN
       NULLIFY (matrix_s)
       CALL get_qs_env(qs_env=qs_env,matrix_s=matrix_s,error=error)
       ALLOCATE(tempmat(2,dft_control%nspins),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO ispin=1,dft_control%nspins
          ALLOCATE(tempmat(1,ispin)%matrix,tempmat(2,ispin)%matrix,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          CALL cp_dbcsr_init(tempmat(1,ispin)%matrix,error=error)
          CALL cp_dbcsr_init(tempmat(2,ispin)%matrix, error=error)
          CALL cp_dbcsr_copy(tempmat(1,ispin)%matrix,matrix_s(1)%matrix,'TEMPMAT',error=error)
          CALL cp_dbcsr_copy(tempmat(2,ispin)%matrix,matrix_s(1)%matrix,'TEMPMAT',error=error)
          CALL cp_dbcsr_set(tempmat(1,ispin)%matrix,0.0_dp,error=error)
          CALL cp_dbcsr_set(tempmat(2,ispin)%matrix,0.0_dp,error=error)
       END DO
       ! integration
       CALL get_qs_kind_set(qs_kind_set,maxco=ldab,maxsgf=lsab,error=error)
       ALLOCATE(cosab(ldab,ldab),sinab(ldab,ldab),work(ldab,ldab),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(dcosab(ldab,ldab,3),dsinab(ldab,ldab,3),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       lsab = MAX(lsab,ldab)
       DO i=1,3
          ALLOCATE(dcost(i,1)%block(lsab,lsab),dsint(i,1)%block(lsab,lsab),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(dcost(i,2)%block(lsab,lsab),dsint(i,2)%block(lsab,lsab),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END DO
    END IF

    !Start the MO derivative calculation
    !loop over all cell vectors
    DO idir=1,3
       zi(idir) = zero
       cosmat => efield%cosmat(idir)%matrix
       sinmat => efield%sinmat(idir)%matrix
       !evaluate the expression needed for the derivative (S_berry * C  and [C^T S_berry C]^-1)
       !first step S_berry * C  and C^T S_berry C
       DO ispin=1, dft_control%nspins   ! spin
          IF(mos(ispin)%mo_set%use_mo_coeff_b)THEN
             CALL get_mo_set(mo_set=mos(ispin)%mo_set,nao=nao,mo_coeff_b=mo_coeff_b,nmo=nmo)
             CALL copy_dbcsr_to_fm(mo_coeff_b,mo_coeff_tmp(ispin)%matrix,error=error)
          ELSE
             CALL get_mo_set(mo_set=mos(ispin)%mo_set,nao=nao,mo_coeff=mo_coeff_tmp(ispin)%matrix,nmo=nmo)
          END IF
          CALL cp_dbcsr_sm_fm_multiply(cosmat,mo_coeff_tmp(ispin)%matrix,opvec(1,ispin)%matrix,ncol=nmo,error=error)
          CALL cp_gemm("T","N",nmo,nmo,nao,1.0_dp,mo_coeff_tmp(ispin)%matrix,opvec(1,ispin)%matrix,0.0_dp,&
               op_fm_set(1,ispin)%matrix,error=error)
          CALL cp_dbcsr_sm_fm_multiply(sinmat,mo_coeff_tmp(ispin)%matrix,opvec(2,ispin)%matrix,ncol=nmo,error=error)
          CALL cp_gemm("T","N",nmo,nmo,nao,1.0_dp,mo_coeff_tmp(ispin)%matrix,opvec(2,ispin)%matrix,0.0_dp,&
               op_fm_set(2,ispin)%matrix,error=error)
       ENDDO
       !second step invert C^T S_berry C
       zdet = one
       DO ispin = 1, dft_control%nspins
          CALL cp_cfm_add_fm(zero,eigrmat(ispin)%matrix,one,op_fm_set(1,ispin)%matrix,error)
          CALL cp_cfm_add_fm(one,eigrmat(ispin)%matrix,-zone,op_fm_set(2,ispin)%matrix,error)
          CALL cp_cfm_set_all(inv_mat(ispin)%matrix,zero,one,error)
          CALL cp_cfm_solve(eigrmat(ispin)%matrix,inv_mat(ispin)%matrix,zdeta,error)
          zdet = zdet*zdeta
       END DO
       zi(idir) = zdet**occ
       zlog(idir) = AIMAG(LOG(zi(idir)))

       IF(.NOT.just_energy) THEN
          !compute the orbital derivative
          DO ispin=1,dft_control%nspins
             inv_work(1,ispin)%matrix%local_data(:,:)=REAL(inv_mat(ispin)%matrix%local_data(:,:),dp)
             inv_work(2,ispin)%matrix%local_data(:,:)=AIMAG(inv_mat(ispin)%matrix%local_data(:,:))
             CALL get_mo_set(mo_set=mos(ispin)%mo_set,nao=nao,nmo=nmo)
             DO i=1,3
                focc = hmat(idir,i)
                CALL cp_gemm("N","N",nao,nmo,nmo,focc,opvec(1,ispin)%matrix,inv_work(2,ispin)%matrix,&
                     1.0_dp,mo_derivs_tmp(idir,ispin)%matrix,error)
                CALL cp_gemm("N","N",nao,nmo,nmo,-focc,opvec(2,ispin)%matrix,inv_work(1,ispin)%matrix,&
                     1.0_dp,mo_derivs_tmp(idir,ispin)%matrix,error)
             END DO
          END DO
       END IF

       !compute nuclear forces
       IF(calculate_forces)THEN
          nkind = SIZE(qs_kind_set)
          natom = SIZE(particle_set)
          kvec(:) = twopi*cell%h_inv(idir,:)

          ! calculate: C [C^T S_berry C]^(-1) C^T
          ! Store this matrix in DBCSR form (only S overlap blocks)
          DO ispin=1,dft_control%nspins
             CALL cp_dbcsr_set(tempmat(1,ispin)%matrix,0.0_dp,error=error)
             CALL cp_dbcsr_set(tempmat(2,ispin)%matrix,0.0_dp,error=error)
             CALL get_mo_set(mo_set=mos(ispin)%mo_set,nao=nao,nmo=nmo)
             CALL cp_gemm("N","N",nao,nmo,nmo,1.0_dp,mo_coeff_tmp(ispin)%matrix,inv_work(1,ispin)%matrix,0.0_dp,&
                  opvec(1,ispin)%matrix,error=error)
             CALL cp_gemm("N","N",nao,nmo,nmo,1.0_dp,mo_coeff_tmp(ispin)%matrix,inv_work(2,ispin)%matrix,0.0_dp,&
                  opvec(2,ispin)%matrix,error=error)
             CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=tempmat(1,ispin)%matrix,&
                  matrix_v=opvec(1,ispin)%matrix,matrix_g=mo_coeff_tmp(ispin)%matrix,ncol=nmo,error=error)
             CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=tempmat(2,ispin)%matrix,&
                  matrix_v=opvec(2,ispin)%matrix,matrix_g=mo_coeff_tmp(ispin)%matrix,ncol=nmo,error=error)
          END DO

          ! Calculation of derivative integrals (da|eikr|b) and (a|eikr|db)
          ALLOCATE (basis_set_list(nkind),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          DO ikind=1,nkind
             qs_kind => qs_kind_set(ikind)
             CALL get_qs_kind(qs_kind=qs_kind,basis_set=basis_set_a,error=error)
             IF (ASSOCIATED(basis_set_a)) THEN
                basis_set_list(ikind)%gto_basis_set => basis_set_a
             ELSE
                NULLIFY(basis_set_list(ikind)%gto_basis_set)
             END IF
          END DO
          !
          CALL neighbor_list_iterator_create(nl_iterator,sab_orb)
          DO WHILE (neighbor_list_iterate(nl_iterator)==0)
             CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,inode=inode,&
                                    iatom=iatom,jatom=jatom,r=rab)
             basis_set_a => basis_set_list(ikind)%gto_basis_set
             IF (.NOT.ASSOCIATED(basis_set_a)) CYCLE
             basis_set_b => basis_set_list(jkind)%gto_basis_set
             IF (.NOT.ASSOCIATED(basis_set_b)) CYCLE
             ! basis ikind
             first_sgfa   =>  basis_set_a%first_sgf
             la_max       =>  basis_set_a%lmax
             la_min       =>  basis_set_a%lmin
             npgfa        =>  basis_set_a%npgf
             nseta        =   basis_set_a%nset
             nsgfa        =>  basis_set_a%nsgf_set
             rpgfa        =>  basis_set_a%pgf_radius
             set_radius_a =>  basis_set_a%set_radius
             sphi_a       =>  basis_set_a%sphi
             zeta         =>  basis_set_a%zet
             ! basis jkind
             first_sgfb   =>  basis_set_b%first_sgf
             lb_max       =>  basis_set_b%lmax
             lb_min       =>  basis_set_b%lmin
             npgfb        =>  basis_set_b%npgf
             nsetb        =   basis_set_b%nset
             nsgfb        =>  basis_set_b%nsgf_set
             rpgfb        =>  basis_set_b%pgf_radius
             set_radius_b =>  basis_set_b%set_radius
             sphi_b       =>  basis_set_b%sphi
             zetb         =>  basis_set_b%zet

             ldsa = SIZE(sphi_a,1)
             ldsb = SIZE(sphi_b,1)
             ra(:) = pbc(particle_set(iatom)%r(:),cell)
             rb(:) = ra+rab
             dab = SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))

             IF (iatom <= jatom) THEN
                irow = iatom
                icol = jatom
             ELSE
                irow = jatom
                icol = iatom
             END IF

             IF (iatom == jatom .AND. dab < 1.e-10_dp) THEN
                fab = 1.0_dp * occ
             ELSE
                fab = 2.0_dp * occ
             END IF

             DO i=1,3
                dcost(i,1)%block = 0.0_dp
                dsint(i,1)%block = 0.0_dp
                dcost(i,2)%block = 0.0_dp
                dsint(i,2)%block = 0.0_dp
             END DO

             DO iset=1,nseta
                ncoa = npgfa(iset)*ncoset(la_max(iset))
                sgfa = first_sgfa(1,iset)
                DO jset=1,nsetb
                   IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE
                   ncob = npgfb(jset)*ncoset(lb_max(jset))
                   sgfb = first_sgfb(1,jset)
                   ! Calculate the primitive integrals (da|b)
                   CALL cossin(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                        lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                        ra,rb,kvec,cosab,sinab,dcosab,dsinab)
                   DO i=1,3
                      CALL contract_all(dcost(i,1)%block,dsint(i,1)%block,&
                           ncoa,nsgfa(iset),sgfa,sphi_a,ldsa,&
                           ncob,nsgfb(jset),sgfb,sphi_b,ldsb,&
                           dcosab(:,:,i),dsinab(:,:,i),ldab,work,ldab)
                   END DO
                   ! Calculate the primitive integrals (a|db)
                   CALL cossin(lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                        la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                        rb,ra,kvec,cosab,sinab,dcosab,dsinab)
                   DO i=1,3
                      dcosab(1:ncoa,1:ncob,i) = TRANSPOSE(dcosab(1:ncob,1:ncoa,i))
                      dsinab(1:ncoa,1:ncob,i) = TRANSPOSE(dsinab(1:ncob,1:ncoa,i))
                      CALL contract_all(dcost(i,2)%block,dsint(i,2)%block,&
                           ncoa,nsgfa(iset),sgfa,sphi_a,ldsa,&
                           ncob,nsgfb(jset),sgfb,sphi_b,ldsb,&
                           dcosab(:,:,i),dsinab(:,:,i),ldab,work,ldab)
                   END DO
                END DO
             END DO
             forcea = 0.0_dp
             forceb = 0.0_dp
             DO ispin=1,dft_control%nspins
                NULLIFY (rblock,iblock)
                CALL cp_dbcsr_get_block_p(matrix=tempmat(1,ispin)%matrix,&
                     row=irow,col=icol,BLOCK=rblock,found=found)
                CPPostcondition(found,cp_failure_level,routineP,error,failure)
                CALL cp_dbcsr_get_block_p(matrix=tempmat(2,ispin)%matrix,&
                     row=irow,col=icol,BLOCK=iblock,found=found)
                CPPostcondition(found,cp_failure_level,routineP,error,failure)
                n1 = SIZE(rblock,1)
                n2 = SIZE(rblock,2)
                CPPostcondition(SIZE(iblock,1)==n1,cp_failure_level,routineP,error,failure)
                CPPostcondition(SIZE(iblock,2)==n2,cp_failure_level,routineP,error,failure)
                CPPostcondition(lsab>=n1,cp_failure_level,routineP,error,failure)
                CPPostcondition(lsab>=n2,cp_failure_level,routineP,error,failure)
                IF (iatom <= jatom) THEN
                   DO i=1,3
                      forcea(i) = forcea(i) + SUM(rblock(1:n1,1:n2)*dsint(i,1)%block(1:n1,1:n2)) &
                                            - SUM(iblock(1:n1,1:n2)*dcost(i,1)%block(1:n1,1:n2))
                      forceb(i) = forceb(i) + SUM(rblock(1:n1,1:n2)*dsint(i,2)%block(1:n1,1:n2)) &
                                            - SUM(iblock(1:n1,1:n2)*dcost(i,2)%block(1:n1,1:n2))
                   END DO
                ELSE
                   DO i=1,3
                      forcea(i) = forcea(i) + SUM(TRANSPOSE(rblock(1:n1,1:n2))*dsint(i,1)%block(1:n2,1:n1)) &
                                            - SUM(TRANSPOSE(iblock(1:n1,1:n2))*dcost(i,1)%block(1:n2,1:n1))
                      forceb(i) = forceb(i) + SUM(TRANSPOSE(rblock(1:n1,1:n2))*dsint(i,2)%block(1:n2,1:n1)) &
                                            - SUM(TRANSPOSE(iblock(1:n1,1:n2))*dcost(i,2)%block(1:n2,1:n1))
                    END DO
                END IF
             END DO
             DO i=1,3
                force_tmp(iatom,:,i)=force_tmp(iatom,:,i)-fab*hmat(i,idir)*forcea(:)
                force_tmp(jatom,:,i)=force_tmp(jatom,:,i)-fab*hmat(i,idir)*forceb(:)
             END DO
          END DO
          CALL neighbor_list_iterator_release(nl_iterator)
          DEALLOCATE (basis_set_list,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END IF
    END DO

    ! make sure the total normalized polarization is within [-1:1]
    DO idir=1,3
       IF( rlog(idir)+zlog(idir) > pi ) zlog(idir) = zlog(idir) - twopi
       IF( rlog(idir)+zlog(idir) < -pi ) zlog(idir) = zlog(idir) + twopi
    END DO
    DO idir=1,3
       qi(idir) = 0.0_dp
       ci(idir) = 0.0_dp
       DO i=1,3
          qi(idir) = qi(idir) + hmat(idir,i) * rlog(i)
          ci(idir) = ci(idir) + hmat(idir,i) * zlog(i)
       END DO
    END DO

    ! Polarization Energy
    ener_field=SUM(qi+ci)
    ! check for smoothness of energy surface
    IF(efield%field_energy == 0.0_dp) THEN
       efield%field_energy = ener_field
    ELSEIF( ABS(efield%field_energy-ener_field) > pi*ABS(SUM(hmat)) ) THEN
       CALL cp_error_message(cp_note_level, routineP, &
        "Large change of e-field energy. Non-smooth energy surface", error)
    ELSE
       efield%field_energy = ener_field
    END IF

    ! Energy
    ener_field = 0.0_dp
    DO i=1,3
       ener_field = ener_field + (fieldpol(i) - 2._dp*twopi*(ci(i)+qi(i)))**2
    END DO
    energy%efield = 0.25_dp*omega/twopi * ener_field

    IF(.NOT.just_energy) THEN
       DO i=1,3
          di(i) = -omega * (fieldpol(i) - 2._dp*twopi*(ci(i)+qi(i)))
       END DO
       ! Add the result to mo_derivativs
       DO ispin=1, dft_control%nspins
          CALL copy_dbcsr_to_fm(mo_derivs(ispin)%matrix,mo_coeff_tmp(ispin)%matrix,error=error)
          DO idir=1,3
             CALL cp_fm_scale_and_add(1.0_dp,mo_coeff_tmp(ispin)%matrix,di(idir),&
                  mo_derivs_tmp(idir,ispin)%matrix,error)
          END DO
       END DO
       DO ispin=1, dft_control%nspins
          CALL copy_fm_to_dbcsr(mo_coeff_tmp(ispin)%matrix,mo_derivs(ispin)%matrix,error=error)
       END DO
    END IF

    IF(calculate_forces)THEN
       DO i=1,3
          DO ia = 1,natom
             CALL get_atomic_kind(particle_set(ia)%atomic_kind,kind_number=ikind)
             iatom = atom_of_kind(ia)
             force(ikind)%efield(1:3,iatom)=force(ikind)%efield(1:3,iatom) + di(i)*force_tmp(ia,1:3,i)
          END DO
       END DO
    END IF

    DO ispin = 1, dft_control%nspins
       CALL cp_cfm_release(eigrmat(ispin)%matrix,error=error)
       CALL cp_cfm_release(inv_mat(ispin)%matrix,error=error)
       CALL cp_fm_release(mo_coeff_tmp(ispin)%matrix,error=error)
       DO i=1,3
          CALL cp_fm_release(mo_derivs_tmp(i,ispin)%matrix,error=error)
       END DO
       DO i = 1, SIZE ( op_fm_set, 1 )
          CALL cp_fm_release(opvec(i,ispin)%matrix,error=error)
          CALL cp_fm_release(op_fm_set(i,ispin)%matrix,error=error)
          CALL cp_fm_release(inv_work(i,ispin)%matrix,error=error)
       END DO
    END DO
    DEALLOCATE(inv_mat,inv_work,op_fm_set,opvec,eigrmat,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(mo_coeff_tmp,mo_derivs_tmp,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    IF(calculate_forces)THEN
       DO ikind=1,SIZE(atomic_kind_set)
          CALL mp_sum(force(ikind)%efield,para_env%group)
       END DO
       DEALLOCATE(atom_of_kind,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (force_tmp,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(cosab,sinab,work,dcosab,dsinab,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DO i=1,3
          DEALLOCATE(dcost(i,1)%block,dsint(i,1)%block,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(dcost(i,2)%block,dsint(i,2)%block,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END DO
       CALL cp_dbcsr_deallocate_matrix_set(tempmat,error)
    END IF
    CALL timestop(handle)

  END SUBROUTINE qs_dispfield_derivatives

! *****************************************************************************
!> \brief ...
!> \param cos_block ...
!> \param sin_block ...
!> \param ncoa ...
!> \param nsgfa ...
!> \param sgfa ...
!> \param sphi_a ...
!> \param ldsa ...
!> \param ncob ...
!> \param nsgfb ...
!> \param sgfb ...
!> \param sphi_b ...
!> \param ldsb ...
!> \param cosab ...
!> \param sinab ...
!> \param ldab ...
!> \param work ...
!> \param ldwork ...
! *****************************************************************************
  SUBROUTINE contract_all(cos_block, sin_block,&
       ncoa,nsgfa,sgfa,sphi_a,ldsa,&
       ncob,nsgfb,sgfb,sphi_b,ldsb,&
       cosab,sinab,ldab,work,ldwork)

    REAL(dp), DIMENSION(:, :), POINTER       :: cos_block, sin_block
    INTEGER, INTENT(IN)                      :: ncoa, nsgfa, sgfa
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: sphi_a
    INTEGER, INTENT(IN)                      :: ldsa, ncob, nsgfb, sgfb
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: sphi_b
    INTEGER, INTENT(IN)                      :: ldsb
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: cosab, sinab
    INTEGER, INTENT(IN)                      :: ldab
    REAL(dp), DIMENSION(:, :)                :: work
    INTEGER, INTENT(IN)                      :: ldwork

! Calculate cosine

    CALL dgemm("N","N",ncoa,nsgfb,ncob,1.0_dp,cosab(1,1),ldab,&
         sphi_b(1,sgfb),ldsb,0.0_dp,work(1,1),ldwork)

    CALL dgemm("T","N",nsgfa,nsgfb,ncoa,1.0_dp,sphi_a(1,sgfa),ldsa,&
        work(1,1),ldwork,1.0_dp,cos_block(sgfa,sgfb),SIZE(cos_block,1))

    ! Calculate sine
    CALL dgemm("N","N",ncoa,nsgfb,ncob,1.0_dp,sinab(1,1),ldab,&
         sphi_b(1,sgfb),ldsb,0.0_dp,work(1,1),ldwork)

    CALL dgemm("T","N",nsgfa,nsgfb,ncoa,1.0_dp,sphi_a(1,sgfa),ldsa,&
         work(1,1),ldwork,1.0_dp,sin_block(sgfa,sgfb),SIZE(sin_block,1))

  END SUBROUTINE contract_all

END MODULE qs_efield_berry
