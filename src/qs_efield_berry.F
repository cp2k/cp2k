!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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
  USE cell_methods,                    ONLY: cell_type,&
                                             pbc
  USE cp_cfm_basic_linalg,             ONLY: cp_cfm_solve
  USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                             cp_cfm_get_info,&
                                             cp_cfm_p_type,&
                                             cp_cfm_release,&
                                             cp_cfm_set_all
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_allocate_matrix_set, cp_dbcsr_create, cp_dbcsr_finalize, &
       cp_dbcsr_get_block_p, cp_dbcsr_init, cp_dbcsr_p_type, cp_dbcsr_set, &
       cp_dbcsr_type, dbcsr_type_no_symmetry
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             copy_fm_to_dbcsr,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_transpose
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: twopi
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: indco,&
                                             ncoset
  USE orbital_symbols,                 ONLY: cgf_symbol
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE qs_period_efield_types,          ONLY: efield_berry_type,&
                                             init_efield_matrices,&
                                             set_efield_matrices
  USE string_utilities,                ONLY: compress,&
                                             uppercase
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_efield_berry'

  ! *** Public subroutines ***

  PUBLIC :: qs_efield_berry_phase

CONTAINS

  ! *****************************************************************************

  SUBROUTINE qs_efield_berry_phase(ks_env,qs_env,calculate_forces,error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL                                  :: calculate_forces
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_efield_berry_phase', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)
    IF(qs_env%dft_control%apply_period_efield)THEN
       IF(calculate_forces.OR.ks_env%s_mstruct_changed)&
            CALL qs_efield_integrals (qs_env,calculate_forces,error)
       CALL qs_efield_mo_derivatives(qs_env,calculate_forces,error)
    END IF
    CALL timestop(handle)
  END SUBROUTINE qs_efield_berry_phase
  ! *****************************************************************************
  SUBROUTINE qs_efield_integrals (qs_env,calculate_forces,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL                                  :: calculate_forces
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_efield_integrals', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, stat
    LOGICAL                                  :: failure
    REAL(dp), DIMENSION(3)                   :: kvec
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: cosmat, dcosmat, dsinmat, &
                                                matrix_s, sinmat
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(efield_berry_type), POINTER         :: efield
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_all

    CALL timeset(routineN,handle)
    failure = .FALSE.
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN
       CALL get_qs_env(qs_env=qs_env,dft_control=dft_control, error=error)
       NULLIFY (matrix_s)
       CALL get_qs_env(qs_env=qs_env,efield=efield,cell=cell, matrix_s=matrix_s,sab_all=sab_all, error=error)
       CALL init_efield_matrices(efield,error)
       ALLOCATE(cosmat(3),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(sinmat(3),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF(calculate_forces)THEN
          ALLOCATE(efield%dcosmat(3),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(efield%dsinmat(3),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       DO i=1,3
          ALLOCATE(cosmat(i)%matrix,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(sinmat(i)%matrix,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          CALL cp_dbcsr_init(cosmat(i)%matrix,error=error)
          CALL cp_dbcsr_init(sinmat(i)%matrix, error=error)
          CALL cp_dbcsr_create(matrix=sinmat(i)%matrix, template=matrix_s(1)%matrix, &
               name="sinmat", matrix_type=dbcsr_type_no_symmetry, error=error)
          CALL cp_dbcsr_alloc_block_from_nbl(sinmat(i)%matrix,sab_all,error=error)
          CALL cp_dbcsr_create(matrix=cosmat(i)%matrix, template=matrix_s(1)%matrix, &
               name="cosmat", matrix_type=dbcsr_type_no_symmetry, error=error)
          CALL cp_dbcsr_alloc_block_from_nbl(cosmat(i)%matrix,sab_all,error=error)
          CALL cp_dbcsr_set(cosmat(i)%matrix,0.0_dp,error=error)
          CALL cp_dbcsr_set(sinmat(i)%matrix,0.0_dp,error=error)
          kvec(:) = twopi*cell%h_inv(i,:)
          IF(calculate_forces)THEN
             CALL berry_mat_derivs(qs_env,cosmat(i)%matrix,sinmat(i)%matrix,kvec,calculate_forces,&
                  dcosmat=dcosmat,dsinmat=dsinmat,error=error)
             CALL set_efield_matrices(efield=efield,cosmat=cosmat,sinmat=sinmat,dcosmat=dcosmat,dsinmat=dsinmat,fielddir=i)
          ELSE
             CALL berry_mat_derivs(qs_env,cosmat(i)%matrix,sinmat(i)%matrix,kvec,calculate_forces,error=error)
             CALL set_efield_matrices(efield=efield,cosmat=cosmat,sinmat=sinmat)
          END IF
       END DO
       CALL set_qs_env(qs_env=qs_env,efield=efield,error=error)
    END IF
    CALL timestop(handle)

  END SUBROUTINE qs_efield_integrals

  ! *****************************************************************************

  SUBROUTINE qs_efield_mo_derivatives(qs_env,calculate_forces,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL                                  :: calculate_forces
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_efield_mo_derivatives', &
      routineP = moduleN//':'//routineN
    COMPLEX(KIND=dp), PARAMETER              :: one = (1.0_dp,0.0_dp) , &
                                                zero = (0.0_dp,0.0_dp)

    COMPLEX(dp)                              :: xphase(3), zdet, zdeta, &
                                                zi(3), zphase(3)
    INTEGER :: fdir, handle, i, ia, iatom, icol_atom, icol_global, &
      icol_local, idim, idir, ikind, irow_atom, irow_global, irow_local, &
      ispin, istat, kind_atom, kkk, nao, natom, ncol_local, nkind, nmo, &
      nmotot, nrow_local, stat, tmp_dim
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, col_atom_index, &
                                                first_sgf, kind_of, last_sgf, &
                                                row_atom_index
    INTEGER, DIMENSION(:), POINTER           :: atom_list, col_indices, &
                                                row_indices
    LOGICAL                                  :: failure, uniform
    REAL(dp)                                 :: charge, ci(3), dd, &
                                                ener_field, fieldfac(3), &
                                                fieldpol(3), myfac, occ
    REAL(dp), DIMENSION(3)                   :: kvec, rcc, ria
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_cfm_p_type), DIMENSION(:), &
      POINTER                                :: eigrmat, inv_mat
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s, mo_derivs
    TYPE(cp_dbcsr_type), POINTER             :: cosmat, dcosmat, dsinmat, &
                                                mo_coeff_b, sinmat
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_coeff_tmp, mo_derivs_tmp
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: inv_work, op_fm_set, opvec
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: ao_work, dcos_full, &
                                                dsin_full, mo_coeff, tmp_im, &
                                                tmp_re
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(efield_berry_type), POINTER         :: efield
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force

    CALL timeset(routineN,handle)
    NULLIFY ( dft_control, cell, particle_set )
    CALL get_qs_env ( qs_env, dft_control=dft_control, cell=cell,&
         particle_set=particle_set, force=force, error=error)
    NULLIFY ( local_particles,atom_list )
    CALL get_qs_env ( qs_env, local_particles=local_particles, error=error)
    NULLIFY ( matrix_s, mos )
    CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s, mos=mos,para_env=para_env, error=error)

    fieldpol=dft_control%period_efield%polarisation
    fieldpol=fieldpol/SQRT(DOT_PRODUCT(fieldpol,fieldpol))
    fieldpol=fieldpol*dft_control%period_efield%strength

    rcc=(/0.0_dp,0.0_dp,0.0_dp/)

    ! nuclear contribution
    ria = twopi * MATMUL ( cell%h_inv, rcc )
    zphase = CMPLX ( COS(ria), SIN(ria), dp )
    DO idir=1,3
       fieldfac(idir)=fieldpol(idir) * SQRT(DOT_PRODUCT(cell%hmat(idir,:),cell%hmat(idir,:)))/twopi
    END DO

    zi = 0._dp
    NULLIFY ( matrix_s, mos,dft_control )
    CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s, mos=mos,&
         dft_control=dft_control,efield=efield,cell=cell,&
         mo_derivs=mo_derivs,energy=energy,para_env=para_env, error=error)

    zi(:) = CMPLX ( 1._dp, 0._dp, dp )
    DO ia = 1,SIZE(particle_set)
       atomic_kind => particle_set(ia)%atomic_kind
       CALL get_atomic_kind(atomic_kind=atomic_kind,core_charge=charge)
       ria = particle_set(ia)%r
       ria = pbc(ria,cell)
       DO i = 1, 3
          kvec(:) = twopi*cell%h_inv(i,:)
          dd = SUM ( kvec(:) * ria(:) )
          zdeta = CMPLX(COS(dd),SIN(dd),KIND=dp)**charge
          zi(i) = zi(i) * zdeta
       END DO
    END DO
    zi = zi * zphase
    ci = AIMAG(LOG(zi))
    ener_field=0.0_dp
    DO idir=1,3
       ener_field=ener_field+ci(idir)*fieldfac(idir)
    END DO

    IF(calculate_forces)THEN
       CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,force=force,error=error)
       nkind=SIZE(atomic_kind_set)
       kkk=1
       DO ikind=1,SIZE(atomic_kind_set)
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               atom_list=atom_list,&
               natom=natom,&
               zeff=charge)
          natom = SIZE(atom_list)
          DO iatom=1,natom
             DO idir=1,3

                IF(para_env%mepos==0)&
                   force(ikind)%efield(idir,iatom)=force(ikind)%efield(idir,iatom)+0.5_dp*fieldfac(idir)*charge
             END DO
          END DO
       END DO
    END IF
    ! occupation
    DO ispin = 1, dft_control%nspins
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,maxocc=occ,uniform_occupation=uniform)
       IF (.NOT.uniform) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Berry phase moments for non uniform MOs' occupation numbers not implemented", &
               error=error, error_level=cp_failure_level)
       END IF
    END DO

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

    NULLIFY(mo_derivs_tmp)!dbcsr->fm
    ALLOCATE(mo_derivs_tmp(SIZE(mo_derivs)),stat=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(mo_coeff_tmp)!dbcsr->fm
    ALLOCATE(mo_coeff_tmp(SIZE(mo_derivs)),stat=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    nmotot = 0

    ! A bit to allocate for the wavefunction
    DO ispin = 1, dft_control%nspins
       NULLIFY(tmp_fm_struct,mo_coeff)
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff,nao=nao,nmo=nmo)
       nmotot = nmotot + nmo
       CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nmo,&
            ncol_global=nmo,para_env=para_env,context=mo_coeff%matrix_struct%context,&
            error=error)
       CALL cp_fm_create (mo_derivs_tmp(ispin)%matrix , mo_coeff%matrix_struct ,error=error)
       CALL cp_fm_create (mo_coeff_tmp(ispin)%matrix , mo_coeff%matrix_struct ,error=error)
       CALL copy_dbcsr_to_fm(mo_derivs(ispin)%matrix,mo_derivs_tmp(ispin)%matrix,error=error)
       DO i = 1, SIZE ( op_fm_set, 1 )
          CALL cp_fm_create (opvec(i,ispin)%matrix , mo_coeff%matrix_struct ,error=error)
          NULLIFY(op_fm_set(i,ispin)%matrix)
          CALL cp_fm_create (op_fm_set(i,ispin)%matrix , tmp_fm_struct ,error=error)
          CALL cp_fm_create ( inv_work(i,ispin)%matrix, op_fm_set(i,ispin)%matrix%matrix_struct ,&
               error=error)
       END DO
       CALL cp_cfm_create ( eigrmat(ispin)%matrix, op_fm_set(1,ispin)%matrix%matrix_struct ,&
            error=error)
       CALL cp_cfm_create ( inv_mat(ispin)%matrix, op_fm_set(1,ispin)%matrix%matrix_struct ,&
            error=error)
       CALL cp_fm_struct_release(tmp_fm_struct,error=error)
    END DO

    ! A lot to allocate for the forces
    IF(calculate_forces)THEN
       CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
            ncol_global=nao,para_env=para_env,context=mo_coeff%matrix_struct%context,&
            error=error)
       NULLIFY (ao_work,dcos_full,dsin_full,tmp_re,tmp_im)
       CALL cp_fm_create ( ao_work, tmp_fm_struct,error=error)
       CALL cp_fm_create ( dcos_full, tmp_fm_struct,error=error)
       CALL cp_fm_create ( dsin_full, tmp_fm_struct,error=error)
       CALL cp_fm_create ( tmp_re, tmp_fm_struct,error=error)
       CALL cp_fm_create ( tmp_im, tmp_fm_struct,error=error)
       CALL cp_fm_struct_release(tmp_fm_struct,error=error)
       natom = SIZE(particle_set)
       ALLOCATE (first_sgf(natom),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (last_sgf(natom),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       CALL get_particle_set(particle_set=particle_set,first_sgf=first_sgf,last_sgf=last_sgf,&
            error=error)
       CALL cp_fm_get_info(dcos_full,ncol_local=ncol_local,nrow_local=nrow_local,&
            row_indices=row_indices, col_indices=col_indices, error=error)
       ALLOCATE (atom_of_kind(natom),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (kind_of(natom),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(row_atom_index(SIZE(row_indices)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind,kind_of=kind_of)
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
    END IF

    !Start the MO derivative calculation
    ria = twopi * REAL(nmotot,dp) * MATMUL ( cell%h_inv, rcc )
    xphase = CMPLX ( COS(ria), SIN(ria), dp )

    zi=0.0_dp

    ! loop over all cell vectors
    fieldfac=fieldfac*occ
    DO idir=1,3
       cosmat=>efield%cosmat(idir)%matrix
       sinmat=>efield%sinmat(idir)%matrix
       IF(ABS(fieldfac(idir)).GT.1.0E-10_dp)THEN
          !evaluate the expression needed for the derivative (S_berry * C  and [C^T S_berry C]^-1)
          !first step S_berry * C  and C^T S_berry C
          DO ispin=1, dft_control%nspins   ! spin
             IF(mos(ispin)%mo_set%use_mo_coeff_b)THEN
                CALL get_mo_set(mo_set=mos(ispin)%mo_set,nao=nao,mo_coeff_b=mo_coeff_b,nmo=nmo)
                CALL copy_dbcsr_to_fm(mo_coeff_b,mo_coeff_tmp(ispin)%matrix,error=error)
             ELSE
                CALL get_mo_set(mo_set=mos(ispin)%mo_set,nao=nao,mo_coeff=mo_coeff_tmp(ispin)%matrix,nmo=nmo)
             END IF
             CALL cp_dbcsr_sm_fm_multiply(cosmat, mo_coeff_tmp(ispin)%matrix, opvec(1,ispin)%matrix, ncol=nmo, error=error)
             CALL cp_fm_gemm("T","N",nmo,nmo,nao,1.0_dp,mo_coeff_tmp(ispin)%matrix,opvec(1,ispin)%matrix,0.0_dp,&
                  op_fm_set(1,ispin)%matrix,error=error)
             CALL cp_dbcsr_sm_fm_multiply(sinmat, mo_coeff_tmp(ispin)%matrix, opvec(2,ispin)%matrix, ncol=nmo, error=error)
             CALL cp_fm_gemm("T","N",nmo,nmo,nao,1.0_dp,mo_coeff_tmp(ispin)%matrix,opvec(2,ispin)%matrix,0.0_dp,&
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

             CALL cp_cfm_set_all (inv_mat(ispin)%matrix,zero,one,error)
             CALL cp_cfm_solve ( eigrmat(ispin)%matrix,inv_mat(ispin)%matrix, zdeta,error )
             zdet=zdet*zdeta
             zi(idir) = zdet
          END DO
          zi(idir)=zi(idir)*xphase(idir)
          ci(idir)=AIMAG(LOG(zi(idir)))
          ener_field=ener_field+ci(idir)*fieldfac(idir)


          !compute the derivative and add the result to mo_derivativs
          DO ispin=1,dft_control%nspins
             CALL cp_cfm_get_info(eigrmat(ispin)%matrix,ncol_local=tmp_dim,error=error)
             CALL get_mo_set(mo_set=mos(ispin)%mo_set,nao=nao,nmo=nmo)
             DO icol_local=1,tmp_dim
                inv_work(1,ispin)%matrix%local_data(:,icol_local)=REAL(inv_mat(ispin)%matrix%local_data(:,icol_local),dp)
                inv_work(2,ispin)%matrix%local_data(:,icol_local)=AIMAG(inv_mat(ispin)%matrix%local_data(:,icol_local))
             END DO
             CALL cp_fm_gemm("N","N",nao,nmo,nmo,fieldfac(idir)/(occ),opvec(1,ispin)%matrix,inv_work(2,ispin)%matrix,&
                  1.0_dp,mo_derivs_tmp(ispin)%matrix,error)
             CALL cp_fm_gemm("N","N",nao,nmo,nmo,-fieldfac(idir)/(occ),opvec(2,ispin)%matrix,inv_work(1,ispin)%matrix,&
                  1.0_dp,mo_derivs_tmp(ispin)%matrix,error)

          END DO

          !now states the funny part, get the derivative w.r.t. the nuclei

          IF(calculate_forces)THEN
             DO ispin=1,dft_control%nspins
                CALL cp_fm_get_info(dcos_full,ncol_local=ncol_local,nrow_local=nrow_local,&
                      row_indices=row_indices, col_indices=col_indices,error=error)
                CALL get_mo_set(mo_set=mos(ispin)%mo_set,nao=nao,nmo=nmo)
                CALL cp_fm_gemm("N","N",nao,nmo,nmo,1.0_dp,mo_coeff_tmp(ispin)%matrix,inv_work(1,ispin)%matrix,0.0_dp,&
                     opvec(1,ispin)%matrix,error=error)
                CALL cp_fm_gemm("N","N",nao,nmo,nmo,1.0_dp,mo_coeff_tmp(ispin)%matrix,inv_work(2,ispin)%matrix,0.0_dp,&
                     opvec(2,ispin)%matrix,error=error)
                CALL cp_fm_gemm("N","T",nao,nao,nmo,1.0_dp,opvec(1,ispin)%matrix,mo_coeff_tmp(ispin)%matrix,0.0_dp,&
                     tmp_re,error=error)
                CALL cp_fm_gemm("N","T",nao,nao,nmo,1.0_dp,opvec(2,ispin)%matrix,mo_coeff_tmp(ispin)%matrix,0.0_dp,&
                     tmp_im,error=error)
                myfac=-0.25_dp*SQRT(DOT_PRODUCT(cell%hmat(idir,:),cell%hmat(idir,:)))/twopi*fieldfac(idir)
                DO fdir=1,3
                   dsinmat=>efield%dsinmat(idir)%deriv(fdir)%matrix
                   dcosmat=>efield%dcosmat(idir)%deriv(fdir)%matrix
                   CALL copy_dbcsr_to_fm(dcosmat,dcos_full,error=error)
                   CALL copy_dbcsr_to_fm(dsinmat,dsin_full,error=error)
                   DO icol_local=1,ncol_local
                      DO irow_local=1,nrow_local
                         icol_atom=col_atom_index(icol_local)
                         irow_atom=row_atom_index(irow_local)
                         ikind=kind_of(irow_atom)
                         kind_atom=atom_of_kind(irow_atom)
                         force(ikind)%efield(fdir,kind_atom)=force(ikind)%efield(fdir,kind_atom) - &
                              tmp_im%local_data(irow_local,icol_local) * &
                              dcos_full%local_data(irow_local,icol_local)*myfac+ &
                              tmp_re%local_data(irow_local,icol_local) * &
                              dsin_full%local_data(irow_local,icol_local)*myfac
                       END DO
                   END DO
                   CALL cp_fm_transpose(dcos_full,ao_work,error)
                   CALL cp_fm_transpose(dsin_full,ao_work,error)
                   DO icol_local=1,ncol_local
                      DO irow_local=1,nrow_local
                         icol_atom=col_atom_index(icol_local)
                         irow_atom=row_atom_index(irow_local)
                         ikind=kind_of(icol_atom)
                         kind_atom=atom_of_kind(icol_atom)
                         force(ikind)%efield(fdir,kind_atom)=force(ikind)%efield(fdir,kind_atom) - &
                              tmp_im%local_data(irow_local,icol_local) * &
                              dcos_full%local_data(irow_local,icol_local)*myfac+ &
                              tmp_re%local_data(irow_local,icol_local) * &
                              dsin_full%local_data(irow_local,icol_local)*myfac
                      END DO
                   END DO
                END DO
             END DO

          END IF
       END IF

    END DO

    DO ispin=1, dft_control%nspins
       CALL copy_fm_to_dbcsr(mo_derivs_tmp(ispin)%matrix,mo_derivs(ispin)%matrix,error=error)
    END DO

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
    energy%efield=ener_field
    DEALLOCATE ( inv_mat, STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE ( mo_coeff_tmp,STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE ( inv_work, STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE ( mo_derivs_tmp, STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE ( op_fm_set, STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE ( opvec, STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE ( eigrmat, STAT = istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    IF(calculate_forces)THEN
       DO ikind=1,SIZE(atomic_kind_set)
          CALL mp_sum(force(ikind)%efield,para_env%group)
       END DO
       CALL cp_fm_release(dcos_full,error=error)
       CALL cp_fm_release(dsin_full,error=error)
       CALL cp_fm_release(tmp_re,error=error)
       CALL cp_fm_release(tmp_im,error=error)
       CALL cp_fm_release(ao_work,error=error)
       DEALLOCATE(atom_of_kind,stat=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(kind_of,stat=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle)

  END SUBROUTINE qs_efield_mo_derivatives



  SUBROUTINE berry_mat_derivs(qs_env,cosmat,sinmat,kvec,calculate_forces,dcosmat,dsinmat,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_type), POINTER             :: cosmat, sinmat
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: kvec
    LOGICAL                                  :: calculate_forces
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: dcosmat, dsinmat
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'berry_mat_derivs', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=12)                        :: cgfsym
    CHARACTER(LEN=default_string_length)     :: name
    INTEGER :: handle, iatom, icol, ikind, inode, irow, iset, istat, j, &
      jatom, jkind, jset, ldab, ldsa, ldsb, ldwork, m, natom, ncoa, ncob, &
      nkind, nseta, nsetb, sgfa, sgfb, stat
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, lb_max, &
                                                lb_min, npgfa, npgfb, nsgfa, &
                                                nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb
    LOGICAL                                  :: failure, found
    REAL(dp), DIMENSION(:, :), POINTER       :: cblock, cosab, sblock, sinab, &
                                                work
    REAL(dp), DIMENSION(:, :, :), POINTER    :: dcosab, dsinab
    REAL(KIND=dp)                            :: dab
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab, rb
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rpgfb, sphi_a, sphi_b, &
                                                zeta, zetb
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(block_p_type), DIMENSION(3)         :: cost, sint
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_all, sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY (atomic_kind_set, particle_set, sab_orb, sab_all, cell)
    CALL get_qs_env(qs_env=qs_env,&
         atomic_kind_set=atomic_kind_set,&
         particle_set=particle_set,cell=cell,&
         sab_orb=sab_orb,sab_all=sab_all,matrix_s=matrix_s, error=error)

    CALL cp_dbcsr_set(sinmat,0.0_dp,error=error)
    CALL cp_dbcsr_set(cosmat,0.0_dp,error=error)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, maxco=ldwork )
    ldab = ldwork
    ALLOCATE(cosab(ldab,ldab),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(sinab(ldab,ldab),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(work(ldwork,ldwork),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    IF(calculate_forces)THEN
       NULLIFY(dcosmat,dsinmat)
       ALLOCATE(dcosab(ldab,ldab,3),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(dsinab(ldab,ldab,3),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    ALLOCATE (basis_set_list(nkind),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    IF(calculate_forces)THEN
       CALL cp_dbcsr_allocate_matrix_set(dcosmat,3,error=error)
       CALL cp_dbcsr_allocate_matrix_set(dsinmat,3,error=error)

       DO m=1,3
          cgfsym = cgf_symbol(1,indco(1:3,m+1))
          name = TRIM(cgfsym(4:))//" DERIVATIVE OF SINMAT"
          CALL compress(name)
          CALL uppercase(name)
          ALLOCATE(dsinmat(m)%matrix)
          CALL cp_dbcsr_init(dsinmat(m)%matrix,error=error)
          CALL cp_dbcsr_create(matrix=dsinmat(m)%matrix, template=matrix_s(1)%matrix, &
               name=TRIM(name), matrix_type=dbcsr_type_no_symmetry, error=error)
          CALL cp_dbcsr_alloc_block_from_nbl(dsinmat(m)%matrix,sab_all,error=error)
          CALL cp_dbcsr_set(dsinmat(m)%matrix,0.0_dp,error=error)


          name = TRIM(cgfsym(4:))//" DERIVATIVE OF COSMAT"
          CALL compress(name)
          CALL uppercase(name)
          ALLOCATE(dcosmat(m)%matrix)
          CALL cp_dbcsr_init(dcosmat(m)%matrix,error=error)
          CALL cp_dbcsr_create(matrix=dcosmat(m)%matrix, template=matrix_s(1)%matrix, &
               name=TRIM(name), matrix_type=dbcsr_type_no_symmetry, error=error)
          CALL cp_dbcsr_alloc_block_from_nbl(dcosmat(m)%matrix,sab_all,error=error)
          CALL cp_dbcsr_set(dcosmat(m)%matrix,0.0_dp,error=error)

       ENDDO
    END IF

    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,orb_basis_set=basis_set_a)
       IF (ASSOCIATED(basis_set_a)) THEN
          basis_set_list(ikind)%gto_basis_set => basis_set_a
       ELSE
          NULLIFY(basis_set_list(ikind)%gto_basis_set)
       END IF
    END DO

    CALL neighbor_list_iterator_create(nl_iterator,sab_all)
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
       irow = iatom
       icol = jatom

       NULLIFY (cblock)
       CALL cp_dbcsr_get_block_p(matrix=cosmat,&
            row=irow,col=icol,BLOCK=cblock,found=found)
       NULLIFY (sblock)
       CALL cp_dbcsr_get_block_p(matrix=sinmat,&
            row=irow,col=icol,BLOCK=sblock,found=found)

       IF (calculate_forces) THEN
          irow = iatom
          icol = jatom
          DO m=1,3
             NULLIFY (sint(m)%block)
             CALL cp_dbcsr_get_block_p(matrix=dsinmat(m)%matrix,&
                  row=irow,col=icol,BLOCK=sint(m)%block,found=found)
             CPPostcondition(found,cp_failure_level,routineP,error,failure)
             NULLIFY (cost(m)%block)
             CALL cp_dbcsr_get_block_p(matrix=dcosmat(m)%matrix,&
                  row=irow,col=icol,BLOCK=cost(m)%block,found=found)
             CPPostcondition(found,cp_failure_level,routineP,error,failure)
          ENDDO
       ENDIF

       IF(ASSOCIATED(cblock).AND..NOT.ASSOCIATED(sblock).OR.&
            .NOT.ASSOCIATED(cblock).AND.ASSOCIATED(sblock)) THEN
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       ENDIF

          ra(:) = pbc(particle_set(iatom)%r(:),cell)
          rb(:) = ra+rab
          dab = SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))

          DO iset=1,nseta

             ncoa = npgfa(iset)*ncoset(la_max(iset))
             sgfa = first_sgfa(1,iset)

             DO jset=1,nsetb

                IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

                ncob = npgfb(jset)*ncoset(lb_max(jset))
                sgfb = first_sgfb(1,jset)

                !               *** Calculate the primitive integrals ***
                IF(calculate_forces)THEN
                   CALL cossin(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                        lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                        ra,rb,kvec,cosab,sinab,dcosab,dsinab)

                ELSE
                   CALL cossin(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                        lb_max(jset),npgfb(jset),zetb(:,jset),rpgfb(:,jset),lb_min(jset),&
                        ra,rb,kvec,cosab,sinab)
                END IF
                IF(ASSOCIATED(cblock).AND.ASSOCIATED(sblock)) THEN
                   CALL contract_all(cblock,sblock,&
                        iatom,ncoa,nsgfa(iset),sgfa,sphi_a,ldsa,&
                        jatom,ncob,nsgfb(jset),sgfb,sphi_b,ldsb,&
                        cosab,sinab,ldab,work,ldwork)
                ENDIF
                IF(calculate_forces)THEN
                   DO m=1,3
                      CALL contract_all(cost(m)%block,sint(m)%block,&
                           iatom,ncoa,nsgfa(iset),sgfa,sphi_a,ldsa,&
                           jatom,ncob,nsgfb(jset),sgfb,sphi_b,ldsb,&
                           dcosab(:,:,m),dsinab(:,:,m),ldab,work,ldwork)

                   END DO
                END IF
             END DO
          END DO

    END DO
    CALL neighbor_list_iterator_release(nl_iterator)

    IF(calculate_forces)THEN
       DO j = 1,3
          CALL cp_dbcsr_finalize(dcosmat(j)%matrix,error=error)
          CALL cp_dbcsr_finalize(dsinmat(j)%matrix,error=error)
       ENDDO
    END IF
    DEALLOCATE (basis_set_list,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(cosab,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(sinab,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(work,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    IF(calculate_forces)THEN
       DEALLOCATE(dcosab,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(dsinab,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

    CALL timestop(handle)

  END SUBROUTINE berry_mat_derivs

  SUBROUTINE contract_all(cos_block, sin_block,&
       iatom,ncoa,nsgfa,sgfa,sphi_a,ldsa,&
       jatom,ncob,nsgfb,sgfb,sphi_b,ldsb,&
       cosab,sinab,ldab,work,ldwork)

    REAL(dp), DIMENSION(:, :), POINTER       :: cos_block, sin_block
    INTEGER, INTENT(IN)                      :: iatom, ncoa, nsgfa, sgfa
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: sphi_a
    INTEGER, INTENT(IN)                      :: ldsa, jatom, ncob, nsgfb, sgfb
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: sphi_b
    INTEGER, INTENT(IN)                      :: ldsb
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: cosab, sinab
    INTEGER, INTENT(IN)                      :: ldab
    REAL(dp), DIMENSION(:, :)                :: work
    INTEGER, INTENT(IN)                      :: ldwork

! Calculate cosine

    CALL dgemm("N","N",ncoa,nsgfb,ncob,&
         1.0_dp,cosab(1,1),ldab,&
         sphi_b(1,sgfb),ldsb,&
         0.0_dp,work(1,1),ldwork)

    CALL dgemm("T","N",nsgfa,nsgfb,ncoa,&
        1.0_dp,sphi_a(1,sgfa),ldsa,&
        work(1,1),ldwork,&
        1.0_dp,cos_block(sgfa,sgfb),&
        SIZE(cos_block,1))

    ! Calculate sine
    CALL dgemm("N","N",ncoa,nsgfb,ncob,&
         1.0_dp,sinab(1,1),ldab,&
         sphi_b(1,sgfb),ldsb,&
         0.0_dp,work(1,1),ldwork)

    CALL dgemm("T","N",nsgfa,nsgfb,ncoa,&
         1.0_dp,sphi_a(1,sgfa),ldsa,&
         work(1,1),ldwork,&
         1.0_dp,sin_block(sgfa,sgfb),&
         SIZE(sin_block,1))

  END SUBROUTINE contract_all

END MODULE qs_efield_berry
