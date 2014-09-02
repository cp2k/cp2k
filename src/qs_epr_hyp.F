!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculates hyperfine values
!> \par History
!>      created 04-2006 [RD]
!>      adapted 02-2007 [JGH]
!> \author R. Declerck (Reinout.Declerck@UGent.be)
! *****************************************************************************
MODULE qs_epr_hyp
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_methods,                    ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_output_handling,              ONLY: cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fourpi
  USE message_passing,                 ONLY: mp_sum,&
                                             mp_sync
  USE particle_types,                  ONLY: particle_type
  USE periodic_table,                  ONLY: ptable
  USE physcon,                         ONLY: a_bohr,&
                                             a_fine,&
                                             e_charge,&
                                             e_gfactor,&
                                             e_mass,&
                                             h_bar,&
                                             mu_perm
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_dr2_gg,&
                                             pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE qs_harmonics_atom,               ONLY: harmonics_atom_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                             rho_atom_coeff,&
                                             rho_atom_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE util,                            ONLY: get_limit
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: qs_epr_hyp_calc

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_epr_hyp'

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE qs_epr_hyp_calc(qs_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_epr_hyp_calc', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    INTEGER :: bo(2), group, ia, iat, iatom, idir1, idir2, ig, ikind, ir, &
      iso, istat, jatom, mepos, natom, natomkind, nkind, num_pe, output_unit, &
      z
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, lsd, paw_atom
    REAL(dp)                                 :: arg, hard_radius, &
                                                hypanisotemp, hypfactor, &
                                                int_radius, rab2, rtemp, sum
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: hypiso, hypiso_one
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: hypaniso
    REAL(KIND=dp), DIMENSION(3)              :: ra, rab
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_p_type), POINTER                 :: hypaniso_gspace, &
                                                rhototspin_elec_gspace
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: rho_rad_h, rho_rad_s
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set
    TYPE(rho_atom_type), POINTER             :: rho_atom
    TYPE(section_vals_type), POINTER         :: dft_section

    NULLIFY(pw_env, cell, atomic_kind_set, qs_kind_set, auxbas_pw_pool, dft_control, &
         logger, dft_section, para_env, particle_set, rho, rho_atom, &
         rho_atom_set, rhototspin_elec_gspace, hypaniso_gspace)

    logger => cp_error_get_logger(error)
    dft_section => section_vals_get_subs_vals(qs_env%input,"DFT",error=error)
    output_unit = cp_print_key_unit_nr(logger,dft_section,&
                  "PRINT%HYPERFINE_COUPLING_TENSOR",&
                  extension=".eprhyp",log_filename=.FALSE.,error=error)
    CALL section_vals_val_get(dft_section,&
                "PRINT%HYPERFINE_COUPLING_TENSOR%INTERACTION_RADIUS",&
                r_val=int_radius,error=error)
    CALL section_vals_val_get(dft_section,"LSD",l_val=lsd,error=error)

    IF ( .NOT. lsd ) THEN
      ! EPR calculation only for LSD
      IF (output_unit>0) THEN
        WRITE (UNIT=output_unit,FMT="(/,T2,A)") &
            "Calculation of EPR hyperfine coupling tensors only for LSD"
      END IF
      NULLIFY(logger,dft_section)
      RETURN
    END IF

    hypfactor = -1.0_dp * mu_perm * e_charge * h_bar * e_gfactor / ( 2.0_dp * e_mass * a_bohr ** 3 )

    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,cell=cell,&
                    rho=rho,atomic_kind_set=atomic_kind_set,qs_kind_set=qs_kind_set,&
                    rho_atom_set=rho_atom_set,pw_env=pw_env,&
                    particle_set=particle_set,para_env=para_env,&
                    error=error)
    group = para_env%group

    IF (output_unit>0) THEN
      WRITE (UNIT=output_unit,FMT="(/,T2,A,/,T2,A)") &
            "Calculation of EPR hyperfine coupling tensors",&
            REPEAT("-",79)
    END IF

    ! allocate hyperfine matrices
    natom = SIZE(particle_set,1)
    ALLOCATE (hypaniso(3,3,natom),STAT=istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (hypiso(natom),STAT=istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (hypiso_one(natom),STAT=istat )
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ! set the matrices to zero
    hypiso = 0.0_dp
    hypiso_one = 0.0_dp
    hypaniso = 0.0_dp

    nkind = SIZE(atomic_kind_set) ! nkind = number of atom types

    DO ikind = 1,nkind ! loop over atom types
      NULLIFY(atom_list, grid_atom, harmonics)
      CALL get_atomic_kind(atomic_kind_set(ikind),&
                           atom_list=atom_list,natom=natomkind,z=z)

      CALL get_qs_kind(qs_kind_set(ikind),harmonics=harmonics,&
                       grid_atom=grid_atom,paw_atom=paw_atom,hard_radius=hard_radius)

      IF(.NOT. paw_atom) CYCLE ! skip the rest and go to next atom type

      num_pe = para_env%num_pe
      mepos  = para_env%mepos
      bo = get_limit( natomkind, num_pe, mepos )

      DO iat = bo(1),bo(2) ! natomkind = # atoms for ikind
        iatom = atom_list(iat)
        rho_atom => rho_atom_set(iatom)
        NULLIFY(rho_rad_h,rho_rad_s)
        CALL get_rho_atom(rho_atom=rho_atom, rho_rad_h=rho_rad_h,&
                          rho_rad_s=rho_rad_s)
        ! Non-relativistic isotropic hyperfine value (hypiso_one)
        DO ia = 1,grid_atom%ng_sphere
           DO iso = 1,harmonics%max_iso_not0
                hypiso_one(iatom) = hypiso_one(iatom) + &
                          (rho_rad_h(1)%r_coef(grid_atom%nr,iso) - &
                          rho_rad_h(2)%r_coef(grid_atom%nr,iso))*&
                          harmonics%slm(ia,iso)*grid_atom%wa(ia)/fourpi
           END DO
        END DO
        ! First calculate hard-soft contributions for the own nucleus
        ! + scalar relativistic isotropic hyperfine value (hypiso)
        DO ir = 1,grid_atom%nr
          IF (grid_atom%rad(ir) <= hard_radius) THEN
            DO ia = 1,grid_atom%ng_sphere
              hypanisotemp = 0.0_dp
              DO iso = 1,harmonics%max_iso_not0
                hypiso(iatom) = hypiso(iatom) + &
                    (rho_rad_h(1)%r_coef(ir,iso) - rho_rad_h(2)%r_coef(ir,iso))*&
                    harmonics%slm(ia,iso)*grid_atom%wr(ir)*grid_atom%wa(ia)*&
                    2._dp/(REAL(z, KIND=dp)*a_fine**2*(1._dp+2._dp*grid_atom%rad(ir)/&
                    (REAL(z, KIND=dp)*a_fine**2))**2*fourpi*grid_atom%rad(ir)**2)
                hypanisotemp = hypanisotemp + &
                    (rho_rad_h(1)%r_coef(ir,iso) - rho_rad_h(2)%r_coef(ir,iso)&
                     -(rho_rad_s(1)%r_coef(ir,iso) - rho_rad_s(2)%r_coef(ir,iso)))*&
                    harmonics%slm(ia,iso)*grid_atom%wr(ir)*grid_atom%wa(ia)/&
                    grid_atom%rad(ir)**3
              END DO  ! iso
              hypaniso(1,1,iatom) = hypaniso(1,1,iatom) + hypanisotemp*&
                       (3._dp*grid_atom%sin_pol(ia)*grid_atom%cos_azi(ia)*&
                       grid_atom%sin_pol(ia)*grid_atom%cos_azi(ia) - 1._dp)
              hypaniso(1,2,iatom) = hypaniso(1,2,iatom) + hypanisotemp*&
                       (3._dp*grid_atom%sin_pol(ia)*grid_atom%cos_azi(ia)*&
                       grid_atom%sin_pol(ia)*grid_atom%sin_azi(ia) - 0._dp)
              hypaniso(1,3,iatom) = hypaniso(1,3,iatom) + hypanisotemp*&
                       (3._dp*grid_atom%sin_pol(ia)*grid_atom%cos_azi(ia)*&
                       grid_atom%cos_pol(ia) - 0._dp)
              hypaniso(2,2,iatom) = hypaniso(2,2,iatom) + hypanisotemp*&
                       (3._dp*grid_atom%sin_pol(ia)*grid_atom%sin_azi(ia)*&
                       grid_atom%sin_pol(ia)*grid_atom%sin_azi(ia) - 1._dp)
              hypaniso(2,3,iatom) = hypaniso(2,3,iatom) + hypanisotemp*&
                       (3._dp*grid_atom%sin_pol(ia)*grid_atom%sin_azi(ia)*&
                       grid_atom%cos_pol(ia) - 0._dp)
              hypaniso(3,3,iatom) = hypaniso(3,3,iatom) + hypanisotemp*&
                       (3._dp*grid_atom%cos_pol(ia)*&
                       grid_atom%cos_pol(ia) - 1._dp)
            END DO  ! ia
          END IF ! hard_radius
        END DO  ! ir

        ! Now calculate hard-soft anisotropic contributions for the other nuclei
        DO jatom=1,natom
          IF (jatom .EQ. iatom) CYCLE ! iatom already done
          rab = pbc(particle_set(iatom)%r,particle_set(jatom)%r,cell)
          rab2 = DOT_PRODUCT(rab,rab)
          IF (SQRT(rab2) <= int_radius) THEN
            DO ir = 1,grid_atom%nr
              IF (grid_atom%rad(ir) <= hard_radius) THEN
                DO ia = 1,grid_atom%ng_sphere
                  hypanisotemp = 0.0_dp
                  rtemp = SQRT(rab2 + grid_atom%rad(ir)**2 + 2.0_dp*grid_atom%rad(ir)*&
                            (rab(1)*grid_atom%sin_pol(ia)*grid_atom%cos_azi(ia) + &
                             rab(2)*grid_atom%sin_pol(ia)*grid_atom%sin_azi(ia) + &
                             rab(3)*grid_atom%cos_pol(ia) ))
                  DO iso = 1,harmonics%max_iso_not0
                    hypanisotemp = hypanisotemp + &
                        (rho_rad_h(1)%r_coef(ir,iso) - rho_rad_h(2)%r_coef(ir,iso)&
                        -(rho_rad_s(1)%r_coef(ir,iso) - rho_rad_s(2)%r_coef(ir,iso)))*&
                        harmonics%slm(ia,iso)*grid_atom%wr(ir)*grid_atom%wa(ia)/&
                        rtemp**5
                  END DO  ! iso
                  hypaniso(1,1,jatom) = hypaniso(1,1,jatom) + hypanisotemp*&
                     (3._dp*(rab(1) + grid_atom%rad(ir)*grid_atom%sin_pol(ia)*grid_atom%cos_azi(ia))*&
                     (rab(1) + grid_atom%rad(ir)*grid_atom%sin_pol(ia)*grid_atom%cos_azi(ia)) - rtemp**2)
                  hypaniso(1,2,jatom) = hypaniso(1,2,jatom) + hypanisotemp*&
                     (3._dp*(rab(1) + grid_atom%rad(ir)*grid_atom%sin_pol(ia)*grid_atom%cos_azi(ia))*&
                     (rab(2) + grid_atom%rad(ir)*grid_atom%sin_pol(ia)*grid_atom%sin_azi(ia)) - 0._dp)
                  hypaniso(1,3,jatom) = hypaniso(1,3,jatom) + hypanisotemp*&
                     (3._dp*(rab(1) + grid_atom%rad(ir)*grid_atom%sin_pol(ia)*grid_atom%cos_azi(ia))*&
                     (rab(3) + grid_atom%rad(ir)*grid_atom%cos_pol(ia)) - 0._dp)
                  hypaniso(2,2,jatom) = hypaniso(2,2,jatom) + hypanisotemp*&
                     (3._dp*(rab(2) + grid_atom%rad(ir)*grid_atom%sin_pol(ia)*grid_atom%sin_azi(ia))*&
                     (rab(2) + grid_atom%rad(ir)*grid_atom%sin_pol(ia)*grid_atom%sin_azi(ia)) - rtemp**2)
                  hypaniso(2,3,jatom) = hypaniso(2,3,jatom) + hypanisotemp*&
                     (3._dp*(rab(2) + grid_atom%rad(ir)*grid_atom%sin_pol(ia)*grid_atom%sin_azi(ia))*&
                     (rab(3) + grid_atom%rad(ir)*grid_atom%cos_pol(ia)) - 0._dp)
                  hypaniso(3,3,jatom) = hypaniso(3,3,jatom) + hypanisotemp*&
                     (3._dp*(rab(3) + grid_atom%rad(ir)*grid_atom%cos_pol(ia))*&
                     (rab(3) + grid_atom%rad(ir)*grid_atom%cos_pol(ia)) - rtemp**2)
                END DO  ! ia
              END IF ! hard_radius
            END DO  ! ir
          END IF ! rab2
        END DO ! jatom
      END DO ! iat
    END DO ! ikind

    ! Now calculate the soft electronic spin density in reciprocal space (g-space)
    ! Plane waves grid to assemble the soft electronic spin density
    CALL pw_env_get(pw_env=pw_env,&
                    auxbas_pw_pool=auxbas_pw_pool,&
                    error=error)

    ALLOCATE (rhototspin_elec_gspace,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL pw_pool_create_pw(auxbas_pw_pool,&
                            rhototspin_elec_gspace%pw,&
                            use_data=COMPLEXDATA1D,&
                            in_space=RECIPROCALSPACE,&
                            error=error)
    CALL pw_zero(rhototspin_elec_gspace%pw,error=error)

    pw_grid => rhototspin_elec_gspace%pw%pw_grid

    ! Load the contribution of the soft electronic density
    CPPrecondition(SIZE(rho%rho_g)>1,cp_failure_level,routineP,error,failure)
    CALL pw_axpy(rho%rho_g(1)%pw,rhototspin_elec_gspace%pw,error=error)
    CALL pw_axpy(rho%rho_g(2)%pw,rhototspin_elec_gspace%pw,alpha=-1._dp,error=error)
    ! grid to assemble anisotropic hyperfine terms
    ALLOCATE (hypaniso_gspace,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    CALL pw_pool_create_pw(auxbas_pw_pool,&
                            hypaniso_gspace%pw,&
                            use_data=COMPLEXDATA1D,&
                            in_space=RECIPROCALSPACE,&
                            error=error)

    DO idir1=1,3
      DO idir2=idir1,3 ! tensor symmetry
        CALL pw_zero(hypaniso_gspace%pw,error=error)
        CALL pw_dr2_gg ( rhototspin_elec_gspace%pw, hypaniso_gspace%pw, &
                         idir1, idir2,error=error)
        DO iatom=1,natom
          sum = 0.0_dp
          ra(:) = pbc(particle_set(iatom)%r,cell)
          DO ig = 1, SIZE ( hypaniso_gspace%pw % cc )
            arg = DOT_PRODUCT(pw_grid % g (:,ig),ra)
            sum = sum + COS(arg)*REAL(hypaniso_gspace%pw%cc(ig),dp) &
                      - SIN(arg)*AIMAG(hypaniso_gspace%pw%cc(ig))
          END DO
          ! Actually, we need -1.0 * fourpi * hypaniso_gspace
          sum = sum * fourpi * (-1.0_dp)
          hypaniso(idir1,idir2,iatom) = hypaniso(idir1,idir2,iatom) + sum
        END DO
      END DO ! idir2
    END DO ! idir1

    CALL pw_pool_give_back_pw(auxbas_pw_pool,rhototspin_elec_gspace%pw,error=error)
    DEALLOCATE (rhototspin_elec_gspace,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL pw_pool_give_back_pw(auxbas_pw_pool,hypaniso_gspace%pw,error=error)
    DEALLOCATE (hypaniso_gspace,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ! Multiply hyperfine matrices with constant*gyromagnetic ratio's
    ! to have it in units of Mhz.

    DO iatom=1,natom
      CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                           z=z)
      hypiso(iatom) = hypiso(iatom) *&
                      2.0_dp/3.0_dp*hypfactor * ptable(z)%gyrom_ratio
      hypiso_one(iatom) = hypiso_one(iatom) *&
                      2.0_dp/3.0_dp*hypfactor * ptable(z)%gyrom_ratio
      DO idir1=1,3
        DO idir2=idir1,3
          hypaniso(idir1,idir2,iatom) = hypaniso(idir1,idir2,iatom) *&
                                        hypfactor/fourpi * ptable(z)%gyrom_ratio
          IF (idir1 /= idir2) THEN
            hypaniso(idir2,idir1,iatom) = hypaniso(idir1,idir2,iatom)
          END IF
        END DO
      END DO
    END DO

    ! Global sum
    CALL mp_sync(group)
    CALL mp_sum(hypaniso,group)
    CALL mp_sum(hypiso,group)
    CALL mp_sum(hypiso_one,group)

    ! Print hyperfine matrices
    IF (output_unit>0) THEN
      DO iatom=1,natom
        CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                             element_symbol=element_symbol,z=z)
        WRITE (UNIT=output_unit,FMT="(T1,I5,T7,A,T10,I3,T14,F16.10,T31,A,T60,F20.10)")&
          iatom,element_symbol, ptable(z)%gyrom_ratio_isotope, ptable(z)%gyrom_ratio,&
          "[Mhz/T]  Sca-Rel A_iso [Mhz]", hypiso(iatom)
        WRITE (UNIT=output_unit,FMT="(T31,A,T60,F20.10)")&
          "         Non-Rel A_iso [Mhz]", hypiso_one(iatom)
        WRITE (UNIT=output_unit,FMT="(T4,A,T18,F20.10,1X,F20.10,1X,F20.10)")&
          "             ", hypaniso(1,1,iatom), hypaniso(1,2,iatom), hypaniso(1,3,iatom),&
          "  A_ani [Mhz]", hypaniso(2,1,iatom), hypaniso(2,2,iatom), hypaniso(2,3,iatom),&
          "             ", hypaniso(3,1,iatom), hypaniso(3,2,iatom), hypaniso(3,3,iatom)
      ENDDO
    END IF

    ! Deallocate the remainings ...
    DEALLOCATE (hypiso,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (hypiso_one,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (hypaniso,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE qs_epr_hyp_calc

END MODULE qs_epr_hyp

