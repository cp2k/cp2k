!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines to handle an external electrostatic field
!>        The external field can be generic and is provided by user input
! *****************************************************************************
MODULE qs_external_potential
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_methods,                    ONLY: cell_type,&
                                             pbc
  USE cp_subsys_types,                 ONLY: cp_subsys_get
  USE f77_blas
  USE force_fields_util,               ONLY: get_generic_info
  USE fparser,                         ONLY: evalf,&
                                             evalfd,&
                                             finalizef,&
                                             initf,&
                                             parsef
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE pw_types,                        ONLY: pw_p_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE string_utilities,                ONLY: compress
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_external_potential'

! *** Public subroutines ***
  PUBLIC :: external_e_potential,&
            external_c_potential

CONTAINS

! *****************************************************************************
!> \brief  Computes the external potential on the grid
!> \param qs_env ...
!> \param error ...
!> \date   12.2009
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  SUBROUTINE external_e_potential(qs_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'external_e_potential', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, k
    INTEGER, DIMENSION(2, 3)                 :: bo_global, bo_local
    LOGICAL                                  :: failure, static_potential
    REAL(kind=dp)                            :: dvol, efunc
    REAL(kind=dp), DIMENSION(3)              :: dr, grid_p
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(pw_p_type), POINTER                 :: v_ee
    TYPE(section_vals_type), POINTER         :: ext_pot_section, input

    failure=.FALSE.
    CALL timeset(routineN,handle)
    IF(qs_env%dft_control%apply_external_potential)THEN
       NULLIFY(cell, particles, v_ee, input, ext_pot_section)
       CALL cp_subsys_get(qs_env%subsys,particles=particles,error=error)
       CALL get_qs_env(qs_env,cell=cell,vee=v_ee,input=input,error=error)
       ext_pot_section => section_vals_get_subs_vals(input,"DFT%EXTERNAL_POTENTIAL",error=error)
       CALL section_vals_val_get(ext_pot_section,"STATIC",l_val=static_potential,error=error)

       IF ((.NOT.static_potential).OR.qs_env%dft_control%eval_external_potential) THEN
          dr=v_ee%pw%pw_grid%dr
          dvol=v_ee%pw%pw_grid%dvol
          v_ee%pw%cr3d=0.0_dp

          bo_local=v_ee%pw%pw_grid%bounds_local
          bo_global=v_ee%pw%pw_grid%bounds

          DO k=bo_local(1,3),bo_local(2,3)
             DO j=bo_local(1,2),bo_local(2,2)
                DO i=bo_local(1,1),bo_local(2,1)
                   grid_p(1)=(i-bo_global(1,1))*dr(1)
                   grid_p(2)=(j-bo_global(1,2))*dr(2)
                   grid_p(3)=(k-bo_global(1,3))*dr(3)
                   CALL get_external_potential(grid_p, ext_pot_section, func=efunc, error=error)
                   v_ee%pw%cr3d(i,j,k)= v_ee%pw%cr3d(i,j,k)+efunc
                END DO
             END DO
          END DO
          qs_env%dft_control%eval_external_potential = .FALSE.
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE external_e_potential

! *****************************************************************************
!> \brief  Computes the force and the energy due to the external potential on the cores
!> \param qs_env ...
!> \param calculate_forces ...
!> \param error ...
!> \date   12.2009
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  SUBROUTINE external_c_potential(qs_env,calculate_forces,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, OPTIONAL                        :: calculate_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'external_c_potential', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom_a, handle, iatom, ikind, &
                                                natom, nkind
    INTEGER, DIMENSION(:), POINTER           :: list
    LOGICAL                                  :: failure, my_force
    REAL(KIND=dp)                            :: ee_core_ener, efunc, zeff
    REAL(KIND=dp), DIMENSION(3)              :: dfunc, r
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(section_vals_type), POINTER         :: ext_pot_section, input

    failure=.FALSE.
    CALL timeset(routineN,handle)
    IF(qs_env%dft_control%apply_external_potential)THEN
       my_force=.FALSE.
       IF(PRESENT(calculate_forces)) my_force=calculate_forces

       CALL get_qs_env(qs_env=qs_env,&
                       atomic_kind_set=atomic_kind_set,&
                       qs_kind_set=qs_kind_set,&
                       energy=energy,&
                       particle_set=particle_set,&
                       input=input,&
                       cell=cell,&
                       error=error)
       ext_pot_section => section_vals_get_subs_vals(input,"DFT%EXTERNAL_POTENTIAL",error=error)
       ee_core_ener=0.0_dp
       nkind = SIZE(atomic_kind_set)

       DO ikind=1,SIZE(atomic_kind_set)
          CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=list,natom=natom)
          CALL get_qs_kind(qs_kind_set(ikind), zeff=zeff)

          natom = SIZE(list)
          DO iatom=1,natom
             atom_a = list(iatom)
             r(:) = pbc(particle_set(atom_a)%r(:),cell)
             CALL get_external_potential(r, ext_pot_section, func=efunc, &
                   dfunc=dfunc, calc_derivatives=my_force, error=error)
             ee_core_ener=ee_core_ener+zeff*efunc
             IF(my_force)THEN
                CALL get_qs_env(qs_env=qs_env,force=force,error=error)
                force(ikind)%eev(:,iatom)=dfunc*zeff
             END IF
          END DO
       END DO
       energy%ee_core=ee_core_ener
    END IF
    CALL timestop(handle)
  END SUBROUTINE external_c_potential

! *****************************************************************************
!> \brief  Low level function for computing the potential and the derivatives
!> \param r ...
!> \param ext_pot_section ...
!> \param func ...
!> \param dfunc ...
!> \param calc_derivatives ...
!> \param error ...
!> \date   12.2009
!> \author Teodoro Laino [tlaino]
! *****************************************************************************
  SUBROUTINE get_external_potential(r, ext_pot_section, func, dfunc, calc_derivatives, error)
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: r
    TYPE(section_vals_type), POINTER         :: ext_pot_section
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: func, dfunc(3)
    LOGICAL, INTENT(IN), OPTIONAL            :: calc_derivatives
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_external_potential', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: coupling_function
    CHARACTER(LEN=default_string_length)     :: def_error, this_error
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: my_par
    INTEGER                                  :: handle, j, stat
    LOGICAL                                  :: check, failure, my_force
    REAL(KIND=dp)                            :: dedf, dx, err, lerr
    REAL(KIND=dp), DIMENSION(:), POINTER     :: my_val

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(my_par, my_val)
    my_force = .FALSE.
    IF (PRESENT(calc_derivatives)) my_force=calc_derivatives
    check=PRESENT(dfunc).EQV.PRESENT(calc_derivatives)
    CPPrecondition(check,cp_failure_level,routineP,error,failure)
    CALL section_vals_val_get(ext_pot_section,"DX",r_val=dx,error=error)
    CALL section_vals_val_get(ext_pot_section,"ERROR_LIMIT",r_val=lerr,error=error)
    CALL get_generic_info(ext_pot_section, "FUNCTION", coupling_function, my_par, my_val,&
         input_variables=(/"X","Y","Z"/), i_rep_sec=1,error=error)
    CALL initf(1)
    CALL parsef(1,TRIM(coupling_function),my_par)

    my_val(1)=r(1)
    my_val(2)=r(2)
    my_val(3)=r(3)

    IF (PRESENT(func)) func = evalf(1,my_val)
    IF (my_force) THEN
       DO j = 1, 3
          dedf = evalfd(1,j,my_val,dx,err)
          IF (ABS(err)>lerr) THEN
             WRITE(this_error,"(A,G12.6,A)")"(",err,")"
             WRITE(def_error,"(A,G12.6,A)")"(",lerr,")"
             CALL compress(this_error,.TRUE.)
             CALL compress(def_error,.TRUE.)
             CALL cp_assert(.FALSE.,cp_warning_level,-300,routineP,&
                  'ASSERTION (cond) failed at line '//cp_to_string(__LINE__)//&
                  ' Error '//TRIM(this_error)//' in computing numerical derivatives larger then'//&
                  TRIM(def_error)//' .',error=error,only_ionode=.TRUE.)
          END IF
          dfunc(j)=dedf
       END DO
    END IF
    DEALLOCATE(my_par,stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(my_val,stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL finalizef()
    CALL timestop(handle)
  END SUBROUTINE get_external_potential

END MODULE qs_external_potential
