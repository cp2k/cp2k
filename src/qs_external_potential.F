!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines to handle an external electrostatic field
!>        The external field can be generic and is provided by user input
! *****************************************************************************
MODULE qs_external_potential
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
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
  USE message_passing,                 ONLY: mp_bcast
  USE particle_types,                  ONLY: particle_type
  USE pw_grid_types,                   ONLY: PW_MODE_LOCAL
  USE pw_types,                        ONLY: pw_p_type,&
                                             pw_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE realspace_grid_cube,             ONLY: cube_to_pw
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
    LOGICAL                                  :: failure, read_from_cube, &
                                                static_potential
    REAL(kind=dp)                            :: dvol, efunc, scaling_factor
    REAL(kind=dp), DIMENSION(3)              :: dr, grid_p
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_p_type), POINTER                 :: v_ee
    TYPE(section_vals_type), POINTER         :: ext_pot_section, input

    failure=.FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(v_ee, input, ext_pot_section, dft_control)
    CALL get_qs_env(qs_env,&
                    vee=v_ee,&
                    input=input,&
                    dft_control=dft_control,&
                    error=error)
    IF(dft_control%apply_external_potential)THEN
       ext_pot_section => section_vals_get_subs_vals(input,"DFT%EXTERNAL_POTENTIAL",error=error)
       CALL section_vals_val_get(ext_pot_section,"STATIC",l_val=static_potential,error=error)
       CALL section_vals_val_get(ext_pot_section,"READ_FROM_CUBE",l_val=read_from_cube,error=error)

       IF ((.NOT.static_potential).OR.dft_control%eval_external_potential) THEN
          IF (read_from_cube) THEN
             CALL section_vals_val_get(ext_pot_section,"SCALING_FACTOR",r_val=scaling_factor,error=error)
             CALL cube_to_pw(v_ee%pw, 'pot.cube', scaling_factor, error)
             dft_control%eval_external_potential = .FALSE.
          ELSE
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
             dft_control%eval_external_potential = .FALSE.
          END IF
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
    LOGICAL                                  :: failure, my_force, pot_on_grid
    REAL(KIND=dp)                            :: ee_core_ener, efunc, zeff
    REAL(KIND=dp), DIMENSION(3)              :: dfunc, r
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_p_type), POINTER                 :: v_ee
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(section_vals_type), POINTER         :: ext_pot_section, input

    failure=.FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(dft_control)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set,&
                    energy=energy,&
                    particle_set=particle_set,&
                    input=input,&
                    cell=cell,&
                    dft_control=dft_control,&
                    error=error)

    IF(dft_control%apply_external_potential)THEN
       IF (dft_control%eval_external_potential) THEN !ensure that external potential is loaded to grid
          CALL external_e_potential(qs_env,error)
       END IF
       my_force=.FALSE.
       IF(PRESENT(calculate_forces)) my_force=calculate_forces
       ext_pot_section => section_vals_get_subs_vals(input,"DFT%EXTERNAL_POTENTIAL",error=error)
       ee_core_ener=0.0_dp
       nkind = SIZE(atomic_kind_set)
       
       !check if external potential on grid has been loaded from a file instead of giving a function
       CALL section_vals_val_get(ext_pot_section,"READ_FROM_CUBE",l_val=pot_on_grid,error=error)
       IF(pot_on_grid) CALL get_qs_env(qs_env,vee=v_ee,input=input,error=error)

       DO ikind=1,SIZE(atomic_kind_set)
          CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=list,natom=natom)
          CALL get_qs_kind(qs_kind_set(ikind), zeff=zeff)

          natom = SIZE(list)
          DO iatom=1,natom
             atom_a = list(iatom)
             !pbc returns r(i) in range [-cell%hmat(i,i)/2, cell%hmat(i,i)/2]
             !for periodic dimensions (assuming the cell is orthorombic).
             !This is not consistent with the potential on grid, where r(i) is
             !in range [0, cell%hmat(i,i)]
             !Use new pbc function with switch positive_range=.TRUE.
             r(:) = pbc(particle_set(atom_a)%r(:), cell, positive_range=.TRUE.)
             !if potential is on grid, interpolate the value at r,
             !otherwise evaluate the given function
             IF(pot_on_grid) THEN
                CALL interpolate_external_potential(r, v_ee%pw, func=efunc, &
                      dfunc=dfunc, calc_derivatives=my_force, error=error)
             ELSE
                CALL get_external_potential(r, ext_pot_section, func=efunc, &
                      dfunc=dfunc, calc_derivatives=my_force, error=error)
             END IF
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
!> \param r                position in realspace
!> \param ext_pot_section ...
!> \param func             external potential at r
!> \param dfunc            derivative of the external potential at r
!> \param calc_derivatives Whether to calulate dfunc
!> \param error ...
!> \date   12.2009
!> \par History
!>      12.2009            created [tlaino] 
!>      11.2014            reading external cube file added [Juha Ritala & Matt Watkins]
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
    
! *****************************************************************************
!> \brief                  subroutine that interpolates the value of the external 
!>                         potential at position r based on the values on the realspace grid
!> \param r                 ...
!> \param grid             external potential pw grid, vee
!> \param func             value of vee at r
!> \param dfunc            derivatives of vee at r
!> \param calc_derivatives calc dfunc
!> \param error ...
! *****************************************************************************
  SUBROUTINE interpolate_external_potential(r, grid, func, dfunc, calc_derivatives, error)
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: r
    TYPE(pw_type), POINTER                   :: grid
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: func, dfunc(3)
    LOGICAL, INTENT(IN), OPTIONAL            :: calc_derivatives
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'interpolate_external_potential', &
      routineP = moduleN//':'//routineN

    INTEGER :: buffer_i, buffer_j, buffer_k, data_source, fd_extra_point, &
      gid, handle, i, i_pbc, ip, j, j_pbc, k, k_pbc, my_rank, num_pe, tag
    INTEGER, DIMENSION(3)                    :: lbounds, lbounds_local, &
                                                lower_inds, ubounds, &
                                                ubounds_local, upper_inds
    LOGICAL                                  :: check, failure, my_force
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: bcast_buffer
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: grid_buffer
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: dgrid
    REAL(KIND=dp), DIMENSION(3)              :: dr, subgrid_origin

    failure = .FALSE.
    CALL timeset(routineN,handle)
    my_force = .FALSE.
    IF (PRESENT(calc_derivatives)) my_force=calc_derivatives
    check=PRESENT(dfunc).EQV.PRESENT(calc_derivatives)
    CPPrecondition(check,cp_failure_level,routineP,error,failure)
    
    IF (my_force) THEN
      ALLOCATE(grid_buffer(0:3,0:3,0:3))
      ALLOCATE(bcast_buffer(0:3))
      ALLOCATE(dgrid(1:2,1:2,1:2,3))
      fd_extra_point = 1
    ELSE
      ALLOCATE(grid_buffer(1:2,1:2,1:2))
      ALLOCATE(bcast_buffer(1:2))
      fd_extra_point = 0
    END IF
    
    ! The values of external potential on grid are distributed among the
    ! processes, so first we have to gather them up
    gid     = grid%pw_grid%para%group
    my_rank = grid%pw_grid%para%my_pos
    num_pe  = grid%pw_grid%para%group_size   
    tag     = 1
    
    dr = grid%pw_grid%dr
    lbounds = grid%pw_grid%bounds(1,:)
    ubounds = grid%pw_grid%bounds(2,:)
    lbounds_local = grid%pw_grid%bounds_local(1,:)
    ubounds_local = grid%pw_grid%bounds_local(2,:)

    ! Determine the indices of grid points that are needed
    lower_inds = lbounds+FLOOR(r/dr)-fd_extra_point
    upper_inds = lower_inds+1+2*fd_extra_point
    
    DO i = lower_inds(1), upper_inds(1)
      ! If index is out of global bounds, assume periodic boundary conditions
      i_pbc = pbc_index(i,lbounds(1),ubounds(1))
      buffer_i = i-lower_inds(1)+1-fd_extra_point
      DO j = lower_inds(2), upper_inds(2)
        j_pbc = pbc_index(j,lbounds(2),ubounds(2))
        buffer_j = j-lower_inds(2)+1-fd_extra_point
        
        ! Find the process that has the data for indices i_pbc and j_pbc 
        ! and store the data to bcast_buffer. Assuming that each process has full z data
        IF (grid%pw_grid%para%mode .NE. PW_MODE_LOCAL ) THEN
           DO ip = 0, num_pe-1
              IF (grid%pw_grid%para%bo (1,1,ip,1)<=i_pbc-lbounds(1)+1 .AND. &
                  grid%pw_grid%para%bo (2,1,ip,1)>=i_pbc-lbounds(1)+1 .AND. &
                  grid%pw_grid%para%bo (1,2,ip,1)<=j_pbc-lbounds(2)+1 .AND. &
                  grid%pw_grid%para%bo (2,2,ip,1)>=j_pbc-lbounds(2)+1 ) THEN
                 data_source = ip
                 EXIT
              END IF
           END DO
           IF (my_rank == data_source) THEN
              IF (lower_inds(3) >= lbounds(3) .AND. upper_inds(3) <= ubounds(3)) THEN
                 bcast_buffer(:) = &
                      grid%cr3d(i_pbc,j_pbc,lower_inds(3):upper_inds(3))
              ELSE
                 DO k = lower_inds(3), upper_inds(3)
                    k_pbc = pbc_index(k,lbounds(3),ubounds(3))
                    buffer_k = k-lower_inds(3)+1-fd_extra_point
                    bcast_buffer(buffer_k) = &
                         grid%cr3d(i_pbc,j_pbc,k_pbc)
                 END DO
              END IF
           END IF
           ! data_source sends data to everyone
           CALL mp_bcast(bcast_buffer,data_source,gid)
           grid_buffer(buffer_i,buffer_j,:) = bcast_buffer
        ELSE
           grid_buffer(buffer_i,buffer_j,:) = grid%cr3d(i_pbc,j_pbc,lower_inds(3):upper_inds(3))
        END IF
      END DO
    END DO
    
    ! Now that all the processes have local external potential data around r,
    ! interpolate the value at r
    subgrid_origin = (lower_inds-lbounds+fd_extra_point)*dr
    func = trilinear_interpolation(r,grid_buffer(1:2,1:2,1:2),subgrid_origin,dr)
    
    ! If the derivative of the potential is needed, approximate the derivative at grid
    ! points using finite differences, and then interpolate the value at r
    IF (my_force) THEN
      CALL d_finite_difference(grid_buffer,dr,dgrid)
      DO i = 1, 3
        dfunc(i) = trilinear_interpolation(r,dgrid(:,:,:,i),subgrid_origin,dr)
      END DO
      DEALLOCATE(dgrid)
    END IF
    
    DEALLOCATE(grid_buffer)
    CALL timestop(handle)
  END SUBROUTINE interpolate_external_potential
  
! *****************************************************************************
!> \brief       subroutine that uses finite differences to approximate the partial
!>              derivatives of the potential based on the given values on grid
!> \param grid  tiny bit of external potential vee
!> \param dr    step size for finite difference
!> \param dgrid derivatives of grid
! *****************************************************************************
  SUBROUTINE d_finite_difference(grid,dr,dgrid)
    REAL(KIND=dp), DIMENSION(0:, 0:, 0:), &
      INTENT(IN)                             :: grid
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: dr
    REAL(KIND=dp), &
      DIMENSION(1:, 1:, 1:, :), INTENT(OUT)  :: dgrid

    INTEGER                                  :: i, j, k

    DO i = 1, SIZE(dgrid,1)
      DO j = 1, SIZE(dgrid,2)
        DO k = 1, SIZE(dgrid,3)
          dgrid(i,j,k,1) = 0.5*(grid(i+1,j,k)-grid(i-1,j,k))/dr(1)
          dgrid(i,j,k,2) = 0.5*(grid(i,j+1,k)-grid(i,j-1,k))/dr(2)
          dgrid(i,j,k,3) = 0.5*(grid(i,j,k+1)-grid(i,j,k-1))/dr(3)
        END DO
      END DO
    END DO
  END SUBROUTINE d_finite_difference
  
! *****************************************************************************
!> \brief             trilinear interpolation function that interpolates value at r based 
!>                    on 2x2x2 grid points around r in subgrid
!> \param r           where to interpolate to 
!> \param subgrid     part of external potential on a grid
!> \param origin      center of grid
!> \param dr          step size
!> \retval value_at_r interpolated value of external potential
! *****************************************************************************
  FUNCTION trilinear_interpolation(r,subgrid,origin,dr) RESULT(value_at_r)
    REAL(KIND=dp), DIMENSION(3)              :: r
    REAL(KIND=dp), DIMENSION(:, :, :)        :: subgrid
    REAL(KIND=dp), DIMENSION(3)              :: origin, dr
    REAL(KIND=dp)                            :: value_at_r

    REAL(KIND=dp), DIMENSION(3)              :: norm_r, norm_r_rev

    norm_r = (r-origin)/dr
    norm_r_rev = 1-norm_r
    value_at_r = subgrid(1,1,1)*PRODUCT(norm_r_rev) + &
                  subgrid(2,1,1)*norm_r(1)*norm_r_rev(2)*norm_r_rev(3) + &
                  subgrid(1,2,1)*norm_r_rev(1)*norm_r(2)*norm_r_rev(3) + &
                  subgrid(1,1,2)*norm_r_rev(1)*norm_r_rev(2)*norm_r(3) + &
                  subgrid(1,2,2)*norm_r_rev(1)*norm_r(2)*norm_r(3) + &
                  subgrid(2,1,2)*norm_r(1)*norm_r_rev(2)*norm_r(3) + &
                  subgrid(2,2,1)*norm_r(1)*norm_r(2)*norm_r_rev(3) + &
                  subgrid(2,2,2)*PRODUCT(norm_r)
  END FUNCTION trilinear_interpolation
  
! *****************************************************************************
!> \brief          get a correct value for possible out of bounds index using periodic
!>                  boundary conditions
!> \param i ...
!> \param lowbound ...
!> \param upbound ...
!> \retval pbc_index ...
! *****************************************************************************
  FUNCTION pbc_index(i, lowbound, upbound)
    INTEGER                                  :: i, lowbound, upbound, &
                                                pbc_index

    IF (i < lowbound) THEN
        pbc_index = upbound+i-lowbound+1
      ELSE IF (i > upbound) THEN
        pbc_index = lowbound+i-upbound-1
      ELSE
        pbc_index = i
      END IF
  END FUNCTION pbc_index

END MODULE qs_external_potential
