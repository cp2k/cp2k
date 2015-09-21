!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!> \author
! *****************************************************************************
MODULE gle_system_dynamics

  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE extended_system_types,           ONLY: map_info_type
  USE gle_system_types,                ONLY: gle_thermo_create,&
                                             gle_type
  USE input_constants,                 ONLY: &
       do_thermo_communication, do_thermo_no_communication, isokin_ensemble, &
       langevin_ensemble, npe_f_ensemble, npe_i_ensemble, &
       nph_uniaxial_damped_ensemble, nph_uniaxial_ensemble, npt_f_ensemble, &
       npt_i_ensemble, nve_ensemble, nvt_ensemble, reftraj_ensemble
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_remove_values,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: global_constraint_type,&
                                             molecule_type
  USE parallel_rng_types,              ONLY: next_random_number,&
                                             read_rng_stream,&
                                             rng_record_length,&
                                             rng_stream_type
  USE particle_types,                  ONLY: particle_type
  USE simpar_types,                    ONLY: simpar_type
  USE termination,                     ONLY: stop_program
  USE thermostat_mapping,              ONLY: thermostat_mapping_region
  USE thermostat_types,                ONLY: thermostat_info_type
  USE thermostat_utils,                ONLY: ke_region_particles,&
                                             momentum_region_particles,&
                                             vel_rescale_particles
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: gle_particles,&
            initialize_gle_part,&
            gle_cholesky_stab, &
            gle_matrix_exp, &
            restart_gle

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gle_system_dynamics'

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param gle ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param particle_set ...
!> \param local_molecules ...
!> \param group ...
!> \param shell_adiabatic ...
!> \param shell_particle_set ...
!> \param core_particle_set ...
!> \param vel ...
!> \param shell_vel ...
!> \param core_vel ...
!> \date
!> \par History
! *****************************************************************************
  SUBROUTINE gle_particles( gle, molecule_kind_set, molecule_set,  particle_set, local_molecules, &
       group, shell_adiabatic, shell_particle_set, core_particle_set, vel, shell_vel, core_vel)

    TYPE(gle_type), POINTER                  :: gle
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    INTEGER, INTENT(IN)                      :: group
    LOGICAL, INTENT(IN), OPTIONAL            :: shell_adiabatic
    TYPE(particle_type), OPTIONAL, POINTER   :: shell_particle_set(:), &
                                                core_particle_set(:)
    REAL(KIND=dp), INTENT(INOUT), OPTIONAL   :: vel(:,:), shell_vel(:,:), &
                                                core_vel(:,:)

    CHARACTER(len=*), PARAMETER :: routineN = 'gle_particles', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iadd, ideg, imap, &
                                                ndim, num
    LOGICAL                                  :: failure, my_shell_adiabatic, &
                                                present_vel
    REAL(dp)                                 :: alpha, beta, rr
    REAL(dp), DIMENSION(:, :), POINTER       :: a_mat, e_tmp, h_tmp, s_tmp
    TYPE(map_info_type), POINTER             :: map_info
    TYPE(rng_stream_type), POINTER           :: rng_stream

    CALL timeset(routineN,handle)
    my_shell_adiabatic = .FALSE.
    failure = .FALSE.
    IF(PRESENT(shell_adiabatic)) my_shell_adiabatic = shell_adiabatic
    NULLIFY( rng_stream)
    present_vel = PRESENT(vel)
    ndim = gle%ndim
    ALLOCATE (s_tmp(ndim,gle%loc_num_gle))
    s_tmp = 0.0_dp
    ALLOCATE (e_tmp(ndim,gle%loc_num_gle))
    ALLOCATE (h_tmp(ndim,gle%loc_num_gle))


    map_info => gle%map_info
    CALL ke_region_particles(map_info, particle_set, molecule_kind_set,&
         local_molecules, molecule_set, group, vel)
    DO ideg = 1, gle%loc_num_gle
       imap = gle%map_info%map_index(ideg)
       gle%nvt(ideg)%kin_energy = map_info%s_kin(imap)
    END DO

    CALL momentum_region_particles(map_info, particle_set, molecule_kind_set,&
         local_molecules, molecule_set, group, vel)

    DO ideg = 1,gle%loc_num_gle
       imap = gle%map_info%map_index(ideg)
       IF (gle%nvt(ideg)%nkt==0.0_dp) CYCLE
       gle%nvt(ideg)%s(1)= map_info%s_kin(imap)
       s_tmp(1,imap) = map_info%s_kin(imap)
       rng_stream => gle%nvt(ideg)%gaussian_rng_stream
       rr =  next_random_number( rng_stream)
       e_tmp(1,imap)= rr
       DO iadd = 2,ndim
          s_tmp(iadd,imap) = gle%nvt(ideg)%s(iadd)
          rr =  next_random_number( rng_stream)
          e_tmp(iadd,imap) = rr
       END DO
    END DO
    num =  gle%loc_num_gle
    a_mat => gle%gle_s
    alpha = 1.0_dp
    beta = 0.0_dp
    CALL DGEMM('N','N',ndim,num,ndim,alpha,a_mat(1,1),ndim,e_tmp(1,1),ndim,beta,h_tmp(1,1),ndim)
!
    a_mat => gle%gle_t
    beta = 1.0_dp
    CALL dgemm("N","N",ndim,num,ndim,alpha,a_mat(1,1),ndim,s_tmp(1,1),ndim,beta,h_tmp(1,1),ndim )


    DO ideg = 1,gle%loc_num_gle
       imap = gle%map_info%map_index(ideg)
       IF (gle%nvt(ideg)%nkt==0.0_dp) CYCLE

       map_info%v_scale(imap) =  h_tmp(1,imap)/s_tmp(1,imap)
       DO iadd = 2,ndim
          gle%nvt(ideg)%s(iadd) = h_tmp(iadd,ideg)
       END DO
    END DO

    CALL vel_rescale_particles(map_info,  molecule_kind_set, molecule_set, particle_set,&
         local_molecules, my_shell_adiabatic, shell_particle_set,core_particle_set,&
         vel, shell_vel, core_vel)

    CALL ke_region_particles(map_info, particle_set, molecule_kind_set,&
         local_molecules, molecule_set, group, vel)
    DO ideg = 1, gle%loc_num_gle
       imap = gle%map_info%map_index(ideg)
       gle%nvt(ideg)%thermostat_energy = gle%nvt(ideg)%thermostat_energy + &
            0.5_dp*( gle%nvt(ideg)%kin_energy - map_info%s_kin(imap))
    END DO

    DEALLOCATE(e_tmp,s_tmp,h_tmp)

    CALL timestop(handle)
  END SUBROUTINE gle_particles


! *****************************************************************************
!> \brief ...
!> \param thermostat_info ...
!> \param simpar ...
!> \param local_molecules ...
!> \param molecule ...
!> \param molecule_kind_set ...
!> \param para_env ...
!> \param gle ...
!> \param gle_section ...
!> \param gci ...
!> \param save_mem ...
!> \author
! *****************************************************************************
  SUBROUTINE initialize_gle_part ( thermostat_info, simpar, local_molecules,&
       molecule, molecule_kind_set, para_env,  gle, gle_section,&
       gci, save_mem)

    TYPE(thermostat_info_type), POINTER      :: thermostat_info
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(gle_type), POINTER                  :: gle
    TYPE(section_vals_type), POINTER         :: gle_section
    TYPE(global_constraint_type), POINTER    :: gci
    LOGICAL, INTENT(IN)                      :: save_mem

    CHARACTER(len=*), PARAMETER :: routineN = 'initialize_gle_part', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, restart
    REAL(dp)                                 :: Mtmp(gle%ndim,gle%ndim)

    failure = .FALSE.
    restart = .FALSE.

    CALL gle_to_particle_mapping ( thermostat_info, simpar, local_molecules,&
         molecule, molecule_kind_set, gle, para_env, gci)

    IF ( gle % ndim /= 0 ) THEN
            CALL init_gle_variables ( gle)
    END IF
    CALL restart_gle(gle,gle_section,save_mem,restart)

    ! here we should have read a_mat and c_mat; whe can therefore compute S and T
    ! deterministic part of the propagator
    CALL gle_matrix_exp((-simpar%dt*0.5_dp)*gle%a_mat, gle%ndim,15,15,gle%gle_t)
    ! stochastic part
    Mtmp=gle%c_mat-MATMUL(gle%gle_t,MATMUL(gle%c_mat,TRANSPOSE(gle%gle_t)))
    CALL gle_cholesky_stab(Mtmp, gle%gle_s, gle%ndim)

  END SUBROUTINE initialize_gle_part

! *****************************************************************************
!> \brief ...
!> \param M ...
!> \param n ...
!> \param j ...
!> \param k ...
!> \param EM ...
!> \author Michele
!> routine to compute matrix exponential via scale & square
! *****************************************************************************
  SUBROUTINE gle_matrix_exp(M, n, j, k, EM)

    INTEGER, INTENT(in)                      :: n
    REAL(dp), INTENT(in)                     :: M(n,n)
    INTEGER, INTENT(in)                      :: j, k
    REAL(dp), INTENT(out)                    :: EM(n,n)

    INTEGER                                  :: i, p
    REAL(dp)                                 :: SM(n,n), tc(j+1)

    tc(1)=1._dp
    DO i=1,j
       tc(i+1)=tc(i)/REAL(i,KIND=dp)
    ENDDO

    !scale
    SM=M*(1._dp/2._dp**k)
    EM=0._dp
    DO i=1,n
       EM(i,i)=tc(j+1)
    ENDDO

    !taylor exp of scaled matrix
    DO p=j,1,-1
       EM=MATMUL(SM,EM)
       DO i=1,n
          EM(i,i)=EM(i,i)+tc(p)
       ENDDO
    ENDDO

    !square
    DO p=1,k
       EM=MATMUL(EM,EM)
    ENDDO
  END SUBROUTINE gle_matrix_exp

! *****************************************************************************
!> \brief ...
!> \param SST ...
!> \param S ...
!> \param n ...
!> \author Michele
!>  "stabilized" cholesky to deal with small & negative diagonal elements,
!>  in practice a LDL^T decomposition is computed, and negative els. are zeroed
! *****************************************************************************
  SUBROUTINE gle_cholesky_stab(SST, S, n)
    INTEGER, INTENT(in)                      :: n
    REAL(dp), INTENT(out)                    :: S(n,n)
    REAL(dp), INTENT(in)                     :: SST(n,n)

    INTEGER                                  :: i, j, k
    REAL(dp)                                 :: D(n), L(n,n)

    D=0._dp
    L=0._dp
    S=0._dp
    DO i=1,n
       L(i,i) = 1.0_dp
       D(i) = SST(i,i)
       DO j=1,i-1
          L(i,j)=SST(i,j);
          DO k=1,j-1
             L(i,j)=L(i,j)-L(i,k)*L(j,k)*D(k)
          ENDDO
          IF(ABS(D(j))>EPSILON(1.0_dp)) L(i,j)=L(i,j)/D(j)
       ENDDO
       DO k=1,i-1
          D(i)=D(i)-L(i,k)*L(i,k)*D(k)
       END DO
    ENDDO
    DO i=1,n
      DO j=1,i
         IF (ABS(D(j))>EPSILON(1.0_dp)) S(i,j)=S(i,j)+L(i,j)*SQRT(D(j))
      ENDDO
    ENDDO

  END SUBROUTINE gle_cholesky_stab

! *****************************************************************************
!> \brief ...
!> \param thermostat_info ...
!> \param simpar ...
!> \param local_molecules ...
!> \param molecule_set ...
!> \param molecule_kind_set ...
!> \param gle ...
!> \param para_env ...
!> \param gci ...
!> \author
! *****************************************************************************
  SUBROUTINE gle_to_particle_mapping (thermostat_info, simpar, local_molecules,&
         molecule_set, molecule_kind_set, gle, para_env, gci)

    TYPE(thermostat_info_type), POINTER      :: thermostat_info
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(gle_type), POINTER                  :: gle
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_constraint_type), POINTER    :: gci

    CHARACTER(LEN=*), PARAMETER :: routineN = 'gle_to_particle_mapping', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, imap, j, mal_size, &
                                                natoms_local, nkind, number, &
                                                region, sum_of_thermostats
    INTEGER, DIMENSION(:), POINTER           :: deg_of_freedom, &
                                                massive_atom_list
    LOGICAL                                  :: do_shell, failure
    REAL(KIND=dp)                            :: fac
    TYPE(map_info_type), POINTER             :: map_info

    failure = .FALSE.
    do_shell = .FALSE.
    NULLIFY (massive_atom_list, deg_of_freedom)
    SELECT CASE ( simpar % ensemble )
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,"Unknown ensemble!",para_env)
    CASE ( nve_ensemble, isokin_ensemble, npe_f_ensemble, npe_i_ensemble, nph_uniaxial_ensemble,&
           nph_uniaxial_damped_ensemble, reftraj_ensemble, langevin_ensemble )
       CALL stop_program(routineN,moduleN,__LINE__,"Never reach this point!",para_env)
    CASE ( nvt_ensemble, npt_i_ensemble, npt_f_ensemble)

      map_info => gle%map_info
      nkind = SIZE(molecule_kind_set)
      sum_of_thermostats = thermostat_info%sum_of_thermostats
      map_info%dis_type  = thermostat_info%dis_type
      number             = thermostat_info%number_of_thermostats
      region             = gle%region

      CALL thermostat_mapping_region(map_info, deg_of_freedom, massive_atom_list,&
         molecule_kind_set, local_molecules, molecule_set, para_env, natoms_local,&
         simpar, number, region, gci, do_shell, thermostat_info%map_loc_thermo_gen,&
         sum_of_thermostats)

      ! This is the local number of available thermostats
      gle%loc_num_gle  = number
      gle%glob_num_gle = sum_of_thermostats
      mal_size = SIZE(massive_atom_list)
      CPASSERT(mal_size/=0)
      CALL gle_thermo_create(gle,  mal_size)
      gle%mal(1:mal_size) = massive_atom_list(1:mal_size)

       ! Sum up the number of degrees of freedom on each thermostat.
       ! first: initialize the target
       map_info%s_kin = 0.0_dp
       DO i = 1, 3
          DO j = 1, natoms_local
             map_info%p_kin(i,j)%point = map_info%p_kin(i,j)%point + 1
          END DO
       END DO

       ! If thermostats are replicated but molecules distributed, we have to
       ! sum s_kin over all processors
       IF (map_info%dis_type == do_thermo_communication ) CALL mp_sum (map_info%s_kin, para_env%group )

       ! We know the total number of system thermostats.
       IF ( (sum_of_thermostats == 1).AND.(map_info%dis_type /= do_thermo_no_communication)) THEN
          fac = map_info%s_kin(1)-deg_of_freedom(1)-simpar%nfree_rot_transl
          IF ( fac == 0.0_dp ) THEN
             CALL stop_program(routineN,moduleN,__LINE__,&
                               "Zero degrees of freedom. Nothing to thermalize!",para_env)
          END IF
          gle%nvt(1)%nkt = simpar%temp_ext*fac
          gle%nvt(1)%degrees_of_freedom = FLOOR(fac)
       ELSE
          DO i = 1, gle%loc_num_gle
             imap = map_info%map_index(i)
             fac = (map_info%s_kin(imap) - deg_of_freedom(i) )
             gle%nvt(i)%nkt = simpar%temp_ext * fac
             gle%nvt(i)%degrees_of_freedom = FLOOR(fac)
          END DO
       END IF
       DEALLOCATE (deg_of_freedom)
       DEALLOCATE ( massive_atom_list)
    END SELECT

  END SUBROUTINE gle_to_particle_mapping

! *****************************************************************************
!> \brief ...
!> \param gle ...
!> \param gle_section ...
!> \param save_mem ...
!> \param restart ...
! *****************************************************************************
  SUBROUTINE restart_gle(gle, gle_section, save_mem, restart)

    TYPE(gle_type), POINTER                  :: gle
    TYPE(section_vals_type), POINTER         :: gle_section
    LOGICAL, INTENT(IN)                      :: save_mem
    LOGICAL, INTENT(OUT)                     :: restart

    CHARACTER(len=*), PARAMETER :: routineN = 'restart_gle', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=rng_record_length)         :: rng_record
    INTEGER                                  :: glob_num, i, ind, j, loc_num, &
                                                n_rep
    LOGICAL                                  :: explicit
    REAL(KIND=dp), DIMENSION(:), POINTER     :: buffer
    TYPE(map_info_type), POINTER             :: map_info
    TYPE(section_vals_type), POINTER         :: work_section

    NULLIFY(buffer, work_section)


    explicit = .FALSE.
    restart = .FALSE.

    IF (ASSOCIATED(gle_section)) THEN
       work_section => section_vals_get_subs_vals(gle_section,"S")
       CALL section_vals_get(work_section, explicit=explicit)
       restart=explicit
    END IF

    IF (restart) THEN
       map_info => gle%map_info
       CALL section_vals_val_get(gle_section,"S%_DEFAULT_KEYWORD_",r_vals=buffer)
       DO i = 1, SIZE ( gle % nvt, 1)
          ind = map_info%index(i)
          ind = (ind-1) * (gle% ndim)
          DO j = 1, SIZE(gle%nvt(i)%s, 1)
             ind = ind + 1
             gle % nvt(i)%s(j)  = buffer ( ind )
          END DO
       END DO

       IF(save_mem) THEN
          NULLIFY(work_section)
          work_section => section_vals_get_subs_vals(gle_section,"S")
          CALL section_vals_remove_values(work_section)
       END IF

      ! Possibly restart the initial thermostat energy
       work_section => section_vals_get_subs_vals(section_vals=gle_section,&
           subsection_name="THERMOSTAT_ENERGY")
       CALL section_vals_get(work_section,explicit=explicit)
       IF (explicit) THEN
         CALL section_vals_val_get(section_vals=work_section,keyword_name="_DEFAULT_KEYWORD_",&
              n_rep_val=n_rep)
         IF (n_rep==gle%glob_num_gle) THEN
            DO i = 1, gle%loc_num_gle
               ind = map_info%index(i)
               CALL section_vals_val_get(section_vals=work_section,keyword_name="_DEFAULT_KEYWORD_",&
                    i_rep_val=ind,r_val=gle%nvt(i)%thermostat_energy)
            END DO
         ELSE
            CALL stop_program(routineN,moduleN,__LINE__,&
                              "Number of thermostat energies not equal to the number of "//&
                              "total thermostats!")
         END IF
       END IF

      ! Possibly restart the random number generators for the different thermostats
       work_section => section_vals_get_subs_vals(section_vals=gle_section,&
           subsection_name="RNG_INIT")

       CALL section_vals_get(work_section,explicit=explicit)
       IF (explicit) THEN
         CALL section_vals_val_get(section_vals=work_section,keyword_name="_DEFAULT_KEYWORD_",&
              n_rep_val=n_rep)

         glob_num =gle%glob_num_gle
         loc_num =gle%loc_num_gle
         IF (n_rep==glob_num) THEN
           DO i = 1, loc_num
              ind =  map_info%index(i)
              CALL section_vals_val_get(section_vals=work_section,keyword_name="_DEFAULT_KEYWORD_",&
                   i_rep_val=ind,c_val=rng_record)
               CALL read_rng_stream(rng_stream=gle%nvt(i)%gaussian_rng_stream,&
                    rng_record=rng_record)
           END DO
         ELSE
           CALL stop_program(routineN,moduleN,__LINE__,&
                             "Number pf restartable stream not equal to the number of "//&
                             "total thermostats!")
         END IF
       END IF
    END IF

  END SUBROUTINE restart_gle

! *****************************************************************************
!> \brief ...
!> \param gle ...
! *****************************************************************************
  SUBROUTINE init_gle_variables (gle)

    TYPE(gle_type), POINTER                  :: gle

    CHARACTER(len=*), PARAMETER :: routineN = 'init_gle_variables', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j
    LOGICAL                                  :: failure
    REAL(dp)                                 :: rr(gle%ndim), &
                                                cc(gle%ndim,gle%ndim)
    TYPE(rng_stream_type), POINTER           :: rng_stream

    failure = .FALSE.

    CALL gle_cholesky_stab(gle%c_mat, cc, gle%ndim)
    DO i = 1, gle%loc_num_gle
      rng_stream => gle%nvt(i)%gaussian_rng_stream
      DO j = 1, gle%ndim
        ! here s should be properly initialized, when it is not read from restart
        rr(j) =  next_random_number( rng_stream)
      END DO
      gle%nvt(i)%s = MATMUL(cc,rr)
    END DO

  END SUBROUTINE init_gle_variables

END MODULE gle_system_dynamics
