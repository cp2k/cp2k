!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_scf [1.0] *
!!
!!   NAME
!!     qs_scf
!!
!!   FUNCTION
!!     Routines for the Quickstep SCF run.
!!
!!   AUTHOR
!!     Matthias Krack (30.04.2001)
!!
!!   MODIFICATION HISTORY
!!     - Joost VandeVondele (02.2002)
!!          added code for: incremental (pab and gvg) update
!!                           initialisation (init_cube, l_info)
!!     - Joost VandeVondele (02.2002)
!!          called the poisson code of the classical part
!!          this takes into account the spherical cutoff and allows for
!!          isolated systems
!!     - Joost VandeVondele (02.2002)
!!          added multiple grid feature
!!          changed to spherical cutoff consistently (?)
!!          therefore removed the gradient correct functionals
!!     - updated with the new QS data structures (10.04.02,MK)
!!     - copy_matrix replaced by transfer_matrix (11.04.02,MK)
!!     - nrebuild_rho and nrebuild_gvg unified (12.04.02,MK)
!!     - set_mo_occupation for smearing of the MO occupation numbers
!!       (17.04.02,MK)
!!     - MO level shifting added (22.04.02,MK)
!!     - Usage of TYPE mo_set_p_type
!!     - Joost VandeVondele (05.2002)
!!           added cholesky based diagonalisation
!!
!!   SOURCE
!******************************************************************************

MODULE qs_scf

! *****************************************************************************

  USE kinds, ONLY: int_size,&
                   wp => dp

  USE atomic_kind_types,        ONLY: atomic_kind_type,&
                                      get_atomic_kind,&
                                      get_atomic_kind_set
  USE basis_set_types,          ONLY: get_gto_basis_set,&
                                      gto_basis_set_type
  USE coefficient_types,        ONLY: coeff_allocate,&
                                      coeff_copy,&
                                      coeff_transform_space,&
                                      coeff_type,&
                                      coeff_zero
  USE cube_utils,               ONLY: cube_info_type,&
                                      destroy_cube_info,&
                                      init_cube_info
  USE dft_types,                ONLY: dft_control_type
  USE external_potential_types, ONLY: all_potential_type,&
                                      get_potential,&
                                      gth_potential_type
  USE functionals
  USE gaussian_gridlevels,      ONLY: gridlevel_info_type,&
                                      init_gaussian_gridlevel,& 
                                      destroy_gaussian_gridlevel
  USE global_types,             ONLY: DEBUG,&
                                      global_environment_type
  USE greens_fn,                ONLY: pw_green_fn_init
  USE hartree,                  ONLY: calculate_hartree
  USE l_utils,                  ONLY: destroy_l_info,&
                                      init_l_info,&
                                      l_info_type
  USE mathconstants,            ONLY: fourpi,twopi
  USE memory_utilities,         ONLY: reallocate
  USE message_passing,          ONLY: mp_max,mp_sum,mp_sync,&
                                      mp_range,mp_allgather
  USE particle_types,           ONLY: particle_type
  USE pw_grid_types,            ONLY: HALFSPACE,pw_grid_type
  USE pw_grids,                 ONLY: pw_find_cutoff,&
                                      pw_grid_construct,&
                                      pw_grid_setup
  USE pw_types,                 ONLY: COMPLEXDATA1D,COMPLEXDATA3D,&
                                      REALDATA3D,REALSPACE,&
                                      RECIPROCALSPACE,&
                                      pw_fft_wrap,pw_type
  USE pws,                      ONLY: init_pw_poisson_solver,&
                                      pw_poisson_solver
  USE qs_blacs,                 ONLY: allocate_blacs_matrix,&
                                      blacs_add_to_element,&
                                      blacs_gemm,&
                                      blacs_get_element,&
                                      blacs_matrix_type,&
                                      blacs_set_all,&
                                      blacs_set_element,&
                                      blacs_syevx,&
                                      blacs_symm,&
                                      copy_blacs_to_blacs_matrix,&
                                      copy_blacs_to_sparse_matrix,&
                                      copy_sparse_to_blacs_matrix,&
                                      finish_blacs,&
                                      get_blacs_matrix_info,&
                                      power_blacs_matrix,&
                                      replicate_blacs_matrix,&
                                      start_blacs,&
                                      symmetrise_blacs_matrix, &
                                      blacs_cholesky_decompose, &
                                      blacs_cholesky_reduce, &
                                      blacs_cholesky_restore
  USE qs_collocate_density,     ONLY: calculate_rho_core,&
                                      calculate_rho_elec,&
                                      calculate_total_rho
  USE qs_core_energies,         ONLY: calculate_ecore,&
                                      calculate_ecore_overlap,&
                                      calculate_ecore_self
  USE qs_diis,                  ONLY: eps_diis,max_diis,scf_diis
  USE qs_environment_types,     ONLY: get_qs_env,&
                                      qs_environment_type,&
                                      set_qs_env
  USE qs_integrate_potential,   ONLY: integrate_v_rspace
  USE qs_mo_types,              ONLY: allocate_mo_set,&
                                      calculate_density_matrix,&
                                      correct_mo_eigenvalues,&
                                      get_mo_set,&
                                      mo_set_p_type,&
                                      mo_set_type,&
                                      read_mo_set,&
                                      set_mo_occupation,&
                                      write_mo_set
  USE qs_overlap,               ONLY: write_blacs_matrix,&
                                      write_sparse_matrix
  USE qs_parser,                ONLY: finish_parser,&
                                      read_object,&
                                      start_parser,&
                                      test_object
  USE simulation_cell,          ONLY: cell_type,&
                                      get_cell
  USE sparse_matrix_types,      ONLY: add_matrices,&
                                      allocate_matrix,&
                                      deallocate_matrix,&
                                      first_block_node,&
                                      get_block_node,&
                                      get_matrix_info,&
                                      next_block_node,&
                                      real_block_node_type,&
                                      real_matrix_p_type,&
                                      replicate_matrix_structure,&
                                      set_matrix,&
                                      symmetrise_diagonal_blocks,&
                                      transfer_matrix
  USE termination,              ONLY: stop_memory,&
                                      stop_program
  USE timings,                  ONLY: timeset,&
                                      timestop
  USE timesl,                   ONLY: cputime

  IMPLICIT NONE

  PRIVATE

  TYPE(blacs_matrix_type), POINTER  :: ortho,scf_work1,scf_work2
  TYPE(cell_type), POINTER          :: cell
  TYPE(dft_control_type), POINTER   :: dft_control
  TYPE(real_matrix_p_type), POINTER :: h,s

  TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
  TYPE(coeff_type), DIMENSION(:), POINTER       :: rho_gspace,rho_rspace,&
                                                   rho_rspace_old,v_gspace,&
                                                   v_rspace
  TYPE(mo_set_p_type), DIMENSION(:), POINTER    :: c
  TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set
  TYPE(pw_grid_type), DIMENSION(:), POINTER     :: pw_grid

  TYPE(coeff_type)          :: drho_rspace,v_rspace_old,vxcg_rspace,&
                               work_gspace,work_rspace
  TYPE(gridlevel_info_type) :: gridlevel_info
  TYPE(l_info_type)         :: l_info
  TYPE(real_matrix_p_type)  :: ks,p,p_old,v

  CHARACTER(LEN=60) :: c_fun,x_fun,xc_fun
  CHARACTER(LEN=10) :: density_guess
  REAL(wp)          :: rel_cutoff,ec,ecore,ecore_overlap,ecore_self,ehartree,&
                       eps_eigval,eps_scf,etotal,ex,level_shift,p_mix,smear,&
                       total_rho_core_rspace,total_rho_elec_rspace,&
                       total_rho_gspace,total_rho_rspace,work_syevx
  INTEGER           :: charge,max_scf,nelectron,nrebuild,nrow_block,&
                       ncol_block,nprow,npcol,maxl,ngrid_level,nkind
  LOGICAL           :: gradient_functional,gth_potential_present,orthorhombic,&
                       rebuild,use_cholesky

  TYPE(coeff_type), DIMENSION(3) :: nabla_rho_rspace

  TYPE(cube_info_type), DIMENSION(:), POINTER :: cube_info
  REAL(wp), DIMENSION(:), POINTER             :: cutoff

! *** Public variables ***

  PUBLIC :: nelectron

! *** Public subroutines ***

  PUBLIC :: read_scf_parameters,&
            scf,&
            write_scf_parameters

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE scf(qs_env,globenv)

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(qs_environment_type), INTENT(INOUT)  :: qs_env

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE scf (MODULE qs_scf)"

!   *** Local parameters ***

    INTEGER, PARAMETER :: maxgridpoints = 100

!   *** Local variables ***

    TYPE(mo_set_type), POINTER :: mo_set

    REAL(wp) :: delta,diis_error,maxradius,t1,t2
    INTEGER  :: cmax,cmaxl,context,handle,homo,igrid_level,imo,iscf,istat,nmo,&
                output_unit
    LOGICAL  :: diis_step,do_level_shift,ionode

    TYPE(blacs_matrix_type), POINTER :: mo_eigenvectors

!   ---------------------------------------------------------------------------



!   *** Quick return, if no SCF iteration is requested ***

    IF (max_scf < 1) RETURN

    CALL timeset("scf","I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    h=h,&
                    s=s,&
                    particle_set=particle_set)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             gth_potential_present=gth_potential_present,&
                             maxl=maxl,&
                             nelectron=nelectron)

    nelectron = nelectron - charge

    CALL get_cell(cell=cell,orthorhombic=orthorhombic)

    gradient_functional = dft_control%gradient_functional
    rel_cutoff = dft_control%qs_control%relative_cutoff
    cutoff => dft_control%qs_control%e_cutoff
    ngrid_level = SIZE(cutoff)

    x_fun = dft_control%exchange_functional
    c_fun = dft_control%correlation_functional
    xc_fun = dft_control%functional

    nkind = SIZE(atomic_kind_set)


    CALL init_functionals()

    CALL init_gaussian_gridlevel(gridlevel_info,ngrid_level,cutoff,rel_cutoff)

    CALL init_grid(globenv)

    ALLOCATE (cube_info(ngrid_level),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"cube_info",0)

    cmax = 0
    DO igrid_level=1,ngrid_level
      maxradius = maxgridpoints*MAXVAL(rho_rspace(igrid_level)%pw%pw_grid%dr)
      CALL init_cube_info(cube_info(igrid_level),&
                          rho_rspace(igrid_level)%pw%pw_grid%dr(:),&
                          maxradius,cmaxl)
      cmax = MAX(cmax,cmaxl)
    END DO

    ! maxl+1 to allow for the forces
    CALL init_l_info(l_info,maxl+1,cmax)

    IF (ionode.AND.globenv%print%scf) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
        "SCF WAVEFUNCTION OPTIMIZATION"
    END IF

    CALL init_scf_run(mo_set,context,qs_env,globenv)

    etotal = 0.0_wp
    iscf = 0
    diis_step = .FALSE.
    rebuild = .TRUE.

    IF (ionode.AND.globenv%print%scf) THEN
      WRITE (UNIT=output_unit,&
             FMT="(/,T3,A,T9,A,T34,A,T49,A,T68,A,/,T3,A)")&
        "Step","Update method","Time","Convergence","Total energy",&
        REPEAT("-",77)
    END IF

!   *** SCF loop ***

    scf_loop: DO

      iscf = iscf + 1

      t1 = cputime()

      IF (globenv%print%density_matrix.AND.globenv%print%each_scf_step) THEN
        CALL write_sparse_matrix(p%matrix,4,6,qs_env,globenv)
      END IF

      CALL build_kohn_sham_matrix(qs_env,globenv)

      rebuild = (MODULO(iscf,nrebuild) == 0)

      IF (globenv%print%kohn_sham_matrix.AND.globenv%print%each_scf_step) THEN
        CALL write_sparse_matrix(ks%matrix,4,6,qs_env,globenv)
      END IF

      CALL copy_sparse_to_blacs_matrix(ks%matrix,scf_work1,context,globenv)

      IF (iscf > 1) THEN
        CALL scf_diis(mo_set,scf_work1,scf_work2,delta,diis_error,diis_step,&
                      context,qs_env,globenv)
      END IF

!     *** Perform MO level shifting, if requested ***

      do_level_shift = ((level_shift /= 0.0_wp).AND.&
                        ((density_guess == "CORE").OR.(iscf > 1)))

      IF (do_level_shift) THEN

        IF (use_cholesky) THEN
           stop 'do_level_shift .and. use_cholesky not checked'
        ENDIF

        CALL get_mo_set(mo_set=mo_set,&
                        homo=homo,&
                        nmo=nmo,&
                        eigenvectors=mo_eigenvectors)

!       *** Orthogonalisation matrix: MO eigenvectors ***

        CALL orthogonalise_matrix(ortho,scf_work1,scf_work2,context,globenv)

!       *** Shift the Kohn-Sham matrix diagonal ***
!       *** elements of the unoccupied MOs      ***

        DO imo=homo+1,nmo
          CALL blacs_add_to_element(scf_work1,imo,imo,level_shift,context,&
                                    globenv)
        END DO

      ELSE

!       *** Orthogonalisation matrix: S**(-1/2) ***
        CALL orthogonalise_matrix(ortho,scf_work1,scf_work2,context,globenv)

      END IF

      CALL eigensolver(scf_work1,mo_set,ortho,scf_work2,context,globenv)

      IF (do_level_shift) THEN

!       *** Unshift the MO eigenvalues, if level-shifting is used ***

        CALL correct_mo_eigenvalues(mo_set,level_shift)

!       *** Store new MOs as orthogonalisation matrix ***

        CALL copy_blacs_to_blacs_matrix(mo_eigenvectors,ortho)

      END IF

!     *** Update the MO occupation numbers, if smearing is used ***

      IF (smear /= 0.0_wp) CALL set_mo_occupation(mo_set,smear)

      IF (globenv%print%each_scf_step) THEN
        CALL write_mo_set(mo_set,atomic_kind_set,particle_set,4,6,context,&
                          globenv)
      END IF

      CALL calculate_density_matrix(mo_set,scf_work1,context,globenv)

      CALL copy_blacs_to_sparse_matrix(scf_work1,ks%matrix,context,globenv)
      CALL symmetrise_diagonal_blocks(ks%matrix)

      t2 = cputime()

      IF (diis_step) THEN
        CALL density_mixing(ks,p,1.0_wp,delta,globenv)
        IF (ionode.AND.globenv%print%scf) THEN
          WRITE (UNIT=output_unit,&
                 FMT="(T2,I5,2X,A,T15,E10.2,T30,F8.2,T40,2F20.10)")&
            iscf,"DIIS",diis_error,t2 - t1,delta,etotal
        END IF
      ELSE
        CALL density_mixing(ks,p,p_mix,delta,globenv)
        IF (ionode.AND.globenv%print%scf) THEN
          WRITE (UNIT=output_unit,&
                 FMT="(T2,I5,2X,A,T15,F6.2,T30,F8.2,T40,2F20.10)")&
            iscf,"Mixing",p_mix,t2 - t1,delta,etotal
        END IF
      END IF

      IF (delta < eps_scf) THEN
        IF (ionode.AND.globenv%print%scf) THEN
          WRITE(UNIT=output_unit,FMT="(/,T3,A,/)")&
            "*** SCF run converged ***"
        END IF
        EXIT scf_loop
      ELSE IF (iscf == max_scf) THEN
        IF (ionode.AND.globenv%print%scf) THEN
          WRITE(UNIT=output_unit,FMT="(/,T3,A,/)")&
            "*** SCF run NOT converged ***"
        END IF
        EXIT scf_loop
      END IF

    END DO scf_loop

    IF (ionode.AND.globenv%print%scf) THEN
      WRITE (UNIT=output_unit,FMT="(/,(T3,A,T40,2F20.10))")&
        "Total electronic density (r-space): ",&
        total_rho_elec_rspace,total_rho_elec_rspace + REAL(nelectron,wp),&
        "Total core charge density (r-space):",&
        total_rho_core_rspace,total_rho_core_rspace - REAL(nelectron,wp)
      WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
        "Total charge density (r-space):     ",total_rho_rspace,&
        "Total charge density (g-space):     ",total_rho_gspace
      WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
        "Overlap energy of the core charge distribution:",ecore_overlap,&
        "Self energy of the core charge distribution:   ",ecore_self,&
        "Core Hamiltonian energy:                       ",ecore,&
        "Hartree energy:                                ",ehartree,&
        "Exchange-correlation energy:                   ",ex + ec,&
        "Total energy:                                  ",etotal
    END IF

!   *** Undo mixing of the density matrix (restore original density) ***

    IF (.NOT.diis_step) CALL transfer_matrix(ks%matrix,p%matrix)

    CALL write_mo_set(mo_set,atomic_kind_set,particle_set,4,6,context,globenv)

    IF (globenv%print%density_matrix) THEN
      CALL write_sparse_matrix(p%matrix,4,6,qs_env,globenv)
    END IF

    IF (globenv%print%kohn_sham_matrix) THEN
      CALL build_kohn_sham_matrix(qs_env,globenv)
      CALL write_sparse_matrix(ks%matrix,4,6,qs_env,globenv)
    END IF

!   *** Write restart file ***

    CALL write_mo_set(mo_set,context,globenv)

    DO igrid_level=1,ngrid_level
      CALL destroy_cube_info(cube_info(igrid_level))
    END DO

    DEALLOCATE (cube_info,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"cube_info")

    CALL destroy_l_info(l_info)
    CALL destroy_gaussian_gridlevel(gridlevel_info,globenv)

    CALL finish_blacs(context,globenv)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE scf

! *****************************************************************************

  SUBROUTINE build_kohn_sham_matrix(qs_env,globenv)

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(qs_environment_type), INTENT(INOUT)  :: qs_env

!   *** Local variables ***

    INTEGER :: handle,output_unit
    LOGICAL :: ionode

!   ---------------------------------------------------------------------------

    CALL timeset("build_kohn_sham_matrix","I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr

    CALL calculate_ecore(h,p,ecore,globenv)

    CALL calculate_rho_elec(p,p_old,rebuild,rho_rspace,rho_rspace_old,&
                            rho_gspace,gridlevel_info,cube_info,l_info,&
                            total_rho_elec_rspace,qs_env,globenv)

    CALL coeff_zero(v_rspace(1))

    IF (gradient_functional) THEN
! we have to wait for the general xc functions since we changed
! to a complex1D data structure for the density
      CALL stop_program("build_kohn","feature has been removed")
!MK      CALL calculate_nabla_rho_elec(globenv)
!MK      CALL coeff_zero(vxcg_rspace)
    END IF

    CALL calculate_xc_potential(globenv)

    CALL calculate_rho_core(rho_rspace(1),cube_info(1),l_info,&
                            total_rho_rspace,qs_env,globenv)

    total_rho_core_rspace = total_rho_rspace - total_rho_elec_rspace

    CALL pw_fft_wrap(rho_rspace(1)%pw,rho_gspace(1)%pw,&
                     debug=(globenv%print%level == DEBUG))

    total_rho_gspace = calculate_total_rho(rho_gspace(1))

    IF (ionode.AND.globenv%print%total_densities) THEN
      WRITE (UNIT=output_unit,FMT="(/,(T3,A,T40,2F20.10))")&
        "Total electronic density (r-space): ",&
        total_rho_elec_rspace,total_rho_elec_rspace + REAL(nelectron,wp),&
        "Total core charge density (r-space):",&
        total_rho_core_rspace,total_rho_core_rspace - REAL(nelectron,wp)
      WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
        "Total charge density (r-space):     ",total_rho_rspace,&
        "Total charge density (g-space):     ",total_rho_gspace
    END IF

    CALL calculate_hartree(rho_gspace(1),ehartree,work_gspace)

    etotal = ecore_overlap + ecore_self + ecore + ehartree + ex + ec

    IF (ionode.AND.globenv%print%scf_energies) THEN
      WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
        "Core Hamiltonian energy:                       ",ecore,&
        "Hartree energy:                                ",ehartree,&
        "Exchange-correlation energy:                   ",ex + ec
    END IF

    CALL pw_fft_wrap(work_gspace%pw,work_rspace%pw,&
                     debug=(globenv%print%level == DEBUG))

    v_rspace(1)%pw%cr3d(:,:,:) = v_rspace(1)%pw%pw_grid%dvol*&
                                 (v_rspace(1)%pw%cr3d(:,:,:) +&
                                  work_rspace%pw%cr3d(:,:,:))

    IF (nrebuild == 1) THEN
      CALL transfer_matrix(h%matrix,ks%matrix)
      CALL integrate_v_rspace(v_rspace,v_gspace,gridlevel_info,v_rspace_old,&
                              cube_info,l_info,rebuild,p,ks,qs_env,globenv)
    ELSE
      IF (rebuild) CALL set_matrix(v%matrix,0.0_wp)
      CALL integrate_v_rspace(v_rspace,v_gspace,gridlevel_info,v_rspace_old,&
                              cube_info,l_info,rebuild,p,v,qs_env,globenv)
      CALL add_matrices(ks%matrix,1.0_wp,h%matrix,1.0_wp,v%matrix)
    END IF

!   *** Keep a copy of the old density matrix ***

    CALL transfer_matrix(p%matrix,p_old%matrix)

    CALL timestop(0.0_wp,handle)

 END SUBROUTINE build_kohn_sham_matrix

! *****************************************************************************

  SUBROUTINE init_grid(globenv)

    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE init_grid (MODULE qs_scf)"

!   *** Local variables ***

    INTEGER :: handle,i,igrid_level,istat,output_unit
    LOGICAL :: ionode

!   ---------------------------------------------------------------------------

    CALL timeset("init_grid","I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr

!   *** Allocate the multi-grid arrays ***

    ALLOCATE (pw_grid(ngrid_level),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"pw_grid",0)
    ALLOCATE (rho_rspace(ngrid_level),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"rho_rspace",0)
    ALLOCATE (rho_gspace(ngrid_level),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"rho_gspace",0)
    ALLOCATE (v_rspace(ngrid_level),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"v_rspace",0)
    ALLOCATE (v_gspace(ngrid_level),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"v_gspace",0)

    ALLOCATE (rho_rspace_old(ngrid_level),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"rho_rspace_old",0)

!   *** Construct the various coeff_type grids ***

    DO igrid_level=1,ngrid_level
      CALL pw_grid_construct(pw_grid(igrid_level))
      pw_grid(igrid_level)%grid_span = HALFSPACE
      IF (igrid_level == 1) THEN
        IF (ionode.AND.globenv%print%pw_grid_information) THEN
          CALL pw_grid_setup(cell,pw_grid(1),cutoff(1),&
                             info=output_unit,&
                             orthorhombic=orthorhombic)
         ELSE
           CALL pw_grid_setup(cell,pw_grid(1),cutoff(1),&
                              orthorhombic=orthorhombic)
         END IF
      ELSE
         IF (ionode.AND.globenv%print%pw_grid_information) THEN
           CALL pw_grid_setup(cell,pw_grid(igrid_level),cutoff(igrid_level),&
                              info=output_unit,&
                              orthorhombic=orthorhombic,&
                              ref_grid=pw_grid(1))
         ELSE
           CALL pw_grid_setup(cell,pw_grid(igrid_level),cutoff(igrid_level),&
                              orthorhombic=orthorhombic,&
                              ref_grid=pw_grid(1))
           END IF
       END IF

       CALL coeff_allocate(rho_rspace(igrid_level),pw_grid(igrid_level),&
                           REALDATA3D)
       rho_rspace(igrid_level)%pw%in_space = REALSPACE
       CALL coeff_allocate(rho_gspace(igrid_level),pw_grid(igrid_level),&
                           COMPLEXDATA1D)
       rho_gspace(igrid_level)%pw%in_space = RECIPROCALSPACE
       CALL coeff_allocate(rho_rspace_old(igrid_level),pw_grid(igrid_level),&
                           REALDATA3D)
       rho_rspace_old%pw%in_space = REALSPACE
       CALL coeff_allocate(v_rspace(igrid_level),pw_grid(igrid_level),&
                           REALDATA3D)
       v_rspace(igrid_level)%pw%in_space = REALSPACE
       CALL coeff_allocate(v_gspace(igrid_level),pw_grid(igrid_level),&
                           COMPLEXDATA1D)
       v_gspace(igrid_level)%pw%in_space = RECIPROCALSPACE

    END DO

!   *** Following quantities are only needed for the highest cutoff grid ***

    CALL pw_green_fn_init(cell%green,pw_grid(1))

    CALL init_pw_poisson_solver(cell%green)

    CALL coeff_allocate(work_rspace,pw_grid(1),REALDATA3D)
    work_rspace%pw%in_space = REALSPACE
    CALL coeff_allocate(work_gspace,pw_grid(1),COMPLEXDATA1D)
    work_gspace%pw%in_space = RECIPROCALSPACE
    CALL coeff_allocate(v_rspace_old,pw_grid(1),REALDATA3D)
    v_rspace_old%pw%in_space = REALSPACE

    IF (gradient_functional) THEN
      DO i=1,3
        CALL coeff_allocate(nabla_rho_rspace(i),pw_grid(1),REALDATA3D)
        nabla_rho_rspace(i)%pw%in_space = REALSPACE
      END DO
      CALL coeff_allocate(drho_rspace,pw_grid(1),REALDATA3D)
      drho_rspace%pw%in_space = REALSPACE
      CALL coeff_allocate(vxcg_rspace,pw_grid(1),REALDATA3D)
      vxcg_rspace%pw%in_space = REALSPACE
    END IF

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE init_grid

! *****************************************************************************

  SUBROUTINE init_scf_run(mo_set,context,qs_env,globenv)

!   Purpose: Initialise a SCF run.

!   History: - Creation (30.04.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), POINTER                :: mo_set
    TYPE(qs_environment_type), INTENT(INOUT)  :: qs_env
    INTEGER, INTENT(OUT)                      :: context

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE init_scf_run (MODULE qs_scf)"

!   *** Local variables ***

    TYPE(blacs_matrix_type), POINTER :: mo_eigenvectors

    INTEGER :: handle,homo,ikind,istat,nao,ndep,nmo,output_unit
    LOGICAL :: ionode

    INTEGER, DIMENSION(0:globenv%num_pe-1) :: nblock_pe,nelement_pe

!   ---------------------------------------------------------------------------

    CALL timeset("init_scf_run","I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr

!   *** Initialise BLACS ***

    CALL start_blacs(nprow,npcol,context,globenv)

!   *** Get the dimension of the full SCF matrices, ***
!   *** i.e. the total number of atomic orbitals    ***

    CALL get_matrix_info(matrix=h%matrix,nrow=nao)

    IF (globenv%print%mo_eigenvalues.OR.&
        (level_shift /= 0.0_wp).OR.&
        (smear /= 0.0_wp)) THEN
      nmo = nao
    ELSE
      nmo = nelectron/2
    END IF

!   *** Allocate a MO set ***

    ALLOCATE (c(1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"c",int_size)
    NULLIFY (c(1)%mo_set)

!   *** Allocate the distributed MO eigenvectors ***

    CALL allocate_mo_set(c(1)%mo_set,nao,nmo,nelectron,2.0_wp,nrow_block,&
                         ncol_block,context,globenv)

!   *** Put the MO set in the QS environment ***

    CALL set_qs_env(qs_env=qs_env,c=c)

!   *** Get the HOMO and the MO eigenvectors ***

    mo_set => c(1)%mo_set

    CALL get_mo_set(mo_set=mo_set,&
                    homo=homo,&
                    eigenvectors=mo_eigenvectors)

!   *** Get BLACS block size of the MO eigenvector matrix      ***
!   *** which has to fit to the other distributed SCF matrices ***

    CALL get_blacs_matrix_info(matrix=mo_eigenvectors,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block)

!   *** Allocate the distributed SCF matrices ***

    CALL allocate_blacs_matrix(new_matrix=ortho,&
                               nrow_global=nao,&
                               ncol_global=nao,&
                               nrow_block=nrow_block,&
                               ncol_block=ncol_block,&
                               name="ORTHOGONALIZATION MATRIX",&
                               context=context,&
                               globenv=globenv)

    CALL replicate_blacs_matrix(prototype_matrix=ortho,&
                                new_matrix=scf_work1,&
                                name="SCF WORK MATRIX 1")

    CALL replicate_blacs_matrix(prototype_matrix=ortho,&
                                new_matrix=scf_work2,&
                                name="SCF WORK MATRIX 2")

    CALL copy_sparse_to_blacs_matrix(h%matrix,scf_work1,context,globenv)

!   *** Redistribute the core Hamiltonian matrix in ***
!   *** order to eliminate redundant atomic blocks  ***

    IF ((globenv%num_pe > 1).AND.gth_potential_present) THEN
      CALL deallocate_matrix(h%matrix)
      CALL replicate_matrix_structure(s%matrix,h%matrix,&
                                      "CORE HAMILTONIAN MATRIX")
      CALL copy_blacs_to_sparse_matrix(scf_work1,h%matrix,context,globenv)
    END IF

    NULLIFY (ks%matrix)
    CALL replicate_matrix_structure(h%matrix,ks%matrix,"KOHN-SHAM MATRIX")

    NULLIFY (p%matrix)
    CALL replicate_matrix_structure(h%matrix,p%matrix,"DENSITY MATRIX")

    NULLIFY (p_old%matrix)
    CALL replicate_matrix_structure(h%matrix,p_old%matrix,"OLD DENSITY MATRIX")

    NULLIFY (v%matrix)
    CALL replicate_matrix_structure(h%matrix,v%matrix,"POTENTIAL MATRIX")

    CALL calculate_ecore_self(ecore_self,qs_env)
    CALL calculate_ecore_overlap(ecore_overlap,qs_env,globenv)

    IF (ionode.AND.globenv%print%scf_energies) THEN
      WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
        "Overlap energy of the core charge distribution:",ecore_overlap,&
        "Self energy of the core charge distribution:   ",ecore_self
    END IF

    if (use_cholesky) then
       ndep=0
       CALL copy_sparse_to_blacs_matrix(s%matrix,ortho,context,globenv)
       CALL blacs_cholesky_decompose(ortho,context,globenv)
    else
       CALL calculate_ortho_matrix(ortho,scf_work2,ndep,context,globenv)
    endif

    IF (ionode.AND.globenv%print%scf) THEN
      WRITE (UNIT=output_unit,FMT="(/,(T3,A,I10))")&
        "Number of electrons:                    ",nelectron,&
        "Number of occupied orbitals:            ",homo,&
        "Number of orbital functions:            ",nao,&
        "Number of independent orbital functions:",nao - ndep
    END IF

    IF (globenv%print%ortho_matrix) THEN
      CALL write_blacs_matrix(ortho,4,6,context,qs_env,globenv)
    END IF

    CALL calculate_first_density_matrix(ortho,mo_set,p,scf_work1,scf_work2,&
                                        context,globenv)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE init_scf_run

! *****************************************************************************

  SUBROUTINE calculate_ortho_matrix(ortho,work,ndep,context,globenv)

!   Purpose: Calculate the orthogonalization matrix (S**(-1/2))

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(blacs_matrix_type), POINTER          :: ortho,work
    INTEGER, INTENT(IN)                       :: context
    INTEGER, INTENT(OUT)                      :: ndep

!   *** Local variables ***

    INTEGER :: handle,nao

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_ortho_matrix","I","",handle)

    CALL copy_sparse_to_blacs_matrix(s%matrix,ortho,context,globenv)
    CALL power_blacs_matrix(ortho,work,-0.5_wp,eps_eigval,ndep,work_syevx,&
                            context,globenv)
    CALL symmetrise_blacs_matrix(ortho,work,context,globenv)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_ortho_matrix

! *****************************************************************************

  SUBROUTINE orthogonalise_matrix(ortho,ks,work,context,globenv)

!   Purpose: Orthogonalisation matrix (S**(-1/2))

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: ks,ortho,work
    TYPE(global_environment_type), INTENT(IN) :: globenv
    INTEGER, INTENT(IN)                       :: context

!   *** Local variables ***

    INTEGER :: handle,nao

!   ---------------------------------------------------------------------------

    CALL timeset("orthogonalise_matrix","I","",handle)

    if (use_cholesky) then

       CALL blacs_cholesky_reduce(ks,ortho,context,globenv)

    else

       CALL get_blacs_matrix_info(matrix=ks,nrow_global=nao)
       CALL blacs_symm("L","U",nao,nao,1.0_wp,ks,ortho,0.0_wp,work,context,&
                    globenv)
       CALL blacs_gemm("T","N",nao,nao,nao,1.0_wp,ortho,work,0.0_wp,ks,context,&
                    globenv)
    endif

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE orthogonalise_matrix

! *****************************************************************************

  SUBROUTINE eigensolver(ks,mo_set,ortho,work,context,globenv)

!   Purpose: Diagonalise the Kohn-Sham matrix to get a new set of MO eigen-
!            vectors and MO eigenvalues.

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), POINTER          :: ks,ortho,work
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), POINTER                :: mo_set
    INTEGER, INTENT(IN)                       :: context

!   *** Local variables ***

    INTEGER :: handle,imo,nao,nmo

    TYPE(blacs_matrix_type), POINTER :: mo_eigenvectors
    REAL(wp), DIMENSION(:), POINTER  :: mo_eigenvalues

!   ---------------------------------------------------------------------------

    CALL timeset("eigensolver","I","",handle)

    CALL get_mo_set(mo_set=mo_set,&
                    nao=nao,&
                    nmo=nmo,&
                    eigenvalues=mo_eigenvalues,&
                    eigenvectors=mo_eigenvectors)

!   *** Diagonalise the Kohn-Sham matrix ***

    IF (use_cholesky) THEN

       CALL blacs_syevx(ks,work,mo_eigenvalues,nmo,work_syevx,context,globenv)
       CALL blacs_cholesky_restore(work,nmo,ortho,context,globenv)
       CALL copy_blacs_to_blacs_matrix(work,mo_eigenvectors)

    ELSE

       CALL blacs_syevx(ks,work,mo_eigenvalues,nmo,work_syevx,context,globenv)

       CALL blacs_gemm("N","N",nao,nmo,nao,1.0_wp,ortho,work,0.0_wp,&
                    mo_eigenvectors,context,globenv)
   
    ENDIF

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE eigensolver

! *****************************************************************************

  SUBROUTINE density_mixing(new_density,old_density,p_mix,delta,globenv)

!   Purpose: Perform a density mixing of the old (last SCF iteration) and the
!            new (current) density matrix.

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(real_matrix_p_type), INTENT(INOUT)   :: new_density,old_density
    REAL(wp), INTENT(IN)                      :: p_mix
    REAL(wp), INTENT(OUT)                     :: delta

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: block_node

    REAL(wp) :: r_mix
    INTEGER  :: handle,i,iblock_col,iblock_row,j,nblock_row

    REAL(wp), DIMENSION(:,:), POINTER :: p_new_block,p_old_block

!   ---------------------------------------------------------------------------

    CALL timeset("density_mixing","I","",handle)

    delta = 0.0_wp

    r_mix = 1.0_wp - p_mix

    CALL get_matrix_info(matrix=old_density%matrix,nblock_row=nblock_row)

    DO iblock_row=1,nblock_row

      block_node => first_block_node(matrix=new_density%matrix,&
                                     block_row=iblock_row)

      DO WHILE (ASSOCIATED(block_node))

        CALL get_block_node(block_node=block_node,&
                            block_col=iblock_col,&
                            block=p_new_block)

        CALL get_block_node(matrix=old_density%matrix,&
                            block_row=iblock_row,&
                            block_col=iblock_col,&
                            block=p_old_block)

        DO j=1,SIZE(p_new_block,2)
          DO i=1,SIZE(p_new_block,1)
            p_new_block(i,j) = p_new_block(i,j) - p_old_block(i,j)
            delta = MAX(delta,ABS(p_new_block(i,j)))
            p_old_block(i,j) = p_old_block(i,j) + p_mix*p_new_block(i,j)
            p_new_block(i,j) = p_old_block(i,j) + r_mix*p_new_block(i,j)
          END DO
        END DO

        block_node => next_block_node(block_node)

      END DO

    END DO

    CALL mp_max(delta,globenv%group)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE density_mixing

! *****************************************************************************

  SUBROUTINE calculate_xc_potential(globenv)

    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    INTEGER :: handle,kg

    INTEGER :: local_planes,local_offset

    INTEGER, DIMENSION(3) :: ng

    REAL(wp), DIMENSION(:,:,:), POINTER :: drho,rho,vxc,vxcg

!
!   Joost VdV ( 05-2002 ) modified to use the mp_allgather function
!   each pe computes only part of the grid and this is broadcasted to all 
!   instead of summed. 
!   This scales significantly better (e.g. factor 3 on 12 cpus 32 H2O)
!

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_xc_potential","I","",handle)

    ng(:) = rho_rspace(1)%pw%pw_grid%npts(3)

    rho => rho_rspace(1)%pw%cr3d(:,:,:)
    vxc => v_rspace(1)%pw%cr3d(:,:,:)

    IF (gradient_functional) THEN
      drho => drho_rspace%pw%cr3d(:,:,:)
      vxcg => vxcg_rspace%pw%cr3d(:,:,:)
    END IF

    ex = 0.0_wp
    ec = 0.0_wp

    CALL mp_range(globenv%mepos,globenv%num_pe,ng(3),local_offset,local_planes)

    DO kg=local_offset+1,local_offset+local_planes

!     *** Exchange functionals ***

      SELECT CASE (TRIM(x_fun))
      CASE ("Becke88")
        CALL becke88(rho(:,:,kg),drho(:,:,kg),ex,vxc(:,:,kg),vxcg(:,:,kg))
      CASE ("PBE")
        CALL pbe_x(rho(:,:,kg),drho(:,:,kg),ex,vxc(:,:,kg),vxcg(:,:,kg))
      CASE ("Perdew86")
        CALL perdew86_x(rho(:,:,kg),drho(:,:,kg),ex,vxc(:,:,kg),vxcg(:,:,kg))
      CASE ("Slater")
        CALL slater(rho(:,:,kg),ex,vxc(:,:,kg))
      CASE ("VWN")
        CALL vwn_x(rho(:,:,kg),ex,vxc(:,:,kg))
      END SELECT

!     *** Correlation functionals ***

      SELECT CASE (TRIM(c_fun))
      CASE ("LYP")
        CALL lyp(rho(:,:,kg),drho(:,:,kg),ec,vxc(:,:,kg),vxcg(:,:,kg))
      CASE ("PBE")
        CALL pbe_c(rho(:,:,kg),drho(:,:,kg),ec,vxc(:,:,kg),vxcg(:,:,kg))
      CASE ("Perdew86")
        CALL perdew86_c(rho(:,:,kg),drho(:,:,kg),ec,vxc(:,:,kg),vxcg(:,:,kg))
      CASE ("VWN")
        CALL vwn_c(rho(:,:,kg),ec,vxc(:,:,kg))
      END SELECT

!     *** Combined exchange-correlation functionals ***

      SELECT CASE (TRIM(xc_fun))
      CASE ("HCTH/93")
        CALL hcth(93,rho(:,:,kg),drho(:,:,kg),ex,vxc(:,:,kg),vxcg(:,:,kg))
      CASE ("HCTH/120")
        CALL hcth(120,rho(:,:,kg),drho(:,:,kg),ex,vxc(:,:,kg),vxcg(:,:,kg))
      CASE ("HCTH/147")
        CALL hcth(147,rho(:,:,kg),drho(:,:,kg),ex,vxc(:,:,kg),vxcg(:,:,kg))
      CASE ("HCTH/407")
        CALL hcth(407,rho(:,:,kg),drho(:,:,kg),ex,vxc(:,:,kg),vxcg(:,:,kg))
      CASE ("Pade")
        CALL pade(rho(:,:,kg),ex,vxc(:,:,kg))
      END SELECT

    END DO

    CALL mp_allgather(globenv%mepos,globenv%num_pe,ng(3),vxc,globenv%group)

    IF (gradient_functional) CALL & 
        mp_allgather(globenv%mepos,globenv%num_pe,ng(3),vxcg,globenv%group)

    CALL mp_sum(ex,globenv%group)
    CALL mp_sum(ec,globenv%group)

    ex = rho_rspace(1)%pw%pw_grid%dvol*ex
    ec = rho_rspace(1)%pw%pw_grid%dvol*ec

!MK IF (gradient_functional) CALL calculate_nabla_vxcg(globenv)

    CALL timestop(0.0_wp,handle)


  END SUBROUTINE calculate_xc_potential

! *****************************************************************************

  SUBROUTINE calculate_first_density_matrix(ortho,mo_set,p,work1,work2,&
                                            context,globenv)

    TYPE(blacs_matrix_type), POINTER          :: ortho,work1,work2
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), POINTER                :: mo_set
    TYPE(real_matrix_p_type), INTENT(INOUT)   :: p
    INTEGER, INTENT(IN)                       :: context

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE calculate_first_density_matrix (MODULE qs_scf)"

!   *** Local variables ***

    TYPE(all_potential_type), POINTER :: all_potential
    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(blacs_matrix_type), POINTER  :: mo_eigenvectors
    TYPE(gth_potential_type), POINTER :: gth_potential
    TYPE(gto_basis_set_type), POINTER :: orb_basis_set

    REAL(wp) :: paa,saa,sum
    INTEGER  :: atom_a,first_sgfa,handle,iatom,ikind,iset,isgf,isgfa,ishell,&
                istat,la,last_sgfa,maxl,natom,nelec,nset,restart_unit

    INTEGER, DIMENSION(:), ALLOCATABLE :: econf

    INTEGER, DIMENSION(:), POINTER   :: atom_list,elec_conf,nshell
    INTEGER, DIMENSION(:,:), POINTER :: first_sgf,l,last_sgf

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_first_density_matrix","I","",handle)

    IF (density_guess == "RESTART") THEN

      CALL read_mo_set(mo_set,context,globenv)
      IF (level_shift /= 0.0_wp) THEN
        CALL get_mo_set(mo_set=mo_set,&
                        eigenvectors=mo_eigenvectors)
        CALL copy_blacs_to_blacs_matrix(mo_eigenvectors,ortho)
      END IF
      CALL calculate_density_matrix(mo_set,work1,context,globenv)
      CALL copy_blacs_to_sparse_matrix(work1,p%matrix,context,globenv)
      CALL symmetrise_diagonal_blocks(p%matrix)

     ELSE IF (density_guess == "CORE") THEN

!     *** It is assumed that work1 holds the upper     ***
!     *** triangle part of the core Hamiltonian matrix ***
      CALL orthogonalise_matrix(ortho,work1,work2,context,globenv)
      CALL eigensolver(work1,mo_set,ortho,work2,context,globenv)
      IF (level_shift /= 0.0_wp) THEN
        CALL get_mo_set(mo_set=mo_set,&
                        eigenvectors=mo_eigenvectors)
        CALL copy_blacs_to_blacs_matrix(mo_eigenvectors,ortho)
      END IF
      CALL calculate_density_matrix(mo_set,work1,context,globenv)
      CALL copy_blacs_to_sparse_matrix(work1,p%matrix,context,globenv)
      CALL symmetrise_diagonal_blocks(p%matrix)

    ELSE IF (density_guess == "ATOMIC") THEN

      CALL blacs_set_all(work1,0.0_wp,context,globenv)

      sum = 0.0_wp

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,maxl=maxl)

      ALLOCATE (econf(0:maxl),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"econf",(maxl + 1)*int_size)

      CALL copy_sparse_to_blacs_matrix(s%matrix,work2,context,globenv)

      DO ikind=1,nkind

        atomic_kind => atomic_kind_set(ikind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             natom=natom,&
                             atom_list=atom_list,&
                             all_potential=all_potential,&
                             gth_potential=gth_potential,&
                             orb_basis_set=orb_basis_set)

        IF (ASSOCIATED(all_potential)) THEN
          CALL get_potential(potential=all_potential,elec_conf=elec_conf)
        ELSE IF (ASSOCIATED(gth_potential)) THEN
          CALL get_potential(potential=gth_potential,elec_conf=elec_conf)
        ELSE
          CYCLE
        END IF

        maxl = SIZE(elec_conf) - 1
        econf(:) = 0
        econf(0:maxl) = elec_conf(0:maxl)

        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               nset=nset,&
                               nshell=nshell,&
                               l=l,&
                               first_sgf=first_sgf,&
                               last_sgf=last_sgf)

        DO iset=1,nset
          DO ishell=1,nshell(iset)
            la = l(ishell,iset)
            nelec = 4*la + 2
            IF (econf(la) > 0) THEN
              IF (econf(la) >= nelec) THEN
                paa = 2.0_wp
                econf(la) = econf(la) - nelec
              ELSE
                paa = REAL(2*econf(la),wp)/REAL(nelec,wp)
                econf(la) = 0
              END IF
              first_sgfa = first_sgf(ishell,iset)
              last_sgfa = last_sgf(ishell,iset)
              DO isgfa=first_sgfa,last_sgfa
                DO iatom=1,natom
                  atom_a = atom_list(iatom)
                  isgf = particle_set(atom_a)%first_sgf + isgfa - 1
                  CALL blacs_set_element(work1,isgf,isgf,paa,context,globenv)
                  CALL blacs_get_element(work2,isgf,isgf,saa,context,globenv)
                  sum = sum + paa*saa
                END DO
              END DO
            END IF
          END DO
        END DO

      END DO

      DEALLOCATE (econf,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"econf")

      CALL mp_sum(sum,globenv%group)

      saa = REAL(nelectron,wp)/sum

      DO ikind=1,nkind

        atomic_kind => atomic_kind_set(ikind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             natom=natom,&
                             atom_list=atom_list,&
                             orb_basis_set=orb_basis_set)

        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               nset=nset,&
                               nshell=nshell,&
                               l=l,&
                               first_sgf=first_sgf,&
                               last_sgf=last_sgf)

        DO iset=1,nset
          DO ishell=1,nshell(iset)
            first_sgfa = first_sgf(ishell,iset)
            last_sgfa = last_sgf(ishell,iset)
            DO isgfa=first_sgfa,last_sgfa
              DO iatom=1,natom
                atom_a = atom_list(iatom)
                isgf = particle_set(atom_a)%first_sgf + isgfa - 1
                CALL blacs_get_element(work1,isgf,isgf,paa,context,globenv)
                paa = paa*saa
                CALL blacs_set_element(work1,isgf,isgf,paa,context,globenv)
              END DO
            END DO
          END DO
        END DO
      END DO

      CALL copy_blacs_to_sparse_matrix(work1,p%matrix,context,globenv)

    ELSE

      CALL stop_program("SUBROUTINE calculate_first_density_matrix "//&
                        "(MODULE qs_scf)",&
                        "An invalid keyword for the initial density "//&
                        "guess was specified")

    END IF

  END SUBROUTINE calculate_first_density_matrix

! *****************************************************************************

  SUBROUTINE read_scf_parameters(globenv)

!   Purpose: Read the parameters for the SCF run.

!   History: - Creation (15.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    CHARACTER(LEN=40) :: keyword

!   ---------------------------------------------------------------------------

!   *** Load the default values ***

    charge = 0
    density_guess = "ATOMIC"
    eps_eigval = 1.0E-5_wp
    eps_scf = 1.0E-5_wp
    eps_diis = 0.1_wp
    level_shift = 0.0_wp
    max_diis = 4
    max_scf = 50
    nprow = 0
    npcol = 0
    nrebuild = 1
    nrow_block = 32
    ncol_block = 32
    p_mix = 0.4_wp
    smear = 0.0_wp
    work_syevx = 1.0_wp
    use_cholesky = .true.

    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label="SCF")

    DO WHILE (test_object(newline=.TRUE.) /= "EOS")
      CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
      SELECT CASE (TRIM(keyword))
      CASE ("CHARGE")
        CALL read_object(charge)
      CASE ("DENSITY_GUESS","SCF_GUESS","GUESS")
        CALL read_object(density_guess,lower_to_upper=.TRUE.)
      CASE ("DENSITY_MIXING","MIXING")
        CALL read_object(p_mix)
      CASE ("EPS_DIIS")
        CALL read_object(eps_diis)
      CASE ("EPS_EIGVAL")
        CALL read_object(eps_eigval)
        use_cholesky = .false.
      CASE ("CHOLESKY_ON")
        use_cholesky = .true.
      CASE ("CHOLESKY_OFF")
        use_cholesky = .false.
      CASE ("EPS_SCF")
        CALL read_object(eps_scf)
      CASE ("LEVEL_SHIFT","LSHIFT")
        CALL read_object(level_shift)
      CASE ("MAX_DIIS")
        CALL read_object(max_diis)
      CASE ("MAX_SCF")
        CALL read_object(max_scf)
      CASE ("NREBUILD")
        CALL read_object(nrebuild)
        nrebuild = MAX(1,nrebuild)
      CASE ("BLOCKSIZE")
        CALL read_object(nrow_block)
        ncol_block = nrow_block
      CASE ("PROCESS_GRID")
        CALL read_object(nprow)
        CALL read_object(npcol)
      CASE ("SMEAR")
        CALL read_object(smear)
      CASE ("WORK_SYEVX")
        CALL read_object(work_syevx)
        work_syevx = MIN(MAX(0.0_wp,work_syevx),1.0_wp)
      END SELECT
    END DO

    CALL finish_parser()

  END SUBROUTINE read_scf_parameters

! *****************************************************************************

  SUBROUTINE write_scf_parameters(globenv)

!   Purpose: Write the parameters for the SCF run to the output unit.

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    INTEGER :: output_unit

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (max_scf < 1) RETURN

    IF (globenv%print%scf) THEN
      output_unit = globenv%scr
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A,/)") "SCF PARAMETERS"
      WRITE (UNIT=output_unit,&
             FMT="(T3,A,/,4(/,T3,A,I5),3(/,T3,A,ES9.2),4(/,T3,A,F5.2))")&
        "Density guess: "//TRIM(density_guess),&
        "charge:            ",charge,&
        "nrebuild:          ",nrebuild,&
        "max_scf:           ",max_scf,&
        "max_diis:          ",max_diis,&
        "eps_scf:           ",eps_scf,&
        "eps_diis:          ",eps_diis,&
        "eps_eigval:        ",eps_eigval,&
        "p_mix:             ",p_mix,&
        "work_syevx:        ",work_syevx,&
        "level_shift [a.u.]:",level_shift,&
        "smear [a.u.]:      ",smear
    END IF

  END SUBROUTINE write_scf_parameters

! *****************************************************************************

END MODULE qs_scf
