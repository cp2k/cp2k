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
!!     1) Joost VandeVondele (02.2002)
!!           added code for : incremental (pab and gvg) update
!!                            initialisation (init_cube, l_info)
!!     2) Joost VandeVondele (02.2002)
!!           called the poisson code of the classical part
!!           this takes into account the spherical cutoff and allows for
!!           isolated systems
!!
!!   SOURCE
!******************************************************************************

MODULE qs_scf

! *****************************************************************************

  USE kinds, ONLY: wp => dp

  USE cell_parameters,       ONLY: abc,cell,cutoff
  USE atomic_kinds,          ONLY: gth_potential_present,kind_info,nkind
  USE atoms,                 ONLY: atom_info
  USE basis_set_types,       ONLY: gto_basis_set_type,maxlcgf
  USE blacs,                 ONLY: allocate_blacs_matrix,&
                                   blacs_gemm,&
                                   blacs_matrix_type,&
                                   blacs_set_all,&
                                   blacs_set_element,&
                                   blacs_get_element,&
                                   blacs_syevx,&
                                   blacs_symm,&
                                   blacs_syrk,&
                                   copy_blacs_to_sparse_matrix,&
                                   copy_sparse_to_blacs_matrix,&
                                   finish_blacs,&
                                   get_blacs_matrix_info,&
                                   power_blacs_matrix,&
                                   replicate_blacs_matrix,&
                                   start_blacs,&
                                   symmetrise_blacs_matrix
  USE coefficient_types,     ONLY: coeff_allocate,&
                                   coeff_type,&
                                   coeff_zero,&
                                   coeff_copy
  USE core_energies,         ONLY: calculate_ecore,&
                                   calculate_ecore_overlap,&
                                   calculate_ecore_self
  USE core_hamiltonian,      ONLY: build_core_hamiltonian_matrix,h,s
  USE diis,                  ONLY: eps_diis,max_diis,scf_diis
  USE functionals
  USE global_types,          ONLY: global_environment_type
  USE input_utilities,       ONLY: finish_input_session,&
                                   read_object,&
                                   search,&
                                   start_input_session
  USE matrix_types,          ONLY: allocate_matrix,&
                                   copy_matrix,&
                                   deallocate_matrix,&
                                   first_block_node,&
                                   get_block_node,&
                                   get_matrix_info,&
                                   next_block_node,&
                                   real_block_node_type,&
                                   real_matrix_set_type,&
                                   replicate_matrix_structure,&
                                   symmetrise_diagonal_blocks
  USE memory_utilities,      ONLY: reallocate
  USE message_passing,       ONLY: mp_max,mp_sum,mp_sync
  USE method_specifications, ONLY: xc_info
  USE mo_types,              ONLY: allocate_mo_set,&
                                   mo_set_type,&
                                   read_mo_set,&
                                   write_mo_set
  USE om_utilities,          ONLY: write_blacs_matrix,&
                                   write_spherical_matrix
  USE pw_grid_types,         ONLY: HALFSPACE,pw_grid_type
  USE pw_grids,              ONLY: pw_find_cutoff,&
                                   pw_grid_construct,&
                                   pw_grid_setup
  USE pws,                   ONLY: init_pw_poisson_solver,pw_poisson_solver
  USE greens_fn,             ONLY: pw_green_fn_init 
  USE hartree,               ONLY: calculate_hartree
  USE pw_types,              ONLY: COMPLEXDATA3D,REALDATA3D,REALSPACE,&
                                   RECIPROCALSPACE,pw_fft_wrap,pw_type
  USE string_utilities,      ONLY: uppercase
  USE termination,           ONLY: stop_program
  USE timings,               ONLY: timeset,timestop
  USE timesl,                ONLY: cputime
  USE cube_utils,            ONLY: cube_info_type,init_cube_info, & 
                                   destroy_cube_info
  USE l_utils,               ONLY: l_info_type,init_l_info, & 
                                   destroy_l_info

  IMPLICIT NONE

  PRIVATE

  TYPE(blacs_matrix_type)    :: ortho,scf_work1,scf_work2
  TYPE(coeff_type)           :: drho_rspace,rho_gspace,rho_rspace,v_rspace,&
                                vxcg_rspace,work_gspace,work_rspace,&
                                rho_rspace_old,v_rspace_old
  TYPE(mo_set_type)          :: alpha_mo
  TYPE(real_matrix_set_type) :: ks,p,p_old,gvg
  TYPE(pw_grid_type)         :: pw_grid
  TYPE(cube_info_type)       :: cube_info
  TYPE(l_info_type)          :: l_info

  CHARACTER(LEN=10) :: density_guess = "CORE"
  REAL(wp)          :: ec = 0.0_wp,&
                       ecore = 0.0_wp,&
                       ecore_overlap = 0.0_wp,&
                       ecore_self = 0.0_wp,&
                       ehartree = 0.0_wp,&
                       eps_eigval = 1.0E-5_wp,&
                       eps_scf = 1.0E-5_wp,&
                       etotal = 0.0_wp,&
                       ex = 0.0_wp,&
                       p_mix = 0.4_wp,&
                       total_rho_core_rspace = 0.0_wp,&
                       total_rho_elec_rspace = 0.0_wp,&
                       total_rho_gspace = 0.0_wp,&
                       total_rho_rspace = 0.0_wp,&
                       work_syevx = 0.0_wp
  INTEGER           :: max_scf = 50,&
                       nrebuild_rho = 5, &
                       nrebuild_gvg = 5, &
                       nelectron = 0,&
                       nrow_block = 32,&
                       ncol_block = 32,&
                       nprow = 0,&
                       npcol = 0,&
                       additional_charge = 0

  logical           :: rebuild_rho = .true.
  logical           :: rebuild_gvg = .true.

  TYPE(coeff_type), DIMENSION(3) :: nabla_rho_rspace

  REAL(wp), DIMENSION(:,:), POINTER :: g_vector

! *** Public variables ***

  PUBLIC :: nelectron

! *** Public subroutines ***

  PUBLIC :: read_scf_parameters,&
            scf,&
            write_scf_parameters

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE scf(globenv)

    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE scf (MODULE qs_scf)"

!   *** Local variables ***

    REAL(wp) :: delta,diis_error,t1,t2
    REAL(wp) :: maxradius
    PARAMETER (maxradius=40.0D0)
    integer  :: maxl
    PARAMETER (maxl=4)
    INTEGER  :: context,handle,iscf,output_unit
    LOGICAL  :: diis_step,ionode

!   ---------------------------------------------------------------------------

!   *** Quick return, if no SCF iteration is requested ***

    IF (max_scf < 1) RETURN

    CALL timeset("scf","I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr

    CALL init_functionals()

    CALL init_grid(globenv)
  
    CALL init_cube_info(cube_info,rho_rspace%pw%pw_grid%dr(:),maxradius)

    CALL init_l_info(l_info,maxl)

    IF (ionode.AND.globenv%print%scf) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
        "SCF WAVEFUNCTION OPTIMIZATION"
    END IF

    CALL init_scf_run(alpha_mo,context,globenv)

    etotal = 0.0_wp
    iscf = 0
    diis_step = .FALSE.
    rebuild_rho = .true.
    rebuild_gvg = .true.

    IF (ionode.AND.globenv%print%scf) THEN
      WRITE (UNIT=output_unit,&
             FMT="(/,T3,A,T9,A,T34,A,T49,A,T68,A,/,T3,A)")&
        "Step","Update method","Time","Convergence","Total energy",&
        REPEAT("-",77)
    END IF

    DO

      iscf = iscf + 1

      t1 = cputime()

      IF (globenv%print%density_matrix.AND.globenv%print%each_scf_step) THEN
        CALL write_spherical_matrix(p%matrix,4,6,globenv)
      END IF

      CALL build_kohn_sham_matrix(globenv)

      if (MODULO(iscf,nrebuild_rho).eq.0) then
         rebuild_rho=.true.
      else
         rebuild_rho=.false.
      endif
      if (MODULO(iscf,nrebuild_gvg).eq.0) then
         rebuild_gvg=.true.
      else
         rebuild_gvg=.false.
      endif

      IF (globenv%print%kohn_sham_matrix.AND.globenv%print%each_scf_step) THEN
        CALL write_spherical_matrix(ks%matrix,4,6,globenv)
      END IF

      CALL copy_sparse_to_blacs_matrix(ks%matrix,scf_work1,context,globenv)
!MK if only gemm is used
!     CALL symmetrise_blacs_matrix(scf_work1,scf_work2,context,globenv)

      IF (iscf > 1) THEN
        CALL scf_diis(alpha_mo,scf_work1,scf_work2,delta,diis_error,diis_step,&
                      context,globenv)
      END IF

      CALL orthogonalise_matrix(ortho,scf_work1,scf_work2,context,globenv)

      CALL eigensolver(scf_work1,alpha_mo,ortho,scf_work2,context,globenv)

      IF (globenv%print%each_scf_step) THEN
        CALL write_mo_set(alpha_mo,4,6,context,globenv)
      END IF

      CALL calculate_density_matrix(alpha_mo,scf_work1,context,globenv)

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
        EXIT
      ELSE IF (iscf == max_scf) THEN
        IF (ionode.AND.globenv%print%scf) THEN
          WRITE(UNIT=output_unit,FMT="(/,T3,A,/)")&
            "*** SCF run NOT converged ***"
        END IF
        EXIT
      END IF

    END DO

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

    IF (.NOT.diis_step) CALL copy_matrix(ks%matrix,p%matrix)

    CALL write_mo_set(alpha_mo,4,6,context,globenv)

    IF (globenv%print%density_matrix) THEN
      CALL write_spherical_matrix(p%matrix,4,6,globenv)
    END IF

    IF (globenv%print%kohn_sham_matrix) THEN
      CALL build_kohn_sham_matrix(globenv)
      CALL write_spherical_matrix(ks%matrix,4,6,globenv)
    END IF

!   *** Write restart file ***

    CALL write_mo_set(alpha_mo,context,globenv)

    CALL destroy_cube_info(cube_info)
    CALL destroy_l_info(l_info)

    CALL finish_blacs(context,globenv)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE scf

! *****************************************************************************

  SUBROUTINE build_kohn_sham_matrix(globenv)

    USE collocate_density,   ONLY: calculate_rho_core,&
                                   calculate_rho_elec,&
                                   calculate_total_rho
    USE integrate_potential, ONLY: integrate_v_rspace
    USE print_keys,          ONLY: DEBUG

    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    INTEGER :: handle,nproduct,output_unit
    LOGICAL :: ionode

!   ---------------------------------------------------------------------------

    CALL timeset("build_kohn_sham_matrix","I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr

    CALL calculate_ecore(h,p,ecore,globenv)

    if (rebuild_rho) then
      CALL coeff_zero(rho_rspace)
    else
      CALL coeff_copy(rho_rspace_old,rho_rspace)
    endif

    CALL calculate_rho_elec(p,p_old,rebuild_rho,rho_rspace,cube_info,l_info, &
                  total_rho_elec_rspace,nproduct,globenv)

    CALL coeff_copy(rho_rspace,rho_rspace_old) ! keep a copy of the old density

    IF (ionode.AND.globenv%print%total_densities) THEN
      WRITE (UNIT=output_unit,FMT="(/,T3,A,I10)")&
        "Number of collocated products (r-space):",nproduct
    END IF

    IF (xc_info%gradient_functional) THEN
      CALL calculate_nabla_rho_elec(globenv)
      CALL coeff_zero(vxcg_rspace)
    END IF

    CALL coeff_zero(v_rspace)

    CALL calculate_xc_potential(globenv)

    CALL calculate_rho_core(rho_rspace,cube_info,l_info,total_rho_rspace,globenv)

    total_rho_core_rspace = total_rho_rspace - total_rho_elec_rspace

    CALL pw_fft_wrap(rho_rspace%pw,rho_gspace%pw,&
                     debug=(globenv%print%level == DEBUG))

    total_rho_gspace = calculate_total_rho(rho_gspace)

    IF (ionode.AND.globenv%print%total_densities) THEN
      WRITE (UNIT=output_unit,FMT="(/,(T3,A,T40,2F20.10))")&
        "Total electronic density (r-space): ",&
        total_rho_elec_rspace,total_rho_elec_rspace + REAL(nelectron,wp),&
        "Total core charge density (r-space):",&
        total_rho_core_rspace,total_rho_core_rspace - REAL(nelectron+additional_charge,wp)
      WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
        "Total charge density (r-space):     ",total_rho_rspace,&
        "Total charge density (g-space):     ",total_rho_gspace
    END IF

    CALL calculate_hartree(rho_gspace,ehartree,work_gspace)
    ! ehartree is returned using a slightly different def. of the volume or so
    ehartree=ehartree*rho_gspace%pw%pw_grid%vol/rho_gspace%pw%pw_grid%dvol

    etotal = ecore_overlap + ecore_self + ecore + ehartree + ex + ec

    IF (ionode.AND.globenv%print%scf_energies) THEN
      WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
        "Overlap energy of the core charge distribution:",ecore_overlap,&
        "Self energy of the core charge distribution:   ",ecore_self,&
        "Core Hamiltonian energy:                       ",ecore,&
        "Hartree energy:                                ",ehartree,&
        "Exchange-correlation energy:                   ",ex + ec
    END IF

    !CALL pw_fft_wrap(rho_gspace%pw,work_rspace%pw,&
    !                 debug=(globenv%print%level == DEBUG))
    CALL pw_fft_wrap(work_gspace%pw,work_rspace%pw,&
                     debug=(globenv%print%level == DEBUG))

    v_rspace%pw%cr3d(:,:,:) = v_rspace%pw%pw_grid%dvol*&
                              (v_rspace%pw%cr3d(:,:,:) +&
                               work_rspace%pw%cr3d(:,:,:))

    CALL copy_matrix(h%matrix,ks%matrix)

    CALL integrate_v_rspace(p,v_rspace,v_rspace_old,cube_info,l_info,rebuild_gvg,&
                             gvg,ks,nproduct,globenv)

    IF (ionode.AND.globenv%print%total_densities) THEN
      WRITE (UNIT=output_unit,FMT="(/,T3,A,I10)")&
        "Number of integrated products (r-space):",nproduct
    END IF

    CALL copy_matrix(p%matrix,p_old%matrix) ! keep a copy of the old density matrix

    CALL timestop(0.0_wp,handle)

 END SUBROUTINE build_kohn_sham_matrix

! *****************************************************************************

  SUBROUTINE init_grid(globenv)

    USE mathconstants,   ONLY: twopi

    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    INTEGER :: handle,i,ig,jg,output_unit
    LOGICAL :: ionode

    REAL(wp), DIMENSION(3) :: dg,dr
    INTEGER, DIMENSION(3)  :: lb_grid,ng,ub_grid

!   ---------------------------------------------------------------------------

    CALL timeset("init_grid","I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr

    CALL pw_grid_construct(pw_grid)

    pw_grid%grid_span = HALFSPACE

    IF (ionode.AND.globenv%print%pw_grid_information) THEN
      CALL pw_grid_setup(cell,pw_grid,cutoff,&
                         info=output_unit,&
                         orthorhombic=.TRUE.)
    ELSE
      CALL pw_grid_setup(cell,pw_grid,cutoff,&
                         orthorhombic=.TRUE.)
    END IF

    CALL pw_green_fn_init ( cell % green, pw_grid  )
    CALL init_pw_poisson_solver ( cell % green )

    CALL coeff_allocate(work_rspace,pw_grid,REALDATA3D)
    work_rspace%pw%in_space = REALSPACE

    CALL coeff_allocate(work_gspace,pw_grid,COMPLEXDATA3D)
    work_gspace%pw%in_space = RECIPROCALSPACE


    CALL coeff_allocate(rho_rspace,pw_grid,REALDATA3D)
    rho_rspace%pw%in_space = REALSPACE
    CALL coeff_allocate(rho_rspace_old,pw_grid,REALDATA3D)
    rho_rspace_old%pw%in_space = REALSPACE

    IF (xc_info%gradient_functional) THEN
      DO i=1,3
        CALL coeff_allocate(nabla_rho_rspace(i),pw_grid,REALDATA3D)
        nabla_rho_rspace(i)%pw%in_space = REALSPACE
      END DO
      CALL coeff_allocate(drho_rspace,pw_grid,REALDATA3D)
      drho_rspace%pw%in_space = REALSPACE
      CALL coeff_allocate(vxcg_rspace,pw_grid,REALDATA3D)
      vxcg_rspace%pw%in_space = REALSPACE
    END IF

    CALL coeff_allocate(rho_gspace,pw_grid,COMPLEXDATA3D)
    rho_gspace%pw%in_space = RECIPROCALSPACE

    CALL coeff_allocate(v_rspace,pw_grid,REALDATA3D)
    v_rspace%pw%in_space = REALSPACE
    CALL coeff_allocate(v_rspace_old,pw_grid,REALDATA3D)
    v_rspace_old%pw%in_space = REALSPACE

!   *** Build g vectors ***

    dr(:) = rho_rspace%pw%pw_grid%dr(:)
    ng(:) = rho_rspace%pw%pw_grid%npts(:)
    dg(:) = twopi/(ng(:)*dr(:))
    lb_grid(:) = rho_rspace%pw%pw_grid%bounds(1,:)
    ub_grid(:) = rho_rspace%pw%pw_grid%bounds(2,:)

    g_vector => reallocate(g_vector,1,MAXVAL(ng(:)),1,3)

    DO i=1,3
      DO ig=lb_grid(i),ub_grid(i)
        IF (ig < 0) THEN
          jg = ig + 1 + ng(i)
        ELSE
          jg = ig + 1
        END IF
        g_vector(jg,i) = REAL(ig,wp)*dg(i)
      END DO
    END DO

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE init_grid

! *****************************************************************************

  SUBROUTINE init_scf_run(mo,context,globenv)

!   Purpose: Initialise a SCF run.

!   History: - Creation (30.04.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), INTENT(OUT)            :: mo
    INTEGER, INTENT(OUT)                      :: context

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE init_scf_run (MODULE qs_scf)"

!   *** Local variables ***

    INTEGER :: handle,homo,ikind,nao,ndep,nmo,output_unit
    LOGICAL :: ionode

    INTEGER, DIMENSION(0:globenv%num_pe-1) :: nblock_pe,nelement_pe

!   ---------------------------------------------------------------------------

    CALL timeset("init_scf_run","I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr

!   *** Initialise BLACS ***

    CALL start_blacs(nprow,npcol,context,globenv)

!   *** Calculate the number of electrons ***

    nelectron = - additional_charge

    DO ikind=1,nkind
      nelectron = nelectron + kind_info(ikind)%natom*kind_info(ikind)%zeff
    END DO

    IF (nelectron.lt.1) THEN
      CALL stop_program(routine,"too few electrons hanging around")
    ENDIF

    IF (MODULO(nelectron,2) == 0) THEN
      homo = nelectron/2
    ELSE
      CALL stop_program(routine,"Odd number of electrons")
    END IF

!   *** Get the dimension of the full SCF matrices, ***
!   *** i.e. the total number of atomic orbitals    ***

    CALL get_matrix_info(matrix=h%matrix,nrow=nao)

    IF (globenv%print%mo_eigenvectors) THEN
      nmo = nao
    ELSE
      nmo = homo
    END IF

!   *** Allocate the distributed MO eigenvectors ***

    CALL allocate_mo_set(mo,nao,nmo,nrow_block,ncol_block,context,globenv)

    mo%homo = homo
    mo%occupation_numbers(1:homo) = 2.0_wp

!   *** Get BLACS block size of the MO eigenvector matrix      ***
!   *** which has to fit to the other distributed SCF matrices ***

    CALL get_blacs_matrix_info(matrix=mo%eigenvectors,&
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
    NULLIFY (gvg%matrix)
    CALL replicate_matrix_structure(h%matrix,gvg%matrix,"POTENTIAL MATRIX")

    CALL calculate_ecore_self(ecore_self)
    CALL calculate_ecore_overlap(ecore_overlap,globenv)

    IF (ionode.AND.globenv%print%scf_energies) THEN
      WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
        "Self energy of the core charge distribution:   ",ecore_self,&
        "Overlap energy of the core charge distribution:",ecore_overlap
    END IF

    CALL calculate_ortho_matrix(ortho,scf_work2,ndep,context,globenv)

    IF (ionode.AND.globenv%print%scf) THEN
      WRITE (UNIT=output_unit,FMT="(/,(T3,A,I10))")&
        "Number of electrons:                    ",nelectron,&
        "Number of occupied orbitals:            ",homo,&
        "Number of orbital functions:            ",nao,&
        "Number of independent orbital functions:",nao - ndep
    END IF

    IF (globenv%print%ortho_matrix) THEN
      CALL write_blacs_matrix(ortho,4,6,context,globenv)
    END IF

    CALL calculate_first_density_matrix(ortho,mo,p,scf_work1,scf_work2,&
                                        context,globenv)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE init_scf_run

! *****************************************************************************

  SUBROUTINE calculate_ortho_matrix(ortho,work,ndep,context,globenv)

!   Purpose: Calculate the orthogonalization matrix (S**(-1/2))

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(blacs_matrix_type), INTENT(OUT)      :: ortho
    TYPE(blacs_matrix_type), INTENT(INOUT)    :: work
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

!   Purpose: Orthogonaliseation matrix (S**(-1/2))

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), INTENT(IN)       :: ortho
    TYPE(blacs_matrix_type), INTENT(INOUT)    :: ks,work
    TYPE(global_environment_type), INTENT(IN) :: globenv
    INTEGER, INTENT(IN)                       :: context

!   *** Local variables ***

    INTEGER :: handle,nao

!   ---------------------------------------------------------------------------

    CALL timeset("orthogonalise_matrix","I","",handle)

    CALL get_blacs_matrix_info(matrix=ks,nrow_global=nao)
    CALL blacs_symm("L","U",nao,nao,1.0_wp,ks,ortho,0.0_wp,work,context,&
                    globenv)
!MK CALL blacs_gemm("T","N",nao,nao,nao,1.0_wp,ks,ortho,0.0_wp,work,context,&
!MK                 globenv)
    CALL blacs_gemm("N","N",nao,nao,nao,1.0_wp,ortho,work,0.0_wp,ks,context,&
                    globenv)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE orthogonalise_matrix

! *****************************************************************************

  SUBROUTINE eigensolver(ks,mo,ortho,work,context,globenv)

!   Purpose: Diagonalise the Kohn-Sham matrix to get a new set of MO eigen-
!            vectors and MO eigenvalues.

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), INTENT(IN)       :: ortho
    TYPE(blacs_matrix_type), INTENT(INOUT)    :: ks,work
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), INTENT(INOUT)          :: mo
    INTEGER, INTENT(IN)                       :: context

!   *** Local variables ***

    INTEGER :: handle,nao,nmo

!   ---------------------------------------------------------------------------

    CALL timeset("eigensolver","I","",handle)

    CALL get_blacs_matrix_info(matrix=mo%eigenvectors,&
                               nrow_global=nao,&
                               ncol_global=nmo)

    CALL blacs_syevx(ks,work,mo%eigenvalues,nmo,work_syevx,context,globenv)

!MK CALL blacs_gemm("T","N",nao,nmo,nao,1.0_wp,ortho,work,0.0_wp,&
!MK                 mo%eigenvectors,context,globenv)
    CALL blacs_symm("L","U",nao,nmo,1.0_wp,ortho,work,0.0_wp,mo%eigenvectors,&
                    context,globenv)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE eigensolver

! *****************************************************************************

  SUBROUTINE calculate_density_matrix(mo,density,context,globenv)

!   Purpose: Calculate the density matrix from the MO eigenvectors and the
!            MO occupation numbers.

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(blacs_matrix_type), INTENT(OUT)      :: density
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), INTENT(IN)             :: mo
    INTEGER, INTENT(IN)                       :: context

!   *** Local variables ***

    INTEGER :: handle

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_density_matrix","I","",handle)

    CALL blacs_syrk("U","N",mo%homo,2.0_wp,mo%eigenvectors,0.0_wp,density,&
                    context,globenv)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_density_matrix

! *****************************************************************************

  SUBROUTINE density_mixing(new_density,old_density,p_mix,delta,globenv)

!   Purpose: Perform a density mixing of the old (last SCF iteration) and the
!            new (current) density matrix.

!   History: - Creation (01.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(real_matrix_set_type), INTENT(INOUT) :: new_density,old_density
    REAL(wp), INTENT(IN)                      :: p_mix
    REAL(wp), INTENT(OUT)                     :: delta

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: new_p_block_node,old_p_block_node
    REAL(wp)                            :: r_mix
    INTEGER                             :: handle,i,iblock_row,j,nblock_row

    REAL(wp), DIMENSION(:,:), POINTER :: new_p_block,old_p_block

!   ---------------------------------------------------------------------------

    CALL timeset("density_mixing","I","",handle)

    delta = 0.0_wp

    r_mix = 1.0_wp - p_mix

    CALL get_matrix_info(matrix=old_density%matrix,nblock_row=nblock_row)

    DO iblock_row=1,nblock_row

      new_p_block_node => first_block_node(matrix=new_density%matrix,&
                                           block_row=iblock_row)
      old_p_block_node => first_block_node(matrix=old_density%matrix,&
                                           block_row=iblock_row)

      DO WHILE (ASSOCIATED(old_p_block_node))

        CALL get_block_node(block_node=new_p_block_node,&
                            block=new_p_block)
        CALL get_block_node(block_node=old_p_block_node,&
                            block=old_p_block)

        DO j=1,SIZE(new_p_block,2)
          DO i=1,SIZE(new_p_block,1)
            new_p_block(i,j) = new_p_block(i,j) - old_p_block(i,j)
            delta = MAX(delta,ABS(new_p_block(i,j)))
            old_p_block(i,j) = old_p_block(i,j) + p_mix*new_p_block(i,j)
            new_p_block(i,j) = old_p_block(i,j) + r_mix*new_p_block(i,j)
          END DO
        END DO

        new_p_block_node => next_block_node(new_p_block_node)
        old_p_block_node => next_block_node(old_p_block_node)

      END DO

    END DO

    CALL mp_max(delta,globenv%group)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE density_mixing

! *****************************************************************************

  SUBROUTINE calculate_nabla_rho_elec(globenv)

    USE print_keys, ONLY: DEBUG

    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    COMPLEX(wp) :: c
    INTEGER     :: handle,i,ig,jg,kg

    INTEGER, DIMENSION(3) :: ng

    COMPLEX(wp), DIMENSION(:,:,:), POINTER :: rho,work

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_drho_elec","I","",handle)

    CALL pw_fft_wrap(rho_rspace%pw,rho_gspace%pw,&
                     debug=(globenv%print%level == DEBUG))

    ng(:) = rho_rspace%pw%pw_grid%npts(:)

    drho_rspace%pw%cr3d(:,:,:) = 0.0_wp

    rho => rho_gspace%pw%cc3d(:,:,:)
    work => work_gspace%pw%cc3d(:,:,:)

    DO i=1,3

      SELECT CASE (i)
      CASE (1)
        DO kg=1,ng(3)
          DO jg=1,ng(2)
            DO ig=1,ng(1)
              c = CMPLX(0.0_wp,g_vector(ig,1),wp)
              work(ig,jg,kg) = c*rho(ig,jg,kg)
            END DO
          END DO
        END DO
      CASE (2)
        DO kg=1,ng(3)
          DO jg=1,ng(2)
            c = CMPLX(0.0_wp,g_vector(jg,2),wp)
            DO ig=1,ng(1)
              work(ig,jg,kg) = c*rho(ig,jg,kg)
            END DO
          END DO
        END DO
      CASE (3)
        DO kg=1,ng(3)
          c = CMPLX(0.0_wp,g_vector(kg,3),wp)
          DO jg=1,ng(2)
            DO ig=1,ng(1)
              work(ig,jg,kg) = c*rho(ig,jg,kg)
            END DO
          END DO
        END DO
      END SELECT

      CALL pw_fft_wrap(work_gspace%pw,nabla_rho_rspace(i)%pw,&
                       debug=(globenv%print%level == DEBUG))

      drho_rspace%pw%cr3d(:,:,:) = drho_rspace%pw%cr3d(:,:,:) +&
                                   nabla_rho_rspace(i)%pw%cr3d(:,:,:)*&
                                   nabla_rho_rspace(i)%pw%cr3d(:,:,:)

    END DO

    drho_rspace%pw%cr3d(:,:,:) = SQRT(drho_rspace%pw%cr3d(:,:,:))

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_nabla_rho_elec

! *****************************************************************************

  SUBROUTINE calculate_xc_potential(globenv)

    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    INTEGER :: handle,kg

    INTEGER, DIMENSION(3) :: ng

    REAL(wp), DIMENSION(:,:,:), POINTER :: drho,rho,vxc,vxcg

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_xc_potential","I","",handle)

    ng(:) = rho_rspace%pw%pw_grid%npts(3)

    rho => rho_rspace%pw%cr3d(:,:,:)
    vxc => v_rspace%pw%cr3d(:,:,:)

    IF (xc_info%gradient_functional) THEN
      drho => drho_rspace%pw%cr3d(:,:,:)
      vxcg => vxcg_rspace%pw%cr3d(:,:,:)
    END IF

    ex = 0.0_wp
    ec = 0.0_wp

    DO kg=1,ng(3)

      IF (globenv%mepos /= MODULO(kg,globenv%num_pe)) CYCLE

!     *** Exchange functionals ***

      SELECT CASE (xc_info%x_functional)
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

      SELECT CASE (xc_info%c_functional)
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

      SELECT CASE (xc_info%x_functional)
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

    CALL mp_sum(vxc,globenv%group)
    IF (xc_info%gradient_functional) CALL mp_sum(vxcg,globenv%group)

    CALL mp_sum(ex,globenv%group)
    CALL mp_sum(ec,globenv%group)

    ex = rho_rspace%pw%pw_grid%dvol*ex
    ec = rho_rspace%pw%pw_grid%dvol*ec

    IF (xc_info%gradient_functional) THEN
      CALL calculate_nabla_vxcg(globenv)
    END IF

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_xc_potential

! *****************************************************************************

  SUBROUTINE calculate_nabla_vxcg(globenv)

    USE print_keys, ONLY: DEBUG

    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    COMPLEX(wp) :: c
    INTEGER     :: handle,i,ig,jg,kg

    INTEGER, DIMENSION(3) :: ng

    COMPLEX(wp), DIMENSION(:,:,:), POINTER :: work

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_nabla_vxcg","I","",handle)

    ng(:) = rho_rspace%pw%pw_grid%npts(:)

    work => work_gspace%pw%cc3d(:,:,:)

    DO i=1,3

      work_rspace%pw%cr3d(:,:,:) = vxcg_rspace%pw%cr3d(:,:,:)*&
                                   nabla_rho_rspace(i)%pw%cr3d(:,:,:)

      CALL pw_fft_wrap(work_rspace%pw,work_gspace%pw,&
                       debug=(globenv%print%level == DEBUG))

      SELECT CASE (i)
      CASE (1)
        DO kg=1,ng(3)
          DO jg=1,ng(2)
            DO ig=1,ng(1)
              c = CMPLX(0.0_wp,g_vector(ig,1),wp)
              work(ig,jg,kg) = c*work(ig,jg,kg)
            END DO
          END DO
        END DO
      CASE (2)
        DO kg=1,ng(3)
          DO jg=1,ng(2)
            c = CMPLX(0.0_wp,g_vector(jg,2),wp)
            DO ig=1,ng(1)
              work(ig,jg,kg) = c*work(ig,jg,kg)
            END DO
          END DO
        END DO
      CASE (3)
        DO kg=1,ng(3)
          c = CMPLX(0.0_wp,g_vector(kg,3),wp)
          DO jg=1,ng(2)
            DO ig=1,ng(1)
              work(ig,jg,kg) = c*work(ig,jg,kg)
            END DO
          END DO
        END DO
      END SELECT

      CALL pw_fft_wrap(work_gspace%pw,work_rspace%pw,&
                       debug=(globenv%print%level == DEBUG))

      v_rspace%pw%cr3d(:,:,:) = v_rspace%pw%cr3d(:,:,:) -&
                                work_rspace%pw%cr3d(:,:,:)

    END DO

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_nabla_vxcg

! *****************************************************************************

  SUBROUTINE calculate_first_density_matrix(ortho,mo,p,work1,work2,context,& 
                                            globenv) 
 
    TYPE(blacs_matrix_type), INTENT(IN)       :: ortho 
    TYPE(blacs_matrix_type), INTENT(INOUT)    :: work1,work2 
    TYPE(global_environment_type), INTENT(IN) :: globenv 
    TYPE(mo_set_type), INTENT(INOUT)          :: mo 
    TYPE(real_matrix_set_type), INTENT(INOUT) :: p 
    INTEGER, INTENT(IN)                       :: context 
 
!   *** Local variables *** 
 
    TYPE(gto_basis_set_type), POINTER :: basis_set_a 
 
    REAL(wp) :: paa,saa,sum 
    INTEGER  :: atoma,first_sgfa,handle,iatom,ikind,iset,isgf,isgfa,ishell,la,& 
                last_sgfa,nelec,restart_unit 
 
    INTEGER, DIMENSION(:), POINTER :: elec_conf 
 
!   --------------------------------------------------------------------------- 
 
    CALL timeset("calculate_first_density_matrix","I","",handle) 
 
    CALL blacs_set_all(work1,0.0_wp,context,globenv) 
 
    IF (density_guess == "RESTART") THEN 
 
      CALL read_mo_set(mo,context,globenv) 
      CALL calculate_density_matrix(mo,work1,context,globenv) 
      CALL copy_blacs_to_sparse_matrix(work1,p%matrix,context,globenv) 
      CALL symmetrise_diagonal_blocks(p%matrix) 
 
    ELSE IF (density_guess == "CORE") THEN 
 
!     *** It is assumed that work1 holds the upper     *** 
!     *** triangle part of the core Hamiltonian matrix *** 
 
!MK if only gemm is used 
!     CALL symmetrise_blacs_matrix(work1,work2,context,globenv) 
      CALL orthogonalise_matrix(ortho,work1,work2,context,globenv) 
      CALL eigensolver(work1,mo,ortho,work2,context,globenv) 
      CALL calculate_density_matrix(mo,work1,context,globenv) 
      CALL copy_blacs_to_sparse_matrix(work1,p%matrix,context,globenv) 
      CALL symmetrise_diagonal_blocks(p%matrix) 
 
    ELSE IF (density_guess == "ATOMIC") THEN 
 
      sum = 0.0_wp 
 
      CALL copy_sparse_to_blacs_matrix(s%matrix,work2,context,globenv) 
 
      NULLIFY (elec_conf) 
      elec_conf => reallocate(elec_conf,0,maxlcgf) 
 
      DO ikind=1,nkind 
 
        IF (ASSOCIATED(kind_info(ikind)%all_potential)) THEN 
          DO la=0,SIZE(kind_info(ikind)%elec_conf)-1 
            elec_conf(la) = kind_info(ikind)%elec_conf(la) 
          END DO 
        ELSE IF (ASSOCIATED(kind_info(ikind)%gth_potential)) THEN 
          DO la=0,SIZE(kind_info(ikind)%gth_potential%elec_conf)-1 
            elec_conf(la) = kind_info(ikind)%gth_potential%elec_conf(la) 
          END DO 
        END IF 
 
        basis_set_a => kind_info(ikind)%orb_basis_set 
 
        DO iset=1,basis_set_a%nset 
          DO ishell=1,basis_set_a%nshell(iset) 
            la = basis_set_a%l(ishell,iset) 
            nelec = 4*la + 2 
            IF (elec_conf(la) > 0) THEN 
              IF (elec_conf(la) >= nelec) THEN 
                paa = 2.0_wp 
                elec_conf(la) = elec_conf(la) - nelec 
              ELSE 
                paa = REAL(2*elec_conf(la),wp)/REAL(nelec,wp) 
                elec_conf(la) = 0 
              END IF 
              first_sgfa = basis_set_a%first_sgf(ishell,iset) 
              last_sgfa = basis_set_a%last_sgf(ishell,iset) 
              DO isgfa=first_sgfa,last_sgfa 
                DO iatom=1,kind_info(ikind)%natom 
                  atoma = kind_info(ikind)%atom_list(iatom) 
                  isgf = atom_info(atoma)%first_sgf + isgfa - 1 
                  CALL blacs_set_element(work1,isgf,isgf,paa,context,globenv) 
                  CALL blacs_get_element(work2,isgf,isgf,saa,context,globenv) 
                  sum = sum + paa*saa 
                END DO 
              END DO 
            END IF 
          END DO 
        END DO 
 
      END DO 
 
      saa = REAL(nelectron,wp)/sum 
 
      DO ikind=1,nkind 
        basis_set_a => kind_info(ikind)%orb_basis_set 
        DO iset=1,basis_set_a%nset 
          DO ishell=1,basis_set_a%nshell(iset) 
            first_sgfa = basis_set_a%first_sgf(ishell,iset) 
            last_sgfa = basis_set_a%last_sgf(ishell,iset) 
            DO isgfa=first_sgfa,last_sgfa 
              DO iatom=1,kind_info(ikind)%natom 
                atoma = kind_info(ikind)%atom_list(iatom) 
                isgf = atom_info(atoma)%first_sgf + isgfa - 1 
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

  SUBROUTINE read_scf_parameters(start_section,end_section,globenv)

!   Purpose: Read the parameters for the SCF run.

!   History: - Creation (15.05.2001, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(INOUT) :: globenv
    CHARACTER(LEN=*), INTENT(IN)                 :: start_section,end_section

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE read_scf_parameters (MODULE qs_scf)"

!   *** Local variables ***

    CHARACTER(LEN=40) :: keyword
    LOGICAL           :: found

!   ---------------------------------------------------------------------------

!   *** Load the default values ***

    density_guess = "CORE"
    eps_eigval = 1.0E-5_wp
    eps_scf = 1.0E-5_wp
    eps_diis = 0.1_wp
    max_diis = 4
    max_scf = 50
    nrebuild_rho = 5
    nrebuild_gvg = 5
    nprow = 0
    npcol = 0
    nrow_block = 32
    ncol_block = 32
    p_mix = 0.4_wp
    work_syevx = 0.0_wp
    additional_charge = 0

    globenv%restart_file_name = "RESTART"

    CALL start_input_session(globenv%input_file_name,globenv)

!   *** Search for the requested input section ***

    CALL search(start_section,.TRUE.,found)

    IF (found) THEN

      DO

        CALL read_object(keyword,skip_lines=1)

        CALL uppercase(keyword)

        SELECT CASE (TRIM(keyword))
        CASE ("CHARGE")
          CALL read_object(additional_charge)
        CASE ("DENSITY_GUESS","SCF_GUESS","GUESS")
          CALL read_object(density_guess)
          CALL uppercase(density_guess)
        CASE ("DENSITY_MIXING","MIXING")
          CALL read_object(p_mix)
        CASE ("EPS_DIIS")
          CALL read_object(eps_diis)
        CASE ("EPS_EIGVAL")
          CALL read_object(eps_eigval)
        CASE ("EPS_SCF")
          CALL read_object(eps_scf)
        CASE ("MAX_DIIS")
          CALL read_object(max_diis)
        CASE ("MAX_SCF")
          CALL read_object(max_scf)
        CASE ("NREBUILD_RHO")
          CALL read_object(nrebuild_rho)
        CASE ("NREBUILD_GVG")
          CALL read_object(nrebuild_gvg)
        CASE ("BLOCKSIZE")
          CALL read_object(nrow_block)
          ncol_block = nrow_block
        CASE ("PROCESS_GRID")
          CALL read_object(nprow)
          CALL read_object(npcol)
        CASE ("RESTART_FILE_NAME","RESTART_FILE","RESTART")
          CALL read_object(globenv%restart_file_name)
        CASE ("WORK_SYEVX")
          CALL read_object(work_syevx)
          work_syevx = MIN(MAX(0.0_wp,work_syevx),1.0_wp)
        CASE DEFAULT
          IF (keyword == end_section) THEN
            EXIT
          ELSE
            CALL stop_program(routine,&
                              "Invalid keyword <"//TRIM(keyword)//&
                              "> found in the input section <"//&
                              TRIM(start_section)//">")
          END IF
        END SELECT

      END DO

    END IF

    CALL finish_input_session

  END SUBROUTINE read_scf_parameters

! *****************************************************************************

  SUBROUTINE write_scf_parameters(lunit)

!   Purpose: Write the parameters for the SCF run to the logical unit number
!            "lunit".

!   ***************************************************************************

    INTEGER, INTENT(IN) :: lunit

!   ---------------------------------------------------------------------------

    IF (max_scf < 1) RETURN

    WRITE (UNIT=lunit,FMT="(/,/,T2,A,/)") "SCF PARAMETERS"

    WRITE (UNIT=lunit,FMT="(T3,A,/,T3,A,I5,2(/,T3,A,ES9.2),/,T3,A,F5.2,/,T3,A,I5,/,T3,A,I5)")&
      "density guess: "//TRIM(density_guess),&
      "max_scf:      ",max_scf,&
      "eps_scf:      ",eps_scf,&
      "eps_eigval:   ",eps_eigval,&
      "p_mix:        ",p_mix,&
      "nrebuild_rho: ",nrebuild_rho,&
      "nrebuild_gvg: ",nrebuild_gvg

    WRITE(UNIT=lunit,FMT="(T3,A,I5)") "CHARGE:     ",additional_charge

    IF (max_diis > 0) THEN
      WRITE (UNIT=lunit,FMT="(T3,A,I5,/,T3,A,ES9.2)")&
        "max_diis:     ",max_diis,&
        "eps_diis:     ",eps_diis
    END IF

  END SUBROUTINE write_scf_parameters

! *****************************************************************************

END MODULE qs_scf
