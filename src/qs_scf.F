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
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE qs_scf

! *****************************************************************************

  USE kinds, ONLY: wp => dp

  USE atomic_kinds,          ONLY: kind_info,nkind
  USE atoms,                 ONLY: atom_info
  USE basis_set_types,       ONLY: gto_basis_set_type,maxlcgf
  USE coefficient_types,     ONLY: coeff_allocate,&
                                   coeff_type,&
                                   coeff_zero
  USE core_hamiltonian,      ONLY: build_core_hamiltonian_matrix,h,s
  USE diis,                  ONLY: eps_diis,max_diis,scf_diis
  USE functionals,           ONLY: init_functionals,vwn_c,vwn_x
  USE global_types,          ONLY: global_environment_type
  USE input_utilities,       ONLY: finish_input_session,&
                                   read_object,&
                                   search,&
                                   start_input_session
  USE mathlib,               ONLY: diagonalize_matrix,&
                                   power_matrix,&
                                   symmetrize_matrix
  USE matrix_types,          ONLY: allocate_matrix,&
                                   copy_matrix,&
                                   first_block_node,&
                                   get_block_node,&
                                   get_matrix_info,&
                                   next_block_node,&
                                   real_block_node_type,&
                                   real_matrix_set_type,&
                                   replicate_matrix_structure
  USE memory_utilities,      ONLY: reallocate
  USE mo_types,              ONLY: allocate_mo_set,&
                                   mo_set_type,&
                                   read_mo_set,&
                                   write_mo_set
  USE om_utilities,          ONLY: write_spherical_matrix
  USE pw_grid_types,         ONLY: HALFSPACE,pw_grid_type
  USE pw_grids,              ONLY: pw_find_cutoff,&
                                   pw_grid_construct,&
                                   pw_grid_setup
  USE pw_types,              ONLY: COMPLEXDATA3D,REALDATA3D,REALSPACE,&
                                   RECIPROCALSPACE,pw_fft_wrap,pw_type
  USE string_utilities,      ONLY: uppercase
  USE termination,           ONLY: stop_program
  USE timings,               ONLY: timeset,timestop
  USE timesl,                ONLY: cputime

  IMPLICIT NONE

  PRIVATE

  TYPE(coeff_type)           :: rho_gspace,rho_rspace,v_rspace
  TYPE(mo_set_type)          :: alpha_mo
  TYPE(real_matrix_set_type) :: ks,ortho,p,scf_work1,scf_work2
  TYPE(pw_grid_type)         :: pw_grid

  CHARACTER(LEN=10) :: density_guess = "CORE"
  REAL(wp)          :: ec = 0.0_wp,&
                       ecore = 0.0_wp,&
                       ecore_overlap = 0.0_wp,&
                       ecore_self = 0.0_wp,&
                       ehartree = 0.0_wp,&
                       eps_scf = 1.0E-5_wp,&
                       etotal = 0.0_wp,&
                       ex = 0.0_wp,&
                       p_mix = 0.4_wp,&
                       total_rho_core_rspace = 0.0_wp,&
                       total_rho_elec_rspace = 0.0_wp,&
                       total_rho_gspace = 0.0_wp,&
                       total_rho_rspace = 0.0_wp
  INTEGER           :: max_scf = 50,&
                       nelectron = 0

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
    INTEGER  :: handle,iscf,output_unit
    LOGICAL  :: diis_step

    REAL(wp), DIMENSION(:,:), POINTER :: scf_work1_matrix

!   ---------------------------------------------------------------------------

    CALL timeset("scf","I","",handle)

    output_unit = globenv%scr

!   *** Quick return, if no SCF iteration is requested ***

    IF (max_scf < 1) RETURN

    CALL init_functionals()

    CALL init_grid(globenv)

    IF (globenv%print%scf) THEN
      WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
        "SCF WAVEFUNCTION OPTIMIZATION"
    END IF

    CALL init_scf_run(globenv,alpha_mo)

    CALL get_block_node(matrix=scf_work1%matrix,&
                        block_row=1,&
                        block_col=1,&
                        block=scf_work1_matrix)

    etotal = 0.0_wp
    iscf = 0
    diis_step = .FALSE.

    IF (globenv%print%scf) THEN
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

      IF (globenv%print%kohn_sham_matrix.AND.globenv%print%each_scf_step) THEN
        CALL write_spherical_matrix(ks%matrix,4,6,globenv)
      END IF

      scf_work1_matrix(:,:) = 0.0_wp
      CALL copy_matrix(ks%matrix,scf_work1%matrix)
      CALL symmetrize_matrix(scf_work1_matrix,"upper_to_lower")

      IF (iscf > 1) THEN
        CALL scf_diis(globenv,alpha_mo,scf_work1,scf_work2,delta,diis_error,&
                      diis_step)
      END IF

      CALL orthogonalize_matrix(ortho,scf_work1,scf_work2)

      CALL eigensolver(scf_work1,alpha_mo,ortho)

      IF (globenv%print%each_scf_step) THEN
        CALL write_mo_set(alpha_mo,4,6,globenv)
      END IF

      CALL calculate_density_matrix(alpha_mo,scf_work1)

      CALL copy_matrix(scf_work1%matrix,ks%matrix)

      t2 = cputime()

      IF (diis_step) THEN
        CALL density_mixing(ks,p,1.0_wp,delta)
        IF (globenv%print%scf) THEN
          WRITE (UNIT=output_unit,&
                 FMT="(T2,I5,2X,A,T15,E10.2,T30,F8.2,T40,2F20.10)")&
            iscf,"DIIS",diis_error,t2 - t1,delta,etotal
        END IF
      ELSE
        CALL density_mixing(ks,p,p_mix,delta)
        IF (globenv%print%scf) THEN
          WRITE (UNIT=output_unit,&
                 FMT="(T2,I5,2X,A,T15,F6.2,T30,F8.2,T40,2F20.10)")&
            iscf,"Mixing",p_mix,t2 - t1,delta,etotal
        END IF
      END IF

      IF (delta < eps_scf) THEN
        IF (globenv%print%scf) THEN
          WRITE(UNIT=output_unit,FMT="(/,T3,A,/)")&
            "*** SCF run converged ***"
        END IF
        IF (.NOT.diis_step) CALL copy_matrix(ks%matrix,p%matrix)
        EXIT
      ELSE IF (iscf == max_scf) THEN
        IF (globenv%print%scf) THEN
          WRITE(UNIT=output_unit,FMT="(/,T3,A,/)")&
            "*** SCF run NOT converged ***"
        END IF
        EXIT
      END IF

    END DO

    IF (globenv%print%scf) THEN
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

    CALL write_mo_set(alpha_mo,4,6,globenv)

    IF (globenv%print%density_matrix) THEN
      CALL write_spherical_matrix(p%matrix,4,6,globenv)
    END IF

    IF (globenv%print%kohn_sham_matrix) THEN
      CALL build_kohn_sham_matrix(globenv)
      CALL write_spherical_matrix(ks%matrix,4,6,globenv)
    END IF

!   *** Write restart file ***

    CALL write_mo_set(alpha_mo,globenv)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE scf

! *****************************************************************************

  SUBROUTINE build_kohn_sham_matrix(globenv)

    USE collocate_density,   ONLY: calculate_rho_core,&
                                   calculate_rho_elec,&
                                   calculate_total_rho
    USE core_energies,       ONLY: calculate_ecore
    USE integrate_potential, ONLY: integrate_v_rspace
    USE print_keys,          ONLY: DEBUG

    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    INTEGER :: handle,nproduct,output_unit

!   ---------------------------------------------------------------------------

    CALL timeset("build_kohn_sham_matrix","I","",handle)

    output_unit = globenv%scr

    CALL calculate_ecore(h,p,ecore)

    CALL coeff_zero(rho_rspace)

    CALL calculate_rho_elec(p,rho_rspace,total_rho_elec_rspace,nproduct)

    IF (globenv%print%total_densities) THEN
      WRITE (UNIT=output_unit,FMT="(/,T3,A,I10)")&
        "Number of collocated products (r-space):",nproduct
    END IF

    CALL coeff_zero(v_rspace)

    CALL calculate_xc_potential(rho_rspace%pw,v_rspace%pw)

    CALL calculate_rho_core(rho_rspace,total_rho_rspace)

    total_rho_core_rspace = total_rho_rspace - total_rho_elec_rspace

    CALL pw_fft_wrap(rho_rspace%pw,rho_gspace%pw,&
                     debug=(globenv%print%level == DEBUG))

    total_rho_gspace = calculate_total_rho(rho_gspace)

    IF (globenv%print%total_densities) THEN
      WRITE (UNIT=output_unit,FMT="(/,(T3,A,T40,2F20.10))")&
        "Total electronic density (r-space): ",&
        total_rho_elec_rspace,total_rho_elec_rspace + REAL(nelectron,wp),&
        "Total core charge density (r-space):",&
        total_rho_core_rspace,total_rho_core_rspace - REAL(nelectron,wp)
      WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
        "Total charge density (r-space):     ",total_rho_rspace,&
        "Total charge density (g-space):     ",total_rho_gspace
    END IF

    CALL calculate_hartree_potential(rho_gspace)

    etotal = ecore_overlap + ecore_self + ecore + ehartree + ex + ec

    IF (globenv%print%scf_energies) THEN
      WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
        "Overlap energy of the core charge distribution:",ecore_overlap,&
        "Self energy of the core charge distribution:   ",ecore_self,&
        "Core Hamiltonian energy:                       ",ecore,&
        "Hartree energy:                                ",ehartree,&
        "Exchange-correlation energy:                   ",ex + ec
    END IF

    CALL pw_fft_wrap(rho_gspace%pw,rho_rspace%pw,&
                     debug=(globenv%print%level == DEBUG))

    v_rspace%pw%cr3d(:,:,:) = v_rspace%pw%pw_grid%dvol*&
                              (v_rspace%pw%cr3d(:,:,:) +&
                               rho_rspace%pw%cr3d(:,:,:))

    CALL copy_matrix(h%matrix,ks%matrix)

    CALL integrate_v_rspace(v_rspace,ks,nproduct)

    IF (globenv%print%total_densities) THEN
      WRITE (UNIT=output_unit,FMT="(/,T3,A,I10)")&
        "Number of integrated products (r-space):",nproduct
    END IF

    CALL timestop(0.0_wp,handle)

 END SUBROUTINE build_kohn_sham_matrix

! *****************************************************************************

  SUBROUTINE init_grid(globenv)

    USE cell_parameters, ONLY: abc,cell,cutoff

    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    INTEGER :: handle,output_unit

!   ---------------------------------------------------------------------------

    CALL timeset("init_grid","I","",handle)

    output_unit = globenv%scr

    CALL pw_grid_construct(pw_grid)

    pw_grid%grid_span = HALFSPACE

    IF (globenv%print%pw_grid_information) THEN
      CALL pw_grid_setup(cell,pw_grid,cutoff,info=output_unit,&
                         orthorhombic=.TRUE.)
    ELSE
      CALL pw_grid_setup(cell,pw_grid,cutoff,orthorhombic=.TRUE.)
    END IF

    CALL coeff_allocate(rho_rspace,pw_grid,REALDATA3D)
    rho_rspace%pw%in_space = REALSPACE

    CALL coeff_allocate(rho_gspace,pw_grid,COMPLEXDATA3D)
    rho_gspace%pw%in_space = RECIPROCALSPACE

    CALL coeff_allocate(v_rspace,pw_grid,REALDATA3D)
    v_rspace%pw%in_space = REALSPACE

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE init_grid

! *****************************************************************************

  SUBROUTINE init_scf_run(globenv,mo)

!   Purpose: Initialize a SCF run.

!   History: - Creation (30.04.2001, Matthias Krack)

!   ***************************************************************************

    USE core_energies, ONLY: calculate_ecore_overlap,&
                             calculate_ecore_self

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(mo_set_type), INTENT(OUT)            :: mo

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE init_scf_run (MODULE qs_scf)"

!   *** Local variables ***

    CHARACTER(LEN=40) :: symmetry
    INTEGER           :: handle,homo,ikind,nao,ndep,nmo,output_unit

!   ---------------------------------------------------------------------------

    CALL timeset("init_scf_run","I","",handle)

    output_unit = globenv%scr

    nelectron = 0

    DO ikind=1,nkind
      nelectron = nelectron + kind_info(ikind)%natom*kind_info(ikind)%zeff
    END DO

    IF (MODULO(nelectron,2) == 0) THEN
      homo = nelectron/2
    ELSE
      CALL stop_program(routine,"Odd number of electrons")
    END IF

!   *** Get the dimension of the full SCF matrices, ***
!   *** i.e. the total number of atomic orbitals    ***

    CALL get_matrix_info(matrix=h%matrix,&
                         nrow=nao,&
                         matrix_symmetry=symmetry)

    nmo = nao
!MK nmo = homo

    CALL allocate_mo_set(mo,nao,nmo)

    mo%homo = homo
    mo%occupation_numbers(1:homo) = 2.0_wp

!   *** Allocate the full SCF matrices ***

    NULLIFY (ortho%matrix)
    CALL allocate_matrix(matrix=ortho%matrix,&
                         nrow=nao,&
                         ncol=nao,&
                         matrix_name="ORTHOGONALIZATION MATRIX",&
                         matrix_symmetry=symmetry)

    NULLIFY (scf_work1%matrix)
    CALL allocate_matrix(matrix=scf_work1%matrix,&
                         nrow=nao,&
                         ncol=nao,&
                         matrix_name="SCF WORK MATRIX 1",&
                         matrix_symmetry="none")

    NULLIFY (scf_work2%matrix)
    CALL allocate_matrix(matrix=scf_work2%matrix,&
                         nrow=nao,&
                         ncol=nao,&
                         matrix_name="SCF WORK MATRIX 2",&
                         matrix_symmetry="none")

    NULLIFY (ks%matrix)
    CALL replicate_matrix_structure(h%matrix,ks%matrix,"KOHN-SHAM MATRIX")
    NULLIFY (p%matrix)
    CALL replicate_matrix_structure(h%matrix,p%matrix,"DENSITY MATRIX")

    CALL calculate_ecore_self(ecore_self)
    CALL calculate_ecore_overlap(ecore_overlap)

    IF (globenv%print%scf_energies) THEN
      WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
        "Self energy of the core charge distribution:   ",ecore_self,&
        "Overlap energy of the core charge distribution:",ecore_overlap
    END IF

    CALL calculate_ortho_matrix(ortho,scf_work1,ndep)

    IF (globenv%print%scf) THEN
      WRITE (UNIT=output_unit,FMT="(/,(T3,A,I10))")&
        "Number of electrons:                    ",nelectron,&
        "Number of occupied orbitals:            ",homo,&
        "Number of orbital functions:            ",nao,&
        "Number of independent orbital functions:",nao - ndep
    END IF

    IF (globenv%print%ortho_matrix) THEN
      CALL write_spherical_matrix(ortho%matrix,4,6,globenv)
    END IF

    CALL calculate_first_density_matrix(ortho,mo,p,scf_work1,scf_work2,globenv)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE init_scf_run

! *****************************************************************************

  SUBROUTINE calculate_ortho_matrix(ortho,work,ndep)

    TYPE(real_matrix_set_type), INTENT(OUT)   :: ortho
    TYPE(real_matrix_set_type), INTENT(INOUT) :: work
    INTEGER, INTENT(OUT)                      :: ndep

!   *** Local variables ***

    INTEGER :: handle,nao

    REAL(wp), DIMENSION(:,:), POINTER :: ortho_matrix,work_matrix

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_ortho_matrix","I","",handle)

    CALL get_matrix_info(matrix=ortho%matrix,nrow=nao)

    CALL get_block_node(matrix=ortho%matrix,&
                        block_row=1,&
                        block_col=1,&
                        block=ortho_matrix)

    CALL get_block_node(matrix=work%matrix,&
                        block_row=1,&
                        block_col=1,&
                        block=work_matrix)

!   *** Calculate the orthogonalization matrix (S**(-1/2)) ***

    work_matrix(:,:) = 0.0_wp

    CALL copy_matrix(s%matrix,work%matrix)

    CALL symmetrize_matrix(work_matrix,"upper_to_lower")

    CALL power_matrix(work_matrix,ortho_matrix,-0.5_wp,1.0E-5_wp,ndep,.FALSE.)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_ortho_matrix

! *****************************************************************************

  SUBROUTINE orthogonalize_matrix(ortho,ks,work)

    TYPE(real_matrix_set_type), INTENT(IN)    :: ortho
    TYPE(real_matrix_set_type), INTENT(INOUT) :: ks,work

!   *** Local variables ***

    INTEGER :: handle,nao

    REAL(wp), DIMENSION(:,:), POINTER :: ks_matrix,ortho_matrix,work_matrix

!   ---------------------------------------------------------------------------

    CALL timeset("orthogonalize_matrix","I","",handle)

    CALL get_matrix_info(matrix=ks%matrix,nrow=nao)

    CALL get_block_node(matrix=ks%matrix,&
                        block_row=1,&
                        block_col=1,&
                        block=ks_matrix)

    CALL get_block_node(matrix=ortho%matrix,&
                        block_row=1,&
                        block_col=1,&
                        block=ortho_matrix)

    CALL get_block_node(matrix=work%matrix,&
                        block_row=1,&
                        block_col=1,&
                        block=work_matrix)

    CALL dgemm("T","N",nao,nao,nao,1.0_wp,ks_matrix,nao,ortho_matrix,nao,&
               0.0_wp,work_matrix,nao)

    CALL dgemm("N","N",nao,nao,nao,1.0_wp,ortho_matrix,nao,work_matrix,nao,&
               0.0_wp,ks_matrix,nao)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE orthogonalize_matrix

! *****************************************************************************

  SUBROUTINE eigensolver(ks,mo,ortho)

    TYPE(real_matrix_set_type), INTENT(IN)    :: ortho
    TYPE(real_matrix_set_type), INTENT(INOUT) :: ks
    TYPE(mo_set_type), INTENT(INOUT)          :: mo

!   *** Local variables ***

    INTEGER :: handle,nao,nmo

    REAL(wp), DIMENSION(:,:), POINTER :: ks_matrix,mo_eigenvectors,ortho_matrix

!   ---------------------------------------------------------------------------

    CALL timeset("eigensolver","I","",handle)

    CALL get_matrix_info(matrix=mo%eigenvectors%matrix,nrow=nao,ncol=nmo)

    CALL get_block_node(matrix=ks%matrix,&
                        block_row=1,&
                        block_col=1,&
                        block=ks_matrix)

    CALL get_block_node(matrix=mo%eigenvectors%matrix,&
                        block_row=1,&
                        block_col=1,&
                        block=mo_eigenvectors)

    CALL get_block_node(matrix=ortho%matrix,&
                        block_row=1,&
                        block_col=1,&
                        block=ortho_matrix)

!   CALL diagonalize_matrix(ks_matrix,mo_eigenvectors,mo%eigenvalues)
!   CALL dgemm("N","N",nao,nmo,nao,1.0_wp,ortho_matrix,nao,mo_eigenvectors,&
!              nao,0.0_wp,ks_matrix,nao)
!   mo_eigenvectors(:,:) = ks_matrix(:,:)

    CALL diagonalize_matrix(ks_matrix,mo%eigenvalues,.FALSE.)

    CALL dgemm("N","N",nao,nmo,nao,1.0_wp,ortho_matrix,nao,ks_matrix,&
               nao,0.0_wp,mo_eigenvectors,nao)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE eigensolver

! *****************************************************************************

  SUBROUTINE calculate_density_matrix(mo,density)

    TYPE(real_matrix_set_type), INTENT(OUT) :: density
    TYPE(mo_set_type), INTENT(IN)           :: mo

!   *** Local variables ***

    INTEGER :: handle,nao

    REAL(wp), DIMENSION(:,:), POINTER :: mo_eigenvectors,density_matrix

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_density_matrix","I","",handle)

    CALL get_matrix_info(matrix=mo%eigenvectors%matrix,nrow=nao)

    CALL get_block_node(matrix=mo%eigenvectors%matrix,&
                        block_row=1,&
                        block_col=1,&
                        block=mo_eigenvectors)

    CALL get_block_node(matrix=density%matrix,&
                        block_row=1,&
                        block_col=1,&
                        block=density_matrix)

    CALL dsyrk("U","N",nao,mo%homo,2.0_wp,mo_eigenvectors,nao,0.0_wp,&
               density_matrix,nao)

    CALL symmetrize_matrix(density_matrix,"upper_to_lower")

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_density_matrix

! *****************************************************************************

  SUBROUTINE density_mixing(new_density,old_density,p_mix,delta)

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

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE density_mixing

! *****************************************************************************

  SUBROUTINE calculate_hartree_potential(rho)

    USE mathconstants, ONLY: fourpi,twopi

    TYPE(coeff_type), TARGET, INTENT(INOUT) :: rho

!   *** Local variables ***

    COMPLEX(wp) :: vhartree
    REAL(wp)    :: g2,gx,gy,gz
    INTEGER     :: handle,i,ig,j,jg,k,kg

    REAL(wp), DIMENSION(3) :: dg

    REAL(wp), DIMENSION(:), POINTER :: dr
    INTEGER, DIMENSION(:), POINTER  :: lb_grid,ng,ub_grid

    COMPLEX(wp), DIMENSION(:,:,:), POINTER :: grid

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_hartree_potential","I","",handle)

    dr => rho%pw%pw_grid%dr(:)
    ng => rho%pw%pw_grid%npts(:)
    lb_grid => rho%pw%pw_grid%bounds(1,:)
    ub_grid => rho%pw%pw_grid%bounds(2,:)

    dg(:) = twopi/(ng(:)*dr(:))

    grid => rho%pw%cc3d(:,:,:)

    ehartree = 0.0_wp

    DO kg=lb_grid(3),ub_grid(3)
      gz = REAL(kg,wp)*dg(3)
      IF (kg < 0) THEN
        k = kg + 1 + ng(3)
      ELSE
        k = kg + 1
      END IF
      DO jg=lb_grid(2),ub_grid(2)
        gy = REAL(jg,wp)*dg(2)
        IF (jg < 0) THEN
          j = jg + 1 + ng(2)
        ELSE
          j = jg + 1
        END IF
        DO ig=lb_grid(1),ub_grid(1)
          gx = REAL(ig,wp)*dg(1)
          IF (ig < 0) THEN
            i = ig + 1 + ng(1)
          ELSE
            i = ig + 1
          END IF
          g2 = gx*gx + gy*gy + gz*gz
          IF (g2 > 1.0E-12_wp) THEN
!!!write (*,"(3I6,2F15.6)") i,j,k,grid(i,j,k)
            vhartree = fourpi*grid(i,j,k)/g2
            ehartree = ehartree + REAL(grid(i,j,k))*REAL(vhartree) +&
                                  AIMAG(grid(i,j,k))*AIMAG(vhartree)
            grid(i,j,k) = vhartree
          END IF
        END DO
      END DO
    END DO

    ehartree = 0.5_wp*rho%pw%pw_grid%vol*ehartree

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_hartree_potential

! *****************************************************************************

  SUBROUTINE calculate_xc_potential(rho,vxc)

    TYPE(pw_type), INTENT(IN)    :: rho
    TYPE(pw_type), INTENT(INOUT) :: vxc

!   *** Local variables ***

    INTEGER :: handle

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_xc_potential","I","",handle)

    ex = 0.0_wp
    ec = 0.0_wp

    CALL vwn_x(rho%cr3d,ex,vxc%cr3d)
    CALL vwn_c(rho%cr3d,ec,vxc%cr3d)

    ex = rho%pw_grid%dvol*ex
    ec = rho%pw_grid%dvol*ec

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_xc_potential

! *****************************************************************************

  SUBROUTINE calculate_first_density_matrix(ortho,mo,p,work1,work2,globenv)

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(real_matrix_set_type), INTENT(IN)    :: ortho
    TYPE(real_matrix_set_type), INTENT(INOUT) :: p,work1,work2
    TYPE(mo_set_type), INTENT(INOUT)          :: mo

!   *** Local variables ***

    TYPE(gto_basis_set_type), POINTER :: basis_set_a

    REAL(wp) :: paa
    INTEGER  :: atoma,first_sgfa,handle,iatom,ikind,iset,isgf,isgfa,ishell,la,&
                last_sgfa,nelec,restart_unit

    INTEGER, DIMENSION(:), POINTER :: elec_conf

    REAL(wp), DIMENSION(:,:), POINTER :: work1_matrix,work2_matrix

!   ---------------------------------------------------------------------------

    CALL timeset("calculate_first_density_matrix","I","",handle)

!   *** Calculate an initial density guess (core Hamiltonian guess) ***

    CALL get_block_node(matrix=work1%matrix,&
                        block_row=1,&
                        block_col=1,&
                        block=work1_matrix)

      work1_matrix(:,:) = 0.0_wp

    IF (density_guess == "RESTART") THEN

      CALL read_mo_set(mo,globenv)
      CALL calculate_density_matrix(mo,work1)

    ELSE IF (density_guess == "CORE") THEN

      CALL copy_matrix(h%matrix,work1%matrix)
      CALL symmetrize_matrix(work1_matrix,"upper_to_lower")
      CALL orthogonalize_matrix(ortho,work1,work2)
      CALL eigensolver(work1,mo,ortho)
      CALL calculate_density_matrix(mo,work1)

    ELSE IF (density_guess == "ATOMIC") THEN

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
                  work1_matrix(isgf,isgf) = paa
                END DO
              END DO
            END IF
          END DO
        END DO

      END DO

    ELSE

      STOP "qs_scf: wrong keyword for guess"

    END IF

    CALL copy_matrix(work1%matrix,p%matrix)

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
    eps_scf = 1.0E-5_wp
    eps_diis = 0.1_wp
    max_diis = 4
    max_scf = 50
    p_mix = 0.4_wp
    globenv%restart_file_name = "RESTART"

    CALL start_input_session(globenv%input_file_name,globenv)

!   *** Search for the requested input section ***

    CALL search(start_section,.TRUE.,found)

    IF (found) THEN

      DO

        CALL read_object(keyword,skip_lines=1)

        CALL uppercase(keyword)

        SELECT CASE (TRIM(keyword))
        CASE ("DENSITY_GUESS","SCF_GUESS","GUESS")
          CALL read_object(density_guess)
          CALL uppercase(density_guess)
        CASE ("DENSITY_MIXING","MIXING")
          CALL read_object(p_mix)
        CASE ("EPS_DIIS")
          CALL read_object(eps_diis)
        CASE ("EPS_SCF")
          CALL read_object(eps_scf)
        CASE ("MAX_DIIS")
          CALL read_object(max_diis)
        CASE ("MAX_SCF")
          CALL read_object(max_scf)
        CASE ("RESTART_FILE_NAME","RESTART_FILE","RESTART")
          CALL read_object(globenv%restart_file_name)
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

    WRITE (lunit,"(/,/,T2,A,/)") "SCF PARAMETERS"

    WRITE (lunit,"(T3,A,/,T3,A,I5,/,T3,A,ES9.2,/,T3,A,F5.2)")&
      "density guess: "//TRIM(density_guess),&
      "max_scf:      ",max_scf,&
      "eps_scf:      ",eps_scf,&
      "p_mix:        ",p_mix

    IF (max_diis > 0) THEN
      WRITE (lunit,"(T3,A,I5,/,T3,A,ES9.2)")&
        "max_diis:     ",max_diis,&
        "eps_diis:     ",eps_diis
    END IF

  END SUBROUTINE write_scf_parameters

! *****************************************************************************

END MODULE qs_scf
