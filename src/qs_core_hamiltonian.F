!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_core_hamiltonian [1.1] *
!!
!!   FUNCTION
!!     Calculation of the core Hamiltonian integral matrix <a|H|b> over
!!     Cartesian Gaussian-type functions.
!!
!!     <a|H|b> = <a|T|b> + <a|V|b>
!!
!!
!!     Kinetic energy:
!!
!!     <a|T|b> = <a|-nabla**2/2|b>
!!               \_______________/
!!                       |
!!                    kinetic
!!
!!     Nuclear potential energy:
!!
!!     a) Allelectron calculation:
!!
!!                         erfc(r)
!!        <a|V|b> = -Z*<a|---------|b>
!!                            r
!!
!!                         1 - erf(r)
!!                = -Z*<a|------------|b>
!!                             r
!!
!!                          1           erf(r)
!!                = -Z*(<a|---|b> - <a|--------|b>)
!!                          r             r
!!
!!                          1
!!                = -Z*(<a|---|b> - N*<ab||c>)
!!                          r
!!
!!                     -Z
!!                = <a|---|b> + Z*N*<ab||c>
!!                      r
!!                  \_______/       \_____/
!!                      |              |
!!                   nuclear        coulomb
!!
!!     b) Pseudopotential calculation (Goedecker, Teter and Hutter; GTH):
!!
!!        <a|V|b> = <a|(V(local) + V(non-local))|b>
!!
!!                = <a|(V(local)|b> + <a|V(non-local))|b>
!!
!!        <a|V(local)|b> = <a|-Z(eff)*erf(SQRT(2)*alpha*r)/r +
!!                            (C1 + C2*(alpha*r)**2 + C3*(alpha*r)**4 +
!!                             C4*(alpha*r)**6)*exp(-(alpha*r)**2/2))|b>
!!
!!        <a|V(non-local)|b> = <a|p(l,i)>*h(i,j)*<p(l,j)|b>
!!
!!   AUTHOR
!!     Matthias Krack (14.09.2000,21.03.02)
!!
!!   LITERATURE
!!     S. Goedecker, M. Teter and J. Hutter, Phys. Rev. B 54, 1703 (1996)
!!     C. Hartwigsen, S. Goedecker and J. Hutter, Phys. Rev. B 58, 3641 (1998)
!!     M. Krack and M. Parrinello, Phys. Chem. Chem. Phys. 2, 2105 (2000)
!!     S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!!
!!   MODIFICATION HISTORY
!!     Joost VandeVondele (April 2003) : added LSD forces
!!     - Non-redundant calculation of the non-local part of the GTH PP
!!       (22.05.2003,MK)
!!     - New parallelization scheme (27.06.2003,MK)
!!
!!   SOURCE
!******************************************************************************

MODULE qs_core_hamiltonian

  USE ai_kinetic,                      ONLY: kinetic
  USE ai_overlap_new,                  ONLY: overlap
  USE ai_overlap_ppl,                  ONLY: overlap_ppl
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE external_potential_types,        ONLY: get_potential,&
                                             gth_potential_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: int_size,&
                                             wp => dp,&
                                             wp_size => dp_size
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: indco,&
                                             init_orbital_pointers,&
                                             nco,&
                                             ncoset
  USE orbital_symbols,                 ONLY: cgf_symbol
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_matrix_pools,                 ONLY: qs_matrix_pools_type
  USE qs_neighbor_list_types,          ONLY: find_neighbor_list,&
                                             first_list,&
                                             first_node,&
                                             get_neighbor_list,&
                                             get_neighbor_list_set,&
                                             get_neighbor_node,&
                                             neighbor_list_set_p_type,&
                                             neighbor_list_set_type,&
                                             neighbor_list_type,&
                                             neighbor_node_type,&
                                             next
  USE qs_overlap,                      ONLY: write_matrix_distribution,&
                                             write_sparse_matrix
  USE qs_parser,                       ONLY: close_file,&
                                             open_file
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE scf_control_types,               ONLY: scf_control_type
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             add_matrices,&
                                             allocate_matrix,&
                                             allocate_matrix_set,&
                                             cp_sm_scale_and_add,&
                                             deallocate_matrix,&
                                             deallocate_matrix_row,&
                                             deallocate_matrix_set,&
                                             first_block_node,&
                                             get_block_node,&
                                             get_matrix_info,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_p_type,&
                                             real_matrix_type,&
                                             replicate_matrix,&
                                             replicate_matrix_structure,&
                                             set_matrix,&
                                             write_blocks_maxabsval
  USE string_utilities,                ONLY: compress,&
                                             uppercase
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_core_hamiltonian"

! *** Public subroutines ***

  PUBLIC :: build_core_hamiltonian_matrix

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE build_core_hamiltonian_matrix(qs_env,globenv,calculate_forces)

!   Purpose: Calculation of the overlap matrix over Cartesian Gaussian
!            functions.

!   History: - Creation (11.03.2002,MK)
!            - Non-redundant calculation of the non-local part of the GTH PP
!              (22.05.2003,MK)
!            - New parallelization scheme (27.06.2003,MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(qs_environment_type), POINTER        :: qs_env
    LOGICAL, INTENT(IN)                       :: calculate_forces

!   *** Local types ***

    TYPE block_p_type
      REAL(wp), DIMENSION(:,:), POINTER :: block
    END TYPE block_p_type

    TYPE vppl_type
      REAL(wp)                          :: alpha_ppl,ppl_radius
      REAL(wp), DIMENSION(:), POINTER   :: cexp_ppl,r2
      INTEGER, DIMENSION(:), POINTER    :: neighbor
      REAL(wp), DIMENSION(:,:), POINTER :: r
    END TYPE vppl_type

    TYPE vppnl_type
      INTEGER                                   :: lppnl,nppnl
      REAL(wp)                                  :: ppnl_radius
      TYPE(block_p_type), DIMENSION(:), POINTER :: sac
      REAL(wp), DIMENSION(:), POINTER           :: alpha_ppnl,r2
      INTEGER, DIMENSION(:), POINTER            :: nprj_ppnl
      INTEGER, DIMENSION(:,:), POINTER          :: neighbor
      REAL(wp), DIMENSION(:,:), POINTER         :: cprj,r,vprj_ppnl
    END TYPE vppnl_type

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "core_hamiltonian"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(gth_potential_type), POINTER :: gth_potential
    TYPE(gto_basis_set_type), POINTER :: orb_basis_set
    TYPE(neighbor_list_type), POINTER :: sab_orb_neighbor_list,&
                                         sac_ppl_neighbor_list,&
                                         sac_ppnl_neighbor_list,&
                                         sbc_ppnl_neighbor_list
    TYPE(neighbor_node_type), POINTER :: sab_orb_neighbor_node,&
                                         sac_ppl_neighbor_node,&
                                         sac_ppnl_neighbor_node,&
                                         sbc_ppnl_neighbor_node
    TYPE(qs_rho_type), POINTER        :: rho
    TYPE(scf_control_type), POINTER   :: scf_control

    CHARACTER(LEN=6)   :: extension
    CHARACTER(LEN=12)  :: cgfsym
    CHARACTER(LEN=80)  :: name
    CHARACTER(LEN=200) :: file_name
    REAL(wp)           :: dab,dac,dbc,f,f0,rab2,rac2,rbc2
    INTEGER            :: atom_a,atom_b,atom_c,first_col,first_prjc,handle,i,&
                          iab,iac,ibc,iatom,ico,icol,ikind,ilist,inode,ipgf,&
                          irow,iset,ishell,ispin,istat,j,jatom,jco,&
                          jkind,jpgf,jset,&
                          jshell,katom,kkind,kneighbor,knode,l,last_jatom,&
                          lc_max,lc_min,ldsab,ldwork,maxblock,maxco,maxdco,&
                          maxder,maxl,maxlgto,maxlppl,maxlppnl,maxppnl,maxsgf,&
                          natom,ncoa,ncoa_sum,ncob,ncoc,ncol,&
                          nder,nkind,nlist,nneighbor,nnode,nnode_c,nprjc,nrow,&
                          nseta,nsetb,nsgf,prjc,sgfa,sgfb,unit_number
    LOGICAL            :: build_kinetic_energy_matrix,gth_potential_present,&
                          new_atom_b,ppl_present,ppnl_present,&
                          return_s_derivatives,return_t_derivatives

    REAL(wp), DIMENSION(1) :: rprjc,zetc
    REAL(wp), DIMENSION(3) :: force_a,force_b,rab,rac,rbc
    INTEGER, DIMENSION(1)  :: first_row,last_row
    INTEGER, DIMENSION(3)  :: cell_c

    TYPE(block_p_type), DIMENSION(:), ALLOCATABLE :: sint,tint
    TYPE(vppl_type), DIMENSION(:), ALLOCATABLE    :: vppl
    TYPE(vppnl_type), DIMENSION(:), ALLOCATABLE   :: vppnl

    REAL(wp), DIMENSION(:,:), ALLOCATABLE   :: hab,pab,sab,vdab,vadb,work
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: sdab
    INTEGER, DIMENSION(:), ALLOCATABLE      :: atom_of_kind,&
                                               first_prj,&
                                               first_sgf,&
                                               last_prj,&
                                               last_sgf,&
                                               kind_of

    TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb,sac_ppl,&
                                                             sac_ppnl,sbc_ppnl
    TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set
    TYPE(qs_force_type), DIMENSION(:), POINTER            :: force
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER       :: h,p,s,t,v,w

    REAL(wp), DIMENSION(:), POINTER   :: set_radius_a,set_radius_b
    REAL(wp), DIMENSION(:,:), POINTER :: h_block,p_block,rpgfa,rpgfb,sphi_a,&
                                         sphi_b,w_block,zeta,zetb
    INTEGER, DIMENSION(:), POINTER    :: la_max,la_min,lb_max,lb_min,npgfa,&
                                         npgfb,nsgfa,nsgfb
    INTEGER, DIMENSION(:,:), POINTER  :: first_sgfa,first_sgfb
    TYPE(distribution_2d_type), POINTER :: distribution_2d

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    NULLIFY (atomic_kind_set)
    NULLIFY (force)
    NULLIFY (h)
    NULLIFY (p)
    NULLIFY (particle_set)
    NULLIFY (rho)
    NULLIFY (s)
    NULLIFY (sab_orb)
    NULLIFY (sac_ppl)
    NULLIFY (sac_ppnl)
    NULLIFY (sbc_ppnl)
    NULLIFY (scf_control)
    NULLIFY (t)
    NULLIFY (w)
    NULLIFY (distribution_2d)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    h=h,&
                    kinetic=t,&
                    particle_set=particle_set,&
                    s=s,&
                    sab_orb=sab_orb,&
                    sac_ppl=sac_ppl,&
                    sac_ppnl=sac_ppnl,&
                    sbc_ppnl=sbc_ppnl,&
                    scf_control=scf_control,&
                    distribution_2d=distribution_2d)

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    IF (calculate_forces) THEN

      CALL timeset(routine_name//" (forces)","I","",handle)

      nder = 1

      ALLOCATE (atom_of_kind(natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "atom_of_kind",natom*int_size)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               atom_of_kind=atom_of_kind)

      CALL get_qs_env(qs_env=qs_env,&
                      force=force,&
                      rho=rho,&
                      w=w)

      p => rho%rho_ao

!     *** If LSD, then combine alpha density and beta density to
!     *** total density: alpha <- alpha + beta   and
!     *** spin density:   beta <- alpha - beta

      IF (SIZE(p) == 2) THEN
        CALL cp_sm_scale_and_add(p(1)%matrix, 1.0_wp,p(2)%matrix,1.0_wp)
        CALL cp_sm_scale_and_add(p(2)%matrix,-2.0_wp,p(1)%matrix,1.0_wp)
        CALL cp_sm_scale_and_add(w(1)%matrix, 1.0_wp,w(2)%matrix,1.0_wp)
        CALL cp_sm_scale_and_add(w(2)%matrix,-2.0_wp,w(1)%matrix,1.0_wp)
      END IF

    ELSE

      CALL timeset(routine_name,"I","",handle)

      IF (globenv%print%derivatives.OR.globenv%print%optical_conductivity) THEN
        nder = 1
      ELSE
        nder = 0
      END IF

    END IF

    IF ((globenv%print%overlap_matrix.AND.globenv%print%derivatives).OR.&
        globenv%print%optical_conductivity) THEN
      return_s_derivatives = .TRUE.
    ELSE
      return_s_derivatives = .FALSE.
    END IF

    maxder = ncoset(nder)

    CALL allocate_matrix_set(s,maxder)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,nsgf=nsgf)

    ALLOCATE (first_sgf(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "first_sgf",natom*int_size)
    ALLOCATE (last_sgf(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "last_sgf",natom*int_size)

    CALL get_particle_set(particle_set=particle_set,&
                          first_sgf=first_sgf,&
                          last_sgf=last_sgf)

    ! assert here that the S-matrix is of the type distribution_2d if this one is non-null

    CALL allocate_matrix(matrix=s(1)%matrix,&
                         nblock_row=natom,&
                         nblock_col=natom,&
                         nrow=nsgf,&
                         ncol=nsgf,&
                         first_row=first_sgf,&
                         last_row=last_sgf,&
                         first_col=first_sgf,&
                         last_col=last_sgf,&
                         matrix_name="OVERLAP MATRIX",&
                         matrix_symmetry="symmetric",&
                         distribution_2d=distribution_2d)

    DEALLOCATE (first_sgf,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "first_sgf")
    DEALLOCATE (last_sgf,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "last_sgf")

    DO i=2,maxder
      cgfsym = cgf_symbol(1,indco(1:3,i))
      name = TRIM(cgfsym(4:))//" DERIVATIVE OF THE OVERLAP MATRIX "//&
             "W.R.T. THE NUCLEAR COORDINATES"
      CALL compress(name)
      CALL uppercase(name)
      CALL replicate_matrix_structure(source=s(1)%matrix,&
                                      TARGET=s(i)%matrix,&
                                      target_name=TRIM(name),&
                                      target_symmetry="antisymmetric")
    END DO

!   *** Allocate the kinetic energy integral ***
!   *** matrix (only needed for printing or  ***
!   *** the OT preconditioner)               ***

    build_kinetic_energy_matrix = (scf_control%use_ot.OR.&
                                   globenv%print%kinetic_energy_matrix.OR.&
                                   globenv%print%potential_energy_matrix.OR.&
                                   globenv%run_type=="EP")

    IF (build_kinetic_energy_matrix) THEN
      IF (globenv%print%derivatives) THEN
        CALL allocate_matrix_set(t,maxder)
        return_t_derivatives = .TRUE.
      ELSE
        CALL allocate_matrix_set(t,1)
        return_t_derivatives = .FALSE.
      END IF
      CALL replicate_matrix_structure(source=s(1)%matrix,&
                                      TARGET=t(1)%matrix,&
                                      target_name="KINETIC ENERGY MATRIX")
      DO i=2,SIZE(t)
        cgfsym = cgf_symbol(1,indco(1:3,i))
        name = TRIM(cgfsym(4:))//" DERIVATIVE OF THE KINETIC ENERGY MATRIX "//&
               "W.R.T. THE NUCLEAR COORDINATES"
        CALL compress(name)
        CALL uppercase(name)
        CALL replicate_matrix_structure(source=s(1)%matrix,&
                                        TARGET=t(i)%matrix,&
                                        target_name=name,&
                                        target_symmetry="antisymmetric")
      END DO
    END IF

!   *** Allocate the core Hamiltonian matrix ***

    CALL allocate_matrix_set(h,1)

    CALL replicate_matrix_structure(source=s(1)%matrix,&
                                    TARGET=h(1)%matrix,&
                                    target_name="CORE HAMILTONIAN MATRIX")

!   *** Allocate work storage ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             gth_potential_present=gth_potential_present,&
                             maxco=maxco,&
                             maxlgto=maxlgto,&
                             maxlppl=maxlppl,&
                             maxlppnl=maxlppnl,&
                             maxppnl=maxppnl,&
                             maxsgf=maxsgf)

    maxl = MAX(maxlgto,maxlppl,maxlppnl)

    CALL init_orbital_pointers(maxl+nder+1)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxdco,&
                             maxder=MAX(1,nder))

    ldsab = MAX(maxco,ncoset(maxlppnl),maxsgf,maxppnl)
    maxblock = MAX(SIZE(s),maxder)

    ALLOCATE (hab(ldsab,ldsab*maxder),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,&
                       "hab",ldsab*ldsab*maxder*wp_size)
    END IF
    hab(:,:) = 0.0_wp

    ALLOCATE (sab(ldsab,ldsab*maxblock),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,&
                       "sab",ldsab*ldsab*maxblock*wp_size)
    END IF
    sab(:,:) = 0.0_wp

    ALLOCATE (sdab(maxdco,maxco,4),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,&
                       "sdab",maxdco*maxco*4*wp_size)
    END IF
    sdab(:,:,:) = 0.0_wp

    ALLOCATE (work(ldsab,ldsab*maxder),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,&
                       "work",ldsab*ldsab*maxder*wp_size)
    END IF
    work(:,:) = 0.0_wp

    IF (calculate_forces) THEN
      ALLOCATE (pab(maxco,maxco),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,&
                         "pab",maxco*maxco*wp_size)
      END IF
      pab(:,:) = 0.0_wp
    END IF

    ALLOCATE (sint(SIZE(s)),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,&
                       "sint",SIZE(s)*int_size)
    END IF
    DO i=1,SIZE(sint)
      NULLIFY (sint(i)%block)
    END DO

    IF (build_kinetic_energy_matrix) THEN
      ALLOCATE (tint(SIZE(t)),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "tint",SIZE(t)*int_size)
      DO i=1,SIZE(tint)
        NULLIFY (tint(i)%block)
      END DO
    END IF

!   *** Load GTH pseudo potential data (local part -> PPL) ***

    ppl_present = ASSOCIATED(sac_ppl)

    IF (ppl_present) THEN

      ALLOCATE (vppl(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "vppl",nkind*int_size)

      DO kkind=1,nkind

        atomic_kind => atomic_kind_set(kkind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             gth_potential=gth_potential)

        NULLIFY (vppl(kkind)%cexp_ppl)

        IF (ASSOCIATED(gth_potential)) THEN
          CALL get_potential(potential=gth_potential,&
                             alpha_ppl=vppl(kkind)%alpha_ppl,&
                             cexp_ppl=vppl(kkind)%cexp_ppl,&
                             ppl_radius=vppl(kkind)%ppl_radius)
        END IF

        NULLIFY (vppl(kkind)%r2)
        NULLIFY (vppl(kkind)%neighbor)
        NULLIFY (vppl(kkind)%r)

      END DO

    END IF

!   *** Load GTH pseudo potential data (non-local part -> PPNL) ***

    ppnl_present = ASSOCIATED(sac_ppnl)

    IF (ppnl_present) THEN

      ALLOCATE (vppnl(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "vppnl",nkind*int_size)

      DO kkind=1,nkind

        atomic_kind => atomic_kind_set(kkind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             gth_potential=gth_potential)

        NULLIFY (vppnl(kkind)%alpha_ppnl)
        NULLIFY (vppnl(kkind)%cprj)
        NULLIFY (vppnl(kkind)%nprj_ppnl)
        NULLIFY (vppnl(kkind)%vprj_ppnl)

        CALL get_potential(potential=gth_potential,&
                           alpha_ppnl=vppnl(kkind)%alpha_ppnl,&
                           cprj=vppnl(kkind)%cprj,&
                           lppnl=vppnl(kkind)%lppnl,&
                           nppnl=vppnl(kkind)%nppnl,&
                           nprj_ppnl=vppnl(kkind)%nprj_ppnl,&
                           ppnl_radius=vppnl(kkind)%ppnl_radius,&
                           vprj_ppnl=vppnl(kkind)%vprj_ppnl)

        NULLIFY (vppnl(kkind)%r2)
        NULLIFY (vppnl(kkind)%neighbor)
        NULLIFY (vppnl(kkind)%r)

        ALLOCATE (vppnl(kkind)%sac(maxder),STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine_name,module_name,__LINE__,&
                           "vppnl(kkind)%sac",maxder*int_size)
        END IF

        DO i=1,maxder
          NULLIFY (vppnl(kkind)%sac(i)%block)
        END DO

      END DO

    END IF

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           orb_basis_set=orb_basis_set)

      IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

      CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                             first_sgf=first_sgfa,&
                             lmax=la_max,&
                             lmin=la_min,&
                             nco_sum=ncoa_sum,&
                             npgf=npgfa,&
                             nset=nseta,&
                             nsgf_set=nsgfa,&
                             pgf_radius=rpgfa,&
                             set_radius=set_radius_a,&
                             sphi=sphi_a,&
                             zet=zeta)

      DO jkind=1,nkind

        atomic_kind => atomic_kind_set(jkind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             orb_basis_set=orb_basis_set)

        IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               first_sgf=first_sgfb,&
                               lmax=lb_max,&
                               lmin=lb_min,&
                               npgf=npgfb,&
                               nset=nsetb,&
                               nsgf_set=nsgfb,&
                               pgf_radius=rpgfb,&
                               set_radius=set_radius_b,&
                               sphi=sphi_b,&
                               zet=zetb)

        iab = ikind + nkind*(jkind - 1)

        IF (.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE

        sab_orb_neighbor_list => first_list(sab_orb(iab)%neighbor_list_set)

        CALL get_neighbor_list_set(neighbor_list_set=&
                                     sab_orb(iab)%neighbor_list_set,&
                                   nlist=nlist)

        DO ilist=1,nlist

          CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list,&
                                 atom=iatom,&
                                 nnode=nnode)

          IF (calculate_forces) atom_a = atom_of_kind(iatom)

          IF (ppl_present) THEN

!           *** Retrieve the data of the SAC_PPL neighbors of atom "iatom" ***

            DO kkind=1,nkind

              iac = ikind + nkind*(kkind - 1)

              IF (.NOT.ASSOCIATED(sac_ppl(iac)%neighbor_list_set)) CYCLE

              sac_ppl_neighbor_list =>&
                find_neighbor_list(neighbor_list_set=&
                                     sac_ppl(iac)%neighbor_list_set,&
                                   atom=iatom)

              CALL get_neighbor_list(neighbor_list=sac_ppl_neighbor_list,&
                                     nnode=nneighbor)

              ALLOCATE (vppl(kkind)%r2(nneighbor),STAT=istat)
              IF (istat /= 0) THEN
                CALL stop_memory(routine_name,module_name,__LINE__,&
                                 "vppl(kkind)%r2",nneighbor*wp_size)
              END IF

              ALLOCATE (vppl(kkind)%neighbor(nneighbor),STAT=istat)
              IF (istat /= 0) THEN
                CALL stop_memory(routine_name,module_name,__LINE__,&
                                 "vppl(kkind)%neighbor",nneighbor*int_size)
              END IF

              ALLOCATE (vppl(kkind)%r(3,nneighbor),STAT=istat)
              IF (istat /= 0) THEN
                CALL stop_memory(routine_name,module_name,__LINE__,&
                                 "vppl(kkind)%r",3*nneighbor*wp_size)
              END IF

              sac_ppl_neighbor_node => first_node(sac_ppl_neighbor_list)

              DO kneighbor=1,nneighbor
                CALL get_neighbor_node(neighbor_node=sac_ppl_neighbor_node,&
                                       neighbor=vppl(kkind)%neighbor(kneighbor),&
                                       r=vppl(kkind)%r(:,kneighbor))
                vppl(kkind)%r2(kneighbor) =&
                  vppl(kkind)%r(1,kneighbor)*vppl(kkind)%r(1,kneighbor) +&
                  vppl(kkind)%r(2,kneighbor)*vppl(kkind)%r(2,kneighbor) +&
                  vppl(kkind)%r(3,kneighbor)*vppl(kkind)%r(3,kneighbor)
                sac_ppl_neighbor_node => next(sac_ppl_neighbor_node)
              END DO

            END DO

          END IF

!         *** Loop over the GTH pseudo potential atoms (non-local part) ***

          IF (ppnl_present) THEN

            DO kkind=1,nkind

              iac = ikind + nkind*(kkind - 1)

              IF (.NOT.ASSOCIATED(sac_ppnl(iac)%neighbor_list_set)) CYCLE

              rprjc(1) = vppnl(kkind)%ppnl_radius

              sac_ppnl_neighbor_list =>&
                find_neighbor_list(neighbor_list_set=&
                                     sac_ppnl(iac)%neighbor_list_set,&
                                   atom=iatom)

              CALL get_neighbor_list(neighbor_list=sac_ppnl_neighbor_list,&
                                     nnode=nneighbor)

              ALLOCATE (vppnl(kkind)%r2(nneighbor),STAT=istat)
              IF (istat /= 0) THEN
                CALL stop_memory(routine_name,module_name,__LINE__,&
                                 "vppnl(kkind)%r2",nneighbor*wp_size)
              END IF

              ALLOCATE (vppnl(kkind)%neighbor(4,nneighbor),STAT=istat)
              IF (istat /= 0) THEN
                CALL stop_memory(routine_name,module_name,__LINE__,&
                                 "vppnl(kkind)%neighbor",4*nneighbor*int_size)
              END IF

              ALLOCATE (vppnl(kkind)%r(3,nneighbor),STAT=istat)
              IF (istat /= 0) THEN
                CALL stop_memory(routine_name,module_name,__LINE__,&
                                 "vppnl(kkind)%r",3*nneighbor*wp_size)
              END IF

              ncol = nneighbor*vppnl(kkind)%nppnl

              DO i=1,SIZE(vppnl(kkind)%sac)
                ALLOCATE (vppnl(kkind)%sac(i)%block(ncoa_sum,ncol),STAT=istat)
                IF (istat /= 0) THEN
                  CALL stop_memory(routine_name,module_name,__LINE__,&
                                   "vppnl(kkind)%sac(i)%block",&
                                   ncoa_sum*ncol*wp_size)
                END IF
              END DO

              ncol = 1

              sac_ppnl_neighbor_node => first_node(sac_ppnl_neighbor_list)

              DO kneighbor=1,nneighbor

                CALL get_neighbor_node(neighbor_node=sac_ppnl_neighbor_node,&
                                       neighbor=vppnl(kkind)%neighbor(1,kneighbor),&
                                       cell=vppnl(kkind)%neighbor(2:4,kneighbor),&
                                       r=rac)

                rac2 = rac(1)*rac(1) + rac(2)*rac(2) + rac(3)*rac(3)
                dac = SQRT(rac2)

                vppnl(kkind)%r(:,kneighbor) = rac(:)
                vppnl(kkind)%r2(kneighbor) = rac2

                nrow = 1

                DO iset=1,nseta

                  ncoa = npgfa(iset)*ncoset(la_max(iset))

                  prjc = 1

                  DO l=0,vppnl(kkind)%lppnl

                    nprjc = vppnl(kkind)%nprj_ppnl(l)*nco(l)

                    IF (nprjc == 0) CYCLE

                    IF (set_radius_a(iset) + rprjc(1) < dac) CYCLE

                    lc_max = l + 2*(vppnl(kkind)%nprj_ppnl(l) - 1)
                    lc_min = l
                    zetc(1) = vppnl(kkind)%alpha_ppnl(l)
                    ncoc = ncoset(lc_max)

!                   *** Calculate the primitive overlap integrals ***

                    CALL overlap(la_max(iset),zeta(1:npgfa(iset),iset),&
                                 rpgfa(1:npgfa(iset),iset),la_min(iset),&
                                 lc_max,zetc,rprjc,lc_min,&
                                 rac,rac2,dac,sab,nder,.TRUE.)

!                   *** Contraction step (projector functions) ***

                    DO i=1,maxder
                      first_col = (i - 1)*SIZE(work,1)
                      CALL dgemm("N","N",ncoa,nprjc,ncoc,&
                                 1.0_wp,sab(1,first_col+1),SIZE(sab,1),&
                                 vppnl(kkind)%cprj(1,prjc),&
                                 SIZE(vppnl(kkind)%cprj,1),&
                                 0.0_wp,work(1,first_col+prjc),SIZE(work,1))
                    END DO

                    prjc = prjc + nprjc

                  END DO

                  DO i=1,maxder
                    first_col = (i - 1)*SIZE(work,1) + 1
                    CALL dgemm("N","N",ncoa,vppnl(kkind)%nppnl,vppnl(kkind)%nppnl,&
                               1.0_wp,work(1,first_col),SIZE(work,1),&
                               vppnl(kkind)%vprj_ppnl(1,1),&
                               SIZE(vppnl(kkind)%vprj_ppnl,1),&
                               0.0_wp,vppnl(kkind)%sac(i)%block(nrow,ncol),&
                               SIZE(vppnl(kkind)%sac(i)%block,1))
                  END DO

                  nrow = nrow + ncoa

                END DO

                ncol = ncol + vppnl(kkind)%nppnl

                sac_ppnl_neighbor_node => next(sac_ppnl_neighbor_node)

              END DO ! kneighbor

            END DO ! kkind

          END IF

          last_jatom = 0

          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list)

          DO inode=1,nnode

            CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                   neighbor=jatom,&
                                   r=rab)

            IF (jatom /= last_jatom) THEN
              new_atom_b = .TRUE.
              last_jatom = jatom
            ELSE
              new_atom_b = .FALSE.
            END IF

            IF (calculate_forces) atom_b = atom_of_kind(jatom)

!           *** Use the symmetry of the first derivatives ***

            IF (iatom == jatom) THEN
              f0 = 1.0_wp
            ELSE
              f0 = 2.0_wp
            END IF

!           *** Create matrix blocks for a new matrix block column ***

            IF (new_atom_b) THEN
              IF (iatom <= jatom) THEN
                irow = iatom
                icol = jatom
              ELSE
                irow = jatom
                icol = iatom
              END IF
              DO i=1,SIZE(s)
                NULLIFY (sint(i)%block)
                CALL add_block_node(matrix=s(i)%matrix,&
                                    block_row=irow,&
                                    block_col=icol,&
                                    block=sint(i)%block)
              END DO
              IF (build_kinetic_energy_matrix) THEN
                DO i=1,SIZE(t)
                  NULLIFY (tint(i)%block)
                  CALL add_block_node(matrix=t(i)%matrix,&
                                      block_row=irow,&
                                      block_col=icol,&
                                      block=tint(i)%block)
                END DO
              END IF
              NULLIFY (h_block)
              CALL add_block_node(matrix=h(1)%matrix,&
                                  block_row=irow,&
                                  block_col=icol,&
                                  block=h_block)
              IF (calculate_forces) THEN
                CALL get_block_node(matrix=p(1)%matrix,&
                                    block_row=irow,&
                                    block_col=icol,&
                                    block=p_block)
                IF (.NOT.ASSOCIATED(p_block)) THEN
                  CALL stop_program(routine_name,module_name,__LINE__,&
                                    "A density matrix block is missing",&
                                    globenv)
                END IF
                CALL get_block_node(matrix=w(1)%matrix,&
                                    block_row=irow,&
                                    block_col=icol,&
                                    block=w_block)
                IF (.NOT.ASSOCIATED(w_block)) THEN
                  CALL stop_program(routine_name,module_name,__LINE__,&
                                    "A weighted density matrix block is missing",&
                                    globenv)
                END IF
              END IF
            END IF

            rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
            dab = SQRT(rab2)

            nrow = 1

            DO iset=1,nseta

              ncoa = npgfa(iset)*ncoset(la_max(iset))
              sgfa = first_sgfa(1,iset)

              DO jset=1,nsetb

                ncob = npgfb(jset)*ncoset(lb_max(jset))
                sgfb = first_sgfb(1,jset)

                IF (set_radius_a(iset) + set_radius_b(jset) >= dab) THEN

                  IF (calculate_forces.AND.(iatom /= jatom)) THEN

!                   *** Decontract W matrix block ***

                    IF (iatom <= jatom) THEN
                      CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
                                 1.0_wp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 w_block(sgfa,sgfb),SIZE(w_block,1),&
                                 0.0_wp,work(1,1),SIZE(work,1))
                    ELSE
                      CALL dgemm("N","T",ncoa,nsgfb(jset),nsgfa(iset),&
                                 1.0_wp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 w_block(sgfb,sgfa),SIZE(w_block,1),&
                                 0.0_wp,work(1,1),SIZE(work,1))
                    END IF
                    CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
                               1.0_wp,work(1,1),SIZE(work,1),&
                               sphi_b(1,sgfb),SIZE(sphi_b,1),&
                               0.0_wp,pab(1,1),SIZE(pab,1))

!                   *** Calculate the primitive overlap integrals ***
!                   *** and the corresponding force contribution  ***

                    CALL overlap(la_max(iset),zeta(1:npgfa(iset),iset),&
                                 rpgfa(1:npgfa(iset),iset),la_min(iset),&
                                 lb_max(jset),zetb(1:npgfb(jset),jset),&
                                 rpgfb(1:npgfb(jset),jset),lb_min(jset),&
                                 rab,rab2,dab,sab,nder,return_s_derivatives,&
                                 sdab,pab,force_a)

                    force(ikind)%overlap(:,atom_a) =&
                      force(ikind)%overlap(:,atom_a) - 2.0_wp*force_a(:)
                    force(jkind)%overlap(:,atom_b) =&
                      force(jkind)%overlap(:,atom_b) + 2.0_wp*force_a(:)

                  ELSE

!                   *** Calculate the primitive overlap integrals ***

                    CALL overlap(la_max(iset),zeta(1:npgfa(iset),iset),&
                                 rpgfa(1:npgfa(iset),iset),la_min(iset),&
                                 lb_max(jset),zetb(1:npgfb(jset),jset),&
                                 rpgfb(1:npgfb(jset),jset),lb_min(jset),&
                                 rab,rab2,dab,sab,nder,return_s_derivatives,&
                                 sdab)

                  END IF

!                 *** Contraction step (overlap matrix and its derivatives) ***

                  DO i=1,SIZE(s)
                    first_col = (i - 1)*SIZE(sab,1) + 1
                    CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                               1.0_wp,sab(1,first_col),SIZE(sab,1),&
                               sphi_b(1,sgfb),SIZE(sphi_b,1),&
                               0.0_wp,work(1,1),SIZE(work,1))
                    IF (iatom <= jatom) THEN
                      CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                                 1.0_wp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 work(1,1),SIZE(work,1),&
                                 1.0_wp,sint(i)%block(sgfa,sgfb),&
                                 SIZE(sint(i)%block,1))
                    ELSE
                      IF (i > 1) THEN
                        f = -1.0_wp
                      ELSE
                        f = 1.0_wp
                      END IF
                      CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                                 f,work(1,1),SIZE(work,1),&
                                 sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 1.0_wp,sint(i)%block(sgfb,sgfa),&
                                 SIZE(sint(i)%block,1))
                    END IF
                  END DO

!                 *** Decontract density matrix block ***

                  IF (calculate_forces) THEN
                    IF (iatom <= jatom) THEN
                      CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
                                 1.0_wp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 p_block(sgfa,sgfb),SIZE(p_block,1),&
                                 0.0_wp,work(1,1),SIZE(work,1))
                    ELSE
                      CALL dgemm("N","T",ncoa,nsgfb(jset),nsgfa(iset),&
                                 1.0_wp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 p_block(sgfb,sgfa),SIZE(p_block,1),&
                                 0.0_wp,work(1,1),SIZE(work,1))
                    END IF
                    CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
                               1.0_wp,work(1,1),SIZE(work,1),&
                               sphi_b(1,sgfb),SIZE(sphi_b,1),&
                               0.0_wp,pab(1,1),SIZE(pab,1))
                  END IF

                  IF (calculate_forces.AND.(iatom /= jatom)) THEN

!                   *** Calculate the primitive kinetic energy integrals ***
!                   *** and the corresponding force contribution         ***

                    CALL kinetic(la_max(iset),zeta(1:npgfa(iset),iset),&
                                 rpgfa(1:npgfa(iset),iset),la_min(iset),&
                                 lb_max(jset),zetb(1:npgfb(jset),jset),&
                                 rpgfb(1:npgfb(jset),jset),lb_min(jset),&
                                 rab,rab2,dab,sdab,hab,nder,&
                                 return_t_derivatives,pab,force_a)

                    force(ikind)%kinetic(:,atom_a) =&
                      force(ikind)%kinetic(:,atom_a) + 2.0_wp*force_a(:)
                    force(jkind)%kinetic(:,atom_b) =&
                      force(jkind)%kinetic(:,atom_b) - 2.0_wp*force_a(:)

                  ELSE

!                   *** Calculate the primitive kinetic energy integrals ***

                    CALL kinetic(la_max(iset),zeta(1:npgfa(iset),iset),&
                                 rpgfa(1:npgfa(iset),iset),la_min(iset),&
                                 lb_max(jset),zetb(1:npgfb(jset),jset),&
                                 rpgfb(1:npgfb(jset),jset),lb_min(jset),&
                                 rab,rab2,dab,sdab,hab,nder,&
                                 return_t_derivatives)

                  END IF

!                 *** Contraction step (kinetic energy integral ***
!                 *** matrix and its derivatives)               ***

                  IF (build_kinetic_energy_matrix) THEN
                    DO i=1,SIZE(t)
                      first_col = (i - 1)*SIZE(hab,1) + 1
                      CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                                 1.0_wp,hab(1,first_col),SIZE(hab,1),&
                                 sphi_b(1,sgfb),SIZE(sphi_b,1),&
                                 0.0_wp,work(1,1),SIZE(work,1))
                      IF (iatom <= jatom) THEN
                        CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                                   1.0_wp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                   work(1,1),SIZE(work,1),&
                                   1.0_wp,tint(i)%block(sgfa,sgfb),&
                                   SIZE(tint(i)%block,1))
                      ELSE
                        IF (i > 1) THEN
                          f = -1.0_wp
                        ELSE
                          f = 1.0_wp
                        END IF
                        CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                                   f,work(1,1),SIZE(work,1),&
                                   sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                   1.0_wp,tint(i)%block(sgfb,sgfa),&
                                   SIZE(tint(i)%block,1))
                      END IF
                    END DO
                  END IF

!                 *** Loop over the GTH pseudo potential atoms (local part) ***

                  IF (ppl_present) THEN

                    DO kkind=1,nkind

                      IF (.NOT.ASSOCIATED(vppl(kkind)%neighbor)) CYCLE

                      DO kneighbor=1,SIZE(vppl(kkind)%neighbor)

                        katom = vppl(kkind)%neighbor(kneighbor)
                        rac(:) = vppl(kkind)%r(:,kneighbor)
                        rac2 = vppl(kkind)%r2(kneighbor)
                        dac = SQRT(rac2)
                        IF (set_radius_a(iset) + vppl(kkind)%ppl_radius < dac) CYCLE

                        rbc(:) = rac(:) - rab(:)
                        rbc2 = rbc(1)*rbc(1) + rbc(2)*rbc(2) + rbc(3)*rbc(3)
                        dbc = SQRT(rbc2)
                        IF (set_radius_b(jset) + vppl(kkind)%ppl_radius < dbc) CYCLE

!                       *** Calculate the GTH pseudo potential forces ***

                        IF (calculate_forces) THEN

                          CALL overlap_ppl(&
                            la_max(iset),zeta(1:npgfa(iset),iset),&
                            rpgfa(1:npgfa(iset),iset),la_min(iset),&
                            lb_max(jset),zetb(1:npgfb(jset),jset),&
                            rpgfb(1:npgfb(jset),jset),lb_min(jset),&
                            vppl(kkind)%cexp_ppl,vppl(kkind)%alpha_ppl,&
                            vppl(kkind)%ppl_radius,&
                            rab,rab2,dab,&
                            rac,rac2,dac,&
                            rbc,rbc2,dbc,&
                            hab,nder,nder,.FALSE.,&
                            pab,force_a,force_b)

!                         *** The derivatives w.r.t. atomic center c are    ***
!                         *** calculated using the translational invariance ***
!                         *** of the first derivatives                      ***

                          atom_c = atom_of_kind(katom)

                          force(ikind)%gth_ppl(:,atom_a) =&
                            force(ikind)%gth_ppl(:,atom_a) + f0*force_a(:)
                          force(kkind)%gth_ppl(:,atom_c) =&
                            force(kkind)%gth_ppl(:,atom_c) - f0*force_a(:)

                          force(jkind)%gth_ppl(:,atom_b) =&
                            force(jkind)%gth_ppl(:,atom_b) + f0*force_b(:)
                          force(kkind)%gth_ppl(:,atom_c) =&
                            force(kkind)%gth_ppl(:,atom_c) - f0*force_b(:)

                        ELSE

                          CALL overlap_ppl(&
                            la_max(iset),zeta(1:npgfa(iset),iset),&
                            rpgfa(1:npgfa(iset),iset),la_min(iset),&
                            lb_max(jset),zetb(1:npgfb(jset),jset),&
                            rpgfb(1:npgfb(jset),jset),lb_min(jset),&
                            vppl(kkind)%cexp_ppl,vppl(kkind)%alpha_ppl,&
                            vppl(kkind)%ppl_radius,&
                            rab,rab2,dab,&
                            rac,rac2,dac,&
                            rbc,rbc2,dbc,&
                            hab,0,0,.FALSE.)

                        END IF

                      END DO

                    END DO

                  END IF

                ELSE

                  DO jco=1,ncob
                    DO ico=1,ncoa
                      hab(ico,jco) = 0.0_wp
                    END DO
                  END DO

                  IF (calculate_forces) THEN

!                   *** Decontract density matrix block ***

                    IF (iatom <= jatom) THEN
                      CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
                                 1.0_wp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 p_block(sgfa,sgfb),SIZE(p_block,1),&
                                 0.0_wp,work(1,1),SIZE(work,1))
                    ELSE
                      CALL dgemm("N","T",ncoa,nsgfb(jset),nsgfa(iset),&
                                 1.0_wp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 p_block(sgfb,sgfa),SIZE(p_block,1),&
                                 0.0_wp,work(1,1),SIZE(work,1))
                    END IF
                    CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
                               1.0_wp,work(1,1),SIZE(work,1),&
                               sphi_b(1,sgfb),SIZE(sphi_b,1),&
                               0.0_wp,pab(1,1),SIZE(pab,1))

                  END IF

                END IF

!               *** Loop over the GTH pseudo potential atoms (non-local part) ***

                IF (ppnl_present.AND.new_atom_b) THEN

                  DO kkind=1,nkind

                    IF (.NOT.ASSOCIATED(vppnl(kkind)%neighbor)) CYCLE

                    rprjc(1) = vppnl(kkind)%ppnl_radius

                    ibc = jkind + nkind*(kkind - 1)

                    IF (.NOT.ASSOCIATED(sbc_ppnl(ibc)%neighbor_list_set)) CYCLE

                    sbc_ppnl_neighbor_list =>&
                      find_neighbor_list(neighbor_list_set=&
                                         sbc_ppnl(ibc)%neighbor_list_set,&
                                         atom=jatom)

                    CALL get_neighbor_list(neighbor_list=sbc_ppnl_neighbor_list,&
                                           nnode=nnode_c)

                    sbc_ppnl_neighbor_node => first_node(sbc_ppnl_neighbor_list)

                    DO knode=1,nnode_c

                      CALL get_neighbor_node(neighbor_node=sbc_ppnl_neighbor_node,&
                                             neighbor=katom,&
                                             cell=cell_c,&
                                             r=rbc)

!                     *** Search operator atom in AC list ***

                      search: DO kneighbor=1,SIZE(vppnl(kkind)%r2)

                        IF ((vppnl(kkind)%neighbor(1,kneighbor) == katom).AND.&
                            (vppnl(kkind)%neighbor(2,kneighbor) == cell_c(1)).AND.&
                            (vppnl(kkind)%neighbor(3,kneighbor) == cell_c(2)).AND.&
                            (vppnl(kkind)%neighbor(4,kneighbor) == cell_c(3))) THEN

                          rac(:) = vppnl(kkind)%r(:,kneighbor)
                          rac2 = vppnl(kkind)%r2(kneighbor)
                          dac = SQRT(rac2)
                          IF (set_radius_a(iset) + rprjc(1) < dac) EXIT search

                          rbc2 = rbc(1)*rbc(1) + rbc(2)*rbc(2) + rbc(3)*rbc(3)
                          dbc = SQRT(rbc2)
                          IF (set_radius_b(jset) + rprjc(1) < dbc) EXIT search

                          prjc = 1

                          DO l=0,vppnl(kkind)%lppnl

                            nprjc = vppnl(kkind)%nprj_ppnl(l)*nco(l)

                            IF (nprjc == 0) CYCLE

                            lc_max = l + 2*(vppnl(kkind)%nprj_ppnl(l) - 1)
                            lc_min = l
                            zetc(1) = vppnl(kkind)%alpha_ppnl(l)
                            ncoc = ncoset(lc_max)

!                           *** Calculate the primitive overlap integrals ***

                            CALL overlap(lb_max(jset),zetb(1:npgfb(jset),jset),&
                                         rpgfb(1:npgfb(jset),jset),lb_min(jset),&
                                         lc_max,zetc,rprjc,lc_min,&
                                         rbc,rbc2,dbc,sab,nder,.TRUE.)

!                           *** Contraction step (projector functions) ***

                            DO i=1,maxder
                              first_col = (i - 1)*SIZE(work,1)
                              CALL dgemm("N","N",ncob,nprjc,ncoc,&
                                         1.0_wp,sab(1,first_col+1),SIZE(sab,1),&
                                         vppnl(kkind)%cprj(1,prjc),&
                                         SIZE(vppnl(kkind)%cprj,1),&
                                         0.0_wp,work(1,first_col+prjc),&
                                         SIZE(work,1))
                            END DO

                            prjc = prjc + nprjc

                          END DO

                          ncol = (kneighbor - 1)*vppnl(kkind)%nppnl + 1

                          CALL dgemm("N","T",ncoa,ncob,vppnl(kkind)%nppnl,&
                                     1.0_wp,vppnl(kkind)%sac(1)%block(nrow,ncol),&
                                     SIZE(vppnl(kkind)%sac(1)%block,1),&
                                     work(1,1),SIZE(work,1),&
                                     1.0_wp,hab(1,1),SIZE(hab,1))

                          IF (calculate_forces) THEN

                            atom_c = atom_of_kind(katom)

                            DO i=1,3

                              first_col = i*SIZE(work,1) + 1

                              CALL dgemm("N","T",ncoa,ncob,vppnl(kkind)%nppnl,&
                                         1.0_wp,vppnl(kkind)%sac(i+1)%block(nrow,ncol),&
                                         SIZE(vppnl(kkind)%sac(i+1)%block,1),&
                                         work(1,1),SIZE(work,1),&
                                         0.0_wp,hab(1,first_col),SIZE(hab,1))

                              f = 0.0_wp

                              DO jco=1,ncob
                                j = first_col + jco - 1
                                DO ico=1,ncoa
                                  f = f + pab(ico,jco)*hab(ico,j)
                                END DO
                              END DO

                              force(ikind)%gth_ppnl(i,atom_a) =&
                                force(ikind)%gth_ppnl(i,atom_a) + f0*f

                              force(kkind)%gth_ppnl(i,atom_c) =&
                                force(kkind)%gth_ppnl(i,atom_c) - f0*f

                              CALL dgemm("N","T",ncoa,ncob,vppnl(kkind)%nppnl,&
                                         1.0_wp,&
                                         vppnl(kkind)%sac(1)%block(nrow,ncol),&
                                         SIZE(vppnl(kkind)%sac(1)%block,1),&
                                         work(1,first_col),SIZE(work,1),&
                                         0.0_wp,hab(1,first_col),SIZE(hab,1))

                              f = 0.0_wp

                              DO jco=1,ncob
                                j = first_col + jco - 1
                                DO ico=1,ncoa
                                  f = f + pab(ico,jco)*hab(ico,j)
                                END DO
                              END DO

                              force(jkind)%gth_ppnl(i,atom_b) =&
                                force(jkind)%gth_ppnl(i,atom_b) + f0*f

                              force(kkind)%gth_ppnl(i,atom_c) =&
                                force(kkind)%gth_ppnl(i,atom_c) - f0*f

                            END DO

                          END IF

                          EXIT search

                        END IF

                      END DO search

                      sbc_ppnl_neighbor_node => next(sbc_ppnl_neighbor_node)

                    END DO ! knode

                  END DO ! kkind

                END IF ! ppnl

!               *** Contraction step (core Hamiltonian matrix) ***

                CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                           1.0_wp,hab(1,1),SIZE(hab,1),&
                           sphi_b(1,sgfb),SIZE(sphi_b,1),&
                           0.0_wp,work(1,1),SIZE(work,1))

                IF (iatom <= jatom) THEN
                  CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                             1.0_wp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                             work(1,1),SIZE(work,1),&
                             1.0_wp,h_block(sgfa,sgfb),SIZE(h_block,1))
                ELSE
                  CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                             1.0_wp,work(1,1),SIZE(work,1),&
                             sphi_a(1,sgfa),SIZE(sphi_a,1),&
                             1.0_wp,h_block(sgfb,sgfa),SIZE(h_block,1))
                END IF

              END DO ! jset

              nrow = nrow + ncoa

            END DO ! iset

            sab_orb_neighbor_node => next(sab_orb_neighbor_node)

          END DO ! jatom => atom B

          IF (ppl_present) THEN
            DO kkind=1,nkind
              IF (ASSOCIATED(vppl(kkind)%r2)) THEN
                DEALLOCATE (vppl(kkind)%r2,STAT=istat)
                IF (istat /= 0) THEN
                  CALL stop_memory(routine_name,module_name,__LINE__,&
                                   "vppl(kkind)%r2")
                END IF
                DEALLOCATE (vppl(kkind)%neighbor,STAT=istat)
                IF (istat /= 0) THEN
                  CALL stop_memory(routine_name,module_name,__LINE__,&
                                   "vppl(kkind)%neighbor")
                END IF
                DEALLOCATE (vppl(kkind)%r,STAT=istat)
                IF (istat /= 0) THEN
                  CALL stop_memory(routine_name,module_name,__LINE__,&
                                   "vppl(kkind)%r")
                END IF
              END IF
            END DO
          END IF

          IF (ppnl_present) THEN
            DO kkind=1,nkind
              IF (ASSOCIATED(vppnl(kkind)%r2)) THEN
                DEALLOCATE (vppnl(kkind)%r2,STAT=istat)
                IF (istat /= 0) THEN
                  CALL stop_memory(routine_name,module_name,__LINE__,&
                                   "vppnl(kkind)%r2")
                END IF
                DEALLOCATE (vppnl(kkind)%neighbor,STAT=istat)
                IF (istat /= 0) THEN
                  CALL stop_memory(routine_name,module_name,__LINE__,&
                                   "vppnl(kkind)%neighbor")
                END IF
                DEALLOCATE (vppnl(kkind)%r,STAT=istat)
                IF (istat /= 0) THEN
                  CALL stop_memory(routine_name,module_name,__LINE__,&
                                   "vppnl(kkind)%r")
                END IF
                DO i=1,SIZE(vppnl(kkind)%sac)
                  DEALLOCATE (vppnl(kkind)%sac(i)%block,STAT=istat)
                  IF (istat /= 0) THEN
                    CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "vppnl(kkind)%sac(i)%block")
                  END IF
                END DO
              END IF
            END DO
          END IF

          sab_orb_neighbor_list => next(sab_orb_neighbor_list)

        END DO ! iatom => atom A

      END DO ! jkind

    END DO ! ikind

!   *** Put the core Hamiltonian matrix in the QS environment ***

    CALL set_qs_env(qs_env=qs_env,&
                    s=s,&
                    kinetic=t,&
                    h=h)

!   *** Release work storage ***

    IF (ALLOCATED(vppl)) THEN
      DO kkind=1,nkind
        NULLIFY (vppl(kkind)%cexp_ppl)
      END DO
      DEALLOCATE (vppl,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "vppl")
    END IF

    IF (ALLOCATED(vppnl)) THEN
      DO kkind=1,nkind
        NULLIFY (vppnl(kkind)%alpha_ppnl)
        NULLIFY (vppnl(kkind)%cprj)
        NULLIFY (vppnl(kkind)%nprj_ppnl)
        NULLIFY (vppnl(kkind)%vprj_ppnl)
        DEALLOCATE (vppnl(kkind)%sac,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "vppnl(kkind)%sac")
      END DO
      DEALLOCATE (vppnl,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "vppnl")
    END IF

    DEALLOCATE (hab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "hab")

    DEALLOCATE (sab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "sab")

    DEALLOCATE (sdab,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "sdab")

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "work")

    IF (calculate_forces) THEN

      DEALLOCATE (atom_of_kind,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "atom_of_kind")

      DEALLOCATE (pab,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "pab")

!     *** If LSD, then recover alpha density and beta density     ***
!     *** from the total density (1) and the spin density (2)     ***
!     *** The W matrix is neglected, since it will be destroyed   ***
!     *** in the calling force routine after leaving this routine ***

      IF (SIZE(p) == 2) THEN
        CALL cp_sm_scale_and_add(p(1)%matrix, 0.5_wp,p(2)%matrix,0.5_wp)
        CALL cp_sm_scale_and_add(p(2)%matrix,-1.0_wp,p(1)%matrix,1.0_wp)
      END IF

    END IF

    DO i=1,SIZE(sint)
      NULLIFY (sint(i)%block)
    END DO
    DEALLOCATE (sint,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "sint")

    IF (build_kinetic_energy_matrix) THEN
      DO i=1,SIZE(tint)
        NULLIFY (tint(i)%block)
      END DO
      DEALLOCATE (tint,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "tint")
    END IF

    CALL timestop(0.0_wp,handle)

!   *** Print the distribution of the overlap matrix blocks ***

    CALL write_matrix_distribution(s(1)%matrix,globenv)

!   *** Print the overlap integral matrix, if requested ***

    IF (globenv%print%overlap_matrix_magnitude) THEN
      IF (globenv%num_pe > 1) THEN
        WRITE (UNIT=extension,FMT="(I6)") globenv%mepos
        file_name = "S_MAGNITUDE"//"."//ADJUSTL(extension)
      ELSE
        file_name = "S_MAGNITUDE"
      END IF
      CALL open_file(file_name=file_name,&
                     file_action="WRITE",&
                     file_form="FORMATTED",&
                     file_status="REPLACE",&
                     unit_number=unit_number)
      CALL write_blocks_maxabsval(s(1)%matrix,unit_number)
      CALL close_file(unit_number=unit_number)
    END IF

    IF (globenv%print%overlap_matrix) THEN
      CALL write_sparse_matrix(s(1)%matrix,4,6,qs_env,globenv)
      IF (globenv%print%derivatives) THEN
        DO i=2,SIZE(s)
          CALL write_sparse_matrix(s(i)%matrix,4,6,qs_env,globenv)
        END DO
      END IF
    END IF

!   *** Print the kinetic energy integral matrix, if requested ***

    IF (globenv%print%kinetic_energy_matrix) THEN
      CALL write_sparse_matrix(t(1)%matrix,4,6,qs_env,globenv)
      IF (globenv%print%derivatives) THEN
        DO i=2,SIZE(t)
          CALL write_sparse_matrix(t(i)%matrix,4,6,qs_env,globenv)
        END DO
      END IF
    END IF

!   *** Print the potential energy matrix, if requested ***

    IF (globenv%print%potential_energy_matrix) THEN
      CALL allocate_matrix_set(v,1)
      CALL replicate_matrix_structure(source=h(1)%matrix,&
                                      TARGET=v(1)%matrix,&
                                      target_name="POTENTIAL ENERGY MATRIX")
      CALL add_matrices(v(1)%matrix,1.0_wp,h(1)%matrix,-1.0_wp,t(1)%matrix)
      CALL write_sparse_matrix(v(1)%matrix,4,6,qs_env,globenv)
      CALL deallocate_matrix_set(v)
    END IF

!   *** Print the core Hamiltonian matrix, if requested ***

    IF (globenv%print%core_hamiltonian_matrix) THEN
      CALL write_sparse_matrix(h(1)%matrix,4,6,qs_env,globenv)
    END IF

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE build_core_hamiltonian_matrix

! *****************************************************************************

END MODULE qs_core_hamiltonian
