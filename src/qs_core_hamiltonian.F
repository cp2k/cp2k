!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2004  CP2K developers group
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
!!     - Joost VandeVondele (April 2003) : added LSD forces
!!     - Non-redundant calculation of the non-local part of the GTH PP
!!       (22.05.2003,MK)
!!     - New parallelization scheme (27.06.2003,MK)
!!     - OpenMP version (07.12.2003,JGH)
!!     - Binary search loop for VPPNL operators (09.01.2004,JGH,MK)
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
  USE external_potential_types,        ONLY: all_potential_type,&
                                             get_potential,&
                                             gth_potential_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: indco,&
                                             init_orbital_pointers,&
                                             nco,&
                                             ncoset
  USE orbital_symbols,                 ONLY: cgf_symbol
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE qs_all_potential,                 ONLY: all_integrals,&
                                             vall_type 
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_neighbor_list_types,          ONLY: find_neighbor_list,&
                                             first_list,&
                                             first_node,&
                                             get_neighbor_list,&
                                             get_neighbor_list_set,&
                                             get_neighbor_node,&
                                             qlist_type,&
                                             block_p_type,& 
                                             neighbor_list_set_p_type,&
                                             neighbor_list_type,&
                                             neighbor_node_type,&
                                             next, &
                                             reduced_3c_list_type
  USE qs_oce_methods,                  ONLY: set_up_oce
  USE qs_oce_types,                    ONLY: allocate_oce_set,&
                                             allocate_vtriple,&
                                             create_oce_set,&
                                             oce_matrix_type, &
                                             prepare_oce_coeff,&
                                             build_reduced_3c_lists,&
                                             retrieve_sac_list,&
                                             vtriple_type
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
                                             deallocate_matrix_set,&
                                             get_block_node,&
                                             real_matrix_p_type,&
                                             replicate_matrix_structure,&
                                             write_blocks_maxabsval
  USE string_utilities,                ONLY: compress,&
                                             uppercase
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: locate,&
                                             sort

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_core_hamiltonian"

!   *** Local types ***

    TYPE vppl_type
      REAL(KIND = dp)                          :: alpha_ppl,ppl_radius
      REAL(KIND = dp), DIMENSION(:), POINTER   :: cexp_ppl,r2
      REAL(KIND = dp), DIMENSION(:,:), POINTER :: r
      INTEGER, DIMENSION(:), POINTER           :: neighbor
    END TYPE vppl_type

    TYPE vppnl_type
      TYPE(qlist_type), DIMENSION(:,:,:), POINTER :: neighbor
      REAL(KIND = dp)                            :: ppnl_radius
      REAL(KIND = dp), DIMENSION(:), POINTER     :: alpha_ppnl
      REAL(KIND = dp), DIMENSION(:,:), POINTER   :: cprj,vprj_ppnl
      INTEGER                                    :: lppnl,nppnl
      INTEGER, DIMENSION(:), POINTER             :: nprj_ppnl
    END TYPE vppnl_type


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

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    LOGICAL, INTENT(IN)                      :: calculate_forces

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "core_hamiltonian"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER     :: atomic_kind
    TYPE(distribution_2d_type), POINTER :: distribution_2d
    TYPE(gth_potential_type), POINTER   :: gth_potential
    TYPE(all_potential_type), POINTER   :: all_potential
    TYPE(gto_basis_set_type), POINTER   :: orb_basis_set
    TYPE(qlist_type), POINTER            :: sac_ppnl_neighbor
    TYPE(neighbor_list_type), POINTER   :: sab_orb_neighbor_list,&
                                           sab_orb_neighbor_list_local,&
                                           sac_ppl_neighbor_list,&
                                           sac_ppnl_neighbor_list,&
                                           sbc_ppnl_neighbor_list
    TYPE(neighbor_node_type), POINTER   :: sab_orb_neighbor_node,&
                                           sac_ppl_neighbor_node,&
                                           sac_ppnl_neighbor_node,&
                                           sbc_ppnl_neighbor_node
    TYPE(qs_rho_type), POINTER          :: rho
    TYPE(scf_control_type), POINTER     :: scf_control

    CHARACTER(LEN=6)   :: extension
    CHARACTER(LEN=12)  :: cgfsym
    CHARACTER(LEN=80)  :: name
    CHARACTER(LEN=200) :: file_name
    REAL(KIND = dp)    :: dab,dac,dbc,f,f0,rab2,rac2,rbc2
    INTEGER            :: atom_a,atom_b,atom_c,first_col,handle,i,iab,&
                          iac,iat,iatom,ibc,icell,ico,icol,idx,ikind,ilist,&
                          inode,ipgf,irow,iset,istat,j,jatom,jcell,jco,&
                          jkind,jset,katom,kcell,kkind,kneighbor,knode,&
                          l,last_jatom,lc_max,lc_min,ldai,ldsab,lppnl,maxblock,&
                          maxco,maxdco,maxder,maxl,maxlgto,maxlppl,maxlppnl,&
                          maxppnl,maxsgf,n,nat,natom,ncoa,ncoa_sum,ncob,ncoc,ncol,&
                          nder,nkind,nlist,nneighbor,nnode,nnode_c,nppnl,nprjc,&
                          nrow,nseta,nsetb,nsgf,prjc,sgfa,sgfb,unit_number
    LOGICAL            :: all_potential_present,&
                          build_kinetic_energy_matrix,gth_potential_present,&
                          new_atom_b,ppl_present,ppnl_present,&
                          return_s_derivatives,return_t_derivatives,&
                          oce_present, rho0_present

    REAL(KIND = dp), DIMENSION(1) :: rprjc,zetc
    REAL(KIND = dp), DIMENSION(3) :: force_a,force_b,rab,rac,rbc
    INTEGER, DIMENSION(3)         :: cell_c,cell_c_max,cell_c_min

    TYPE(block_p_type), DIMENSION(:), POINTER    :: oceh,oces
    TYPE(block_p_type), DIMENSION(:,:), POINTER  :: oceht,ocest
    TYPE(block_p_type), DIMENSION(:), POINTER    :: sint,tint
    TYPE(block_p_type), DIMENSION(:,:), POINTER  :: sintt,tintt
    TYPE(vppl_type), DIMENSION(:), POINTER       :: vppl
    TYPE(vppnl_type), DIMENSION(:), POINTER      :: vppnl
    TYPE(vall_type), DIMENSION(:), POINTER       :: vall
    TYPE(vppl_type), DIMENSION(:,:), POINTER     :: vpplt
    TYPE(vppnl_type), DIMENSION(:,:), POINTER    :: vppnlt
    TYPE(vall_type), DIMENSION(:,:), POINTER     :: vallt
    TYPE(qlist_type), DIMENSION(:,:,:), POINTER   :: neighbor_att

    REAL(KIND = dp), DIMENSION(:,:), POINTER     :: hab,pab,sab,work
    REAL(KIND = dp), DIMENSION(:,:,:), POINTER   :: sdab
    REAL(KIND = dp), DIMENSION(:,:,:), POINTER   :: habt,pabt,sabt,workt
    REAL(KIND = dp), DIMENSION(:,:,:,:), POINTER :: sdabt

    REAL(KIND = dp), DIMENSION(:,:,:,:), ALLOCATABLE :: ai_work

    INTEGER, DIMENSION(:), ALLOCATABLE :: atom_of_kind,&
                                          first_sgf,&
                                          last_prj,&
                                          last_sgf

    TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb,sac_ppl,&
                                                             sac_ppnl,sbc_ppnl,&
                                                             sac_oce,sbc_oce,&
                                                             sac_3c,sbc_3c
    TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set
    TYPE(qs_force_type), DIMENSION(:), POINTER            :: force
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER       :: matrix_h,matrix_p,matrix_s, &
                                                             matrix_t,matrix_v,matrix_w
    TYPE(oce_matrix_type), POINTER                        :: oce
    TYPE(vtriple_type), DIMENSION(:,:),POINTER            :: VH_3c_list
    TYPE(reduced_3c_list_type), DIMENSION(:), POINTER     :: reduced_3c_oce,&
                                                             reduced_3c_rho0 

    REAL(KIND = dp), DIMENSION(:), POINTER   :: set_radius_a,set_radius_b
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: h_block,p_block,rpgfa,rpgfb,&
                                                sphi_a,sphi_b,w_block,zeta,zetb

    INTEGER, DIMENSION(:), POINTER    :: atom_list,la_max,la_min,lb_max,lb_min,&
                                         npgfa,npgfb,nsgfa,nsgfb
    INTEGER, DIMENSION(:,:), POINTER  :: first_sgfa,first_sgfb

    INTEGER :: ithread,nthread
    INTEGER :: omp_get_max_threads,omp_get_thread_num

#if defined(__AIX) && defined(__INIT_WITH_NANS)
    REAL(dp) :: nans_dp
    nans_dp = z"7FF0000000000001"
#endif

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    NULLIFY (atomic_kind_set)
    NULLIFY (distribution_2d)
    NULLIFY (force)
    NULLIFY (matrix_h)
    NULLIFY (habt)
    NULLIFY (oce)
    NULLIFY (VH_3c_list)
    NULLIFY (matrix_p)
    NULLIFY (pabt)
    NULLIFY (particle_set)
    NULLIFY (rho)
    NULLIFY (matrix_s)
    NULLIFY (sab_orb)
    NULLIFY (sabt)
    NULLIFY (sac_ppl)
    NULLIFY (sac_ppnl)
    NULLIFY (sbc_ppnl)
    NULLIFY (sac_oce)
    NULLIFY (sbc_oce)
    NULLIFY (sac_3c)
    NULLIFY (sbc_3c)
    NULLIFY (scf_control)
    NULLIFY (sdabt)
    NULLIFY (matrix_t)
    NULLIFY (matrix_w)
    NULLIFY (workt)
! Default
    oce_present = .FALSE.
    rho0_present = .FALSE.

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    matrix_h=matrix_h,&
                    kinetic=matrix_t,&
                    oce=oce, &
                    particle_set=particle_set,&
                    matrix_s=matrix_s,&
                    sab_orb=sab_orb,&
                    sac_ppl=sac_ppl,&
                    sac_ppnl=sac_ppnl,&
                    sbc_ppnl=sbc_ppnl,&
                    sac_oce=sac_oce,&
                    sbc_oce=sbc_oce,&
                    VH_3c_list=VH_3c_list,&
                    sac_3c=sac_3c,&
                    sbc_3c=sbc_3c,&
                    scf_control=scf_control,&
                    distribution_2d=distribution_2d)

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    IF (calculate_forces) THEN

      CALL timeset(routine_name//" (forces)","I"," ",handle)

      nder = 1

      ALLOCATE (atom_of_kind(natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "atom_of_kind",natom*int_size)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               atom_of_kind=atom_of_kind)

      CALL get_qs_env(qs_env=qs_env,&
                      force=force,&
                      rho=rho,&
                      matrix_w=matrix_w)

      matrix_p => rho%rho_ao

!     *** If LSD, then combine alpha density and beta density to
!     *** total density: alpha <- alpha + beta   and
!     *** spin density:   beta <- alpha - beta
!     (since all things can be computed based on the sum of these matrices anyway)
!     (matrix_p is restored at the end of the run, matrix_w is left in its modified state
!     (as it should not be needed afterwards)

      IF (SIZE(matrix_p) == 2) THEN
        CALL cp_sm_scale_and_add(matrix_p(1)%matrix, 1.0_dp,matrix_p(2)%matrix,1.0_dp)
        CALL cp_sm_scale_and_add(matrix_p(2)%matrix,-2.0_dp,matrix_p(1)%matrix,1.0_dp)
        CALL cp_sm_scale_and_add(matrix_w(1)%matrix, 1.0_dp,matrix_w(2)%matrix,1.0_dp)
        CALL cp_sm_scale_and_add(matrix_w(2)%matrix,-2.0_dp,matrix_w(1)%matrix,1.0_dp)
      END IF

    ELSE

      CALL timeset(routine_name,"I"," ",handle)

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

    CALL allocate_matrix_set(matrix_s,maxder)

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

    CALL allocate_matrix(matrix=matrix_s(1)%matrix,&
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

!   *** Allocate the matrix of coefficients for one center expansions
    IF(qs_env%dft_control%qs_control%gapw) THEN
      CALL create_oce_set(oce)
      CALL allocate_oce_set(oce,natom,nkind,maxder)
      CALL prepare_oce_coeff(particle_set,oce,natom,nsgf,&
                             first_sgf,last_sgf,distribution_2d)
    ENDIF

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
      CALL replicate_matrix_structure(source=matrix_s(1)%matrix,&
                                      target=matrix_s(i)%matrix,&
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
    ! assume false by default
    return_t_derivatives = .FALSE.

    IF (build_kinetic_energy_matrix) THEN
      IF (globenv%print%derivatives) THEN
        CALL allocate_matrix_set(matrix_t,maxder)
        return_t_derivatives = .TRUE.
      ELSE
        CALL allocate_matrix_set(matrix_t,1)
        return_t_derivatives = .FALSE.
      END IF
      CALL replicate_matrix_structure(source=matrix_s(1)%matrix,&
                                      target=matrix_t(1)%matrix,&
                                      target_name="KINETIC ENERGY MATRIX")
      DO i=2,SIZE(matrix_t)
        cgfsym = cgf_symbol(1,indco(1:3,i))
        name = TRIM(cgfsym(4:))//" DERIVATIVE OF THE KINETIC ENERGY MATRIX "//&
               "W.R.T. THE NUCLEAR COORDINATES"
        CALL compress(name)
        CALL uppercase(name)
        CALL replicate_matrix_structure(source=matrix_s(1)%matrix,&
                                        target=matrix_t(i)%matrix,&
                                        target_name=name,&
                                        target_symmetry="antisymmetric")
      END DO
    END IF

!   *** Allocate the core Hamiltonian matrix ***

    CALL allocate_matrix_set(matrix_h,1)

    CALL replicate_matrix_structure(source=matrix_s(1)%matrix,&
                                    target=matrix_h(1)%matrix,&
                                    target_name="CORE HAMILTONIAN MATRIX")

!   *** Allocate work storage ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             gth_potential_present=gth_potential_present,&
                             all_potential_present=all_potential_present, &
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
    maxblock = MAX(SIZE(matrix_s),maxder)
    nthread = 1
!$  nthread = omp_get_max_threads()
!MK
    ldai = ncoset(maxl+nder+1)

    CALL reallocate(habt,1,ldsab,1,ldsab*maxder,0,nthread-1)
    CALL reallocate(sabt,1,ldsab,1,ldsab*maxblock,0,nthread-1)
    CALL reallocate(sdabt,1,maxdco,1,maxco,1,4,0,nthread-1)
    CALL reallocate(workt,1,ldsab,1,ldsab*maxder,0,nthread-1)
    IF (calculate_forces) THEN
      CALL reallocate(pabt,1,maxco,1,maxco,0,nthread-1)
    END IF

    ALLOCATE (sintt(SIZE(matrix_s),0:nthread-1),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,&
                       "sintt",SIZE(matrix_s)*int_size*nthread)
    END IF
    DO i=1,SIZE(sintt,1)
      DO j=0,nthread-1
        NULLIFY (sintt(i,j)%block)
      END DO
    END DO

    IF (build_kinetic_energy_matrix) THEN
      ALLOCATE (tintt(SIZE(matrix_t),0:nthread-1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "tintt",SIZE(matrix_t)*int_size*nthread)
      DO i=1,SIZE(tintt,1)
        DO j=0,nthread-1
          NULLIFY (tintt(i,j)%block)
        END DO
      END DO
    END IF

!   ***  Initialize blocks for the one center expansione ***
    IF(qs_env%dft_control%qs_control%gapw) THEN
       oce_present = ASSOCIATED(sac_oce)
       rho0_present = ASSOCIATED(sac_3c)
       IF(rho0_present) THEN
         CALL allocate_vtriple(VH_3c_list,nkind,natom)
       END IF

       ALLOCATE (oceht(SIZE(matrix_s),0:nthread-1),STAT=istat)
       IF (istat /= 0) THEN
         CALL stop_memory(routine_name,module_name,__LINE__,&
                       "oceht",SIZE(matrix_s)*int_size*nthread)
       END IF
       ALLOCATE (ocest(SIZE(matrix_s),0:nthread-1),STAT=istat)
       IF (istat /= 0) THEN
         CALL stop_memory(routine_name,module_name,__LINE__,&
                       "ocest",SIZE(matrix_s)*int_size*nthread)
       END IF

       DO i=1,SIZE(oceht,1)
         DO j=0,nthread-1
           NULLIFY (oceht(i,j)%block)
           NULLIFY (ocest(i,j)%block)
         END DO
       END DO
     ENDIF

!   *** Load GTH pseudo potential data (local part -> PPL) ***

    ppl_present = ASSOCIATED(sac_ppl)

    IF (ppl_present) THEN

      ALLOCATE (vpplt(nkind,0:nthread-1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "vpplt",nkind*int_size)

      DO kkind=1,nkind

        atomic_kind => atomic_kind_set(kkind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             gth_potential=gth_potential)

        DO j=0,nthread-1
          NULLIFY (vpplt(kkind,j)%cexp_ppl)
          IF (ASSOCIATED(gth_potential)) THEN
            CALL get_potential(potential=gth_potential,&
                               alpha_ppl=vpplt(kkind,j)%alpha_ppl,&
                               cexp_ppl=vpplt(kkind,j)%cexp_ppl,&
                               ppl_radius=vpplt(kkind,j)%ppl_radius)
          END IF
          NULLIFY (vpplt(kkind,j)%r2)
          NULLIFY (vpplt(kkind,j)%neighbor)
          NULLIFY (vpplt(kkind,j)%r)
        END DO

      END DO

    END IF

!   *** Load GTH pseudo potential data (non-local part -> PPNL) ***

    ppnl_present = ASSOCIATED(sac_ppnl)

    IF (ppnl_present) THEN

      ALLOCATE (vppnlt(nkind,0:nthread-1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "vppnlt",nkind*int_size)

      DO kkind=1,nkind

        atomic_kind => atomic_kind_set(kkind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             gth_potential=gth_potential)

        DO j=0,nthread-1
          NULLIFY (vppnlt(kkind,j)%alpha_ppnl)
          NULLIFY (vppnlt(kkind,j)%cprj)
          NULLIFY (vppnlt(kkind,j)%neighbor)
          NULLIFY (vppnlt(kkind,j)%nprj_ppnl)
          NULLIFY (vppnlt(kkind,j)%vprj_ppnl)
          IF (ASSOCIATED(gth_potential)) THEN
            CALL get_potential(potential=gth_potential,&
                               alpha_ppnl=vppnlt(kkind,j)%alpha_ppnl,&
                               cprj=vppnlt(kkind,j)%cprj,&
                               lppnl=vppnlt(kkind,j)%lppnl,&
                               nppnl=vppnlt(kkind,j)%nppnl,&
                               nprj_ppnl=vppnlt(kkind,j)%nprj_ppnl,&
                               ppnl_radius=vppnlt(kkind,j)%ppnl_radius,&
                               vprj_ppnl=vppnlt(kkind,j)%vprj_ppnl)
          END IF
        END DO

      END DO

    END IF

    IF (all_potential_present) THEN

      ALLOCATE(vallt(nkind,0:nthread-1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "vallt",nkind*int_size)

      DO kkind=1,nkind

        atomic_kind => atomic_kind_set(kkind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             all_potential=all_potential)

        DO j=0,nthread-1
          NULLIFY(vallt(kkind,j)%neighbor)
          NULLIFY(vallt(kkind,j)%force_c)
          NULLIFY(vallt(kkind,j)%index_c)
          IF (ASSOCIATED(all_potential)) THEN
            CALL get_potential(potential=all_potential,&
              alpha_core_charge=vallt(kkind,j)%alpha_c,&
              zeff=vallt(kkind,j)%zeta_c,&
              ccore_charge=vallt(kkind,j)%core_charge,&
              core_charge_radius=vallt(kkind,j)%core_radius)
          ELSE
            vallt(kkind,j)%alpha_c = 0.0_dp
            vallt(kkind,j)%zeta_c = 0.0_dp
            vallt(kkind,j)%core_charge = 0.0_dp
            vallt(kkind,j)%core_radius = 0.0_dp
          END IF 
        END DO  ! j
      END DO  ! kkind

    END IF


!    *** Retrieve data of SAC_OCE and of SAC_3C neighbors of atom "iatom" ***

    DO ikind = 1,nkind
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           natom=nat,atom_list=atom_list)
      DO iat = 1,nat
        iatom = atom_list(iat)
        IF(oce_present) THEN
           CALL retrieve_sac_list(oce%voce,sac_oce,iatom,ikind,nkind,globenv)
        END IF
        IF(rho0_present) THEN
          CALL retrieve_sac_list(VH_3c_list,sac_3c,iatom,ikind,nkind,globenv)
        ENDIF
      END DO
    END DO

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

        CALL get_neighbor_list_set(neighbor_list_set=&
                                   sab_orb(iab)%neighbor_list_set,&
                                   nlist=nlist)

        NULLIFY ( sab_orb_neighbor_list )

!$OMP parallel &
!$OMP default(none) &
!$OMP private(ithread,hab,sab,sdab,work,pab,ilist,iatom,nnode,atom_a,kkind,iac) &
!$OMP private(sac_ppl_neighbor_list,nneighbor,istat,sac_ppl_neighbor_node,kneighbor) &
!$OMP private(sac_ppnl_neighbor_list,ncol,sac_ppnl_neighbor_node,rac,rac2,dac,nrow) &
!$OMP private(iset,ncoa,prjc,l,nprjc,lc_max,lc_min,zetc,ncoc,i,first_col,last_jatom) &
!$OMP private(sab_orb_neighbor_node,inode,jatom,rab,new_atom_b,atom_b,f0,irow,icol) &
!$OMP private(reduced_3c_oce,reduced_3c_rho0)&
!$OMP private(rab2,dab,sgfa,jset,ncob,sgfb,f,katom,rbc,rbc2,dbc,atom_c,jco,ico,rprjc) &
!$OMP private(ibc,sbc_ppnl_neighbor_list,nnode_c,sbc_ppnl_neighbor_node,knode) &
!$OMP private(cell_c,h_block,p_block,w_block,j,force_a,force_b,neighbor_att) &
!$OMP private(sab_orb_neighbor_list_local,sint,tint,vppl,vppnl,vall,ai_work,oceh,oces) &
!$OMP private(icell,jcell,kcell,idx,n,lppnl,nppnl,cell_c_max,cell_c_min,sac_ppnl_neighbor) &
!$OMP shared(sbc_ppnl,ikind,jkind,force,return_t_derivatives,return_s_derivatives) &
!$OMP shared(lb_min,rpgfb,zetb,sphi_b,sphi_a,nsgfa,nsgfb,set_radius_b,ldai,maxlppl) &
!$OMP shared(sab_orb,first_sgfb,lb_max,npgfb,nsetb,first_sgfa,matrix_w,globenv,matrix_p,matrix_h,matrix_t,iab) &
!$OMP shared(oce,sac_oce,sbc_oce,sac_3c,sbc_3c,VH_3c_list,oce_present)&
!$OMP shared(build_kinetic_energy_matrix,matrix_s,maxder,nder,la_min,rpgfa,zeta) &
!$OMP shared(set_radius_a,nco,ncoset,npgfa,nseta,ncoa_sum,vppnlt,vallt,sac_ppnl,ppnl_present,all_potential_present) &
!$OMP shared(vpplt,sac_ppl,nkind,ppl_present,atom_of_kind,sab_orb_neighbor_list,rho0_present) &
!$OMP shared(nlist,calculate_forces,pabt,workt,sdabt,sabt,habt,la_max,sintt,tintt,oceht,ocest,qs_env)
        ithread = 0
!$      ithread = omp_get_thread_num()
        hab => habt(:,:,ithread)
        sab => sabt(:,:,ithread)
        sdab => sdabt(:,:,:,ithread)
        work => workt(:,:,ithread)
        IF (calculate_forces) pab => pabt(:,:,ithread)
        sint => sintt(:,ithread)
        IF (build_kinetic_energy_matrix) THEN
          tint => tintt(:,ithread)
        END IF
        IF (ppl_present) THEN
          vppl => vpplt(:,ithread)
        END IF
        IF (ppnl_present) THEN
          vppnl => vppnlt(:,ithread)
        END IF
        IF(qs_env%dft_control%qs_control%gapw) THEN
          oceh => oceht(:,ithread)
          oces => ocest(:,ithread)
        END IF
        IF (all_potential_present) THEN
          vall => vallt(:,ithread)
        END IF
        ALLOCATE (ai_work(ldai,ldai,MAX(1,ncoset(maxlppl)),ncoset(nder+1)),STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine_name,module_name,__LINE__,&
                           "ai_work",ldai*ldai*ncoset(maxlppl)*ncoset(nder+1)*dp_size)
        END IF
#if defined(__AIX) && defined(__INIT_WITH_NANS)
        ai_work = nans_dp
#endif

!$OMP do
        DO ilist=1,nlist

!$OMP critical (qs_core_neighbor_list)
          IF ( .NOT. ASSOCIATED(sab_orb_neighbor_list) ) THEN
            sab_orb_neighbor_list => first_list(sab_orb(iab)%neighbor_list_set)
          ELSE
            sab_orb_neighbor_list => next(sab_orb_neighbor_list)
          END IF
          sab_orb_neighbor_list_local => sab_orb_neighbor_list
!$OMP end critical (qs_core_neighbor_list)
          CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list_local,&
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
                                 "vppl(kkind)%r2",nneighbor*dp_size)
              END IF

              ALLOCATE (vppl(kkind)%neighbor(nneighbor),STAT=istat)
              IF (istat /= 0) THEN
                CALL stop_memory(routine_name,module_name,__LINE__,&
                                 "vppl(kkind)%neighbor",nneighbor*int_size)
              END IF

              ALLOCATE (vppl(kkind)%r(3,nneighbor),STAT=istat)
              IF (istat /= 0) THEN
                CALL stop_memory(routine_name,module_name,__LINE__,&
                                 "vppl(kkind)%r",3*nneighbor*dp_size)
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

              lppnl = vppnl(kkind)%lppnl
              nppnl = vppnl(kkind)%nppnl
              rprjc(1) = vppnl(kkind)%ppnl_radius

              sac_ppnl_neighbor_list =>&
                find_neighbor_list(neighbor_list_set=&
                                     sac_ppnl(iac)%neighbor_list_set,&
                                   atom=iatom)

              CALL get_neighbor_list(neighbor_list=sac_ppnl_neighbor_list,&
                                     nnode=nneighbor)

!             *** Find the proper cell index ranges ***

              cell_c_max(:) = 0
              cell_c_min(:) = 0

              sac_ppnl_neighbor_node => first_node(sac_ppnl_neighbor_list)

              DO kneighbor=1,nneighbor
                CALL get_neighbor_node(neighbor_node=sac_ppnl_neighbor_node,&
                                       cell=cell_c)
                DO i=1,3
                  cell_c_max(i) = MAX(cell_c_max(i),cell_c(i))
                  cell_c_min(i) = MIN(cell_c_min(i),cell_c(i))
                END DO
                sac_ppnl_neighbor_node => next(sac_ppnl_neighbor_node)
              END DO

!             *** Allocate and initialize the sac_ppnl neighbor lists ***

              ALLOCATE (vppnl(kkind)%neighbor(cell_c_min(1):cell_c_max(1),&
                                              cell_c_min(2):cell_c_max(2),&
                                              cell_c_min(3):cell_c_max(3)),&
                        STAT=istat)
              IF (istat /= 0 ) THEN
                CALL stop_memory(routine_name,module_name,__LINE__,&
                                 "vppnl(kkind)%neighbor",&
                                 (cell_c_max(1) - cell_c_min(1) + 1)*&
                                 (cell_c_max(2) - cell_c_min(2) + 1)*&
                                 (cell_c_max(3) - cell_c_min(3) + 1)*&
                                 int_size)
              END IF
              vppnl(kkind)%neighbor(:,:,:)%n = 0

!             *** Find proper array size for each cell ***

              sac_ppnl_neighbor_node => first_node(sac_ppnl_neighbor_list)

              DO kneighbor=1,nneighbor
                CALL get_neighbor_node(neighbor_node=sac_ppnl_neighbor_node,&
                                       cell=cell_c)
                vppnl(kkind)%neighbor(cell_c(1),cell_c(2),cell_c(3))%n =&
                  vppnl(kkind)%neighbor(cell_c(1),cell_c(2),cell_c(3))%n + 1
                sac_ppnl_neighbor_node => next(sac_ppnl_neighbor_node)
              END DO

!             *** Allocate the sac_ppnl neighbor list ***

              DO kcell=cell_c_min(3),cell_c_max(3)
                DO jcell=cell_c_min(2),cell_c_max(2)
                  DO icell=cell_c_min(1),cell_c_max(1)
                    sac_ppnl_neighbor => vppnl(kkind)%neighbor(icell,jcell,kcell)
                    NULLIFY (sac_ppnl_neighbor%sac)
                    NULLIFY (sac_ppnl_neighbor%r2)
                    NULLIFY (sac_ppnl_neighbor%r)
                    NULLIFY (sac_ppnl_neighbor%index_list)
                    NULLIFY (sac_ppnl_neighbor%list)
                    n = sac_ppnl_neighbor%n
                    IF (n > 0) THEN
                      CALL reallocate(sac_ppnl_neighbor%r2,1,n)
                      CALL reallocate(sac_ppnl_neighbor%r,1,3,1,n)
                      CALL reallocate(sac_ppnl_neighbor%index_list,1,n)
                      CALL reallocate(sac_ppnl_neighbor%list,1,n)
                      ALLOCATE (sac_ppnl_neighbor%sac(maxder))
                      IF (istat /= 0) THEN
                        CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "sac_ppnl_neighbor%sac",maxder*int_size)
                      END IF
                      DO i=1,maxder
                        ncol = n*nppnl
                        ALLOCATE (sac_ppnl_neighbor%sac(i)%block(ncoa_sum,ncol),&
                                  STAT=istat)
                        IF (istat /= 0) THEN
                          CALL stop_memory(routine_name,module_name,__LINE__,&
                                           "sac_ppnl_neighbor%sac(i)%block",&
                                           ncoa_sum*ncol*dp_size)
                        END IF
#if defined(__AIX) && defined(__INIT_WITH_NANS)
                        sac_ppnl_neighbor%sac(i)%block = nans_dp
#endif
                      END DO
                    END IF
                  END DO
                END DO
              END DO

!             *** Fill sac_ppnl neighbor lists ***

              vppnl(kkind)%neighbor(:,:,:)%n = 0

              sac_ppnl_neighbor_node => first_node(sac_ppnl_neighbor_list)

              DO kneighbor=1,nneighbor
                CALL get_neighbor_node(neighbor_node=sac_ppnl_neighbor_node,&
                                       neighbor=katom,&
                                       cell=cell_c,&
                                       r=rac)
                sac_ppnl_neighbor => vppnl(kkind)%neighbor(cell_c(1),&
                                                           cell_c(2),&
                                                           cell_c(3))
                sac_ppnl_neighbor%n = sac_ppnl_neighbor%n + 1
                sac_ppnl_neighbor%list(sac_ppnl_neighbor%n) = katom
                rac2 = rac(1)*rac(1) + rac(2)*rac(2) + rac(3)*rac(3)
                sac_ppnl_neighbor%r2(sac_ppnl_neighbor%n) = rac2
                sac_ppnl_neighbor%r(:,sac_ppnl_neighbor%n) = rac(:)
                sac_ppnl_neighbor_node => next(sac_ppnl_neighbor_node)
              END DO

!             *** Build sac overlap integral blocks ***

              DO kcell=cell_c_min(3),cell_c_max(3)
                DO jcell=cell_c_min(2),cell_c_max(2)
                  DO icell=cell_c_min(1),cell_c_max(1)

                    sac_ppnl_neighbor => vppnl(kkind)%neighbor(icell,jcell,kcell)

!                   *** Sort sac_ppnl neighbor lists ***

                    IF (sac_ppnl_neighbor%n > 0) THEN
                      CALL sort(sac_ppnl_neighbor%list,&
                                sac_ppnl_neighbor%n,&
                                sac_ppnl_neighbor%index_list)
                    ELSE
                      CYCLE
                    END IF

                    ncol = 1

                    DO kneighbor=1,sac_ppnl_neighbor%n

                      idx = sac_ppnl_neighbor%index_list(kneighbor)

                      rac(:) = sac_ppnl_neighbor%r(:,idx)
                      rac2 = sac_ppnl_neighbor%r2(idx)
                      dac = SQRT(rac2)

                      nrow = 1

                      DO iset=1,nseta

                        ncoa = npgfa(iset)*ncoset(la_max(iset))

                        prjc = 1

                        DO l=0,lppnl

                          nprjc = vppnl(kkind)%nprj_ppnl(l)*nco(l)

                          IF (nprjc == 0) CYCLE

                          IF (set_radius_a(iset) + rprjc(1) < dac) CYCLE

                          lc_max = l + 2*(vppnl(kkind)%nprj_ppnl(l) - 1)
                          lc_min = l
                          zetc(1) = vppnl(kkind)%alpha_ppnl(l)
                          ncoc = ncoset(lc_max)

!                         *** Calculate the primitive overlap integrals ***

#if defined(__AIX) && defined(__INIT_WITH_NANS)
                          ai_work = nans_dp
                          sab = nans_dp
#endif
                          CALL overlap(la_max(iset),la_min(iset),npgfa(iset),&
                                       rpgfa(:,iset),zeta(:,iset),&
                                       lc_max,lc_min,1,rprjc,zetc,&
                                       rac,dac,sab,nder,.TRUE.,ai_work,ldai)

!                         *** Contraction step (projector functions) ***

                          DO i=1,maxder
                            first_col = (i - 1)*SIZE(work,1)
                            CALL dgemm("N","N",ncoa,nprjc,ncoc,&
                                       1.0_dp,sab(1,first_col+1),SIZE(sab,1),&
                                       vppnl(kkind)%cprj(1,prjc),&
                                       SIZE(vppnl(kkind)%cprj,1),0.0_dp,&
                                       work(1,first_col+prjc),SIZE(work,1))
                          END DO

                          prjc = prjc + nprjc

                        END DO

                        DO i=1,maxder
                          first_col = (i - 1)*SIZE(work,1) + 1
                          CALL dgemm("N","N",ncoa,nppnl,nppnl,&
                                     1.0_dp,work(1,first_col),SIZE(work,1),&
                                     vppnl(kkind)%vprj_ppnl(1,1),&
                                     SIZE(vppnl(kkind)%vprj_ppnl,1),0.0_dp,&
                                     sac_ppnl_neighbor%sac(i)%block(nrow,ncol),&
                                     SIZE(sac_ppnl_neighbor%sac(i)%block,1))
                        END DO

                        nrow = nrow + ncoa

                      END DO

                      ncol = ncol + nppnl

                    END DO ! kneighbor

                  END DO ! icell
                END DO ! jcell
              END DO ! kcell

            END DO ! kkind

          END IF

          IF(all_potential_present) THEN
            DO kkind = 1,nkind
              IF(ASSOCIATED(VH_3c_list(kkind,iatom)%neighbor)) &
                 vall(kkind)%neighbor =>  VH_3c_list(kkind,iatom)%neighbor
            END DO !  kkind
          END IF

          last_jatom = 0

          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)

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

            IF(oce_present) THEN
              CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                    reduced_3c_oce=reduced_3c_oce)
!            *** Build reduced 3-center lists
              CALL build_reduced_3c_lists(atomic_kind_set, iatom, jatom, jkind,&
                                          sbc_oce, oce%voce, reduced_3c_oce)

!            *** Calculate the coefficients of the one center expansions ***
              CALL set_up_oce(atomic_kind_set,reduced_3c_oce,oce,&
                              ikind,iatom,jkind,jatom,rab,&
                              oceh,oces,nder)
            END IF

            IF(rho0_present) THEN
              CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                    reduced_3c_rho0=reduced_3c_rho0)
!            *** Build reduced 3-center lists
              CALL build_reduced_3c_lists(atomic_kind_set, iatom, jatom, jkind,&
                                          sbc_3c, VH_3c_list, reduced_3c_rho0)
            END IF

            IF (calculate_forces) atom_b = atom_of_kind(jatom)

!           *** Use the symmetry of the first derivatives ***

            IF (iatom == jatom) THEN
              f0 = 1.0_dp
            ELSE
              f0 = 2.0_dp
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
              DO i=1,SIZE(matrix_s)
                NULLIFY (sint(i)%block)
                CALL add_block_node(matrix=matrix_s(i)%matrix,&
                                    block_row=irow,&
                                    block_col=icol,&
                                    block=sint(i)%block)
              END DO
              IF (build_kinetic_energy_matrix) THEN
                DO i=1,SIZE(matrix_t)
                  NULLIFY (tint(i)%block)
                  CALL add_block_node(matrix=matrix_t(i)%matrix,&
                                      block_row=irow,&
                                      block_col=icol,&
                                      block=tint(i)%block)
                END DO
              END IF
              NULLIFY (h_block)
              CALL add_block_node(matrix=matrix_h(1)%matrix,&
                                  block_row=irow,&
                                  block_col=icol,&
                                  block=h_block)
              IF (calculate_forces) THEN
                CALL get_block_node(matrix=matrix_p(1)%matrix,&
                                    block_row=irow,&
                                    block_col=icol,&
                                    block=p_block)
                IF (.NOT.ASSOCIATED(p_block)) THEN
                  CALL stop_program(routine_name,module_name,__LINE__,&
                                    "A density matrix block is missing",&
                                    globenv)
                END IF
                CALL get_block_node(matrix=matrix_w(1)%matrix,&
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
                                 1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 w_block(sgfa,sgfb),SIZE(w_block,1),&
                                 0.0_dp,work(1,1),SIZE(work,1))
                    ELSE
                      CALL dgemm("N","T",ncoa,nsgfb(jset),nsgfa(iset),&
                                 1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 w_block(sgfb,sgfa),SIZE(w_block,1),&
                                 0.0_dp,work(1,1),SIZE(work,1))
                    END IF
                    CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
                               1.0_dp,work(1,1),SIZE(work,1),&
                               sphi_b(1,sgfb),SIZE(sphi_b,1),&
                               0.0_dp,pab(1,1),SIZE(pab,1))

!                   *** Calculate the primitive overlap integrals ***
!                   *** and the corresponding force contribution  ***

#if defined(__AIX) && defined(__INIT_WITH_NANS)
                    ai_work = nans_dp
                    sab = nans_dp
                    sdab = nans_dp
#endif
                    CALL overlap(la_max(iset),la_min(iset),npgfa(iset),&
                                 rpgfa(:,iset),zeta(:,iset),&
                                 lb_max(jset),lb_min(jset),npgfb(jset),&
                                 rpgfb(:,jset),zetb(:,jset),&
                                 rab,dab,sab,nder,return_s_derivatives,&
                                 ai_work,ldai,sdab,pab,force_a)
!$OMP critical (qs_core_force_overlap)
                    force(ikind)%overlap(:,atom_a) =&
                      force(ikind)%overlap(:,atom_a) - 2.0_dp*force_a(:)
                    force(jkind)%overlap(:,atom_b) =&
                      force(jkind)%overlap(:,atom_b) + 2.0_dp*force_a(:)
!$OMP end critical (qs_core_force_overlap)

                  ELSE

!                   *** Calculate the primitive overlap integrals ***

#if defined(__AIX) && defined(__INIT_WITH_NANS)
                    ai_work = nans_dp
                    sab = nans_dp
                    sdab = nans_dp
#endif
                    CALL overlap(la_max(iset),la_min(iset),npgfa(iset),&
                                 rpgfa(:,iset),zeta(:,iset),&
                                 lb_max(jset),lb_min(jset),npgfb(jset),&
                                 rpgfb(:,jset),zetb(:,jset),&
                                 rab,dab,sab,nder,return_s_derivatives,&
                                 ai_work,ldai,sdab)

                  END IF

!                 *** Contraction step (overlap matrix and its derivatives) ***

                  DO i=1,SIZE(matrix_s)
                    IF ((i > 1).AND.(.NOT.return_s_derivatives)) CYCLE
                    first_col = (i - 1)*SIZE(sab,1) + 1
                    CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                               1.0_dp,sab(1,first_col),SIZE(sab,1),&
                               sphi_b(1,sgfb),SIZE(sphi_b,1),&
                               0.0_dp,work(1,1),SIZE(work,1))
                    IF (iatom <= jatom) THEN
                      CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                                 1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 work(1,1),SIZE(work,1),&
                                 1.0_dp,sint(i)%block(sgfa,sgfb),&
                                 SIZE(sint(i)%block,1))
                    ELSE
                      IF (i > 1) THEN
                        f = -1.0_dp
                      ELSE
                        f = 1.0_dp
                      END IF
                      CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                                 f,work(1,1),SIZE(work,1),&
                                 sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 1.0_dp,sint(i)%block(sgfb,sgfa),&
                                 SIZE(sint(i)%block,1))
                    END IF
                  END DO

!                 *** Decontract density matrix block ***

                  IF (calculate_forces) THEN
                    IF (iatom <= jatom) THEN
                      CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
                                 1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 p_block(sgfa,sgfb),SIZE(p_block,1),&
                                 0.0_dp,work(1,1),SIZE(work,1))
                    ELSE
                      CALL dgemm("N","T",ncoa,nsgfb(jset),nsgfa(iset),&
                                 1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 p_block(sgfb,sgfa),SIZE(p_block,1),&
                                 0.0_dp,work(1,1),SIZE(work,1))
                    END IF
                    CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
                               1.0_dp,work(1,1),SIZE(work,1),&
                               sphi_b(1,sgfb),SIZE(sphi_b,1),&
                               0.0_dp,pab(1,1),SIZE(pab,1))
                  END IF

                  IF (calculate_forces.AND.(iatom /= jatom)) THEN

!                   *** Calculate the primitive kinetic energy integrals ***
!                   *** and the corresponding force contribution         ***

#if defined(__AIX) && defined(__INIT_WITH_NANS)
                    ai_work = nans_dp
                    hab = nans_dp
#endif
                    CALL kinetic(la_max(iset),la_min(iset),npgfa(iset),&
                                 rpgfa(:,iset),zeta(:,iset),&
                                 lb_max(jset),lb_min(jset),npgfb(jset),&
                                 rpgfb(:,jset),zetb(:,jset),&
                                 rab,dab,sdab,hab,nder,return_t_derivatives,&
                                 ai_work,ldai,pab,force_a)

!$OMP critical (qs_core_force_kinetic)
                    force(ikind)%kinetic(:,atom_a) =&
                      force(ikind)%kinetic(:,atom_a) + 2.0_dp*force_a(:)
                    force(jkind)%kinetic(:,atom_b) =&
                      force(jkind)%kinetic(:,atom_b) - 2.0_dp*force_a(:)
!$OMP end critical (qs_core_force_kinetic)

                  ELSE

!                   *** Calculate the primitive kinetic energy integrals ***

#if defined(__AIX) && defined(__INIT_WITH_NANS)
                    ai_work = nans_dp
                    hab = nans_dp
#endif
                    CALL kinetic(la_max(iset),la_min(iset),npgfa(iset),&
                                 rpgfa(:,iset),zeta(:,iset),&
                                 lb_max(jset),lb_min(jset),npgfb(jset),&
                                 rpgfb(:,jset),zetb(:,jset),&
                                 rab,dab,sdab,hab,nder,return_t_derivatives,&
                                 ai_work,ldai)

                  END IF

!                 *** Contraction step (kinetic energy integral ***
!                 *** matrix and its derivatives)               ***

                  IF (build_kinetic_energy_matrix) THEN
                    DO i=1,SIZE(matrix_t)
                      IF ((i > 1).AND.(.NOT.return_t_derivatives)) CYCLE
                      first_col = (i - 1)*SIZE(hab,1) + 1
                      CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                                 1.0_dp,hab(1,first_col),SIZE(hab,1),&
                                 sphi_b(1,sgfb),SIZE(sphi_b,1),&
                                 0.0_dp,work(1,1),SIZE(work,1))
                      IF (iatom <= jatom) THEN
                        CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                                   1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                   work(1,1),SIZE(work,1),&
                                   1.0_dp,tint(i)%block(sgfa,sgfb),&
                                   SIZE(tint(i)%block,1))
                      ELSE
                        IF (i > 1) THEN
                          f = -1.0_dp
                        ELSE
                          f = 1.0_dp
                        END IF
                        CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                                   f,work(1,1),SIZE(work,1),&
                                   sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                   1.0_dp,tint(i)%block(sgfb,sgfa),&
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

#if defined(__AIX) && defined(__INIT_WITH_NANS)
                    ai_work = nans_dp
#endif
                          CALL overlap_ppl(&
                            la_max(iset),la_min(iset),npgfa(iset),&
                            rpgfa(:,iset),zeta(:,iset),&
                            lb_max(jset),lb_min(jset),npgfb(jset),&
                            rpgfb(:,jset),zetb(:,jset),&
                            vppl(kkind)%cexp_ppl,vppl(kkind)%alpha_ppl,&
                            vppl(kkind)%ppl_radius,&
                            rab,dab,rac,dac,rbc,dbc,&
                            hab,nder,nder,.FALSE.,ai_work,&
                            pab,force_a,force_b)

!                         *** The derivatives w.r.t. atomic center c are    ***
!                         *** calculated using the translational invariance ***
!                         *** of the first derivatives                      ***

                          atom_c = atom_of_kind(katom)

!$OMP critical (qs_core_force_gth_ppl)
                          force(ikind)%gth_ppl(:,atom_a) =&
                            force(ikind)%gth_ppl(:,atom_a) + f0*force_a(:)
                          force(kkind)%gth_ppl(:,atom_c) =&
                            force(kkind)%gth_ppl(:,atom_c) - f0*force_a(:)

                          force(jkind)%gth_ppl(:,atom_b) =&
                            force(jkind)%gth_ppl(:,atom_b) + f0*force_b(:)
                          force(kkind)%gth_ppl(:,atom_c) =&
                            force(kkind)%gth_ppl(:,atom_c) - f0*force_b(:)
!$OMP end critical (qs_core_force_gth_ppl)

                        ELSE

#if defined(__AIX) && defined(__INIT_WITH_NANS)
                    ai_work = nans_dp
#endif
                          CALL overlap_ppl(&
                            la_max(iset),la_min(iset),npgfa(iset),&
                            rpgfa(:,iset),zeta(:,iset),&
                            lb_max(jset),lb_min(jset),npgfb(jset),&
                            rpgfb(:,jset),zetb(:,jset),&
                            vppl(kkind)%cexp_ppl,vppl(kkind)%alpha_ppl,&
                            vppl(kkind)%ppl_radius,&
                            rab,dab,rac,dac,rbc,dbc,&
                            hab,0,0,.FALSE.,ai_work)

                        END IF

                      END DO

                    END DO

                  END IF

                ELSE

                  DO jco=1,ncob
                    DO ico=1,ncoa
                      hab(ico,jco) = 0.0_dp
                    END DO
                  END DO

                  IF (calculate_forces) THEN

!                   *** Decontract density matrix block ***

                    IF (iatom <= jatom) THEN
                      CALL dgemm("N","N",ncoa,nsgfb(jset),nsgfa(iset),&
                                 1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 p_block(sgfa,sgfb),SIZE(p_block,1),&
                                 0.0_dp,work(1,1),SIZE(work,1))
                    ELSE
                      CALL dgemm("N","T",ncoa,nsgfb(jset),nsgfa(iset),&
                                 1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                                 p_block(sgfb,sgfa),SIZE(p_block,1),&
                                 0.0_dp,work(1,1),SIZE(work,1))
                    END IF
                    CALL dgemm("N","T",ncoa,ncob,nsgfb(jset),&
                               1.0_dp,work(1,1),SIZE(work,1),&
                               sphi_b(1,sgfb),SIZE(sphi_b,1),&
                               0.0_dp,pab(1,1),SIZE(pab,1))

                  END IF

                END IF

!               *** Loop over the GTH pseudo potential atoms (non-local part) ***
                
                IF (ppnl_present.AND.new_atom_b) THEN
                   
                   DO kkind=1,nkind
                      
                      IF (.NOT.ASSOCIATED(vppnl(kkind)%neighbor)) CYCLE
                      
                      ibc = jkind + nkind*(kkind - 1)
                      
                      IF (.NOT.ASSOCIATED(sbc_ppnl(ibc)%neighbor_list_set)) CYCLE
                      
                      lppnl = vppnl(kkind)%lppnl
                      nppnl = vppnl(kkind)%nppnl
                      rprjc(1) = vppnl(kkind)%ppnl_radius
                      
                      sbc_ppnl_neighbor_list =>&
                           find_neighbor_list(neighbor_list_set=&
                           sbc_ppnl(ibc)%neighbor_list_set,&
                           atom=jatom)
                      
                      CALL get_neighbor_list(neighbor_list=sbc_ppnl_neighbor_list,&
                           nnode=nnode_c)
                      
                      sbc_ppnl_neighbor_node => first_node(sbc_ppnl_neighbor_list)
                      
                      knodeloop: DO knode=1,nnode_c
                         
                         CALL get_neighbor_node(neighbor_node=sbc_ppnl_neighbor_node,&
                              neighbor=katom,&
                              cell=cell_c,&
                              r=rbc)
                         
                         ! The < and > operators don't work correctly on arrays in
                         ! the openmp version
                         neighbor_att => vppnl(kkind)%neighbor
                         DO i=1, 3
                            IF ( (cell_c(i) < LBOUND(neighbor_att,i)) .OR. &
                                 (cell_c(i) > UBOUND(neighbor_att,i))) THEN
                               sbc_ppnl_neighbor_node => next(sbc_ppnl_neighbor_node)
                               CYCLE knodeloop
                            END IF
                         END DO

                      sac_ppnl_neighbor => vppnl(kkind)%neighbor(cell_c(1),&
                                                                 cell_c(2),&
                                                                 cell_c(3))
                      

                      IF (sac_ppnl_neighbor%n == 0) THEN
                        sbc_ppnl_neighbor_node => next(sbc_ppnl_neighbor_node)
                        CYCLE
                      END IF

!                     *** Locate operator atom in the sac_ppnl neighbor list ***
                      kneighbor = locate(sac_ppnl_neighbor%list,katom)

                      IF (kneighbor == 0) THEN
                        sbc_ppnl_neighbor_node => next(sbc_ppnl_neighbor_node)
                        CYCLE
                      END IF

                      idx = sac_ppnl_neighbor%index_list(kneighbor)
                      rac(:) = sac_ppnl_neighbor%r(:,idx)

                      rac2 = sac_ppnl_neighbor%r2(idx)
                      dac = SQRT(rac2)
                      IF (set_radius_a(iset) + rprjc(1) < dac) THEN
                        sbc_ppnl_neighbor_node => next(sbc_ppnl_neighbor_node)
                        CYCLE
                      END IF

                      rbc2 = rbc(1)*rbc(1) + rbc(2)*rbc(2) + rbc(3)*rbc(3)
                      dbc = SQRT(rbc2)
                      IF (set_radius_b(jset) + rprjc(1) < dbc) THEN
                        sbc_ppnl_neighbor_node => next(sbc_ppnl_neighbor_node)
                        CYCLE
                      END IF

                      prjc = 1

                      DO l=0,lppnl

                        nprjc = vppnl(kkind)%nprj_ppnl(l)*nco(l)

                        IF (nprjc == 0) CYCLE

                        lc_max = l + 2*(vppnl(kkind)%nprj_ppnl(l) - 1)
                        lc_min = l
                        zetc(1) = vppnl(kkind)%alpha_ppnl(l)
                        ncoc = ncoset(lc_max)

!                       *** Calculate the primitive overlap integrals ***

#if defined(__AIX) && defined(__INIT_WITH_NANS)
                        ai_work = nans_dp
                        sab = nans_dp
#endif
                        CALL overlap(lb_max(jset),lb_min(jset),npgfb(jset),&
                                     rpgfb(:,jset),zetb(:,jset),&
                                     lc_max,lc_min,1,rprjc,zetc,&
                                     rbc,dbc,sab,nder,.TRUE.,ai_work,ldai)

!                       *** Contraction step (projector functions) ***

                        DO i=1,maxder
                          first_col = (i - 1)*SIZE(work,1)
                          CALL dgemm("N","N",ncob,nprjc,ncoc,1.0_dp,&
                                     sab(1,first_col+1),SIZE(sab,1),&
                                     vppnl(kkind)%cprj(1,prjc),&
                                     SIZE(vppnl(kkind)%cprj,1),&
                                     0.0_dp,work(1,first_col+prjc),&
                                     SIZE(work,1))
                        END DO

                        prjc = prjc + nprjc

                      END DO

                      ncol = (kneighbor - 1)*nppnl + 1

                      CALL dgemm("N","T",ncoa,ncob,nppnl,1.0_dp,&
                                 sac_ppnl_neighbor%sac(1)%block(nrow,ncol),&
                                 SIZE(sac_ppnl_neighbor%sac(1)%block,1),&
                                 work(1,1),SIZE(work,1),&
                                 1.0_dp,hab(1,1),SIZE(hab,1))

                      IF (calculate_forces) THEN

                        atom_c = atom_of_kind(katom)

                        DO i=1,3

                          first_col = i*SIZE(work,1) + 1

                          CALL dgemm("N","T",ncoa,ncob,nppnl,1.0_dp,&
                                     sac_ppnl_neighbor%sac(i+1)%block(nrow,ncol),&
                                     SIZE(sac_ppnl_neighbor%sac(i+1)%block,1),&
                                     work(1,1),SIZE(work,1),&
                                     0.0_dp,hab(1,first_col),SIZE(hab,1))

                          f = 0.0_dp

                          DO jco=1,ncob
                            j = first_col + jco - 1
                            DO ico=1,ncoa
                              f = f + pab(ico,jco)*hab(ico,j)
                            END DO
                          END DO

!$OMP critical (qs_core_force_gth_ppnl)
                          force(ikind)%gth_ppnl(i,atom_a) =&
                            force(ikind)%gth_ppnl(i,atom_a) + f0*f

                          force(kkind)%gth_ppnl(i,atom_c) =&
                            force(kkind)%gth_ppnl(i,atom_c) - f0*f
!$OMP end critical (qs_core_force_gth_ppnl)

                          CALL dgemm("N","T",ncoa,ncob,nppnl,1.0_dp,&
                                     sac_ppnl_neighbor%sac(1)%block(nrow,ncol),&
                                     SIZE(sac_ppnl_neighbor%sac(1)%block,1),&
                                     work(1,first_col),SIZE(work,1),&
                                     0.0_dp,hab(1,first_col),SIZE(hab,1))

                          f = 0.0_dp

                          DO jco=1,ncob
                            j = first_col + jco - 1
                            DO ico=1,ncoa
                              f = f + pab(ico,jco)*hab(ico,j)
                            END DO
                          END DO

!$OMP critical (qs_core_force_gth_ppnl)
                          force(jkind)%gth_ppnl(i,atom_b) =&
                            force(jkind)%gth_ppnl(i,atom_b) + f0*f

                          force(kkind)%gth_ppnl(i,atom_c) =&
                            force(kkind)%gth_ppnl(i,atom_c) - f0*f
!$OMP end critical (qs_core_force_gth_ppnl)

                        END DO

                      END IF

                      sbc_ppnl_neighbor_node => next(sbc_ppnl_neighbor_node)

                    END DO knodeloop

                  END DO ! kkind

                END IF ! ppnl

                IF (all_potential_present ) THEN

                  IF(calculate_forces) THEN

                    CALL all_integrals(la_max(iset),la_min(iset),npgfa(iset),&
                                 rpgfa(:,iset),zeta(:,iset),&
                                 lb_max(jset),lb_min(jset),npgfb(jset),&
                                 rpgfb(:,jset),zetb(:,jset),&
                                 rab,dab,jkind,jatom,nkind,sbc_3c,vall,&
                                 hab,nder,pab,force_a,force_b)

!$OMP critical (qs_core_force_all_potential)
                    force(ikind)%all_potential(:,atom_a) = &
                        force(ikind)%all_potential(:,atom_a) + f0*force_a(:)
                    force(jkind)%all_potential(:,atom_b) = & 
                        force(jkind)%all_potential(:,atom_b) + f0*force_b(:)
!$OMP end critical (qs_core_force_all_potential)

                    DO kkind = 1,nkind
                      IF (.NOT. ASSOCIATED(vall(kkind)%force_c)) CYCLE
                      nnode_c = SIZE(vall(kkind)%index_c,1)
                      DO knode = 1,nnode_c
                        katom = vall(kkind)%index_c(knode)
                        atom_c = atom_of_kind(katom)
!$OMP critical (qs_core_force_all_potential)
                        force(kkind)%all_potential(1:3,atom_c) = & 
                          force(kkind)%all_potential(1:3,atom_c) + &
                          f0*vall(kkind)%force_c(1:3,knode)
!$OMP end critical (qs_core_force_all_potential)
                      END DO  ! knode
                      DEALLOCATE(vall(kkind)%index_c,vall(kkind)%force_c,STAT=istat)
                    END DO  ! kkind

                  ELSE

                    CALL all_integrals(la_max(iset),la_min(iset),npgfa(iset),&
                                 rpgfa(:,iset),zeta(:,iset),&
                                 lb_max(jset),lb_min(jset),npgfb(jset),&
                                 rpgfb(:,jset),zetb(:,jset),&
                                 rab,dab,jkind,jatom,nkind,sbc_3c,vall,&
                                 hab,nder)
                    DO kkind = 1,nkind
                      IF (.NOT. ASSOCIATED(vall(kkind)%force_c)) CYCLE
                      DEALLOCATE(vall(kkind)%index_c,vall(kkind)%force_c,STAT=istat)
                    END DO  ! kkind
                  END IF

                END IF  ! all_potential

!               *** Contraction step (core Hamiltonian matrix) ***

                CALL dgemm("N","N",ncoa,nsgfb(jset),ncob,&
                           1.0_dp,hab(1,1),SIZE(hab,1),&
                           sphi_b(1,sgfb),SIZE(sphi_b,1),&
                           0.0_dp,work(1,1),SIZE(work,1))

                IF (iatom <= jatom) THEN
                  CALL dgemm("T","N",nsgfa(iset),nsgfb(jset),ncoa,&
                             1.0_dp,sphi_a(1,sgfa),SIZE(sphi_a,1),&
                             work(1,1),SIZE(work,1),&
                             1.0_dp,h_block(sgfa,sgfb),SIZE(h_block,1))
                ELSE
                  CALL dgemm("T","N",nsgfb(jset),nsgfa(iset),ncoa,&
                             1.0_dp,work(1,1),SIZE(work,1),&
                             sphi_a(1,sgfa),SIZE(sphi_a,1),&
                             1.0_dp,h_block(sgfb,sgfa),SIZE(h_block,1))
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
              IF (ASSOCIATED(vppnl(kkind)%neighbor)) THEN
                DO kcell=LBOUND(vppnl(kkind)%neighbor,3),&
                         UBOUND(vppnl(kkind)%neighbor,3)
                  DO jcell=LBOUND(vppnl(kkind)%neighbor,2),&
                           UBOUND(vppnl(kkind)%neighbor,2)
                    DO icell=LBOUND(vppnl(kkind)%neighbor,1),&
                             UBOUND(vppnl(kkind)%neighbor,1)
                      sac_ppnl_neighbor => vppnl(kkind)%neighbor(icell,&
                                                                 jcell,&
                                                                 kcell)
                      IF (sac_ppnl_neighbor%n > 0) THEN
                        DEALLOCATE (sac_ppnl_neighbor%r2,STAT=istat)
                        IF (istat /= 0) THEN
                          CALL stop_memory(routine_name,module_name,__LINE__,&
                                           "sac_ppnl_neighbor%r2")
                        END IF
                        DEALLOCATE (sac_ppnl_neighbor%r,STAT=istat)
                        IF (istat /= 0) THEN
                          CALL stop_memory(routine_name,module_name,__LINE__,&
                                           "sac_ppnl_neighbor%r")
                        END IF
                        DEALLOCATE (sac_ppnl_neighbor%list,STAT=istat)
                        IF (istat /= 0) THEN
                          CALL stop_memory(routine_name,module_name,__LINE__,&
                                           "sac_ppnl_neighbor%list")
                        END IF
                        DEALLOCATE (sac_ppnl_neighbor%index_list,STAT=istat)
                        IF (istat /= 0) THEN
                          CALL stop_memory(routine_name,module_name,__LINE__,&
                                           "sac_ppnl_neighbor%index_list")
                        END IF
                        DO i=1,maxder
                          DEALLOCATE (sac_ppnl_neighbor%sac(i)%block,STAT=istat)
                          IF (istat /= 0) THEN
                            CALL stop_memory(routine_name,module_name,__LINE__,&
                                             "sac_ppnl_neighbor%sac(i)%block")
                          END IF
                        END DO
                        DEALLOCATE (sac_ppnl_neighbor%sac,STAT=istat)
                        IF (istat /= 0) THEN
                          CALL stop_memory(routine_name,module_name,__LINE__,&
                                           "sac_ppnl_neighbor%sac")
                        END IF
                      END IF
                    END DO
                  END DO
                END DO
                DEALLOCATE (vppnl(kkind)%neighbor,STAT=istat)
                IF (istat /= 0) THEN
                  CALL stop_memory(routine_name,module_name,__LINE__,&
                                   "vppnl(kkind)%neighbor")
                END IF
              END IF
            END DO
          END IF

          IF (all_potential_present) THEN
            DO kkind=1,nkind
              IF (ASSOCIATED(vall(kkind)%neighbor)) THEN
                NULLIFY (vall(kkind)%neighbor)
                NULLIFY (vall(kkind)%force_c)
                NULLIFY (vall(kkind)%index_c)
              END IF
            END DO 
          END IF

        END DO ! iatom => atom A

        DEALLOCATE (ai_work,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                         "ai_work")

!$OMP end parallel

      END DO ! jkind

    END DO ! ikind
!   *** Put the core Hamiltonian matrix in the QS environment ***

    CALL set_qs_env(qs_env=qs_env,&
                    matrix_s=matrix_s,&
                    kinetic=matrix_t,&
                    matrix_h=matrix_h,&
                    oce=oce,&
                    VH_3c_list=VH_3c_list)
!   *** Release work storage ***

    IF (ppl_present) THEN
      DO kkind=1,nkind
        DO j=0,nthread-1
          NULLIFY (vpplt(kkind,j)%cexp_ppl)
        END DO
      END DO
      DEALLOCATE (vpplt,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "vpplt")
    END IF

    IF (ppnl_present) THEN
      DO kkind=1,nkind
        DO j=0,nthread-1
          NULLIFY (vppnlt(kkind,j)%alpha_ppnl)
          NULLIFY (vppnlt(kkind,j)%cprj)
          NULLIFY (vppnlt(kkind,j)%nprj_ppnl)
          NULLIFY (vppnlt(kkind,j)%vprj_ppnl)
        END DO
      END DO
      DEALLOCATE (vppnlt,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "vppnlt")
    END IF

    IF (all_potential_present) THEN
      DO kkind=1,nkind
        DO j =0,nthread-1
          NULLIFY (vallt(kkind,j)%neighbor)
          NULLIFY (vallt(kkind,j)%force_c)
          NULLIFY (vallt(kkind,j)%index_c)
        END DO
      END DO
      DEALLOCATE (vallt,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "vallt")
    END IF

    DEALLOCATE (habt,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "hab")
    DEALLOCATE (sabt,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "sab")
    DEALLOCATE (sdabt,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "sdab")
    DEALLOCATE (workt,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "work")
    IF (calculate_forces) THEN

      DEALLOCATE (atom_of_kind,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "atom_of_kind")
      DEALLOCATE (pabt,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "pab")
!     *** If LSD, then recover alpha density and beta density     ***
!     *** from the total density (1) and the spin density (2)     ***
!     *** The W matrix is neglected, since it will be destroyed   ***
!     *** in the calling force routine after leaving this routine ***

      IF (SIZE(matrix_p) == 2) THEN
        CALL cp_sm_scale_and_add(matrix_p(1)%matrix, 0.5_dp,matrix_p(2)%matrix,0.5_dp)
        CALL cp_sm_scale_and_add(matrix_p(2)%matrix,-1.0_dp,matrix_p(1)%matrix,1.0_dp)
      END IF

    END IF

    DO i=1,SIZE(sintt,1)
      DO j=0,nthread-1
        NULLIFY (sintt(i,j)%block)
      END DO
    END DO
    DEALLOCATE (sintt,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "sintt")

    IF (build_kinetic_energy_matrix) THEN
      DO i=1,SIZE(tintt)
        DO j=0,nthread-1
          NULLIFY (tintt(i,j)%block)
        END DO
      END DO
      DEALLOCATE (tintt,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "tintt")
    END IF

    IF(qs_env%dft_control%qs_control%gapw) THEN
       DO i=1,SIZE(oceht,1)
         DO j=0,nthread-1
           NULLIFY (oceht(i,j)%block)
           NULLIFY (ocest(i,j)%block)
         END DO
       END DO
       DEALLOCATE (oceht,STAT=istat)
       IF (istat /= 0) THEN
         CALL stop_memory(routine_name,module_name,__LINE__,"oceht")
       END IF
       DEALLOCATE (ocest,STAT=istat)
       IF (istat /= 0) THEN
         CALL stop_memory(routine_name,module_name,__LINE__,"ocest")
       END IF
    ENDIF

    CALL timestop(0.0_dp,handle)

!   *** Print the distribution of the overlap matrix blocks ***

    CALL write_matrix_distribution(matrix_s(1)%matrix,globenv)

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
      CALL write_blocks_maxabsval(matrix_s(1)%matrix,unit_number)
      CALL close_file(unit_number=unit_number)
    END IF

    IF (globenv%print%overlap_matrix) THEN
      CALL write_sparse_matrix(matrix_s(1)%matrix,4,6,qs_env,globenv)
      IF (globenv%print%derivatives) THEN
        DO i=2,SIZE(matrix_s)
          CALL write_sparse_matrix(matrix_s(i)%matrix,4,6,qs_env,globenv)
        END DO
      END IF
    END IF

!   *** Print the kinetic energy integral matrix, if requested ***

    IF (globenv%print%kinetic_energy_matrix) THEN
      CALL write_sparse_matrix(matrix_t(1)%matrix,4,6,qs_env,globenv)
      IF (globenv%print%derivatives) THEN
        DO i=2,SIZE(matrix_t)
          CALL write_sparse_matrix(matrix_t(i)%matrix,4,6,qs_env,globenv)
        END DO
      END IF
    END IF

!   *** Print the potential energy matrix, if requested ***

    IF (globenv%print%potential_energy_matrix) THEN
      CALL allocate_matrix_set(matrix_v,1)
      CALL replicate_matrix_structure(source=matrix_h(1)%matrix,&
                                      target=matrix_v(1)%matrix,&
                                      target_name="POTENTIAL ENERGY MATRIX")
      CALL add_matrices(matrix_v(1)%matrix,1.0_dp,matrix_h(1)%matrix,-1.0_dp,matrix_t(1)%matrix)
      CALL write_sparse_matrix(matrix_v(1)%matrix,4,6,qs_env,globenv)
      CALL deallocate_matrix_set(matrix_v)
    END IF

!   *** Print the core Hamiltonian matrix, if requested ***

    IF (globenv%print%core_hamiltonian_matrix) THEN
      CALL write_sparse_matrix(matrix_h(1)%matrix,4,6,qs_env,globenv)
    END IF

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE build_core_hamiltonian_matrix

! *****************************************************************************

END MODULE qs_core_hamiltonian
