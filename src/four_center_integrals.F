!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!***** cp2k/four_center_integrals [1.0] *
!!
!!   NAME
!!     four_center_integrals
!!
!!   FUNCTION
!!     routines to calculate Hartree-Fock Exchange and the corresponding energy
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Manuel Guidon 
!!
!!   MODIFICATION HISTORY
!!     11.2006 created [Manuel Guidon]
!!
!****************************************************************************
MODULE four_center_integrals
  
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE erf_fn,                          ONLY: erfc
  USE input_constants,                 ONLY: do_hf_screening_none,&
                                             do_hf_screening_short
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_8
  USE libint_interface,                ONLY: evaluate_deriv_eri,&
                                             evaluate_eri,&
                                             evaluate_eri_screen,&
                                             initialize_libderiv,&
                                             initialize_libint,&
                                             lib_deriv,&
                                             lib_int,&
                                             pbc
  USE mathconstants
  USE mathlib,                         ONLY: symmetrize_matrix
  USE message_passing,                 ONLY: mp_max,&
                                             mp_sum,&
                                             mp_sum_i8,&
                                             mp_sync
  USE orbital_pointers
  USE particle_types,                  ONLY: particle_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE sparse_matrix_types,             ONLY: copy_local_sm_to_replicated_fm,&
                                             first_block_node,&
                                             get_block_node,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE x_hf,                            ONLY: x_hf_box_scheme,&
                                             x_hf_cell_type,&
                                             x_hf_distribution,&
                                             x_hf_screening_type,&
                                             x_hf_set_distr_energy,&
                                             x_hf_set_distr_forces,&
                                             x_hf_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
   
  PUBLIC  integrate_four_center, derivatives_four_center, compute_memory_needs, &
          integrate_four_center_symm, derivatives_four_center_symm

  TYPE(lib_int), SAVE      :: lib
  TYPE(lib_deriv), SAVE    :: deriv
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'four_center_integrals'
  REAL(KIND=dp), PARAMETER :: epsilon = 0.00001_dp
   INTEGER :: iteration_count = 0

!***
!****************************************************************************
  
  CONTAINS

!!***** four_center_integrals/coulomb4 [1.0] *
!!
!!   NAME
!!     coulomb4
!!
!!   FUNCTION
!!     calculates two-electron integrals of a quartet/shell using the library 
!!     lib_int
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - r: position
!!     - z: zeta
!!     - l: angular momentum
!!     - npgf: number of primitive cartesian gaussian in actual shell
!!     - primitive_integrals: array of primitive_integrals
!!     - screening_parameter: contains info for libint 
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     11.2006 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE coulomb4(ra, rb, rc, rd, zeta, zetb, zetc, zetd,&
                      la_min, la_max, lb_min, lb_max,&
                      lc_min, lc_max, ld_min, ld_max,&
                      npgfa, npgfb, npgfc, npgfd,&
                      ncoa,ncob,ncoc,ncod,&
                      primitive_integrals,&
                      screening_parameter)
  
    REAL(dp), INTENT(IN)                     :: ra(3), rb(3), rc(3), rd(3)
    INTEGER, INTENT(IN)                      :: npgfa, ld_max, ld_min, &
                                                lc_max, lc_min, lb_max, &
                                                lb_min, la_max, la_min, &
                                                npgfd, npgfc, npgfb
    REAL(dp), DIMENSION(1:npgfa), INTENT(IN) :: zeta
    REAL(dp), DIMENSION(1:npgfb), INTENT(IN) :: zetb
    REAL(dp), DIMENSION(1:npgfc), INTENT(IN) :: zetc
    REAL(dp), DIMENSION(1:npgfd), INTENT(IN) :: zetd
    INTEGER, INTENT(IN)                      :: ncoa, ncob, ncoc, ncod
    REAL(dp), &
      DIMENSION(ncoa, ncob, ncoc, ncod)      :: primitive_integrals
    TYPE(x_hf_screening_type)                :: screening_parameter

    INTEGER :: ipgf, jpgf, kpgf, la, lb, lc, ld, lpgf, offset_a, offset_a1, &
      offset_b, offset_b1, offset_c, offset_c1, offset_d, offset_d1

    primitive_integrals = 0.0_dp

    DO lpgf = 1,npgfd
      offset_d1 = (lpgf-1)*ncoset(ld_max)
      DO ld = ld_min, ld_max
        offset_d = offset_d1 + ncoset(ld-1)
        DO kpgf = 1,npgfc
          offset_c1 = (kpgf-1)*ncoset(lc_max)
          DO lc = lc_min, lc_max
            offset_c = offset_c1 + ncoset(lc-1)
            DO jpgf = 1,npgfb
              offset_b1 = (jpgf-1)*ncoset(lb_max)
              DO lb = lb_min, lb_max
                offset_b = offset_b1 + ncoset(lb-1)
                DO ipgf = 1,npgfa
                  offset_a1 = (ipgf-1)*ncoset(la_max)
                  DO la = la_min,la_max
                    offset_a = offset_a1 + ncoset(la-1)
                    !Build primitives
                    CALL evaluate_eri(lib, ra, rb, rc, rd,&
                                      zeta(ipgf), zetb(jpgf), zetc(kpgf), zetd(lpgf),&
                                      la, lb, lc ,ld,&
                                      ncoa,ncob,ncoc,ncod,&
                                      offset_a,offset_b,offset_c,offset_d, &
                                      primitive_integrals, &
                                      screening_parameter)
                  END DO !la
                END DO !ipgf
              END DO !lb
            END DO !jpgf
          END DO !lc
        END DO !kpgf
      END DO !ld
    END DO !lpgf
  END SUBROUTINE coulomb4
!***************************************************************************

!!***** four_center_integrals/screen4 [1.0] *
!!
!!   NAME
!!     screen4
!!
!!   FUNCTION
!!     calculates max values of two-electron integrals in a quartet/shell
!!     w.r.t. different zetas using the library lib_int
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - r: position
!!     - z: zeta
!!     - l: angular momentum
!!     - npgf: number of primitive cartesian gaussian in actual shell
!!     - max_val: schwarz screening value
!!     - screening_parameter: contains info for libint
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!    03.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE screen4(ra, rb, zeta, zetb, &
                      la_min, la_max, lb_min, lb_max,&
                      npgfa, npgfb, &
                      max_val, screening_parameter)

    REAL(dp), INTENT(IN)                     :: ra(3), rb(3)
    REAL(dp), DIMENSION(:), INTENT(IN)       :: zeta, zetb
    INTEGER, INTENT(IN)                      :: la_min, la_max, lb_min, &
                                                lb_max, npgfa, npgfb
    REAL(dp), INTENT(INOUT)                  :: max_val
    TYPE(x_hf_screening_type)                :: screening_parameter

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'screen4'

    INTEGER                                  :: ipgf, jpgf, la, lb
    REAL(dp)                                 :: max_val_temp

    max_val_temp = max_val
    DO ipgf = 1,npgfa
      DO jpgf = 1,npgfb
        DO la = la_min,la_max
          DO lb = lb_min, lb_max
            !Build primitives
            CALL evaluate_eri_screen(lib, ra, rb, ra, rb,&
                                     zeta(ipgf), zetb(jpgf), zeta(ipgf), zetb(jpgf),&
                                     la, lb, la ,lb,&
                                     max_val_temp, screening_parameter)
            max_val = MAX(max_val, max_val_temp)
          END DO !lb
        END DO !la
      END DO !jpgf
    END DO !ipgf

  END SUBROUTINE screen4


!!***** four_center_integrals/compute_screening_matrices [1.0] *
!!
!!   NAME
!!     compute_screening_matrices
!!
!!   FUNCTION
!!     calculates the screening matrices 
!!
!!   NOTES
!!     - already computes the square root
!!     - also takes the primitive to sgf transformation into account
!!     - in this way all stored contracted integrals are smaller or equal to their estimate
!!   INPUTS
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!    05.2007 refactored [Joost VandeVondele]
!!
!!****************************************************************************

  SUBROUTINE compute_screening_matrices(screen_set,screen_atom,screen_max,particle_set,&
                                        atomic_kind_set,screening_parameter,read_from_memory,buffer,error)

    REAL(kind=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: screen_set
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: screen_atom
    REAL(KIND=dp)                            :: screen_max
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(x_hf_screening_type)                :: screening_parameter
    LOGICAL                                  :: read_from_memory
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      OPTIONAL, POINTER                      :: buffer
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_screening_matrices', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, jatom, jset, k, &
                                                katom, kset, max_set, natom, &
                                                ncob, nsetb, nsetc, sgfb, stat
    INTEGER, DIMENSION(:), POINTER           :: lb_max, lb_min, lc_max, &
                                                lc_min, npgfb, npgfc, nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfb
    LOGICAL                                  :: failure
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: contraction_max
    REAL(dp), DIMENSION(:, :), POINTER       :: sphi_b, zetb, zetc
    REAL(KIND=dp)                            :: max_val, p2s_b, p2s_c, rb(3), &
                                                rc(3)
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_b, orb_basis_c

    CALL timeset(routineN,handle)

    failure = .FALSE.

    natom=SIZE(particle_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, maxnset=max_set)
    ALLOCATE(screen_set(natom,natom,max_set,max_set),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(screen_atom(natom,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(contraction_max(max_set,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    screen_set=0.0_dp
    screen_atom=0.0_dp
    contraction_max=0.0_dp

    DO jatom=1,natom

      atom_kind => particle_set(jatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_b)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                             lmax=lb_max,&
                             npgf=npgfb,&
                             nset=nsetb,&
                             first_sgf=first_sgfb,&
                             sphi=sphi_b,&
                             nsgf_set=nsgfb)
      DO jset = 1,nsetb
        ! takes the primitive to contracted transformation into account
        ncob = npgfb(jset)*ncoset(lb_max(jset))
        sgfb = first_sgfb(1,jset)
        ! if the primitives are assumed to be all of max_val2, max_val2*p2s_b becomes
        ! the maximum value after multiplication with sphi_b
        contraction_max(jset,jatom) = MAXVAL((/(SUM(ABS(sphi_b(1:ncob,k))),k=sgfb,sgfb+nsgfb(jset)-1)/))
      ENDDO
    ENDDO

    DO jatom=1,natom

      atom_kind => particle_set(jatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_b)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                             lmax=lb_max,&
                             lmin=lb_min,&
                             npgf=npgfb,&
                             nset=nsetb,&
                             zet=zetb)
      rb=particle_set(jatom)%r(:)
      DO katom = 1,natom
        atom_kind => particle_set(katom)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atom_kind,&
                             orb_basis_set=orb_basis_c)
        CALL get_gto_basis_set(gto_basis_set=orb_basis_c,&
                               lmax=lc_max,&
                               lmin=lc_min,&
                               npgf=npgfc,&
                               nset=nsetc,&
                               zet=zetc)
        rc=particle_set(katom)%r(:)
        DO jset = 1,nsetb
          p2s_b = contraction_max(jset,jatom)
          DO kset = 1,nsetc
            p2s_c = contraction_max(kset,katom)
            max_val = 0.0_dp
            IF( read_from_memory ) THEN
              screen_set(jatom,katom,jset,kset) = buffer(jatom,katom,jset,kset)
            ELSE
              CALL screen4(rb, rc, &
                            zetb(:,jset), zetc(:,kset), &
                            lb_min(jset), lb_max(jset), lc_min(kset), lc_max(kset),&
                            npgfb(jset), npgfc(kset), &
                            max_val, screening_parameter)
              screen_set(jatom,katom,jset,kset) = SQRT(max_val)*p2s_b*p2s_c
              IF( PRESENT(buffer) ) buffer(jatom,katom,jset,kset) = screen_set(jatom,katom,jset,kset)
            END IF
          END DO
        END DO
        screen_atom(jatom,katom)=MAXVAL(screen_set(jatom,katom,1:nsetb,1:nsetc))
      END DO
    END DO
    screen_max=MAXVAL(screen_atom)
    DEALLOCATE(contraction_max, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle)
  END SUBROUTINE compute_screening_matrices

!!***** four_center_integrals/compute_pmax_matrices [1.0] *
!!
!!   NAME
!!     compute_pmax_matrices
!!
!!   FUNCTION
!!     calculates the maximum of the density matrix in compressed form for screening purposes 
!!
!!   NOTES
!!     - computes for each pair of shells the maximum absolute value of p
!!
!!   INPUTS
!!     - screen_pmat_forces: if false just returns a matrix of 1.0_dp
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!    05.2007 created [Joost VandeVondele]
!!
!****************************************************************************
  
  SUBROUTINE compute_pmax_matrices(pmax_set,pmax_atom,pmax_max, &
              full_density,full_density_beta,screen_pmat_forces,particle_set,atomic_kind_set,error)

    REAL(kind=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: pmax_set
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: pmax_atom
    REAL(KIND=dp)                            :: pmax_max
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: full_density, &
                                                full_density_beta
    LOGICAL, INTENT(IN)                      :: screen_pmat_forces
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_pmax_matrices', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, jatom, jset, katom, kset, max_set, mb, mc, natom, &
      nsetb, nsetc, nsgfb_total, nsgfc_total, pb, pb1, pb2, pc, pc1, pc2, &
      sgfb, sgfc, stat
    INTEGER, DIMENSION(:), POINTER           :: last_sgf_global, nsgfb, nsgfc
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfb, first_sgfc
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: pmax_tmp
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_b, orb_basis_c

    CALL timeset(routineN,handle)
    failure = .FALSE.

    natom=SIZE(particle_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, maxnset=max_set)

    ! create a mapping array of easy access to the density matrix elements
    ALLOCATE(last_sgf_global(0:natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    last_sgf_global(0)=0
    DO jatom=1,natom
      atom_kind => particle_set(jatom)%atomic_kind
 
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_b)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                             nsgf=nsgfb_total)
      last_sgf_global(jatom) = last_sgf_global(jatom-1)+nsgfb_total
    END DO

    ! allocate and compute pmax
    ALLOCATE(pmax_set(max_set,max_set,natom,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(pmax_atom(natom,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF (screen_pmat_forces) THEN

       pmax_set=0.0_dp
       pmax_atom=0.0_dp

       DO jatom=1,natom

         atom_kind => particle_set(jatom)%atomic_kind
         CALL get_atomic_kind(atomic_kind=atom_kind,&
                              orb_basis_set=orb_basis_b)
         CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                                nsgf=nsgfb_total,&
                                nset=nsetb,& 
                                first_sgf=first_sgfb,&
                                nsgf_set=nsgfb)
         pb1 = last_sgf_global(jatom) - nsgfb_total

         DO katom = 1,natom
           atom_kind => particle_set(katom)%atomic_kind
           CALL get_atomic_kind(atomic_kind=atom_kind,&
                                orb_basis_set=orb_basis_c)
           CALL get_gto_basis_set(gto_basis_set=orb_basis_c,&
                                  nsgf=nsgfc_total,&
                                  nset=nsetc,&
                                  first_sgf=first_sgfc,&
                                  nsgf_set=nsgfc)
           pc1 = last_sgf_global(katom) - nsgfc_total

           DO jset = 1,nsetb
             sgfb = first_sgfb(1,jset)
             DO kset = 1,nsetc
               sgfc = first_sgfc(1,kset)
               pmax_tmp = 0.0_dp
               DO mc = 1,nsgfc(kset)
                  pc2 = mc+sgfc-1
                  pc = pc1 + pc2
                  DO mb = 1,nsgfb(jset)
                     pb2 = mb+sgfb-1
                     pb = pb1 + pb2
                     pmax_tmp=MAX(pmax_tmp,ABS(full_density(pb,pc)))
                     IF( ASSOCIATED(full_density_beta) ) pmax_tmp=MAX(pmax_tmp,ABS(full_density_beta(pb,pc)))
                  ENDDO
               ENDDO
               pmax_set(jset,kset,jatom,katom) = pmax_tmp
             END DO
           END DO
           pmax_atom(jatom,katom)=MAXVAL(pmax_set(1:nsetb,1:nsetc,jatom,katom))
         END DO
       END DO

    ELSE

       pmax_set=1.0_dp
       pmax_atom=1.0_dp

    ENDIF

    pmax_max=MAXVAL(pmax_atom)

    DEALLOCATE(last_sgf_global,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE compute_pmax_matrices
!!***** four_center_integrals/compute_memory_needs [1.0] *
!!
!!   NAME
!!     compute_memory_needs
!!
!!   FUNCTION
!!     calculates the needed memory to store the matrix elements, taking the screening into account
!!
!!   NOTES
!!     based on the row/col distribution given, computes for each local block_ab the needed memory to store the corresponding integrals
!!     this nees to remain in sync with the integrate_four_center routine.
!!     This distribution is then used to load balance HFX.
!!
!!   INPUTS
!!     memory block ab : output : each element corresponds to the number of integrals that need
!!                                to be stored if this block of H_ks is local to the CPU.
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!    05.2007 created [Joost VandeVondele]
!!
!!****************************************************************************

  SUBROUTINE compute_memory_needs(memory_blocks_ab,eps_schwarz,particle_set,atomic_kind_set,para_env,&
                                  screening_parameter,error)
    INTEGER, DIMENSION(:, :), INTENT(OUT)    :: memory_blocks_ab
    REAL(KIND=dp), INTENT(IN)                :: eps_schwarz
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(x_hf_screening_type)                :: screening_parameter
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_memory_needs', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iatom, iblock, iset, jatom, jset, katom, kset, latom, &
      lset, max_am, natom, nseta, nsetb, nsetc, nsetd, stat
    INTEGER, DIMENSION(:), POINTER           :: nsgfa, nsgfb, nsgfc, nsgfd
    LOGICAL                                  :: failure, include_ab
    REAL(KIND=dp)                            :: max_val1, max_val2, screen_max
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: screen_atom
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: screen_set
    TYPE(gto_basis_set_type), POINTER        :: orb_basis

    CALL timeset(routineN,handle)
    failure=.FALSE.

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, maxlgto=max_am)
    CALL initialize_libint(lib, max_am, error)

    memory_blocks_ab = 0

    ! set up screening matrices
    CALL compute_screening_matrices(screen_set,screen_atom,screen_max,particle_set,&
                                    atomic_kind_set,screening_parameter,.FALSE.,error=error)

    ! loop over this CPU's share of the atoms and compute the memory cost of the atomic blocks
    ! do not compute this for blocks that will not be associated
    natom=SIZE(particle_set)
    iblock=0
    DO iatom=1,natom
      DO jatom=1,natom
        iblock=iblock+1
        IF (MOD(iblock,para_env%num_pe) .NE. para_env%mepos) CYCLE
        CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, orb_basis_set=orb_basis)
        CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nseta, nsgf_set=nsgfa)
        CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind, orb_basis_set=orb_basis)
        CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetb, nsgf_set=nsgfb)
        DO katom=1,natom
          CALL get_atomic_kind(atomic_kind=particle_set(katom)%atomic_kind, orb_basis_set=orb_basis)
          CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetc, nsgf_set=nsgfc)
          DO latom=1,natom
            ! IF( ((iatom <= jatom).AND.(katom<=latom)).AND.(katom+latom<=iatom+jatom))  THEN
            IF( ((.TRUE.).AND.(katom<=latom)).AND.(katom+latom<=iatom+jatom))  THEN
              IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<MIN(iatom,jatom))) CYCLE
    
            
              ! take the symmetric storage form into account
               IF (iatom > jatom) THEN
                 include_ab = (MODULO(iatom + jatom,2) /= 0)
               ELSE
                 include_ab = (MODULO(iatom + jatom,2) == 0)
               END IF
               IF (.NOT. include_ab) CYCLE

              CALL get_atomic_kind(atomic_kind=particle_set(latom)%atomic_kind, orb_basis_set=orb_basis)
              CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetd, nsgf_set=nsgfd)

              DO iset = 1, nseta
                DO jset = 1,nsetb
                  max_val1 = screen_set(iatom,jatom,iset,jset)
                  DO kset = 1,nsetc
                    DO lset = 1,nsetd
                      max_val2 = max_val1 * screen_set(katom,latom,kset,lset)
                      IF(max_val2<eps_schwarz) CYCLE
                      memory_blocks_ab(iatom,jatom)=memory_blocks_ab(iatom,jatom) + &
                                                    nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            END IF
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DEALLOCATE(screen_set,screen_atom,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL mp_sum(memory_blocks_ab,para_env%group)

    CALL timestop(handle)

  END SUBROUTINE compute_memory_needs

!!***** four_center_integrals/integrate_four_center [1.0] *
!!
!!   NAME
!!     integrate_four_center 
!!
!!   FUNCTION
!!     computes four center integrals for a full basis set and updates the
!!     Kohn-Sham-Matrix and energy
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_env 
!!     - ks_matrix
!!     - energy
!!     - rho: density matrix
!!     - para_env
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     11.2006 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE integrate_four_center(qs_env,ks_matrix,energy,rho,hf_section,para_env,& 
                                   geometry_did_change,do_periodic, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: hf_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL                                  :: geometry_did_change, &
                                                do_periodic
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'integrate_four_center', &
      routineP = moduleN//':'//routineN

    INTEGER :: current_counter, handle, handle_getP, handle_main, iatom, &
      integral_counter, iset, iw, jatom, jset, katom, kset, latom, lset, ma, &
      max_am, max_k, max_set, maxsgf, mb, mc, md, natom, nco_max, ncoa, ncob, &
      ncoc, ncod, nseta, nsetb, nsetc, nsetd, nsgfa_total, nsgfb_total, &
      nsgfc_total, nsgfd_total, nspins, number_of_precalc_integrals, p_index, &
      pa, pa1, pa2, pb, pb1, pb2, pc, pc1, pc2, pd, pd1, pd2, sgfa, sgfb, &
      sgfc, sgfd, stat, total_integrals_to_store, &
      total_integrals_to_store_aver, total_integrals_to_store_max
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: last_sgf_global
    INTEGER, DIMENSION(:), POINTER :: la_max, la_min, lb_max, lb_min, lc_max, &
      lc_min, ld_max, ld_min, npgfa, npgfb, npgfc, npgfd, nsgfa, nsgfb, &
      nsgfc, nsgfd
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb, &
                                                first_sgfc, first_sgfd
    LOGICAL                                  :: buffer_overflow, failure
    REAL(dp) :: box_size(3), ene_x_aa, ene_x_bb, eps_schwarz, fac, &
      hf_fraction, max_val1, max_val2, memory_usage, ra(3), rb(3), rc(3), &
      rd(3), screen_max, tmp1, tmp11, tmp2, tmp21, tmp3, tmp31
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: primitive_integrals, T1
    REAL(dp), DIMENSION(:, :), POINTER       :: screen_atom, sphi_a, sphi_b, &
                                                sphi_c, sphi_d, zeta, zetb, &
                                                zetc, zetd
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: screen_set
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: full_density, &
                                                full_density_beta, &
                                                sparse_block, &
                                                sparse_block_beta
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), POINTER      :: screen_periodic
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_a, orb_basis_b, &
                                                orb_basis_c, orb_basis_d
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(real_block_node_type), POINTER      :: block_node, block_node_beta
    TYPE(x_hf_screening_type)                :: screening_parameter

! many timings somewhat more clear

    CALL mp_sync(para_env%group)
    CALL timeset(routineN,"I",'',handle)

    CALL section_vals_val_get(hf_section,"FRACTION",r_val=hf_fraction,error=error)
    CALL section_vals_val_get(hf_section,"EPS_SCHWARZ",r_val=eps_schwarz,error=error)
   
    number_of_precalc_integrals = qs_env%x_data%number_of_elements
    screening_parameter = qs_env%x_data%screening_parameter
    integral_counter = 1
    memory_usage = 0.0_dp
    buffer_overflow = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    error=error)
 
    natom = SIZE(particle_set,1)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxsgf=maxsgf,&
                             maxnset=max_set,&
                             maxlgto=max_am)
    
    CALL initialize_libint(lib, max_am, error)
    
    nspins = qs_env%dft_control%nspins
 
    ALLOCATE(last_sgf_global(0:natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    last_sgf_global(0)=0
    DO iatom=1,natom
      atom_kind => particle_set(iatom)%atomic_kind
 
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_a)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_a,&
                             nsgf=nsgfa_total)
      last_sgf_global(iatom) = last_sgf_global(iatom-1)+nsgfa_total
    END DO
 
    ene_x_aa = 0.0_dp; ene_x_bb = 0.0_dp

    CALL timeset(routineN//"_getP",handle_getP)

    ALLOCATE(full_density(1,1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL copy_local_sm_to_replicated_fm(rho%rho_ao(1)%matrix,full_density)
    CALL symmetrize_matrix(full_density,"upper_to_lower")
    CALL mp_sum(full_density,para_env%group)
    IF(nspins == 2) THEN
      ALLOCATE(full_density_beta(1,1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      CALL copy_local_sm_to_replicated_fm(rho%rho_ao(2)%matrix,full_density_beta)
      CALL symmetrize_matrix(full_density_beta,"upper_to_lower")
      CALL mp_sum(full_density_beta,para_env%group)
    END IF

    CALL timestop(handle_getP)


    CALL timeset(routineN//"_main",handle_main)
    max_k = qs_env%x_data%number_of_k_loops
    IF( do_periodic ) THEN
      box_size(1) = qs_env%cell%hmat(1,1)
      box_size(2) = qs_env%cell%hmat(2,2)
      box_size(3) = qs_env%cell%hmat(3,3)
      max_k = MIN(max_k, SIZE(qs_env%x_data%neighbor_cells))
      ALLOCATE(screen_periodic(natom,natom,max_set,max_set,SIZE(qs_env%x_data%neighbor_cells)))
      CALL comp_screen_mat_per(screen_set,screen_atom,screen_max,screen_periodic,particle_set,atomic_kind_set,&
                               screening_parameter,.NOT.geometry_did_change,&
                               qs_env%x_data%x_screening_buffer_periodic,&
                               qs_env%x_data%neighbor_cells,box_size,max_k,error)
    ELSE
      CALL compute_screening_matrices(screen_set,screen_atom,screen_max,particle_set,atomic_kind_set,&
                                    screening_parameter,.NOT.geometry_did_change,qs_env%x_data%x_screening_buffer, error)
    END IF
    ! precompute maximum nco and allocate scratch
    nco_max=0
    DO iatom=1,natom
      atom_kind => particle_set(iatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_a)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_a, lmax=la_max, lmin=la_min,&
                             npgf=npgfa, nset=nseta, zet=zeta, nsgf_set=nsgfa,&
                             first_sgf=first_sgfa, sphi=sphi_a, nsgf=nsgfa_total)
      DO iset = 1, nseta
         nco_max = MAX(nco_max,npgfa(iset)*ncoset(la_max(iset)))
      ENDDO
    ENDDO
    ALLOCATE(primitive_integrals(nco_max**4),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    primitive_integrals=0.0_dp
    ALLOCATE(T1(nco_max**4),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    T1=0.0_dp


    total_integrals_to_store=0
    memory_usage = 0.0_dp
    ene_x_aa = 0.0_dp
    ene_x_bb = 0.0_dp
    IF(nspins == 1) THEN
      fac = 0.5_dp * hf_fraction
    ELSE
      fac = 1.0_dp * hf_fraction
    END IF

    tmp11 = 0.0_dp
    DO iatom=1,ks_matrix(1)%matrix%nblock_row

      block_node => first_block_node(ks_matrix(1)%matrix,iatom)
      IF(nspins == 2) block_node_beta => first_block_node(ks_matrix(2)%matrix,iatom)

      ! parallel distribution might mean that there won't be an associated latom
      IF (.NOT. ASSOCIATED(block_node)) CYCLE
 
      atom_kind => particle_set(iatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_a)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_a, lmax=la_max, lmin=la_min,&
                             npgf=npgfa, nset=nseta, zet=zeta, nsgf_set=nsgfa,&
                             first_sgf=first_sgfa, sphi=sphi_a, nsgf=nsgfa_total)
      ra=particle_set(iatom)%r(:)
      pa1 = last_sgf_global(iatom) - nsgfa_total

      DO jatom = 1,natom

        IF (screen_atom(iatom,jatom) * screen_max < eps_schwarz) CYCLE

        atom_kind => particle_set(jatom)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atom_kind,&
                             orb_basis_set=orb_basis_b)
        CALL get_gto_basis_set(gto_basis_set=orb_basis_b, lmax=lb_max, lmin=lb_min,&
                               npgf=npgfb, nset=nsetb, zet=zetb, nsgf_set=nsgfb,&
                               first_sgf=first_sgfb, sphi=sphi_b, nsgf=nsgfb_total)
        rb=particle_set(jatom)%r(:)
        pb1 = last_sgf_global(jatom) - nsgfb_total

        DO katom = 1,natom
          atom_kind => particle_set(katom)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atom_kind,&
                               orb_basis_set=orb_basis_c)
          CALL get_gto_basis_set(gto_basis_set=orb_basis_c, lmax=lc_max, lmin=lc_min,&
                                 npgf=npgfc, nset=nsetc, zet=zetc, nsgf_set=nsgfc,&
                                 first_sgf=first_sgfc, sphi=sphi_c, nsgf=nsgfc_total)
          rc=particle_set(katom)%r(:)
          pc1 = last_sgf_global(katom) - nsgfc_total

          block_node => first_block_node(ks_matrix(1)%matrix,iatom)
          IF(nspins == 2) block_node_beta => first_block_node(ks_matrix(2)%matrix,iatom)
          DO WHILE (ASSOCIATED(block_node))
           CALL get_block_node(block_node=block_node,&
                                block_col=latom,&
                                BLOCK=sparse_block)

            block_node => next_block_node(block_node)      
            IF(nspins == 2) THEN
              CALL get_block_node(block_node=block_node_beta,&
                                  block_col=latom,&
                                  BLOCK=sparse_block_beta)
              block_node_beta => next_block_node(block_node_beta)
            END IF

            IF (screen_atom(iatom,jatom) * screen_atom(katom,latom) < eps_schwarz) CYCLE

            atom_kind => particle_set(latom)%atomic_kind
            CALL get_atomic_kind(atomic_kind=atom_kind,&
                                 orb_basis_set=orb_basis_d)
            CALL get_gto_basis_set(gto_basis_set=orb_basis_d, lmax=ld_max, lmin=ld_min,&
                                   npgf=npgfd, nset=nsetd, zet=zetd, nsgf_set=nsgfd,&
                                   first_sgf=first_sgfd, sphi=sphi_d, nsgf=nsgfd_total)
            rd=particle_set(latom)%r(:)
            pd1 = last_sgf_global(latom) - nsgfd_total
           
            IF(iatom<latom) THEN
              tmp3 = 2.0_dp
              tmp31 = 2.0_dp
            ELSE
              tmp3 = 1.0_dp
              tmp31 = 1.0_dp
            END IF

            DO iset = 1, nseta
              ncoa = npgfa(iset)*ncoset(la_max(iset))
              sgfa = first_sgfa(1,iset)
              DO jset = 1,nsetb
                ncob = npgfb(jset)*ncoset(lb_max(jset))
                sgfb = first_sgfb(1,jset)
                max_val1 = screen_set(iatom,jatom,iset,jset)
                IF (max_val1 * screen_atom(katom,latom) < eps_schwarz) CYCLE
                DO kset = 1,nsetc
                  ncoc = npgfc(kset)*ncoset(lc_max(kset))
                  sgfc = first_sgfc(1,kset)
                  DO lset = 1,nsetd
                    ncod = npgfd(lset)*ncoset(ld_max(lset)) 
                    sgfd = first_sgfd(1,lset)
                    max_val2 = max_val1 * screen_set(katom,latom,kset,lset)
                    IF(max_val2<eps_schwarz) CYCLE
                    
                    current_counter = nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)
                    total_integrals_to_store=total_integrals_to_store+current_counter
                    buffer_overflow = integral_counter + current_counter - 1 > number_of_precalc_integrals
                    IF(.NOT.buffer_overflow .AND. .NOT.geometry_did_change) THEN 
                      DO md = 1,nsgfd(lset)
                        pd2 = md+sgfd-1
                        pd = pd1 + pd2
                        DO mc = 1,nsgfc(kset)
                          pc2 = mc+sgfc-1
                          pc = pc1 + pc2
                          DO mb = 1,nsgfb(jset)
                            pb2 = mb+sgfb-1
                            pb = pb1 + pb2
                            tmp1 = fac*full_density(pb,pc)
                            IF( nspins == 2 ) tmp11 = fac*full_density_beta(pb,pc)
                            DO ma = 1,nsgfa(iset)
                              pa2 = ma+sgfa-1
                              pa = pa1 + pa2
                              tmp2 = tmp1*qs_env%x_data%x_buffer(integral_counter)
                              sparse_block(pa2,pd2) = sparse_block(pa2,pd2) - tmp2
                              tmp2 = tmp3*full_density(pa,pd)*tmp2
                              ene_x_aa = ene_x_aa - tmp2
                              IF( nspins == 2 ) THEN
                                tmp21 = tmp11*qs_env%x_data%x_buffer(integral_counter)
                                sparse_block_beta(pa2,pd2) = sparse_block_beta(pa2,pd2) - tmp21
                                tmp21 = tmp31*full_density_beta(pa,pd)*tmp21
                                ene_x_bb = ene_x_bb -tmp21
                              END IF
                              integral_counter = integral_counter + 1
                            END DO
                          END DO
                        END DO
                      END DO
                    END IF
                    IF(geometry_did_change .OR. buffer_overflow) THEN
                      IF( .NOT. do_periodic ) THEN
                        CALL coulomb4(ra, rb, rc, rd,&
                                      zeta(:,iset), zetb(:,jset), zetc(:,kset), zetd(:,lset),&
                                      la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                                      lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                                      npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset),&
                                      ncoa,ncob,ncoc,ncod,&
                                      primitive_integrals,&
                                      screening_parameter)
                      ELSE
                        CALL coulomb4_periodic(ra, rb, rc, rd,&
                                               zeta(1:npgfa(iset),iset), zetb(1:npgfb(jset),jset), &
                                               zetc(1:npgfc(kset),kset), zetd(1:npgfd(lset),lset),&
                                               la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                                               lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                                               npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset),&
                                               ncoa,ncob,ncoc,ncod,&
                                               primitive_integrals(1),&
                                               screening_parameter, box_size, &
                                               qs_env%x_data%neighbor_cells, screen_periodic(iatom,jatom,iset,jset,:),&
                                               screen_periodic(katom,latom,kset,lset,:), eps_schwarz, max_k)
                      END IF


                      CALL dgemm("T","N",ncob*ncoc*ncod,nsgfa(iset),ncoa,&
                                 1.0_dp, primitive_integrals(1),ncoa,&
                                 sphi_a(1,sgfa), SIZE(sphi_a,1),&
                                 0.0_dp, T1(1),ncob*ncoc*ncod)

                      CALL dgemm("T","N",nsgfa(iset)*ncoc*ncod,nsgfb(jset),ncob,&
                                 1.0_dp, T1(1),ncob,&
                                 sphi_b(1,sgfb), SIZE(sphi_b,1),&
                                 0.0_dp, primitive_integrals(1),nsgfa(iset)*ncoc*ncod)

                      CALL dgemm("T","N",nsgfa(iset)*nsgfb(jset)*ncod,nsgfc(kset),ncoc,&
                                 1.0_dp, primitive_integrals(1),ncoc,&
                                 sphi_c(1,sgfc), SIZE(sphi_c,1),&
                                 0.0_dp, T1(1),nsgfa(iset)*nsgfb(jset)*ncod)

                      CALL dgemm("T","N",nsgfa(iset)*nsgfb(jset)*nsgfc(kset),nsgfd(lset),ncod,&
                                 1.0_dp, T1(1),ncod,&
                                 sphi_d(1,sgfd), SIZE(sphi_d,1),&
                                 0.0_dp, primitive_integrals(1),nsgfa(iset)*nsgfb(jset)*nsgfc(kset))


                      p_index=0
                      DO md = 1,nsgfd(lset)
                        pd2 = md+sgfd-1
                        pd = pd1 + pd2
                        DO mc = 1,nsgfc(kset)
                          pc2 = mc+sgfc-1
                          pc = pc1 + pc2
                          DO mb = 1,nsgfb(jset)
                            pb2 = mb+sgfb-1
                            pb = pb1 + pb2
                            tmp1 = fac*full_density(pb,pc)
                            IF( nspins == 2 ) tmp11 = fac*full_density_beta(pb,pc)
                            DO ma = 1,nsgfa(iset)
                              pa2 = ma+sgfa-1
                              pa = pa1 + pa2
                              p_index=p_index+1
                              tmp2 = tmp1*primitive_integrals(p_index)
                              sparse_block(pa2,pd2) = sparse_block(pa2,pd2) - tmp2
                              tmp2 = tmp3*full_density(pa,pd)*tmp2
                              ene_x_aa = ene_x_aa - tmp2
                              IF( nspins == 2 ) THEN
                                tmp21 = tmp11*primitive_integrals(p_index)
                                sparse_block_beta(pa2,pd2) = sparse_block_beta(pa2,pd2) - tmp21
                                tmp21 = tmp31*full_density_beta(pa,pd)*tmp21
                                ene_x_bb = ene_x_bb -tmp21
                              END IF
                              IF(geometry_did_change .AND. (.NOT.buffer_overflow)) THEN
                                qs_env%x_data%x_buffer(integral_counter) = &
                                   primitive_integrals(p_index)
                                integral_counter = integral_counter + 1
                              END IF
                            END DO
                          END DO
                        END DO
                      END DO
                    ENDIF
                  END DO !lset
                END DO !kset
              END DO !jset
            END DO !set

          END DO !latom
        END DO !katom
      END DO !jatom
    END DO !iatom

    DEALLOCATE(primitive_integrals,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(T1,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL mp_sum(ene_x_aa,para_env%group)
    IF(nspins==2) CALL mp_sum(ene_x_bb,para_env%group)
    energy%ex = 0.5_dp*(ene_x_aa+ene_x_bb)

    CALL timestop(handle_main)
    
    iw = cp_print_key_unit_nr(logger,hf_section,"HF_INFO",&
             extension=".scfLog",error=error)
    total_integrals_to_store_max=total_integrals_to_store
    total_integrals_to_store_aver=(total_integrals_to_store+para_env%num_pe-1)/(para_env%num_pe)
    CALL mp_max(total_integrals_to_store_max,para_env%group)
    CALL mp_sum(total_integrals_to_store_aver,para_env%group)

    IF (iw>0) THEN
      IF(memory_usage < integral_counter * dp_size /1024_dp/1024_dp) THEN
        memory_usage = integral_counter * dp_size /1024_dp/1024_dp
      END IF
      WRITE (UNIT=iw,FMT="(/,(T3,A,T60,F20.10))")&
               "Memory used for x_buffer [MB]:                 ",memory_usage
      WRITE (UNIT=iw,FMT="((T3,A,T60,I20))")&
               "Number of integrals stored in buffer:          ",integral_counter-1
      WRITE (UNIT=iw,FMT="((T3,A,T60,I20))")&
               "Maximum number of elements to store:          ",total_integrals_to_store_max
      WRITE (UNIT=iw,FMT="((T3,A,T60,I20),/)")&
               "Average number of elements to store:          ",total_integrals_to_store_aver
    END IF
    CALL cp_print_key_finished_output(iw,logger,hf_section,&
                                      "HF_INFO", error=error)
    
    IF( do_periodic ) THEN
      DEALLOCATE(screen_periodic, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    DEALLOCATE(screen_set, screen_atom, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(last_sgf_global,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(full_density,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(nspins==2) THEN
      DEALLOCATE(full_density_beta,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    !
    ! this wraps to free_libint, but currently causes a segfault
    ! as a result, we don't call it, but some memory remains allocated
    !  
    ! CALL terminate_libint(lib)

    ! many timings somewhat more clear
    CALL mp_sync(para_env%group)

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE integrate_four_center

!***************************************************************************

!!***** four_center_integrals/forces4 [1.0] *
!!
!!   NAME
!!     forces4
!!
!!   FUNCTION
!!     calculates two-electron derivatived of a quartet/shell using the library 
!!     lib_deriv
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - r: position
!!     - z: zeta
!!     - l: angular momentum
!!     - npgf: number of primitive cartesian gaussian in actual shell
!!     - primitive_forces: array of primitive_forces
!!     - work_forces: work_storage
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     03.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE forces4(ra, rb, rc, rd, zeta, zetb, zetc, zetd,&
                    la_min, la_max, lb_min, lb_max,&
                    lc_min, lc_max, ld_min, ld_max,&
                    npgfa, npgfb, npgfc, npgfd,&
                    ncoa, ncob, ncoc, ncod,&
                    primitive_forces, work_forces, eps_schwarz,&
                    screening_parameter )
  
    REAL(dp), INTENT(IN)                     :: ra(3), rb(3), rc(3), rd(3)
    REAL(dp), DIMENSION(:), INTENT(IN)       :: zeta, zetb, zetc, zetd
    INTEGER, INTENT(IN) :: la_min, la_max, lb_min, lb_max, lc_min, lc_max, &
      ld_min, ld_max, npgfa, npgfb, npgfc, npgfd, ncoa, ncob, ncoc, ncod
    REAL(dp), &
      DIMENSION(ncoa, ncob, ncoc, ncod, 12)  :: primitive_forces
    REAL(dp), DIMENSION(:, :), POINTER       :: work_forces
    REAL(dp), INTENT(IN)                     :: eps_schwarz
    TYPE(x_hf_screening_type)                :: screening_parameter

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'forces4'

    INTEGER :: ipgf, jpgf, kpgf, la, lb, lc, ld, lpgf, offset_a, offset_a1, &
      offset_b, offset_b1, offset_c, offset_c1, offset_d, offset_d1

    primitive_forces = 0.0_dp
    DO ipgf = 1,npgfa
      offset_a1 = (ipgf-1)*ncoset(la_max)
      DO jpgf = 1,npgfb
        offset_b1 = (jpgf-1)*ncoset(lb_max)
        DO kpgf = 1,npgfc
          offset_c1 = (kpgf-1)*ncoset(lc_max)
          DO lpgf = 1,npgfd
            offset_d1 = (lpgf-1)*ncoset(ld_max)
            DO la = la_min,la_max
              offset_a = offset_a1 + ncoset(la-1)
              DO lb = lb_min, lb_max
                offset_b = offset_b1 + ncoset(lb-1)
                DO lc = lc_min, lc_max
                  offset_c = offset_c1 + ncoset(lc-1)
                  DO ld = ld_min, ld_max
                    offset_d = offset_d1 + ncoset(ld-1) 
                    !Build primitives
                    CALL evaluate_deriv_eri(deriv, ra, rb, rc, rd,&
                                           zeta(ipgf), zetb(jpgf), zetc(kpgf), zetd(lpgf),&
                                           la, lb, lc ,ld,&
                                           work_forces,&
                                           ncoa, ncob, ncoc, ncod, &
                                           primitive_forces, &
                                           offset_a, offset_b, offset_c, offset_d,&
                                           screening_parameter)

                  END DO !ld
                END DO !lc
              END DO !lb
            END DO !la

          END DO !lpgf
        END DO !kpgf
      END DO !jpgf
    END DO !ipgf
  
  END SUBROUTINE forces4
!***************************************************************************

!!***** four_center_integrals/derivatives_four_center [1.0] *
!!
!!   NAME
!!     derivatives_four_center 
!!
!!   FUNCTION
!!     computes four center derivatives for a full basis set and updates the
!!     forces%fock_4c arrays
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_env 
!!     - ks_matrix
!!     - energy
!!     - rho: density matrix
!!     - para_env
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     03.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE derivatives_four_center(qs_env,ks_matrix,energy,rho,hf_section,para_env,& 
                                     do_periodic,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: hf_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(IN)                      :: do_periodic
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'derivatives_four_center', &
      routineP = moduleN//':'//routineN

    INTEGER :: coord, handle, i_atom, iatom, ikind, iset, j_atom, jatom, &
      jkind, jset, k_atom, katom, kkind, kset, l_atom, latom, lkind, lset, &
      ma, max_am, max_k, max_set, maxsgf, mb, mc, md, natom, ncoa, ncob, &
      ncoc, ncod, nseta, nsetb, nsetc, nsetd, nsgfa_total, nsgfb_total, &
      nsgfc_total, nsgfd_total, nspins, pa, pa1, pa2, pb, pb1, pb2, pc, pc1, &
      pc2, pd, pd1, pd2, sgfa, sgfb, sgfc, sgfd, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of, &
                                                last_sgf_global
    INTEGER, DIMENSION(:), POINTER :: la_max, la_min, lb_max, lb_min, lc_max, &
      lc_min, ld_max, ld_min, ncgfa, ncgfb, ncgfc, ncgfd, npgfa, npgfb, &
      npgfc, npgfd, nsgfa, nsgfb, nsgfc, nsgfd
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb, &
                                                first_sgfc, first_sgfd
    LOGICAL                                  :: failure, screen_pmat_forces
    REAL(dp) :: box_size(3), eps_schwarz, fac, hf_fraction, max_val1, &
      max_val2, ra(3), rb(3), rc(3), rd(3), screen_max, temp1, temp1b, temp2
    REAL(dp), DIMENSION(:, :), POINTER       :: sphi_a, sphi_b, sphi_c, &
                                                sphi_d, work_forces, zeta, &
                                                zetb, zetc, zetd
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: screen_set, T1
    REAL(dp), DIMENSION(:, :, :, :, :), &
      POINTER                                :: primitive_forces
    REAL(KIND=dp)                            :: pmax_max
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: full_density, &
                                                full_density_beta, pmax_atom, &
                                                screen_atom, sparse_block
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: pmax_set
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), POINTER      :: screen_periodic
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_a, orb_basis_b, &
                                                orb_basis_c, orb_basis_d
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(real_block_node_type), POINTER      :: block_node
    TYPE(x_hf_screening_type)                :: screening_parameter

! many timings somewhat more clear

    CALL mp_sync(para_env%group)
    CALL timeset(routineN,"I",'',handle)

    CALL section_vals_val_get(hf_section,"FRACTION",r_val=hf_fraction,error=error)
    CALL section_vals_val_get(hf_section,"EPS_SCHWARZ_FORCES",r_val=eps_schwarz,error=error)
    CALL section_vals_val_get(hf_section,"SCREEN_PMAT_FORCES",l_val=screen_pmat_forces,error=error)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    error=error)
 
    natom = SIZE(particle_set,1)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxsgf=maxsgf,&
                             maxnset=max_set,&
                             maxlgto=max_am)

   
    CALL initialize_libderiv(deriv, max_am, error)

    screening_parameter = qs_env%x_data%screening_parameter

    ALLOCATE(work_forces(((max_am+1)*(max_am+2)/2)**4,12),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env=qs_env, force=force, error=error)

    ALLOCATE (atom_of_kind(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (kind_of(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure) 
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             atom_of_kind=atom_of_kind,&
                             kind_of=kind_of)
    DO ikind = 1,SIZE(force)
      force(ikind)%fock_4c(:,:) = 0.0_dp
    END DO

    nspins = qs_env%dft_control%nspins
 
    ALLOCATE(last_sgf_global(0:natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    last_sgf_global(0)=0
    DO iatom=1,natom
      atom_kind => particle_set(iatom)%atomic_kind
 
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_a)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_a,&
                             nsgf=nsgfa_total)
      last_sgf_global(iatom) = last_sgf_global(iatom-1)+nsgfa_total
    END DO
 
    ALLOCATE(full_density(1,1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL copy_local_sm_to_replicated_fm(rho%rho_ao(1)%matrix,full_density)
    CALL symmetrize_matrix(full_density,"upper_to_lower")
    CALL mp_sum(full_density,para_env%group)
    IF(nspins == 2) THEN
      ALLOCATE(full_density_beta(1,1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      CALL copy_local_sm_to_replicated_fm(rho%rho_ao(2)%matrix,full_density_beta)
      CALL symmetrize_matrix(full_density_beta,"upper_to_lower")
      CALL mp_sum(full_density_beta,para_env%group)
    ELSE
      NULLIFY(full_density_beta)
    END IF

    CALL compute_pmax_matrices(pmax_set,pmax_atom,pmax_max, &
                  full_density,full_density_beta,screen_pmat_forces,particle_set,atomic_kind_set,error)

    ! 
    ! notice that we use eps_schwarz == EPS_SCHWARZ_FORCES
    !
    max_k = qs_env%x_data%number_of_k_loops
    IF( do_periodic ) THEN
      box_size(1) = qs_env%cell%hmat(1,1)
      box_size(2) = qs_env%cell%hmat(2,2)
      box_size(3) = qs_env%cell%hmat(3,3)
      max_k = MIN(max_k, SIZE(qs_env%x_data%neighbor_cells))
      ALLOCATE(screen_periodic(natom,natom,max_set,max_set,SIZE(qs_env%x_data%neighbor_cells)))
      CALL comp_screen_mat_per(screen_set,screen_atom,screen_max,screen_periodic,particle_set,atomic_kind_set,&
                               screening_parameter,.FALSE.,&
                               neighbor_cells=qs_env%x_data%neighbor_cells,&
                               box_size=box_size,&
                               max_k=max_k,error=error)
    ELSE
      CALL compute_screening_matrices(screen_set,screen_atom,screen_max,particle_set,atomic_kind_set,&
                                      screening_parameter,.FALSE.,error=error)
    END IF

    DO iatom=1,ks_matrix(1)%matrix%nblock_row
      block_node => first_block_node(ks_matrix(1)%matrix,iatom)
      atom_kind => particle_set(iatom)%atomic_kind

      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_a)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_a,&
                             lmax=la_max,&
                             lmin=la_min,&
                             npgf=npgfa,&
                             nset=nseta,&
                             zet=zeta,&
                             nsgf_set=nsgfa,&
                             ncgf_set=ncgfa,&
                             first_sgf=first_sgfa,&
                             sphi=sphi_a,&
                             nsgf=nsgfa_total)
      ra=particle_set(iatom)%r(:)
      
      i_atom = atom_of_kind(iatom)
      ikind = kind_of(iatom)
      pa1 = last_sgf_global(iatom) - nsgfa_total

      DO jatom = 1,natom
        atom_kind => particle_set(jatom)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atom_kind,&
                             orb_basis_set=orb_basis_b)
        CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                               lmax=lb_max,&
                               lmin=lb_min,&
                               npgf=npgfb,&
                               nset=nsetb,&
                               zet=zetb,&
                               nsgf_set=nsgfb,&
                               ncgf_set=ncgfb,&
                               first_sgf=first_sgfb,&
                               sphi=sphi_b,&
                               nsgf=nsgfb_total)
        rb=particle_set(jatom)%r(:)

        j_atom = atom_of_kind(jatom)
        jkind = kind_of(jatom)
        pb1 = last_sgf_global(jatom) - nsgfb_total
     
        DO katom = 1,natom
          atom_kind => particle_set(katom)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atom_kind,&
                               orb_basis_set=orb_basis_c)
          CALL get_gto_basis_set(gto_basis_set=orb_basis_c,&
                                 lmax=lc_max,&
                                 lmin=lc_min,&
                                 npgf=npgfc,&
                                 nset=nsetc,&
                                 zet=zetc,&
                                 nsgf_set=nsgfc,&
                                 ncgf_set=ncgfc,&
                                 first_sgf=first_sgfc,&
                                 sphi=sphi_c,&
                                 nsgf=nsgfc_total)
          rc=particle_set(katom)%r(:)
          
          k_atom = atom_of_kind(katom)
          kkind = kind_of(katom)
          pc1 = last_sgf_global(katom) - nsgfc_total
           
          block_node => first_block_node(ks_matrix(1)%matrix,iatom)
          DO WHILE (ASSOCIATED(block_node))
            CALL get_block_node(block_node=block_node,&
                                block_col=latom,&
                                BLOCK=sparse_block)

            block_node => next_block_node(block_node)
          
            IF((iatom==jatom .AND. iatom==katom .AND. iatom==latom)) CYCLE 

            atom_kind => particle_set(latom)%atomic_kind
            CALL get_atomic_kind(atomic_kind=atom_kind,&
                                 orb_basis_set=orb_basis_d)
            CALL get_gto_basis_set(gto_basis_set=orb_basis_d,&
                                   lmax=ld_max,&
                                   lmin=ld_min,&
                                   npgf=npgfd,&
                                   nset=nsetd,&
                                   zet=zetd,&
                                   nsgf_set=nsgfd,&
                                   ncgf_set=ncgfd,&
                                   first_sgf=first_sgfd,&
                                   sphi=sphi_d,&
                                   nsgf=nsgfd_total)
 
            rd=particle_set(latom)%r(:)
            
            l_atom = atom_of_kind(latom)
            lkind = kind_of(latom)
            pd1 = last_sgf_global(latom) - nsgfd_total
             
            IF(nspins == 1) THEN
              fac = 0.25_dp * hf_fraction
            ELSE
              fac = 0.5_dp * hf_fraction
            END IF
            IF(iatom<latom) fac = fac * 2.0_dp
            DO iset = 1, nseta
              ncoa = npgfa(iset)*ncoset(la_max(iset))
              sgfa = first_sgfa(1,iset)
              DO jset = 1,nsetb
                ncob = npgfb(jset)*ncoset(lb_max(jset))
                sgfb = first_sgfb(1,jset)
                max_val1 = screen_set(iatom,jatom,iset,jset)
                DO kset = 1,nsetc
                  ncoc = npgfc(kset)*ncoset(lc_max(kset))
                  sgfc = first_sgfc(1,kset)
                  DO lset = 1,nsetd
                    ncod = npgfd(lset)*ncoset(ld_max(lset)) 
                    sgfd = first_sgfd(1,lset)

                    max_val2 = max_val1 * screen_set(katom,latom,kset,lset) * &
                               pmax_set(jset,kset,jatom,katom) * pmax_set(iset,lset,iatom,latom)
                    IF(max_val2<eps_schwarz)  CYCLE
                    
                    ALLOCATE(primitive_forces(ncoa,ncob,ncoc,ncod,12),STAT=stat)
                    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                    primitive_forces = 0.0_dp
                    
                    ALLOCATE(T1(ncoa,ncob,ncoc,ncod),STAT=stat)
                    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                    IF( .NOT. do_periodic ) THEN 
                      CALL forces4(ra, rb, rc, rd,&
                                   zeta(:,iset), zetb(:,jset), zetc(:,kset), zetd(:,lset),&
                                   la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                                   lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                                   npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset),&
                                   ncoa, ncob, ncoc, ncod,&
                                   primitive_forces,&
                                   work_forces,& 
                                   eps_schwarz,&
                                   screening_parameter)
                    ELSE
                      CALL forces4_periodic(ra, rb, rc, rd,&
                                            zeta(:,iset), zetb(:,jset), zetc(:,kset), zetd(:,lset),&
                                            la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                                            lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                                            npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset),&
                                            ncoa, ncob, ncoc, ncod,&
                                            primitive_forces,&
                                            work_forces,&
                                            screening_parameter, box_size, &
                                            qs_env%x_data%neighbor_cells, screen_periodic(iatom,jatom,iset,jset,:),&
                                            screen_periodic(katom,latom,kset,lset,:), eps_schwarz, max_k)

                    END IF
                    DO coord = 1,12
                      CALL dgemm("T","N",ncob*ncoc*ncod,nsgfa(iset),ncoa,&
                                 1.0_dp, primitive_forces(1,1,1,1,coord),ncoa,&
                                 sphi_a(1,sgfa), SIZE(sphi_a,1),&
                                 0.0_dp, T1(1,1,1,1),ncob*ncoc*ncod)

  
                      CALL dgemm("T","N",ncoa*ncoc*ncod,nsgfb(jset),ncob,&
                                 1.0_dp, T1(1,1,1,1),ncob,&
                                 sphi_b(1,sgfb), SIZE(sphi_b,1),&
                                 0.0_dp, primitive_forces(1,1,1,1,coord),ncoa*ncoc*ncod)

  
                      CALL dgemm("T","N",ncoa*ncob*ncod,nsgfc(kset),ncoc,&
                                 1.0_dp, primitive_forces(1,1,1,1,coord),ncoc,&
                                 sphi_c(1,sgfc), SIZE(sphi_c,1),&
                                 0.0_dp, T1(1,1,1,1),ncoa*ncob*ncod)

  
                      CALL dgemm("T","N",ncoa*ncob*ncoc,nsgfd(lset),ncod,&
                                 1.0_dp, T1(1,1,1,1),ncod,&
                                 sphi_d(1,sgfd), SIZE(sphi_d,1),&
                                 0.0_dp, primitive_forces(1,1,1,1,coord),ncoa*ncob*ncoc)
                    END DO              
                        

                    DO mb = 1, nsgfb(jset)
                      pb2 = mb+sgfb-1
                      pb = pb1 + pb2
                      DO mc = 1, nsgfc(kset)
                        pc2 = mc+sgfc-1
                        pc = pc1 + pc2
                        temp1 = full_density(pc,pb) * fac
                        IF( nspins == 2 ) temp1b = full_density_beta(pc,pb) * fac
                        DO md = 1,nsgfd(lset)
                          pd2 = md+sgfd-1
                          pd = pd1 + pd2
                          DO ma = 1,nsgfa(iset) 
                            pa2 = ma+sgfa-1
                            pa = pa1 + pa2
                            temp2 = temp1 * full_density(pa,pd)
                    
                            force(ikind)%fock_4c(1,i_atom) = force(ikind)%fock_4c(1,i_atom) - &
                                                             temp2 * primitive_forces(ma,mb,mc,md,1)
                            force(ikind)%fock_4c(2,i_atom) = force(ikind)%fock_4c(2,i_atom) - &
                                                             temp2 * primitive_forces(ma,mb,mc,md,2) 
                            force(ikind)%fock_4c(3,i_atom) = force(ikind)%fock_4c(3,i_atom) - &
                                                             temp2 * primitive_forces(ma,mb,mc,md,3)
        
                            force(jkind)%fock_4c(1,j_atom) = force(jkind)%fock_4c(1,j_atom) - &
                                                             temp2 * primitive_forces(ma,mb,mc,md,4) 
                            force(jkind)%fock_4c(2,j_atom) = force(jkind)%fock_4c(2,j_atom) - &
                                                             temp2 * primitive_forces(ma,mb,mc,md,5)
                            force(jkind)%fock_4c(3,j_atom) = force(jkind)%fock_4c(3,j_atom) - &
                                                             temp2 * primitive_forces(ma,mb,mc,md,6)

                            force(kkind)%fock_4c(1,k_atom) = force(kkind)%fock_4c(1,k_atom) - &
                                                             temp2 * primitive_forces(ma,mb,mc,md,7)
                            force(kkind)%fock_4c(2,k_atom) = force(kkind)%fock_4c(2,k_atom) - &
                                                             temp2 * primitive_forces(ma,mb,mc,md,8)
                            force(kkind)%fock_4c(3,k_atom) = force(kkind)%fock_4c(3,k_atom) - &
                                                             temp2 * primitive_forces(ma,mb,mc,md,9)

                            force(lkind)%fock_4c(1,l_atom) = force(lkind)%fock_4c(1,l_atom) - &
                                                             temp2 * primitive_forces(ma,mb,mc,md,10)
                            force(lkind)%fock_4c(2,l_atom) = force(lkind)%fock_4c(2,l_atom) - &
                                                             temp2 * primitive_forces(ma,mb,mc,md,11)
                            force(lkind)%fock_4c(3,l_atom) = force(lkind)%fock_4c(3,l_atom) - &
                                                             temp2 * primitive_forces(ma,mb,mc,md,12)
                            IF(nspins == 2) THEN
                              temp2 = temp1b*full_density_beta(pa,pd)
                              force(ikind)%fock_4c(1,i_atom) = force(ikind)%fock_4c(1,i_atom) - &
                                                               temp2 * primitive_forces(ma,mb,mc,md,1)
                              force(ikind)%fock_4c(2,i_atom) = force(ikind)%fock_4c(2,i_atom) - &
                                                               temp2 * primitive_forces(ma,mb,mc,md,2) 
                              force(ikind)%fock_4c(3,i_atom) = force(ikind)%fock_4c(3,i_atom) - &
                                                               temp2 * primitive_forces(ma,mb,mc,md,3)

                              force(jkind)%fock_4c(1,j_atom) = force(jkind)%fock_4c(1,j_atom) - &
                                                               temp2 * primitive_forces(ma,mb,mc,md,4) 
                              force(jkind)%fock_4c(2,j_atom) = force(jkind)%fock_4c(2,j_atom) - &
                                                               temp2 * primitive_forces(ma,mb,mc,md,5)
                              force(jkind)%fock_4c(3,j_atom) = force(jkind)%fock_4c(3,j_atom) - &
                                                               temp2 * primitive_forces(ma,mb,mc,md,6)

                              force(kkind)%fock_4c(1,k_atom) = force(kkind)%fock_4c(1,k_atom) - &
                                                               temp2 * primitive_forces(ma,mb,mc,md,7)
                              force(kkind)%fock_4c(2,k_atom) = force(kkind)%fock_4c(2,k_atom) - &
                                                               temp2 * primitive_forces(ma,mb,mc,md,8)
                              force(kkind)%fock_4c(3,k_atom) = force(kkind)%fock_4c(3,k_atom) - &
                                                               temp2 * primitive_forces(ma,mb,mc,md,9)

                              force(lkind)%fock_4c(1,l_atom) = force(lkind)%fock_4c(1,l_atom) - &
                                                               temp2 * primitive_forces(ma,mb,mc,md,10)
                              force(lkind)%fock_4c(2,l_atom) = force(lkind)%fock_4c(2,l_atom) - &
                                                               temp2 * primitive_forces(ma,mb,mc,md,11)
                              force(lkind)%fock_4c(3,l_atom) = force(lkind)%fock_4c(3,l_atom) - &
                                                               temp2 * primitive_forces(ma,mb,mc,md,12)
                            END IF 
 
                          END DO !ma
                        END DO !mb
                      END DO !mc
                    END DO !md

                    DEALLOCATE(primitive_forces,T1,STAT=stat)
                    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                  END DO !lset
                END DO !ksetm
               END DO !jset
            END DO !iset
         
          END DO !latom
        END DO !katom
      END DO !jatom
    END DO !iatom

    IF( do_periodic ) THEN
      DEALLOCATE(screen_periodic, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    DEALLOCATE(last_sgf_global,screen_set, screen_atom, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(full_density,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(nspins==2) THEN
      DEALLOCATE(full_density_beta,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE(work_forces, atom_of_kind,kind_of,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(pmax_set,pmax_atom,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    !
    ! this wraps to free_libint, but currently causes a segfault
    ! as a result, we don't call it, but some memory remains allocated
    !  
    ! CALL terminate_libderiv(deriv)

    CALL mp_sync(para_env%group)

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE derivatives_four_center
!***************************************************************************

!!***** four_center_integrals/integrate_four_center_symm [1.0] *
!!
!!   NAME
!!     integrate_four_center_symm
!!
!!   FUNCTION
!!     computes four center integrals for a full basis set and updates the
!!     Kohn-Sham-Matrix and energy. Uses all 8 eri symmetries
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_env 
!!     - ks_matrix
!!     - energy
!!     - rho: density matrix
!!     - para_env
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     06.2007 created [Manuel Guidon]
!!     08.2007 optimized load balance [Manuel Guidon] 
!!
!!****************************************************************************

  SUBROUTINE integrate_four_center_symm(qs_env,ks_matrix,energy,rho,hf_section,para_env,& 
                                   geometry_did_change,do_periodic,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: hf_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL                                  :: geometry_did_change, &
                                                do_periodic
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'integrate_four_center_symm', &
      routineP = moduleN//':'//routineN

    INTEGER :: bin, box_idx1, box_idx2, box_idx3, box_idx4, current_counter, &
      handle, handle_dist_ks, handle_getP, handle_main, iatom, &
      integral_counter, iset, iw, jatom, jset, katom, kset, latom, lset, ma, &
      max_am, max_k, max_set, maxsgf, mb, mc, md, my_bin_size, my_istart, &
      my_jstart, my_kstart, my_lstart, natom, nco_max, ncoa, ncob, ncoc, &
      ncod, nseta, nsetb, nsetc, nsetd, nsgfa_total, nsgfb_total, &
      nsgfc_total, nsgfd_total, nspins, number_of_precalc_integrals, p_index, &
      pa, pa1, pa2, pb, pb1, pb2, pc, pc1, pc2, pd, pd1, pd2, pgf_radius_box, &
      sgfa, sgfb, sgfc, sgfd, stat, total_integrals_to_store, &
      total_integrals_to_store_aver, total_integrals_to_store_max
    INTEGER(KIND=int_8)                      :: my_current_counter
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: last_sgf_global
    INTEGER, DIMENSION(:), POINTER :: la_max, la_min, lb_max, lb_min, lc_max, &
      lc_min, ld_max, ld_min, npgfa, npgfb, npgfc, npgfd, nsgfa, nsgfb, &
      nsgfc, nsgfd
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb, &
                                                first_sgfc, first_sgfd, &
                                                is_assoc_atomic_block
    LOGICAL                                  :: buffer_overflow, &
                                                do_p_screening, failure, &
                                                first_j_loop, first_k_loop, &
                                                first_l_loop
    REAL(dp) :: box_size(3), ene_x_aa, ene_x_bb, eps_schwarz, fac, &
      hf_fraction, K, Kbox, max_val1, max_val2, memory_usage, omega_box, &
      ra(3), rb(3), Rbox, rc(3), rd(3), screen_max, symm_fac, tmp
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: primitive_integrals, T1
    REAL(dp), DIMENSION(:, :), POINTER       :: full_ks, full_ks_beta, &
                                                screen_atom, sphi_a, sphi_b, &
                                                sphi_c, sphi_d, zeta, zetb, &
                                                zetc, zetd
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: screen_set
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: full_density, &
                                                full_density_beta, &
                                                sparse_block, &
                                                sparse_block_beta
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), POINTER      :: screen_periodic
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_a, orb_basis_b, &
                                                orb_basis_c, orb_basis_d
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(real_block_node_type), POINTER      :: block_node, block_node_beta
    TYPE(x_hf_screening_type)                :: screening_parameter

! many timings somewhat more clear

    CALL mp_sync(para_env%group)
    CALL timeset(routineN,"I",'',handle)

    CALL section_vals_val_get(hf_section,"FRACTION",r_val=hf_fraction,error=error)
    CALL section_vals_val_get(hf_section,"EPS_SCHWARZ",r_val=eps_schwarz,error=error)
   
    number_of_precalc_integrals = qs_env%x_data%number_of_elements
    screening_parameter = qs_env%x_data%screening_parameter
    integral_counter = 1
    memory_usage = 0.0_dp
    buffer_overflow = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
   
    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    error=error)
 
    natom = SIZE(particle_set,1)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxsgf=maxsgf,&
                             maxnset=max_set,&
                             maxlgto=max_am)
    
    CALL initialize_libint(lib, max_am, error)
    
    nspins = qs_env%dft_control%nspins
 
    ALLOCATE(last_sgf_global(0:natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    last_sgf_global(0)=0
    DO iatom=1,natom
      atom_kind => particle_set(iatom)%atomic_kind
 
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_a)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_a,&
                             nsgf=nsgfa_total)
      last_sgf_global(iatom) = last_sgf_global(iatom-1)+nsgfa_total
    END DO
 
    CALL timeset(routineN//"_getP",handle_getP)

    ALLOCATE(full_density(1,1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL copy_local_sm_to_replicated_fm(rho%rho_ao(1)%matrix,full_density)
    CALL symmetrize_matrix(full_density,"upper_to_lower")
    CALL mp_sum(full_density,para_env%group)
    IF(nspins == 2) THEN
      ALLOCATE(full_density_beta(1,1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      CALL copy_local_sm_to_replicated_fm(rho%rho_ao(2)%matrix,full_density_beta)
      CALL symmetrize_matrix(full_density_beta,"upper_to_lower")
      CALL mp_sum(full_density_beta,para_env%group)
    END IF

    CALL timestop(handle_getP)

    CALL section_vals_val_get(hf_section,"SCREEN_ON_INITIAL_P",l_val=do_p_screening,error=error)
    IF( do_p_screening) THEN
      IF(geometry_did_change ) THEN
        CALL update_pmax_matrix(qs_env%x_data%initial_p, &
                    full_density,full_density_beta,&
                    particle_set,atomic_kind_set,error)
      END IF
    END IF

    ALLOCATE(full_ks(SIZE(full_density,1),SIZE(full_density,2)))
    full_ks = 0.0_dp
    IF(nspins==2) THEN
      ALLOCATE(full_ks_beta(SIZE(full_density,1),SIZE(full_density,2)))
      full_ks_beta = 0.0_dp
    END IF
   
    max_k = qs_env%x_data%number_of_k_loops
    IF( do_periodic ) THEN
      box_size(1) = qs_env%cell%hmat(1,1)
      box_size(2) = qs_env%cell%hmat(2,2)
      box_size(3) = qs_env%cell%hmat(3,3)
      max_k = MIN(max_k, SIZE(qs_env%x_data%neighbor_cells))
      ALLOCATE(screen_periodic(natom,natom,max_set,max_set,SIZE(qs_env%x_data%neighbor_cells)))
      CALL comp_screen_mat_per(screen_set,screen_atom,screen_max,screen_periodic,particle_set,atomic_kind_set,&
                               screening_parameter,.NOT.geometry_did_change,&
                               qs_env%x_data%x_screening_buffer_periodic,&
                               qs_env%x_data%neighbor_cells,box_size,max_k,error)
      IF(geometry_did_change) THEN
        CALL screen_farfield_per(qs_env%x_data%far_field_boxes, &
                                 particle_set, atomic_kind_set, screening_parameter, &
                                 qs_env%x_data%neighbor_cells, box_size, error)
       END IF

    ELSE
      CALL compute_screening_matrices(screen_set,screen_atom,screen_max,particle_set,atomic_kind_set,&
                                    screening_parameter,.NOT.geometry_did_change,qs_env%x_data%x_screening_buffer, error)
      IF(geometry_did_change) THEN
        CALL screen_farfield(qs_env%x_data%far_field_boxes,&
                            particle_set, atomic_kind_set, screening_parameter,&
                            error)
      END IF
    END IF
    CALL timeset(routineN//"_main",handle_main)

    total_integrals_to_store=0
    memory_usage = 0.0_dp
    IF(nspins == 1) THEN
      fac = 0.5_dp * hf_fraction
    ELSE
      fac = 1.0_dp * hf_fraction
    END IF
    
    ALLOCATE(is_assoc_atomic_block(natom,natom))
    is_assoc_atomic_block = 0
    DO iatom=1,ks_matrix(1)%matrix%nblock_row
      block_node => first_block_node(ks_matrix(1)%matrix,iatom)
      DO WHILE (ASSOCIATED(block_node))
        CALL get_block_node(block_node=block_node,&
                            block_col=jatom,&
                            BLOCK=sparse_block)

        is_assoc_atomic_block(iatom,jatom) = 1
        is_assoc_atomic_block(jatom,iatom) = 1
        block_node => next_block_node(block_node)
      END DO
    END DO
    
    CALL mp_sum(is_assoc_atomic_block,para_env%group)

    IF( geometry_did_change ) THEN
      IF( do_p_screening) THEN
        CALL load_balance_energy(qs_env%x_data,eps_schwarz,particle_set,atomic_kind_set,para_env,&
                                 screening_parameter,screen_set,screen_atom,screen_max,&
                                 is_assoc_atomic_block,ks_matrix,do_periodic,pmax_set=qs_env%x_data%initial_p,error=error)
      ELSE
        CALL load_balance_energy(qs_env%x_data,eps_schwarz,particle_set,atomic_kind_set,para_env,&
                                 screening_parameter,screen_set,screen_atom,screen_max,&
                                 is_assoc_atomic_block,ks_matrix,do_periodic,error=error)
      END IF
    END IF
    ! precompute maximum nco and allocate scratch
    nco_max=0
    DO iatom=1,natom
      atom_kind => particle_set(iatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_a)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_a, lmax=la_max, lmin=la_min,&
                             npgf=npgfa, nset=nseta, zet=zeta, nsgf_set=nsgfa,&
                             first_sgf=first_sgfa, sphi=sphi_a, nsgf=nsgfa_total)
      DO iset = 1, nseta
         nco_max = MAX(nco_max,npgfa(iset)*ncoset(la_max(iset)))
      ENDDO
    ENDDO
    ALLOCATE(primitive_integrals(nco_max**4),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    primitive_integrals=0.0_dp
    ALLOCATE(T1(nco_max**4),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    T1=0.0_dp

    ! the full loop
    my_bin_size = SIZE(qs_env%x_data%distribution_energy)
    first_j_loop = .TRUE.
    first_l_loop = .TRUE.
    first_k_loop = .TRUE.
    outer: DO bin=1,my_bin_size
      my_istart = qs_env%x_data%distribution_energy(bin)%istart
      my_jstart = qs_env%x_data%distribution_energy(bin)%jstart
      my_kstart = qs_env%x_data%distribution_energy(bin)%kstart
      my_lstart = qs_env%x_data%distribution_energy(bin)%lstart
      my_current_counter = 0
      IF(my_istart == 0 .OR. my_jstart == 0 .OR. my_kstart == 0 .OR. my_lstart == 0 .OR. &
         qs_env%x_data%distribution_energy(bin)%number_of_integrals == 0) CYCLE
      DO iatom=my_istart,natom
        atom_kind => particle_set(iatom)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atom_kind,&
                             orb_basis_set=orb_basis_a)
        CALL get_gto_basis_set(gto_basis_set=orb_basis_a, lmax=la_max, lmin=la_min,&
                               npgf=npgfa, nset=nseta, zet=zeta, nsgf_set=nsgfa,&
                               first_sgf=first_sgfa, sphi=sphi_a, nsgf=nsgfa_total)
        ra=particle_set(iatom)%r(:)
        pa1 = last_sgf_global(iatom) - nsgfa_total
  
        IF(first_j_loop) jatom=my_jstart-1
        IF(.NOT. first_j_loop) jatom = iatom-1
        DO WHILE(jatom+1<=natom)
          jatom = jatom + 1
          IF (screen_atom(iatom,jatom) * screen_max < eps_schwarz) CYCLE
          atom_kind => particle_set(jatom)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atom_kind,&
                               orb_basis_set=orb_basis_b)
          CALL get_gto_basis_set(gto_basis_set=orb_basis_b, lmax=lb_max, lmin=lb_min,&
                                 npgf=npgfb, nset=nsetb, zet=zetb, nsgf_set=nsgfb,&
                                 first_sgf=first_sgfb, sphi=sphi_b, nsgf=nsgfb_total)
          rb=particle_set(jatom)%r(:)
          pb1 = last_sgf_global(jatom) - nsgfb_total

        
          IF(first_k_loop) katom=my_kstart-1
          IF(.NOT. first_k_loop) katom = 0
          DO WHILE(katom+1<=natom)
            katom = katom + 1
            atom_kind => particle_set(katom)%atomic_kind
            CALL get_atomic_kind(atomic_kind=atom_kind,&
                                 orb_basis_set=orb_basis_c)
            CALL get_gto_basis_set(gto_basis_set=orb_basis_c, lmax=lc_max, lmin=lc_min,&
                                   npgf=npgfc, nset=nsetc, zet=zetc, nsgf_set=nsgfc,&
                                   first_sgf=first_sgfc, sphi=sphi_c, nsgf=nsgfc_total)
            rc=particle_set(katom)%r(:)
            pc1 = last_sgf_global(katom) - nsgfc_total

             IF(first_l_loop) THEN
                latom=my_lstart-1
             ELSE
                latom = katom-1
             END IF
             DO WHILE(latom+1<=natom)
              latom = latom + 1
              IF(.NOT. (is_assoc_atomic_block(iatom,latom)>=1 .AND. &
                        is_assoc_atomic_block(iatom,katom)>=1 .AND. &
                        is_assoc_atomic_block(jatom,katom)>=1 .AND. &
                        is_assoc_atomic_block(jatom,latom)>=1 ) ) CYCLE
              IF(katom+latom<=iatom+jatom)  THEN
                IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE
                IF (screen_atom(iatom,jatom) * screen_atom(katom,latom) < eps_schwarz) CYCLE
  
                 my_current_counter = my_current_counter + 1
                 IF(my_current_counter > qs_env%x_data%distribution_energy(bin)%number_of_integrals) THEN
                   first_j_loop = .TRUE.
                   first_l_loop = .TRUE.
                   first_k_loop = .TRUE.
                   CYCLE outer
                 END IF

                !calculate symmetry_factor
                symm_fac = 0.5_dp
                IF(iatom==jatom) symm_fac = symm_fac*2.0_dp
                IF(katom==latom) symm_fac = symm_fac*2.0_dp
                IF(iatom==katom .AND. jatom==latom .AND. iatom/=jatom .AND. katom/=latom) symm_fac = symm_fac*2.0_dp
                IF(iatom==katom .AND. iatom==jatom .AND. katom==latom) symm_fac = symm_fac*2.0_dp
                symm_fac = 1.0_dp / symm_fac

                atom_kind => particle_set(latom)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atom_kind,&
                                     orb_basis_set=orb_basis_d)
                CALL get_gto_basis_set(gto_basis_set=orb_basis_d, lmax=ld_max, lmin=ld_min,&
                                       npgf=npgfd, nset=nsetd, zet=zetd, nsgf_set=nsgfd,&
                                       first_sgf=first_sgfd, sphi=sphi_d, nsgf=nsgfd_total)
                rd=particle_set(latom)%r(:)
                pd1 = last_sgf_global(latom) - nsgfd_total
           
                DO iset = 1, nseta
                  ncoa = npgfa(iset)*ncoset(la_max(iset))
                  sgfa = first_sgfa(1,iset)
                  DO jset = 1,nsetb
                    ncob = npgfb(jset)*ncoset(lb_max(jset))
                    sgfb = first_sgfb(1,jset)
                    max_val1 = screen_set(iatom,jatom,iset,jset)
                    IF (max_val1 * screen_atom(katom,latom) < eps_schwarz) CYCLE
                    DO kset = 1,nsetc
                      ncoc = npgfc(kset)*ncoset(lc_max(kset))
                      sgfc = first_sgfc(1,kset)
                      DO lset = 1,nsetd
                        ncod = npgfd(lset)*ncoset(ld_max(lset)) 
                        sgfd = first_sgfd(1,lset)
                        max_val2 = max_val1 * screen_set(katom,latom,kset,lset)
                        IF( do_p_screening) THEN
                          max_val2 = max_val2 * MAX(qs_env%x_data%initial_p(iset,kset,iatom,katom), &
                                                    qs_env%x_data%initial_p(jset,lset,jatom,latom), &
                                                    qs_env%x_data%initial_p(iset,lset,iatom,latom), &
                                                    qs_env%x_data%initial_p(kset,jset,katom,jatom))
                        END IF
                        IF(max_val2<eps_schwarz) CYCLE
                        box_idx1 = (iatom-1)*natom+jatom
                        box_idx2 = (iset-1)*max_set+jset
                        box_idx3 = (katom-1)*natom+latom
                        box_idx4 = (kset-1)*max_set+lset
                        Rbox = (qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%center(1) - &
                                qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%center(1))**2 + &
                               (qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%center(2) - &
                                qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%center(2))**2 + &
                               (qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%center(3) - &
                                qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%center(3))**2
                        Rbox = SQRT(Rbox)
                        Rbox = Rbox - qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%R - &
                                      qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%R
                        pgf_radius_box = INT(1.0_dp/ &
                              (SQRT(2.0_dp*qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%alpha)*erfc(epsilon)))
                        pgf_radius_box = pgf_radius_box + INT(1.0_dp/ & 
                              (SQRT(2.0_dp*qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%alpha)*erfc(epsilon))) + 2
                        IF(Rbox>pgf_radius_box) THEN
                          max_val2 = 0.0_dp
                          SELECT CASE(screening_parameter%screening_type)
                            CASE(do_hf_screening_none)
                              Kbox = qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                     qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%K
                              max_val2 = 1.0_dp/Rbox*Kbox
                              IF( do_p_screening) THEN
                                max_val2 = max_val2 * MAX(qs_env%x_data%initial_p(iset,kset,iatom,katom), &
                                                          qs_env%x_data%initial_p(jset,lset,jatom,latom), &
                                                          qs_env%x_data%initial_p(iset,lset,iatom,latom), &
                                                          qs_env%x_data%initial_p(kset,jset,katom,jatom))
                              END IF
                            CASE(do_hf_screening_short)
                              omega_box = 1.0_dp/qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%alpha + &
                                          1.0_dp/qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%alpha + &
                                          1.0_dp/(screening_parameter%omega**2)
                              omega_box = 1.0_dp/omega_box
                              Kbox = qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                     qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%K
                              max_val2 = erfc(omega_box*Rbox)/Rbox*Kbox
                              IF( do_p_screening) THEN
                                max_val2 = max_val2 * MAX(qs_env%x_data%initial_p(iset,kset,iatom,katom), &
                                                          qs_env%x_data%initial_p(jset,lset,jatom,latom), &
                                                          qs_env%x_data%initial_p(iset,lset,iatom,latom), &
                                                          qs_env%x_data%initial_p(kset,jset,katom,jatom))
                              END IF
                          END SELECT
                          IF(max_val2<eps_schwarz) CYCLE
                        END IF

                        current_counter = nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)
                        
                        total_integrals_to_store=total_integrals_to_store+current_counter
                        buffer_overflow = integral_counter + current_counter - 1 > number_of_precalc_integrals
                        IF(.NOT.buffer_overflow .AND. .NOT.geometry_did_change) THEN 
                          DO md = 1,nsgfd(lset)
                            pd2 = md+sgfd-1
                            pd = pd1 + pd2
                            DO mc = 1,nsgfc(kset)
                              pc2 = mc+sgfc-1
                              pc = pc1 + pc2
                              DO mb = 1,nsgfb(jset)
                                pb2 = mb+sgfb-1
                                pb = pb1 + pb2
                                DO ma = 1,nsgfa(iset)
                                  pa2 = ma+sgfa-1
                                  pa = pa1 + pa2
                                  tmp = fac*symm_fac*qs_env%x_data%x_buffer(integral_counter)
                                  full_ks(pa,pd) = full_ks(pa,pd)-tmp*full_density(pb,pc)
                                  full_ks(pa,pc) = full_ks(pa,pc)-tmp*full_density(pb,pd)
                                  full_ks(pb,pc) = full_ks(pb,pc)-tmp*full_density(pa,pd)
                                  full_ks(pb,pd) = full_ks(pb,pd)-tmp*full_density(pa,pc)
                                  IF( nspins == 2 ) THEN
                                    full_ks_beta(pa,pd) = full_ks_beta(pa,pd)-tmp*full_density_beta(pb,pc)
                                    full_ks_beta(pa,pc) = full_ks_beta(pa,pc)-tmp*full_density_beta(pb,pd)
                                    full_ks_beta(pb,pc) = full_ks_beta(pb,pc)-tmp*full_density_beta(pa,pd)
                                    full_ks_beta(pb,pd) = full_ks_beta(pb,pd)-tmp*full_density_beta(pa,pc)
                                  END IF
                                  integral_counter = integral_counter + 1
                                END DO
                              END DO
                            END DO
                          END DO
                        END IF
  
                        IF(geometry_did_change .OR. buffer_overflow) THEN
                          IF( .NOT.do_periodic ) THEN 
                            CALL coulomb4(ra, rb, rc, rd,&
                                          zeta(1:npgfa(iset),iset), zetb(1:npgfb(jset),jset), &
                                          zetc(1:npgfc(kset),kset), zetd(1:npgfd(lset),lset),&
                                          la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                                          lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                                          npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset),&
                                          ncoa,ncob,ncoc,ncod,&
                                          primitive_integrals(1),&
                                          screening_parameter)
                          ELSE
                            CALL coulomb4_periodic(ra, rb, rc, rd,&
                                                   zeta(1:npgfa(iset),iset), zetb(1:npgfb(jset),jset), &
                                                   zetc(1:npgfc(kset),kset), zetd(1:npgfd(lset),lset),&
                                                   la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                                                   lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                                                   npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset),&
                                                   ncoa,ncob,ncoc,ncod,&
                                                   primitive_integrals(1),&
                                                   screening_parameter, box_size, &
                                                   qs_env%x_data%neighbor_cells, screen_periodic(iatom,jatom,iset,jset,:),&
                                                   screen_periodic(katom,latom,kset,lset,:), eps_schwarz, max_k)
                          END IF
  
                          CALL dgemm("T","N",ncob*ncoc*ncod,nsgfa(iset),ncoa,&
                                     1.0_dp, primitive_integrals(1),ncoa,&
                                     sphi_a(1,sgfa), SIZE(sphi_a,1),&
                                     0.0_dp, T1(1),ncob*ncoc*ncod)
  
                          CALL dgemm("T","N",nsgfa(iset)*ncoc*ncod,nsgfb(jset),ncob,&
                                     1.0_dp, T1(1),ncob,&
                                     sphi_b(1,sgfb), SIZE(sphi_b,1),&
                                     0.0_dp, primitive_integrals(1),nsgfa(iset)*ncoc*ncod)
  
                          CALL dgemm("T","N",nsgfa(iset)*nsgfb(jset)*ncod,nsgfc(kset),ncoc,&
                                     1.0_dp, primitive_integrals(1),ncoc,&
                                     sphi_c(1,sgfc), SIZE(sphi_c,1),&
                                     0.0_dp, T1(1),nsgfa(iset)*nsgfb(jset)*ncod)
  
                          CALL dgemm("T","N",nsgfa(iset)*nsgfb(jset)*nsgfc(kset),nsgfd(lset),ncod,&
                                     1.0_dp, T1(1),ncod,&
                                     sphi_d(1,sgfd), SIZE(sphi_d,1),&
                                     0.0_dp, primitive_integrals(1),nsgfa(iset)*nsgfb(jset)*nsgfc(kset))

                          p_index=0 
                          DO md = 1,nsgfd(lset)
                            pd2 = md+sgfd-1
                            pd = pd1 + pd2
                            DO mc = 1,nsgfc(kset)
                              pc2 = mc+sgfc-1
                              pc = pc1 + pc2
                              DO mb = 1,nsgfb(jset)
                                pb2 = mb+sgfb-1
                                pb = pb1 + pb2
                                DO ma = 1,nsgfa(iset)
                                  pa2 = ma+sgfa-1
                                  pa = pa1 + pa2
                                  p_index=p_index+1
                                  tmp = fac*symm_fac*primitive_integrals(p_index)
                                  full_ks(pa,pd) = full_ks(pa,pd)-tmp*full_density(pb,pc)
                                  full_ks(pa,pc) = full_ks(pa,pc)-tmp*full_density(pb,pd)
                                  full_ks(pb,pc) = full_ks(pb,pc)-tmp*full_density(pa,pd)
                                  full_ks(pb,pd) = full_ks(pb,pd)-tmp*full_density(pa,pc)
                                  IF( nspins == 2 ) THEN
                                    full_ks_beta(pa,pd) = full_ks_beta(pa,pd)-tmp*full_density_beta(pb,pc)
                                    full_ks_beta(pa,pc) = full_ks_beta(pa,pc)-tmp*full_density_beta(pb,pd)
                                    full_ks_beta(pb,pc) = full_ks_beta(pb,pc)-tmp*full_density_beta(pa,pd)
                                    full_ks_beta(pb,pd) = full_ks_beta(pb,pd)-tmp*full_density_beta(pa,pc)
                                  END IF
                                  IF(geometry_did_change .AND. (.NOT.buffer_overflow)) THEN
                                    qs_env%x_data%x_buffer(integral_counter) = &
                                    primitive_integrals(p_index)
                                    integral_counter = integral_counter + 1
                                  END IF
                                END DO
                              END DO
                            END DO
                          END DO
                        END IF
                      END DO !lset
                    END DO !kset
                  END DO !jset
                END DO !iset

              END IF !symmetry cycle 
            END DO !latom
            first_l_loop = .FALSE.
          END DO !katom
          first_k_loop = .FALSE.
        END DO !jatom
        first_j_loop = .FALSE.
      END DO !iatom
    END DO outer !bin


    DEALLOCATE(primitive_integrals,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(T1,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    ene_x_aa = 0.0_dp
    ene_x_bb = 0.0_dp
   
    ! symmetrize ks_matrix, calculate energies 
    DO pa=1,SIZE(full_ks,1)
      DO pd=pa,SIZE(full_ks,2)
        IF(pa/=pd) THEN
          full_ks(pa,pd) = (full_ks(pa,pd) + full_ks(pd,pa))*0.5_dp
          full_ks(pd,pa) = full_ks(pa,pd)
          ene_x_aa = ene_x_aa + 2.0_dp * full_density(pa,pd)*full_ks(pa,pd)
          IF(nspins==2) THEN
            full_ks_beta(pa,pd) = (full_ks_beta(pa,pd) + &
                                   full_ks_beta(pd,pa)) * 0.5_dp
            full_ks_beta(pd,pa) = full_ks_beta(pa,pd) 
            ene_x_bb = ene_x_bb + 2.0_dp * full_density_beta(pa,pd)*full_ks_beta(pa,pd)
          END IF
        ELSE
          ene_x_aa = ene_x_aa + full_density(pa,pd)*full_ks(pa,pd)
          IF(nspins==2) ene_x_bb = ene_x_bb + full_density_beta(pa,pd)*full_ks_beta(pa,pd)
        END IF
      END DO
    END DO
    
    CALL mp_sum(full_ks,para_env%group)
    IF(nspins==2) CALL mp_sum(full_ks_beta,para_env%group)
    CALL mp_sum(ene_x_aa, para_env%group)
    IF(nspins==2) CALL mp_sum(ene_x_bb,para_env%group)

    CALL timeset(routineN//"_dist_KS",handle_dist_ks)
    DO iatom=1,ks_matrix(1)%matrix%nblock_row
      block_node => first_block_node(ks_matrix(1)%matrix,iatom)
      IF(nspins == 2) block_node_beta => first_block_node(ks_matrix(2)%matrix,iatom)
      ! parallel distribution might mean that there won't be an associated latom
      IF (.NOT. ASSOCIATED(block_node)) CYCLE
      atom_kind => particle_set(iatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_a)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_a,&
                             nsgf=nsgfa_total,nset=nseta,&
                             first_sgf=first_sgfa,&
                             nsgf_set=nsgfa)
      pa1 = last_sgf_global(iatom) - nsgfa_total
      block_node => first_block_node(ks_matrix(1)%matrix,iatom)
      IF(nspins == 2) block_node_beta => first_block_node(ks_matrix(2)%matrix,iatom)
      DO WHILE (ASSOCIATED(block_node))
        CALL get_block_node(block_node=block_node,&
                            block_col=latom,&
                            BLOCK=sparse_block)
        block_node => next_block_node(block_node)
        IF(nspins == 2) THEN
          CALL get_block_node(block_node=block_node_beta,&
                              block_col=latom,&
                              BLOCK=sparse_block_beta)
          block_node_beta => next_block_node(block_node_beta)
        END IF
        atom_kind => particle_set(latom)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atom_kind,&
                             orb_basis_set=orb_basis_d)
        CALL get_gto_basis_set(gto_basis_set=orb_basis_d,&
                               nsgf=nsgfd_total,nset=nsetd,&
                               first_sgf=first_sgfd,&
                               nsgf_set=nsgfd)
        pd1 = last_sgf_global(latom) - nsgfd_total
        DO lset=1,nsetd
          sgfd = first_sgfd(1,lset)
          DO iset=1,nseta
            sgfa = first_sgfa(1,iset)
            DO md = 1,nsgfd(lset)
              pd2 = md+sgfd-1
              pd = pd1 + pd2
              DO ma = 1,nsgfa(iset)
                pa2 = ma+sgfa-1
                pa = pa1 + pa2
                sparse_block(pa2,pd2) = full_ks(pa,pd)
                IF(nspins==2) sparse_block_beta(pa2,pd2) = full_ks_beta(pa,pd)
              END DO 
            END DO
          END DO    
        END DO
      END DO
    END DO
    
    CALL timestop(handle_dist_ks)
    energy%ex = 0.5_dp*(ene_x_aa+ene_x_bb)

    CALL timestop(handle_main)
    
    iw = cp_print_key_unit_nr(logger,hf_section,"HF_INFO",&
             extension=".scfLog",error=error)
    total_integrals_to_store_max=total_integrals_to_store
    total_integrals_to_store_aver=(total_integrals_to_store+para_env%num_pe-1)/(para_env%num_pe)
    CALL mp_max(total_integrals_to_store_max,para_env%group)
    CALL mp_sum(total_integrals_to_store_aver,para_env%group)

    IF (iw>0) THEN
      IF(memory_usage < integral_counter * dp_size /1024_dp/1024_dp) THEN
        memory_usage = integral_counter * dp_size /1024_dp/1024_dp
      END IF
      WRITE (UNIT=iw,FMT="(/,(T3,A,T60,F20.10))")&
               "Memory used for x_buffer [MB]:                 ",memory_usage
      WRITE (UNIT=iw,FMT="((T3,A,T60,I20))")&
               "Number of integrals stored in buffer:          ",integral_counter-1

      WRITE (UNIT=iw,FMT="((T3,A,T60,I20))")&
               "Maximum number of elements to store:          ",total_integrals_to_store_max

      WRITE (UNIT=iw,FMT="((T3,A,T60,I20),/)")&
               "Average number of elements to store:          ",total_integrals_to_store_aver
    END IF

    CALL cp_print_key_finished_output(iw,logger,hf_section,&
                                      "HF_INFO", error=error)

    IF( do_periodic ) THEN
      DEALLOCATE(screen_periodic, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    DEALLOCATE(screen_set, screen_atom, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(full_ks, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(last_sgf_global,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(full_density,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(is_assoc_atomic_block,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(nspins==2) THEN
      DEALLOCATE(full_density_beta,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(full_ks_beta, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    !
    ! this wraps to free_libint, but currently causes a segfault
    ! as a result, we don't call it, but some memory remains allocated
    !  
    ! CALL terminate_libint(lib)

    ! many timings somewhat more clear
    CALL mp_sync(para_env%group)

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE integrate_four_center_symm
!***************************************************************************

!!***** four_center_integrals/derivatives_four_center_symm [1.0] *
!!
!!   NAME
!!     derivatives_four_center_symm
!!
!!   FUNCTION
!!     computes four center derivatives for a full basis set and updates the
!!     forces%fock_4c arrays. Uses all 8 eri symmetries
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_env 
!!     - ks_matrix
!!     - energy
!!     - rho: density matrix
!!     - para_env
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     06.2007 created [Manuel Guidon]
!!     08.2007 optimized load balance [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE derivatives_four_center_symm(qs_env,ks_matrix,energy,rho,hf_section,para_env,& 
                                    do_periodic, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: hf_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(IN)                      :: do_periodic
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'derivatives_four_center_symm', &
      routineP = moduleN//':'//routineN

    INTEGER :: bin, box_idx1, box_idx2, box_idx3, box_idx4, coord, &
      forces_map(4,2), handle, i_atom, iatom, ikind, iset, j_atom, jatom, &
      jkind, jset, k_atom, katom, kkind, kset, l_atom, latom, lkind, lset, &
      ma, max_am, max_k, max_set, maxsgf, mb, mc, md, my_bin_size, my_istart, &
      my_jstart, my_kstart, my_lstart, natom, nco_max, ncoa, ncob, ncoc, &
      ncod, nseta, nsetb, nsetc, nsetd, nsgfa_total, nsgfb_total, &
      nsgfc_total, nsgfd_total, nspins, p_index, pa, pa1, pa2, pb, pb1, pb2, &
      pc, pc1, pc2, pd, pd1, pd2, pgf_radius_box, sgfa, sgfb, sgfc, sgfd, stat
    INTEGER(KIND=int_8)                      :: my_current_counter
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of, &
                                                last_sgf_global
    INTEGER, DIMENSION(:), POINTER :: la_max, la_min, lb_max, lb_min, lc_max, &
      lc_min, ld_max, ld_min, ncgfa, ncgfb, ncgfc, ncgfd, npgfa, npgfb, &
      npgfc, npgfd, nsgfa, nsgfb, nsgfc, nsgfd
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb, &
                                                first_sgfc, first_sgfd, &
                                                is_assoc_atomic_block
    LOGICAL                                  :: failure, first_j_loop, &
                                                first_k_loop, first_l_loop, &
                                                screen_pmat_forces
    REAL(dp)                                 :: Kbox, omega_box, Rbox
    REAL(KIND=dp) :: box_size(3), eps_schwarz, fac, hf_fraction, max_val1, &
      max_val2, pmax_max, ra(3), rb(3), rc(3), rd(3), screen_max, symm_fac, &
      temp1, temp1b, temp2, temp3, temp3b
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:), TARGET                   :: primitive_forces
    REAL(KIND=dp), DIMENSION(:), POINTER     :: T1, T2
    REAL(KIND=dp), DIMENSION(:, :), POINTER :: full_density, &
      full_density_beta, pmax_atom, screen_atom, sparse_block, sphi_a, &
      sphi_b, sphi_c, sphi_d, work_forces, zeta, zetb, zetc, zetd
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: pmax_set, screen_set
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), POINTER      :: screen_periodic
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_a, orb_basis_b, &
                                                orb_basis_c, orb_basis_d
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(real_block_node_type), POINTER      :: block_node
    TYPE(x_hf_screening_type)                :: screening_parameter

! many timings somewhat more clear

    CALL mp_sync(para_env%group)
    CALL timeset(routineN,"I",'',handle)
    
    CALL get_qs_env(qs_env=qs_env, force=force, error=error)
    DO ikind = 1,SIZE(force)
      force(ikind)%fock_4c(:,:) = 0.0_dp
    END DO

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    error=error)
 
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxsgf=maxsgf,&
                             maxnset=max_set,&
                             maxlgto=max_am)
   
    natom = SIZE(particle_set,1)
    !one atom systems have no contribution to forces
    IF( natom == 1) THEN
      CALL timestop(0.0_dp,handle)
      RETURN
    END IF


    CALL section_vals_val_get(hf_section,"FRACTION",r_val=hf_fraction,error=error)
    CALL section_vals_val_get(hf_section,"EPS_SCHWARZ_FORCES",r_val=eps_schwarz,error=error)
    CALL section_vals_val_get(hf_section,"SCREEN_PMAT_FORCES",l_val=screen_pmat_forces,error=error)


    CALL initialize_libderiv(deriv, max_am, error)

    screening_parameter = qs_env%x_data%screening_parameter

    ALLOCATE(work_forces(((max_am+1)*(max_am+2)/2)**4,12),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE (atom_of_kind(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (kind_of(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure) 
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             atom_of_kind=atom_of_kind,&
                             kind_of=kind_of)

    nspins = qs_env%dft_control%nspins
 
    ALLOCATE(last_sgf_global(0:natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    last_sgf_global(0)=0
    DO iatom=1,natom
      atom_kind => particle_set(iatom)%atomic_kind
 
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_a)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_a,&
                             nsgf=nsgfa_total)
      last_sgf_global(iatom) = last_sgf_global(iatom-1)+nsgfa_total
    END DO
 
    ALLOCATE(full_density(1,1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL copy_local_sm_to_replicated_fm(rho%rho_ao(1)%matrix,full_density)
    CALL symmetrize_matrix(full_density,"upper_to_lower")
    CALL mp_sum(full_density,para_env%group)
    IF(nspins == 2) THEN
      ALLOCATE(full_density_beta(1,1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      CALL copy_local_sm_to_replicated_fm(rho%rho_ao(2)%matrix,full_density_beta)
      CALL symmetrize_matrix(full_density_beta,"upper_to_lower")
      CALL mp_sum(full_density_beta,para_env%group)
    ELSE
      NULLIFY(full_density_beta)
    END IF

    CALL compute_pmax_matrices(pmax_set,pmax_atom,pmax_max, &
                  full_density,full_density_beta,screen_pmat_forces,particle_set,atomic_kind_set,error)

    ! 
    ! notice that we use eps_schwarz == EPS_SCHWARZ_FORCES
    !
    max_k = qs_env%x_data%number_of_k_loops
    IF( do_periodic ) THEN
      box_size(1) = qs_env%cell%hmat(1,1)
      box_size(2) = qs_env%cell%hmat(2,2)
      box_size(3) = qs_env%cell%hmat(3,3)
      max_k = MIN(max_k, SIZE(qs_env%x_data%neighbor_cells))
      ALLOCATE(screen_periodic(natom,natom,max_set,max_set,SIZE(qs_env%x_data%neighbor_cells)))
      CALL comp_screen_mat_per(screen_set,screen_atom,screen_max,screen_periodic,particle_set,atomic_kind_set,&
                                               screening_parameter,.FALSE.,&
                                               neighbor_cells=qs_env%x_data%neighbor_cells,&
                                               box_size=box_size,&
                                               max_k=max_k,error=error)
    ELSE
      CALL compute_screening_matrices(screen_set,screen_atom,screen_max,particle_set,atomic_kind_set,&
                                      screening_parameter,.FALSE.,error=error)
      !farfield screening matrix already there
    END IF

    ALLOCATE(is_assoc_atomic_block(natom,natom))
    is_assoc_atomic_block = 0
    DO iatom=1,ks_matrix(1)%matrix%nblock_row
      block_node => first_block_node(ks_matrix(1)%matrix,iatom)
      DO WHILE (ASSOCIATED(block_node))
        CALL get_block_node(block_node=block_node,&
                            block_col=jatom,&
                            BLOCK=sparse_block)

        is_assoc_atomic_block(iatom,jatom) = 1
        is_assoc_atomic_block(jatom,iatom) = 1
        block_node => next_block_node(block_node)
      END DO
    END DO

    CALL mp_sum(is_assoc_atomic_block,para_env%group)

    CALL load_balance_forces(qs_env%x_data,eps_schwarz,particle_set,atomic_kind_set,para_env,&
                             screening_parameter,screen_set,screen_atom,screen_max,&
                             pmax_set,is_assoc_atomic_block,ks_matrix,do_periodic,error)

    ! precompute maximum nco and allocate scratch
    nco_max=0
    DO iatom=1,natom
      atom_kind => particle_set(iatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_a)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_a, lmax=la_max, lmin=la_min,&
                             npgf=npgfa, nset=nseta, zet=zeta, nsgf_set=nsgfa,&
                             first_sgf=first_sgfa, sphi=sphi_a, nsgf=nsgfa_total)
      DO iset = 1, nseta
         nco_max = MAX(nco_max,npgfa(iset)*ncoset(la_max(iset)))
      ENDDO
    ENDDO
    ALLOCATE(primitive_forces(12*nco_max**4),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    primitive_forces=0.0_dp
    ALLOCATE(T1(nco_max**4),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    T1=0.0_dp

    my_bin_size = SIZE(qs_env%x_data%distribution_forces)
    first_j_loop = .TRUE.
    first_l_loop = .TRUE.
    first_k_loop = .TRUE.
    outer: DO bin=1,my_bin_size
      my_istart = qs_env%x_data%distribution_forces(bin)%istart
      my_jstart = qs_env%x_data%distribution_forces(bin)%jstart
      my_kstart = qs_env%x_data%distribution_forces(bin)%kstart
      my_lstart = qs_env%x_data%distribution_forces(bin)%lstart
      my_current_counter = 0
      IF(my_istart == 0 .OR. my_jstart == 0 .OR. my_kstart == 0 .OR. my_lstart == 0 .OR. &
         qs_env%x_data%distribution_forces(bin)%number_of_integrals == 0) CYCLE
      DO iatom=my_istart,natom
        atom_kind => particle_set(iatom)%atomic_kind

        CALL get_atomic_kind(atomic_kind=atom_kind,&
                             orb_basis_set=orb_basis_a)
        CALL get_gto_basis_set(gto_basis_set=orb_basis_a,&
                               lmax=la_max,&
                               lmin=la_min,&
                               npgf=npgfa,&
                               nset=nseta,&
                               zet=zeta,&
                               nsgf_set=nsgfa,&
                               ncgf_set=ncgfa,&
                               first_sgf=first_sgfa,&
                               sphi=sphi_a,&
                               nsgf=nsgfa_total)
        ra=particle_set(iatom)%r(:)
        i_atom = atom_of_kind(iatom)
        ikind = kind_of(iatom)
        pa1 = last_sgf_global(iatom) - nsgfa_total
        forces_map(1,1) = ikind
        forces_map(1,2) = i_atom
        IF(first_j_loop) jatom=my_jstart-1
        IF(.NOT. first_j_loop) jatom = iatom-1
        DO WHILE(jatom+1<=natom)
          jatom = jatom + 1

          IF (screen_atom(iatom,jatom) * screen_max < eps_schwarz) CYCLE
          atom_kind => particle_set(jatom)%atomic_kind
          CALL get_atomic_kind(atomic_kind=atom_kind,&
                               orb_basis_set=orb_basis_b)
          CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                                 lmax=lb_max,&
                                 lmin=lb_min,&
                                 npgf=npgfb,&
                                 nset=nsetb,&
                                 zet=zetb,&
                                 nsgf_set=nsgfb,&
                                 ncgf_set=ncgfb,&
                                 first_sgf=first_sgfb,&
                                 sphi=sphi_b,&
                                 nsgf=nsgfb_total)
          rb=particle_set(jatom)%r(:)
          j_atom = atom_of_kind(jatom)
          jkind = kind_of(jatom)
          pb1 = last_sgf_global(jatom) - nsgfb_total
          forces_map(2,1) = jkind
          forces_map(2,2) = j_atom
          IF(first_k_loop) katom=my_kstart-1
            IF(.NOT. first_k_loop) katom = 0
            DO WHILE(katom+1<=natom)
              katom = katom + 1
            atom_kind => particle_set(katom)%atomic_kind
            CALL get_atomic_kind(atomic_kind=atom_kind,&
                                 orb_basis_set=orb_basis_c)
            CALL get_gto_basis_set(gto_basis_set=orb_basis_c,&
                                   lmax=lc_max,&
                                   lmin=lc_min,&
                                   npgf=npgfc,&
                                   nset=nsetc,&
                                   zet=zetc,&
                                   nsgf_set=nsgfc,&
                                   ncgf_set=ncgfc,&
                                   first_sgf=first_sgfc,&
                                   sphi=sphi_c,&
                                   nsgf=nsgfc_total)
            rc=particle_set(katom)%r(:)
            k_atom = atom_of_kind(katom)
            kkind = kind_of(katom)
            pc1 = last_sgf_global(katom) - nsgfc_total
            forces_map(3,1) = kkind
            forces_map(3,2) = k_atom
            IF(first_l_loop) THEN
                  latom=my_lstart-1
               ELSE
                  latom = katom-1
               END IF
               DO WHILE(latom+1<=natom)
                latom = latom + 1
   
              !All four centers equivalent => zero-contribution
              IF((iatom==jatom .AND. iatom==katom .AND. iatom==latom)) CYCLE 
              IF (screen_atom(iatom,jatom) * screen_atom(katom,latom) < eps_schwarz) CYCLE
              IF(.NOT. (is_assoc_atomic_block(iatom,latom)>=1 .AND. &
                      is_assoc_atomic_block(iatom,katom)>=1 .AND. &
                      is_assoc_atomic_block(jatom,katom)>=1 .AND. &
                      is_assoc_atomic_block(jatom,latom)>=1 ) ) CYCLE
  
              IF(katom+latom<=iatom+jatom)  THEN
                IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE
  
                my_current_counter = my_current_counter + 1
                   IF(my_current_counter > qs_env%x_data%distribution_forces(bin)%number_of_integrals) THEN
                     first_j_loop = .TRUE.
                     first_l_loop = .TRUE.
                     first_k_loop = .TRUE.
                     CYCLE outer
                   END IF
  
                atom_kind => particle_set(latom)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atom_kind,&
                                     orb_basis_set=orb_basis_d)
                CALL get_gto_basis_set(gto_basis_set=orb_basis_d,&
                                       lmax=ld_max,&
                                       lmin=ld_min,&
                                       npgf=npgfd,&
                                       nset=nsetd,&
                                       zet=zetd,&
                                       nsgf_set=nsgfd,&
                                       ncgf_set=ncgfd,&
                                       first_sgf=first_sgfd,&
                                       sphi=sphi_d,&
                                       nsgf=nsgfd_total)
 
                rd=particle_set(latom)%r(:)
              
                l_atom = atom_of_kind(latom)
                lkind = kind_of(latom)
                pd1 = last_sgf_global(latom) - nsgfd_total
                forces_map(4,1) = lkind
                forces_map(4,2) = l_atom
               
                IF(nspins == 1) THEN
                  fac = 0.25_dp * hf_fraction
                ELSE
                  fac = 0.5_dp * hf_fraction
                END IF
                !calculate symmetry_factor
                symm_fac = 0.25_dp
                IF(iatom==jatom) symm_fac = symm_fac*2.0_dp
                IF(katom==latom) symm_fac = symm_fac*2.0_dp
                IF(iatom==katom .AND. jatom==latom .AND. iatom/=jatom .AND. katom/=latom) symm_fac = symm_fac*2.0_dp
                IF(iatom==katom .AND. iatom==jatom .AND. katom==latom) symm_fac = symm_fac*2.0_dp
  
                symm_fac = 1.0_dp/symm_fac
                fac = fac * symm_fac             
  
                DO iset = 1, nseta
                  ncoa = npgfa(iset)*ncoset(la_max(iset))
                  sgfa = first_sgfa(1,iset)
                  DO jset = 1,nsetb
                    ncob = npgfb(jset)*ncoset(lb_max(jset))
                    sgfb = first_sgfb(1,jset)
                    max_val1 = screen_set(iatom,jatom,iset,jset)
                    IF (max_val1 * screen_atom(katom,latom) < eps_schwarz) CYCLE
                    DO kset = 1,nsetc
                      ncoc = npgfc(kset)*ncoset(lc_max(kset))
                      sgfc = first_sgfc(1,kset)
                      DO lset = 1,nsetd
                        ncod = npgfd(lset)*ncoset(ld_max(lset)) 
                        sgfd = first_sgfd(1,lset)
   
                        max_val2 = max_val1 * screen_set(katom,latom,kset,lset) * &
                                   2.0_dp*MAX(pmax_set(iset,kset,iatom,katom) * pmax_set(jset,lset,jatom,latom),&
                                              pmax_set(iset,lset,iatom,latom) * pmax_set(kset,jset,katom,jatom))
                        IF(max_val2<eps_schwarz)  CYCLE
                        box_idx1 = (iatom-1)*natom+jatom
                        box_idx2 = (iset-1)*max_set+jset
                        box_idx3 = (katom-1)*natom+latom
                        box_idx4 = (kset-1)*max_set+lset
                        Rbox = (qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%center(1) - &
                                qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%center(1))**2 + &
                               (qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%center(2) - &
                                qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%center(2))**2 + &
                               (qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%center(3) - &
                                qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%center(3))**2
                        Rbox = SQRT(Rbox)
                        Rbox = Rbox - qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%R - &
                                      qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%R
                        pgf_radius_box = INT(1.0_dp/ &
                              (SQRT(2.0_dp*qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%alpha)*erfc(epsilon)))
                        pgf_radius_box = pgf_radius_box + INT(1.0_dp/ &
                              (SQRT(2.0_dp*qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%alpha)*erfc(epsilon))) + 2
                        IF(Rbox>pgf_radius_box) THEN
                          max_val2 = 0.0_dp
                          SELECT CASE(screening_parameter%screening_type)
                            CASE(do_hf_screening_none)
                              Kbox = qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                     qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%K
                              max_val2 = 1.0_dp/Rbox*Kbox* &
                                         2.0_dp*MAX(pmax_set(iset,kset,iatom,katom) * pmax_set(jset,lset,jatom,latom),&
                                         pmax_set(iset,lset,iatom,latom) * pmax_set(kset,jset,katom,jatom))
                           CASE(do_hf_screening_short)
                              omega_box = 1.0_dp/qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%alpha + &
                                          1.0_dp/qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%alpha + &
                                          1.0_dp/(screening_parameter%omega**2)
                              omega_box = 1.0_dp/omega_box
                              Kbox = qs_env%x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                     qs_env%x_data%far_field_boxes(box_idx3, box_idx4)%K
                              max_val2 = erfc(omega_box*Rbox)/Rbox*Kbox * &
                                         2.0_dp*MAX(pmax_set(iset,kset,iatom,katom) * pmax_set(jset,lset,jatom,latom),&
                                         pmax_set(iset,lset,iatom,latom) * pmax_set(kset,jset,katom,jatom))
                          END SELECT
                          IF(max_val2<eps_schwarz) CYCLE
                        END IF
                        IF( .NOT. do_periodic ) THEN
                          CALL forces4(ra, rb, rc, rd,&
                                       zeta(:,iset), zetb(:,jset), zetc(:,kset), zetd(:,lset),&
                                       la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                                       lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                                       npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset),&
                                       ncoa, ncob, ncoc, ncod,&
                                       primitive_forces,&
                                       work_forces,&
                                       eps_schwarz,&
                                       screening_parameter)
                        ELSE
                          CALL forces4_periodic(ra, rb, rc, rd,&
                                                zeta(:,iset), zetb(:,jset), zetc(:,kset), zetd(:,lset),&
                                                la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                                                lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                                                npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset),&
                                                ncoa, ncob, ncoc, ncod,&
                                                primitive_forces,&
                                                work_forces,&
                                                screening_parameter, box_size, &
                                                qs_env%x_data%neighbor_cells, screen_periodic(iatom,jatom,iset,jset,:),&
                                                screen_periodic(katom,latom,kset,lset,:), eps_schwarz, max_k)
                        END IF
 
                        DO coord = 1,12
                          T2 => primitive_forces((coord-1)*ncoa*ncob*ncoc*ncod+1:coord*ncoa*ncob*ncoc*ncod)
                          CALL dgemm("T","N",ncob*ncoc*ncod,nsgfa(iset),ncoa,&
                                     1.0_dp, T2(1), ncoa,&
                                     sphi_a(1,sgfa), SIZE(sphi_a,1),&
                                     0.0_dp, T1(1),ncob*ncoc*ncod)

                          CALL dgemm("T","N",nsgfa(iset)*ncoc*ncod,nsgfb(jset),ncob,&
                                     1.0_dp, T1(1),ncob,&
                                     sphi_b(1,sgfb), SIZE(sphi_b,1),&
                                     0.0_dp, T2(1),&
                                     nsgfa(iset)*ncoc*ncod)
  
                          CALL dgemm("T","N",nsgfa(iset)*nsgfb(jset)*ncod,nsgfc(kset),ncoc,&
                                     1.0_dp, T2(1),&
                                     ncoc,&
                                     sphi_c(1,sgfc), SIZE(sphi_c,1),&
                                     0.0_dp, T1(1),nsgfa(iset)*nsgfb(jset)*ncod)
  
                          CALL dgemm("T","N",nsgfa(iset)*nsgfb(jset)*nsgfc(kset),nsgfd(lset),ncod,&
                                     1.0_dp, T1(1),ncod,&
                                     sphi_d(1,sgfd), SIZE(sphi_d,1),&
                                     0.0_dp, T2(1),&
                                     nsgfa(iset)*nsgfb(jset)*nsgfc(kset))
  
  
                          p_index = 0
                          DO md = 1,nsgfd(lset)
                            pd2 = md+sgfd-1
                            pd = pd1 + pd2
                            DO mc = 1, nsgfc(kset)
                              pc2 = mc+sgfc-1
                              pc = pc1 + pc2
                              DO mb = 1, nsgfb(jset)
                                pb2 = mb+sgfb-1
                                pb = pb1 + pb2
                                temp1 = full_density(pc,pb) * fac
                                IF( nspins == 2 ) temp1b = full_density_beta(pc,pb) * fac
                                DO ma = 1,nsgfa(iset) 
                                  pa2 = ma+sgfa-1
                                  pa = pa1 + pa2
                                  temp3 = full_density(pb,pd) * fac
                                  IF( nspins == 2 ) temp3b = full_density_beta(pb,pd) * fac
                                  p_index = p_index + 1 
                                  temp2 = temp1 * full_density(pa,pd) + temp3 * full_density(pa,pc)
                                  force(forces_map((coord-1)/3 + 1,1))%fock_4c(MOD(coord-1,3) + 1,forces_map((coord-1)/3 + 1,2)) = &
                                  force(forces_map((coord-1)/3 + 1,1))%fock_4c(MOD(coord-1,3) + 1,forces_map((coord-1)/3 + 1,2)) - &
                                  temp2 * T2(p_index) !primitive_forces((coord-1)*ncoa*ncob*ncoc*ncod+p_index)
                                  IF(nspins == 2) THEN
                                    temp2 = temp1b*full_density_beta(pa,pd) + temp3b*full_density_beta(pa,pc)
                                    force(forces_map((coord-1)/3 + 1,1))%fock_4c(MOD(coord-1,3) + 1,&
                                          forces_map((coord-1)/3 + 1,2)) = &
                                    force(forces_map((coord-1)/3 + 1,1))%fock_4c(MOD(coord-1,3) + 1,&
                                          forces_map((coord-1)/3 + 1,2)) - &
                                    temp2 * T2(p_index)!primitive_forces((coord-1)*ncoa*ncob*ncoc*ncod+p_index)
                                  END IF 
                                END DO
                              END DO !ma
                            END DO !mb
                          END DO !mc
                        END DO !md
    
                      END DO !lset
                    END DO !ksetm
                  END DO !jset
                END DO !iset
              END IF !cycle
            END DO !latom
            first_l_loop = .FALSE.
          END DO !katom
          first_k_loop = .FALSE.
        END DO !jatom
        first_j_loop = .FALSE.
      END DO !iatom
    END DO outer !bin
    
    IF( do_periodic ) THEN
      DEALLOCATE(screen_periodic, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE(T1,primitive_forces, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(last_sgf_global,screen_set, screen_atom, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(full_density,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(is_assoc_atomic_block,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(nspins==2) THEN
      DEALLOCATE(full_density_beta,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE(work_forces, atom_of_kind,kind_of,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(pmax_set,pmax_atom,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    !
    ! this wraps to free_libint, but currently causes a segfault
    ! as a result, we don't call it, but some memory remains allocated
    !  
    ! CALL terminate_libderiv(deriv)
    
    CALL mp_sync(para_env%group)

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE derivatives_four_center_symm

!***************************************************************************

!!***** four_center_integrals/load_balance_energy [1.0] *
!!
!!   NAME
!!     load_balance_energy 
!!
!!   FUNCTION
!!     Distributes the computation of eri's to all available processes.
!!
!!   NOTES
!!     The optimization is done via a simple binning procedure:
!!     In a first step the total amount of integrals in the system is calculated,
!!     taking into account the sparsity of the KS-matrix and the screening based
!!     on the Schwarz-Inequality. Then, the procedure uses different bin-sizes 
!!     in order to optimize the load balance.
!!     At the end each process owns an unique array of inidices-ranges 
!!     that are used to decide whether a process has to calculate a certain
!!     bunch of integrals.
!!     
!!
!!   INPUTS
!!     - x_data: Object that stores the indices array
!!     - eps_schwarz: screening parameter
!!     - particle_set, atomic_kind_set, para_env
!!     - screening_parameter: short/longrange screening information 
!!     - screen_set, screen_atom, screen_max: screening matrices
!!     - is_assoc_atomic_kind: KS-matrix sparsity
!!     - ks_matrix: used for parallel computation of total number of integrals
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     06.2007 created [Manuel Guidon]
!!     08.2007 new parallel scheme [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE load_balance_energy(x_data, eps_schwarz, particle_set,atomic_kind_set,para_env,&
                                 screening_parameter, screen_set, screen_atom, screen_max,&
                                 is_assoc_atomic_block_global,ks_matrix,do_periodic,error, pmax_set)
    TYPE(x_hf_type), POINTER                 :: x_data
    REAL(KIND=dp), INTENT(IN)                :: eps_schwarz
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(x_hf_screening_type), INTENT(IN)    :: screening_parameter
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: screen_set
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: screen_atom
    REAL(KIND=dp), INTENT(IN)                :: screen_max
    INTEGER, DIMENSION(:, :), POINTER        :: is_assoc_atomic_block_global
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    LOGICAL                                  :: do_periodic
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(dp), DIMENSION(:, :, :, :), &
      OPTIONAL, POINTER                      :: pmax_set

    CHARACTER(LEN=*), PARAMETER :: routineN = 'load_balance_energy', &
      routineP = moduleN//':'//routineN

    INTEGER :: box_idx1, box_idx2, box_idx3, box_idx4, first_from, first_to, &
      handle, i, iatom, icpu, ipos, iset, jatom, jset, katom, kset, latom, &
      lset, max_set, mepos, natom, nseta, nsetb, nsetc, nsetd, &
      pgf_radius_box, restore_iatom, restore_jatom, restore_katom, &
      restore_latom, stat
    INTEGER(KIND=int_8) :: current_counter, distribution_counter_end, &
      distribution_counter_start, int_per_cpu, max_dev, min_dev, nbins, ncpu, &
      opt_max, opt_nbin, restore_counter, total_integrals, total_integrals2
    INTEGER(KIND=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: cost_per_cpu
    INTEGER(KIND=int_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: cost_array
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: tmp_pos
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: first_dist
    INTEGER, DIMENSION(:), POINTER           :: npgfa, npgfb, npgfc, npgfd, &
                                                nsgfa, nsgfb, nsgfc, nsgfd
    LOGICAL                                  :: do_p_screening, failure
    REAL(dp)                                 :: K, Kbox, omega_box, ra(3), &
                                                rb(3), Rbox, rc(3), rd(3)
    REAL(KIND=dp)                            :: max_val1, max_val2
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: zeta, zetb, zetc, zetd
    TYPE(gto_basis_set_type), POINTER        :: orb_basis
    TYPE(x_hf_distribution), DIMENSION(:), &
      POINTER                                :: ptr_to_tmp_dist
    TYPE(x_hf_distribution), &
      DIMENSION(:, :), POINTER               :: tmp_dist

    do_p_screening = .FALSE.
    IF(PRESENT(pmax_set)) do_p_screening = .TRUE.
    

    natom=SIZE(particle_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, maxnset=max_set) 
    CALL timeset(routineN,handle)
    failure=.FALSE.

    ncpu = para_env%num_pe
    IF(ncpu == 1) THEN
      ALLOCATE(tmp_dist(1,1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      tmp_dist(1,1)%number_of_integrals = HUGE(tmp_dist(1,1)%number_of_integrals)
      tmp_dist(1,1)%istart = 1
      tmp_dist(1,1)%jstart = 1
      tmp_dist(1,1)%kstart = 1
      tmp_dist(1,1)%lstart = 1
      ptr_to_tmp_dist => tmp_dist(1,:)
      CALL x_hf_set_distr_energy(ptr_to_tmp_dist,x_data,error=error)
      DEALLOCATE(tmp_dist,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
      !Simple distribution on first loop variable
      ALLOCATE(first_dist(ncpu,2), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)      
      first_dist = -1
      first_from = 1
      first_to = 1
      outer: DO icpu=1,ncpu
        IF(first_from > natom ) EXIT outer
        first_to = first_from + natom/ncpu + 1
        first_to = MIN(first_to, natom)
        first_dist(icpu,1) = first_from
        first_dist(icpu,2) = first_to
        IF(first_to == natom) EXIT outer
        first_from = first_to + 1
       END DO outer

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set)
      !Calculate total numbers of integrals that have to be calculated (wrt screening and symmetry)
      total_integrals = 0
      IF(first_dist(para_env%mepos+1,1) /= -1 ) THEN
        DO iatom=first_dist(para_env%mepos+1,1),first_dist(para_env%mepos+1,2)!1,natom
          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, orb_basis_set=orb_basis)
          ra=particle_set(iatom)%r(:)
          CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nseta, nsgf_set=nsgfa, zet=zeta, npgf=npgfa)
          DO jatom=iatom,natom
            IF (screen_atom(iatom,jatom) * screen_max < eps_schwarz) CYCLE
            CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind, orb_basis_set=orb_basis)
            rb=particle_set(jatom)%r(:)
            CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetb, nsgf_set=nsgfb, zet=zetb, npgf=npgfb)
            DO katom=1,natom
              CALL get_atomic_kind(atomic_kind=particle_set(katom)%atomic_kind, orb_basis_set=orb_basis)
              rc=particle_set(katom)%r(:)
              CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetc, nsgf_set=nsgfc, zet=zetc, npgf=npgfc)
              DO latom=katom,natom
                IF (screen_atom(iatom,jatom) * screen_atom(katom,latom) < eps_schwarz) CYCLE
                IF(.NOT. (is_assoc_atomic_block_global(iatom,latom)  >=1 .AND. &
                          is_assoc_atomic_block_global(iatom,katom)  >=1 .AND. &
                          is_assoc_atomic_block_global(jatom,katom)  >=1 .AND. &
                          is_assoc_atomic_block_global(jatom,latom)  >=1 ) ) CYCLE 
 
                IF((katom+latom<=iatom+jatom))  THEN
                  IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE
                 
                  CALL get_atomic_kind(atomic_kind=particle_set(latom)%atomic_kind, orb_basis_set=orb_basis)
                  rd=particle_set(latom)%r(:)
                  CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetd, nsgf_set=nsgfd, zet=zetd, npgf=npgfd)

                  DO iset = 1, nseta
                    DO jset = 1,nsetb
                      max_val1 = screen_set(iatom,jatom,iset,jset)
                      IF (max_val1 * screen_atom(katom,latom) < eps_schwarz) CYCLE
                      DO kset = 1,nsetc
                        DO lset = 1,nsetd
                          max_val2 = max_val1 * screen_set(katom,latom,kset,lset)
                          IF( do_p_screening) THEN
                            max_val2 = max_val2 * MAX(pmax_set(iset,kset,iatom,katom), pmax_set(jset,lset,jatom,latom),&
                                                  pmax_set(iset,lset,iatom,latom) , pmax_set(kset,jset,katom,jatom))
                          END IF
                          IF(max_val2<eps_schwarz) CYCLE
                          box_idx1 = (iatom-1)*natom+jatom
                          box_idx2 = (iset-1)*max_set+jset
                          box_idx3 = (katom-1)*natom+latom
                          box_idx4 = (kset-1)*max_set+lset
                          Rbox = (x_data%far_field_boxes(box_idx1,box_idx2)%center(1) - &
                                  x_data%far_field_boxes(box_idx3,box_idx4)%center(1))**2 + &
                                 (x_data%far_field_boxes(box_idx1,box_idx2)%center(2) - &
                                  x_data%far_field_boxes(box_idx3,box_idx4)%center(2))**2 + &
                                 (x_data%far_field_boxes(box_idx1,box_idx2)%center(3) - &
                                  x_data%far_field_boxes(box_idx3,box_idx4)%center(3))**2
                          Rbox = SQRT(Rbox)
                          Rbox = Rbox - x_data%far_field_boxes(box_idx1, box_idx2)%R - &
                                        x_data%far_field_boxes(box_idx3, box_idx4)%R
                          pgf_radius_box = INT(1.0_dp/ &
                                (SQRT(2.0_dp*x_data%far_field_boxes(box_idx1, box_idx2)%alpha)*erfc(epsilon)))
                          pgf_radius_box = pgf_radius_box + INT(1.0_dp/ &
                                (SQRT(2.0_dp*x_data%far_field_boxes(box_idx3, box_idx4)%alpha)*erfc(epsilon))) + 2
                          IF(Rbox>pgf_radius_box) THEN
                            max_val2 = 0.0_dp
                            SELECT CASE(screening_parameter%screening_type)
                              CASE(do_hf_screening_none)
                                Kbox = x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                       x_data%far_field_boxes(box_idx3, box_idx4)%K
                                max_val2 = 1.0_dp/Rbox*Kbox
                                IF( do_p_screening) THEN
                                  max_val2 = max_val2 * MAX(pmax_set(iset,kset,iatom,katom), pmax_set(jset,lset,jatom,latom),&
                                             pmax_set(iset,lset,iatom,latom) , pmax_set(kset,jset,katom,jatom))
                                END IF
                              CASE(do_hf_screening_short)
                                omega_box = 1.0_dp/x_data%far_field_boxes(box_idx1, box_idx2)%alpha + &
                                            1.0_dp/x_data%far_field_boxes(box_idx3, box_idx4)%alpha + &
                                            1.0_dp/(screening_parameter%omega**2)
                                omega_box = 1.0_dp/omega_box
                                Kbox = x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                       x_data%far_field_boxes(box_idx3, box_idx4)%K
                                max_val2 = erfc(omega_box*Rbox)/Rbox*Kbox
                                IF( do_p_screening) THEN
                                  max_val2 = max_val2 * MAX(pmax_set(iset,kset,iatom,katom), pmax_set(jset,lset,jatom,latom),&
                                                 pmax_set(iset,lset,iatom,latom) , pmax_set(kset,jset,katom,jatom))
                                END IF
                            END SELECT
                            IF(max_val2<eps_schwarz) CYCLE
                          END IF
                          total_integrals = total_integrals + nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                END IF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      END IF
      CALL mp_sum_i8(total_integrals, para_env%group)
      opt_max = HUGE(opt_max)  
      ALLOCATE(cost_array(11,11))
      DO ipos = 1,11
        cost_array(ipos,1) = 2**(ipos-1)
        cost_array(ipos,2) = 0 
      END DO
     
      DO i = 1,11
        nbins = cost_array(i,1)
        IF(MOD(i,para_env%num_pe) /= para_env%mepos ) CYCLE
        !floor value
        total_integrals2 = total_integrals +(ncpu*nbins- MOD(total_integrals,ncpu*nbins))
        int_per_cpu = total_integrals2 / ncpu
        int_per_cpu = int_per_cpu/nbins
        ALLOCATE(cost_per_cpu(ncpu))
        cost_per_cpu = 0
        current_counter = 0
        restore_counter = 0
        ipos = 1
        mepos = 1
        distribution_counter_start = 1
        distribution_counter_end = 0
        DO iatom=1,natom
          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, orb_basis_set=orb_basis)
          ra=particle_set(iatom)%r(:)
          CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nseta, nsgf_set=nsgfa, zet=zeta, npgf=npgfa)
          DO jatom=iatom,natom
            IF (screen_atom(iatom,jatom) * screen_max < eps_schwarz) CYCLE
            CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind, orb_basis_set=orb_basis)
            rb=particle_set(jatom)%r(:)
            CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetb, nsgf_set=nsgfb, zet=zetb, npgf=npgfb)
            DO katom=1,natom
              CALL get_atomic_kind(atomic_kind=particle_set(katom)%atomic_kind, orb_basis_set=orb_basis)
              rc=particle_set(katom)%r(:)
              CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetc, nsgf_set=nsgfc, zet=zetc, npgf=npgfc)
              DO latom=katom,natom
                IF (screen_atom(iatom,jatom) * screen_atom(katom,latom) < eps_schwarz) CYCLE
                IF(.NOT. (is_assoc_atomic_block_global(iatom,latom)  >=1 .AND. &
                          is_assoc_atomic_block_global(iatom,katom)  >=1 .AND. &
                          is_assoc_atomic_block_global(jatom,katom)  >=1 .AND. &
                          is_assoc_atomic_block_global(jatom,latom)  >=1 ) ) CYCLE
  
                IF((katom+latom<=iatom+jatom))  THEN
                  IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE

                  distribution_counter_end = distribution_counter_end + 1
  
                  CALL get_atomic_kind(atomic_kind=particle_set(latom)%atomic_kind, orb_basis_set=orb_basis)
                  rd=particle_set(latom)%r(:)
                  CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetd, nsgf_set=nsgfd, zet=zetd, npgf=npgfd)
  
                  DO iset = 1, nseta
                    DO jset = 1,nsetb
                      max_val1 = screen_set(iatom,jatom,iset,jset)
                      IF (max_val1 * screen_atom(katom,latom) < eps_schwarz) CYCLE
                      DO kset = 1,nsetc
                        DO lset = 1,nsetd
                          max_val2 = max_val1 * screen_set(katom,latom,kset,lset)
                          IF( do_p_screening) THEN
                            max_val2 = max_val2 * MAX(pmax_set(iset,kset,iatom,katom), pmax_set(jset,lset,jatom,latom),&
                                       pmax_set(iset,lset,iatom,latom) , pmax_set(kset,jset,katom,jatom))
                          END IF

                          IF(max_val2<eps_schwarz) CYCLE
                          box_idx1 = (iatom-1)*natom+jatom
                          box_idx2 = (iset-1)*max_set+jset
                          box_idx3 = (katom-1)*natom+latom
                          box_idx4 = (kset-1)*max_set+lset
                          Rbox = (x_data%far_field_boxes(box_idx1, box_idx2)%center(1) - &
                                  x_data%far_field_boxes(box_idx3, box_idx4)%center(1))**2 + &
                                 (x_data%far_field_boxes(box_idx1, box_idx2)%center(2) - &
                                  x_data%far_field_boxes(box_idx3, box_idx4)%center(2))**2 + &
                                 (x_data%far_field_boxes(box_idx1, box_idx2)%center(3) - &
                                  x_data%far_field_boxes(box_idx3, box_idx4)%center(3))**2
                          Rbox = SQRT(Rbox)
                          Rbox = Rbox - x_data%far_field_boxes(box_idx1, box_idx2)%R - &
                                        x_data%far_field_boxes(box_idx3, box_idx4)%R
                          pgf_radius_box = INT(1.0_dp/ &
                                 (SQRT(2.0_dp*x_data%far_field_boxes(box_idx1, box_idx2)%alpha)*erfc(epsilon)))
                          pgf_radius_box = pgf_radius_box + INT(1.0_dp/ &
                                 (SQRT(2.0_dp*x_data%far_field_boxes(box_idx3, box_idx4)%alpha)*erfc(epsilon))) + 2
   
                          IF(Rbox>pgf_radius_box) THEN
                            max_val2 = 0.0_dp
                            SELECT CASE(screening_parameter%screening_type)
                              CASE(do_hf_screening_none)
                                Kbox = x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                       x_data%far_field_boxes(box_idx3, box_idx4)%K
                                max_val2 = 1.0_dp/Rbox*Kbox
                                IF( do_p_screening) THEN
                                  max_val2 = max_val2 * MAX(pmax_set(iset,kset,iatom,katom), pmax_set(jset,lset,jatom,latom),&
                                             pmax_set(iset,lset,iatom,latom) , pmax_set(kset,jset,katom,jatom))
                                END IF
                              CASE(do_hf_screening_short)
                                omega_box = 1.0_dp/x_data%far_field_boxes(box_idx1, box_idx2)%alpha + &
                                            1.0_dp/x_data%far_field_boxes(box_idx3, box_idx4)%alpha + &
                                            1.0_dp/(screening_parameter%omega**2)
                                omega_box = 1.0_dp/omega_box
                                Kbox = x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                       x_data%far_field_boxes(box_idx3, box_idx4)%K
                                max_val2 = erfc(omega_box*Rbox)/Rbox*Kbox
                                IF( do_p_screening) THEN
                                  max_val2 = max_val2 * MAX(pmax_set(iset,kset,iatom,katom), pmax_set(jset,lset,jatom,latom),&
                                             pmax_set(iset,lset,iatom,latom) , pmax_set(kset,jset,katom,jatom))
                                END IF
                            END SELECT
                            IF(max_val2<eps_schwarz) CYCLE
                          END IF
                          current_counter = current_counter + nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                  IF( current_counter >= int_per_cpu) THEN
                    cost_per_cpu(mepos) = cost_per_cpu(mepos) + current_counter
                    restore_counter = current_counter
                    current_counter = 0
                    mepos = mepos + 1
                    IF( mepos>ncpu) THEN
                      mepos = mepos - ncpu
                      ipos = ipos + 2
                    END IF
                  END IF
                END IF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
   
        !fill the last proc
        IF(ncpu > 0 ) THEN
          IF(current_counter == 0 ) THEN
            cost_per_cpu(ncpu) = cost_per_cpu(ncpu) + restore_counter
          ELSE
            cost_per_cpu(ncpu) = cost_per_cpu(ncpu) + current_counter
          END IF
        END IF
       
        max_dev = 0
        min_dev = 0
        DO mepos = 1,ncpu 
          IF(total_integrals/ncpu>cost_per_cpu(mepos)) THEN
            IF(ABS(total_integrals/ncpu-cost_per_cpu(mepos))>min_dev) min_dev = ABS(total_integrals/ncpu-cost_per_cpu(mepos))
          ELSE
            IF(ABS(total_integrals/ncpu-cost_per_cpu(mepos))>max_dev) max_dev = ABS(total_integrals/ncpu-cost_per_cpu(mepos))
          END IF
        END DO
        cost_array(i,2) = max_dev 
  
        DEALLOCATE(cost_per_cpu,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      END DO

      CALL mp_sum_i8(cost_array(:,2), para_env%group)
      opt_max = HUGE(opt_max)
      DO i = 1,11
        IF(1000*cost_array(i,2)/(total_integrals/ncpu)<opt_max) THEN
          opt_max = 1000*cost_array(i,2)/(total_integrals/ncpu)
          opt_nbin = cost_array(i,1)
        END IF
      END DO
      !Call routine again for nbins = opt_nbin
      nbins = opt_nbin
      !floor value
      int_per_cpu = CEILING(REAL(total_integrals)/REAL(ncpu)/REAL(nbins),int_8)
      ALLOCATE(tmp_dist(ncpu,nbins),STAT=stat)
      tmp_dist(:,:)%istart = 0
      tmp_dist(:,:)%jstart = 0
      tmp_dist(:,:)%kstart = 0
      tmp_dist(:,:)%lstart = 0
      tmp_dist(:,:)%number_of_integrals = 0
      tmp_dist(1,1)%istart = 1
      tmp_dist(1,1)%jstart = 1
      tmp_dist(1,1)%kstart = 1
      tmp_dist(1,1)%lstart = 1
      ALLOCATE(tmp_pos(ncpu))
      tmp_pos = 1

      current_counter = 0
      restore_counter = 0
      ipos = 1
      mepos = 1
      distribution_counter_start = 1
      distribution_counter_end = 0

      DO iatom=1,natom
        CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, orb_basis_set=orb_basis)
        ra=particle_set(iatom)%r(:)
        CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nseta, nsgf_set=nsgfa, zet=zeta, npgf=npgfa)
        DO jatom=iatom,natom
          IF (screen_atom(iatom,jatom) * screen_max < eps_schwarz) CYCLE
          CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind, orb_basis_set=orb_basis)
          rb=particle_set(jatom)%r(:)
          CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetb, nsgf_set=nsgfb, zet=zetb, npgf=npgfb)
          DO katom=1,natom
            CALL get_atomic_kind(atomic_kind=particle_set(katom)%atomic_kind, orb_basis_set=orb_basis)
            rc=particle_set(katom)%r(:)
            CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetc, nsgf_set=nsgfc, zet=zetc, npgf=npgfc)
            DO latom=katom,natom
              IF (screen_atom(iatom,jatom) * screen_atom(katom,latom) < eps_schwarz) CYCLE
              IF(.NOT. (is_assoc_atomic_block_global(iatom,latom)  >=1 .AND. &
                        is_assoc_atomic_block_global(iatom,katom)  >=1 .AND. &
                        is_assoc_atomic_block_global(jatom,katom)  >=1 .AND. &
                        is_assoc_atomic_block_global(jatom,latom)  >=1 ) ) CYCLE
 
              IF((katom+latom<=iatom+jatom))  THEN
                IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE
                distribution_counter_end = distribution_counter_end + 1
                CALL get_atomic_kind(atomic_kind=particle_set(latom)%atomic_kind, orb_basis_set=orb_basis)
                rd=particle_set(latom)%r(:)
                CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetd, nsgf_set=nsgfd, zet=zetd, npgf=npgfd)
  
                DO iset = 1, nseta
                  DO jset = 1,nsetb
                    max_val1 = screen_set(iatom,jatom,iset,jset)
                    IF (max_val1 * screen_atom(katom,latom) < eps_schwarz) CYCLE
                    DO kset = 1,nsetc
                      DO lset = 1,nsetd
                        max_val2 = max_val1 * screen_set(katom,latom,kset,lset)
                        IF( do_p_screening) THEN
                          max_val2 = max_val2 * MAX(pmax_set(iset,kset,iatom,katom), pmax_set(jset,lset,jatom,latom),&
                                     pmax_set(iset,lset,iatom,latom) , pmax_set(kset,jset,katom,jatom))
                        END IF
                        IF(max_val2<eps_schwarz) CYCLE
                        box_idx1 = (iatom-1)*natom+jatom
                        box_idx2 = (iset-1)*max_set+jset
                        box_idx3 = (katom-1)*natom+latom
                        box_idx4 = (kset-1)*max_set+lset
                        Rbox = (x_data%far_field_boxes(box_idx1, box_idx2)%center(1) - &
                                x_data%far_field_boxes(box_idx3, box_idx4)%center(1))**2 + &
                               (x_data%far_field_boxes(box_idx1, box_idx2)%center(2) - &
                                x_data%far_field_boxes(box_idx3, box_idx4)%center(2))**2 + &
                               (x_data%far_field_boxes(box_idx1, box_idx2)%center(3) - &
                                x_data%far_field_boxes(box_idx3, box_idx4)%center(3))**2
                        Rbox = SQRT(Rbox)
                        Rbox = Rbox - x_data%far_field_boxes(box_idx1, box_idx2)%R - &
                                      x_data%far_field_boxes(box_idx3, box_idx4)%R
                        pgf_radius_box = INT(1.0_dp/&
                                (SQRT(2.0_dp*x_data%far_field_boxes(box_idx1, box_idx2)%alpha)*erfc(epsilon)))
                        pgf_radius_box = pgf_radius_box + INT(1.0_dp/&
                                (SQRT(2.0_dp*x_data%far_field_boxes(box_idx3, box_idx4)%alpha)*erfc(epsilon))) + 2
                        IF(Rbox>pgf_radius_box) THEN
                          max_val2 = 0.0_dp
                          SELECT CASE(screening_parameter%screening_type)
                            CASE(do_hf_screening_none)
                              Kbox = x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                     x_data%far_field_boxes(box_idx3, box_idx4)%K
                              max_val2 = 1.0_dp/Rbox*Kbox
                              IF( do_p_screening) THEN
                                max_val2 = max_val2 * MAX(pmax_set(iset,kset,iatom,katom), pmax_set(jset,lset,jatom,latom),&
                                           pmax_set(iset,lset,iatom,latom) , pmax_set(kset,jset,katom,jatom))
                              END IF
                            CASE(do_hf_screening_short)
                              omega_box = 1.0_dp/x_data%far_field_boxes(box_idx1, box_idx2)%alpha + &
                                          1.0_dp/x_data%far_field_boxes(box_idx3, box_idx4)%alpha + &
                                          1.0_dp/(screening_parameter%omega**2)
                              omega_box = 1.0_dp/omega_box
                              Kbox = x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                     x_data%far_field_boxes(box_idx3, box_idx4)%K
                              max_val2 = erfc(omega_box*Rbox)/Rbox*Kbox
                              IF( do_p_screening) THEN
                                max_val2 = max_val2 * MAX(pmax_set(iset,kset,iatom,katom), pmax_set(jset,lset,jatom,latom),&
                                           pmax_set(iset,lset,iatom,latom) , pmax_set(kset,jset,katom,jatom))
                              END IF
                          END SELECT
                          IF(max_val2<eps_schwarz) CYCLE
                        END IF
                        current_counter = current_counter + nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
                IF( current_counter >= int_per_cpu) THEN
                  tmp_dist(mepos,tmp_pos(mepos))%number_of_integrals = distribution_counter_end - distribution_counter_start +1
                  distribution_counter_start = distribution_counter_end + 1
                  restore_counter = current_counter
                  current_counter = 0
                  tmp_pos(mepos) = tmp_pos(mepos) + 1
                  mepos = mepos + 1
                  restore_iatom = iatom
                  restore_jatom = jatom
                  restore_katom = katom
                  restore_latom = latom  
                  IF( mepos>ncpu) THEN
                    mepos = mepos - ncpu
                  END IF
                 
                  tmp_dist(mepos,tmp_pos(mepos))%istart= iatom
                  tmp_dist(mepos,tmp_pos(mepos))%jstart= jatom
                  tmp_dist(mepos,tmp_pos(mepos))%kstart= katom
                  tmp_dist(mepos,tmp_pos(mepos))%lstart= latom + 1 !overflow should be recorvering automatically
                END IF
              END IF
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      !fill the last proc
      IF(ncpu > 1 ) THEN
        tmp_dist(ncpu,tmp_pos(ncpu))%istart= restore_iatom
        tmp_dist(ncpu,tmp_pos(ncpu))%jstart= restore_jatom
        tmp_dist(ncpu,tmp_pos(ncpu))%kstart= restore_katom
        tmp_dist(ncpu,tmp_pos(ncpu))%lstart= restore_latom + 1 !overflow should be recorvering automatically
      END IF
      IF(ncpu > 0 ) THEN
        IF(current_counter == 0 ) THEN
          tmp_dist(ncpu,tmp_pos(ncpu))%number_of_integrals = distribution_counter_end-distribution_counter_start+1 
        ELSE
          tmp_dist(ncpu,tmp_pos(ncpu))%number_of_integrals = distribution_counter_end-distribution_counter_start+1 
        END IF
      END IF
      DO mepos=1,ncpu
        ptr_to_tmp_dist => tmp_dist(mepos,:)
        IF(para_env%mepos == mepos-1)&
           CALL x_hf_set_distr_energy(ptr_to_tmp_dist,x_data,error=error)
      END DO
      DEALLOCATE(tmp_dist,tmp_pos,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle)
  END SUBROUTINE load_balance_energy
!***************************************************************************

!!***** four_center_integrals/load_balance_forces [1.0] *
!!
!!   NAME
!!     load_balance_forces
!!
!!   FUNCTION
!!     Distributes the computation of eri-derivatives to all available processes.
!!
!!   NOTES
!!     The optimization is done via a simple binning procedure:
!!     In a first step the total amount of integrals in the system is calculated,
!!     taking into account the sparsity of the KS-matrix, the screening based
!!     on the Schwarz-Inequality and the density-matrix. Then, the procedure uses 
!!     different bin-sizes in order to optimize the load balance.
!!     At the end each process owns an unique array of inidices-ranges 
!!     that are used to decide whether a process has to calculate a certain
!!     bunch of integrals.
!!     
!!
!!   INPUTS
!!     - x_data: Object that stores the indices array
!!     - eps_schwarz: screening parameter
!!     - particle_set, atomic_kind_set, para_env
!!     - screening_parameter: short/longrange screening information 
!!     - screen_set, screen_atom, screen_max, pmax_set: screening matrices
!!     - is_assoc_atomic_kind: KS-matrix sparsity
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     06.2007 created [Manuel Guidon]
!!     08.2007 new parallel scheme [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE load_balance_forces(x_data, eps_schwarz, particle_set,atomic_kind_set,para_env,&
                                 screening_parameter, screen_set, screen_atom, screen_max,&
                                 pmax_set, is_assoc_atomic_block_global, ks_matrix, do_periodic, &
                                 error)
    TYPE(x_hf_type), POINTER                 :: x_data
    REAL(KIND=dp), INTENT(IN)                :: eps_schwarz
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(x_hf_screening_type), INTENT(IN)    :: screening_parameter
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: screen_set
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: screen_atom
    REAL(KIND=dp), INTENT(IN)                :: screen_max
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: pmax_set
    INTEGER, DIMENSION(:, :), POINTER        :: is_assoc_atomic_block_global
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    LOGICAL                                  :: do_periodic
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'load_balance_forces', &
      routineP = moduleN//':'//routineN

    INTEGER :: box_idx1, box_idx2, box_idx3, box_idx4, first_from, first_to, &
      handle, i, iatom, icpu, ipos, iset, jatom, jset, katom, kset, latom, &
      lset, max_set, mepos, natom, nseta, nsetb, nsetc, nsetd, &
      pgf_radius_box, restore_iatom, restore_jatom, restore_katom, &
      restore_latom, stat
    INTEGER(KIND=int_8) :: current_counter, distribution_counter_end, &
      distribution_counter_start, int_per_cpu, max_dev, min_dev, nbins, ncpu, &
      opt_max, opt_nbin, restore_counter, total_integrals, total_integrals2
    INTEGER(KIND=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: cost_per_cpu
    INTEGER(KIND=int_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: cost_array
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: tmp_pos
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: first_dist
    INTEGER, DIMENSION(:), POINTER           :: nsgfa, nsgfb, nsgfc, nsgfd
    LOGICAL                                  :: failure
    REAL(dp)                                 :: Kbox, omega_box, Rbox
    REAL(KIND=dp)                            :: max_val1, max_val2
    TYPE(gto_basis_set_type), POINTER        :: orb_basis
    TYPE(x_hf_distribution), DIMENSION(:), &
      POINTER                                :: ptr_to_tmp_dist
    TYPE(x_hf_distribution), &
      DIMENSION(:, :), POINTER               :: tmp_dist

    natom=SIZE(particle_set)

    CALL timeset(routineN,handle)
    failure=.FALSE.
    ncpu = para_env%num_pe
    IF(ncpu == 1) THEN
      ALLOCATE(tmp_dist(1,1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      tmp_dist(1,1)%number_of_integrals = HUGE(tmp_dist(1,1)%number_of_integrals)
      tmp_dist(1,1)%istart = 1
      tmp_dist(1,1)%jstart = 1
      tmp_dist(1,1)%kstart = 1
      tmp_dist(1,1)%lstart = 1
      ptr_to_tmp_dist => tmp_dist(1,:)
      CALL x_hf_set_distr_forces(ptr_to_tmp_dist,x_data,error=error)
      DEALLOCATE(tmp_dist,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
      !Simple distribution on first loop variable
      ALLOCATE(first_dist(ncpu,2), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      first_dist = -1
      first_from = 1
      first_to = 1
      outer: DO icpu=1,ncpu
        IF(first_from > natom ) EXIT outer
        first_to = first_from + natom/ncpu + 1
        first_to = MIN(first_to, natom)
        first_dist(icpu,1) = first_from
        first_dist(icpu,2) = first_to
        IF(first_to == natom) EXIT outer
        first_from = first_to + 1
       END DO outer
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,maxnset=max_set)
 
      !Calculate total numbers of integrals that have to be calculated (wrt screening and symmetry)
      total_integrals = 0
      distribution_counter_end = 0
      IF(first_dist(para_env%mepos+1,1) /= -1 ) THEN
        DO iatom=first_dist(para_env%mepos+1,1),first_dist(para_env%mepos+1,2)
          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, orb_basis_set=orb_basis)
          CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nseta, nsgf_set=nsgfa)
          DO jatom=iatom,natom 
            IF (screen_atom(iatom,jatom) * screen_max < eps_schwarz) CYCLE
            CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind, orb_basis_set=orb_basis)
            CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetb, nsgf_set=nsgfb)
            DO katom=1,natom
              CALL get_atomic_kind(atomic_kind=particle_set(katom)%atomic_kind, orb_basis_set=orb_basis)
              CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetc, nsgf_set=nsgfc)
              DO latom=katom,natom
                !All four centers equivalent => zero-contribution
                IF((iatom==jatom .AND. iatom==katom .AND. iatom==latom)) CYCLE
  
                IF (screen_atom(iatom,jatom) * screen_atom(katom,latom) < eps_schwarz) CYCLE
                IF(.NOT. (is_assoc_atomic_block_global(iatom,latom)  >=1 .AND. &
                          is_assoc_atomic_block_global(iatom,katom)  >=1 .AND. &
                          is_assoc_atomic_block_global(jatom,katom)  >=1 .AND. &
                          is_assoc_atomic_block_global(jatom,latom)  >=1 ) ) CYCLE
  
                IF((katom+latom<=iatom+jatom))  THEN
                  IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE
                
                  CALL get_atomic_kind(atomic_kind=particle_set(latom)%atomic_kind, orb_basis_set=orb_basis)
                  CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetd, nsgf_set=nsgfd)
  
                  DO iset = 1, nseta
                    DO jset = 1,nsetb
                      max_val1 = screen_set(iatom,jatom,iset,jset)
                      IF (max_val1 * screen_atom(katom,latom) < eps_schwarz) CYCLE
                      DO kset = 1,nsetc
                        DO lset = 1,nsetd
                          max_val2 = max_val1 * screen_set(katom,latom,kset,lset) * &
                                     2.0_dp*MAX(pmax_set(iset,kset,iatom,katom) * pmax_set(jset,lset,jatom,latom),&
                                                pmax_set(iset,lset,iatom,latom) * pmax_set(kset,jset,katom,jatom))
                          IF(max_val2<eps_schwarz) CYCLE
                          IF(.NOT. do_periodic) THEN
                            box_idx1 = (iatom-1)*natom+jatom
                            box_idx2 = (iset-1)*max_set+jset
                            box_idx3 = (katom-1)*natom+latom
                            box_idx4 = (kset-1)*max_set+lset
                            Rbox = (x_data%far_field_boxes(box_idx1,box_idx2)%center(1) - &
                                    x_data%far_field_boxes(box_idx3,box_idx4)%center(1))**2 + &
                                   (x_data%far_field_boxes(box_idx1,box_idx2)%center(2) - &
                                    x_data%far_field_boxes(box_idx3,box_idx4)%center(2))**2 + &
                                   (x_data%far_field_boxes(box_idx1,box_idx2)%center(3) - &
                                    x_data%far_field_boxes(box_idx3,box_idx4)%center(3))**2
                            Rbox = SQRT(Rbox)
                            Rbox = Rbox - x_data%far_field_boxes(box_idx1, box_idx2)%R - &
                                          x_data%far_field_boxes(box_idx3, box_idx4)%R
                            pgf_radius_box = INT(1.0_dp/ &
                                  (SQRT(2.0_dp*x_data%far_field_boxes(box_idx1, box_idx2)%alpha)*erfc(epsilon)))
                            pgf_radius_box = pgf_radius_box + INT(1.0_dp/ &
                                  (SQRT(2.0_dp*x_data%far_field_boxes(box_idx3, box_idx4)%alpha)*erfc(epsilon))) + 2
                            IF(Rbox>pgf_radius_box) THEN
                              max_val2 = 0.0_dp
                              SELECT CASE(screening_parameter%screening_type)
                                CASE(do_hf_screening_none)
                                  Kbox = x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                         x_data%far_field_boxes(box_idx3, box_idx4)%K
                                  max_val2 = 1.0_dp/Rbox*Kbox*&
                                             2.0_dp*MAX(pmax_set(iset,kset,iatom,katom) * pmax_set(jset,lset,jatom,latom),&
                                             pmax_set(iset,lset,iatom,latom) * pmax_set(kset,jset,katom,jatom))
                                CASE(do_hf_screening_short)
                                  omega_box = 1.0_dp/x_data%far_field_boxes(box_idx1, box_idx2)%alpha + &
                                              1.0_dp/x_data%far_field_boxes(box_idx3, box_idx4)%alpha + &
                                              1.0_dp/(screening_parameter%omega**2)
                                  omega_box = 1.0_dp/omega_box
                                  Kbox = x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                         x_data%far_field_boxes(box_idx3, box_idx4)%K
                                  max_val2 = erfc(omega_box*Rbox)/Rbox*Kbox* &
                                             2.0_dp*MAX(pmax_set(iset,kset,iatom,katom) * pmax_set(jset,lset,jatom,latom),&
                                             pmax_set(iset,lset,iatom,latom) * pmax_set(kset,jset,katom,jatom))
                              END SELECT
                              IF(max_val2<eps_schwarz) CYCLE
                            END IF
                          END IF
  
                          total_integrals = total_integrals + nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                END IF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      END IF
      CALL mp_sum_i8(total_integrals,para_env%group)
  
      opt_max = HUGE(opt_max)
      ALLOCATE(cost_array(11,11))
      DO ipos = 1,11
        cost_array(ipos,1) = 2**(ipos-1)
        cost_array(ipos,2) = 0
      END DO
      DO i = 1,11
        nbins = cost_array(i,1)
        IF(MOD(i,para_env%num_pe) /= para_env%mepos ) CYCLE
        !floor value
        total_integrals2 = total_integrals +(ncpu*nbins- MOD(total_integrals,ncpu*nbins))
        int_per_cpu = total_integrals2 / ncpu
        int_per_cpu = int_per_cpu/nbins
        ALLOCATE(cost_per_cpu(ncpu))
        cost_per_cpu = 0
        current_counter = 0
        restore_counter = 0
        ipos = 1
        mepos = 1
        distribution_counter_start = 1
        distribution_counter_end = 0
        DO iatom=1,natom
          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, orb_basis_set=orb_basis)
          CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nseta, nsgf_set=nsgfa)
          DO jatom=iatom,natom
            IF (screen_atom(iatom,jatom) * screen_max < eps_schwarz) CYCLE
            CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind, orb_basis_set=orb_basis)
            CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetb, nsgf_set=nsgfb)
            DO katom=1,natom
              CALL get_atomic_kind(atomic_kind=particle_set(katom)%atomic_kind, orb_basis_set=orb_basis)
              CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetc, nsgf_set=nsgfc)
              DO latom=katom,natom
                !All four centers equivalent => zero-contribution
                IF((iatom==jatom .AND. iatom==katom .AND. iatom==latom)) CYCLE
  
                IF (screen_atom(iatom,jatom) * screen_atom(katom,latom) < eps_schwarz) CYCLE
                IF(.NOT. (is_assoc_atomic_block_global(iatom,latom)  >=1 .AND. &
                          is_assoc_atomic_block_global(iatom,katom)  >=1 .AND. &
                          is_assoc_atomic_block_global(jatom,katom)  >=1 .AND. &
                          is_assoc_atomic_block_global(jatom,latom)  >=1 ) ) CYCLE
  
                IF((katom+latom<=iatom+jatom))  THEN
                  IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE
                  distribution_counter_end = distribution_counter_end + 1
  
                  CALL get_atomic_kind(atomic_kind=particle_set(latom)%atomic_kind, orb_basis_set=orb_basis)
                  CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetd, nsgf_set=nsgfd)
  
                  DO iset = 1, nseta
                    DO jset = 1,nsetb
                      max_val1 = screen_set(iatom,jatom,iset,jset)
                      IF (max_val1 * screen_atom(katom,latom) < eps_schwarz) CYCLE
                      DO kset = 1,nsetc
                        DO lset = 1,nsetd
                          max_val2 = max_val1 * screen_set(katom,latom,kset,lset) * &
                                   2.0_dp*MAX(pmax_set(iset,kset,iatom,katom) * pmax_set(jset,lset,jatom,latom),&
                                              pmax_set(iset,lset,iatom,latom) * pmax_set(kset,jset,katom,jatom))
                          IF(max_val2<eps_schwarz) CYCLE
                          IF(.NOT. do_periodic) THEN
                            box_idx1 = (iatom-1)*natom+jatom
                            box_idx2 = (iset-1)*max_set+jset
                            box_idx3 = (katom-1)*natom+latom
                            box_idx4 = (kset-1)*max_set+lset
                            Rbox = (x_data%far_field_boxes(box_idx1,box_idx2)%center(1) - &
                                    x_data%far_field_boxes(box_idx3,box_idx4)%center(1))**2 + &
                                   (x_data%far_field_boxes(box_idx1,box_idx2)%center(2) - &
                                    x_data%far_field_boxes(box_idx3,box_idx4)%center(2))**2 + &
                                   (x_data%far_field_boxes(box_idx1,box_idx2)%center(3) - &
                                   x_data%far_field_boxes(box_idx3,box_idx4)%center(3))**2
                            Rbox = SQRT(Rbox)
                            Rbox = Rbox - x_data%far_field_boxes(box_idx1, box_idx2)%R - &
                                          x_data%far_field_boxes(box_idx3, box_idx4)%R
                            pgf_radius_box = INT(1.0_dp/ &
                                  (SQRT(2.0_dp*x_data%far_field_boxes(box_idx1, box_idx2)%alpha)*erfc(epsilon)))
                            pgf_radius_box = pgf_radius_box + INT(1.0_dp/ &
                                  (SQRT(2.0_dp*x_data%far_field_boxes(box_idx3, box_idx4)%alpha)*erfc(epsilon))) + 2
                             IF(Rbox>pgf_radius_box) THEN
                              max_val2 = 0.0_dp
                              SELECT CASE(screening_parameter%screening_type)
                                CASE(do_hf_screening_none)
                                  Kbox = x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                         x_data%far_field_boxes(box_idx3, box_idx4)%K
                                  max_val2 = 1.0_dp/Rbox*Kbox * &
                                             2.0_dp*MAX(pmax_set(iset,kset,iatom,katom) * pmax_set(jset,lset,jatom,latom),&
                                             pmax_set(iset,lset,iatom,latom) * pmax_set(kset,jset,katom,jatom))
                                CASE(do_hf_screening_short)
                                  omega_box = 1.0_dp/x_data%far_field_boxes(box_idx1, box_idx2)%alpha + &
                                              1.0_dp/x_data%far_field_boxes(box_idx3, box_idx4)%alpha + &
                                              1.0_dp/(screening_parameter%omega**2)
                                  omega_box = 1.0_dp/omega_box
                                  Kbox = x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                         x_data%far_field_boxes(box_idx3, box_idx4)%K
                                  max_val2 = erfc(omega_box*Rbox)/Rbox*Kbox * &
                                             2.0_dp*MAX(pmax_set(iset,kset,iatom,katom) * pmax_set(jset,lset,jatom,latom),&
                                             pmax_set(iset,lset,iatom,latom) * pmax_set(kset,jset,katom,jatom))
                              END SELECT
                              IF(max_val2<eps_schwarz) CYCLE
                            END IF
                          END IF
                          current_counter = current_counter + nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                  IF( current_counter >= int_per_cpu) THEN
                    cost_per_cpu(mepos) = cost_per_cpu(mepos) + current_counter
                    distribution_counter_start = distribution_counter_end + 1
                    restore_counter = current_counter
                    current_counter = 0
                    mepos = mepos + 1
                    IF( mepos>ncpu) THEN
                      mepos = mepos - ncpu
                      ipos = ipos + 2
                    END IF
                  END IF
                END IF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        !fill the last proc
        IF(ncpu > 0 ) THEN
          IF(current_counter == 0 ) THEN
            cost_per_cpu(ncpu) = cost_per_cpu(ncpu) + restore_counter
          ELSE
            cost_per_cpu(ncpu) = cost_per_cpu(ncpu) + current_counter
          END IF
        END IF
       
        max_dev = 0
        min_dev = 0
        DO mepos = 1,ncpu 
          IF(total_integrals/ncpu>cost_per_cpu(mepos)) THEN
            IF(ABS(total_integrals/ncpu-cost_per_cpu(mepos))>min_dev) min_dev = ABS(total_integrals/ncpu-cost_per_cpu(mepos))
          ELSE
            IF(ABS(total_integrals/ncpu-cost_per_cpu(mepos))>max_dev) max_dev = ABS(total_integrals/ncpu-cost_per_cpu(mepos))
          END IF
        END DO
        cost_array(i,2) = max_dev
  
        IF(1000*max_dev/(total_integrals/ncpu)<opt_max) THEN
          opt_max = 1000*max_dev/(total_integrals/ncpu) 
          opt_nbin = nbins
        END IF
        DEALLOCATE(cost_per_cpu, STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      END DO
     
      CALL mp_sum_i8(cost_array(:,2), para_env%group)
      opt_max = HUGE(opt_max)
      DO i = 1,11
        IF(1000*cost_array(i,2)/(total_integrals/ncpu)<opt_max) THEN
          opt_max = 1000*cost_array(i,2)/(total_integrals/ncpu)
          opt_nbin = cost_array(i,1)
        END IF
      END DO
      !Call routine again for nbins = opt_nbin
      nbins = opt_nbin
      !floor value
      int_per_cpu = CEILING(REAL(total_integrals)/REAL(ncpu)/REAL(nbins),int_8)
      ALLOCATE(tmp_dist(ncpu,nbins),STAT=stat)
      tmp_dist(:,:)%istart = 0
      tmp_dist(:,:)%jstart = 0
      tmp_dist(:,:)%kstart = 0
      tmp_dist(:,:)%lstart = 0
      tmp_dist(:,:)%number_of_integrals = 0
      tmp_dist(1,1)%istart = 1
      tmp_dist(1,1)%jstart = 1
      tmp_dist(1,1)%kstart = 1
      tmp_dist(1,1)%lstart = 1
      ALLOCATE(tmp_pos(ncpu))
      tmp_pos = 1
  
      current_counter = 0
      restore_counter = 0
      ipos = 1
      mepos = 1
      distribution_counter_start = 1
      distribution_counter_end = 0
      DO iatom=1,natom
        CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, orb_basis_set=orb_basis)
        CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nseta, nsgf_set=nsgfa)
        DO jatom=iatom,natom
          IF (screen_atom(iatom,jatom) * screen_max < eps_schwarz) CYCLE
          CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind, orb_basis_set=orb_basis)
          CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetb, nsgf_set=nsgfb)
          DO katom=1,natom
            CALL get_atomic_kind(atomic_kind=particle_set(katom)%atomic_kind, orb_basis_set=orb_basis)
            CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetc, nsgf_set=nsgfc)
            DO latom=katom,natom
            !All four centers equivalent => zero-contribution
              IF((iatom==jatom .AND. iatom==katom .AND. iatom==latom)) CYCLE
  
              IF (screen_atom(iatom,jatom) * screen_atom(katom,latom) < eps_schwarz) CYCLE
              IF(.NOT. (is_assoc_atomic_block_global(iatom,latom)  >=1 .AND. &
                        is_assoc_atomic_block_global(iatom,katom)  >=1 .AND. &
                        is_assoc_atomic_block_global(jatom,katom)  >=1 .AND. &
                        is_assoc_atomic_block_global(jatom,latom)  >=1 ) ) CYCLE
  
              IF((katom+latom<=iatom+jatom))  THEN
                IF( ((iatom+jatom).EQ.(katom+latom) ) .AND.(katom<iatom)) CYCLE
                distribution_counter_end = distribution_counter_end + 1
  
                CALL get_atomic_kind(atomic_kind=particle_set(latom)%atomic_kind, orb_basis_set=orb_basis)
                CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetd, nsgf_set=nsgfd)
  
                DO iset = 1, nseta
                  DO jset = 1,nsetb
                    max_val1 = screen_set(iatom,jatom,iset,jset)
                    IF (max_val1 * screen_atom(katom,latom) < eps_schwarz) CYCLE
                    DO kset = 1,nsetc
                      DO lset = 1,nsetd
                        max_val2 = max_val1 * screen_set(katom,latom,kset,lset) * &
                                   2.0_dp*MAX(pmax_set(iset,kset,iatom,katom) * pmax_set(jset,lset,jatom,latom),&
                                              pmax_set(iset,lset,iatom,latom) * pmax_set(kset,jset,katom,jatom))
                        IF(max_val2<eps_schwarz) CYCLE
                        IF(.NOT. do_periodic) THEN
                          box_idx1 = (iatom-1)*natom+jatom
                          box_idx2 = (iset-1)*max_set+jset
                          box_idx3 = (katom-1)*natom+latom
                          box_idx4 = (kset-1)*max_set+lset
                          Rbox = (x_data%far_field_boxes(box_idx1,box_idx2)%center(1) - &
                                  x_data%far_field_boxes(box_idx3,box_idx4)%center(1))**2 + &
                                 (x_data%far_field_boxes(box_idx1,box_idx2)%center(2) - &
                                  x_data%far_field_boxes(box_idx3,box_idx4)%center(2))**2 + &
                                 (x_data%far_field_boxes(box_idx1,box_idx2)%center(3) - &
                                  x_data%far_field_boxes(box_idx3,box_idx4)%center(3))**2
                          Rbox = SQRT(Rbox)
                          Rbox = Rbox - x_data%far_field_boxes(box_idx1, box_idx2)%R - &
                                        x_data%far_field_boxes(box_idx3, box_idx4)%R
                          pgf_radius_box = INT(1.0_dp/ &
                                (SQRT(2.0_dp*x_data%far_field_boxes(box_idx1, box_idx2)%alpha)*erfc(epsilon)))
                          pgf_radius_box = pgf_radius_box + INT(1.0_dp/ &
                                (SQRT(2.0_dp*x_data%far_field_boxes(box_idx3, box_idx4)%alpha)*erfc(epsilon))) + 2
                          IF(Rbox>pgf_radius_box) THEN
                            max_val2 = 0.0_dp
                            SELECT CASE(screening_parameter%screening_type)
                              CASE(do_hf_screening_none)
                                Kbox = x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                       x_data%far_field_boxes(box_idx3, box_idx4)%K
                                max_val2 = 1.0_dp/Rbox*Kbox * &
                                           2.0_dp*MAX(pmax_set(iset,kset,iatom,katom) * pmax_set(jset,lset,jatom,latom),&
                                           pmax_set(iset,lset,iatom,latom) * pmax_set(kset,jset,katom,jatom))
                              CASE(do_hf_screening_short)
                                omega_box = 1.0_dp/x_data%far_field_boxes(box_idx1, box_idx2)%alpha + &
                                            1.0_dp/x_data%far_field_boxes(box_idx3, box_idx4)%alpha + &
                                            1.0_dp/(screening_parameter%omega**2)
                                omega_box = 1.0_dp/omega_box
                                Kbox = x_data%far_field_boxes(box_idx1, box_idx2)%K * &
                                       x_data%far_field_boxes(box_idx3, box_idx4)%K
                                max_val2 = erfc(omega_box*Rbox)/Rbox*Kbox * &
                                           2.0_dp*MAX(pmax_set(iset,kset,iatom,katom) * pmax_set(jset,lset,jatom,latom),&
                                           pmax_set(iset,lset,iatom,latom) * pmax_set(kset,jset,katom,jatom))
                            END SELECT
                            IF(max_val2<eps_schwarz) CYCLE
                          END IF
                        END IF
                        current_counter = current_counter + nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)
                        ENDDO
                    ENDDO
                  ENDDO
                ENDDO
                IF( current_counter >= int_per_cpu) THEN
                  tmp_dist(mepos,tmp_pos(mepos))%number_of_integrals = distribution_counter_end - distribution_counter_start +1
                  distribution_counter_start = distribution_counter_end + 1
                  restore_counter = current_counter
                  current_counter = 0
                  tmp_pos(mepos) = tmp_pos(mepos) + 1
                  mepos = mepos + 1
                  restore_iatom = iatom
                  restore_jatom = jatom
                  restore_katom = katom
                  restore_latom = latom
                  IF( mepos>ncpu) THEN
                    mepos = mepos - ncpu
                    ipos = ipos + 2
                  END IF
                  tmp_dist(mepos,tmp_pos(mepos))%istart= iatom
                  tmp_dist(mepos,tmp_pos(mepos))%jstart= jatom
                  tmp_dist(mepos,tmp_pos(mepos))%kstart= katom
                  tmp_dist(mepos,tmp_pos(mepos))%lstart= latom + 1 !overflow should be recorvering automatically
                END IF
              END IF
          ENDDO
          ENDDO
        ENDDO
      ENDDO
      tmp_dist(ncpu,tmp_pos(ncpu))%istart= restore_iatom
      tmp_dist(ncpu,tmp_pos(ncpu))%jstart= restore_jatom
      tmp_dist(ncpu,tmp_pos(ncpu))%kstart= restore_katom
      tmp_dist(ncpu,tmp_pos(ncpu))%lstart= restore_latom + 1 !overflow should be recorvering automatically
      IF(current_counter == 0 ) THEN
        tmp_dist(ncpu,tmp_pos(ncpu))%number_of_integrals = distribution_counter_end-distribution_counter_start+1
      ELSE
        tmp_dist(ncpu,tmp_pos(ncpu))%number_of_integrals = distribution_counter_end-distribution_counter_start+1
      END IF
      DO mepos=1,ncpu
        ptr_to_tmp_dist => tmp_dist(mepos,:)
        IF(para_env%mepos == mepos-1)&
          CALL x_hf_set_distr_forces(ptr_to_tmp_dist,x_data,error=error)
      END DO
      DEALLOCATE(tmp_dist,tmp_pos,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle)
  END SUBROUTINE load_balance_forces
!***************************************************************************

!!***** four_center_integrals/coulomb4_periodic [1.0] *
!!
!!   NAME
!!     coulomb4_periodic
!!
!!   FUNCTION
!!     calculates two-electron integrals of a quartet/shell using the library 
!!     lib_int in the case of periodic boundary conditions
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - r: position
!!     - z: zeta
!!     - l: angular momentum
!!     - npgf: number of primitive cartesian gaussian in actual shell
!!     - primitive_integrals: array of primitive_integrals
!!     - box_size, neighbor_cells: contains info for pbc
!!     - screening_parameter: contains info for libint 
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     06.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE coulomb4_periodic(ra, rb, rc, rd, zeta, zetb, zetc, zetd,&
                               la_min, la_max, lb_min, lb_max,&
                               lc_min, lc_max, ld_min, ld_max,&
                               npgfa, npgfb, npgfc, npgfd,&
                               ncoa,ncob,ncoc,ncod,&
                               primitive_integrals,&
                               screening_parameter, box_size, neighbor_cells, &
                               screen1, screen2, eps_schwarz, max_k)
  
    REAL(dp), INTENT(IN)                     :: ra(3), rb(3), rc(3), rd(3)
    INTEGER, INTENT(IN)                      :: npgfa, ld_max, ld_min, &
                                                lc_max, lc_min, lb_max, &
                                                lb_min, la_max, la_min, &
                                                npgfd, npgfc, npgfb
    REAL(dp), DIMENSION(1:npgfa), INTENT(IN) :: zeta
    REAL(dp), DIMENSION(1:npgfb), INTENT(IN) :: zetb
    REAL(dp), DIMENSION(1:npgfc), INTENT(IN) :: zetc
    REAL(dp), DIMENSION(1:npgfd), INTENT(IN) :: zetd
    INTEGER, INTENT(IN)                      :: ncoa, ncob, ncoc, ncod
    REAL(dp), &
      DIMENSION(ncoa, ncob, ncoc, ncod)      :: primitive_integrals
    TYPE(x_hf_screening_type)                :: screening_parameter
    REAL(dp), INTENT(IN)                     :: box_size(3)
    TYPE(x_hf_cell_type), DIMENSION(:), &
      POINTER                                :: neighbor_cells
    REAL(dp), DIMENSION(:)                   :: screen1, screen2
    REAL(dp), INTENT(IN)                     :: eps_schwarz
    INTEGER, INTENT(IN)                      :: max_k

    INTEGER :: i, ipgf, j, jpgf, k, kpgf, la, lb, lc, ld, lpgf, offset_a, &
      offset_a1, offset_b, offset_b1, offset_c, offset_c1, offset_d, offset_d1
    REAL(dp) :: B11(3), C11(3), D11(3), Eta, im_B(3), im_D(3), P(3), &
      pbc_B(3), pbc_D(3), PQ(3), Q(3), shift(3), shift1(3), temp(3), &
      tmp_D(3), Zeta1, Zeta_A, Zeta_B, Zeta_C, Zeta_D

    primitive_integrals = 0.0_dp
    temp = rb -ra
    CALL pbc(temp, box_size, pbc_B)
    B11 = ra + pbc_B
    temp = rd-rc
    CALL pbc(temp, box_size, pbc_D)
    D11 = rc + pbc_D
    DO i = 1,SIZE(neighbor_cells)
      im_B(1) = B11(1) + neighbor_cells(i)%cell(1) * box_size(1)
      im_B(2) = B11(2) + neighbor_cells(i)%cell(2) * box_size(2)
      im_B(3) = B11(3) + neighbor_cells(i)%cell(3) * box_size(3)
      DO j = 1,SIZE(neighbor_cells)
        IF(screen1(i)*screen2(j)<eps_schwarz) CYCLE
        im_D(1) = D11(1) + neighbor_cells(j)%cell(1) * box_size(1)
        im_D(2) = D11(2) + neighbor_cells(j)%cell(2) * box_size(2)
        im_D(3) = D11(3) + neighbor_cells(j)%cell(3) * box_size(3)
        DO lpgf = 1,npgfd
          offset_d1 = (lpgf-1)*ncoset(ld_max)
          Zeta_D = zetd(lpgf)
          DO ld = ld_min, ld_max
            offset_d = offset_d1 + ncoset(ld-1)
            DO kpgf = 1,npgfc
              offset_c1 = (kpgf-1)*ncoset(lc_max)
              Zeta_C = zetc(kpgf)
              Eta  =  Zeta_C + Zeta_D
              Q = (Zeta_C*rc + Zeta_D*im_D)/Eta
              DO lc = lc_min, lc_max
                offset_c = offset_c1 + ncoset(lc-1)
                DO jpgf = 1,npgfb
                  offset_b1 = (jpgf-1)*ncoset(lb_max)
                  Zeta_B = zetb(jpgf)
                  DO lb = lb_min, lb_max
                    offset_b = offset_b1 + ncoset(lb-1)
                    DO ipgf = 1,npgfa
                      offset_a1 = (ipgf-1)*ncoset(la_max)
                      Zeta_A = zeta(ipgf)
                      Zeta1 = Zeta_A + Zeta_B
                      P = (Zeta_A*ra + Zeta_B*im_B)/Zeta1
                      temp = P-Q
                      CALL pbc(temp,box_size,PQ)
                      shift = -PQ + temp
                      DO k=1,max_k!size(neighbor_cells)
                        shift1(1) =  neighbor_cells(k)%cell(1) * box_size(1)
                        shift1(2) =  neighbor_cells(k)%cell(2) * box_size(2)
                        shift1(3) =  neighbor_cells(k)%cell(3) * box_size(3)
                        C11 = rc + shift + shift1
                        tmp_D = im_D + shift + shift1
                        DO la = la_min,la_max
                          offset_a = offset_a1 + ncoset(la-1)
                          !Build primitives
                          CALL evaluate_eri(lib, ra, im_B, C11, tmp_D,&
                                            zeta(ipgf), zetb(jpgf), zetc(kpgf), zetd(lpgf),&
                                            la, lb, lc ,ld,&
                                            ncoa,ncob,ncoc,ncod,&
                                            offset_a,offset_b,offset_c,offset_d, &
                                            primitive_integrals, &
                                            screening_parameter)
                        END DO !la
                      END DO
                    END DO !ipgf
                  END DO !lb
                END DO !jpgf
              END DO !lc
            END DO !kpgf
          END DO !ld
        END DO !lpgf
      END DO
    END DO
  END SUBROUTINE coulomb4_periodic
!***************************************************************************

!!***** four_center_integrals/comp_screen_mat_per [1.0] *
!!
!!   NAME
!!     comp_screen_mat_per
!!
!!   FUNCTION
!!     calculates the screening matrices in the case of periodic boundary
!!     conditions
!!
!!   NOTES
!!     - already computes the square root
!!     - also takes the primitive to sgf transformation into account
!!     - in this way all stored contracted integrals are smaller or equal to their estimate
!!
!!   INPUTS
!!     - screen_xxx: screening matrices
!!     - read_from_memory: Flag to control reading (when true, data come from buffer)
!!     - buffer: stores the matrices
!!     - neighbor_cells, box_size: contain information fort pbc
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!    06.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE comp_screen_mat_per(screen_set, screen_atom, screen_max,screen_periodic,particle_set,&
                                 atomic_kind_set,screening_parameter, read_from_memory, buffer,&
                                 neighbor_cells, box_size, max_k, error)
    REAL(kind=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: screen_set
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: screen_atom
    REAL(KIND=dp)                            :: screen_max
    REAL(kind=dp), &
      DIMENSION(:, :, :, :, :), POINTER      :: screen_periodic
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(x_hf_screening_type)                :: screening_parameter
    LOGICAL                                  :: read_from_memory
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), OPTIONAL, &
      POINTER                                :: buffer
    TYPE(x_hf_cell_type), DIMENSION(:), &
      POINTER                                :: neighbor_cells
    REAL(KIND=dp)                            :: box_size(3)
    INTEGER, INTENT(IN)                      :: max_k
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'comp_screen_mat_per', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, jatom, jset, k, &
                                                katom, kset, max_set, natom, &
                                                ncob, nsetb, nsetc, sgfb, stat
    INTEGER, DIMENSION(:), POINTER           :: lb_max, lb_min, lc_max, &
                                                lc_min, npgfb, npgfc, nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfb
    LOGICAL                                  :: failure
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: max_val
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: contraction_max
    REAL(dp), DIMENSION(:, :), POINTER       :: sphi_b, zetb, zetc
    REAL(KIND=dp)                            :: p2s_b, p2s_c, rb(3), rc(3)
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_b, orb_basis_c

    CALL timeset(routineN,handle)

    failure = .FALSE.

    natom=SIZE(particle_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, maxnset=max_set)
    ALLOCATE(screen_set(natom,natom,max_set,max_set),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(screen_atom(natom,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(contraction_max(max_set,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    screen_set=0.0_dp
    screen_atom=0.0_dp
    contraction_max=0.0_dp
    screen_periodic = 0.0_dp
    ALLOCATE(max_val(SIZE(neighbor_cells)))
    max_val = 0.0_dp
    
    DO jatom=1,natom

      atom_kind => particle_set(jatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_b)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                             lmax=lb_max,&
                             npgf=npgfb,&
                             nset=nsetb,&
                             first_sgf=first_sgfb,&
                             sphi=sphi_b,&
                             nsgf_set=nsgfb)
      DO jset = 1,nsetb
        ! takes the primitive to contracted transformation into account
        ncob = npgfb(jset)*ncoset(lb_max(jset))
        sgfb = first_sgfb(1,jset)
        ! if the primitives are assumed to be all of max_val2, max_val2*p2s_b becomes
        ! the maximum value after multiplication with sphi_b
        contraction_max(jset,jatom) = MAXVAL((/(SUM(ABS(sphi_b(1:ncob,k))),k=sgfb,sgfb+nsgfb(jset)-1)/))
      ENDDO
    ENDDO

    DO jatom=1,natom

      atom_kind => particle_set(jatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_b)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                             lmax=lb_max,&
                             lmin=lb_min,&
                             npgf=npgfb,&
                             nset=nsetb,&
                             zet=zetb)
      rb=particle_set(jatom)%r(:)
      DO katom = 1,natom
        atom_kind => particle_set(katom)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atom_kind,&
                             orb_basis_set=orb_basis_c)
        CALL get_gto_basis_set(gto_basis_set=orb_basis_c,&
                               lmax=lc_max,&
                               lmin=lc_min,&
                               npgf=npgfc,&
                               nset=nsetc,&
                               zet=zetc)
        rc=particle_set(katom)%r(:)
        DO jset = 1,nsetb
          p2s_b = contraction_max(jset,jatom)
          DO kset = 1,nsetc
            p2s_c = contraction_max(kset,katom)
            IF( read_from_memory ) THEN
              screen_periodic(jatom,katom,jset,kset,:) = buffer(jatom,katom,jset,kset,:)
            ELSE
              max_val = 0.0_dp
              CALL screen4_periodic(rb, rc, &
                                    zetb(:,jset), zetc(:,kset), &
                                    lb_min(jset), lb_max(jset), lc_min(kset), lc_max(kset),&
                                    npgfb(jset), npgfc(kset), &
                                    max_val, screening_parameter, neighbor_cells, box_size, max_k)
                screen_periodic(jatom,katom,jset,kset,:) = SQRT(max_val)*p2s_b*p2s_c
                IF( PRESENT(buffer) ) buffer(jatom,katom,jset,kset,:) = screen_periodic(jatom,katom,jset,kset,:)
              END IF
              screen_set(jatom,katom,jset,kset) = MAXVAL(screen_periodic(jatom,katom,jset,kset,:))
          END DO
        END DO
        screen_atom(jatom,katom)=MAXVAL(screen_set(jatom,katom,1:nsetb,1:nsetc))
      END DO
    END DO

    screen_max=MAXVAL(screen_atom)
    DEALLOCATE(contraction_max,max_val, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle)

  END SUBROUTINE comp_screen_mat_per
!***************************************************************************

!!***** four_center_integrals/screen4_periodic [1.0] *
!!
!!   NAME
!!     screen4_periodic
!!
!!   FUNCTION
!!     calculates max values of two-electron integrals in a quartet/shell
!!     w.r.t. different zetas using the library lib_int in the case of 
!!     periodic boundary conditions
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - r: position
!!     - z: zeta
!!     - l: angular momentum
!!     - npgf: number of primitive cartesian gaussian in actual shell
!!     - max_val: schwarz screening value
!!     - screening_parameter: contains info for libint
!!     - neighbor_cells: Array containing the considered cells
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!    06.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE screen4_periodic(ra, rb, zeta, zetb, &
                             la_min, la_max, lb_min, lb_max,&
                             npgfa, npgfb, &
                             max_val, screening_parameter, neighbor_cells, &
                             box_size, max_k)

    REAL(dp), INTENT(IN)                     :: ra(3), rb(3)
    REAL(dp), DIMENSION(:), INTENT(IN)       :: zeta, zetb
    INTEGER, INTENT(IN)                      :: la_min, la_max, lb_min, &
                                                lb_max, npgfa, npgfb
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: max_val
    TYPE(x_hf_screening_type)                :: screening_parameter
    TYPE(x_hf_cell_type), DIMENSION(:), &
      POINTER                                :: neighbor_cells
    REAL(KIND=dp), INTENT(IN)                :: box_size(3)
    INTEGER, INTENT(IN)                      :: max_k

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'screen4_periodic'

    INTEGER                                  :: i, ipgf, jpgf, k, la, lb
    REAL(dp)                                 :: B11(3), im_B(3), &
                                                max_val_temp, pbc_B(3), &
                                                shift1(3), temp(3)

    temp = rb -ra
    CALL pbc(temp, box_size, pbc_B)
    B11 = ra + pbc_B
    DO i = 1,SIZE(neighbor_cells)
      im_B(1) = B11(1) + neighbor_cells(i)%cell(1) * box_size(1)
      im_B(2) = B11(2) + neighbor_cells(i)%cell(2) * box_size(2)
      im_B(3) = B11(3) + neighbor_cells(i)%cell(3) * box_size(3)
      DO k=1,max_k!size(neighbor_cells)
        shift1(1) =  neighbor_cells(k)%cell(1) * box_size(1)
        shift1(2) =  neighbor_cells(k)%cell(2) * box_size(2)
        shift1(3) =  neighbor_cells(k)%cell(3) * box_size(3)
        DO ipgf = 1,npgfa
          DO jpgf = 1,npgfb
            DO la = la_min,la_max
              DO lb = lb_min, lb_max
                !Build primitives
                 max_val_temp = 0.0_dp
                 CALL evaluate_eri_screen(lib, ra, im_B, ra, im_B,&
                                          zeta(ipgf), zetb(jpgf), zeta(ipgf), zetb(jpgf),&
                                          la, lb, la ,lb,&
                                           max_val_temp, screening_parameter)
                 max_val(i) = MAX(max_val(i), max_val_temp)
              END DO !lb
            END DO !la
          END DO !jpgf
        END DO !ipgf
      END DO
    END DO
  END SUBROUTINE screen4_periodic
!***************************************************************************

!!***** four_center_integrals/forces4_periodic [1.0] *
!!
!!   NAME
!!     forces4_periodic
!!
!!   FUNCTION
!!     calculates two-electron derivatived of a quartet/shell using the library 
!!     lib_deriv wrt pbc
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - r: position
!!     - z: zeta
!!     - l: angular momentum
!!     - npgf: number of primitive cartesian gaussian in actual shell
!!     - primitive_forces: array of primitive_forces
!!     - work_forces: work_storage
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     07.2007 created [Manuel Guidon]
!!
!!****************************************************************************

  SUBROUTINE forces4_periodic(ra, rb, rc, rd, zeta, zetb, zetc, zetd,&
                              la_min, la_max, lb_min, lb_max,&
                              lc_min, lc_max, ld_min, ld_max,&
                              npgfa, npgfb, npgfc, npgfd,&
                              ncoa, ncob, ncoc, ncod,&
                              primitive_forces, work_forces,&
                              screening_parameter, box_size, neighbor_cells, screen1, screen2, eps_schwarz, max_k)
  
    REAL(dp), INTENT(IN)                     :: ra(3), rb(3), rc(3), rd(3)
    REAL(dp), DIMENSION(:), INTENT(IN)       :: zeta, zetb, zetc, zetd
    INTEGER, INTENT(IN) :: la_min, la_max, lb_min, lb_max, lc_min, lc_max, &
      ld_min, ld_max, npgfa, npgfb, npgfc, npgfd, ncoa, ncob, ncoc, ncod
    REAL(dp), &
      DIMENSION(ncoa, ncob, ncoc, ncod, 12)  :: primitive_forces
    REAL(dp), DIMENSION(:, :), POINTER       :: work_forces
    TYPE(x_hf_screening_type)                :: screening_parameter
    REAL(KIND=dp), INTENT(IN)                :: box_size(3)
    TYPE(x_hf_cell_type), DIMENSION(:), &
      POINTER                                :: neighbor_cells
    REAL(dp), DIMENSION(:)                   :: screen1, screen2
    REAL(dp), INTENT(IN)                     :: eps_schwarz
    INTEGER, INTENT(IN)                      :: max_k

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'forces4_periodic'

    INTEGER :: i, ipgf, j, jpgf, k, kpgf, la, lb, lc, ld, lpgf, offset_a, &
      offset_a1, offset_b, offset_b1, offset_c, offset_c1, offset_d, offset_d1
    REAL(dp) :: B11(3), C11(3), D11(3), Eta, im_B(3), im_D(3), P(3), &
      pbc_B(3), pbc_D(3), PQ(3), Q(3), shift(3), shift1(3), temp(3), &
      tmp_D(3), Zeta1, Zeta_A, Zeta_B, Zeta_C, Zeta_D

    primitive_forces = 0.0_dp
    temp = rb -ra
    CALL pbc(temp, box_size, pbc_B)
    B11 = ra + pbc_B
    temp = rd-rc
    CALL pbc(temp, box_size, pbc_D)
    D11 = rc + pbc_D
    DO i = 1,SIZE(neighbor_cells)
      im_B(1) = B11(1) + neighbor_cells(i)%cell(1) * box_size(1)
      im_B(2) = B11(2) + neighbor_cells(i)%cell(2) * box_size(2)
      im_B(3) = B11(3) + neighbor_cells(i)%cell(3) * box_size(3)
      DO j = 1,SIZE(neighbor_cells)
        IF(screen1(i)*screen2(j)<eps_schwarz) CYCLE
        im_D(1) = D11(1) + neighbor_cells(j)%cell(1) * box_size(1)
        im_D(2) = D11(2) + neighbor_cells(j)%cell(2) * box_size(2)
        im_D(3) = D11(3) + neighbor_cells(j)%cell(3) * box_size(3)
        DO lpgf = 1,npgfd
          offset_d1 = (lpgf-1)*ncoset(ld_max)
          Zeta_D = zetd(lpgf)
          DO ld = ld_min, ld_max
            offset_d = offset_d1 + ncoset(ld-1)
            DO kpgf = 1,npgfc
              offset_c1 = (kpgf-1)*ncoset(lc_max)
              Zeta_C = zetc(kpgf)
              Eta  =  Zeta_C + Zeta_D
              Q = (Zeta_C*rc + Zeta_D*im_D)/Eta
              DO lc = lc_min, lc_max
                offset_c = offset_c1 + ncoset(lc-1)
                DO jpgf = 1,npgfb
                  offset_b1 = (jpgf-1)*ncoset(lb_max)
                  Zeta_B = zetb(jpgf)
                  DO lb = lb_min, lb_max
                    offset_b = offset_b1 + ncoset(lb-1)
                    DO ipgf = 1,npgfa
                      offset_a1 = (ipgf-1)*ncoset(la_max)
                      Zeta_A = zeta(ipgf)
                      Zeta1 = Zeta_A + Zeta_B
                      P = (Zeta_A*ra + Zeta_B*im_B)/Zeta1
                      temp = P-Q
                      CALL pbc(temp,box_size,PQ)
                      shift = -PQ + temp
                      DO k=1,max_k!size(neighbor_cells)
                        shift1(1) =  neighbor_cells(k)%cell(1) * box_size(1)
                        shift1(2) =  neighbor_cells(k)%cell(2) * box_size(2)
                        shift1(3) =  neighbor_cells(k)%cell(3) * box_size(3)
                     
                        C11 = rc + shift + shift1
                        tmp_D = im_D + shift + shift1
                        DO la = la_min,la_max
                          offset_a = offset_a1 + ncoset(la-1)
                          !Build primitives
                          CALL evaluate_deriv_eri(deriv, ra, im_B, C11, tmp_D,&
                                                  zeta(ipgf), zetb(jpgf), zetc(kpgf), zetd(lpgf),&
                                                  la, lb, lc ,ld,&
                                                  work_forces,&
                                                  ncoa, ncob, ncoc, ncod, &
                                                  primitive_forces, &
                                                  offset_a, offset_b, offset_c, offset_d,&
                                                  screening_parameter)
                        END DO !la
                      END DO
                    END DO !ipgf
                  END DO !lb
                END DO !jpgf
              END DO !lc
            END DO !kpgf
          END DO !ld
        END DO !lpgf
      END DO
    END DO
  END SUBROUTINE forces4_periodic
!***************************************************************************

!!***** four_center_integrals/screen_farfield [1.0] *
!!
!!   NAME
!!     screen_farfield
!!
!!   FUNCTION
!!     Estimates the contribution of an ERI taking 1/R dependence into account.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - boxes: Buffer to store screening matrices
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     08.2007 created [Manuel Guidon]
!!
!!****************************************************************************
 
  SUBROUTINE screen_farfield(boxes, particle_set, atomic_kind_set, screening_parameter, error)
    TYPE(x_hf_box_scheme), DIMENSION(:, :), &
      POINTER                                :: boxes
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(x_hf_screening_type)                :: screening_parameter
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'screen_farfield', &
      routineP = moduleN//':'//routineN

    INTEGER :: box_idx1, box_idx2, handle, jatom, jpgf, jset, k, katom, kpgf, &
      kset, max_set, natom, ncob, nsetb, nsetc, pgf_radius, sgfb, stat
    INTEGER, DIMENSION(:), POINTER           :: lb_max, lb_min, lc_max, &
                                                lc_min, npgfb, npgfc, nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfb
    LOGICAL                                  :: failure
    REAL(dp)                                 :: px_max, px_min, py_max, &
                                                py_min, pz_max, pz_min, rbc2
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: contraction_max
    REAL(dp), DIMENSION(:, :), POINTER       :: rpgfb, rpgfc, sphi_b, zetb, &
                                                zetc
    REAL(KIND=dp) :: alpha, K_max, Kbc, p2s_b, p2s_c, R_max, rb(3), Rbc, &
      rc(3), rp(3), system_size(3), x_max, x_min, y_max, y_min, z_max, z_min, &
      zet1, zet2
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_b, orb_basis_c

    CALL timeset(routineN,handle)

    failure = .FALSE.

    natom=SIZE(particle_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, maxnset=max_set)
   
    ALLOCATE(contraction_max(max_set,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    contraction_max=0.0_dp

    x_min = HUGE(x_min)
    x_max =-HUGE(x_max)
    y_min = HUGE(y_min)
    y_max = -HUGE(y_max)
    z_min = HUGE(z_min)
    z_max = -HUGE(z_max)
    DO jatom=1,natom
      rb=particle_set(jatom)%r(:)
      IF(rb(1)>x_max) x_max=rb(1)
      IF(rb(1)<x_min) x_min=rb(1)
      IF(rb(2)>y_max) y_max=rb(2)
      IF(rb(2)<y_min) y_min=rb(2)
      IF(rb(3)>z_max) z_max=rb(3)
      IF(rb(3)<z_min) z_min=rb(3)
    END DO

    system_size(1) = x_max - x_min
    system_size(2) = y_max - y_min
    system_size(3) = z_max - z_min
    DEALLOCATE(boxes,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(boxes(natom**2,max_set**2), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    boxes(:,:)%center(1) = 0.0_dp
    boxes(:,:)%center(2) = 0.0_dp
    boxes(:,:)%center(3) = 0.0_dp
    boxes(:,:)%K = -HUGE(Kbc)
    boxes(:,:)%R = -HUGE(Rbc)
    boxes(:,:)%alpha = HUGE(alpha)
    DO jatom=1,natom
      atom_kind => particle_set(jatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_b)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                             lmax=lb_max,&
                             npgf=npgfb,&
                             nset=nsetb,&
                             first_sgf=first_sgfb,&
                             sphi=sphi_b,&
                             nsgf_set=nsgfb)
      DO jset = 1,nsetb
        ! takes the primitive to contracted transformation into account
        ncob = npgfb(jset)*ncoset(lb_max(jset))
        sgfb = first_sgfb(1,jset)
        ! if the primitives are assumed to be all of max_val2, max_val2*p2s_b becomes
        ! the maximum value after multiplication with sphi_b
        contraction_max(jset,jatom) = MAXVAL((/(SUM(ABS(sphi_b(1:ncob,k))),k=sgfb,sgfb+nsgfb(jset)-1)/))
      ENDDO
    ENDDO

    DO jatom=1,natom
      atom_kind => particle_set(jatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_b)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                             lmax=lb_max,&
                             lmin=lb_min,&
                             npgf=npgfb,&
                             nset=nsetb,&
                             pgf_radius=rpgfb,&
                             zet=zetb)
      rb=particle_set(jatom)%r(:)
      DO katom = 1,natom
        atom_kind => particle_set(katom)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atom_kind,&
                             orb_basis_set=orb_basis_c)
        CALL get_gto_basis_set(gto_basis_set=orb_basis_c,&
                               lmax=lc_max,&
                               lmin=lc_min,&
                               npgf=npgfc,&
                               nset=nsetc,&
                               pgf_radius=rpgfc,&
                               zet=zetc)
        rc=particle_set(katom)%r(:)
        rbc2 = (rb(1)-rc(1))**2+(rb(2)-rc(2))**2+(rb(3)-rc(3))**2
        box_idx1 = (jatom-1)*natom + katom 
        DO jset = 1,nsetb
          p2s_b = contraction_max(jset,jatom)
          DO kset = 1,nsetc
            p2s_c = contraction_max(kset,katom)
            K_max = TINY(K_max)
            R_max = TINY(R_max)
            px_min = HUGE(px_min)
            px_max = -HUGE(px_max)
            py_min = HUGE(py_min)
            py_max = -HUGE(py_max)
            pz_min = HUGE(pz_min)
            pz_max = -HUGE(pz_max)
            box_idx2 = (jset-1)* max_set + kset
            DO jpgf = 1,npgfb(jset)
              zet1 = zetb(jpgf,jset)
              DO kpgf = 1,npgfc(kset)
                zet2 = zetc(kpgf,kset)
                alpha = zet1+zet2
                rp = (zet1*rb+zet2*rc)/alpha
                IF(rp(1)<px_min) px_min = rp(1)
                IF(rp(2)<py_min) py_min = rp(2)
                IF(rp(3)<pz_min) pz_min = rp(3)
                IF(rp(1)>px_max) px_max = rp(1)
                IF(rp(2)>py_max) py_max = rp(2)
                IF(rp(3)>pz_max) pz_max = rp(3)
                Kbc = 1.0_dp/alpha*EXP(-zet1*zet2/alpha*(rbc2))
                IF( alpha < boxes(box_idx1, box_idx2)%alpha ) boxes(box_idx1, box_idx2)%alpha = alpha
                IF( Kbc > boxes(box_idx1, box_idx2)%K ) boxes(box_idx1, box_idx2)%K = Kbc*p2s_b*p2s_c*&
                    pi**1.25_dp*SQRT(2.0_dp)
              END DO
            END DO
            boxes(box_idx1, box_idx2)%center(1) = (px_max + px_min) / 2.0_dp
            boxes(box_idx1, box_idx2)%center(2) = (py_max + py_min) / 2.0_dp
            boxes(box_idx1, box_idx2)%center(3) = (pz_max + pz_min) / 2.0_dp
            R_max = TINY(R_max)
            DO jpgf = 1,npgfb(jset)
              zet1 = zetb(jpgf,jset)
              DO kpgf = 1,npgfc(kset)
                zet2 = zetc(kpgf,kset)
                alpha = zet1+zet2
                rp = (zet1*rb+zet2*rc)/alpha
                R_max = SQRT((rp(1)-boxes(box_idx1, box_idx2)%center(1))**2 + &
                             (rp(2)-boxes(box_idx1, box_idx2)%center(2))**2 + &
                             (rp(3)-boxes(box_idx1, box_idx2)%center(3))**2)
                IF( R_max > boxes(box_idx1, box_idx2)%R ) boxes(box_idx1, box_idx2)%R = R_max
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
    DEALLOCATE(contraction_max, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle)
  END SUBROUTINE screen_farfield
!***************************************************************************

!!***** four_center_integrals/screen_farfield_per [1.0] *
!!
!!   NAME
!!     screen_farfield_per
!!
!!   FUNCTION
!!     Estimates the contribution of an ERI taking 1/R dependence into account.
!!     Periodic case
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - boxes: Buffer to store screening matrices
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     09.2007 created [Manuel Guidon]
!!
!!****************************************************************************
 
  SUBROUTINE screen_farfield_per(boxes, particle_set, atomic_kind_set, screening_parameter, &
                                 neighbor_cells, box_size, error)
    TYPE(x_hf_box_scheme), DIMENSION(:, :), &
      POINTER                                :: boxes
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(x_hf_screening_type)                :: screening_parameter
    TYPE(x_hf_cell_type), DIMENSION(:), &
      POINTER                                :: neighbor_cells
    REAL(KIND=dp), INTENT(IN)                :: box_size(3)
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'screen_farfield_per', &
      routineP = moduleN//':'//routineN

    INTEGER :: box_idx1, box_idx2, handle, i, jatom, jpgf, jset, k, katom, &
      kpgf, kset, max_set, natom, ncob, nsetb, nsetc, sgfb, stat
    INTEGER, DIMENSION(:), POINTER           :: lb_max, lb_min, lc_max, &
                                                lc_min, npgfb, npgfc, nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfb
    LOGICAL                                  :: failure
    REAL(dp)                                 :: C11(3), im_C(3), pbc_C(3), &
                                                px_max, px_min, py_max, &
                                                py_min, pz_max, pz_min, rbc2, &
                                                temp(3)
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: contraction_max
    REAL(dp), DIMENSION(:, :), POINTER       :: rpgfb, rpgfc, sphi_b, zetb, &
                                                zetc
    REAL(KIND=dp)                            :: alpha, K_max, Kbc, p2s_b, &
                                                p2s_c, R_max, rb(3), Rbc, &
                                                rc(3), rp(3), zet1, zet2
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_b, orb_basis_c

    CALL timeset(routineN,handle)

    failure = .FALSE.

    natom=SIZE(particle_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, maxnset=max_set)
   
    ALLOCATE(contraction_max(max_set,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    contraction_max=0.0_dp

    DEALLOCATE(boxes,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(boxes(natom**2,max_set**2), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    boxes(:,:)%center(1) = 0.0_dp
    boxes(:,:)%center(2) = 0.0_dp
    boxes(:,:)%center(3) = 0.0_dp
    boxes(:,:)%K = -HUGE(Kbc)
    boxes(:,:)%R = -HUGE(Rbc)
    boxes(:,:)%alpha = HUGE(alpha)
    DO jatom=1,natom
      atom_kind => particle_set(jatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_b)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                             lmax=lb_max,&
                             npgf=npgfb,&
                             nset=nsetb,&
                             first_sgf=first_sgfb,&
                             sphi=sphi_b,&
                             nsgf_set=nsgfb)
      DO jset = 1,nsetb
        ! takes the primitive to contracted transformation into account
        ncob = npgfb(jset)*ncoset(lb_max(jset))
        sgfb = first_sgfb(1,jset)
        ! if the primitives are assumed to be all of max_val2, max_val2*p2s_b becomes
        ! the maximum value after multiplication with sphi_b
        contraction_max(jset,jatom) = MAXVAL((/(SUM(ABS(sphi_b(1:ncob,k))),k=sgfb,sgfb+nsgfb(jset)-1)/))
      ENDDO
    ENDDO

    DO jatom=1,natom
      atom_kind => particle_set(jatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_b)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                             lmax=lb_max,&
                             lmin=lb_min,&
                             npgf=npgfb,&
                             nset=nsetb,&
                             pgf_radius=rpgfb,&
                             zet=zetb)
      rb=particle_set(jatom)%r(:)
      DO katom = 1,natom
        atom_kind => particle_set(katom)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atom_kind,&
                             orb_basis_set=orb_basis_c)
        CALL get_gto_basis_set(gto_basis_set=orb_basis_c,&
                               lmax=lc_max,&
                               lmin=lc_min,&
                               npgf=npgfc,&
                               nset=nsetc,&
                               pgf_radius=rpgfc,&
                               zet=zetc)
        rc=particle_set(katom)%r(:)
        box_idx1 = (jatom-1)*natom + katom 
        temp = rc-rb
        CALL pbc(temp, box_size, pbc_C)
        C11 = rb + pbc_C
        DO jset = 1,nsetb
          p2s_b = contraction_max(jset,jatom)
          DO kset = 1,nsetc
            p2s_c = contraction_max(kset,katom)
            K_max = TINY(K_max)
            R_max = TINY(R_max)
            px_min = HUGE(px_min)
            px_max = -HUGE(px_max)
            py_min = HUGE(py_min)
            py_max = -HUGE(py_max)
            pz_min = HUGE(pz_min)
            pz_max = -HUGE(pz_max)
            box_idx2 = (jset-1)* max_set + kset
            DO i = 1,SIZE(neighbor_cells)
              im_C(1) = C11(1) + neighbor_cells(i)%cell(1) * box_size(1)
              im_C(2) = C11(2) + neighbor_cells(i)%cell(2) * box_size(2)
              im_C(3) = C11(3) + neighbor_cells(i)%cell(3) * box_size(3)
              rbc2 = (rb(1)-im_C(1))**2+(rb(2)-im_C(2))**2+(rb(3)-im_C(3))**2
              DO jpgf = 1,npgfb(jset)
                zet1 = zetb(jpgf,jset)
                DO kpgf = 1,npgfc(kset)
                  zet2 = zetc(kpgf,kset)
                  alpha = zet1+zet2
                  rp = (zet1*rb+zet2*im_C)/alpha
                  IF(rp(1)<px_min) px_min = rp(1)
                  IF(rp(2)<py_min) py_min = rp(2)
                  IF(rp(3)<pz_min) pz_min = rp(3)
                  IF(rp(1)>px_max) px_max = rp(1)
                  IF(rp(2)>py_max) py_max = rp(2)
                  IF(rp(3)>pz_max) pz_max = rp(3)
                  Kbc = 1.0_dp/alpha*EXP(-zet1*zet2/alpha*(rbc2))
                  IF( alpha < boxes(box_idx1, box_idx2)%alpha ) boxes(box_idx1, box_idx2)%alpha = alpha
                  IF( Kbc > boxes(box_idx1, box_idx2)%K ) boxes(box_idx1, box_idx2)%K = Kbc*p2s_b*p2s_c*&
                    pi**1.25_dp*SQRT(2.0_dp)
                END DO
              END DO
            END DO
            boxes(box_idx1, box_idx2)%center(1) = (px_max + px_min) / 2.0_dp
            boxes(box_idx1, box_idx2)%center(2) = (py_max + py_min) / 2.0_dp
            boxes(box_idx1, box_idx2)%center(3) = (pz_max + pz_min) / 2.0_dp
            R_max = TINY(R_max)
            DO i = 1,SIZE(neighbor_cells)
              im_C(1) = C11(1) + neighbor_cells(i)%cell(1) * box_size(1)
              im_C(2) = C11(2) + neighbor_cells(i)%cell(2) * box_size(2)
              im_C(3) = C11(3) + neighbor_cells(i)%cell(3) * box_size(3)
              rbc2 = (rb(1)-im_C(1))**2+(rb(2)-im_C(2))**2+(rb(3)-im_C(3))**2
              DO jpgf = 1,npgfb(jset)
                zet1 = zetb(jpgf,jset)
                DO kpgf = 1,npgfc(kset)
                  zet2 = zetc(kpgf,kset)
                  alpha = zet1+zet2
                  rp = (zet1*rb+zet2*im_C)/alpha
                  R_max = SQRT((rp(1)-boxes(box_idx1, box_idx2)%center(1))**2 + &
                               (rp(2)-boxes(box_idx1, box_idx2)%center(2))**2 + &
                               (rp(3)-boxes(box_idx1, box_idx2)%center(3))**2)
                  IF( R_max > boxes(box_idx1, box_idx2)%R ) boxes(box_idx1, box_idx2)%R = R_max
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
    DEALLOCATE(contraction_max, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle)
  END SUBROUTINE screen_farfield_per

!!***** four_center_integrals/update_pmax_matrix [1.0] *
!!
!!   NAME
!!     update_pmax_matrix
!!
!!   FUNCTION
!!     updates the maximum of the density matrix in compressed form for screening purposes 
!!
!!   NOTES
!!     - updates for each pair of shells the maximum absolute value of p
!!
!!   INPUTS
!!    -  pmax_set: buffer to store matrix
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!    09.2007 created [Manuel Guidon]
!!
!****************************************************************************
  
  SUBROUTINE update_pmax_matrix(pmax_set, &
              full_density,full_density_beta,particle_set,atomic_kind_set,error)

    REAL(kind=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: pmax_set
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: full_density, &
                                                full_density_beta
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_pmax_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, jatom, jset, katom, kset, max_set, mb, mc, natom, &
      nsetb, nsetc, nsgfb_total, nsgfc_total, pb, pb1, pb2, pc, pc1, pc2, &
      sgfb, sgfc, stat
    INTEGER, DIMENSION(:), POINTER           :: last_sgf_global, nsgfb, nsgfc
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfb, first_sgfc
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: pmax_tmp
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_b, orb_basis_c

    CALL timeset(routineN,handle)
    failure = .FALSE.

    natom=SIZE(particle_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, maxnset=max_set)

    ! create a mapping array of easy access to the density matrix elements
    ALLOCATE(last_sgf_global(0:natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    last_sgf_global(0)=0
    DO jatom=1,natom
      atom_kind => particle_set(jatom)%atomic_kind
 
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_b)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                             nsgf=nsgfb_total)
      last_sgf_global(jatom) = last_sgf_global(jatom-1)+nsgfb_total
    END DO

    pmax_set=0.0_dp

    DO jatom=1,natom

      atom_kind => particle_set(jatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_b)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                             nsgf=nsgfb_total,&
                             nset=nsetb,& 
                             first_sgf=first_sgfb,&
                             nsgf_set=nsgfb)
      pb1 = last_sgf_global(jatom) - nsgfb_total

      DO katom = 1,natom
        atom_kind => particle_set(katom)%atomic_kind
        CALL get_atomic_kind(atomic_kind=atom_kind,&
                             orb_basis_set=orb_basis_c)
        CALL get_gto_basis_set(gto_basis_set=orb_basis_c,&
                               nsgf=nsgfc_total,&
                               nset=nsetc,&
                               first_sgf=first_sgfc,&
                               nsgf_set=nsgfc)
        pc1 = last_sgf_global(katom) - nsgfc_total

        DO jset = 1,nsetb
          sgfb = first_sgfb(1,jset)
          DO kset = 1,nsetc
            sgfc = first_sgfc(1,kset)
            pmax_tmp = 0.0_dp
            DO mc = 1,nsgfc(kset)
               pc2 = mc+sgfc-1
               pc = pc1 + pc2
               DO mb = 1,nsgfb(jset)
                  pb2 = mb+sgfb-1
                  pb = pb1 + pb2
                  pmax_tmp=MAX(pmax_tmp,ABS(full_density(pb,pc)))
                  IF( ASSOCIATED(full_density_beta) ) pmax_tmp=MAX(pmax_tmp,ABS(full_density_beta(pb,pc)))
               ENDDO
            ENDDO
            pmax_set(jset,kset,jatom,katom) = pmax_tmp
          END DO
        END DO
      END DO
    END DO

    DEALLOCATE(last_sgf_global,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE update_pmax_matrix

END MODULE four_center_integrals
