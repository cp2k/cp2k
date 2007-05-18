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
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE libint_interface,                ONLY: evaluate_deriv_eri,&
                                             evaluate_eri,&
                                             evaluate_eri_screen,&
                                             initialize_libderiv,&
                                             initialize_libint,&
                                             lib_deriv,&
                                             lib_int
  USE mathconstants
  USE mathlib,                         ONLY: symmetrize_matrix
  USE message_passing,                 ONLY: mp_max,&
                                             mp_sum,&
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
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
   
  PUBLIC  integrate_four_center, derivatives_four_center, compute_memory_needs

  TYPE(lib_int), SAVE      :: lib
  TYPE(lib_deriv), SAVE    :: deriv
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'four_center_integrals'
  INTEGER:: file_counter_e = 0, file_counter_f = 100

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
!!     - screeb: precomputed screening values w.r.t. different zetas
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
                      primitive_integrals,screenij,screenkl, eps_schwarz)
  
    REAL(dp), INTENT(IN)                     :: ra(3), rb(3), rc(3), rd(3)
    REAL(dp), DIMENSION(:), INTENT(IN)       :: zeta, zetb, zetc, zetd
    INTEGER, INTENT(IN)                      :: la_min, la_max, lb_min, &
                                                lb_max, lc_min, lc_max, &
                                                ld_min, ld_max, npgfa, npgfb, &
                                                npgfc, npgfd
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: primitive_integrals
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: screenij, screenkl
    REAL(dp), INTENT(IN)                     :: eps_schwarz

    INTEGER :: ipgf, jpgf, kpgf, la, lb, lc, ld, lpgf, offset_a, offset_a1, &
      offset_b, offset_b1, offset_c, offset_c1, offset_d, offset_d1
    REAL(dp)                                 :: max_val1, max_val2

    DO ipgf = 1,npgfa
      offset_a1 = (ipgf-1)*ncoset(la_max)
      DO jpgf = 1,npgfb
        offset_b1 = (jpgf-1)*ncoset(lb_max)
        max_val1 = screenij(ipgf,jpgf)
        DO kpgf = 1,npgfc
          offset_c1 = (kpgf-1)*ncoset(lc_max)
          DO lpgf = 1,npgfd
            offset_d1 = (lpgf-1)*ncoset(ld_max)
            max_val2 = max_val1*screenkl(kpgf,lpgf)
            IF(max_val2<eps_schwarz) CYCLE
            DO la = la_min,la_max
              offset_a = offset_a1 + ncoset(la-1)
              DO lb = lb_min, lb_max
                offset_b = offset_b1 + ncoset(lb-1)
                DO lc = lc_min, lc_max
                  offset_c = offset_c1 + ncoset(lc-1)
                  DO ld = ld_min, ld_max
                    offset_d = offset_d1 + ncoset(ld-1)
                    !Build primitives
                    CALL evaluate_eri(lib, ra, rb, rc, rd,&
                                      zeta(ipgf), zetb(jpgf), zetc(kpgf), zetd(lpgf),&
                                      la, lb, lc ,ld,&
                                      primitive_integrals(offset_a+1:,&
                                                          offset_b+1:,&
                                                          offset_c+1:,&
                                                          offset_d+1:))
                  END DO !ld
                END DO !lc
              END DO !lb
            END DO !la

          END DO !lpgf
        END DO !kpgf
      END DO !jpgf
    END DO !ipgf
  END SUBROUTINE coulomb4
!***************************************************************************

!!***** four_center_integrals/screen4 [1.0] *
!!
!!   NAME
!!     screen4
!!
!!   FUNCTION
!!     calculates max values of two-electron integrals in a quartet/shell
!!     w.r.t. to different sets using the library lib_int
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

  SUBROUTINE screen4(ra, rb, rc, rd, zeta, zetb, zetc, zetd,&
                      la_min, la_max, lb_min, lb_max,&
                      lc_min, lc_max, ld_min, ld_max,&
                      npgfa, npgfb, npgfc, npgfd,&
                      max_val)
  
    REAL(dp), INTENT(IN)                     :: ra(3), rb(3), rc(3), rd(3)
    REAL(dp), DIMENSION(:), INTENT(IN)       :: zeta, zetb, zetc, zetd
    INTEGER, INTENT(IN)                      :: la_min, la_max, lb_min, &
                                                lb_max, lc_min, lc_max, &
                                                ld_min, ld_max, npgfa, npgfb, &
                                                npgfc, npgfd
    REAL(dp), INTENT(INOUT)                  :: max_val

    INTEGER                                  :: ipgf, jpgf, kpgf, la, lb, lc, &
                                                ld, lpgf
    REAL(dp)                                 :: max_val_temp

    max_val_temp = max_val
    DO ipgf = 1,npgfa
      DO jpgf = 1,npgfb
        DO kpgf = 1,npgfc
          DO lpgf = 1,npgfd

            DO la = la_min,la_max
              DO lb = lb_min, lb_max
                DO lc = lc_min, lc_max
                  DO ld = ld_min, ld_max
                    !Build primitives
                    CALL evaluate_eri_screen(lib, ra, rb, rc, rd,&
                                             zeta(ipgf), zetb(jpgf), zetc(kpgf), zetd(lpgf),&
                                             la, lb, lc ,ld,&
                                             max_val_temp)
                    max_val = MAX(max_val, max_val_temp)
                  END DO !ld
                END DO !lc
              END DO !lb
            END DO !la

          END DO !lpgf
        END DO !kpgf
      END DO !jpgf
    END DO !ipgf
  
  END SUBROUTINE screen4
!***************************************************************************

!!***** four_center_integrals/screen41 [1.0] *
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

  SUBROUTINE screen41(ra, rb, rc, rd, zeta, zetb, zetc, zetd,&
                      la_min, la_max, lb_min, lb_max,&
                      lc_min, lc_max, ld_min, ld_max,&
                      max_val)

    REAL(dp), INTENT(IN)                     :: ra(3), rb(3), rc(3), rd(3), &
                                                zeta, zetb, zetc, zetd
    INTEGER, INTENT(IN)                      :: la_min, la_max, lb_min, &
                                                lb_max, lc_min, lc_max, &
                                                ld_min, ld_max
    REAL(dp), INTENT(INOUT)                  :: max_val

    INTEGER                                  :: la, lb, lc, ld
    REAL(dp)                                 :: max_val_temp

    max_val_temp = max_val

    DO la = la_min,la_max
      DO lb = lb_min, lb_max
        DO lc = lc_min, lc_max
          DO ld = ld_min, ld_max
            !Build primitives
            CALL evaluate_eri_screen(lib, ra, rb, rc, rd,&
                                     zeta, zetb, zetc, zetd,&
                                     la, lb, lc ,ld,&
                                     max_val_temp)
            max_val = MAX(max_val, max_val_temp)
          END DO !ld
        END DO !lc
      END DO !lb
    END DO !la

  END SUBROUTINE screen41

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
  SUBROUTINE compute_screening_matrices(screen_pgf,screen_set,screen_atom,particle_set,atomic_kind_set,error)

    REAL(kind=dp), &
      DIMENSION(:, :, :, :, :, :), POINTER   :: screen_pgf
    REAL(kind=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: screen_set
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: screen_atom
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_screening_matrices', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iatom, iset, jatom, jpgf, jset, katom, kpgf, kset, &
      max_pgf, max_set, natom, nseta, nsetb, nsetc, stat, ncob,ncoc,sgfb,sgfc,k
    INTEGER, DIMENSION(:), POINTER           :: lb_max, lb_min, lc_max, &
                                                lc_min, npgfa, npgfb, npgfc,nsgfb,nsgfc
    LOGICAL                                  :: failure
    REAL(dp), DIMENSION(:, :), POINTER       :: zetb, zetc, sphi_b, sphi_c
    REAL(KIND=dp)                            :: max_val, max_val2, rb(3), &
                                                rc(3),p2s_b,p2s_c
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_a, orb_basis_b, &
                                                orb_basis_c
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfb, first_sgfc


    CALL timeset(routineN,handle)

    failure = .FALSE.

    natom=SIZE(particle_set)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, maxnset=max_set)

    max_pgf = 0
    DO iatom=1,natom

      atom_kind => particle_set(iatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_a)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_a,&
                             npgf=npgfa,&
                             nset=nseta)
      DO iset=1,nseta
        max_pgf = MAX(max_pgf,npgfa(iset))
      END DO
    END DO

    ! write(6,*) "screen mem",natom,max_set,max_pgf,natom*natom*max_set*max_set*max_pgf*max_pgf

    ALLOCATE(screen_pgf(natom,natom,max_set,max_set,max_pgf, max_pgf),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(screen_set(natom,natom,max_set,max_set),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(screen_atom(natom,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    screen_pgf=0.0_dp
    screen_set=0.0_dp
    screen_atom=0.0_dp

    DO jatom=1,natom

      atom_kind => particle_set(jatom)%atomic_kind
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           orb_basis_set=orb_basis_b)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_b,&
                             lmax=lb_max,&
                             lmin=lb_min,&
                             npgf=npgfb,&
                             nset=nsetb,& 
                             first_sgf=first_sgfb,&
                             sphi=sphi_b,&
                             nsgf_set=nsgfb,&
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
                               sphi=sphi_c,&
                               first_sgf=first_sgfc,&
                               nsgf_set=nsgfc,&
                               zet=zetc)
        rc=particle_set(katom)%r(:)
        DO jset = 1,nsetb
          ! takes the primitive to contracted transformation into account
          ncob = npgfb(jset)*ncoset(lb_max(jset))
          sgfb = first_sgfb(1,jset)
          ! if the primitives are assumed to be all of max_val2, max_val2*p2s_b becomes
          ! the maximum value after multiplication with sphi_b
          p2s_b = MAXVAL((/(SUM(ABS(sphi_b(1:ncob,k))),k=sgfb,sgfb+nsgfb(jset)-1)/))

          DO kset = 1,nsetc
            ! takes the primitive to contracted transformation into account
            ncoc = npgfc(kset)*ncoset(lc_max(kset))
            sgfc = first_sgfc(1,kset)
            p2s_c = MAXVAL((/(SUM(ABS(sphi_c(1:ncoc,k))),k=sgfc,sgfc+nsgfc(kset)-1)/))

            DO jpgf = 1,npgfb(jset)
              DO kpgf = 1,npgfc(kset)
                max_val2=0.0_dp
                CALL screen41(rb, rc, rb, rc,&
                              zetb(jpgf,jset), zetc(kpgf,kset), zetb(jpgf,jset), zetc(kpgf,kset),&
                              lb_min(jset), lb_max(jset), lc_min(kset), lc_max(kset),&
                              lb_min(jset), lb_max(jset), lc_min(kset), lc_max(kset),&
                              max_val2)
                screen_pgf(jatom,katom,jset,kset,jpgf,kpgf) = SQRT(max_val2)*p2s_b*p2s_c
              END DO
            END DO
            screen_set(jatom,katom,jset,kset) = MAXVAL(screen_pgf(jatom,katom,jset,kset,1:npgfb(jset),1:npgfc(kset)))
          END DO
        END DO
        screen_atom(jatom,katom)=MAXVAL(screen_set(jatom,katom,1:nsetb,1:nsetc))
      END DO
    END DO

    CALL timestop(handle)

  END SUBROUTINE compute_screening_matrices
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
  SUBROUTINE compute_memory_needs(memory_blocks_ab,eps_schwarz,particle_set,atomic_kind_set,para_env,error)
    INTEGER, DIMENSION(:,:), INTENT(OUT) :: memory_blocks_ab
    REAL(KIND=dp), INTENT(IN)            :: eps_schwarz
    TYPE(particle_type), DIMENSION(:), &
       POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
       POINTER                                :: atomic_kind_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_memory_needs', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp), DIMENSION(:,:,:,:,:,:), POINTER :: screening2
    REAL(KIND=dp), DIMENSION(:,:,:,:), POINTER :: schwarz_screening
    REAL(KIND=dp), DIMENSION(:,:), POINTER :: screen_atom
    REAL(KIND=dp)                          :: max_val1,max_val2
    TYPE(gto_basis_set_type), POINTER        :: orb_basis
    INTEGER :: stat,handle,iatom,jatom,katom,latom,natom,iset,jset,kset,lset
    INTEGER :: nseta,nsetb,nsetc,nsetd,iblock, max_am
    INTEGER, DIMENSION(:), POINTER :: nsgfa, nsgfb, nsgfc, nsgfd
    LOGICAL :: failure, include_ab

    CALL timeset(routineN,handle)
    failure=.FALSE.

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, maxlgto=max_am)
    CALL initialize_libint(lib, max_am, error)

    memory_blocks_ab = 0

    ! set up screening matrices
    CALL compute_screening_matrices(screening2,schwarz_screening,screen_atom,particle_set,atomic_kind_set,error)

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

            ! take the symmetric storage form into account
            IF (iatom > latom) THEN
              include_ab = (MODULO(iatom + latom,2) /= 0)
            ELSE
              include_ab = (MODULO(iatom + latom,2) == 0)
            END IF
            IF (.NOT. include_ab) CYCLE

            CALL get_atomic_kind(atomic_kind=particle_set(latom)%atomic_kind, orb_basis_set=orb_basis)
            CALL get_gto_basis_set(gto_basis_set=orb_basis, nset=nsetd, nsgf_set=nsgfd)

            DO iset = 1, nseta
              DO jset = 1,nsetb
                max_val1 = schwarz_screening(iatom,jatom,iset,jset)
                DO kset = 1,nsetc
                  DO lset = 1,nsetd
                    max_val2 = max_val1 * schwarz_screening(katom,latom,kset,lset)
                    IF(max_val2<eps_schwarz) CYCLE
                    memory_blocks_ab(iatom,latom)=memory_blocks_ab(iatom,latom) + &
                                                  nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DEALLOCATE(schwarz_screening,screening2,screen_atom,STAT=stat)
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
                                   geometry_did_change, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: hf_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL                                  :: geometry_did_change
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'integrate_four_center', &
      routineP = moduleN//':'//routineN

    INTEGER :: cgfa_first, cgfa_last, cgfb_first, cgfb_last, cgfc_first, &
      cgfc_last, cgfd_first, cgfd_last, current_counter, handle, handle_getP, &
      handle_main, iatom, integral_counter, ipgf, iset, iw, jatom, jpgf, &
      jset, katom, kpgf, kset, la, latom, lb, lc, ld, lpgf, lset, ma, max_am, &
      max_set, maxsgf, mb, mc, md, na, natom, nb, nc, ncoa, ncob, ncoc, ncod, &
      nd, nseta, nsetb, nsetc, nsetd, nsgfa_total, nsgfb_total, nsgfc_total, &
      nsgfd_total, nspins, number_of_precalc_integrals, offset_a1, offset_b1, &
      offset_c1, offset_d1, pa, pa1, pa2, pb, pb1, pb2, pc, pc1, pc2, pd, &
      pd1, pd2, pna, pnb, pnc, pnd, sgfa, sgfb, sgfc, sgfd
    INTEGER :: stat, total_integrals_to_store, total_integrals_to_store_aver, &
      total_integrals_to_store_max
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: last_sgf_global
    INTEGER, DIMENSION(:), POINTER :: la_max, la_min, lb_max, lb_min, lc_max, &
      lc_min, ld_max, ld_min, npgfa, npgfb, npgfc, npgfd, nsgfa, nsgfb, &
      nsgfc, nsgfd
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb, &
                                                first_sgfc, first_sgfd
    LOGICAL                                  :: failure
    REAL(dp) :: ene_x_aa, ene_x_bb, eps_schwarz, fac, hf_fraction, max_val1, &
      max_val2, max_val3, max_val4, memory_usage, prim, ra(3), rb(3), rc(3), &
      rd(3), sa, sb, sc, sd, tmp1, tmp11, tmp2, tmp21, tmp3, tmp31, max_int_val
    REAL(dp), DIMENSION(:, :), POINTER       :: screen_atom, sphi_a, sphi_b, &
                                                sphi_c, sphi_d, zeta, zetb, &
                                                zetc, zetd
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: primitive_integrals, &
                                                schwarz_screening, T, T1
    REAL(dp), DIMENSION(:, :, :, :, :, :), &
      POINTER                                :: screening2
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: full_density, &
                                                full_density_beta, &
                                                sparse_block, &
                                                sparse_block_beta
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_a, orb_basis_b, &
                                                orb_basis_c, orb_basis_d
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(real_block_node_type), POINTER      :: block_node, block_node_beta

! many timings somewhat more clear

    CALL mp_sync(para_env%group)
    CALL timeset(routineN,"I",'',handle)

    CALL section_vals_val_get(hf_section,"FRACTION",r_val=hf_fraction,error=error)
    CALL section_vals_val_get(hf_section,"EPS_SCHWARZ",r_val=eps_schwarz,error=error)
   
    number_of_precalc_integrals = qs_env%x_data%number_of_elements
    integral_counter = 1
    memory_usage = 0.0_dp
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

    CALL compute_screening_matrices(screening2,schwarz_screening,screen_atom,particle_set,atomic_kind_set,error)
    
    CALL timeset(routineN//"_main",handle_main)

    total_integrals_to_store=0
    memory_usage = 0.0_dp
    ene_x_aa = 0.0_dp
    ene_x_bb = 0.0_dp
    IF(nspins == 1) THEN
      fac = 0.5_dp * hf_fraction
    ELSE
      fac = 1.0_dp * hf_fraction
    END IF

    ALLOCATE(T(ncoset(max_am),&
               ncoset(max_am),&
               ncoset(max_am),&
               ncoset(max_am)),STAT=stat)
    tmp11 = 0.0_dp
    DO iatom=1,ks_matrix(1)%matrix%nblock_row

      block_node => first_block_node(ks_matrix(1)%matrix,iatom)
      IF(nspins == 2) block_node_beta => first_block_node(ks_matrix(2)%matrix,iatom)
 
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
                             first_sgf=first_sgfa,&
                             sphi=sphi_a,&
                             nsgf=nsgfa_total)
      ra=particle_set(iatom)%r(:)
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
                               first_sgf=first_sgfb,&
                               sphi=sphi_b,&
                               nsgf=nsgfb_total)
        rb=particle_set(jatom)%r(:)
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
                                 first_sgf=first_sgfc,&
                                 sphi=sphi_c,&
                                 nsgf=nsgfc_total)
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
                                   first_sgf=first_sgfd,&
                                   sphi=sphi_d,&
                                   nsgf=nsgfd_total)
 
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
                max_val1 = schwarz_screening(iatom,jatom,iset,jset)
                DO kset = 1,nsetc
                  ncoc = npgfc(kset)*ncoset(lc_max(kset))
                  sgfc = first_sgfc(1,kset)
                  DO lset = 1,nsetd
                    ncod = npgfd(lset)*ncoset(ld_max(lset)) 
                    sgfd = first_sgfd(1,lset)
                    max_val2 = max_val1 * schwarz_screening(katom,latom,kset,lset)
                    IF(max_val2<eps_schwarz) CYCLE
                    
                    current_counter = nsgfa(iset)*nsgfb(jset)*nsgfc(kset)*nsgfd(lset)
                    total_integrals_to_store=total_integrals_to_store+current_counter

                    IF( integral_counter + current_counter<= number_of_precalc_integrals) THEN
                      IF(.NOT.geometry_did_change) THEN
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
                      ELSE
                        ALLOCATE(primitive_integrals(ncoa,ncob,ncoc,ncod),STAT=stat)
                        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                        ALLOCATE(T1(ncoa,ncob,ncoc,ncod),STAT=stat)
                        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

                        !Get primitive array into cache
                        primitive_integrals = 0.0_dp
                        CALL coulomb4(ra, rb, rc, rd,&
                                      zeta(:,iset), zetb(:,jset), zetc(:,kset), zetd(:,lset),&
                                      la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                                      lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                                      npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset),&
                                      primitive_integrals,screening2(iatom,jatom,iset,jset,:,:),&
                                      screening2(katom,latom,kset,lset,:,:),eps_schwarz)

                        CALL dgemm("T","N",ncob*ncoc*ncod,nsgfa(iset),ncoa,&
                                   1.0_dp, primitive_integrals(1,1,1,1),ncoa,&
                                   sphi_a(1,sgfa), SIZE(sphi_a,1),&
                                   0.0_dp, T1(1,1,1,1),ncob*ncoc*ncod)

                        CALL dgemm("T","N",ncoa*ncoc*ncod,nsgfb(jset),ncob,&
                                   1.0_dp, T1(1,1,1,1),ncob,&
                                   sphi_b(1,sgfb), SIZE(sphi_b,1),&
                                   0.0_dp, primitive_integrals(1,1,1,1),ncoa*ncoc*ncod)

                        CALL dgemm("T","N",ncoa*ncob*ncod,nsgfc(kset),ncoc,&
                                   1.0_dp, primitive_integrals(1,1,1,1),ncoc,&
                                   sphi_c(1,sgfc), SIZE(sphi_c,1),&
                                   0.0_dp, T1(1,1,1,1),ncoa*ncob*ncod)

                        CALL dgemm("T","N",ncoa*ncob*ncoc,nsgfd(lset),ncod,&
                                   1.0_dp, T1(1,1,1,1),ncod,&
                                   sphi_d(1,sgfd), SIZE(sphi_d,1),&
                                   0.0_dp, primitive_integrals(1,1,1,1),ncoa*ncob*ncoc)

                        ! max_int_val=MAXVAL(ABS(primitive_integrals(1:nsgfa(iset),1:nsgfb(jset),1:nsgfc(kset),1:nsgfd(lset))))
                        ! max_int_val should be smaller equal max_val2
                        ! write(2000+para_env%mepos,*) max_int_val,max_val2,max_int_val/max_val2

                        DEALLOCATE(T1,STAT=stat)
                        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

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
                                tmp2 = tmp1*primitive_integrals(ma,mb,mc,md)
                                sparse_block(pa2,pd2) = sparse_block(pa2,pd2) - tmp2
                                tmp2 = tmp3*full_density(pa,pd)*tmp2
                                ene_x_aa = ene_x_aa - tmp2
                                IF( nspins == 2 ) THEN
                                  tmp21 = tmp11*primitive_integrals(ma,mb,mc,md)
                                  sparse_block_beta(pa2,pd2) = sparse_block_beta(pa2,pd2) - tmp21
                                  tmp21 = tmp31*full_density_beta(pa,pd)*tmp21
                                  ene_x_bb = ene_x_bb -tmp21
                                END IF
                                qs_env%x_data%x_buffer(integral_counter) = &
                                   primitive_integrals(ma,&
                                                       mb,&
                                                       mc,&
                                                       md)
                                integral_counter = integral_counter + 1
                              END DO
                            END DO
                          END DO
                        END DO

                        DEALLOCATE(primitive_integrals,STAT=stat)
                        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                      ENDIF
                    ELSE
                      DO lpgf = 1,npgfd(lset)
                        offset_d1 = (lpgf-1)*ncoset(ld_max(lset))
                        cgfd_first = offset_d1 + ncoset(ld_min(lset)-1)+1
                        cgfd_last = offset_d1 + ncoset(ld_max(lset))
                        DO kpgf = 1,npgfc(kset)
                         offset_c1 = (kpgf-1)*ncoset(lc_max(kset))
                         cgfc_first = offset_c1 + ncoset(lc_min(kset)-1)+1
                         cgfc_last = offset_c1 + ncoset(lc_max(kset))
                         max_val3 = screening2(katom,latom,kset,lset,kpgf,lpgf)
                          DO jpgf = 1,npgfb(jset)
                            offset_b1 = (jpgf-1)*ncoset(lb_max(jset))
                            cgfb_first = offset_b1 + ncoset(lb_min(jset)-1)+1
                            cgfb_last = offset_b1 + ncoset(lb_max(jset))
                            DO ipgf = 1,npgfa(iset)
                              max_val4 = max_val3* screening2(iatom,jatom,iset,jset,ipgf,jpgf)
                              IF(max_val4<eps_schwarz) CYCLE
                              offset_a1 = (ipgf-1)*ncoset(la_max(iset))
                              cgfa_first = offset_a1 + ncoset(la_min(iset)-1)+1
                              cgfa_last = offset_a1+ncoset(la_max(iset))

                              DO ld = ld_min(lset),ld_max(lset)
                                DO lc = lc_min(kset), lc_max(kset)
                                  DO lb = lb_min(jset), lb_max(jset)
                                    DO la = la_min(iset), la_max(iset)
                                      !Build primitives
                                      CALL evaluate_eri(lib, ra, rb, rc, rd,&
                                                        zeta(ipgf,iset), zetb(jpgf,jset),&
                                                        zetc(kpgf,kset), zetd(lpgf,lset),&
                                                        la, lb, lc ,ld,&
                                                        T(ncoset(la-1)+1:,&
                                                          ncoset(lb-1)+1:,&
                                                          ncoset(lc-1)+1:,&
                                                          ncoset(ld-1)+1:))
                                    END DO !ld
                                  END DO !lc
                                END DO !lb
                              END DO !la
                              
                              DO nd = cgfd_first,cgfd_last
                                pnd = nd-offset_d1
                                DO nc = cgfc_first,cgfc_last
                                  pnc = nc-offset_c1
                                  DO nb = cgfb_first,cgfb_last 
                                    pnb = nb-offset_b1
                                    DO na = cgfa_first,cgfa_last
                                      pna = na - offset_a1
                                      prim = T(pna,pnb,pnc,pnd)
                                      IF(ABS(prim)<eps_schwarz ) CYCLE
                                      DO mc = sgfc,nsgfc(kset)+sgfc-1
                                        sc = prim*sphi_c(nc,mc)
                                        pc2=mc+pc1
                                        IF(ABS(sc) < eps_schwarz) CYCLE
                                        DO mb = sgfb,nsgfb(jset)+sgfb-1
                                          sb = sc * sphi_b(nb,mb)
                                          pb2=mb+pb1
                                          tmp1 = fac*full_density(pb2,pc2)
                                          IF( nspins == 2) tmp11 = fac*full_density_beta(pb2,pc2)
                                          IF(ABS(tmp1)+ABS(tmp11)<eps_schwarz) CYCLE
                                          IF(ABS(sb) <eps_schwarz) CYCLE
                                          DO md = sgfd,nsgfd(lset)+sgfd-1
                                            sd = sb * sphi_d(nd,md)
                                            pd2=md+pd1
                                            IF(ABS(sd) <eps_schwarz) CYCLE
                                            DO ma = sgfa,nsgfa(iset)+sgfa-1
                                              sa = sd * sphi_a(na,ma)
                                              pa2=ma+pa1
                                              IF(ABS(sa) <eps_schwarz) CYCLE
                                              sparse_block(ma,md) = sparse_block(ma,md) - &
                                                                    tmp1*sa
                                              ene_x_aa = ene_x_aa - tmp3*tmp1*full_density(pa2,pd2)*sa
                                              IF(nspins == 2 ) THEN
                                                sparse_block_beta(ma,md) = sparse_block_beta(ma,md) - &
                                                                             tmp11*sa
                                                ene_x_aa = ene_x_aa - tmp3*tmp11*full_density_beta(pa2,pd2)*sa
                                              END IF
                                            END DO !ma
                                          END DO !mb
                                        END DO !mc
                                      END DO !md

                                    END DO !na
                                  END DO !nb
                                END DO !nc
                              END DO !nd

                            END DO !lpgf
                          END DO !kpgf
                        END DO !jpgf
                      END DO !ipgf
                    END IF

                  END DO !lset
                END DO !kset
              END DO !jset
            END DO !set

          END DO !latom
        END DO !katom
      END DO !jatom
    END DO !iatom


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

    DEALLOCATE(last_sgf_global,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(schwarz_screening,screening2,screen_atom,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(full_density,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(T,STAT=stat)
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
                    primitive_forces, work_forces, screenij, screenkl, eps_schwarz )
  
    REAL(dp), INTENT(IN)                     :: ra(3), rb(3), rc(3), rd(3)
    REAL(dp), DIMENSION(:), INTENT(IN)       :: zeta, zetb, zetc, zetd
    INTEGER, INTENT(IN)                      :: la_min, la_max, lb_min, &
                                                lb_max, lc_min, lc_max, &
                                                ld_min, ld_max, npgfa, npgfb, &
                                                npgfc, npgfd
    REAL(dp), DIMENSION(:, :, :, :, :), &
      POINTER                                :: primitive_forces
    REAL(dp), DIMENSION(:, :), POINTER       :: work_forces
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: screenij, screenkl
    REAL(dp), INTENT(IN)                     :: eps_schwarz

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'forces4'

    INTEGER :: ipgf, jpgf, kpgf, la, lb, lc, ld, lpgf, offset_a, offset_a1, &
      offset_b, offset_b1, offset_c, offset_c1, offset_d, offset_d1
    REAL(dp)                                 :: max_val1, max_val2

    DO ipgf = 1,npgfa
      offset_a1 = (ipgf-1)*ncoset(la_max)
      DO jpgf = 1,npgfb
        offset_b1 = (jpgf-1)*ncoset(lb_max)
        max_val1 = screenij(ipgf,jpgf)
        DO kpgf = 1,npgfc
          offset_c1 = (kpgf-1)*ncoset(lc_max)
          DO lpgf = 1,npgfd
            offset_d1 = (lpgf-1)*ncoset(ld_max)
            max_val2 = max_val1*screenkl(kpgf,lpgf)
            IF(max_val2<eps_schwarz) CYCLE
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
                                           work_forces, primitive_forces,&
                                           offset_a, offset_b, offset_c, offset_d)

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
                                     error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: hf_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'derivatives_four_center', &
      routineP = moduleN//':'//routineN

    INTEGER :: coord, handle, i_atom, iatom, ikind, iset, j_atom, jatom, &
      jkind, jset, k_atom, katom, kkind, kset, l_atom, latom, lkind, lset, &
      ma, max_am, max_set, maxsgf, mb, mc, md, natom, ncoa, ncob, ncoc, ncod, &
      nseta, nsetb, nsetc, nsetd, nsgfa_total, nsgfb_total, nsgfc_total, &
      nsgfd_total, nspins, pa, pa1, pa2, pb, pb1, pb2, pc, pc1, pc2, pd, pd1, &
      pd2, sgfa, sgfb, sgfc, sgfd, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of, &
                                                last_sgf_global
    INTEGER, DIMENSION(:), POINTER :: la_max, la_min, lb_max, lb_min, lc_max, &
      lc_min, ld_max, ld_min, ncgfa, ncgfb, ncgfc, ncgfd, npgfa, npgfb, &
      npgfc, npgfd, nsgfa, nsgfb, nsgfc, nsgfd
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb, &
                                                first_sgfc, first_sgfd
    LOGICAL                                  :: failure
    REAL(dp) :: eps_schwarz, fac, hf_fraction, max_val1, max_val2, ra(3), &
      rb(3), rc(3), rd(3), temp1, temp1b, temp2
    REAL(dp), DIMENSION(:, :), POINTER       :: sphi_a, sphi_b, sphi_c, &
                                                sphi_d, work_forces, zeta, &
                                                zetb, zetc, zetd
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: schwarz_screening, T1
    REAL(dp), DIMENSION(:, :, :, :, :), &
      POINTER                                :: primitive_forces
    REAL(dp), DIMENSION(:, :, :, :, :, :), &
      POINTER                                :: screening2
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: full_density, &
                                                full_density_beta, &
                                                screen_atom, sparse_block
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

    CALL timeset(routineN,"I",'',handle)

    CALL section_vals_val_get(hf_section,"FRACTION",r_val=hf_fraction,error=error)
    CALL section_vals_val_get(hf_section,"EPS_SCHWARZ",r_val=eps_schwarz,error=error)

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
    END IF

    CALL compute_screening_matrices(screening2,schwarz_screening,screen_atom,particle_set,atomic_kind_set,error)
  
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
                max_val1 = schwarz_screening(iatom,jatom,iset,jset)
                DO kset = 1,nsetc
                  ncoc = npgfc(kset)*ncoset(lc_max(kset))
                  sgfc = first_sgfc(1,kset)
                  DO lset = 1,nsetd
                    ncod = npgfd(lset)*ncoset(ld_max(lset)) 
                    sgfd = first_sgfd(1,lset)

                    max_val2 = max_val1 * schwarz_screening(katom,latom,kset,lset)
                    IF(max_val2<eps_schwarz) CYCLE
                    
                    ALLOCATE(primitive_forces(ncoa,ncob,ncoc,ncod,12),STAT=stat)
                    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                    primitive_forces = 0.0_dp
                    
                    ALLOCATE(T1(ncoa,ncob,ncoc,ncod),STAT=stat)
                    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                    
                    CALL forces4(ra, rb, rc, rd,&
                                 zeta(:,iset), zetb(:,jset), zetc(:,kset), zetd(:,lset),&
                                 la_min(iset), la_max(iset), lb_min(jset), lb_max(jset),&
                                 lc_min(kset), lc_max(kset), ld_min(lset), ld_max(lset),&
                                 npgfa(iset), npgfb(jset), npgfc(kset), npgfd(lset),&
                                 primitive_forces,&
                                 work_forces, screening2(iatom,jatom,iset,jset,:,:),&
                                 screening2(katom,latom,kset,lset,:,:), eps_schwarz)
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

    DEALLOCATE(last_sgf_global,schwarz_screening, screening2, screen_atom, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(full_density,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF(nspins==2) THEN
      DEALLOCATE(full_density_beta,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE(work_forces, atom_of_kind,kind_of,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    !
    ! this wraps to free_libint, but currently causes a segfault
    ! as a result, we don't call it, but some memory remains allocated
    !  
    ! CALL terminate_libderiv(deriv)

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE derivatives_four_center
!***************************************************************************

END MODULE four_center_integrals
