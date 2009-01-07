!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Several screening methods used in HFX calcualtions 
!> \par History
!>      04.2008 created [Manuel Guidon]
!> \author Manuel Guidon 
! *****************************************************************************
MODULE hfx_screening_methods 
  USE cell_types,                      ONLY: cell_type,&
                                             pbc,&
                                             real_to_scaled,&
                                             scaled_to_real
  USE erf_fn,                          ONLY: erfc
  USE f77_blas
  USE hfx_libint_interface,            ONLY: evaluate_eri_screen
  USE hfx_libint_wrapper,              ONLY: lib_int
  USE hfx_types,                       ONLY: hfx_basis_type,&
                                             hfx_box_scheme,&
                                             hfx_cell_type,&
                                             hfx_potential_type
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE mathconstants
  USE orbital_pointers
  USE particle_types,                  ONLY: particle_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
   
  PUBLIC screen4, comp_screen_mat, comp_pmax_mat, comp_screen_mat_per, &
         screen4_per, screen_farfield, screen_farfield_per, update_pmax_mat
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_screening_methods'
  REAL(dp), PARAMETER      :: epsilon = 0.00001_dp
  REAL(dp), PARAMETER      :: erfc_epsilon = 0.99998871620832942100_dp !erfc(epsilon)

!***
  
  CONTAINS

! *****************************************************************************
!> \brief calculates max values of two-electron integrals in a quartet/shell
!>      w.r.t. different zetas using the library lib_int
!> \param r position
!> \param z zeta
!> \param l angular momentum
!> \param npgf number of primitive cartesian gaussian in actual shell
!> \param max_val schwarz screening value
!> \param potential_parameter contains info for libint
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>     03.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE screen4(lib, ra, rb, zeta, zetb, &
                     la_min, la_max, lb_min, lb_max,&
                     npgfa, npgfb, &
                     max_val, potential_parameter, i_thread)

    TYPE(lib_int)                            :: lib
    REAL(dp), INTENT(IN)                     :: ra(3), rb(3)
    REAL(dp), DIMENSION(:), INTENT(IN)       :: zeta, zetb
    INTEGER, INTENT(IN)                      :: la_min, la_max, lb_min, &
                                                lb_max, npgfa, npgfb
    REAL(dp), INTENT(INOUT)                  :: max_val
    TYPE(hfx_potential_type)                 :: potential_parameter
    INTEGER, INTENT(IN)                      :: i_thread

    CHARACTER(len=*), PARAMETER :: routineN = 'screen4', &
      routineP = moduleN//':'//routineN

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
                                     max_val_temp, potential_parameter)
            max_val = MAX(max_val, max_val_temp)
          END DO !lb
        END DO !la
      END DO !jpgf
    END DO !ipgf

  END SUBROUTINE screen4

! *****************************************************************************
!> \brief calculates the screening matrices 
!> \note
!>      - already computes the square root
!>      - also takes the primitive to sgf transformation into account
!>      - in this way all stored contracted integrals are smaller or equal to their estimate
!> \par History
!>     05.2007 refactored [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE comp_screen_mat(lib, screen_set,screen_atom,screen_max,particle_set,&
                             max_set,potential_parameter,read_from_memory,&
                             kind_of, basis_parameter, i_thread, buffer, error)
    TYPE(lib_int)                            :: lib
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: screen_set
    REAL(dp), DIMENSION(:, :), POINTER       :: screen_atom
    REAL(dp)                                 :: screen_max
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: max_set
    TYPE(hfx_potential_type)                 :: potential_parameter
    LOGICAL                                  :: read_from_memory
    INTEGER                                  :: kind_of(*)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    INTEGER, INTENT(IN)                      :: i_thread
    REAL(dp), DIMENSION(:, :, :, :), &
      OPTIONAL, POINTER                      :: buffer
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'comp_screen_mat', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: jatom, jkind, jset, k, katom, &
                                                kkind, kset, natom, ncob, &
                                                nsetb, nsetc, sgfb, stat
    INTEGER, DIMENSION(:), POINTER           :: lb_max, lb_min, lc_max, &
                                                lc_min, npgfb, npgfc, nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfb
    LOGICAL                                  :: failure
    REAL(dp)                                 :: max_val, p2s_b, p2s_c, rb(3), &
                                                rc(3)
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: contraction_max
    REAL(dp), DIMENSION(:, :), POINTER       :: sphi_b, zetb, zetc

    failure = .FALSE.

    natom=SIZE(particle_set)
    ALLOCATE(screen_set(max_set,max_set,natom,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(screen_atom(natom,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(contraction_max(max_set,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    screen_set=0.0_dp
    screen_atom=0.0_dp
    contraction_max=0.0_dp

    DO jatom=1,natom
      jkind = kind_of(jatom)
      nsetb = basis_parameter(jkind)%nset
      npgfb => basis_parameter(jkind)%npgf
      lb_max => basis_parameter(jkind)%lmax
      first_sgfb => basis_parameter(jkind)%first_sgf
      sphi_b => basis_parameter(jkind)%sphi
      nsgfb => basis_parameter(jkind)%nsgf
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
      jkind = kind_of(jatom)
      lb_max => basis_parameter(jkind)%lmax
      lb_min => basis_parameter(jkind)%lmin
      npgfb => basis_parameter(jkind)%npgf
      nsetb = basis_parameter(jkind)%nset
      zetb => basis_parameter(jkind)%zet
      rb=particle_set(jatom)%r(:)
      DO katom = jatom,natom
        kkind = kind_of(katom)
        lc_max => basis_parameter(kkind)%lmax
        lc_min => basis_parameter(kkind)%lmin
        npgfc => basis_parameter(kkind)%npgf
        nsetc = basis_parameter(kkind)%nset
        zetc => basis_parameter(kkind)%zet
        rc=particle_set(katom)%r(:)
        DO jset = 1,nsetb
          p2s_b = contraction_max(jset,jatom)
          DO kset = 1,nsetc
            p2s_c = contraction_max(kset,katom)
            max_val = 0.0_dp
            IF( read_from_memory ) THEN
              screen_set(jset,kset,jatom,katom) = buffer(jset,kset,jatom,katom)
            ELSE
              CALL screen4(lib, rb, rc, &
                           zetb(:,jset), zetc(:,kset), &
                           lb_min(jset), lb_max(jset), lc_min(kset), lc_max(kset),&
                           npgfb(jset), npgfc(kset), &
                           max_val, potential_parameter, i_thread)
              screen_set(jset,kset,jatom,katom) = SQRT(max_val)*p2s_b*p2s_c
              IF( PRESENT(buffer) ) buffer(jset,kset,jatom,katom) = screen_set(jset,kset,jatom,katom)
            END IF
          END DO
        END DO
        screen_atom(jatom,katom)=MAXVAL(screen_set(1:nsetb,1:nsetc,jatom,katom))
      END DO
    END DO

    !! Symmetrize it
    DO jatom=1,natom
      DO katom=jatom,natom
        screen_atom(katom,jatom) = screen_atom(jatom,katom)
        DO jset=1,max_set
          DO kset=1,max_set
            screen_set(kset,jset,katom,jatom) = screen_set(jset,kset,jatom,katom)
          END DO
        END DO
      END DO
    END DO

    screen_max=MAXVAL(screen_atom)
    DEALLOCATE(contraction_max, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE comp_screen_mat

! *****************************************************************************
!> \brief calculates the maximum of the density matrix in compressed form for screening purposes 
!> \param screen_pmat_forces if false just returns a matrix of 1.0_dp
!> \note
!>      - computes for each pair of shells the maximum absolute value of p
!> \par History
!>     05.2007 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE comp_pmax_mat(pmax_set,pmax_atom,pmax_max, &
                           full_density,full_density_beta,screen_pmat_forces,particle_set,max_set,&
                           kind_of, basis_parameter, error)

    REAL(dp), DIMENSION(:, :, :, :), POINTER :: pmax_set
    REAL(dp), DIMENSION(:, :), POINTER       :: pmax_atom
    REAL(dp)                                 :: pmax_max
    REAL(dp), DIMENSION(:, :), POINTER       :: full_density, &
                                                full_density_beta
    LOGICAL, INTENT(IN)                      :: screen_pmat_forces
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: max_set
    INTEGER                                  :: kind_of(*)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'comp_pmax_mat', &
      routineP = moduleN//':'//routineN

    INTEGER :: jatom, jkind, jset, katom, kkind, kset, mb, mc, natom, nsetb, &
      nsetc, nsgfb_total, nsgfc_total, pb, pb1, pb2, pc, pc1, pc2, sgfb, &
      sgfc, stat
    INTEGER, DIMENSION(:), POINTER           :: last_sgf_global, nsgfb, nsgfc
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfb, first_sgfc
    LOGICAL                                  :: failure
    REAL(dp)                                 :: pmax_tmp

    failure = .FALSE.

    natom=SIZE(particle_set)

    ! create a mapping array of easy access to the density matrix elements
    ALLOCATE(last_sgf_global(0:natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    last_sgf_global(0)=0
    DO jatom=1,natom
      jkind = kind_of(jatom)
      last_sgf_global(jatom) = last_sgf_global(jatom-1)+ basis_parameter(jkind)%nsgf_total
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
         jkind = kind_of(jatom)
         nsgfb_total = basis_parameter(jkind)%nsgf_total
         nsetb = basis_parameter(jkind)%nset
         first_sgfb => basis_parameter(jkind)%first_sgf
         nsgfb => basis_parameter(jkind)%nsgf
         pb1 = last_sgf_global(jatom) - nsgfb_total
         DO katom = 1,natom
           kkind = kind_of(katom)
           nsgfc_total = basis_parameter(kkind)%nsgf_total
           nsetc = basis_parameter(kkind)%nset
           first_sgfc => basis_parameter(kkind)%first_sgf
           nsgfc => basis_parameter(kkind)%nsgf
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


  END SUBROUTINE comp_pmax_mat

! *****************************************************************************
!> \brief calculates the screening matrices in the case of periodic boundary
!>      conditions
!> \param screen_xxx screening matrices
!> \param read_from_memory Flag to control reading (when true, data come from buffer)
!> \param buffer stores the matrices
!> \param neighbor_cells contain information for pbc
!> \note
!>      - already computes the square root
!>      - also takes the primitive to sgf transformation into account
!>      - in this way all stored contracted integrals are smaller or equal to their estimate
!> \par History
!>     06.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE comp_screen_mat_per(lib, screen_set, screen_atom, screen_max,screen_periodic,particle_set,&
                                 max_set,potential_parameter, read_from_memory, kind_of, basis_parameter, buffer,&
                                 neighbor_cells, i_thread, cell, error)
    TYPE(lib_int)                            :: lib
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: screen_set
    REAL(dp), DIMENSION(:, :), POINTER       :: screen_atom
    REAL(dp)                                 :: screen_max
    REAL(dp), DIMENSION(:, :, :, :, :), &
      POINTER                                :: screen_periodic
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: max_set
    TYPE(hfx_potential_type)                 :: potential_parameter
    LOGICAL                                  :: read_from_memory
    INTEGER                                  :: kind_of(*)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    REAL(dp), DIMENSION(:, :, :, :, :), &
      OPTIONAL, POINTER                      :: buffer
    TYPE(hfx_cell_type), DIMENSION(:), &
      POINTER                                :: neighbor_cells
    INTEGER, INTENT(IN)                      :: i_thread
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'comp_screen_mat_per', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: jatom, jkind, jset, k, katom, &
                                                kkind, kset, natom, ncob, &
                                                nsetb, nsetc, sgfb, stat
    INTEGER, DIMENSION(:), POINTER           :: lb_max, lb_min, lc_max, &
                                                lc_min, npgfb, npgfc, nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfb
    LOGICAL                                  :: failure
    REAL(dp)                                 :: p2s_b, p2s_c, rb(3), rc(3)
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: max_val
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: contraction_max
    REAL(dp), DIMENSION(:, :), POINTER       :: sphi_b, zetb, zetc

    failure = .FALSE.

    natom=SIZE(particle_set)
    ALLOCATE(screen_set(max_set,max_set,natom,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(screen_atom(natom,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(contraction_max(max_set,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    screen_set=0.0_dp
    screen_atom=0.0_dp
    contraction_max=0.0_dp
    screen_periodic = 0.0_dp
    ALLOCATE(max_val(SIZE(neighbor_cells)),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    max_val = 0.0_dp
    
    DO jatom=1,natom
      jkind = kind_of(jatom)
      lb_max => basis_parameter(jkind)%lmax
      npgfb => basis_parameter(jkind)%npgf
      nsetb = basis_parameter(jkind)%nset
      first_sgfb => basis_parameter(jkind)%first_sgf
      sphi_b => basis_parameter(jkind)%sphi
      nsgfb => basis_parameter(jkind)%nsgf 
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
      jkind = kind_of(jatom)
      lb_max => basis_parameter(jkind)%lmax
      lb_min => basis_parameter(jkind)%lmin
      npgfb => basis_parameter(jkind)%npgf
      nsetb = basis_parameter(jkind)%nset
      zetb => basis_parameter(jkind)%zet
      rb=particle_set(jatom)%r(:)
      DO katom = 1,natom!jatom,natom
        kkind = kind_of(katom)
        lc_max => basis_parameter(kkind)%lmax
        lc_min => basis_parameter(kkind)%lmin
        npgfc => basis_parameter(kkind)%npgf
        nsetc = basis_parameter(kkind)%nset
        zetc => basis_parameter(kkind)%zet
        rc=particle_set(katom)%r(:)
        DO jset = 1,nsetb
          p2s_b = contraction_max(jset,jatom)
          DO kset = 1,nsetc
            p2s_c = contraction_max(kset,katom)
            IF( read_from_memory ) THEN
              screen_periodic(:,jset,kset,jatom,katom) = buffer(:,jset,kset,jatom,katom)
            ELSE
              max_val = 0.0_dp
              CALL screen4_per(lib, rb, rc, &
                               zetb(:,jset), zetc(:,kset), &
                               lb_min(jset), lb_max(jset), lc_min(kset), lc_max(kset),&
                               npgfb(jset), npgfc(kset), &
                               max_val, potential_parameter, neighbor_cells, &
                               i_thread, cell)
                screen_periodic(:,jset,kset,jatom,katom) = SQRT(max_val)*p2s_b*p2s_c
                IF( PRESENT(buffer) ) buffer(:,jset,kset,jatom,katom) = screen_periodic(:,jset,kset,jatom,katom)
              END IF
              screen_set(jset,kset,jatom,katom) = MAXVAL(screen_periodic(:,jset,kset,jatom,katom))
          END DO
        END DO
        screen_atom(jatom,katom)=MAXVAL(screen_set(1:nsetb,1:nsetc,jatom,katom))
      END DO
    END DO
    !! Symmetrize it
    !! temporarily disable this
!    DO jatom=1,natom
!      DO katom=jatom,natom
!        screen_atom(katom,jatom) = screen_atom(jatom,katom)
!        DO jset=1,max_set
!          DO kset=1,max_set
!            screen_set(kset,jset,katom,jatom) = screen_set(jset,kset,jatom,katom)
!          END DO
!        END DO
!      END DO
!    END DO

    screen_max=MAXVAL(screen_atom)
    DEALLOCATE(contraction_max,max_val, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE comp_screen_mat_per

! *****************************************************************************
!> \brief calculates max values of two-electron integrals in a quartet/shell
!>      w.r.t. different zetas using the library lib_int in the case of 
!>      periodic boundary conditions
!> \param r position
!> \param z zeta
!> \param l angular momentum
!> \param npgf number of primitive cartesian gaussian in actual shell
!> \param max_val schwarz screening value
!> \param potential_parameter contains info for libint
!> \param neighbor_cells Array containing the considered cells
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>     06.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE screen4_per(lib, ra, rb, zeta, zetb, &
                         la_min, la_max, lb_min, lb_max,&
                         npgfa, npgfb, &
                         max_val, potential_parameter, neighbor_cells, &
                         i_thread, cell)

    TYPE(lib_int)                            :: lib
    REAL(dp), INTENT(IN)                     :: ra(3), rb(3)
    REAL(dp), DIMENSION(:), INTENT(IN)       :: zeta, zetb
    INTEGER, INTENT(IN)                      :: la_min, la_max, lb_min, &
                                                lb_max, npgfa, npgfb
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: max_val
    TYPE(hfx_potential_type)                 :: potential_parameter
    TYPE(hfx_cell_type), DIMENSION(:), &
      POINTER                                :: neighbor_cells
    INTEGER, INTENT(IN)                      :: i_thread
    TYPE(cell_type), POINTER                 :: cell

    CHARACTER(len=*), PARAMETER :: routineN = 'screen4_per', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ipgf, jpgf, la, lb
    REAL(dp)                                 :: B11(3), im_B(3), &
                                                max_val_temp, pbc_B(3), s(3), &
                                                temp(3)

    temp = rb -ra
    pbc_B = pbc(temp,cell)
    B11 = ra + pbc_B
    DO i = 1,SIZE(neighbor_cells)
      CALL real_to_scaled(s,B11,cell)
      s = s + neighbor_cells(i)%cell
      CALL scaled_to_real(im_B,s,cell)
      DO ipgf = 1,npgfa
        DO jpgf = 1,npgfb
          DO la = la_min,la_max
            DO lb = lb_min, lb_max
              !Build primitives
               max_val_temp = 0.0_dp
               CALL evaluate_eri_screen(lib, ra, im_B, ra, im_B,&
                                        zeta(ipgf), zetb(jpgf), zeta(ipgf), zetb(jpgf),&
                                        la, lb, la ,lb,&
                                        max_val_temp, potential_parameter)
               max_val(i) = MAX(max_val(i), max_val_temp)
            END DO !lb
          END DO !la
        END DO !jpgf
      END DO !ipgf
   END DO
  END SUBROUTINE screen4_per

! *****************************************************************************
!> \brief Estimates the contribution of an ERI taking 1/R dependence into account.
!> \param boxes Buffer to store screening matrices
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE screen_farfield(boxes, particle_set, max_set, &
                             kind_of, basis_parameter, error)
    TYPE(hfx_box_scheme), DIMENSION(:, :), &
      POINTER                                :: boxes
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: max_set
    INTEGER                                  :: kind_of(*)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'screen_farfield', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: jatom, jkind, jpgf, jset, k, &
                                                katom, kkind, kpgf, kset, &
                                                natom, ncob, nsetb, nsetc, &
                                                sgfb, stat
    INTEGER(int_8)                           :: box_idx1, box_idx2
    INTEGER, DIMENSION(:), POINTER           :: lb_max, lb_min, lc_max, &
                                                lc_min, npgfb, npgfc, nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfb
    LOGICAL                                  :: failure
    REAL(dp) :: alpha, K_max, Kbc, p2s_b, p2s_c, px_max, px_min, py_max, &
      py_min, pz_max, pz_min, R_max, rb(3), Rbc, rbc2, rc(3), rp(3), &
      system_size(3), x_max, x_min, y_max, y_min, z_max, z_min, zet1, zet2
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: contraction_max
    REAL(dp), DIMENSION(:, :), POINTER       :: sphi_b, zetb, zetc

    failure = .FALSE.

    natom=SIZE(particle_set)
   
    ALLOCATE(contraction_max(max_set,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    contraction_max=0.0_dp
    alpha =0.0_dp
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
    ALLOCATE(boxes(max_set**2,natom**2), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    boxes(:,:)%center(1) = 0.0_dp
    boxes(:,:)%center(2) = 0.0_dp
    boxes(:,:)%center(3) = 0.0_dp
    boxes(:,:)%K = -HUGE(Kbc)
    boxes(:,:)%R = -HUGE(Rbc)
    boxes(:,:)%alphaInv = HUGE(alpha)
    DO jatom=1,natom
      jkind = kind_of(jatom)
      lb_max => basis_parameter(jkind)%lmax
      npgfb => basis_parameter(jkind)%npgf
      nsetb = basis_parameter(jkind)%nset
      first_sgfb => basis_parameter(jkind)%first_sgf
      sphi_b => basis_parameter(jkind)%sphi
      nsgfb => basis_parameter(jkind)%nsgf
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
      jkind = kind_of(jatom)
      lb_max => basis_parameter(jkind)%lmax
      lb_min => basis_parameter(jkind)%lmin
      npgfb => basis_parameter(jkind)%npgf
      nsetb = basis_parameter(jkind)%nset
      zetb => basis_parameter(jkind)%zet
      rb=particle_set(jatom)%r(:)
      DO katom = 1,natom
        kkind = kind_of(katom)
        lc_max => basis_parameter(kkind)%lmax
        lc_min => basis_parameter(kkind)%lmin
        npgfc => basis_parameter(kkind)%npgf
        nsetc = basis_parameter(kkind)%nset
        zetc => basis_parameter(kkind)%zet
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
                IF( alpha < boxes(box_idx2, box_idx1)%alphaInv ) boxes(box_idx2, box_idx1)%alphaInv = alpha
                IF( Kbc > boxes(box_idx2, box_idx1)%K ) boxes(box_idx2, box_idx1)%K = Kbc*p2s_b*p2s_c*&
                    pi**1.25_dp*SQRT(2.0_dp)
              END DO
            END DO
            !! Store the inverse alpha
            boxes(box_idx2, box_idx1)%alphaInv = 1.0_dp/alpha
            !! Store the radius
            boxes(box_idx2, box_idx1)%radius = INT(1.0_dp/ &
                      (SQRT(2.0_dp*alpha)*erfc_epsilon))

            boxes(box_idx2, box_idx1)%center(1) = (px_max + px_min) / 2.0_dp
            boxes(box_idx2, box_idx1)%center(2) = (py_max + py_min) / 2.0_dp
            boxes(box_idx2, box_idx1)%center(3) = (pz_max + pz_min) / 2.0_dp
            R_max = TINY(R_max)
            DO jpgf = 1,npgfb(jset)
              zet1 = zetb(jpgf,jset)
              DO kpgf = 1,npgfc(kset)
                zet2 = zetc(kpgf,kset)
                alpha = zet1+zet2
                rp = (zet1*rb+zet2*rc)/alpha
                R_max = SQRT((rp(1)-boxes(box_idx2, box_idx1)%center(1))**2 + &
                             (rp(2)-boxes(box_idx2, box_idx1)%center(2))**2 + &
                             (rp(3)-boxes(box_idx2, box_idx1)%center(3))**2)
                IF( R_max > boxes(box_idx2, box_idx1)%R ) boxes(box_idx2, box_idx1)%R = R_max
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
    DEALLOCATE(contraction_max, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE screen_farfield

! *****************************************************************************
!> \brief Estimates the contribution of an ERI taking 1/R dependence into account.
!>      Periodic case
!> \param boxes Buffer to store screening matrices
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE screen_farfield_per(boxes, particle_set, max_set, &
                                 neighbor_cells, kind_of, &
                                 basis_parameter, cell, error)
    TYPE(hfx_box_scheme), DIMENSION(:, :), &
      POINTER                                :: boxes
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: max_set
    TYPE(hfx_cell_type), DIMENSION(:), &
      POINTER                                :: neighbor_cells
    INTEGER                                  :: kind_of(*)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'screen_farfield_per', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, jatom, jkind, jpgf, jset, &
                                                k, katom, kkind, kpgf, kset, &
                                                natom, ncob, nsetb, nsetc, &
                                                sgfb, stat
    INTEGER(int_8)                           :: box_idx1, box_idx2
    INTEGER, DIMENSION(:), POINTER           :: lb_max, lb_min, lc_max, &
                                                lc_min, npgfb, npgfc, nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfb
    LOGICAL                                  :: failure
    REAL(dp) :: alpha, C11(3), im_C(3), K_max, Kbc, p2s_b, p2s_c, pbc_C(3), &
      px_max, px_min, py_max, py_min, pz_max, pz_min, R_max, rb(3), Rbc, &
      rbc2, rc(3), rp(3), s(3), temp(3), zet1, zet2
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: contraction_max
    REAL(dp), DIMENSION(:, :), POINTER       :: sphi_b, zetb, zetc

    failure = .FALSE.
    alpha = 0.0_dp
    natom=SIZE(particle_set)
   
    ALLOCATE(contraction_max(max_set,natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    contraction_max=0.0_dp

    DEALLOCATE(boxes,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(boxes(max_set**2,natom**2), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    boxes(:,:)%center(1) = 0.0_dp
    boxes(:,:)%center(2) = 0.0_dp
    boxes(:,:)%center(3) = 0.0_dp
    boxes(:,:)%K = -HUGE(Kbc)
    boxes(:,:)%R = -HUGE(Rbc)
    boxes(:,:)%alphaInv = HUGE(alpha)
    DO jatom=1,natom
      jkind = kind_of(jatom)
      lb_max => basis_parameter(jkind)%lmax
      npgfb => basis_parameter(jkind)%npgf
      nsetb = basis_parameter(jkind)%nset
      first_sgfb => basis_parameter(jkind)%first_sgf
      sphi_b => basis_parameter(jkind)%sphi
      nsgfb => basis_parameter(jkind)%nsgf
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
      jkind = kind_of(jatom)
      lb_max => basis_parameter(jkind)%lmax
      lb_min => basis_parameter(jkind)%lmin
      npgfb => basis_parameter(jkind)%npgf
      nsetb = basis_parameter(jkind)%nset
      zetb => basis_parameter(jkind)%zet
      rb=particle_set(jatom)%r(:)
      DO katom = 1,natom
        kkind = kind_of(katom)
        lc_max => basis_parameter(kkind)%lmax
        lc_min => basis_parameter(kkind)%lmin
        npgfc => basis_parameter(kkind)%npgf
        nsetc = basis_parameter(kkind)%nset
        zetc => basis_parameter(kkind)%zet
        rc=particle_set(katom)%r(:)
        box_idx1 = (jatom-1)*natom + katom 
        temp = rc-rb
        pbc_C = pbc(temp,cell)
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
              CALL real_to_scaled(s,C11,cell)
              s = s + neighbor_cells(i)%cell
              CALL scaled_to_real(im_C,s,cell)
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
                  IF( alpha < boxes(box_idx2, box_idx1)%alphaInv ) boxes(box_idx2, box_idx1)%alphaInv = alpha
                  IF( Kbc > boxes(box_idx2, box_idx1)%K ) boxes(box_idx2, box_idx1)%K = Kbc*p2s_b*p2s_c*&
                    pi**1.25_dp*SQRT(2.0_dp)
                END DO
              END DO
            END DO
            !! Store inverse alpha
            boxes(box_idx2, box_idx1)%alphaInv = 1.0_dp/alpha
            !! Store the radius
            boxes(box_idx2, box_idx1)%radius = INT(1.0_dp/ &
                      (SQRT(2.0_dp*alpha)*erfc_epsilon))
            boxes(box_idx2, box_idx1)%center(1) = (px_max + px_min) / 2.0_dp
            boxes(box_idx2, box_idx1)%center(2) = (py_max + py_min) / 2.0_dp
            boxes(box_idx2, box_idx1)%center(3) = (pz_max + pz_min) / 2.0_dp
            R_max = TINY(R_max)
            DO i = 1,SIZE(neighbor_cells)
              CALL real_to_scaled(s,C11,cell)
              s = s + neighbor_cells(i)%cell
              CALL scaled_to_real(im_C,s,cell)

              rbc2 = (rb(1)-im_C(1))**2+(rb(2)-im_C(2))**2+(rb(3)-im_C(3))**2
              DO jpgf = 1,npgfb(jset)
                zet1 = zetb(jpgf,jset)
                DO kpgf = 1,npgfc(kset)
                  zet2 = zetc(kpgf,kset)
                  alpha = zet1+zet2
                  rp = (zet1*rb+zet2*im_C)/alpha
                  R_max = SQRT((rp(1)-boxes(box_idx2, box_idx1)%center(1))**2 + &
                               (rp(2)-boxes(box_idx2, box_idx1)%center(2))**2 + &
                               (rp(3)-boxes(box_idx2, box_idx1)%center(3))**2)
                  IF( R_max > boxes(box_idx2, box_idx1)%R ) boxes(box_idx2, box_idx1)%R = R_max
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
    DEALLOCATE(contraction_max, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE screen_farfield_per

! *****************************************************************************
!> \brief updates the maximum of the density matrix in compressed form for screening purposes 
!> \param pmax_set buffer to store matrix
!> \note
!>      - updates for each pair of shells the maximum absolute value of p
!> \par History
!>     09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE update_pmax_mat(pmax_set, full_density, full_density_beta, particle_set, &
                             max_set ,kind_of,basis_parameter,error)

    REAL(dp), DIMENSION(:, :, :, :), POINTER :: pmax_set
    REAL(dp), DIMENSION(:, :), POINTER       :: full_density, &
                                                full_density_beta
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: max_set
    INTEGER                                  :: kind_of(*)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_pmax_mat', &
      routineP = moduleN//':'//routineN

    INTEGER :: jatom, jkind, jset, katom, kkind, kset, mb, mc, natom, nsetb, &
      nsetc, nsgfb_total, nsgfc_total, pb, pb1, pb2, pc, pc1, pc2, sgfb, &
      sgfc, stat
    INTEGER, DIMENSION(:), POINTER           :: last_sgf_global, nsgfb, nsgfc
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfb, first_sgfc
    LOGICAL                                  :: failure
    REAL(dp)                                 :: pmax_tmp

    failure = .FALSE.

    natom=SIZE(particle_set)

    ! create a mapping array of easy access to the density matrix elements
    ALLOCATE(last_sgf_global(0:natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    last_sgf_global(0)=0
    DO jatom=1,natom
      jkind = kind_of(jatom)
      last_sgf_global(jatom) = last_sgf_global(jatom-1)+basis_parameter(jkind)%nsgf_total
    END DO

    pmax_set = 0.0_dp

    DO jatom=1,natom
      jkind = kind_of(jatom)
      nsgfb_total = basis_parameter(jkind)%nsgf_total
      nsetb = basis_parameter(jkind)%nset
      first_sgfb => basis_parameter(jkind)%first_sgf
      nsgfb => basis_parameter(jkind)%nsgf
      pb1 = last_sgf_global(jatom) - nsgfb_total

      DO katom = jatom,natom
        kkind = kind_of(katom)
        nsgfc_total = basis_parameter(kkind)%nsgf_total
        nsetc = basis_parameter(kkind)%nset
        first_sgfc => basis_parameter(kkind)%first_sgf
        nsgfc => basis_parameter(kkind)%nsgf
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

    !! Symmetrize it
    DO jatom=1,natom
      DO katom=jatom,natom
        DO jset=1,max_set
          DO kset=1,max_set
            pmax_set(kset,jset,katom,jatom) = pmax_set(jset,kset,jatom,katom)
          END DO
        END DO
      END DO
    END DO

    DEALLOCATE(last_sgf_global,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE update_pmax_mat

END MODULE hfx_screening_methods
