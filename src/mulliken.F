!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief compute mulliken charges
!>      we (currently) define them as c_i = 1/2 [ (PS)_{ii} + (SP)_{ii} ]
!> \author Joost VandeVondele March 2003
! *****************************************************************************
MODULE mulliken
  USE atomic_charges,                  ONLY: print_atomic_charges,&
                                             print_multipoles
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp_control_types,                ONLY: mulliken_restraint_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_get_block_p, cp_dbcsr_get_info, cp_dbcsr_iterator, &
       cp_dbcsr_iterator_blocks_left, cp_dbcsr_iterator_next_block, &
       cp_dbcsr_iterator_start, cp_dbcsr_iterator_stop, cp_dbcsr_p_type, &
       cp_dbcsr_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE semi_empirical_int_arrays,       ONLY: se_orbital_pointer
  USE semi_empirical_mpole_types,      ONLY: semi_empirical_mpole_type
  USE semi_empirical_types,            ONLY: get_se_param,&
                                             semi_empirical_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mulliken'

! *** Public subroutines ***

  PUBLIC :: mulliken_charges, mulliken_restraint, multipoles

  INTERFACE mulliken_charges
     MODULE PROCEDURE mulliken_charges_a, mulliken_charges_b, mulliken_charges_c
  END INTERFACE

  INTERFACE multipoles
     MODULE PROCEDURE multipoles_a, multipoles_b
  END INTERFACE

CONTAINS

! *****************************************************************************
!> \brief computes the energy and density matrix derivate of a constraint on the
!>      mulliken charges
!>
!>      optional outputs:
!>      computes energy (added)
!>      contribution to KS matrix (added)
!>      contribution to W  matrix (added)
!> \param mulliken_restraint_control additional parameters needed to control the restraint
!> \param para_env para_env of the matrices
!> \param s_matrix ,p_matrix : containing the respective quantities
!> \param p_matrix ...
!> \param energy ...
!> \param order_p ...
!> \param ks_matrix ...
!> \param w_matrix ...
!> \param error ...
!> \par History
!>      06.2004 created [Joost VandeVondele]
!> \note
!>      contribution to the KS matrix is derivative wrt P
!>      contribution to the W matrix is derivate wrt S (sign?)
!>      needed for orbital and ionic forces respectively
! *****************************************************************************
  SUBROUTINE mulliken_restraint(mulliken_restraint_control,para_env, &
                                s_matrix,p_matrix,energy,order_p,ks_matrix,w_matrix,error)
    TYPE(mulliken_restraint_type), &
      INTENT(IN)                             :: mulliken_restraint_control
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_dbcsr_type), POINTER             :: s_matrix
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: p_matrix
    REAL(KIND=dp), OPTIONAL                  :: energy, order_p
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: ks_matrix, w_matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mulliken_restraint', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, iblock_col, iblock_row, &
                                                ispin, nblock, nspin
    LOGICAL                                  :: found
    REAL(kind=dp)                            :: mult, my_energy, my_order_p
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: charges, charges_deriv, &
                                                ks_block, p_block, s_block, &
                                                w_block
    TYPE(cp_dbcsr_iterator)                  :: iter

! here we get the numbers for charges

    nspin=SIZE(p_matrix)
    CALL cp_dbcsr_get_info(s_matrix,nblkrows_total=nblock)

    ALLOCATE(charges(nblock,nspin))
    ALLOCATE(charges_deriv(nblock,nspin))
    CALL compute_charges(p_matrix,s_matrix,charges,para_env)
    !
    ! this can be used to check the correct implementation of the derivative
    ! CALL rf_deriv_check(mulliken_restraint_control,charges)
    !
    CALL restraint_functional(mulliken_restraint_control, &
                              charges,charges_deriv,my_energy,my_order_p)

    IF (PRESENT(order_p)) THEN
        order_p=my_order_p
    ENDIF
    IF (PRESENT(energy)) THEN
        energy=my_energy
    ENDIF

    IF (PRESENT(ks_matrix)) THEN

       DO ispin=1,nspin
          CALL cp_dbcsr_iterator_start(iter, s_matrix)
          DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
             CALL cp_dbcsr_iterator_next_block(iter, iblock_row, iblock_col, s_block,blk)
             CALL cp_dbcsr_get_block_p(matrix=ks_matrix(ispin)%matrix,&
                  row=iblock_row,col=iblock_col,BLOCK=ks_block,found=found)

             IF ( .NOT. (ASSOCIATED(s_block) .AND. ASSOCIATED(ks_block) ) ) THEN
                CALL stop_program(routineN,moduleN,__LINE__,&
                                  "Unexpected s / ks structure",para_env)
             END IF
                mult = 0.5_dp * charges_deriv(iblock_row,ispin) + &
                       0.5_dp * charges_deriv(iblock_col,ispin)

                ks_block=ks_block+mult*s_block

          ENDDO
          CALL cp_dbcsr_iterator_stop(iter)
       ENDDO

    ENDIF

    IF (PRESENT(w_matrix)) THEN

       DO ispin=1,nspin
          CALL cp_dbcsr_iterator_start(iter, p_matrix(ispin)%matrix)
          DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
             CALL cp_dbcsr_iterator_next_block(iter, iblock_row, iblock_col, p_block, blk)
             CALL cp_dbcsr_get_block_p(matrix=w_matrix(ispin)%matrix,&
                  row=iblock_row,col=iblock_col,BLOCK=w_block,found=found)

             ! we can cycle if a block is not present
             IF ( .NOT. (ASSOCIATED(w_block) .AND. ASSOCIATED(p_block) ) ) CYCLE

             ! minus sign relates to convention for W
             mult = -0.5_dp * charges_deriv(iblock_row,ispin)   &
                    -0.5_dp * charges_deriv(iblock_col,ispin)

             w_block=w_block+mult*p_block

          END DO
          CALL cp_dbcsr_iterator_stop(iter)
       ENDDO

    ENDIF

    DEALLOCATE(charges)
    DEALLOCATE(charges_deriv)

  END SUBROUTINE mulliken_restraint

! *****************************************************************************
!> \brief computes energy and derivatives given a set of charges
!>       this implementation uses the spin density on a number of atoms
!>       as a penalty function
!> \param mulliken_restraint_control ...
!> \param charges (nblock,nspin)
!> \param charges_deriv derivate wrt the corresponding charge entry
!> \param energy ...
!> \param order_p ...
!> \par History
!>      06.2004 created [Joost VandeVondele]
!>      02.2005 added more general form [Joost VandeVondele]
!> \note
!>       should be easy to adapt for other specialized cases
! *****************************************************************************
SUBROUTINE restraint_functional(mulliken_restraint_control,charges, &
                                charges_deriv,energy,order_p)
    TYPE(mulliken_restraint_type), &
      INTENT(IN)                             :: mulliken_restraint_control
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: charges, charges_deriv
    REAL(KIND=dp), INTENT(OUT)               :: energy, order_p

    INTEGER                                  :: I
    REAL(KIND=dp)                            :: dum

    charges_deriv=0.0_dp
    order_p=0.0_dp

    DO I=1,mulliken_restraint_control%natoms
       order_p=order_p+charges(mulliken_restraint_control%atoms(I),1) &
                      -charges(mulliken_restraint_control%atoms(I),2)  ! spin density on the relevant atoms
    ENDDO
    ! energy
    energy=mulliken_restraint_control%strength*(order_p-mulliken_restraint_control%target)**2
    ! derivative
    dum=2*mulliken_restraint_control%strength*(order_p-mulliken_restraint_control%target)
    DO I=1,mulliken_restraint_control%natoms
       charges_deriv(mulliken_restraint_control%atoms(I),1)=  dum
       charges_deriv(mulliken_restraint_control%atoms(I),2)= -dum
    ENDDO
END SUBROUTINE restraint_functional

! *****************************************************************************
!> \brief checks numerically the analytic derivative in restraint_functional
!> \param mulliken_restraint_control ...
!> \param charges (nblock,nspin)
!> \par History
!>      06.2004 created [Joost VandeVondele]
!> \note
!>       for testing only, uncomment the call in mulliken restraint
! *****************************************************************************
SUBROUTINE rf_deriv_check(mulliken_restraint_control,charges)
    TYPE(mulliken_restraint_type), &
      INTENT(IN)                             :: mulliken_restraint_control
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: charges

    REAL(KIND=dp), PARAMETER                 :: stepsize = 0.00001_dp

    INTEGER                                  :: I, J
    REAL(KIND=dp)                            :: dum, energy1, energy2
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: charges_deriv, charges_local, &
                                                deriv_num

    ALLOCATE(charges_local(1:SIZE(charges,1),1:SIZE(charges,2)))
    ALLOCATE(charges_deriv(1:SIZE(charges,1),1:SIZE(charges,2)))
    ALLOCATE(deriv_num(1:SIZE(charges,1),1:SIZE(charges,2)))

    ! finite difference check
    charges_local=charges
    do_i: DO I=1,SIZE(charges,1)
      do_j: DO J=1,SIZE(charges,2)
         ! x+dx
         charges_local(I,J)=charges_local(I,J)+stepsize
         CALL restraint_functional(mulliken_restraint_control,charges_local, &
                                   charges_deriv,energy1,dum)
         charges_local(I,J)=charges_local(I,J)-stepsize ! back to orig
         ! x-dx
         charges_local(I,J)=charges_local(I,J)-stepsize
         CALL restraint_functional(mulliken_restraint_control,charges_local, &
                                   charges_deriv,energy2,dum)
         charges_local(I,J)=charges_local(I,J)+stepsize ! back to orig
         ! num deriv
         deriv_num(I,J)=(energy1-energy2)/(2*stepsize)
      ENDDO do_j
    ENDDO do_i
    ! analytic results vs deriv_num
    CALL restraint_functional(mulliken_restraint_control,charges, &
                                   charges_deriv,energy2,dum)
    WRITE(*,*) "deriv check : max abs diff, max abs", MAXVAL(ABS(charges_deriv-deriv_num)), &
                                                      MAXVAL(ABS(charges_deriv))

    DEALLOCATE(charges_local,charges_deriv,deriv_num)

END SUBROUTINE rf_deriv_check

! *****************************************************************************
!> \brief compute the mulliken charges
!> \param p_matrix , s_matrix, para_env
!> \param s_matrix ...
!> \param charges previously allocated with the right size (natom,nspin)
!> \param para_env ...
!> \par History
!>      06.2004 created [Joost VandeVondele]
!> \note
!>      charges are computed per spin in the LSD case
! *****************************************************************************
  SUBROUTINE compute_charges(p_matrix,s_matrix,charges,para_env)
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: p_matrix
    TYPE(cp_dbcsr_type), POINTER             :: s_matrix
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: charges
    TYPE(cp_para_env_type), POINTER          :: para_env

    INTEGER                                  :: blk, iblock_col, iblock_row, &
                                                ispin, nblock, nspin
    LOGICAL                                  :: found
    REAL(kind=dp)                            :: mult
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_block, s_block
    TYPE(cp_dbcsr_iterator)                  :: iter

! here we get the numbers for charges

    nspin=SIZE(p_matrix)
    CALL cp_dbcsr_get_info(s_matrix,nblkrows_total=nblock)

    charges=0.0_dp
    DO ispin=1,nspin
       CALL cp_dbcsr_iterator_start(iter, s_matrix)
       DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
          CALL cp_dbcsr_iterator_next_block(iter, iblock_row, iblock_col, s_block, blk)
          CALL cp_dbcsr_get_block_p(matrix=p_matrix(ispin)%matrix,&
               row=iblock_row,col=iblock_col,BLOCK=p_block,found=found)

          ! we can cycle if a block is not present
          IF ( .NOT. (ASSOCIATED(s_block) .AND. ASSOCIATED(p_block) ) ) CYCLE

          IF (iblock_row.EQ.iblock_col) THEN
             mult=0.5_dp ! avoid double counting of diagonal blocks
          ELSE
             mult=1.0_dp
          ENDIF
          charges(iblock_row,ispin)=charges(iblock_row,ispin) + &
                                            mult*SUM(p_block*s_block)
          charges(iblock_col,ispin)=charges(iblock_col,ispin) + &
                                            mult*SUM(p_block*s_block)

       ENDDO
       CALL cp_dbcsr_iterator_stop(iter)
    ENDDO
    CALL mp_sum(charges,para_env%group)

  END SUBROUTINE compute_charges

! *****************************************************************************
!> \brief compute the mulliken charge derivatives
!> \param p_matrix , s_matrix, para_env
!> \param s_matrix ...
!> \param charges ...
!> \param dcharges previously allocated with the right size (natom,3)
!> \param para_env ...
!> \par History
!>      01.2012 created [JHU]
! *****************************************************************************
  SUBROUTINE compute_dcharges(p_matrix,s_matrix,charges,dcharges,para_env)
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: p_matrix, s_matrix
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: charges, dcharges
    TYPE(cp_para_env_type), POINTER          :: para_env

    INTEGER                                  :: blk, iblock_col, iblock_row, &
                                                ider, ispin, nblock, nspin
    LOGICAL                                  :: found
    REAL(kind=dp)                            :: mult
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ds_block, p_block, s_block
    TYPE(cp_dbcsr_iterator)                  :: iter

! here we get the numbers for charges

    nspin=SIZE(p_matrix)
    CALL cp_dbcsr_get_info(s_matrix(1)%matrix,nblkrows_total=nblock)

    charges=0.0_dp
    dcharges=0.0_dp
    DO ispin=1,nspin
       CALL cp_dbcsr_iterator_start(iter, s_matrix(1)%matrix)
       DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
          CALL cp_dbcsr_iterator_next_block(iter, iblock_row, iblock_col, s_block, blk)
          CALL cp_dbcsr_get_block_p(matrix=p_matrix(ispin)%matrix,&
               row=iblock_row,col=iblock_col,BLOCK=p_block,found=found)

          ! we can cycle if a block is not present
          IF ( .NOT. (ASSOCIATED(s_block) .AND. ASSOCIATED(p_block) ) ) CYCLE

          IF (iblock_row.EQ.iblock_col) THEN
             mult=0.5_dp ! avoid double counting of diagonal blocks
          ELSE
             mult=1.0_dp
          ENDIF
          charges(iblock_row,ispin)=charges(iblock_row,ispin) + mult*SUM(p_block*s_block)
          charges(iblock_col,ispin)=charges(iblock_col,ispin) + mult*SUM(p_block*s_block)
          DO ider=1,3
             CALL cp_dbcsr_get_block_p(matrix=s_matrix(ider+1)%matrix,&
                  row=iblock_row,col=iblock_col,BLOCK=ds_block,found=found)
             dcharges(iblock_row,ider)=dcharges(iblock_row,ider) + mult*SUM(p_block*ds_block)
             dcharges(iblock_col,ider)=dcharges(iblock_col,ider) + mult*SUM(p_block*ds_block)
          END DO

       ENDDO
       CALL cp_dbcsr_iterator_stop(iter)
    ENDDO
    CALL mp_sum(charges,para_env%group)
    CALL mp_sum(dcharges,para_env%group)

  END SUBROUTINE compute_dcharges


! *****************************************************************************
!> \brief print the mulliken charges to scr on ionode
!> \param p_matrix , s_matrix, para_env
!> \param s_matrix ...
!> \param para_env ...
!> \param particle_set (needed for Z)
!> \param qs_kind_set ...
!> \param scr unit for output
!> \param title ...
!> \param error ...
!> \par History
!>      06.2004 adapted to remove explicit matrix multiply [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE mulliken_charges_a(p_matrix,s_matrix,para_env,particle_set, &
                                qs_kind_set, scr,title,error)

    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: p_matrix
    TYPE(cp_dbcsr_type), POINTER             :: s_matrix
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    INTEGER                                  :: scr
    CHARACTER(LEN=*)                         :: title
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mulliken_charges_a', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nblock, nspin
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: charges

    CALL timeset(routineN,handle)

    ! here we get the numbers for charges
    nspin=SIZE(p_matrix)
    CALL cp_dbcsr_get_info(s_matrix,nblkrows_total=nblock)

    ALLOCATE(charges(nblock,nspin))
    CALL compute_charges(p_matrix,s_matrix,charges,para_env)

    CALL print_atomic_charges(particle_set,qs_kind_set,scr,title,electronic_charges=charges)

    DEALLOCATE(charges)

    CALL timestop(handle)

  END SUBROUTINE mulliken_charges_a

! *****************************************************************************
!> \brief ...
!> \param p_matrix ...
!> \param s_matrix ...
!> \param para_env ...
!> \param mcharge ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE mulliken_charges_b(p_matrix,s_matrix,para_env,mcharge,error)

    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: p_matrix
    TYPE(cp_dbcsr_type), POINTER             :: s_matrix
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: mcharge
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mulliken_charges_b', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    CALL compute_charges(p_matrix,s_matrix,mcharge,para_env)

    CALL timestop(handle)

  END SUBROUTINE mulliken_charges_b

! *****************************************************************************
!> \brief ...
!> \param p_matrix ...
!> \param s_matrix ...
!> \param para_env ...
!> \param mcharge ...
!> \param dmcharge ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE mulliken_charges_c(p_matrix,s_matrix,para_env,mcharge,dmcharge,error)

    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: p_matrix, s_matrix
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: mcharge, dmcharge
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mulliken_charges_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    CALL compute_dcharges(p_matrix,s_matrix,mcharge,dmcharge,para_env)

    CALL timestop(handle)

  END SUBROUTINE mulliken_charges_c

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param p_matrix ...
!> \param pscp ...
!> \param s_matrix ...
!> \param charge ...
!> \param dipole ...
!> \param quadrupole ...
!> \param para_env ...
! *****************************************************************************
  SUBROUTINE compute_multipoles(atomic_kind_set,qs_kind_set,particle_set,&
                                p_matrix,pscp,s_matrix,charge,dipole,quadrupole,para_env)
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: p_matrix
    TYPE(cp_dbcsr_type), POINTER             :: pscp, s_matrix
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charge
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dipole
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: quadrupole
    TYPE(cp_para_env_type), POINTER          :: para_env

    INTEGER :: blk, iatom, iblock_col, iblock_row, ikind, indi, indj, ispin, &
      jint, natorb_a, nblock, nspin, size_1c_int
    LOGICAL                                  :: defined, found
    REAL(kind=dp)                            :: fac
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_block, pscp_block, s_block
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(semi_empirical_mpole_type), POINTER :: mpole
    TYPE(semi_empirical_type), POINTER       :: se_kind_a

    nspin=SIZE(p_matrix)
    CALL cp_dbcsr_get_info(s_matrix,nblkrows_total=nblock)

    charge=     0.0_dp
    dipole=     0.0_dp
    quadrupole= 0.0_dp

    DO ispin=1,nspin
       CALL cp_dbcsr_iterator_start(iter, s_matrix)
       DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
          CALL cp_dbcsr_iterator_next_block(iter, iblock_row, iblock_col, s_block, blk)
          CALL cp_dbcsr_get_block_p(matrix=p_matrix(ispin)%matrix,&
               row=iblock_row,col=iblock_col,BLOCK=p_block,found=found)
          ! we can cycle if a block is not present
          IF ( .NOT. (ASSOCIATED(s_block) .AND. ASSOCIATED(p_block) ) ) CYCLE
          IF (iblock_row.NE.iblock_col) CYCLE

          iatom=iblock_row

          CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
          CALL get_qs_kind(qs_kind_set(ikind), se_parameter=se_kind_a)
          CALL get_se_param(se_kind_a,defined=defined,natorb=natorb_a)

          CALL cp_dbcsr_get_block_p(matrix=p_matrix(ispin)%matrix,&
               row=iatom,col=iatom,BLOCK=p_block,found=found)
          CALL cp_dbcsr_get_block_p(matrix=pscp,&
               row=iatom,col=iatom,BLOCK=pscp_block,found=found)

          charge(iatom)=charge(iatom)+se_kind_a%zeff
          size_1c_int = SIZE(se_kind_a%w_mpole)
          DO jint = 1, size_1c_int

             mpole => se_kind_a%w_mpole(jint)%mpole
             indi  = se_orbital_pointer(mpole%indi)
             indj  = se_orbital_pointer(mpole%indj)
             fac   = 1.0_dp
             IF (indi/=indj) fac = 2.0_dp

             IF(mpole%task(1))THEN
               charge(iatom) = charge(iatom)       + fac*p_block(indi,indj)*mpole%c
               charge(iatom) = charge(iatom)       + fac*pscp_block(indi,indj)*mpole%c
             END IF
             IF(mpole%task(2))THEN
               dipole(:,iatom) = dipole(:,iatom)   + fac*p_block(indi,indj)*mpole%d(:)
               dipole(:,iatom) = dipole(:,iatom)   + fac*pscp_block(indi,indj)*mpole%d(:)
             END IF
             IF(mpole%task(3))THEN
               quadrupole(:,:,iatom) = quadrupole(:,:,iatom) + fac*p_block(indi,indj)*mpole%qc(:,:)/3.0_dp
               quadrupole(:,:,iatom) = quadrupole(:,:,iatom) + fac*pscp_block(indi,indj)*mpole%qc(:,:)/3.0_dp
             END IF

             END DO

       ENDDO
       CALL cp_dbcsr_iterator_stop(iter)
    ENDDO
    CALL mp_sum(charge,para_env%group)
    CALL mp_sum(dipole,para_env%group)
    CALL mp_sum(quadrupole,para_env%group)

  END SUBROUTINE compute_multipoles

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param p_matrix ...
!> \param pscp ...
!> \param s_matrix ...
!> \param para_env ...
!> \param scr ...
!> \param title ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE multipoles_a(atomic_kind_set,qs_kind_set,particle_set,&
                          p_matrix,pscp,s_matrix,para_env,scr,title,error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: p_matrix
    TYPE(cp_dbcsr_type), POINTER             :: pscp, s_matrix
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: scr
    CHARACTER(LEN=*)                         :: title
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'multipoles_a', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nblock, nspin
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charge
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dipole
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: quadrupole

    CALL timeset(routineN,handle)

    ! here we get the numbers for charges
    nspin=SIZE(p_matrix)
    CALL cp_dbcsr_get_info(s_matrix,nblkrows_total=nblock)

    ALLOCATE(charge(nblock))
    ALLOCATE(dipole(3,nblock))
    ALLOCATE(quadrupole(3,3,nblock))

    CALL compute_multipoles(atomic_kind_set,qs_kind_set,particle_set,&
                            p_matrix,pscp,s_matrix,charge,dipole,quadrupole,para_env)

    CALL print_multipoles(particle_set,qs_kind_set,scr,title,charge=charge,&
                          dipole=dipole,quadrupole=quadrupole)

    DEALLOCATE(charge)
    DEALLOCATE(dipole)
    DEALLOCATE(quadrupole)

    CALL timestop(handle)

  END SUBROUTINE multipoles_a

! *****************************************************************************
!> \brief ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param p_matrix ...
!> \param pscp ...
!> \param s_matrix ...
!> \param para_env ...
!> \param mcharge ...
!> \param mdipole ...
!> \param mquadrupole ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE multipoles_b(atomic_kind_set,qs_kind_set,particle_set,&
                          p_matrix,pscp,s_matrix,para_env,&
                          mcharge,mdipole,mquadrupole,error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: p_matrix
    TYPE(cp_dbcsr_type), POINTER             :: pscp, s_matrix
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mcharge
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: mdipole
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: mquadrupole
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'multipoles_b', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    CALL compute_multipoles(atomic_kind_set,qs_kind_set,particle_set,&
                            p_matrix,pscp,s_matrix,mcharge,mdipole,mquadrupole,para_env)

    CALL timestop(handle)

  END SUBROUTINE multipoles_b

END MODULE mulliken
