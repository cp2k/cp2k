!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2003  CP2K group                                            !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****** cp2k/mulliken [1.0] *
!!
!!   NAME
!!     mulliken
!!
!!   FUNCTION
!!     compute mulliken charges
!!     we (currently) define them as c_i = 1/2 [ (PS)_{ii} + (SP)_{ii} ]
!!
!!   AUTHOR
!!     Joost VandeVondele March 2003
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE mulliken
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply,&
                                             cp_sm_plus_fm_fm_t
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE sparse_matrix_types,             ONLY: deallocate_matrix,&
                                             get_block_node,&
                                             get_matrix_info,&
                                             first_block_node, &
                                             next_block_node, &
                                             real_block_node_type,&
                                             real_matrix_p_type,&
                                             real_matrix_type,&
                                             replicate_matrix,&
                                             set_matrix
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: module_name = "mulliken"

! *** Public subroutines ***

  PUBLIC :: mulliken_charges, mulliken_restraint

! *****************************************************************************

CONTAINS
!!****f* mulliken/mulliken_restraint *
!!
!!   NAME
!!     mulliken_restraint
!!
!!   FUNCTION
!!     computes the energy and density matrix derivate of a constraint on the 
!!     mulliken charges
!!
!!     optional outputs:
!!     computes energy (added)
!!     contribution to KS matrix (added)
!!     contribution to W  matrix (added)
!!
!!   NOTES
!!     contribution to the KS matrix is derivative wrt P
!!     contribution to the W matrix is derivate wrt S (sign?)
!!     needed for orbital and ionic forces respectively
!!
!!   INPUTS
!!    - para_env                    : para_env of the matrices 
!!    - s_matrix,p_matrix : containing the respective quantities
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [Joost VandeVondele]
!!
!!   SOURCE
!!*** **********************************************************************
  SUBROUTINE mulliken_restraint(para_env,s_matrix,p_matrix,energy,ks_matrix,w_matrix)
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(real_matrix_type), POINTER          :: s_matrix
    TYPE(real_matrix_p_type), DIMENSION(:), &
                              POINTER        :: p_matrix

    TYPE(real_matrix_p_type), DIMENSION(:), &
                    POINTER, OPTIONAL        :: ks_matrix, w_matrix
    REAL(KIND=dp), OPTIONAL                  :: energy

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: charges
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: charges_deriv
    INTEGER :: iblock_row, iblock_col, ispin, nspin, nblock
    REAL(kind=dp)                            :: mult, my_energy, strength
    TYPE(real_block_node_type), POINTER      :: block_node
    REAL(kind=dp), DIMENSION(:,:), POINTER   :: s_block, ks_block, p_block, w_block

    ! here we get the numbers for charges
    nspin=SIZE(p_matrix)
    CALL get_matrix_info(s_matrix,nblock_row=nblock)

    ALLOCATE(charges(nblock,nspin))
    ALLOCATE(charges_deriv(nblock,nspin))
    CALL compute_charges(p_matrix,s_matrix,charges,para_env)
    CALL restraint_functional(charges,charges_deriv,my_energy)

    IF (PRESENT(energy)) THEN 
        energy=energy+my_energy
    ENDIF

    IF (PRESENT(ks_matrix)) THEN
       DO ispin=1,nspin
          DO iblock_row=1,nblock
             block_node => first_block_node(s_matrix,iblock_row)
             DO WHILE (ASSOCIATED(block_node))
                CALL get_block_node(block_node=block_node,&
                                    block_col=iblock_col,&
                                    BLOCK=s_block)
                CALL get_block_node(matrix=ks_matrix(ispin)%matrix,&
                                    block_row=iblock_row,&
                                    block_col=iblock_col,&
                                    BLOCK=ks_block)
    
                mult = 0.5_dp * charges_deriv(iblock_row,ispin) + &
                       0.5_dp * charges_deriv(iblock_col,ispin)
   
                ks_block=ks_block+mult*s_block
   
                block_node => next_block_node(block_node)
             END DO
          END DO
       ENDDO
    ENDIF

    IF (PRESENT(w_matrix)) THEN
       DO ispin=1,nspin
          DO iblock_row=1,nblock
             block_node => first_block_node(p_matrix(ispin)%matrix,iblock_row)
             DO WHILE (ASSOCIATED(block_node))
                CALL get_block_node(block_node=block_node,&
                                    block_col=iblock_col,&
                                    BLOCK=p_block)
                CALL get_block_node(matrix=w_matrix(ispin)%matrix,&
                                    block_row=iblock_row,&
                                    block_col=iblock_col,&
                                    BLOCK=w_block)
                ! minus sign relates to convention for W 
                mult = -0.5_dp * charges_deriv(iblock_row,ispin)   &
                       -0.5_dp * charges_deriv(iblock_col,ispin)
   
                w_block=w_block+mult*p_block
   
                block_node => next_block_node(block_node)
             END DO
          END DO
       ENDDO
    ENDIF

    DEALLOCATE(charges)
    DEALLOCATE(charges_deriv)
  END SUBROUTINE mulliken_restraint

!!****f* mulliken/restraint_functional *
!!
!!   NAME
!!      restraint_functional
!!   FUNCTION
!!      computes energy and derivatives given a set of charges
!!     
!!   NOTES
!!      still hard coded right now
!!   INPUTS
!!    - charges(nblock,nspin)
!!    - output :
!!        - energy
!!        - charges_deriv : derivate wrt the corresponding charge entry
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [Joost VandeVondele]
!!
!!   SOURCE
!!*** **********************************************************************
SUBROUTINE restraint_functional(charges,charges_deriv,energy)
    REAL(KIND=dp), DIMENSION(:,:), POINTER :: charges, charges_deriv
    REAL(KIND=dp)                          :: energy
    REAL(KIND=dp)                          :: dum,strength

    charges_deriv=0.0_dp
    energy=0.0_dp
    strength=0.4_dp

    dum=(charges(2,1)-charges(2,2)) - 1.0_dp
    energy=energy+strength*dum**2
    charges_deriv(2,1)=  strength*2*dum
    charges_deriv(2,2)= -strength*2*dum
END SUBROUTINE restraint_functional

!!****f* mulliken/compute_charges *
!!
!!   NAME
!!     compute_charges
!!
!!   FUNCTION
!!     compute the mulliken charges
!!
!!   NOTES
!!     charges are computed per spin in the LSD case
!!
!!   INPUTS
!!    - p_matrix, s_matrix, para_env
!!    - charges : previously allocated with the right size (natom,nspin)
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [Joost VandeVondele]
!!
!!   SOURCE
!!*** **********************************************************************
  SUBROUTINE compute_charges(p_matrix,s_matrix,charges,para_env)
    TYPE(real_matrix_p_type), DIMENSION(:), &
                              POINTER        :: p_matrix
    TYPE(real_matrix_type), POINTER          :: s_matrix
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: charges
    TYPE(cp_para_env_type), POINTER          :: para_env

    INTEGER :: iblock_row, iblock_col, ispin, nspin, nblock
    REAL(kind=dp)                            :: mult
    TYPE(real_block_node_type), POINTER      :: block_node
    REAL(kind=dp), DIMENSION(:,:), POINTER   :: s_block, p_block

    ! here we get the numbers for charges
    nspin=SIZE(p_matrix)
    CALL get_matrix_info(s_matrix,nblock_row=nblock)

    charges=0.0_dp
    DO ispin=1,nspin
       DO iblock_row=1,nblock
          block_node => first_block_node(s_matrix,iblock_row)
          DO WHILE (ASSOCIATED(block_node))
             CALL get_block_node(block_node=block_node,&
                                 block_col=iblock_col,&
                                 BLOCK=s_block)
             CALL get_block_node(matrix=p_matrix(ispin)%matrix,&
                                 block_row=iblock_row,&
                                 block_col=iblock_col,&
                                 BLOCK=p_block)

             IF (iblock_row.EQ.iblock_col) THEN
                mult=0.5_dp ! avoid double counting of diagonal blocks
             ELSE
                mult=1.0_dp
             ENDIF
             charges(iblock_row,ispin)=charges(iblock_row,ispin) + &
                                               mult*SUM(p_block*s_block)
             charges(iblock_col,ispin)=charges(iblock_col,ispin) + &
                                               mult*SUM(p_block*s_block)

             block_node => next_block_node(block_node)
          END DO
       END DO
    ENDDO
    CALL mp_sum(charges,para_env%group)

  END SUBROUTINE compute_charges


!!****f* mulliken/mulliken_charges *
!!
!!   NAME
!!     mulliken_charges
!!
!!   FUNCTION
!!     print the mulliken charges to scr on ionode
!!
!!   NOTES
!!
!!   INPUTS
!!    - p_matrix, s_matrix, para_env
!!    - particle_set (needed for Z)
!!    - scr : unit for output
!!    - ionode : true if scr is connected, and output is needed
!!
!!   MODIFICATION HISTORY
!!     06.2004 adapted to remove explicit matrix multiply [Joost VandeVondele]
!!
!!   SOURCE
!!*** **********************************************************************
  SUBROUTINE mulliken_charges(p_matrix,s_matrix,para_env,particle_set,scr,ionode)

    TYPE(real_matrix_p_type), DIMENSION(:), &
                              POINTER        :: p_matrix
    TYPE(real_matrix_type), POINTER          :: s_matrix
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER                                  :: scr
    LOGICAL                                  :: ionode

    CHARACTER(LEN=2)                         :: element_symbol
    INTEGER                                  :: handle, i, iblock, &
                                                ispin, nblock, nspin, zeff
    REAL(KIND=dp)                            :: total_charge
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: charges

! GO

    CALL timeset("mulliken_charges","I","",handle)

    ! here we get the numbers for charges
    nspin=SIZE(p_matrix)
    CALL get_matrix_info(s_matrix,nblock_row=nblock)

    ALLOCATE(charges(nblock,nspin))
    CALL compute_charges(p_matrix,s_matrix,charges,para_env)

    ! obviously output is not yet very nice
    IF (ionode) THEN
       WRITE(scr,'(A)') ""
       WRITE(scr,'(A)') "  Mulliken charges "
       SELECT CASE (nspin)
       CASE(1)
       WRITE(scr,'(A)') "  Atom     |    Charge"
       CASE DEFAULT
       WRITE(scr,'(A)') "  Atom     |    Charge | Spin diff charge"
       END SELECT
       total_charge = 0.0_dp
       IF (SIZE(particle_set) .NE. nblock) CALL stop_program("Mulliken","unexpected number of atoms")
       DO iblock=1,nblock
          CALL get_atomic_kind(atomic_kind=particle_set(iblock)%atomic_kind,&
                             element_symbol=element_symbol,&
                             zeff=zeff)

          SELECT CASE (nspin)
          CASE(1)
          WRITE(scr,'(I6,A2,A2,F12.6)') iblock,"  ",element_symbol,zeff-charges(iblock,1)
          total_charge=total_charge+ zeff-charges(iblock,1)
          CASE DEFAULT
          WRITE(scr,'(I6,A2,A2,2F12.6)') iblock,"  ",element_symbol, &
                      zeff-(charges(iblock,1)+charges(iblock,2)),(charges(iblock,1)-charges(iblock,2))
          total_charge=total_charge+ zeff-(charges(iblock,1)+charges(iblock,2))
          END SELECT
       ENDDO
       WRITE(scr,'(A,F10.6)') "  Total     ",total_charge

    ENDIF
    DEALLOCATE(charges)
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE

END MODULE mulliken
