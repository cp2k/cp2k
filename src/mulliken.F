!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2003-2005  CP2K group                                       !
!-----------------------------------------------------------------------------!


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
  USE atomic_charges,                  ONLY: print_atomic_charges
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE cp_control_types,                ONLY: mulliken_restraint_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE sparse_matrix_types,             ONLY: first_block_node,&
                                             get_block_node,&
                                             get_matrix_info,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_p_type,&
                                             real_matrix_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mulliken'

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
!!    - mulliken_restraint_control  : additional parameters needed to control the restraint
!!    - para_env                    : para_env of the matrices 
!!    - s_matrix,p_matrix : containing the respective quantities
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [Joost VandeVondele]
!!
!!   SOURCE
!!*** **********************************************************************
  SUBROUTINE mulliken_restraint(mulliken_restraint_control,para_env, &
                                s_matrix,p_matrix,energy,order_p,ks_matrix,w_matrix)
    TYPE(mulliken_restraint_type), &
      INTENT(IN)                             :: mulliken_restraint_control
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(real_matrix_type), POINTER          :: s_matrix
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: p_matrix
    REAL(KIND=dp), OPTIONAL                  :: energy, order_p
    TYPE(real_matrix_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: ks_matrix, w_matrix

    INTEGER                                  :: iblock_col, iblock_row, &
                                                ispin, nblock, nspin
    REAL(kind=dp)                            :: mult, my_energy, my_order_p
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: charges, charges_deriv, &
                                                ks_block, p_block, s_block, &
                                                w_block
    TYPE(real_block_node_type), POINTER      :: block_node

! here we get the numbers for charges

    nspin=SIZE(p_matrix)
    CALL get_matrix_info(s_matrix,nblock_row=nblock)

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

                IF ( .NOT. (ASSOCIATED(s_block) .AND. ASSOCIATED(ks_block) ) ) THEN
                         CALL stop_program("Mulliken","unexpected s / ks structure")
                ENDIF
    
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
                IF ( .NOT. (ASSOCIATED(w_block) .AND. ASSOCIATED(p_block) ) ) THEN
                         CALL stop_program("Mulliken","unexpected p / w structure")
                ENDIF
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
!!      this implementation uses the spin density on a number of atoms 
!!      as a penalty function
!!     
!!   NOTES
!!      should be easy to adapt for other specialized cases
!!
!!   INPUTS
!!    - charges(nblock,nspin)
!!    - output :
!!        - energy
!!        - charges_deriv : derivate wrt the corresponding charge entry
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [Joost VandeVondele]
!!     02.2005 added more general form [Joost VandeVondele]
!!
!!   SOURCE
!!*** **********************************************************************
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

!!****f* mulliken/rf_deriv_check *
!!
!!   NAME
!!      rf_deriv_check
!!   FUNCTION
!!      checks numerically the analytic derivative in restraint_functional
!!
!!   NOTES
!!      for testing only, uncomment the call in mulliken restraint
!!
!!   INPUTS
!!    - charges(nblock,nspin)
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [Joost VandeVondele]
!!
!!   SOURCE
!!*** **********************************************************************
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
    WRITE(6,*) "deriv check : max abs diff, max abs", MAXVAL(ABS(charges_deriv-deriv_num)), &
                                                      MAXVAL(ABS(charges_deriv))

    DEALLOCATE(charges_local,charges_deriv,deriv_num)

END SUBROUTINE rf_deriv_check

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
      POINTER                                :: p_matrix
    TYPE(real_matrix_type), POINTER          :: s_matrix
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: charges
    TYPE(cp_para_env_type), POINTER          :: para_env

    INTEGER                                  :: iblock_col, iblock_row, &
                                                ispin, nblock, nspin
    REAL(kind=dp)                            :: mult
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_block, s_block
    TYPE(real_block_node_type), POINTER      :: block_node

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

             IF ( .NOT. (ASSOCIATED(s_block) .AND. ASSOCIATED(p_block))) THEN
                         CALL stop_program("Mulliken","unexpected p / w structure")
             ENDIF

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
  SUBROUTINE mulliken_charges(p_matrix,s_matrix,para_env,particle_set, &
                              scr,title)

    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: p_matrix
    TYPE(real_matrix_type), POINTER          :: s_matrix
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER                                  :: scr
    CHARACTER(LEN=*)                         :: title

    CHARACTER(LEN=2)                         :: element_symbol
    INTEGER                                  :: handle, i, iblock, nblock, &
                                                nspin
    REAL(KIND=dp)                            :: total_charge, zeff
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: charges

! GO

    CALL timeset("mulliken_charges","I","",handle)

    ! here we get the numbers for charges
    nspin=SIZE(p_matrix)
    CALL get_matrix_info(s_matrix,nblock_row=nblock)

    ALLOCATE(charges(nblock,nspin))
    CALL compute_charges(p_matrix,s_matrix,charges,para_env)

    CALL print_atomic_charges(particle_set,scr,title,electronic_charges=charges)

    DEALLOCATE(charges)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE

END MODULE mulliken
