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
  USE cp_fm_basic_linalg,              ONLY: cp_sm_fm_multiply, &
                                             cp_sm_plus_fm_fm_t
  USE cp_fm_types,                     ONLY: cp_fm_type, &
                                             cp_fm_create, &
                                             cp_fm_release
  USE kinds,                           ONLY: wp => dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE qs_mo_types,                     ONLY: mo_set_p_type, &
                                             get_mo_set
  USE sparse_matrix_types,             ONLY: real_matrix_type, &
                                             real_matrix_p_type,&
                                             set_matrix, &
                                             replicate_matrix, &
                                             deallocate_matrix, &
                                             get_block_node, &
                                             get_matrix_info
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: module_name = "mulliken"

! *** Public subroutines ***

  PUBLIC :: mulliken_charges

! *****************************************************************************

CONTAINS

! computes mulliken charges, avoiding full nao x nao matrices
! we give a breakup in alpha and beta spin
! still needs generalization for smear

  SUBROUTINE mulliken_charges(mo_array,s_matrix,particle_set,scr,ionode)
    TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mo_array
    TYPE(real_matrix_type), POINTER                    :: s_matrix
    TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
    LOGICAL                                            :: ionode
    INTEGER                                            :: scr

    ! locals

    INTEGER :: nspin,ispin,handle,homo
    TYPE(real_matrix_p_type), POINTER, DIMENSION(:)    :: charge_matrix
    TYPE(cp_fm_type), POINTER                 :: fm_temp,eigenvectors
    REAL(KIND=wp), DIMENSION(:,:), POINTER             :: block
    REAL(KIND=wp), DIMENSION(:,:), POINTER              :: charges
        CHARACTER(LEN=2) :: element_symbol
        integer :: zeff

    INTEGER :: nblock,iblock,i
    REAL(KIND=wp)       :: total_charge

    ! GO
    CALL timeset("mulliken_charges","I","",handle)
    NULLIFY(charge_matrix,fm_temp,block,eigenvectors)

    nspin=SIZE(mo_array)
    ALLOCATE(charge_matrix(nspin)) 
    DO ispin=1,nspin
       NULLIFY(charge_matrix(ispin)%matrix)
    ENDDO

    ! compute the sparse SP avoiding full matrix multiplies (we really only need the diag, but well..)
    DO ispin=1,nspin
      CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
                      homo=homo,&
                      eigenvectors=eigenvectors)
       CALL cp_fm_create(fm_temp,eigenvectors%matrix_struct)
       CALL cp_sm_fm_multiply(s_matrix,eigenvectors,fm_temp, &
                              homo,eigenvectors%matrix_struct%para_env)
       
       CALL replicate_matrix(s_matrix,charge_matrix(ispin)%matrix,"charge_matrix")
       CALL set_matrix(charge_matrix(ispin)%matrix,0.0_wp)
       ! still need an 'occ' number here.
       CALL cp_sm_plus_fm_fm_t(sparse_matrix=charge_matrix(ispin)%matrix,&
                               matrix_v=fm_temp,&
                               matrix_g=eigenvectors,&
                               first_col=1,&
                               last_col=homo)
       CALL cp_fm_release(fm_temp)
    ENDDO

    ! here we get the numbers for charges

    CALL get_matrix_info(charge_matrix(1)%matrix,nblock_row=nblock)
    ALLOCATE(charges(nblock,nspin))
    charges=0.0_wp

    DO ispin=1,nspin
      DO iblock=1,nblock
          CALL get_block_node(charge_matrix(ispin)%matrix,iblock,iblock,block=block)
          IF (ASSOCIATED(block)) THEN
            DO i=1,size(block,1)
               charges(iblock,ispin)=charges(iblock,ispin)+block(i,i)
            ENDDO
          ENDIF
      ENDDO
    ENDDO

    CALL mp_sum(charges,eigenvectors%matrix_struct%para_env%group)

    ! obviously output is not yet very nice

    IF (ionode) THEN
       write(scr,'(A)') ""
       write(scr,'(A)') "  Mulliken charges "
       SELECT CASE (nspin)
       CASE(1)
       write(scr,'(A)') "  Atom     |    Charge"
       CASE DEFAULT
       write(scr,'(A)') "  Atom     |    Charge | Spin diff charge"
       END SELECT
       total_charge = 0.0_wp
       if (size(particle_set) .NE. nblock) call stop_program("Mulliken","unexpected number of atoms")
       DO iblock=1,nblock
          CALL get_atomic_kind(atomic_kind=particle_set(iblock)%atomic_kind,&
                             element_symbol=element_symbol,&
                             zeff=zeff)

          SELECT CASE (nspin)
          CASE(1)
          write(scr,'(I6,A2,A2,F12.6)') iblock,"  ",element_symbol,zeff-2.0*charges(iblock,1)
          total_charge=total_charge+ zeff-2.0*charges(iblock,1)
          CASE DEFAULT
          write(scr,'(I6,A2,A2,2F12.6)') iblock,"  ",element_symbol, &
                      zeff-(charges(iblock,1)+charges(iblock,2)),(charges(iblock,1)-charges(iblock,2))
          total_charge=total_charge+ zeff-(charges(iblock,1)+charges(iblock,2))
          END SELECT
       ENDDO
       write(scr,'(A,F10.6)') "  Total     ",total_charge

    ENDIF

    ! give back data
 
    DO ispin=1,nspin
       CALL deallocate_matrix(charge_matrix(ispin)%matrix)
    ENDDO

    DEALLOCATE(charge_matrix)
    CALL timestop(0.0_wp,handle)

  END SUBROUTINE

END MODULE mulliken
