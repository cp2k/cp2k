!------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations          !
!   Copyright (C) 2000  CP2K developers group                                  !
!------------------------------------------------------------------------------!
!!****** cp2k/tbmd_initialize [1.0] *
!!
!!   NAME
!!     tbmd_initialize
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
!
    MODULE tbmd_initialize

  USE kinds,                           ONLY: dbl
  USE particle_types,                  ONLY: particle_type
  USE tbmd_types,                      ONLY: tb_hopping_type,&
                                             tbatom_param_type
  USE termination,                     ONLY: stop_memory

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: tbmd_init, tb_get_numel, tb_get_basize
!
!!*****
!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
!!****** tbmd_initialize/tbmd_init [1.0] *
!!
!!   NAME
!!     tbmd_init
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

      SUBROUTINE tbmd_init(part,tbatom,tbhop)
    TYPE(particle_type), INTENT(IN)          :: part(:)
    TYPE(tbatom_param_type), DIMENSION(:), &
      INTENT(INOUT)                          :: tbatom
    TYPE(tb_hopping_type), DIMENSION(:, :), &
      INTENT(INOUT)                          :: tbhop

    INTEGER                                  :: i, iat, ii, ios, is, j, jat, &
                                                js, ju, k, n, natom_type, ni, &
                                                nj

!
! count number of atoms per atom type

        natom_type = size(tbatom)
        DO iat = 1, natom_type
          tbatom(iat) %natoms = 0
        END DO
        DO iat = 1, size(part)
          i = part(iat) %prop%ptype
          tbatom(i) %natoms = tbatom(i) %natoms + 1
        END DO
! count number of basis function per atom type
! set up diagonal part of hamiltonian
        DO iat = 1, natom_type
          n = tbatom(iat) %bsf(0) + 3*tbatom(iat) %bsf(1) + &
            5*tbatom(iat) %bsf(2) + 7*tbatom(iat) %bsf(3)
          tbatom(iat) %nbsf = n
          ALLOCATE (tbatom(iat)%h_diag(n),STAT=ios)
          IF (ios/=0) CALL stop_memory('tbmd_initialize', &
            'tbatom(iat)%h_diag',n)
          ii = 0
          DO i = 0, 3
            DO j = 1, tbatom(iat) %bsf(i)
              tbatom(iat) %h_diag(ii+1:ii+2*i+1) = tbatom(iat) &
                %orbital_energy(i,j)
              ii = ii + 2*i + 1
            END DO
          END DO
        END DO
! set up shell pointers
        DO iat = 1, natom_type
          n = tbatom(iat) %bsf(0) + tbatom(iat) %bsf(1) + &
            tbatom(iat) %bsf(2) + tbatom(iat) %bsf(3)
          ALLOCATE (tbatom(iat)%spoint(n,2),STAT=ios)
          IF (ios/=0) CALL stop_memory('tbmd_initialize', &
            'tbatom(iat)%spoint',2*n)
          ii = 0
          DO i = 0, 3
            DO j = 1, tbatom(iat) %bsf(i)
              ii = ii + 1
              tbatom(iat) %spoint(ii,2) = i
            END DO
          END DO
          tbatom(iat) %spoint(1,1) = 0
          DO j = 2, n
            ii = tbatom(iat) %spoint(j-1,2)
            tbatom(iat) %spoint(j,1) = tbatom(iat) %spoint(j-1,1) + 2*ii + 1
          END DO
        END DO
! pointer: shell pair to hopping potential
        DO iat = 1, natom_type
          DO jat = 1, natom_type
            ni = size(tbatom(iat)%spoint(:,1))
            nj = size(tbatom(jat)%spoint(:,1))
            ALLOCATE (tbhop(iat,jat)%llpoint(ni,nj),STAT=ios)
            IF (ios/=0) CALL stop_memory('tbmd_initialize', &
              'tbatom(iat)%llpoint',ni*nj)
            k = 0
            DO i = 1, ni
              is = tbatom(iat) %spoint(i,2)
              ju = 1
              IF (iat==jat) ju = i
              DO j = ju, nj
                js = tbatom(jat) %spoint(j,2)
                tbhop(iat,jat) %llpoint(i,j) = k + 1
                k = k + min(is,js)
              END DO
            END DO
          END DO
        END DO
!
      END SUBROUTINE tbmd_init
!!*****
!------------------------------------------------------------------------------!
      FUNCTION tb_get_basize(tbatom) RESULT (basis_size)
    TYPE(tbatom_param_type), DIMENSION(:), &
      INTENT(IN)                             :: tbatom
    INTEGER                                  :: basis_size

    INTEGER                                  :: i

        basis_size = 0
        DO i = 1, size(tbatom)
          basis_size = basis_size + tbatom(i) %nbsf*tbatom(i) %natoms
        END DO
      END FUNCTION tb_get_basize
!------------------------------------------------------------------------------!
      FUNCTION tb_get_numel(tbatom,charge) RESULT (number_of_electrons)
    TYPE(tbatom_param_type), DIMENSION(:), &
      INTENT(IN)                             :: tbatom
    REAL(dbl), INTENT(IN)                    :: charge
    INTEGER                                  :: number_of_electrons

    INTEGER                                  :: i

        number_of_electrons = 0
        DO i = 1, size(tbatom)
          number_of_electrons = number_of_electrons + &
            tbatom(i) %valence_charge*tbatom(i) %natoms
        END DO
        number_of_electrons = number_of_electrons + charge
      END FUNCTION tb_get_numel
!------------------------------------------------------------------------------!
    END MODULE tbmd_initialize
!------------------------------------------------------------------------------!
