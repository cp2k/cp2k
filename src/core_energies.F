!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/core_energies [1.0] *
!!
!!   NAME
!!     core_energies
!!
!!   FUNCTION
!!     Calculation of the energies concerning the core charge distribution
!!
!!   AUTHOR
!!     Matthias Krack (27.04.2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE core_energies

! *****************************************************************************

  USE kinds, ONLY: wp => dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: calculate_ecore,&
            calculate_ecore_overlap,&
            calculate_ecore_self

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE calculate_ecore(h,p,ecore)

!   Purpose: Calculate the core Hamiltonian energy which includes the kinetic
!            and the potential energy of the electrons.

!   ***************************************************************************

    USE matrix_types, ONLY: first_block_node,&
                            get_block_node,&
                            get_matrix_info,&
                            next_block_node,&
                            real_block_node_type,&
                            real_matrix_set_type

    TYPE(real_matrix_set_type), INTENT(IN) :: h,p
    REAL(wp), INTENT(OUT)                  :: ecore

!   *** Local variables ***

    TYPE(real_block_node_type), POINTER :: h_block_node
    REAL(wp)                            :: ecore_block
    INTEGER                             :: i,iblock_col,iblock_row,j,nblock_row

    REAL(wp), DIMENSION(:,:), POINTER :: h_block,p_block

!   ---------------------------------------------------------------------------

    ecore = 0.0_wp

    CALL get_matrix_info(matrix=h%matrix,nblock_row=nblock_row)

    DO iblock_row=1,nblock_row

      h_block_node => first_block_node(matrix=h%matrix,block_row=iblock_row)

      DO WHILE (ASSOCIATED(h_block_node))

        CALL get_block_node(block_node=h_block_node,&
                            block_col=iblock_col,&
                            block=h_block)

        CALL get_block_node(matrix=p%matrix,&
                            block_row=iblock_row,&
                            block_col=iblock_col,&
                            block=p_block)

        ecore_block = 0.0_wp

        DO j=1,SIZE(h_block,2)
          DO i=1,SIZE(h_block,1)
            ecore_block = ecore_block + p_block(i,j)*h_block(i,j)
          END DO
        END DO

        IF (iblock_row == iblock_col) THEN
          ecore = ecore + ecore_block
        ELSE
          ecore = ecore + 2.0_wp*ecore_block
        END IF

        h_block_node => next_block_node(h_block_node)

      END DO

    END DO

  END SUBROUTINE calculate_ecore

! *****************************************************************************

  SUBROUTINE calculate_ecore_overlap(ecore_overlap)

!   Purpose: Calculate the overlap energy of the core charge distribution.

!   ***************************************************************************

    USE atomic_kinds,    ONLY: kind_info
    USE atoms,           ONLY: atom_info,natom
    USE cell_parameters, ONLY: abc

    REAL(wp), INTENT(OUT) :: ecore_overlap

!   *** Local variables ***

    REAL(wp) :: alpha_a,alpha_ab,alpha_b,dab,dabmax,radius_a,radius_b,za,zb,zab
    INTEGER  :: atom_a,atom_b,iatom,icell,ikind,jatom,jcell,jkind,kcell

    REAL(wp), DIMENSION(3) :: ra,rab,rb
    INTEGER, DIMENSION(3)  :: ncell

    REAL(wp), EXTERNAL :: erfc

!   ---------------------------------------------------------------------------

    ecore_overlap = 0.0_wp

    DO iatom=1,natom

      ikind = atom_info(iatom)%kind

      alpha_a = kind_info(ikind)%alpha_core_charge
      radius_a = kind_info(ikind)%core_charge_radius
      za = REAL(kind_info(ikind)%zeff,wp)

      ra(:) = atom_info(iatom)%r_pbc(:)

      DO jatom=iatom+1,natom

        jkind = atom_info(jatom)%kind

        alpha_b = kind_info(jkind)%alpha_core_charge
        radius_b = kind_info(jkind)%core_charge_radius
        zb = REAL(kind_info(jkind)%zeff,wp)

        rb(:) = atom_info(jatom)%r_pbc(:)

        dabmax = radius_a + radius_b

        ncell(:) = CEILING(dabmax/abc(:))

        alpha_ab = SQRT(alpha_a*alpha_b/(alpha_a + alpha_b))
        zab = za*zb

        DO icell=-ncell(1),ncell(1)
          DO jcell=-ncell(2),ncell(2)
            DO kcell=-ncell(3),ncell(3)

!             *** Calculate the distance vector between the two atoms ***

              rab(:) = rb(:) - ra(:) + REAL((/icell,jcell,kcell/),wp)*abc(:)

!             *** Calculate the interatomic distance ***

              dab = SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))

              IF (dab < dabmax) THEN
                ecore_overlap = ecore_overlap + zab*erfc(alpha_ab*dab)/dab
              END IF

            END DO
          END DO
        END DO

      END DO

    END DO

  END SUBROUTINE calculate_ecore_overlap

! *****************************************************************************

  SUBROUTINE calculate_ecore_self(ecore_self)

!   Purpose: Calculate the self energy of the core charge distribution.

!   ***************************************************************************

    USE atomic_kinds,  ONLY: kind_info,nkind
    USE mathconstants, ONLY: twopi

    REAL(wp), INTENT(OUT) :: ecore_self

!   *** Local variables ***

    INTEGER :: ikind

!   ---------------------------------------------------------------------------

    ecore_self = 0.0_wp

    DO ikind=1,nkind
      ecore_self = ecore_self - REAL(kind_info(ikind)%natom,wp)*&
                                REAL(kind_info(ikind)%zeff,wp)**2*&
                                SQRT(kind_info(ikind)%alpha_core_charge)
    END DO

    ecore_self = ecore_self/SQRT(twopi)

  END SUBROUTINE calculate_ecore_self

! *****************************************************************************

END MODULE core_energies
