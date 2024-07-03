!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of charge equilibration in xTB
!> \author JGH
! **************************************************************************************************
MODULE xtb_eeq
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE atprop_types,                    ONLY: atprop_array_init,&
                                              atprop_type
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              xtb_control_type
   USE cp_fm_basic_linalg,              ONLY: cp_fm_solve
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_element,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_set_element,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: oorootpi
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xtb_eeq'

   PUBLIC :: xtb_eeq_calculation

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param charges ...
!> \param cnumbers ...
!> \param eeq_energy ...
! **************************************************************************************************
   SUBROUTINE xtb_eeq_calculation(qs_env, charges, cnumbers, eeq_energy)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: charges
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: cnumbers
      REAL(KIND=dp), INTENT(INOUT)                       :: eeq_energy

      CHARACTER(len=*), PARAMETER :: routineN = 'xtb_eeq_calculation'

      INTEGER                                            :: handle, iac, iar, iatom, ic, ikind, ir, &
                                                            jkind, natom, ncloc, nkind, nrloc, ns
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      INTEGER, DIMENSION(:), POINTER                     :: cind, rind
      LOGICAL                                            :: defined, do_ewald
      REAL(KIND=dp)                                      :: ala, alb, dr, gam, grc, kappa, lambda, &
                                                            totalcharge, xi, xr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: adiag, chia
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gab
      REAL(KIND=dp), DIMENSION(3)                        :: ri, rij, rj
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: mat_struct, vec_struct
      TYPE(cp_fm_type)                                   :: eeq_mat, rhs_vec
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(xtb_atom_type), POINTER                       :: xtb_atom_a, xtb_atom_b
      TYPE(xtb_control_type), POINTER                    :: xtb_control

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set, &
                      particle_set=particle_set, &
                      cell=cell, &
                      atprop=atprop, &
                      dft_control=dft_control)
      CALL get_qs_env(qs_env, nkind=nkind, natom=natom)

      xtb_control => dft_control%qs_control%xtb_control

      totalcharge = dft_control%charge

      IF (atprop%energy) THEN
         CALL atprop_array_init(atprop%atecoul, natom)
      END IF

      ! gamma[a,b]
      ALLOCATE (gab(nkind, nkind))
      gab = 0.0_dp
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, defined=defined)
         IF (.NOT. defined) CYCLE
         CALL get_xtb_atom_param(xtb_atom_a, alpg=ala)
         DO jkind = 1, nkind
            CALL get_qs_kind(qs_kind_set(jkind), xtb_parameter=xtb_atom_b)
            CALL get_xtb_atom_param(xtb_atom_b, defined=defined)
            IF (.NOT. defined) CYCLE
            CALL get_xtb_atom_param(xtb_atom_b, alpg=alb)
            !
            gab(ikind, jkind) = SQRT(1._dp/(ala*ala + alb*alb))
            !
         END DO
      END DO

      ! Chi[a,a]
      ALLOCATE (chia(natom), adiag(natom))
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, xi=xi, kappa0=kappa, eta=gam)
         !
         chia(iatom) = xi - kappa*SQRT(cnumbers(iatom))
         adiag(iatom) = gam + 2.0_dp*gab(ikind, ikind)*oorootpi
         !
      END DO

      do_ewald = xtb_control%do_ewald
      IF (do_ewald) THEN
         ! response matrix
         ns = natom
         CALL get_qs_env(qs_env=qs_env, para_env=para_env, blacs_env=blacs_env)
         CALL cp_fm_struct_create(mat_struct, context=blacs_env, para_env=para_env, &
                                  nrow_global=ns, ncol_global=ns)
         CALL cp_fm_create(eeq_mat, mat_struct)
      ELSE
         ns = natom + 1
         CALL get_qs_env(qs_env=qs_env, para_env=para_env, blacs_env=blacs_env)
         CALL cp_fm_struct_create(mat_struct, context=blacs_env, para_env=para_env, &
                                  nrow_global=ns, ncol_global=ns)
         CALL cp_fm_create(eeq_mat, mat_struct)
         CALL cp_fm_struct_create(vec_struct, context=blacs_env, para_env=para_env, &
                                  nrow_global=ns, ncol_global=1)
         CALL cp_fm_create(rhs_vec, vec_struct)
         !
         CALL cp_fm_get_info(eeq_mat, nrow_local=nrloc, ncol_local=ncloc, &
                             row_indices=rind, col_indices=cind)
         ! set up matrix
         CALL cp_fm_set_all(eeq_mat, 1.0_dp, 0.0_dp)
         CALL cp_fm_set_all(rhs_vec, 0.0_dp)
         DO ir = 1, nrloc
            iar = rind(ir)
            IF (iar > natom) CYCLE
            ri(1:3) = particle_set(iar)%r(1:3)
            DO ic = 1, ncloc
               iac = cind(ic)
               IF (iac > natom) CYCLE
               rj(1:3) = particle_set(iac)%r(1:3)
               IF (iar == iac) THEN
                  grc = adiag(iar)
               ELSE
                  ikind = kind_of(iar)
                  jkind = kind_of(iac)
                  rij(1:3) = ri(1:3) - rj(1:3)
                  dr = SQRT(SUM(rij**2))
                  grc = erf(gab(ikind, jkind)*dr)/dr
               END IF
               CALL cp_fm_set_element(eeq_mat, iar, iac, grc)
            END DO
         END DO
         ! set up rhs vector
         DO ir = 1, nrloc
            iar = rind(ir)
            IF (iar > natom) THEN
               xr = totalcharge
            ELSE
               xr = -chia(iar)
            END IF
            CALL cp_fm_set_element(rhs_vec, iar, 1, xr)
         END DO
         !
         ! fm matrices are overwritten!
         CALL cp_fm_solve(eeq_mat, rhs_vec)
         !
         charges = 0.0_dp
         lambda = 0.0_dp
         DO ir = 1, nrloc
            iar = rind(ir)
            IF (iar <= natom) THEN
               CALL cp_fm_get_element(rhs_vec, iar, 1, xr)
            ELSE
               CALL cp_fm_get_element(rhs_vec, iar, 1, lambda)
            END IF
            charges(iar) = xr
         END DO
         CALL para_env%sum(lambda)
         CALL para_env%sum(charges)
         !
         ! energy:  -0.5*q^T.X - lambda*totalcharge
         eeq_energy = 0.5*SUM(charges(1:natom)*chia(1:natom)) - lambda*totalcharge
         !
      END IF

      CALL cp_fm_struct_release(mat_struct)
      CALL cp_fm_struct_release(vec_struct)
      CALL cp_fm_release(eeq_mat)
      CALL cp_fm_release(rhs_vec)

      DEALLOCATE (gab, chia, adiag)

      CALL timestop(handle)

   END SUBROUTINE xtb_eeq_calculation

END MODULE xtb_eeq

