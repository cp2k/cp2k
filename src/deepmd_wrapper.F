!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface to the DeePMD-kit or a c++ wrapper.
!> \par History
!>      07.2019 created [Yongbin Zhuang]
!>      06.2021 refactored [Yunpei Liu]
!>      10.2023 adapt to DeePMD-kit C Interface [Yunpei Liu]
!> \author Yongbin Zhuang
! **************************************************************************************************

MODULE deepmd_wrapper

   USE ISO_C_BINDING,                   ONLY: C_CHAR,&
                                              C_DOUBLE,&
                                              C_INT,&
                                              C_LOC,&
                                              C_NULL_CHAR,&
                                              C_NULL_PTR,&
                                              C_PTR
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: dp_deep_pot, dp_deep_pot_compute, deepmd_model_release

#if(__DEEPMD)
   INTERFACE
      FUNCTION dp_deep_pot_c(c_model) BIND(C, name="DP_NewDeepPot")
         USE ISO_C_BINDING, ONLY: C_PTR, &
                                  C_CHAR, &
                                  C_DOUBLE, &
                                  C_INT, &
                                  C_NULL_CHAR, &
                                  C_LOC
      CHARACTER(LEN=1, KIND=C_CHAR)                      :: c_model(*)
      TYPE(C_PTR)                                        :: dp_deep_pot_c

      END FUNCTION

      SUBROUTINE dp_deep_pot_compute_c(dp, natom, &
                                       coord, atype, cell, &
                                       energy, force, virial, &
                                       atomic_energy, atomic_virial) BIND(C, name="DP_DeepPotCompute")
         USE ISO_C_BINDING, ONLY: C_PTR, &
                                  C_INT, &
                                  C_CHAR, &
                                  C_DOUBLE
      TYPE(C_PTR), VALUE                                 :: dp
      INTEGER(C_INT), VALUE                              :: natom
      REAL(C_DOUBLE), DIMENSION(natom, 3), INTENT(IN)    :: coord
      INTEGER(C_INT), DIMENSION(natom), INTENT(IN)       :: atype
      REAL(C_DOUBLE), DIMENSION(9), INTENT(IN)           :: cell
      REAL(C_DOUBLE), INTENT(OUT)                        :: energy
      REAL(C_DOUBLE), DIMENSION(natom, 3), INTENT(OUT)   :: force
      REAL(C_DOUBLE), DIMENSION(9), INTENT(OUT)          :: virial
      REAL(C_DOUBLE), DIMENSION(natom), INTENT(OUT)      :: atomic_energy
      REAL(C_DOUBLE), DIMENSION(natom, 9), INTENT(OUT)   :: atomic_virial

      END SUBROUTINE
   END INTERFACE
#endif

CONTAINS

! **************************************************************************************************
!> \brief Load DP from a model file.
!> \param model Path to the model file.
!> \return Pointer to the DP model.
! **************************************************************************************************

   FUNCTION dp_deep_pot(model) RESULT(pot)
      CHARACTER(len=*), INTENT(INOUT)                    :: model
      TYPE(C_PTR)                                        :: pot

#if(__DEEPMD)
      pot = dp_deep_pot_c(c_model=TRIM(model)//C_NULL_CHAR)
#else
      MARK_USED(model)
      pot = C_NULL_PTR
#endif
   END FUNCTION dp_deep_pot

! **************************************************************************************************
!> \brief Compute energy, force and virial from DP.
!> \param dp Pointer to the DP model.
!> \param natom Number of atoms.
!> \param coord Coordinates of the atoms.
!> \param atype Atom types.
!> \param cell Cell vectors.
!> \param energy Potential energy.
!> \param force Forces.
!> \param virial Virial tensor.
!> \param atomic_energy Atomic energies.
!> \param atomic_virial Atomic virial tensors.
! **************************************************************************************************

   SUBROUTINE dp_deep_pot_compute(dp, natom, coord, atype, cell, energy, force, virial, atomic_energy, atomic_virial)
      USE, INTRINSIC :: ISO_C_BINDING
      TYPE(C_PTR), VALUE :: dp
      INTEGER(C_INT), VALUE :: natom
      REAL(C_DOUBLE), DIMENSION(natom, 3), INTENT(IN) :: coord
      INTEGER(C_INT), DIMENSION(natom), INTENT(IN) :: atype
      REAL(C_DOUBLE), DIMENSION(9), INTENT(IN), OPTIONAL :: cell
      REAL(C_DOUBLE), INTENT(OUT), OPTIONAL :: energy
      REAL(C_DOUBLE), DIMENSION(natom, 3), INTENT(OUT), OPTIONAL :: force
      REAL(C_DOUBLE), DIMENSION(9), INTENT(OUT), OPTIONAL :: virial
      REAL(C_DOUBLE), DIMENSION(natom), INTENT(OUT), OPTIONAL :: atomic_energy
      REAL(C_DOUBLE), DIMENSION(natom, 9), INTENT(OUT), OPTIONAL :: atomic_virial

#if(__DEEPMD)
      CALL dp_deep_pot_compute_c(dp, natom, coord, atype, cell, energy, force, virial, atomic_energy, atomic_virial)
#else
      CPABORT("CP2K was compiled without libdeepmd_c library.")
      MARK_USED(dp)
      MARK_USED(natom)
      MARK_USED(coord)
      MARK_USED(atype)
      MARK_USED(cell)
      MARK_USED(energy)
      MARK_USED(force)
      MARK_USED(virial)
      MARK_USED(atomic_energy)
      MARK_USED(atomic_virial)
#endif
   END SUBROUTINE

! **************************************************************************************************
!> \brief Releases a deepmd model and all its ressources.
!> \param model Pointer to the DP model.
! **************************************************************************************************
   SUBROUTINE deepmd_model_release(model)
      USE ISO_C_BINDING, ONLY: C_PTR, &
                               C_NULL_PTR
      TYPE(C_PTR), INTENT(INOUT)                         :: model

      model = C_NULL_PTR
   END SUBROUTINE deepmd_model_release

END MODULE deepmd_wrapper
