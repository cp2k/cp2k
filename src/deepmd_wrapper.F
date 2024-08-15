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
                                              C_NULL_CHAR,&
                                              C_NULL_PTR,&
                                              C_PTR
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: deepmd_model_type, deepmd_model_load, deepmd_model_compute, deepmd_model_release

   TYPE deepmd_model_type
      PRIVATE
      TYPE(C_PTR)                          :: c_ptr = C_NULL_PTR
   END TYPE deepmd_model_type

CONTAINS

! **************************************************************************************************
!> \brief Load DP from a model file.
!> \param filename Path to the model file.
!> \return Pointer to the DP model.
! **************************************************************************************************
   FUNCTION deepmd_model_load(filename) RESULT(model)
      CHARACTER(len=*), INTENT(INOUT)                    :: filename
      TYPE(deepmd_model_type)                            :: model

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'deepmd_model_load'

      INTEGER                                            :: handle
      INTERFACE
         FUNCTION NewDeepPot(filename) BIND(C, name="DP_NewDeepPot")
            IMPORT :: C_PTR, C_CHAR
            CHARACTER(kind=C_CHAR), DIMENSION(*)               :: filename
            TYPE(C_PTR)                                        :: NewDeepPot
         END FUNCTION
      END INTERFACE

      CALL timeset(routineN, handle)

#if defined(__DEEPMD)
      model%c_ptr = NewDeepPot(filename=TRIM(filename)//C_NULL_CHAR)
#else
      CPABORT("CP2K was compiled without libdeepmd_c library.")
      MARK_USED(filename)
      MARK_USED(model)
#endif

      CALL timestop(handle)
   END FUNCTION deepmd_model_load

! **************************************************************************************************
!> \brief Compute energy, force and virial from DP.
!> \param model Pointer to the DP model.
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
   SUBROUTINE deepmd_model_compute(model, natom, coord, atype, cell, energy, force, virial, &
                                   atomic_energy, atomic_virial)
      TYPE(deepmd_model_type)                            :: model
      INTEGER                                            :: natom
      REAL(kind=dp), DIMENSION(natom, 3), INTENT(IN)     :: coord
      INTEGER, DIMENSION(natom), INTENT(IN)              :: atype
      REAL(kind=dp), DIMENSION(9), INTENT(IN)            :: cell
      REAL(kind=dp), INTENT(OUT)                         :: energy
      REAL(kind=dp), DIMENSION(natom, 3), INTENT(OUT)    :: force
      REAL(kind=dp), DIMENSION(9), INTENT(OUT)           :: virial
      REAL(kind=dp), DIMENSION(natom), INTENT(OUT)       :: atomic_energy
      REAL(kind=dp), DIMENSION(natom, 9), INTENT(OUT)    :: atomic_virial

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deepmd_model_compute'

      INTEGER                                            :: handle
      INTERFACE
         SUBROUTINE DeepPotCompute(model, natom, coord, atype, cell, energy, force, virial, &
                                   atomic_energy, atomic_virial) BIND(C, name="DP_DeepPotCompute")
            IMPORT :: C_PTR, C_INT, C_DOUBLE
            TYPE(C_PTR), VALUE                                 :: model
            INTEGER(C_INT), VALUE                              :: natom
            REAL(C_DOUBLE), DIMENSION(natom, 3)                :: coord
            INTEGER(C_INT), DIMENSION(natom)                   :: atype
            REAL(C_DOUBLE), DIMENSION(9)                       :: cell
            REAL(C_DOUBLE)                                     :: energy
            REAL(C_DOUBLE), DIMENSION(natom, 3)                :: force
            REAL(C_DOUBLE), DIMENSION(9)                       :: virial
            REAL(C_DOUBLE), DIMENSION(natom)                   :: atomic_energy
            REAL(C_DOUBLE), DIMENSION(natom, 9)                :: atomic_virial
         END SUBROUTINE
      END INTERFACE

      CALL timeset(routineN, handle)

#if defined(__DEEPMD)
      CALL DeepPotCompute(model=model%c_ptr, &
                          natom=natom, &
                          coord=coord, &
                          atype=atype, &
                          cell=cell, &
                          energy=energy, &
                          force=force, &
                          virial=virial, &
                          atomic_energy=atomic_energy, &
                          atomic_virial=atomic_virial)
#else
      CPABORT("CP2K was compiled without libdeepmd_c library.")
      MARK_USED(model)
      MARK_USED(natom)
      MARK_USED(coord)
      MARK_USED(atype)
      MARK_USED(cell)
      energy = 0.0_dp
      force = 0.0_dp
      virial = 0.0_dp
      atomic_energy = 0.0_dp
      atomic_virial = 0.0_dp
#endif

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Releases a deepmd model and all its ressources.
!> \param model Pointer to the DP model.
! **************************************************************************************************
   SUBROUTINE deepmd_model_release(model)
      TYPE(deepmd_model_type)                            :: model

      model%c_ptr = C_NULL_PTR
   END SUBROUTINE deepmd_model_release

END MODULE deepmd_wrapper
