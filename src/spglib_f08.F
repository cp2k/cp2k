!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface for SPGLIB symmetry routines
!> \par History
!> \author jgh
! **************************************************************************************************
#ifdef __SPGLIB
   MODULE spglib_f08

   USE iso_c_binding,                   ONLY: c_char,&
                                              c_double,&
                                              c_int

      IMPLICIT NONE

      PRIVATE

      INTERFACE

         FUNCTION spg_get_symmetry(rotation, translation, max_size, lattice, &
                                   position, types, num_atom, symprec) bind(c)
            import c_int, c_double
            INTEGER(c_int), INTENT(inout) :: rotation(3, 3, *)
            REAL(c_double), INTENT(inout) :: translation(3, *)
            INTEGER(c_int), INTENT(in), value :: max_size
            REAL(c_double), INTENT(in) :: lattice(3, 3), position(3, *)
            INTEGER(c_int), INTENT(in) :: types(*)
            INTEGER(c_int), INTENT(in), value :: num_atom
            REAL(c_double), INTENT(in), value :: symprec
            INTEGER(c_int) :: spg_get_symmetry
         END FUNCTION spg_get_symmetry

         FUNCTION spg_get_multiplicity(lattice, position, types, num_atom, symprec) bind(c)
            import c_int, c_double
            REAL(c_double), INTENT(in) :: lattice(3, 3), position(3, *)
            INTEGER(c_int), INTENT(in) :: types(*)
            INTEGER(c_int), INTENT(in), value :: num_atom
            REAL(c_double), INTENT(in), value :: symprec
            INTEGER(c_int) :: spg_get_multiplicity
         END FUNCTION spg_get_multiplicity

         FUNCTION spg_get_international(symbol, lattice, position, types, num_atom, symprec) bind(c)
            import c_char, c_int, c_double
            CHARACTER(kind=c_char), INTENT(out) :: symbol(11)
            REAL(c_double), INTENT(in) :: lattice(3, 3), position(3, *)
            INTEGER(c_int), INTENT(in) :: types(*)
            INTEGER(c_int), INTENT(in), value :: num_atom
            REAL(c_double), INTENT(in), value :: symprec
            INTEGER(c_int) :: spg_get_international ! the number corresponding to 'symbol'. 0 on failure
         END FUNCTION spg_get_international

         FUNCTION spg_get_schoenflies(symbol, lattice, position, types, num_atom, symprec) bind(c)
            import c_char, c_int, c_double
            CHARACTER(kind=c_char), INTENT(out) :: symbol(7)
            REAL(c_double), INTENT(in) :: lattice(3, 3), position(3, *)
            INTEGER(c_int), INTENT(in) :: types(*)
            INTEGER(c_int), INTENT(in), value :: num_atom
            REAL(c_double), INTENT(in), value :: symprec
            INTEGER(c_int) :: spg_get_schoenflies ! the number corresponding to 'symbol'. 0 on failure
         END FUNCTION spg_get_schoenflies

         FUNCTION spg_get_pointgroup(symbol, trans_mat, rotations, num_rotations) bind(c)
            import c_char, c_int, c_double
            CHARACTER(kind=c_char) :: symbol(6)
            INTEGER(c_int), INTENT(inout) :: trans_mat(3, 3)
            INTEGER(c_int), INTENT(in) :: rotations(3, 3, *)
            INTEGER(c_int), INTENT(in), value :: num_rotations
            INTEGER(c_int) :: spg_get_pointgroup
         END FUNCTION spg_get_pointgroup

         FUNCTION spg_get_ir_reciprocal_mesh(grid_point, map, mesh, &
                                             is_shift, is_time_reversal, lattice, position, types, num_atom, symprec) bind(c)
            import c_int, c_double
!   Beware the map refers to positions starting at 0
            INTEGER(c_int), INTENT(out) :: grid_point(3, *), map(*) ! size is product(mesh)
            INTEGER(c_int), INTENT(in) :: mesh(3), is_shift(3)
            INTEGER(c_int), INTENT(in), value :: is_time_reversal
            REAL(c_double), INTENT(in) :: lattice(3, 3), position(3, *)
            INTEGER(c_int), INTENT(in) :: types(*)
            INTEGER(c_int), INTENT(in), value :: num_atom
            REAL(c_double), INTENT(in), value :: symprec
            INTEGER(c_int) :: spg_get_ir_reciprocal_mesh ! the number of points in the reduced mesh
         END FUNCTION spg_get_ir_reciprocal_mesh

         FUNCTION spg_get_major_version() bind(c, name='spg_get_major_version')
      INTEGER                                            :: spg_get_major_version

         END FUNCTION spg_get_major_version

         FUNCTION spg_get_minor_version() bind(c, name='spg_get_minor_version')
      INTEGER                                            :: spg_get_minor_version

         END FUNCTION spg_get_minor_version

         FUNCTION spg_get_micro_version() bind(c, name='spg_get_micro_version')
      INTEGER                                            :: spg_get_micro_version

         END FUNCTION spg_get_micro_version

      END INTERFACE

      PUBLIC :: spg_get_symmetry, spg_get_multiplicity, spg_get_international, spg_get_schoenflies, &
                spg_get_pointgroup, spg_get_ir_reciprocal_mesh, &
                spg_get_major_version, spg_get_minor_version, spg_get_micro_version

   END MODULE spglib_f08
#else
! **************************************************************************************************
!> \brief This is a stub for the Interface for SPGLIB symmetry routines
!> \par History
!> \author jgh
! **************************************************************************************************
   MODULE spglib_f08

   USE kinds, ONLY: dp

#include "./base/base_uses.f90"

      IMPLICIT NONE
      PRIVATE
      PUBLIC :: spg_get_symmetry, spg_get_multiplicity, spg_get_international, spg_get_schoenflies, &
                spg_get_pointgroup, spg_get_ir_reciprocal_mesh, &
                spg_get_major_version, spg_get_minor_version, spg_get_micro_version

   CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param rotation ...
!> \param translation ...
!> \param max_size ...
!> \param lattice ...
!> \param position ...
!> \param types ...
!> \param num_atom ...
!> \param symprec ...
!> \retval spg_get_symmetry ...
! **************************************************************************************************
      FUNCTION spg_get_symmetry(rotation, translation, max_size, lattice, &
                                position, types, num_atom, symprec)
      INTEGER, INTENT(inout)                             :: rotation(:, :, :)
      REAL(KIND=dp), INTENT(inout)                       :: translation(:, :)
      INTEGER, INTENT(in)                                :: max_size
      REAL(KIND=dp), INTENT(in)                          :: lattice(:, :), position(:, :)
      INTEGER, INTENT(in)                                :: types(:), num_atom
      REAL(KIND=dp), INTENT(in)                          :: symprec
      INTEGER                                            :: spg_get_symmetry

         MARK_USED(rotation)
         MARK_USED(translation)
         MARK_USED(max_size)
         MARK_USED(lattice)
         MARK_USED(position)
         MARK_USED(types)
         MARK_USED(num_atom)
         MARK_USED(symprec)
         spg_get_symmetry = 0
         CPABORT("Requires linking to the SPGLIB library.")
      END FUNCTION spg_get_symmetry

! **************************************************************************************************
!> \brief ...
!> \param lattice ...
!> \param position ...
!> \param types ...
!> \param num_atom ...
!> \param symprec ...
!> \retval spg_get_multiplicity ...
! **************************************************************************************************
      FUNCTION spg_get_multiplicity(lattice, position, types, num_atom, symprec)
      REAL(KIND=dp), INTENT(in)                          :: lattice(:, :), position(:, :)
      INTEGER, INTENT(in)                                :: types(:), num_atom
      REAL(KIND=dp), INTENT(in)                          :: symprec
      INTEGER                                            :: spg_get_multiplicity

         MARK_USED(lattice)
         MARK_USED(position)
         MARK_USED(types)
         MARK_USED(num_atom)
         MARK_USED(symprec)
         spg_get_multiplicity = 0
         CPABORT("Requires linking to the SPGLIB library.")
      END FUNCTION spg_get_multiplicity

! **************************************************************************************************
!> \brief ...
!> \param symbol ...
!> \param lattice ...
!> \param position ...
!> \param types ...
!> \param num_atom ...
!> \param symprec ...
!> \retval spg_get_international ...
! **************************************************************************************************
      FUNCTION spg_get_international(symbol, lattice, position, types, num_atom, symprec)
      CHARACTER, INTENT(out)                             :: symbol(11)
      REAL(KIND=dp), INTENT(in)                          :: lattice(:, :), position(:, :)
      INTEGER, INTENT(in)                                :: types(:), num_atom
      REAL(KIND=dp), INTENT(in)                          :: symprec
      INTEGER                                            :: spg_get_international

! the number corresponding to 'symbol'. 0 on failure
         MARK_USED(lattice)
         MARK_USED(position)
         MARK_USED(types)
         MARK_USED(num_atom)
         MARK_USED(symprec)
         spg_get_international = 0
         symbol = " "
         CPABORT("Requires linking to the SPGLIB library.")
      END FUNCTION spg_get_international

! **************************************************************************************************
!> \brief ...
!> \param symbol ...
!> \param lattice ...
!> \param position ...
!> \param types ...
!> \param num_atom ...
!> \param symprec ...
!> \retval spg_get_schoenflies ...
! **************************************************************************************************
      FUNCTION spg_get_schoenflies(symbol, lattice, position, types, num_atom, symprec)
      CHARACTER, INTENT(out)                             :: symbol(7)
      REAL(KIND=dp), INTENT(in)                          :: lattice(:, :), position(:, :)
      INTEGER, INTENT(in)                                :: types(:), num_atom
      REAL(KIND=dp), INTENT(in)                          :: symprec
      INTEGER                                            :: spg_get_schoenflies

! the number corresponding to 'symbol'. 0 on failure
         MARK_USED(lattice)
         MARK_USED(position)
         MARK_USED(types)
         MARK_USED(num_atom)
         MARK_USED(symprec)
         spg_get_schoenflies = 0
         symbol = " "
         CPABORT("Requires linking to the SPGLIB library.")
      END FUNCTION spg_get_schoenflies

! **************************************************************************************************
!> \brief ...
!> \param symbol ...
!> \param trans_mat ...
!> \param rotations ...
!> \param num_rotations ...
!> \retval spg_get_pointgroup ...
! **************************************************************************************************
      FUNCTION spg_get_pointgroup(symbol, trans_mat, rotations, num_rotations)
      CHARACTER                                          :: symbol(6)
      INTEGER, INTENT(inout)                             :: trans_mat(:, :)
      INTEGER, INTENT(in)                                :: rotations(:, :, :), num_rotations
      INTEGER                                            :: spg_get_pointgroup

         MARK_USED(trans_mat)
         MARK_USED(rotations)
         MARK_USED(num_rotations)
         spg_get_pointgroup = 0
         symbol = " "
         CPABORT("Requires linking to the SPGLIB library.")
      END FUNCTION spg_get_pointgroup

! **************************************************************************************************
!> \brief ...
!> \param grid_point ...
!> \param map ...
!> \param mesh ...
!> \param is_shift ...
!> \param is_time_reversal ...
!> \param lattice ...
!> \param position ...
!> \param types ...
!> \param num_atom ...
!> \param symprec ...
!> \retval spg_get_ir_reciprocal_mesh ...
! **************************************************************************************************
      FUNCTION spg_get_ir_reciprocal_mesh(grid_point, map, mesh, &
                                          is_shift, is_time_reversal, lattice, position, types, num_atom, symprec)
      INTEGER, INTENT(out)                               :: grid_point(:, :), map(:)
      INTEGER, INTENT(in)                                :: mesh(:), is_shift(:), is_time_reversal
      REAL(KIND=dp), INTENT(in)                          :: lattice(:, :), position(:, :)
      INTEGER, INTENT(in)                                :: types(:), num_atom
      REAL(KIND=dp), INTENT(in)                          :: symprec
      INTEGER                                            :: spg_get_ir_reciprocal_mesh

! size is product(mesh)

! the number of points in the reduced mesh
         MARK_USED(grid_point)
         MARK_USED(map)
         MARK_USED(mesh)
         MARK_USED(is_shift)
         MARK_USED(is_time_reversal)
         MARK_USED(lattice)
         MARK_USED(position)
         MARK_USED(types)
         MARK_USED(num_atom)
         MARK_USED(symprec)
         spg_get_ir_reciprocal_mesh = 0
         CPABORT("Requires linking to the SPGLIB library.")
      END FUNCTION spg_get_ir_reciprocal_mesh

! **************************************************************************************************
!> \brief ...
!> \retval spg_get_major_version ...
! **************************************************************************************************
      FUNCTION spg_get_major_version()
      INTEGER                                            :: spg_get_major_version

         spg_get_major_version = 0
      END FUNCTION spg_get_major_version

! **************************************************************************************************
!> \brief ...
!> \retval spg_get_minor_version ...
! **************************************************************************************************
      FUNCTION spg_get_minor_version()
      INTEGER                                            :: spg_get_minor_version

         spg_get_minor_version = 0
      END FUNCTION spg_get_minor_version

! **************************************************************************************************
!> \brief ...
!> \retval spg_get_micro_version ...
! **************************************************************************************************
      FUNCTION spg_get_micro_version()
      INTEGER                                            :: spg_get_micro_version

         spg_get_micro_version = 0
      END FUNCTION spg_get_micro_version

   END MODULE spglib_f08

#endif
