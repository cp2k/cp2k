MODULE trajana_trajectory_types
   USE trajana_kinds,                   ONLY: dp

   IMPLICIT NONE
   PRIVATE

   TYPE, PUBLIC :: cell_type
      REAL(dp) :: h(3, 3) = 0.0_dp
      LOGICAL :: periodic(3) = .TRUE.
      LOGICAL :: valid = .FALSE.
   END TYPE cell_type

   TYPE, PUBLIC :: frame_type
      INTEGER :: number = 0
      INTEGER :: natoms = 0
      CHARACTER(LEN=:), ALLOCATABLE :: comment
      CHARACTER(LEN=32), ALLOCATABLE :: label(:)
      REAL(dp), ALLOCATABLE :: value(:, :)
      TYPE(cell_type) :: cell
   END TYPE frame_type

END MODULE trajana_trajectory_types
