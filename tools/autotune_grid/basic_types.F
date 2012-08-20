MODULE basic_types

  USE kinds

  TYPE primitive_type
     INTEGER ::  la_max,la_min,lb_max,lb_min,igrid
     REAL(KIND=dp) :: zeta,zetb,ra(3),rab(3),rab2
     REAL(KIND=dp) :: intsum
  END TYPE

  TYPE realspace_grid_type
     INTEGER :: ngpts                                    ! # grid points
     INTEGER, DIMENSION (3) :: npts                      ! # grid points per dimension
     INTEGER, DIMENSION (3) :: lb                        ! lower bounds
     INTEGER, DIMENSION (3) :: ub                        ! upper bounds
     REAL(KIND=dp), DIMENSION ( :, :, : ),POINTER :: r   ! the grid
     REAL(KIND=dp), DIMENSION(3) :: dr                   ! grid spacing
     REAL(KIND=dp), DIMENSION(3,3) :: dh                 ! grid matrix
     INTEGER, DIMENSION (3) :: perd                      ! periodicity enforced
     INTEGER :: border                                   ! border points
     INTEGER :: ngpts_local                              ! local dimensions
     INTEGER, DIMENSION (3) :: npts_local
     INTEGER, DIMENSION (3) :: lb_local
     INTEGER, DIMENSION (3) :: ub_local
  END TYPE realspace_grid_type

END MODULE basic_types
