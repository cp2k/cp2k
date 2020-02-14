!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
MODULE grid_api
   USE grid_collocate,                  ONLY: &
        GRID_FUNC_AB, GRID_FUNC_ADBmDAB_X, GRID_FUNC_ADBmDAB_Y, GRID_FUNC_ADBmDAB_Z, &
        GRID_FUNC_ARDBmDARB_XX, GRID_FUNC_ARDBmDARB_XY, GRID_FUNC_ARDBmDARB_XZ, &
        GRID_FUNC_ARDBmDARB_YX, GRID_FUNC_ARDBmDARB_YY, GRID_FUNC_ARDBmDARB_YZ, &
        GRID_FUNC_ARDBmDARB_ZX, GRID_FUNC_ARDBmDARB_ZY, GRID_FUNC_ARDBmDARB_ZZ, &
        GRID_FUNC_DABpADB_X, GRID_FUNC_DABpADB_Y, GRID_FUNC_DABpADB_Z, GRID_FUNC_DADB, &
        GRID_FUNC_DX, GRID_FUNC_DXDX, GRID_FUNC_DXDY, GRID_FUNC_DY, GRID_FUNC_DYDY, &
        GRID_FUNC_DYDZ, GRID_FUNC_DZ, GRID_FUNC_DZDX, GRID_FUNC_DZDZ, collocate_pgf_product
   USE grid_integrate,                  ONLY: integrate_pgf_product

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'grid_api'

   PUBLIC :: GRID_FUNC_AB
   PUBLIC :: GRID_FUNC_DADB
   PUBLIC :: GRID_FUNC_ADBmDAB_X
   PUBLIC :: GRID_FUNC_ADBmDAB_Y
   PUBLIC :: GRID_FUNC_ADBmDAB_Z
   PUBLIC :: GRID_FUNC_ARDBmDARB_XX
   PUBLIC :: GRID_FUNC_ARDBmDARB_XY
   PUBLIC :: GRID_FUNC_ARDBmDARB_XZ
   PUBLIC :: GRID_FUNC_ARDBmDARB_YX
   PUBLIC :: GRID_FUNC_ARDBmDARB_YY
   PUBLIC :: GRID_FUNC_ARDBmDARB_YZ
   PUBLIC :: GRID_FUNC_ARDBmDARB_ZX
   PUBLIC :: GRID_FUNC_ARDBmDARB_ZY
   PUBLIC :: GRID_FUNC_ARDBmDARB_ZZ
   PUBLIC :: GRID_FUNC_DABpADB_X
   PUBLIC :: GRID_FUNC_DABpADB_Y
   PUBLIC :: GRID_FUNC_DABpADB_Z
   PUBLIC :: GRID_FUNC_DX
   PUBLIC :: GRID_FUNC_DY
   PUBLIC :: GRID_FUNC_DZ
   PUBLIC :: GRID_FUNC_DXDY
   PUBLIC :: GRID_FUNC_DYDZ
   PUBLIC :: GRID_FUNC_DZDX
   PUBLIC :: GRID_FUNC_DXDX
   PUBLIC :: GRID_FUNC_DYDY
   PUBLIC :: GRID_FUNC_DZDZ

   PUBLIC :: collocate_pgf_product
   PUBLIC :: integrate_pgf_product

CONTAINS

END MODULE grid_api
