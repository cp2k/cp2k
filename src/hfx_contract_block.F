!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
!> \brief routines to contract density matrix blocks with the for center
!>        integrals to yield the Kohn-Sham matrix. The specialized routines
!>        are about 1.2-2.0 as fast as the default one.
!> \par History
!>      10.2009 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
MODULE hfx_contract_block
  USE kinds,                           ONLY: dp
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_contract_block'
  PUBLIC :: contract_block
CONTAINS
! *****************************************************************************
!> \brief ...
!> \param ma_max ...
!> \param mb_max ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE contract_block(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: ma_max, mb_max, mc_max, md_max
    REAL(KIND=dp) :: kbd(mb_max*md_max), kbc(mb_max*mc_max), &
      kad(ma_max*md_max), kac(ma_max*mc_max), pbd(mb_max*md_max), &
      pbc(mb_max*mc_max), pad(ma_max*md_max), pac(ma_max*mc_max), &
      prim(ma_max*mb_max*mc_max*md_max), scale

#if defined (__LIBINT)
    SELECT CASE(ma_max)
    CASE(1)
      SELECT CASE(mb_max)
      CASE(1)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_1_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_1_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_1_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_1_1_4(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_1_1_5(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_1_1_6(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_1_1_7(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_1_1_9(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_1_1_10(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_1_1_11(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_1_1_15(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_1_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_1_2_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_1_2_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_1_2_4(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_1_2_5(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_1_2_6(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_1_2_7(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_1_2_9(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(3)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_1_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_1_3_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_1_3_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_1_3_4(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_1_3_5(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_1_3_6(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(4)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_1_4_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_1_4_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_1_4_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_1_4_4(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(5)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_1_5_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_1_5_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_1_5_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(6)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_1_6_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_1_6_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_1_6_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(7)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_1_7_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_1_7_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_1_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_1_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_1_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_1_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_1_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_1_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_1_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_1_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_1_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(9)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_1_9_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_1_9_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_1_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_1_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_1_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_1_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_1_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_1_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_1_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_1_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_1_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(10)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_1_10_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_1_10(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_1_10(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_1_10(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_1_10(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_1_10(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_1_10(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_1_10(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_1_10(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_1_10(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_1_10(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(11)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_1_11_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_1_11(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_1_11(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_1_11(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_1_11(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_1_11(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_1_11(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_1_11(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_1_11(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_1_11(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_1_11(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(15)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_1_15_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_1_15(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_1_15(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_1_15(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_1_15(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_1_15(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_1_15(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_1_15(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_1_15(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_1_15(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_1_15(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(2)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_2_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_2_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_2_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_2_1_4(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_2_1_5(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_2_1_6(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_2_1_7(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_2_1_9(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_2_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_2_2_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_2_2_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_2_2_4(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(3)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_2_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_2_3_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_2_3_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(4)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_2_4_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_2_4_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_2_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_2_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_2_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_2_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_2_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_2_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_2_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_2_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_2_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(5)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_2_5_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_2_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_2_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_2_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_2_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_2_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_2_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_2_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_2_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_2_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_2_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(6)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_2_6_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_2_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_2_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_2_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_2_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_2_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_2_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_2_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_2_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_2_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_2_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(7)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_2_7_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_2_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_2_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_2_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_2_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_2_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_2_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_2_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_2_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_2_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_2_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(9)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_2_9_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_2_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_2_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_2_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_2_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_2_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_2_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_2_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_2_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_2_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_2_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(10)
          CALL block_1_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_1_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_1_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(3)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_3_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_3_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_3_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_3_1_4(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_3_1_5(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_3_1_6(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_3_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_3_2_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_3_2_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(3)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_3_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_3_3_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_3_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_3_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_3_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_3_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_3_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_3_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_3_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_3_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_3_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(4)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_3_4_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_3_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_3_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_3_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_3_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_3_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_3_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_3_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_3_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_3_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_3_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(5)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_3_5_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_3_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_3_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_3_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_3_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_3_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_3_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_3_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_3_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_3_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_3_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(6)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_3_6_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_3_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_3_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_3_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_3_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_3_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_3_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_3_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_3_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_3_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_3_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(7)
          CALL block_1_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_1_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_1_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_1_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_1_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(4)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_4_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_4_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_4_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_4_1_4(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_4_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_4_2_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_4_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_4_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_4_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_4_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_4_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_4_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_4_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_4_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_4_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(3)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_4_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_4_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_4_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_4_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_4_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_4_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_4_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_4_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_4_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_4_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_4_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(4)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_4_4_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_4_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_4_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_4_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_4_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_4_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_4_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_4_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_4_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_4_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_4_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(5)
          CALL block_1_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_1_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_1_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_1_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_1_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_1_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_1_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(5)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_5_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_5_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_5_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_5_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_5_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_5_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_5_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_5_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_5_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_5_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_5_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_5_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_5_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_5_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(3)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_5_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_5_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_5_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_5_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_5_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_5_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_5_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_5_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_5_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_5_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_5_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(4)
          CALL block_1_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_1_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_1_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_1_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_1_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_1_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_1_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_1_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(6)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_6_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_6_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_6_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_6_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_6_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_6_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_6_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_6_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_6_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_6_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_6_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_6_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_6_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_6_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(3)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_6_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_6_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_6_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_6_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_6_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_6_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_6_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_6_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_6_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_6_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_6_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(4)
          CALL block_1_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_1_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_1_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_1_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_1_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_1_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_1_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_1_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(7)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_7_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_7_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_7_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_7_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_7_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_7_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_7_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_7_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_7_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_7_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_7_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_7_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_7_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_7_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_7_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_7_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_7_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_7_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_7_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_7_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_7_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_7_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(3)
          CALL block_1_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_1_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_1_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_1_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_1_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_1_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_1_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_1_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_1_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(9)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_9_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_9_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_9_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_9_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_9_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_9_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_9_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_9_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_9_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_9_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_9_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_9_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_9_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_9_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_9_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_9_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_9_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_9_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_9_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_9_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_9_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_9_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(3)
          CALL block_1_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_1_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_1_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_1_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_1_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_1_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_1_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_1_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_1_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(10)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_10_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_10_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_10_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_10_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_10_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_10_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_10_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_10_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_10_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_10_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_10_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          CALL block_1_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(3)
          CALL block_1_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_1_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_1_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_1_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_1_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_1_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_1_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_1_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_1_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(11)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_11_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_11_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_11_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_11_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_11_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_11_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_11_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_11_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_11_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_11_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_11_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          CALL block_1_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(3)
          CALL block_1_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_1_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_1_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_1_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_1_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_1_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_1_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_1_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_1_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(15)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_1_15_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_1_15_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_1_15_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_1_15_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_1_15_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_1_15_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_1_15_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_1_15_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_1_15_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_1_15_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_1_15_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          CALL block_1_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(3)
          CALL block_1_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_1_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_1_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_1_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_1_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_1_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_1_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_1_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_1_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE DEFAULT
        CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      END SELECT
    CASE(2)
      SELECT CASE(mb_max)
      CASE(1)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_2_1_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_2_1_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_2_1_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_2_1_1_4(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_2_1_1_5(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_2_1_1_6(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_2_1_1_7(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_2_1_1_9(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_2_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_2_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_2_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_2_1_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_2_1_2_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_2_1_2_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_2_1_2_4(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_2_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_2_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_2_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_2_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_2_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_2_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_2_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(3)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_2_1_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_2_1_3_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_2_1_3_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_2_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_2_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_2_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_2_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_2_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_2_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_2_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_2_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(4)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_2_1_4_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_2_1_4_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_2_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_2_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_2_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_2_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_2_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_2_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_2_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_2_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_2_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(5)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_2_1_5_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_2_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_2_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_2_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_2_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_2_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_2_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_2_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_2_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_2_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_2_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(6)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_2_1_6_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_2_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_2_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_2_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_2_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_2_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_2_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_2_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_2_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_2_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_2_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(7)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_2_1_7_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_2_1_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_2_1_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_2_1_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_2_1_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_2_1_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_2_1_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_2_1_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_2_1_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_2_1_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_2_1_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(9)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_2_1_9_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_2_1_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_2_1_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_2_1_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_2_1_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_2_1_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_2_1_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_2_1_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_2_1_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_2_1_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_2_1_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(10)
          CALL block_2_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_2_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_2_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(2)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_2_2_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_2_2_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_2_2_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_2_2_1_4(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_2_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_2_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_2_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_2_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_2_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_2_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_2_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_2_2_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_2_2_2_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_2_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_2_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_2_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_2_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_2_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_2_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_2_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_2_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_2_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(3)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_2_2_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_2_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_2_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_2_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_2_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_2_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_2_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_2_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_2_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_2_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_2_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(4)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_2_2_4_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_2_2_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_2_2_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_2_2_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_2_2_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_2_2_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_2_2_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_2_2_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_2_2_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_2_2_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_2_2_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(5)
          CALL block_2_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_2_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_2_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_2_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_2_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_2_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_2_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(3)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_2_3_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_2_3_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_2_3_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_2_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_2_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_2_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_2_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_2_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_2_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_2_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_2_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_2_3_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_2_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_2_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_2_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_2_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_2_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_2_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_2_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_2_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_2_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_2_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(3)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_2_3_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_2_3_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_2_3_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_2_3_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_2_3_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_2_3_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_2_3_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_2_3_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_2_3_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_2_3_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_2_3_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(4)
          CALL block_2_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_2_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_2_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_2_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_2_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_2_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_2_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_2_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(4)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_2_4_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_2_4_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_2_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_2_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_2_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_2_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_2_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_2_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_2_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_2_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_2_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_2_4_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_2_4_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_2_4_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_2_4_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_2_4_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_2_4_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_2_4_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_2_4_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_2_4_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_2_4_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_2_4_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(3)
          CALL block_2_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_2_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_2_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_2_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_2_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_2_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_2_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_2_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_2_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(5)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_2_5_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_2_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_2_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_2_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_2_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_2_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_2_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_2_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_2_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_2_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_2_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          CALL block_2_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(3)
          CALL block_2_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_2_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_2_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_2_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_2_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_2_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_2_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_2_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_2_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(6)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_2_6_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_2_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_2_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_2_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_2_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_2_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_2_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_2_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_2_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_2_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_2_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          CALL block_2_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(3)
          CALL block_2_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_2_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_2_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_2_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_2_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_2_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_2_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_2_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_2_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(7)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_2_7_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_2_7_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_2_7_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_2_7_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_2_7_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_2_7_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_2_7_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_2_7_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_2_7_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_2_7_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_2_7_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          CALL block_2_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(3)
          CALL block_2_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_2_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_2_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_2_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_2_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_2_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_2_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_2_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_2_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(9)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_2_9_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_2_9_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_2_9_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_2_9_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_2_9_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_2_9_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_2_9_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_2_9_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_2_9_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_2_9_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_2_9_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          CALL block_2_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(3)
          CALL block_2_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_2_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_2_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_2_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_2_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_2_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_2_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_2_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_2_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(10)
        CALL block_2_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(11)
        CALL block_2_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(15)
        CALL block_2_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE DEFAULT
        CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      END SELECT
    CASE(3)
      SELECT CASE(mb_max)
      CASE(1)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_3_1_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_3_1_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_3_1_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_3_1_1_4(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_3_1_1_5(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_3_1_1_6(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_3_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_3_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_3_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_3_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_3_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_3_1_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_3_1_2_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_3_1_2_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_3_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_3_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_3_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_3_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_3_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_3_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_3_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_3_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(3)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_3_1_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_3_1_3_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_3_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_3_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_3_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_3_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_3_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_3_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_3_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_3_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_3_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(4)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_3_1_4_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_3_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_3_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_3_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_3_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_3_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_3_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_3_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_3_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_3_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_3_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(5)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_3_1_5_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_3_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_3_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_3_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_3_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_3_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_3_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_3_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_3_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_3_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_3_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(6)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_3_1_6_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_3_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_3_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_3_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_3_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_3_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_3_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_3_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_3_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_3_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_3_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(7)
          CALL block_3_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_3_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_3_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_3_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_3_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(2)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_3_2_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_3_2_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_3_2_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_3_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_3_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_3_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_3_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_3_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_3_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_3_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_3_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_3_2_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_3_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_3_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_3_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_3_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_3_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_3_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_3_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_3_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_3_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_3_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(3)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_3_2_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_3_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_3_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_3_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_3_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_3_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_3_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_3_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_3_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_3_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_3_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(4)
          CALL block_3_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_3_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_3_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_3_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_3_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_3_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_3_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_3_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(3)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_3_3_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_3_3_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_3_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_3_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_3_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_3_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_3_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_3_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_3_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_3_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_3_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_3_3_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_3_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_3_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_3_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_3_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_3_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_3_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_3_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_3_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_3_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_3_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(3)
          CALL block_3_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_3_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_3_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_3_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_3_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_3_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_3_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_3_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_3_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(4)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_3_4_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_3_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_3_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_3_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_3_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_3_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_3_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_3_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_3_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_3_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_3_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          CALL block_3_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(3)
          CALL block_3_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_3_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_3_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_3_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_3_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_3_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_3_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_3_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_3_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(5)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_3_5_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_3_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_3_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_3_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_3_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_3_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_3_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_3_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_3_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_3_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_3_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          CALL block_3_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(3)
          CALL block_3_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_3_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_3_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_3_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_3_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_3_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_3_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_3_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_3_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(6)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_3_6_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_3_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_3_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_3_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_3_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_3_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_3_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_3_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_3_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_3_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_3_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          CALL block_3_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(3)
          CALL block_3_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_3_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_3_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_3_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_3_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_3_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_3_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_3_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_3_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(7)
        CALL block_3_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(9)
        CALL block_3_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(10)
        CALL block_3_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(11)
        CALL block_3_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(15)
        CALL block_3_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE DEFAULT
        CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      END SELECT
    CASE(4)
      SELECT CASE(mb_max)
      CASE(1)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_4_1_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_4_1_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_4_1_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_4_1_1_4(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_4_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_4_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_4_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_4_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_4_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_4_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_4_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_4_1_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_4_1_2_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_4_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_4_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_4_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_4_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_4_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_4_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_4_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_4_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_4_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(3)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_4_1_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_4_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_4_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_4_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_4_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_4_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_4_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_4_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_4_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_4_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_4_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(4)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_4_1_4_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_4_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_4_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_4_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_4_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_4_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_4_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_4_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_4_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_4_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_4_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(5)
          CALL block_4_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_4_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_4_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_4_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_4_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_4_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_4_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(2)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_4_2_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_4_2_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_4_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_4_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_4_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_4_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_4_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_4_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_4_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_4_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_4_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_4_2_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_4_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_4_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_4_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_4_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_4_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_4_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_4_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_4_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_4_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_4_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(3)
          CALL block_4_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_4_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_4_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_4_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_4_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_4_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_4_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_4_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_4_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(3)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_4_3_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_4_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_4_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_4_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_4_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_4_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_4_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_4_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_4_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_4_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_4_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          CALL block_4_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(3)
          CALL block_4_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_4_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_4_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_4_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_4_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_4_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_4_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_4_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_4_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(4)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_4_4_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_4_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_4_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_4_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_4_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_4_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_4_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_4_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_4_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_4_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_4_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          CALL block_4_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(3)
          CALL block_4_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_4_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_4_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_4_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_4_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_4_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_4_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_4_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_4_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(5)
        CALL block_4_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(6)
        CALL block_4_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(7)
        CALL block_4_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(9)
        CALL block_4_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(10)
        CALL block_4_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(11)
        CALL block_4_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(15)
        CALL block_4_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE DEFAULT
        CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      END SELECT
    CASE(5)
      SELECT CASE(mb_max)
      CASE(1)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_5_1_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_5_1_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_5_1_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_5_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_5_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_5_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_5_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_5_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_5_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_5_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_5_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_5_1_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_5_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_5_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_5_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_5_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_5_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_5_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_5_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_5_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_5_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_5_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(3)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_5_1_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_5_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_5_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_5_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_5_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_5_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_5_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_5_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_5_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_5_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_5_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(4)
          CALL block_5_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_5_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_5_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_5_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_5_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_5_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_5_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_5_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(2)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_5_2_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_5_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_5_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_5_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_5_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_5_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_5_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_5_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_5_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_5_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_5_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          CALL block_5_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(3)
          CALL block_5_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_5_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_5_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_5_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_5_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_5_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_5_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_5_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_5_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(3)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_5_3_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_5_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_5_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_5_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_5_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_5_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_5_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_5_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_5_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_5_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_5_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          CALL block_5_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(3)
          CALL block_5_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_5_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_5_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_5_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_5_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_5_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_5_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_5_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_5_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(4)
        CALL block_5_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(5)
        CALL block_5_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(6)
        CALL block_5_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(7)
        CALL block_5_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(9)
        CALL block_5_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(10)
        CALL block_5_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(11)
        CALL block_5_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(15)
        CALL block_5_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE DEFAULT
        CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      END SELECT
    CASE(6)
      SELECT CASE(mb_max)
      CASE(1)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_6_1_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_6_1_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_6_1_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_6_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_6_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_6_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_6_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_6_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_6_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_6_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_6_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_6_1_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_6_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_6_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_6_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_6_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_6_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_6_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_6_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_6_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_6_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_6_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(3)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_6_1_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_6_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_6_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_6_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_6_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_6_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_6_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_6_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_6_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_6_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_6_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(4)
          CALL block_6_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_6_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_6_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_6_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_6_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_6_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_6_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_6_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(2)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_6_2_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_6_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_6_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_6_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_6_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_6_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_6_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_6_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_6_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_6_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_6_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          CALL block_6_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(3)
          CALL block_6_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_6_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_6_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_6_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_6_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_6_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_6_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_6_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_6_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(3)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_6_3_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_6_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_6_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_6_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_6_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_6_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_6_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_6_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_6_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_6_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_6_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          CALL block_6_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(3)
          CALL block_6_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_6_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_6_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_6_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_6_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_6_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_6_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_6_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_6_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(4)
        CALL block_6_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(5)
        CALL block_6_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(6)
        CALL block_6_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(7)
        CALL block_6_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(9)
        CALL block_6_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(10)
        CALL block_6_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(11)
        CALL block_6_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(15)
        CALL block_6_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE DEFAULT
        CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      END SELECT
    CASE(7)
      SELECT CASE(mb_max)
      CASE(1)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_7_1_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_7_1_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_7_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_7_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_7_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_7_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_7_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_7_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_7_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_7_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_7_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_7_1_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_7_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_7_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_7_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_7_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_7_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_7_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_7_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_7_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_7_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_7_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(3)
          CALL block_7_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_7_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_7_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_7_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_7_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_7_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_7_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_7_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_7_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(2)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_7_2_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_7_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_7_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_7_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_7_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_7_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_7_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_7_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_7_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_7_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_7_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          CALL block_7_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(3)
          CALL block_7_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_7_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_7_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_7_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_7_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_7_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_7_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_7_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_7_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(3)
        CALL block_7_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(4)
        CALL block_7_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(5)
        CALL block_7_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(6)
        CALL block_7_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(7)
        CALL block_7_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(9)
        CALL block_7_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(10)
        CALL block_7_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(11)
        CALL block_7_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(15)
        CALL block_7_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE DEFAULT
        CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      END SELECT
    CASE(9)
      SELECT CASE(mb_max)
      CASE(1)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_9_1_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_9_1_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_9_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_9_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_9_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_9_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_9_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_9_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_9_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_9_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_9_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_9_1_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_9_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_9_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_9_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_9_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_9_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_9_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_9_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_9_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_9_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_9_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(3)
          CALL block_9_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_9_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_9_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_9_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_9_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_9_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_9_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_9_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_9_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(2)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_9_2_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_9_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_9_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_9_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_9_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_9_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_9_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_9_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_9_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_9_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_9_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          CALL block_9_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(3)
          CALL block_9_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_9_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_9_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_9_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_9_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_9_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_9_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_9_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_9_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(3)
        CALL block_9_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(4)
        CALL block_9_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(5)
        CALL block_9_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(6)
        CALL block_9_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(7)
        CALL block_9_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(9)
        CALL block_9_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(10)
        CALL block_9_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(11)
        CALL block_9_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(15)
        CALL block_9_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE DEFAULT
        CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      END SELECT
    CASE(10)
      SELECT CASE(mb_max)
      CASE(1)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_10_1_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_10_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_10_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_10_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_10_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_10_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_10_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_10_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_10_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_10_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_10_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          CALL block_10_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(3)
          CALL block_10_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_10_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_10_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_10_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_10_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_10_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_10_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_10_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_10_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(2)
        CALL block_10_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(3)
        CALL block_10_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(4)
        CALL block_10_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(5)
        CALL block_10_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(6)
        CALL block_10_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(7)
        CALL block_10_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(9)
        CALL block_10_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(10)
        CALL block_10_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(11)
        CALL block_10_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(15)
        CALL block_10_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE DEFAULT
        CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      END SELECT
    CASE(11)
      SELECT CASE(mb_max)
      CASE(1)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_11_1_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_11_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_11_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_11_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_11_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_11_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_11_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_11_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_11_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_11_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_11_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          CALL block_11_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(3)
          CALL block_11_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_11_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_11_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_11_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_11_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_11_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_11_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_11_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_11_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(2)
        CALL block_11_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(3)
        CALL block_11_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(4)
        CALL block_11_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(5)
        CALL block_11_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(6)
        CALL block_11_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(7)
        CALL block_11_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(9)
        CALL block_11_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(10)
        CALL block_11_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(11)
        CALL block_11_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(15)
        CALL block_11_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE DEFAULT
        CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      END SELECT
    CASE(15)
      SELECT CASE(mb_max)
      CASE(1)
        SELECT CASE(mc_max)
        CASE(1)
          SELECT CASE(md_max)
          CASE(1)
            CALL block_15_1_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(2)
            CALL block_15_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(3)
            CALL block_15_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(4)
            CALL block_15_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(5)
            CALL block_15_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(6)
            CALL block_15_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(7)
            CALL block_15_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(9)
            CALL block_15_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(10)
            CALL block_15_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(11)
            CALL block_15_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE(15)
            CALL block_15_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          CASE DEFAULT
            CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
          END SELECT
        CASE(2)
          CALL block_15_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(3)
          CALL block_15_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(4)
          CALL block_15_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(5)
          CALL block_15_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(6)
          CALL block_15_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(7)
          CALL block_15_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(9)
          CALL block_15_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(10)
          CALL block_15_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(11)
          CALL block_15_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE(15)
          CALL block_15_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        CASE DEFAULT
          CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
        END SELECT
      CASE(2)
        CALL block_15_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(3)
        CALL block_15_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(4)
        CALL block_15_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(5)
        CALL block_15_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(6)
        CALL block_15_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(7)
        CALL block_15_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(9)
        CALL block_15_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(10)
        CALL block_15_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(11)
        CALL block_15_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE(15)
        CALL block_15_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      CASE DEFAULT
        CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
      END SELECT
    CASE DEFAULT
      CALL block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    END SELECT
#else
    RETURN
#endif
  END SUBROUTINE contract_block

#if defined (__LIBINT)
! *****************************************************************************
!> \brief ...
!> \param ma_max ...
!> \param mb_max ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_default(ma_max,mb_max,mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: ma_max, mb_max, mc_max, md_max
    REAL(KIND=dp) :: kbd(mb_max*md_max), kbc(mb_max*mc_max), &
      kad(ma_max*md_max), kac(ma_max*mc_max), pbd(mb_max*md_max), &
      pbc(mb_max*mc_max), pad(ma_max*md_max), pac(ma_max*mc_max), &
      prim(ma_max*mb_max*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:mb_max*md_max) = 0.0_dp
      kbc(1:mb_max*mc_max) = 0.0_dp
      kad(1:ma_max*md_max) = 0.0_dp
      kac(1:ma_max*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,mb_max
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*mb_max+mb)
            p_bc = pbc((mc-1)*mb_max+mb)
            DO ma = 1,ma_max
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*ma_max+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*ma_max+ma)
              kad((md-1)*ma_max+ma) =  kad((md-1)*ma_max+ma)-tmp*p_bc
              kac((mc-1)*ma_max+ma) = kac((mc-1)*ma_max+ma)-tmp*p_bd
            END DO
            kbd((md-1)*mb_max+mb) = kbd((md-1)*mb_max+mb) - ks_bd
            kbc((mc-1)*mb_max+mb) = kbc((mc-1)*mb_max+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_default
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*1), kad(1*1), kac(1*1), pbd(1*1), &
      pbc(1*1), pad(1*1), pac(1*1), prim(1*1*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*2), kbc(1*1), kad(1*2), kac(1*1), pbd(1*2), &
      pbc(1*1), pad(1*2), pac(1*1), prim(1*1*1*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*2) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:1*2) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_1_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*3), kbc(1*1), kad(1*3), kac(1*1), pbd(1*3), &
      pbc(1*1), pad(1*3), pac(1*1), prim(1*1*1*3), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*3) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:1*3) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,3
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_1_3
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_1_4(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*4), kbc(1*1), kad(1*4), kac(1*1), pbd(1*4), &
      pbc(1*1), pad(1*4), pac(1*1), prim(1*1*1*4), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*4) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:1*4) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,4
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_1_4
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_1_5(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*5), kbc(1*1), kad(1*5), kac(1*1), pbd(1*5), &
      pbc(1*1), pad(1*5), pac(1*1), prim(1*1*1*5), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*5) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:1*5) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,5
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_1_5
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_1_6(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*6), kbc(1*1), kad(1*6), kac(1*1), pbd(1*6), &
      pbc(1*1), pad(1*6), pac(1*1), prim(1*1*1*6), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*6) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:1*6) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,6
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_1_6
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_1_7(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*7), kbc(1*1), kad(1*7), kac(1*1), pbd(1*7), &
      pbc(1*1), pad(1*7), pac(1*1), prim(1*1*1*7), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*7) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:1*7) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,7
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_1_7
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_1_9(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*9), kbc(1*1), kad(1*9), kac(1*1), pbd(1*9), &
      pbc(1*1), pad(1*9), pac(1*1), prim(1*1*1*9), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*9) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:1*9) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,9
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_1_9
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_1_10(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*10), kbc(1*1), kad(1*10), kac(1*1), pbd(1*10), &
      pbc(1*1), pad(1*10), pac(1*1), prim(1*1*1*10), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*10) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:1*10) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,10
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_1_10
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_1_11(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*11), kbc(1*1), kad(1*11), kac(1*1), pbd(1*11), &
      pbc(1*1), pad(1*11), pac(1*1), prim(1*1*1*11), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*11) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:1*11) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,11
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_1_11
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_1_15(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*15), kbc(1*1), kad(1*15), kac(1*1), pbd(1*15), &
      pbc(1*1), pad(1*15), pac(1*1), prim(1*1*1*15), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*15) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:1*15) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,15
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_1_15
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*2), kad(1*1), kac(1*2), pbd(1*1), &
      pbc(1*2), pad(1*1), pac(1*2), prim(1*1*2*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_2_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_2_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*2), kbc(1*2), kad(1*2), kac(1*2), pbd(1*2), &
      pbc(1*2), pad(1*2), pac(1*2), prim(1*1*2*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*2) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:1*2) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_2_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_2_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*3), kbc(1*2), kad(1*3), kac(1*2), pbd(1*3), &
      pbc(1*2), pad(1*3), pac(1*2), prim(1*1*2*3), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*3) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:1*3) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,3
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_2_3
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_2_4(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*4), kbc(1*2), kad(1*4), kac(1*2), pbd(1*4), &
      pbc(1*2), pad(1*4), pac(1*2), prim(1*1*2*4), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*4) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:1*4) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,4
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_2_4
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_2_5(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*5), kbc(1*2), kad(1*5), kac(1*2), pbd(1*5), &
      pbc(1*2), pad(1*5), pac(1*2), prim(1*1*2*5), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*5) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:1*5) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,5
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_2_5
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_2_6(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*6), kbc(1*2), kad(1*6), kac(1*2), pbd(1*6), &
      pbc(1*2), pad(1*6), pac(1*2), prim(1*1*2*6), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*6) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:1*6) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,6
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_2_6
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_2_7(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*7), kbc(1*2), kad(1*7), kac(1*2), pbd(1*7), &
      pbc(1*2), pad(1*7), pac(1*2), prim(1*1*2*7), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*7) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:1*7) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,7
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_2_7
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_2_9(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*9), kbc(1*2), kad(1*9), kac(1*2), pbd(1*9), &
      pbc(1*2), pad(1*9), pac(1*2), prim(1*1*2*9), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*9) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:1*9) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,9
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_2_9
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*2), kad(1*md_max), kac(1*2), &
      pbd(1*md_max), pbc(1*2), pad(1*md_max), pac(1*2), prim(1*1*2*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*3), kad(1*1), kac(1*3), pbd(1*1), &
      pbc(1*3), pad(1*1), pac(1*3), prim(1*1*3*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*3) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*3) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,3
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_3_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_3_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*2), kbc(1*3), kad(1*2), kac(1*3), pbd(1*2), &
      pbc(1*3), pad(1*2), pac(1*3), prim(1*1*3*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*2) = 0.0_dp
      kbc(1:1*3) = 0.0_dp
      kad(1:1*2) = 0.0_dp
      kac(1:1*3) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,3
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_3_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_3_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*3), kbc(1*3), kad(1*3), kac(1*3), pbd(1*3), &
      pbc(1*3), pad(1*3), pac(1*3), prim(1*1*3*3), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*3) = 0.0_dp
      kbc(1:1*3) = 0.0_dp
      kad(1:1*3) = 0.0_dp
      kac(1:1*3) = 0.0_dp
      p_index = 0
      DO md = 1,3
        DO mc = 1,3
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_3_3
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_3_4(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*4), kbc(1*3), kad(1*4), kac(1*3), pbd(1*4), &
      pbc(1*3), pad(1*4), pac(1*3), prim(1*1*3*4), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*4) = 0.0_dp
      kbc(1:1*3) = 0.0_dp
      kad(1:1*4) = 0.0_dp
      kac(1:1*3) = 0.0_dp
      p_index = 0
      DO md = 1,4
        DO mc = 1,3
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_3_4
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_3_5(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*5), kbc(1*3), kad(1*5), kac(1*3), pbd(1*5), &
      pbc(1*3), pad(1*5), pac(1*3), prim(1*1*3*5), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*5) = 0.0_dp
      kbc(1:1*3) = 0.0_dp
      kad(1:1*5) = 0.0_dp
      kac(1:1*3) = 0.0_dp
      p_index = 0
      DO md = 1,5
        DO mc = 1,3
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_3_5
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_3_6(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*6), kbc(1*3), kad(1*6), kac(1*3), pbd(1*6), &
      pbc(1*3), pad(1*6), pac(1*3), prim(1*1*3*6), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*6) = 0.0_dp
      kbc(1:1*3) = 0.0_dp
      kad(1:1*6) = 0.0_dp
      kac(1:1*3) = 0.0_dp
      p_index = 0
      DO md = 1,6
        DO mc = 1,3
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_3_6
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*3), kad(1*md_max), kac(1*3), &
      pbd(1*md_max), pbc(1*3), pad(1*md_max), pac(1*3), prim(1*1*3*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*3) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*3) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,3
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_3
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_4_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*4), kad(1*1), kac(1*4), pbd(1*1), &
      pbc(1*4), pad(1*1), pac(1*4), prim(1*1*4*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*4) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*4) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,4
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_4_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_4_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*2), kbc(1*4), kad(1*2), kac(1*4), pbd(1*2), &
      pbc(1*4), pad(1*2), pac(1*4), prim(1*1*4*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*2) = 0.0_dp
      kbc(1:1*4) = 0.0_dp
      kad(1:1*2) = 0.0_dp
      kac(1:1*4) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,4
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_4_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_4_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*3), kbc(1*4), kad(1*3), kac(1*4), pbd(1*3), &
      pbc(1*4), pad(1*3), pac(1*4), prim(1*1*4*3), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*3) = 0.0_dp
      kbc(1:1*4) = 0.0_dp
      kad(1:1*3) = 0.0_dp
      kac(1:1*4) = 0.0_dp
      p_index = 0
      DO md = 1,3
        DO mc = 1,4
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_4_3
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_4_4(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*4), kbc(1*4), kad(1*4), kac(1*4), pbd(1*4), &
      pbc(1*4), pad(1*4), pac(1*4), prim(1*1*4*4), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*4) = 0.0_dp
      kbc(1:1*4) = 0.0_dp
      kad(1:1*4) = 0.0_dp
      kac(1:1*4) = 0.0_dp
      p_index = 0
      DO md = 1,4
        DO mc = 1,4
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_4_4
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*4), kad(1*md_max), kac(1*4), &
      pbd(1*md_max), pbc(1*4), pad(1*md_max), pac(1*4), prim(1*1*4*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*4) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*4) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,4
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_4
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_5_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*5), kad(1*1), kac(1*5), pbd(1*1), &
      pbc(1*5), pad(1*1), pac(1*5), prim(1*1*5*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*5) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*5) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,5
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_5_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_5_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*2), kbc(1*5), kad(1*2), kac(1*5), pbd(1*2), &
      pbc(1*5), pad(1*2), pac(1*5), prim(1*1*5*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*2) = 0.0_dp
      kbc(1:1*5) = 0.0_dp
      kad(1:1*2) = 0.0_dp
      kac(1:1*5) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,5
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_5_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_5_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*3), kbc(1*5), kad(1*3), kac(1*5), pbd(1*3), &
      pbc(1*5), pad(1*3), pac(1*5), prim(1*1*5*3), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*3) = 0.0_dp
      kbc(1:1*5) = 0.0_dp
      kad(1:1*3) = 0.0_dp
      kac(1:1*5) = 0.0_dp
      p_index = 0
      DO md = 1,3
        DO mc = 1,5
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_5_3
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*5), kad(1*md_max), kac(1*5), &
      pbd(1*md_max), pbc(1*5), pad(1*md_max), pac(1*5), prim(1*1*5*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*5) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*5) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,5
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_5
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_6_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*6), kad(1*1), kac(1*6), pbd(1*1), &
      pbc(1*6), pad(1*1), pac(1*6), prim(1*1*6*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*6) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*6) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,6
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_6_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_6_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*2), kbc(1*6), kad(1*2), kac(1*6), pbd(1*2), &
      pbc(1*6), pad(1*2), pac(1*6), prim(1*1*6*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*2) = 0.0_dp
      kbc(1:1*6) = 0.0_dp
      kad(1:1*2) = 0.0_dp
      kac(1:1*6) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,6
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_6_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_6_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*3), kbc(1*6), kad(1*3), kac(1*6), pbd(1*3), &
      pbc(1*6), pad(1*3), pac(1*6), prim(1*1*6*3), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*3) = 0.0_dp
      kbc(1:1*6) = 0.0_dp
      kad(1:1*3) = 0.0_dp
      kac(1:1*6) = 0.0_dp
      p_index = 0
      DO md = 1,3
        DO mc = 1,6
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_6_3
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*6), kad(1*md_max), kac(1*6), &
      pbd(1*md_max), pbc(1*6), pad(1*md_max), pac(1*6), prim(1*1*6*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*6) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*6) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,6
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_6
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_7_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*7), kad(1*1), kac(1*7), pbd(1*1), &
      pbc(1*7), pad(1*1), pac(1*7), prim(1*1*7*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*7) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*7) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,7
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_7_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_7_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*2), kbc(1*7), kad(1*2), kac(1*7), pbd(1*2), &
      pbc(1*7), pad(1*2), pac(1*7), prim(1*1*7*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*2) = 0.0_dp
      kbc(1:1*7) = 0.0_dp
      kad(1:1*2) = 0.0_dp
      kac(1:1*7) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,7
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_7_2
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*7), kad(1*md_max), kac(1*7), &
      pbd(1*md_max), pbc(1*7), pad(1*md_max), pac(1*7), prim(1*1*7*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*7) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*7) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,7
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_7
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_9_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*9), kad(1*1), kac(1*9), pbd(1*1), &
      pbc(1*9), pad(1*1), pac(1*9), prim(1*1*9*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*9) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*9) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,9
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_9_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_9_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*2), kbc(1*9), kad(1*2), kac(1*9), pbd(1*2), &
      pbc(1*9), pad(1*2), pac(1*9), prim(1*1*9*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*2) = 0.0_dp
      kbc(1:1*9) = 0.0_dp
      kad(1:1*2) = 0.0_dp
      kac(1:1*9) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,9
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_9_2
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*9), kad(1*md_max), kac(1*9), &
      pbd(1*md_max), pbc(1*9), pad(1*md_max), pac(1*9), prim(1*1*9*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*9) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*9) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,9
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_9
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_10_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*10), kad(1*1), kac(1*10), pbd(1*1), &
      pbc(1*10), pad(1*1), pac(1*10), prim(1*1*10*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*10) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*10) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,10
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_10_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_10(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*10), kad(1*md_max), kac(1*10), &
      pbd(1*md_max), pbc(1*10), pad(1*md_max), pac(1*10), &
      prim(1*1*10*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*10) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*10) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,10
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_10
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_11_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*11), kad(1*1), kac(1*11), pbd(1*1), &
      pbc(1*11), pad(1*1), pac(1*11), prim(1*1*11*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*11) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*11) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,11
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_11_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_11(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*11), kad(1*md_max), kac(1*11), &
      pbd(1*md_max), pbc(1*11), pad(1*md_max), pac(1*11), &
      prim(1*1*11*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*11) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*11) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,11
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_11
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_15_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*15), kad(1*1), kac(1*15), pbd(1*1), &
      pbc(1*15), pad(1*1), pac(1*15), prim(1*1*15*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*15) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*15) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,15
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_15_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_1_15(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*15), kad(1*md_max), kac(1*15), &
      pbd(1*md_max), pbc(1*15), pad(1*md_max), pac(1*15), &
      prim(1*1*15*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*15) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*15) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,15
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_1_15
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*1), kbc(2*1), kad(1*1), kac(1*1), pbd(2*1), &
      pbc(2*1), pad(1*1), pac(1*1), prim(1*2*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*1) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*2), kbc(2*1), kad(1*2), kac(1*1), pbd(2*2), &
      pbc(2*1), pad(1*2), pac(1*1), prim(1*2*1*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*2) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:1*2) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_1_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*3), kbc(2*1), kad(1*3), kac(1*1), pbd(2*3), &
      pbc(2*1), pad(1*3), pac(1*1), prim(1*2*1*3), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*3) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:1*3) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,3
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_1_3
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_1_4(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*4), kbc(2*1), kad(1*4), kac(1*1), pbd(2*4), &
      pbc(2*1), pad(1*4), pac(1*1), prim(1*2*1*4), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*4) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:1*4) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,4
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_1_4
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_1_5(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*5), kbc(2*1), kad(1*5), kac(1*1), pbd(2*5), &
      pbc(2*1), pad(1*5), pac(1*1), prim(1*2*1*5), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*5) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:1*5) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,5
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_1_5
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_1_6(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*6), kbc(2*1), kad(1*6), kac(1*1), pbd(2*6), &
      pbc(2*1), pad(1*6), pac(1*1), prim(1*2*1*6), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*6) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:1*6) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,6
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_1_6
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_1_7(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*7), kbc(2*1), kad(1*7), kac(1*1), pbd(2*7), &
      pbc(2*1), pad(1*7), pac(1*1), prim(1*2*1*7), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*7) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:1*7) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,7
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_1_7
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_1_9(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*9), kbc(2*1), kad(1*9), kac(1*1), pbd(2*9), &
      pbc(2*1), pad(1*9), pac(1*1), prim(1*2*1*9), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*9) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:1*9) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,9
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_1_9
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*1), kad(1*md_max), kac(1*1), &
      pbd(2*md_max), pbc(2*1), pad(1*md_max), pac(1*1), prim(1*2*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*1), kbc(2*2), kad(1*1), kac(1*2), pbd(2*1), &
      pbc(2*2), pad(1*1), pac(1*2), prim(1*2*2*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*1) = 0.0_dp
      kbc(1:2*2) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,2
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_2_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_2_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*2), kbc(2*2), kad(1*2), kac(1*2), pbd(2*2), &
      pbc(2*2), pad(1*2), pac(1*2), prim(1*2*2*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*2) = 0.0_dp
      kbc(1:2*2) = 0.0_dp
      kad(1:1*2) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,2
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_2_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_2_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*3), kbc(2*2), kad(1*3), kac(1*2), pbd(2*3), &
      pbc(2*2), pad(1*3), pac(1*2), prim(1*2*2*3), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*3) = 0.0_dp
      kbc(1:2*2) = 0.0_dp
      kad(1:1*3) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,3
        DO mc = 1,2
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_2_3
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_2_4(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*4), kbc(2*2), kad(1*4), kac(1*2), pbd(2*4), &
      pbc(2*2), pad(1*4), pac(1*2), prim(1*2*2*4), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*4) = 0.0_dp
      kbc(1:2*2) = 0.0_dp
      kad(1:1*4) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,4
        DO mc = 1,2
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_2_4
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*2), kad(1*md_max), kac(1*2), &
      pbd(2*md_max), pbc(2*2), pad(1*md_max), pac(1*2), prim(1*2*2*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*2) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,2
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*1), kbc(2*3), kad(1*1), kac(1*3), pbd(2*1), &
      pbc(2*3), pad(1*1), pac(1*3), prim(1*2*3*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*1) = 0.0_dp
      kbc(1:2*3) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*3) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,3
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_3_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_3_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*2), kbc(2*3), kad(1*2), kac(1*3), pbd(2*2), &
      pbc(2*3), pad(1*2), pac(1*3), prim(1*2*3*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*2) = 0.0_dp
      kbc(1:2*3) = 0.0_dp
      kad(1:1*2) = 0.0_dp
      kac(1:1*3) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,3
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_3_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_3_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*3), kbc(2*3), kad(1*3), kac(1*3), pbd(2*3), &
      pbc(2*3), pad(1*3), pac(1*3), prim(1*2*3*3), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*3) = 0.0_dp
      kbc(1:2*3) = 0.0_dp
      kad(1:1*3) = 0.0_dp
      kac(1:1*3) = 0.0_dp
      p_index = 0
      DO md = 1,3
        DO mc = 1,3
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_3_3
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*3), kad(1*md_max), kac(1*3), &
      pbd(2*md_max), pbc(2*3), pad(1*md_max), pac(1*3), prim(1*2*3*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*3) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*3) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,3
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_3
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_4_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*1), kbc(2*4), kad(1*1), kac(1*4), pbd(2*1), &
      pbc(2*4), pad(1*1), pac(1*4), prim(1*2*4*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*1) = 0.0_dp
      kbc(1:2*4) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*4) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,4
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_4_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_4_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*2), kbc(2*4), kad(1*2), kac(1*4), pbd(2*2), &
      pbc(2*4), pad(1*2), pac(1*4), prim(1*2*4*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*2) = 0.0_dp
      kbc(1:2*4) = 0.0_dp
      kad(1:1*2) = 0.0_dp
      kac(1:1*4) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,4
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_4_2
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*4), kad(1*md_max), kac(1*4), &
      pbd(2*md_max), pbc(2*4), pad(1*md_max), pac(1*4), prim(1*2*4*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*4) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*4) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,4
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_4
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_5_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*1), kbc(2*5), kad(1*1), kac(1*5), pbd(2*1), &
      pbc(2*5), pad(1*1), pac(1*5), prim(1*2*5*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*1) = 0.0_dp
      kbc(1:2*5) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*5) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,5
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_5_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*5), kad(1*md_max), kac(1*5), &
      pbd(2*md_max), pbc(2*5), pad(1*md_max), pac(1*5), prim(1*2*5*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*5) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*5) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,5
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_5
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_6_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*1), kbc(2*6), kad(1*1), kac(1*6), pbd(2*1), &
      pbc(2*6), pad(1*1), pac(1*6), prim(1*2*6*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*1) = 0.0_dp
      kbc(1:2*6) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*6) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,6
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_6_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*6), kad(1*md_max), kac(1*6), &
      pbd(2*md_max), pbc(2*6), pad(1*md_max), pac(1*6), prim(1*2*6*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*6) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*6) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,6
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_6
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_7_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*1), kbc(2*7), kad(1*1), kac(1*7), pbd(2*1), &
      pbc(2*7), pad(1*1), pac(1*7), prim(1*2*7*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*1) = 0.0_dp
      kbc(1:2*7) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*7) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,7
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_7_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*7), kad(1*md_max), kac(1*7), &
      pbd(2*md_max), pbc(2*7), pad(1*md_max), pac(1*7), prim(1*2*7*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*7) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*7) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,7
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_7
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_9_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*1), kbc(2*9), kad(1*1), kac(1*9), pbd(2*1), &
      pbc(2*9), pad(1*1), pac(1*9), prim(1*2*9*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*1) = 0.0_dp
      kbc(1:2*9) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*9) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,9
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_9_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*9), kad(1*md_max), kac(1*9), &
      pbd(2*md_max), pbc(2*9), pad(1*md_max), pac(1*9), prim(1*2*9*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*9) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*9) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,9
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2_9
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*mc_max), kad(1*md_max), &
      kac(1*mc_max), pbd(2*md_max), pbc(2*mc_max), pad(1*md_max), &
      pac(1*mc_max), prim(1*2*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*mc_max) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_3_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(3*1), kbc(3*1), kad(1*1), kac(1*1), pbd(3*1), &
      pbc(3*1), pad(1*1), pac(1*1), prim(1*3*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*1) = 0.0_dp
      kbc(1:3*1) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_3_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_3_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(3*2), kbc(3*1), kad(1*2), kac(1*1), pbd(3*2), &
      pbc(3*1), pad(1*2), pac(1*1), prim(1*3*1*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*2) = 0.0_dp
      kbc(1:3*1) = 0.0_dp
      kad(1:1*2) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,1
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_3_1_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_3_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(3*3), kbc(3*1), kad(1*3), kac(1*1), pbd(3*3), &
      pbc(3*1), pad(1*3), pac(1*1), prim(1*3*1*3), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*3) = 0.0_dp
      kbc(1:3*1) = 0.0_dp
      kad(1:1*3) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,3
        DO mc = 1,1
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_3_1_3
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_3_1_4(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(3*4), kbc(3*1), kad(1*4), kac(1*1), pbd(3*4), &
      pbc(3*1), pad(1*4), pac(1*1), prim(1*3*1*4), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*4) = 0.0_dp
      kbc(1:3*1) = 0.0_dp
      kad(1:1*4) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,4
        DO mc = 1,1
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_3_1_4
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_3_1_5(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(3*5), kbc(3*1), kad(1*5), kac(1*1), pbd(3*5), &
      pbc(3*1), pad(1*5), pac(1*1), prim(1*3*1*5), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*5) = 0.0_dp
      kbc(1:3*1) = 0.0_dp
      kad(1:1*5) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,5
        DO mc = 1,1
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_3_1_5
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_3_1_6(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(3*6), kbc(3*1), kad(1*6), kac(1*1), pbd(3*6), &
      pbc(3*1), pad(1*6), pac(1*1), prim(1*3*1*6), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*6) = 0.0_dp
      kbc(1:3*1) = 0.0_dp
      kad(1:1*6) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,6
        DO mc = 1,1
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_3_1_6
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(3*md_max), kbc(3*1), kad(1*md_max), kac(1*1), &
      pbd(3*md_max), pbc(3*1), pad(1*md_max), pac(1*1), prim(1*3*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*md_max) = 0.0_dp
      kbc(1:3*1) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_3_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_3_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(3*1), kbc(3*2), kad(1*1), kac(1*2), pbd(3*1), &
      pbc(3*2), pad(1*1), pac(1*2), prim(1*3*2*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*1) = 0.0_dp
      kbc(1:3*2) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,2
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_3_2_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_3_2_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(3*2), kbc(3*2), kad(1*2), kac(1*2), pbd(3*2), &
      pbc(3*2), pad(1*2), pac(1*2), prim(1*3*2*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*2) = 0.0_dp
      kbc(1:3*2) = 0.0_dp
      kad(1:1*2) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,2
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_3_2_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_3_2_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(3*3), kbc(3*2), kad(1*3), kac(1*2), pbd(3*3), &
      pbc(3*2), pad(1*3), pac(1*2), prim(1*3*2*3), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*3) = 0.0_dp
      kbc(1:3*2) = 0.0_dp
      kad(1:1*3) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,3
        DO mc = 1,2
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_3_2_3
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(3*md_max), kbc(3*2), kad(1*md_max), kac(1*2), &
      pbd(3*md_max), pbc(3*2), pad(1*md_max), pac(1*2), prim(1*3*2*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*md_max) = 0.0_dp
      kbc(1:3*2) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,2
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_3_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_3_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(3*1), kbc(3*3), kad(1*1), kac(1*3), pbd(3*1), &
      pbc(3*3), pad(1*1), pac(1*3), prim(1*3*3*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*1) = 0.0_dp
      kbc(1:3*3) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*3) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,3
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_3_3_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_3_3_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(3*2), kbc(3*3), kad(1*2), kac(1*3), pbd(3*2), &
      pbc(3*3), pad(1*2), pac(1*3), prim(1*3*3*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*2) = 0.0_dp
      kbc(1:3*3) = 0.0_dp
      kad(1:1*2) = 0.0_dp
      kac(1:1*3) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,3
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_3_3_2
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_3_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(3*md_max), kbc(3*3), kad(1*md_max), kac(1*3), &
      pbd(3*md_max), pbc(3*3), pad(1*md_max), pac(1*3), prim(1*3*3*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*md_max) = 0.0_dp
      kbc(1:3*3) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*3) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,3
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_3_3
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_3_4_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(3*1), kbc(3*4), kad(1*1), kac(1*4), pbd(3*1), &
      pbc(3*4), pad(1*1), pac(1*4), prim(1*3*4*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*1) = 0.0_dp
      kbc(1:3*4) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*4) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,4
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_3_4_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_3_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(3*md_max), kbc(3*4), kad(1*md_max), kac(1*4), &
      pbd(3*md_max), pbc(3*4), pad(1*md_max), pac(1*4), prim(1*3*4*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*md_max) = 0.0_dp
      kbc(1:3*4) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*4) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,4
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_3_4
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_3_5_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(3*1), kbc(3*5), kad(1*1), kac(1*5), pbd(3*1), &
      pbc(3*5), pad(1*1), pac(1*5), prim(1*3*5*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*1) = 0.0_dp
      kbc(1:3*5) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*5) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,5
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_3_5_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_3_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(3*md_max), kbc(3*5), kad(1*md_max), kac(1*5), &
      pbd(3*md_max), pbc(3*5), pad(1*md_max), pac(1*5), prim(1*3*5*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*md_max) = 0.0_dp
      kbc(1:3*5) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*5) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,5
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_3_5
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_3_6_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(3*1), kbc(3*6), kad(1*1), kac(1*6), pbd(3*1), &
      pbc(3*6), pad(1*1), pac(1*6), prim(1*3*6*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*1) = 0.0_dp
      kbc(1:3*6) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*6) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,6
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_3_6_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_3_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(3*md_max), kbc(3*6), kad(1*md_max), kac(1*6), &
      pbd(3*md_max), pbc(3*6), pad(1*md_max), pac(1*6), prim(1*3*6*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*md_max) = 0.0_dp
      kbc(1:3*6) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*6) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,6
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_3_6
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(3*md_max), kbc(3*mc_max), kad(1*md_max), &
      kac(1*mc_max), pbd(3*md_max), pbc(3*mc_max), pad(1*md_max), &
      pac(1*mc_max), prim(1*3*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*md_max) = 0.0_dp
      kbc(1:3*mc_max) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_3
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_4_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(4*1), kbc(4*1), kad(1*1), kac(1*1), pbd(4*1), &
      pbc(4*1), pad(1*1), pac(1*1), prim(1*4*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*1) = 0.0_dp
      kbc(1:4*1) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_4_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_4_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(4*2), kbc(4*1), kad(1*2), kac(1*1), pbd(4*2), &
      pbc(4*1), pad(1*2), pac(1*1), prim(1*4*1*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*2) = 0.0_dp
      kbc(1:4*1) = 0.0_dp
      kad(1:1*2) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,1
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_4_1_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_4_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(4*3), kbc(4*1), kad(1*3), kac(1*1), pbd(4*3), &
      pbc(4*1), pad(1*3), pac(1*1), prim(1*4*1*3), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*3) = 0.0_dp
      kbc(1:4*1) = 0.0_dp
      kad(1:1*3) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,3
        DO mc = 1,1
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_4_1_3
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_4_1_4(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(4*4), kbc(4*1), kad(1*4), kac(1*1), pbd(4*4), &
      pbc(4*1), pad(1*4), pac(1*1), prim(1*4*1*4), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*4) = 0.0_dp
      kbc(1:4*1) = 0.0_dp
      kad(1:1*4) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,4
        DO mc = 1,1
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_4_1_4
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(4*md_max), kbc(4*1), kad(1*md_max), kac(1*1), &
      pbd(4*md_max), pbc(4*1), pad(1*md_max), pac(1*1), prim(1*4*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*md_max) = 0.0_dp
      kbc(1:4*1) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_4_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_4_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(4*1), kbc(4*2), kad(1*1), kac(1*2), pbd(4*1), &
      pbc(4*2), pad(1*1), pac(1*2), prim(1*4*2*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*1) = 0.0_dp
      kbc(1:4*2) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,2
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_4_2_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_4_2_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(4*2), kbc(4*2), kad(1*2), kac(1*2), pbd(4*2), &
      pbc(4*2), pad(1*2), pac(1*2), prim(1*4*2*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*2) = 0.0_dp
      kbc(1:4*2) = 0.0_dp
      kad(1:1*2) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,2
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_4_2_2
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_4_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(4*md_max), kbc(4*2), kad(1*md_max), kac(1*2), &
      pbd(4*md_max), pbc(4*2), pad(1*md_max), pac(1*2), prim(1*4*2*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*md_max) = 0.0_dp
      kbc(1:4*2) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,2
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_4_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_4_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(4*1), kbc(4*3), kad(1*1), kac(1*3), pbd(4*1), &
      pbc(4*3), pad(1*1), pac(1*3), prim(1*4*3*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*1) = 0.0_dp
      kbc(1:4*3) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*3) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,3
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_4_3_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_4_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(4*md_max), kbc(4*3), kad(1*md_max), kac(1*3), &
      pbd(4*md_max), pbc(4*3), pad(1*md_max), pac(1*3), prim(1*4*3*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*md_max) = 0.0_dp
      kbc(1:4*3) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*3) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,3
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_4_3
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_4_4_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(4*1), kbc(4*4), kad(1*1), kac(1*4), pbd(4*1), &
      pbc(4*4), pad(1*1), pac(1*4), prim(1*4*4*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*1) = 0.0_dp
      kbc(1:4*4) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*4) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,4
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_4_4_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_4_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(4*md_max), kbc(4*4), kad(1*md_max), kac(1*4), &
      pbd(4*md_max), pbc(4*4), pad(1*md_max), pac(1*4), prim(1*4*4*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*md_max) = 0.0_dp
      kbc(1:4*4) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*4) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,4
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_4_4
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(4*md_max), kbc(4*mc_max), kad(1*md_max), &
      kac(1*mc_max), pbd(4*md_max), pbc(4*mc_max), pad(1*md_max), &
      pac(1*mc_max), prim(1*4*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*md_max) = 0.0_dp
      kbc(1:4*mc_max) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_4
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_5_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(5*1), kbc(5*1), kad(1*1), kac(1*1), pbd(5*1), &
      pbc(5*1), pad(1*1), pac(1*1), prim(1*5*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:5*1) = 0.0_dp
      kbc(1:5*1) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,5
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*5+mb)
            p_bc = pbc((mc-1)*5+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*5+mb) = kbd((md-1)*5+mb) - ks_bd
            kbc((mc-1)*5+mb) = kbc((mc-1)*5+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_5_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_5_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(5*2), kbc(5*1), kad(1*2), kac(1*1), pbd(5*2), &
      pbc(5*1), pad(1*2), pac(1*1), prim(1*5*1*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:5*2) = 0.0_dp
      kbc(1:5*1) = 0.0_dp
      kad(1:1*2) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,1
          DO mb = 1,5
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*5+mb)
            p_bc = pbc((mc-1)*5+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*5+mb) = kbd((md-1)*5+mb) - ks_bd
            kbc((mc-1)*5+mb) = kbc((mc-1)*5+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_5_1_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_5_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(5*3), kbc(5*1), kad(1*3), kac(1*1), pbd(5*3), &
      pbc(5*1), pad(1*3), pac(1*1), prim(1*5*1*3), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:5*3) = 0.0_dp
      kbc(1:5*1) = 0.0_dp
      kad(1:1*3) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,3
        DO mc = 1,1
          DO mb = 1,5
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*5+mb)
            p_bc = pbc((mc-1)*5+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*5+mb) = kbd((md-1)*5+mb) - ks_bd
            kbc((mc-1)*5+mb) = kbc((mc-1)*5+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_5_1_3
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(5*md_max), kbc(5*1), kad(1*md_max), kac(1*1), &
      pbd(5*md_max), pbc(5*1), pad(1*md_max), pac(1*1), prim(1*5*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:5*md_max) = 0.0_dp
      kbc(1:5*1) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,5
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*5+mb)
            p_bc = pbc((mc-1)*5+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*5+mb) = kbd((md-1)*5+mb) - ks_bd
            kbc((mc-1)*5+mb) = kbc((mc-1)*5+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_5_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_5_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(5*1), kbc(5*2), kad(1*1), kac(1*2), pbd(5*1), &
      pbc(5*2), pad(1*1), pac(1*2), prim(1*5*2*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:5*1) = 0.0_dp
      kbc(1:5*2) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,2
          DO mb = 1,5
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*5+mb)
            p_bc = pbc((mc-1)*5+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*5+mb) = kbd((md-1)*5+mb) - ks_bd
            kbc((mc-1)*5+mb) = kbc((mc-1)*5+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_5_2_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_5_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(5*md_max), kbc(5*2), kad(1*md_max), kac(1*2), &
      pbd(5*md_max), pbc(5*2), pad(1*md_max), pac(1*2), prim(1*5*2*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:5*md_max) = 0.0_dp
      kbc(1:5*2) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,2
          DO mb = 1,5
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*5+mb)
            p_bc = pbc((mc-1)*5+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*5+mb) = kbd((md-1)*5+mb) - ks_bd
            kbc((mc-1)*5+mb) = kbc((mc-1)*5+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_5_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_5_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(5*1), kbc(5*3), kad(1*1), kac(1*3), pbd(5*1), &
      pbc(5*3), pad(1*1), pac(1*3), prim(1*5*3*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:5*1) = 0.0_dp
      kbc(1:5*3) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*3) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,3
          DO mb = 1,5
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*5+mb)
            p_bc = pbc((mc-1)*5+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*5+mb) = kbd((md-1)*5+mb) - ks_bd
            kbc((mc-1)*5+mb) = kbc((mc-1)*5+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_5_3_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_5_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(5*md_max), kbc(5*3), kad(1*md_max), kac(1*3), &
      pbd(5*md_max), pbc(5*3), pad(1*md_max), pac(1*3), prim(1*5*3*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:5*md_max) = 0.0_dp
      kbc(1:5*3) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*3) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,3
          DO mb = 1,5
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*5+mb)
            p_bc = pbc((mc-1)*5+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*5+mb) = kbd((md-1)*5+mb) - ks_bd
            kbc((mc-1)*5+mb) = kbc((mc-1)*5+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_5_3
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(5*md_max), kbc(5*mc_max), kad(1*md_max), &
      kac(1*mc_max), pbd(5*md_max), pbc(5*mc_max), pad(1*md_max), &
      pac(1*mc_max), prim(1*5*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:5*md_max) = 0.0_dp
      kbc(1:5*mc_max) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,5
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*5+mb)
            p_bc = pbc((mc-1)*5+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*5+mb) = kbd((md-1)*5+mb) - ks_bd
            kbc((mc-1)*5+mb) = kbc((mc-1)*5+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_5
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_6_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(6*1), kbc(6*1), kad(1*1), kac(1*1), pbd(6*1), &
      pbc(6*1), pad(1*1), pac(1*1), prim(1*6*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:6*1) = 0.0_dp
      kbc(1:6*1) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,6
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*6+mb)
            p_bc = pbc((mc-1)*6+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*6+mb) = kbd((md-1)*6+mb) - ks_bd
            kbc((mc-1)*6+mb) = kbc((mc-1)*6+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_6_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_6_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(6*2), kbc(6*1), kad(1*2), kac(1*1), pbd(6*2), &
      pbc(6*1), pad(1*2), pac(1*1), prim(1*6*1*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:6*2) = 0.0_dp
      kbc(1:6*1) = 0.0_dp
      kad(1:1*2) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,1
          DO mb = 1,6
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*6+mb)
            p_bc = pbc((mc-1)*6+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*6+mb) = kbd((md-1)*6+mb) - ks_bd
            kbc((mc-1)*6+mb) = kbc((mc-1)*6+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_6_1_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_6_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(6*3), kbc(6*1), kad(1*3), kac(1*1), pbd(6*3), &
      pbc(6*1), pad(1*3), pac(1*1), prim(1*6*1*3), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:6*3) = 0.0_dp
      kbc(1:6*1) = 0.0_dp
      kad(1:1*3) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,3
        DO mc = 1,1
          DO mb = 1,6
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*6+mb)
            p_bc = pbc((mc-1)*6+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*6+mb) = kbd((md-1)*6+mb) - ks_bd
            kbc((mc-1)*6+mb) = kbc((mc-1)*6+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_6_1_3
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(6*md_max), kbc(6*1), kad(1*md_max), kac(1*1), &
      pbd(6*md_max), pbc(6*1), pad(1*md_max), pac(1*1), prim(1*6*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:6*md_max) = 0.0_dp
      kbc(1:6*1) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,6
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*6+mb)
            p_bc = pbc((mc-1)*6+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*6+mb) = kbd((md-1)*6+mb) - ks_bd
            kbc((mc-1)*6+mb) = kbc((mc-1)*6+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_6_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_6_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(6*1), kbc(6*2), kad(1*1), kac(1*2), pbd(6*1), &
      pbc(6*2), pad(1*1), pac(1*2), prim(1*6*2*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:6*1) = 0.0_dp
      kbc(1:6*2) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,2
          DO mb = 1,6
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*6+mb)
            p_bc = pbc((mc-1)*6+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*6+mb) = kbd((md-1)*6+mb) - ks_bd
            kbc((mc-1)*6+mb) = kbc((mc-1)*6+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_6_2_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_6_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(6*md_max), kbc(6*2), kad(1*md_max), kac(1*2), &
      pbd(6*md_max), pbc(6*2), pad(1*md_max), pac(1*2), prim(1*6*2*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:6*md_max) = 0.0_dp
      kbc(1:6*2) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,2
          DO mb = 1,6
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*6+mb)
            p_bc = pbc((mc-1)*6+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*6+mb) = kbd((md-1)*6+mb) - ks_bd
            kbc((mc-1)*6+mb) = kbc((mc-1)*6+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_6_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_6_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(6*1), kbc(6*3), kad(1*1), kac(1*3), pbd(6*1), &
      pbc(6*3), pad(1*1), pac(1*3), prim(1*6*3*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:6*1) = 0.0_dp
      kbc(1:6*3) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*3) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,3
          DO mb = 1,6
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*6+mb)
            p_bc = pbc((mc-1)*6+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*6+mb) = kbd((md-1)*6+mb) - ks_bd
            kbc((mc-1)*6+mb) = kbc((mc-1)*6+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_6_3_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_6_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(6*md_max), kbc(6*3), kad(1*md_max), kac(1*3), &
      pbd(6*md_max), pbc(6*3), pad(1*md_max), pac(1*3), prim(1*6*3*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:6*md_max) = 0.0_dp
      kbc(1:6*3) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*3) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,3
          DO mb = 1,6
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*6+mb)
            p_bc = pbc((mc-1)*6+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*6+mb) = kbd((md-1)*6+mb) - ks_bd
            kbc((mc-1)*6+mb) = kbc((mc-1)*6+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_6_3
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(6*md_max), kbc(6*mc_max), kad(1*md_max), &
      kac(1*mc_max), pbd(6*md_max), pbc(6*mc_max), pad(1*md_max), &
      pac(1*mc_max), prim(1*6*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:6*md_max) = 0.0_dp
      kbc(1:6*mc_max) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,6
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*6+mb)
            p_bc = pbc((mc-1)*6+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*6+mb) = kbd((md-1)*6+mb) - ks_bd
            kbc((mc-1)*6+mb) = kbc((mc-1)*6+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_6
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_7_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(7*1), kbc(7*1), kad(1*1), kac(1*1), pbd(7*1), &
      pbc(7*1), pad(1*1), pac(1*1), prim(1*7*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:7*1) = 0.0_dp
      kbc(1:7*1) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,7
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*7+mb)
            p_bc = pbc((mc-1)*7+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*7+mb) = kbd((md-1)*7+mb) - ks_bd
            kbc((mc-1)*7+mb) = kbc((mc-1)*7+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_7_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_7_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(7*2), kbc(7*1), kad(1*2), kac(1*1), pbd(7*2), &
      pbc(7*1), pad(1*2), pac(1*1), prim(1*7*1*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:7*2) = 0.0_dp
      kbc(1:7*1) = 0.0_dp
      kad(1:1*2) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,1
          DO mb = 1,7
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*7+mb)
            p_bc = pbc((mc-1)*7+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*7+mb) = kbd((md-1)*7+mb) - ks_bd
            kbc((mc-1)*7+mb) = kbc((mc-1)*7+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_7_1_2
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_7_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(7*md_max), kbc(7*1), kad(1*md_max), kac(1*1), &
      pbd(7*md_max), pbc(7*1), pad(1*md_max), pac(1*1), prim(1*7*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:7*md_max) = 0.0_dp
      kbc(1:7*1) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,7
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*7+mb)
            p_bc = pbc((mc-1)*7+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*7+mb) = kbd((md-1)*7+mb) - ks_bd
            kbc((mc-1)*7+mb) = kbc((mc-1)*7+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_7_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_7_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(7*1), kbc(7*2), kad(1*1), kac(1*2), pbd(7*1), &
      pbc(7*2), pad(1*1), pac(1*2), prim(1*7*2*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:7*1) = 0.0_dp
      kbc(1:7*2) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,2
          DO mb = 1,7
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*7+mb)
            p_bc = pbc((mc-1)*7+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*7+mb) = kbd((md-1)*7+mb) - ks_bd
            kbc((mc-1)*7+mb) = kbc((mc-1)*7+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_7_2_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_7_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(7*md_max), kbc(7*2), kad(1*md_max), kac(1*2), &
      pbd(7*md_max), pbc(7*2), pad(1*md_max), pac(1*2), prim(1*7*2*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:7*md_max) = 0.0_dp
      kbc(1:7*2) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,2
          DO mb = 1,7
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*7+mb)
            p_bc = pbc((mc-1)*7+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*7+mb) = kbd((md-1)*7+mb) - ks_bd
            kbc((mc-1)*7+mb) = kbc((mc-1)*7+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_7_2
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(7*md_max), kbc(7*mc_max), kad(1*md_max), &
      kac(1*mc_max), pbd(7*md_max), pbc(7*mc_max), pad(1*md_max), &
      pac(1*mc_max), prim(1*7*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:7*md_max) = 0.0_dp
      kbc(1:7*mc_max) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,7
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*7+mb)
            p_bc = pbc((mc-1)*7+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*7+mb) = kbd((md-1)*7+mb) - ks_bd
            kbc((mc-1)*7+mb) = kbc((mc-1)*7+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_7
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_9_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(9*1), kbc(9*1), kad(1*1), kac(1*1), pbd(9*1), &
      pbc(9*1), pad(1*1), pac(1*1), prim(1*9*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:9*1) = 0.0_dp
      kbc(1:9*1) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,9
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*9+mb)
            p_bc = pbc((mc-1)*9+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*9+mb) = kbd((md-1)*9+mb) - ks_bd
            kbc((mc-1)*9+mb) = kbc((mc-1)*9+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_9_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_9_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(9*2), kbc(9*1), kad(1*2), kac(1*1), pbd(9*2), &
      pbc(9*1), pad(1*2), pac(1*1), prim(1*9*1*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:9*2) = 0.0_dp
      kbc(1:9*1) = 0.0_dp
      kad(1:1*2) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,1
          DO mb = 1,9
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*9+mb)
            p_bc = pbc((mc-1)*9+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*9+mb) = kbd((md-1)*9+mb) - ks_bd
            kbc((mc-1)*9+mb) = kbc((mc-1)*9+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_9_1_2
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_9_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(9*md_max), kbc(9*1), kad(1*md_max), kac(1*1), &
      pbd(9*md_max), pbc(9*1), pad(1*md_max), pac(1*1), prim(1*9*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:9*md_max) = 0.0_dp
      kbc(1:9*1) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,9
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*9+mb)
            p_bc = pbc((mc-1)*9+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*9+mb) = kbd((md-1)*9+mb) - ks_bd
            kbc((mc-1)*9+mb) = kbc((mc-1)*9+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_9_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_9_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(9*1), kbc(9*2), kad(1*1), kac(1*2), pbd(9*1), &
      pbc(9*2), pad(1*1), pac(1*2), prim(1*9*2*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:9*1) = 0.0_dp
      kbc(1:9*2) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,2
          DO mb = 1,9
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*9+mb)
            p_bc = pbc((mc-1)*9+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*9+mb) = kbd((md-1)*9+mb) - ks_bd
            kbc((mc-1)*9+mb) = kbc((mc-1)*9+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_9_2_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_9_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(9*md_max), kbc(9*2), kad(1*md_max), kac(1*2), &
      pbd(9*md_max), pbc(9*2), pad(1*md_max), pac(1*2), prim(1*9*2*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:9*md_max) = 0.0_dp
      kbc(1:9*2) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*2) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,2
          DO mb = 1,9
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*9+mb)
            p_bc = pbc((mc-1)*9+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*9+mb) = kbd((md-1)*9+mb) - ks_bd
            kbc((mc-1)*9+mb) = kbc((mc-1)*9+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_9_2
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(9*md_max), kbc(9*mc_max), kad(1*md_max), &
      kac(1*mc_max), pbd(9*md_max), pbc(9*mc_max), pad(1*md_max), &
      pac(1*mc_max), prim(1*9*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:9*md_max) = 0.0_dp
      kbc(1:9*mc_max) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,9
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*9+mb)
            p_bc = pbc((mc-1)*9+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*9+mb) = kbd((md-1)*9+mb) - ks_bd
            kbc((mc-1)*9+mb) = kbc((mc-1)*9+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_9
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_10_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(10*1), kbc(10*1), kad(1*1), kac(1*1), pbd(10*1), &
      pbc(10*1), pad(1*1), pac(1*1), prim(1*10*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:10*1) = 0.0_dp
      kbc(1:10*1) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,10
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*10+mb)
            p_bc = pbc((mc-1)*10+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*10+mb) = kbd((md-1)*10+mb) - ks_bd
            kbc((mc-1)*10+mb) = kbc((mc-1)*10+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_10_1_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_10_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(10*md_max), kbc(10*1), kad(1*md_max), kac(1*1), &
      pbd(10*md_max), pbc(10*1), pad(1*md_max), pac(1*1), &
      prim(1*10*1*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:10*md_max) = 0.0_dp
      kbc(1:10*1) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,10
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*10+mb)
            p_bc = pbc((mc-1)*10+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*10+mb) = kbd((md-1)*10+mb) - ks_bd
            kbc((mc-1)*10+mb) = kbc((mc-1)*10+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_10_1
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(10*md_max), kbc(10*mc_max), kad(1*md_max), &
      kac(1*mc_max), pbd(10*md_max), pbc(10*mc_max), pad(1*md_max), &
      pac(1*mc_max), prim(1*10*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:10*md_max) = 0.0_dp
      kbc(1:10*mc_max) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,10
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*10+mb)
            p_bc = pbc((mc-1)*10+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*10+mb) = kbd((md-1)*10+mb) - ks_bd
            kbc((mc-1)*10+mb) = kbc((mc-1)*10+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_10
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_11_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(11*1), kbc(11*1), kad(1*1), kac(1*1), pbd(11*1), &
      pbc(11*1), pad(1*1), pac(1*1), prim(1*11*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:11*1) = 0.0_dp
      kbc(1:11*1) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,11
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*11+mb)
            p_bc = pbc((mc-1)*11+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*11+mb) = kbd((md-1)*11+mb) - ks_bd
            kbc((mc-1)*11+mb) = kbc((mc-1)*11+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_11_1_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_11_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(11*md_max), kbc(11*1), kad(1*md_max), kac(1*1), &
      pbd(11*md_max), pbc(11*1), pad(1*md_max), pac(1*1), &
      prim(1*11*1*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:11*md_max) = 0.0_dp
      kbc(1:11*1) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,11
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*11+mb)
            p_bc = pbc((mc-1)*11+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*11+mb) = kbd((md-1)*11+mb) - ks_bd
            kbc((mc-1)*11+mb) = kbc((mc-1)*11+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_11_1
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(11*md_max), kbc(11*mc_max), kad(1*md_max), &
      kac(1*mc_max), pbd(11*md_max), pbc(11*mc_max), pad(1*md_max), &
      pac(1*mc_max), prim(1*11*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:11*md_max) = 0.0_dp
      kbc(1:11*mc_max) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,11
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*11+mb)
            p_bc = pbc((mc-1)*11+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*11+mb) = kbd((md-1)*11+mb) - ks_bd
            kbc((mc-1)*11+mb) = kbc((mc-1)*11+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_11
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_15_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(15*1), kbc(15*1), kad(1*1), kac(1*1), pbd(15*1), &
      pbc(15*1), pad(1*1), pac(1*1), prim(1*15*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:15*1) = 0.0_dp
      kbc(1:15*1) = 0.0_dp
      kad(1:1*1) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,15
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*15+mb)
            p_bc = pbc((mc-1)*15+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*15+mb) = kbd((md-1)*15+mb) - ks_bd
            kbc((mc-1)*15+mb) = kbc((mc-1)*15+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_15_1_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_15_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(15*md_max), kbc(15*1), kad(1*md_max), kac(1*1), &
      pbd(15*md_max), pbc(15*1), pad(1*md_max), pac(1*1), &
      prim(1*15*1*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:15*md_max) = 0.0_dp
      kbc(1:15*1) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,15
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*15+mb)
            p_bc = pbc((mc-1)*15+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*15+mb) = kbd((md-1)*15+mb) - ks_bd
            kbc((mc-1)*15+mb) = kbc((mc-1)*15+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_15_1
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_1_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(15*md_max), kbc(15*mc_max), kad(1*md_max), &
      kac(1*mc_max), pbd(15*md_max), pbc(15*mc_max), pad(1*md_max), &
      pac(1*mc_max), prim(1*15*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:15*md_max) = 0.0_dp
      kbc(1:15*mc_max) = 0.0_dp
      kad(1:1*md_max) = 0.0_dp
      kac(1:1*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,15
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*15+mb)
            p_bc = pbc((mc-1)*15+mb)
            DO ma = 1,1
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*1+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*1+ma)
              kad((md-1)*1+ma) =  kad((md-1)*1+ma)-tmp*p_bc
              kac((mc-1)*1+ma) = kac((mc-1)*1+ma)-tmp*p_bd
            END DO
            kbd((md-1)*15+mb) = kbd((md-1)*15+mb) - ks_bd
            kbc((mc-1)*15+mb) = kbc((mc-1)*15+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_1_15
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*1), kad(2*1), kac(2*1), pbd(1*1), &
      pbc(1*1), pad(2*1), pac(2*1), prim(2*1*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:2*1) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*2), kbc(1*1), kad(2*2), kac(2*1), pbd(1*2), &
      pbc(1*1), pad(2*2), pac(2*1), prim(2*1*1*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*2) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:2*2) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_1_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*3), kbc(1*1), kad(2*3), kac(2*1), pbd(1*3), &
      pbc(1*1), pad(2*3), pac(2*1), prim(2*1*1*3), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*3) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:2*3) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,3
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_1_3
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_1_4(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*4), kbc(1*1), kad(2*4), kac(2*1), pbd(1*4), &
      pbc(1*1), pad(2*4), pac(2*1), prim(2*1*1*4), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*4) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:2*4) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,4
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_1_4
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_1_5(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*5), kbc(1*1), kad(2*5), kac(2*1), pbd(1*5), &
      pbc(1*1), pad(2*5), pac(2*1), prim(2*1*1*5), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*5) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:2*5) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,5
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_1_5
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_1_6(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*6), kbc(1*1), kad(2*6), kac(2*1), pbd(1*6), &
      pbc(1*1), pad(2*6), pac(2*1), prim(2*1*1*6), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*6) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:2*6) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,6
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_1_6
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_1_7(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*7), kbc(1*1), kad(2*7), kac(2*1), pbd(1*7), &
      pbc(1*1), pad(2*7), pac(2*1), prim(2*1*1*7), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*7) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:2*7) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,7
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_1_7
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_1_9(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*9), kbc(1*1), kad(2*9), kac(2*1), pbd(1*9), &
      pbc(1*1), pad(2*9), pac(2*1), prim(2*1*1*9), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*9) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:2*9) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,9
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_1_9
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*1), kad(2*md_max), kac(2*1), &
      pbd(1*md_max), pbc(1*1), pad(2*md_max), pac(2*1), prim(2*1*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*2), kad(2*1), kac(2*2), pbd(1*1), &
      pbc(1*2), pad(2*1), pac(2*2), prim(2*1*2*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:2*1) = 0.0_dp
      kac(1:2*2) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_2_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_2_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*2), kbc(1*2), kad(2*2), kac(2*2), pbd(1*2), &
      pbc(1*2), pad(2*2), pac(2*2), prim(2*1*2*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*2) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:2*2) = 0.0_dp
      kac(1:2*2) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_2_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_2_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*3), kbc(1*2), kad(2*3), kac(2*2), pbd(1*3), &
      pbc(1*2), pad(2*3), pac(2*2), prim(2*1*2*3), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*3) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:2*3) = 0.0_dp
      kac(1:2*2) = 0.0_dp
      p_index = 0
      DO md = 1,3
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_2_3
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_2_4(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*4), kbc(1*2), kad(2*4), kac(2*2), pbd(1*4), &
      pbc(1*2), pad(2*4), pac(2*2), prim(2*1*2*4), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*4) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:2*4) = 0.0_dp
      kac(1:2*2) = 0.0_dp
      p_index = 0
      DO md = 1,4
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_2_4
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*2), kad(2*md_max), kac(2*2), &
      pbd(1*md_max), pbc(1*2), pad(2*md_max), pac(2*2), prim(2*1*2*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*2) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*3), kad(2*1), kac(2*3), pbd(1*1), &
      pbc(1*3), pad(2*1), pac(2*3), prim(2*1*3*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*3) = 0.0_dp
      kad(1:2*1) = 0.0_dp
      kac(1:2*3) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,3
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_3_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_3_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*2), kbc(1*3), kad(2*2), kac(2*3), pbd(1*2), &
      pbc(1*3), pad(2*2), pac(2*3), prim(2*1*3*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*2) = 0.0_dp
      kbc(1:1*3) = 0.0_dp
      kad(1:2*2) = 0.0_dp
      kac(1:2*3) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,3
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_3_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_3_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*3), kbc(1*3), kad(2*3), kac(2*3), pbd(1*3), &
      pbc(1*3), pad(2*3), pac(2*3), prim(2*1*3*3), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*3) = 0.0_dp
      kbc(1:1*3) = 0.0_dp
      kad(1:2*3) = 0.0_dp
      kac(1:2*3) = 0.0_dp
      p_index = 0
      DO md = 1,3
        DO mc = 1,3
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_3_3
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*3), kad(2*md_max), kac(2*3), &
      pbd(1*md_max), pbc(1*3), pad(2*md_max), pac(2*3), prim(2*1*3*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*3) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*3) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,3
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_3
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_4_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*4), kad(2*1), kac(2*4), pbd(1*1), &
      pbc(1*4), pad(2*1), pac(2*4), prim(2*1*4*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*4) = 0.0_dp
      kad(1:2*1) = 0.0_dp
      kac(1:2*4) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,4
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_4_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_4_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*2), kbc(1*4), kad(2*2), kac(2*4), pbd(1*2), &
      pbc(1*4), pad(2*2), pac(2*4), prim(2*1*4*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*2) = 0.0_dp
      kbc(1:1*4) = 0.0_dp
      kad(1:2*2) = 0.0_dp
      kac(1:2*4) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,4
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_4_2
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*4), kad(2*md_max), kac(2*4), &
      pbd(1*md_max), pbc(1*4), pad(2*md_max), pac(2*4), prim(2*1*4*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*4) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*4) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,4
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_4
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_5_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*5), kad(2*1), kac(2*5), pbd(1*1), &
      pbc(1*5), pad(2*1), pac(2*5), prim(2*1*5*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*5) = 0.0_dp
      kad(1:2*1) = 0.0_dp
      kac(1:2*5) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,5
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_5_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*5), kad(2*md_max), kac(2*5), &
      pbd(1*md_max), pbc(1*5), pad(2*md_max), pac(2*5), prim(2*1*5*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*5) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*5) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,5
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_5
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_6_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*6), kad(2*1), kac(2*6), pbd(1*1), &
      pbc(1*6), pad(2*1), pac(2*6), prim(2*1*6*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*6) = 0.0_dp
      kad(1:2*1) = 0.0_dp
      kac(1:2*6) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,6
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_6_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*6), kad(2*md_max), kac(2*6), &
      pbd(1*md_max), pbc(1*6), pad(2*md_max), pac(2*6), prim(2*1*6*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*6) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*6) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,6
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_6
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_7_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*7), kad(2*1), kac(2*7), pbd(1*1), &
      pbc(1*7), pad(2*1), pac(2*7), prim(2*1*7*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*7) = 0.0_dp
      kad(1:2*1) = 0.0_dp
      kac(1:2*7) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,7
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_7_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_7(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*7), kad(2*md_max), kac(2*7), &
      pbd(1*md_max), pbc(1*7), pad(2*md_max), pac(2*7), prim(2*1*7*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*7) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*7) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,7
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_7
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_9_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*9), kad(2*1), kac(2*9), pbd(1*1), &
      pbc(1*9), pad(2*1), pac(2*9), prim(2*1*9*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*9) = 0.0_dp
      kad(1:2*1) = 0.0_dp
      kac(1:2*9) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,9
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_9_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1_9(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*9), kad(2*md_max), kac(2*9), &
      pbd(1*md_max), pbc(1*9), pad(2*md_max), pac(2*9), prim(2*1*9*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*9) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*9) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,9
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1_9
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*mc_max), kad(2*md_max), &
      kac(2*mc_max), pbd(1*md_max), pbc(1*mc_max), pad(2*md_max), &
      pac(2*mc_max), prim(2*1*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*mc_max) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_2_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*1), kbc(2*1), kad(2*1), kac(2*1), pbd(2*1), &
      pbc(2*1), pad(2*1), pac(2*1), prim(2*2*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*1) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:2*1) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_2_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_2_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*2), kbc(2*1), kad(2*2), kac(2*1), pbd(2*2), &
      pbc(2*1), pad(2*2), pac(2*1), prim(2*2*1*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*2) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:2*2) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_2_1_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_2_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*3), kbc(2*1), kad(2*3), kac(2*1), pbd(2*3), &
      pbc(2*1), pad(2*3), pac(2*1), prim(2*2*1*3), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*3) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:2*3) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,3
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_2_1_3
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_2_1_4(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*4), kbc(2*1), kad(2*4), kac(2*1), pbd(2*4), &
      pbc(2*1), pad(2*4), pac(2*1), prim(2*2*1*4), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*4) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:2*4) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,4
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_2_1_4
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*1), kad(2*md_max), kac(2*1), &
      pbd(2*md_max), pbc(2*1), pad(2*md_max), pac(2*1), prim(2*2*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_2_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_2_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*1), kbc(2*2), kad(2*1), kac(2*2), pbd(2*1), &
      pbc(2*2), pad(2*1), pac(2*2), prim(2*2*2*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*1) = 0.0_dp
      kbc(1:2*2) = 0.0_dp
      kad(1:2*1) = 0.0_dp
      kac(1:2*2) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,2
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_2_2_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_2_2_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*2), kbc(2*2), kad(2*2), kac(2*2), pbd(2*2), &
      pbc(2*2), pad(2*2), pac(2*2), prim(2*2*2*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*2) = 0.0_dp
      kbc(1:2*2) = 0.0_dp
      kad(1:2*2) = 0.0_dp
      kac(1:2*2) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,2
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_2_2_2
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*2), kad(2*md_max), kac(2*2), &
      pbd(2*md_max), pbc(2*2), pad(2*md_max), pac(2*2), prim(2*2*2*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*2) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*2) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,2
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_2_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_2_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*1), kbc(2*3), kad(2*1), kac(2*3), pbd(2*1), &
      pbc(2*3), pad(2*1), pac(2*3), prim(2*2*3*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*1) = 0.0_dp
      kbc(1:2*3) = 0.0_dp
      kad(1:2*1) = 0.0_dp
      kac(1:2*3) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,3
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_2_3_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*3), kad(2*md_max), kac(2*3), &
      pbd(2*md_max), pbc(2*3), pad(2*md_max), pac(2*3), prim(2*2*3*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*3) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*3) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,3
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_2_3
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_2_4_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*1), kbc(2*4), kad(2*1), kac(2*4), pbd(2*1), &
      pbc(2*4), pad(2*1), pac(2*4), prim(2*2*4*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*1) = 0.0_dp
      kbc(1:2*4) = 0.0_dp
      kad(1:2*1) = 0.0_dp
      kac(1:2*4) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,4
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_2_4_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_2_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*4), kad(2*md_max), kac(2*4), &
      pbd(2*md_max), pbc(2*4), pad(2*md_max), pac(2*4), prim(2*2*4*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*4) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*4) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,4
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_2_4
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*mc_max), kad(2*md_max), &
      kac(2*mc_max), pbd(2*md_max), pbc(2*mc_max), pad(2*md_max), &
      pac(2*mc_max), prim(2*2*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*mc_max) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_3_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(3*1), kbc(3*1), kad(2*1), kac(2*1), pbd(3*1), &
      pbc(3*1), pad(2*1), pac(2*1), prim(2*3*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*1) = 0.0_dp
      kbc(1:3*1) = 0.0_dp
      kad(1:2*1) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_3_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_3_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(3*2), kbc(3*1), kad(2*2), kac(2*1), pbd(3*2), &
      pbc(3*1), pad(2*2), pac(2*1), prim(2*3*1*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*2) = 0.0_dp
      kbc(1:3*1) = 0.0_dp
      kad(1:2*2) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,1
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_3_1_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_3_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(3*3), kbc(3*1), kad(2*3), kac(2*1), pbd(3*3), &
      pbc(3*1), pad(2*3), pac(2*1), prim(2*3*1*3), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*3) = 0.0_dp
      kbc(1:3*1) = 0.0_dp
      kad(1:2*3) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,3
        DO mc = 1,1
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_3_1_3
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(3*md_max), kbc(3*1), kad(2*md_max), kac(2*1), &
      pbd(3*md_max), pbc(3*1), pad(2*md_max), pac(2*1), prim(2*3*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*md_max) = 0.0_dp
      kbc(1:3*1) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_3_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_3_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(3*1), kbc(3*2), kad(2*1), kac(2*2), pbd(3*1), &
      pbc(3*2), pad(2*1), pac(2*2), prim(2*3*2*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*1) = 0.0_dp
      kbc(1:3*2) = 0.0_dp
      kad(1:2*1) = 0.0_dp
      kac(1:2*2) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,2
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_3_2_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(3*md_max), kbc(3*2), kad(2*md_max), kac(2*2), &
      pbd(3*md_max), pbc(3*2), pad(2*md_max), pac(2*2), prim(2*3*2*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*md_max) = 0.0_dp
      kbc(1:3*2) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*2) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,2
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_3_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_3_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(3*1), kbc(3*3), kad(2*1), kac(2*3), pbd(3*1), &
      pbc(3*3), pad(2*1), pac(2*3), prim(2*3*3*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*1) = 0.0_dp
      kbc(1:3*3) = 0.0_dp
      kad(1:2*1) = 0.0_dp
      kac(1:2*3) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,3
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_3_3_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_3_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(3*md_max), kbc(3*3), kad(2*md_max), kac(2*3), &
      pbd(3*md_max), pbc(3*3), pad(2*md_max), pac(2*3), prim(2*3*3*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*md_max) = 0.0_dp
      kbc(1:3*3) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*3) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,3
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_3_3
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(3*md_max), kbc(3*mc_max), kad(2*md_max), &
      kac(2*mc_max), pbd(3*md_max), pbc(3*mc_max), pad(2*md_max), &
      pac(2*mc_max), prim(2*3*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*md_max) = 0.0_dp
      kbc(1:3*mc_max) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_3
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_4_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(4*1), kbc(4*1), kad(2*1), kac(2*1), pbd(4*1), &
      pbc(4*1), pad(2*1), pac(2*1), prim(2*4*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*1) = 0.0_dp
      kbc(1:4*1) = 0.0_dp
      kad(1:2*1) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_4_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_4_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(4*2), kbc(4*1), kad(2*2), kac(2*1), pbd(4*2), &
      pbc(4*1), pad(2*2), pac(2*1), prim(2*4*1*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*2) = 0.0_dp
      kbc(1:4*1) = 0.0_dp
      kad(1:2*2) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,1
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_4_1_2
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(4*md_max), kbc(4*1), kad(2*md_max), kac(2*1), &
      pbd(4*md_max), pbc(4*1), pad(2*md_max), pac(2*1), prim(2*4*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*md_max) = 0.0_dp
      kbc(1:4*1) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_4_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_4_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(4*1), kbc(4*2), kad(2*1), kac(2*2), pbd(4*1), &
      pbc(4*2), pad(2*1), pac(2*2), prim(2*4*2*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*1) = 0.0_dp
      kbc(1:4*2) = 0.0_dp
      kad(1:2*1) = 0.0_dp
      kac(1:2*2) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,2
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_4_2_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_4_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(4*md_max), kbc(4*2), kad(2*md_max), kac(2*2), &
      pbd(4*md_max), pbc(4*2), pad(2*md_max), pac(2*2), prim(2*4*2*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*md_max) = 0.0_dp
      kbc(1:4*2) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*2) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,2
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_4_2
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(4*md_max), kbc(4*mc_max), kad(2*md_max), &
      kac(2*mc_max), pbd(4*md_max), pbc(4*mc_max), pad(2*md_max), &
      pac(2*mc_max), prim(2*4*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*md_max) = 0.0_dp
      kbc(1:4*mc_max) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_4
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_5_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(5*1), kbc(5*1), kad(2*1), kac(2*1), pbd(5*1), &
      pbc(5*1), pad(2*1), pac(2*1), prim(2*5*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:5*1) = 0.0_dp
      kbc(1:5*1) = 0.0_dp
      kad(1:2*1) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,5
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*5+mb)
            p_bc = pbc((mc-1)*5+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*5+mb) = kbd((md-1)*5+mb) - ks_bd
            kbc((mc-1)*5+mb) = kbc((mc-1)*5+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_5_1_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(5*md_max), kbc(5*1), kad(2*md_max), kac(2*1), &
      pbd(5*md_max), pbc(5*1), pad(2*md_max), pac(2*1), prim(2*5*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:5*md_max) = 0.0_dp
      kbc(1:5*1) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,5
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*5+mb)
            p_bc = pbc((mc-1)*5+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*5+mb) = kbd((md-1)*5+mb) - ks_bd
            kbc((mc-1)*5+mb) = kbc((mc-1)*5+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_5_1
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(5*md_max), kbc(5*mc_max), kad(2*md_max), &
      kac(2*mc_max), pbd(5*md_max), pbc(5*mc_max), pad(2*md_max), &
      pac(2*mc_max), prim(2*5*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:5*md_max) = 0.0_dp
      kbc(1:5*mc_max) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,5
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*5+mb)
            p_bc = pbc((mc-1)*5+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*5+mb) = kbd((md-1)*5+mb) - ks_bd
            kbc((mc-1)*5+mb) = kbc((mc-1)*5+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_5
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_6_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(6*1), kbc(6*1), kad(2*1), kac(2*1), pbd(6*1), &
      pbc(6*1), pad(2*1), pac(2*1), prim(2*6*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:6*1) = 0.0_dp
      kbc(1:6*1) = 0.0_dp
      kad(1:2*1) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,6
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*6+mb)
            p_bc = pbc((mc-1)*6+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*6+mb) = kbd((md-1)*6+mb) - ks_bd
            kbc((mc-1)*6+mb) = kbc((mc-1)*6+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_6_1_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(6*md_max), kbc(6*1), kad(2*md_max), kac(2*1), &
      pbd(6*md_max), pbc(6*1), pad(2*md_max), pac(2*1), prim(2*6*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:6*md_max) = 0.0_dp
      kbc(1:6*1) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,6
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*6+mb)
            p_bc = pbc((mc-1)*6+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*6+mb) = kbd((md-1)*6+mb) - ks_bd
            kbc((mc-1)*6+mb) = kbc((mc-1)*6+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_6_1
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(6*md_max), kbc(6*mc_max), kad(2*md_max), &
      kac(2*mc_max), pbd(6*md_max), pbc(6*mc_max), pad(2*md_max), &
      pac(2*mc_max), prim(2*6*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:6*md_max) = 0.0_dp
      kbc(1:6*mc_max) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,6
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*6+mb)
            p_bc = pbc((mc-1)*6+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*6+mb) = kbd((md-1)*6+mb) - ks_bd
            kbc((mc-1)*6+mb) = kbc((mc-1)*6+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_6
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_7_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(7*1), kbc(7*1), kad(2*1), kac(2*1), pbd(7*1), &
      pbc(7*1), pad(2*1), pac(2*1), prim(2*7*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:7*1) = 0.0_dp
      kbc(1:7*1) = 0.0_dp
      kad(1:2*1) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,7
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*7+mb)
            p_bc = pbc((mc-1)*7+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*7+mb) = kbd((md-1)*7+mb) - ks_bd
            kbc((mc-1)*7+mb) = kbc((mc-1)*7+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_7_1_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_7_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(7*md_max), kbc(7*1), kad(2*md_max), kac(2*1), &
      pbd(7*md_max), pbc(7*1), pad(2*md_max), pac(2*1), prim(2*7*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:7*md_max) = 0.0_dp
      kbc(1:7*1) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,7
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*7+mb)
            p_bc = pbc((mc-1)*7+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*7+mb) = kbd((md-1)*7+mb) - ks_bd
            kbc((mc-1)*7+mb) = kbc((mc-1)*7+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_7_1
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(7*md_max), kbc(7*mc_max), kad(2*md_max), &
      kac(2*mc_max), pbd(7*md_max), pbc(7*mc_max), pad(2*md_max), &
      pac(2*mc_max), prim(2*7*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:7*md_max) = 0.0_dp
      kbc(1:7*mc_max) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,7
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*7+mb)
            p_bc = pbc((mc-1)*7+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*7+mb) = kbd((md-1)*7+mb) - ks_bd
            kbc((mc-1)*7+mb) = kbc((mc-1)*7+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_7
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_9_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(9*1), kbc(9*1), kad(2*1), kac(2*1), pbd(9*1), &
      pbc(9*1), pad(2*1), pac(2*1), prim(2*9*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:9*1) = 0.0_dp
      kbc(1:9*1) = 0.0_dp
      kad(1:2*1) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,9
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*9+mb)
            p_bc = pbc((mc-1)*9+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*9+mb) = kbd((md-1)*9+mb) - ks_bd
            kbc((mc-1)*9+mb) = kbc((mc-1)*9+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_9_1_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_9_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(9*md_max), kbc(9*1), kad(2*md_max), kac(2*1), &
      pbd(9*md_max), pbc(9*1), pad(2*md_max), pac(2*1), prim(2*9*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:9*md_max) = 0.0_dp
      kbc(1:9*1) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,9
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*9+mb)
            p_bc = pbc((mc-1)*9+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*9+mb) = kbd((md-1)*9+mb) - ks_bd
            kbc((mc-1)*9+mb) = kbc((mc-1)*9+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_9_1
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(9*md_max), kbc(9*mc_max), kad(2*md_max), &
      kac(2*mc_max), pbd(9*md_max), pbc(9*mc_max), pad(2*md_max), &
      pac(2*mc_max), prim(2*9*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:9*md_max) = 0.0_dp
      kbc(1:9*mc_max) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,9
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*9+mb)
            p_bc = pbc((mc-1)*9+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*9+mb) = kbd((md-1)*9+mb) - ks_bd
            kbc((mc-1)*9+mb) = kbc((mc-1)*9+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_9
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(10*md_max), kbc(10*mc_max), kad(2*md_max), &
      kac(2*mc_max), pbd(10*md_max), pbc(10*mc_max), pad(2*md_max), &
      pac(2*mc_max), prim(2*10*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:10*md_max) = 0.0_dp
      kbc(1:10*mc_max) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,10
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*10+mb)
            p_bc = pbc((mc-1)*10+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*10+mb) = kbd((md-1)*10+mb) - ks_bd
            kbc((mc-1)*10+mb) = kbc((mc-1)*10+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_10
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(11*md_max), kbc(11*mc_max), kad(2*md_max), &
      kac(2*mc_max), pbd(11*md_max), pbc(11*mc_max), pad(2*md_max), &
      pac(2*mc_max), prim(2*11*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:11*md_max) = 0.0_dp
      kbc(1:11*mc_max) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,11
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*11+mb)
            p_bc = pbc((mc-1)*11+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*11+mb) = kbd((md-1)*11+mb) - ks_bd
            kbc((mc-1)*11+mb) = kbc((mc-1)*11+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_11
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_2_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(15*md_max), kbc(15*mc_max), kad(2*md_max), &
      kac(2*mc_max), pbd(15*md_max), pbc(15*mc_max), pad(2*md_max), &
      pac(2*mc_max), prim(2*15*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:15*md_max) = 0.0_dp
      kbc(1:15*mc_max) = 0.0_dp
      kad(1:2*md_max) = 0.0_dp
      kac(1:2*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,15
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*15+mb)
            p_bc = pbc((mc-1)*15+mb)
            DO ma = 1,2
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*2+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*2+ma)
              kad((md-1)*2+ma) =  kad((md-1)*2+ma)-tmp*p_bc
              kac((mc-1)*2+ma) = kac((mc-1)*2+ma)-tmp*p_bd
            END DO
            kbd((md-1)*15+mb) = kbd((md-1)*15+mb) - ks_bd
            kbc((mc-1)*15+mb) = kbc((mc-1)*15+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_2_15
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_1_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*1), kad(3*1), kac(3*1), pbd(1*1), &
      pbc(1*1), pad(3*1), pac(3*1), prim(3*1*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:3*1) = 0.0_dp
      kac(1:3*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_1_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_1_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*2), kbc(1*1), kad(3*2), kac(3*1), pbd(1*2), &
      pbc(1*1), pad(3*2), pac(3*1), prim(3*1*1*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*2) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:3*2) = 0.0_dp
      kac(1:3*1) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_1_1_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_1_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*3), kbc(1*1), kad(3*3), kac(3*1), pbd(1*3), &
      pbc(1*1), pad(3*3), pac(3*1), prim(3*1*1*3), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*3) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:3*3) = 0.0_dp
      kac(1:3*1) = 0.0_dp
      p_index = 0
      DO md = 1,3
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_1_1_3
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_1_1_4(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*4), kbc(1*1), kad(3*4), kac(3*1), pbd(1*4), &
      pbc(1*1), pad(3*4), pac(3*1), prim(3*1*1*4), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*4) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:3*4) = 0.0_dp
      kac(1:3*1) = 0.0_dp
      p_index = 0
      DO md = 1,4
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_1_1_4
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_1_1_5(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*5), kbc(1*1), kad(3*5), kac(3*1), pbd(1*5), &
      pbc(1*1), pad(3*5), pac(3*1), prim(3*1*1*5), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*5) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:3*5) = 0.0_dp
      kac(1:3*1) = 0.0_dp
      p_index = 0
      DO md = 1,5
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_1_1_5
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_1_1_6(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*6), kbc(1*1), kad(3*6), kac(3*1), pbd(1*6), &
      pbc(1*1), pad(3*6), pac(3*1), prim(3*1*1*6), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*6) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:3*6) = 0.0_dp
      kac(1:3*1) = 0.0_dp
      p_index = 0
      DO md = 1,6
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_1_1_6
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*1), kad(3*md_max), kac(3*1), &
      pbd(1*md_max), pbc(1*1), pad(3*md_max), pac(3*1), prim(3*1*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:3*md_max) = 0.0_dp
      kac(1:3*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_1_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*2), kad(3*1), kac(3*2), pbd(1*1), &
      pbc(1*2), pad(3*1), pac(3*2), prim(3*1*2*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:3*1) = 0.0_dp
      kac(1:3*2) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_1_2_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_1_2_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*2), kbc(1*2), kad(3*2), kac(3*2), pbd(1*2), &
      pbc(1*2), pad(3*2), pac(3*2), prim(3*1*2*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*2) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:3*2) = 0.0_dp
      kac(1:3*2) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_1_2_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_1_2_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*3), kbc(1*2), kad(3*3), kac(3*2), pbd(1*3), &
      pbc(1*2), pad(3*3), pac(3*2), prim(3*1*2*3), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*3) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:3*3) = 0.0_dp
      kac(1:3*2) = 0.0_dp
      p_index = 0
      DO md = 1,3
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_1_2_3
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*2), kad(3*md_max), kac(3*2), &
      pbd(1*md_max), pbc(1*2), pad(3*md_max), pac(3*2), prim(3*1*2*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:3*md_max) = 0.0_dp
      kac(1:3*2) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_1_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_1_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*3), kad(3*1), kac(3*3), pbd(1*1), &
      pbc(1*3), pad(3*1), pac(3*3), prim(3*1*3*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*3) = 0.0_dp
      kad(1:3*1) = 0.0_dp
      kac(1:3*3) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,3
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_1_3_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_1_3_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*2), kbc(1*3), kad(3*2), kac(3*3), pbd(1*2), &
      pbc(1*3), pad(3*2), pac(3*3), prim(3*1*3*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*2) = 0.0_dp
      kbc(1:1*3) = 0.0_dp
      kad(1:3*2) = 0.0_dp
      kac(1:3*3) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,3
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_1_3_2
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*3), kad(3*md_max), kac(3*3), &
      pbd(1*md_max), pbc(1*3), pad(3*md_max), pac(3*3), prim(3*1*3*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*3) = 0.0_dp
      kad(1:3*md_max) = 0.0_dp
      kac(1:3*3) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,3
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_1_3
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_1_4_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*4), kad(3*1), kac(3*4), pbd(1*1), &
      pbc(1*4), pad(3*1), pac(3*4), prim(3*1*4*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*4) = 0.0_dp
      kad(1:3*1) = 0.0_dp
      kac(1:3*4) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,4
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_1_4_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*4), kad(3*md_max), kac(3*4), &
      pbd(1*md_max), pbc(1*4), pad(3*md_max), pac(3*4), prim(3*1*4*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*4) = 0.0_dp
      kad(1:3*md_max) = 0.0_dp
      kac(1:3*4) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,4
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_1_4
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_1_5_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*5), kad(3*1), kac(3*5), pbd(1*1), &
      pbc(1*5), pad(3*1), pac(3*5), prim(3*1*5*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*5) = 0.0_dp
      kad(1:3*1) = 0.0_dp
      kac(1:3*5) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,5
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_1_5_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_1_5(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*5), kad(3*md_max), kac(3*5), &
      pbd(1*md_max), pbc(1*5), pad(3*md_max), pac(3*5), prim(3*1*5*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*5) = 0.0_dp
      kad(1:3*md_max) = 0.0_dp
      kac(1:3*5) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,5
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_1_5
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_1_6_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*6), kad(3*1), kac(3*6), pbd(1*1), &
      pbc(1*6), pad(3*1), pac(3*6), prim(3*1*6*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*6) = 0.0_dp
      kad(1:3*1) = 0.0_dp
      kac(1:3*6) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,6
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_1_6_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_1_6(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*6), kad(3*md_max), kac(3*6), &
      pbd(1*md_max), pbc(1*6), pad(3*md_max), pac(3*6), prim(3*1*6*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*6) = 0.0_dp
      kad(1:3*md_max) = 0.0_dp
      kac(1:3*6) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,6
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_1_6
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*mc_max), kad(3*md_max), &
      kac(3*mc_max), pbd(1*md_max), pbc(1*mc_max), pad(3*md_max), &
      pac(3*mc_max), prim(3*1*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*mc_max) = 0.0_dp
      kad(1:3*md_max) = 0.0_dp
      kac(1:3*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_2_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*1), kbc(2*1), kad(3*1), kac(3*1), pbd(2*1), &
      pbc(2*1), pad(3*1), pac(3*1), prim(3*2*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*1) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:3*1) = 0.0_dp
      kac(1:3*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_2_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_2_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*2), kbc(2*1), kad(3*2), kac(3*1), pbd(2*2), &
      pbc(2*1), pad(3*2), pac(3*1), prim(3*2*1*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*2) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:3*2) = 0.0_dp
      kac(1:3*1) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_2_1_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_2_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*3), kbc(2*1), kad(3*3), kac(3*1), pbd(2*3), &
      pbc(2*1), pad(3*3), pac(3*1), prim(3*2*1*3), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*3) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:3*3) = 0.0_dp
      kac(1:3*1) = 0.0_dp
      p_index = 0
      DO md = 1,3
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_2_1_3
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*1), kad(3*md_max), kac(3*1), &
      pbd(2*md_max), pbc(2*1), pad(3*md_max), pac(3*1), prim(3*2*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:3*md_max) = 0.0_dp
      kac(1:3*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_2_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_2_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*1), kbc(2*2), kad(3*1), kac(3*2), pbd(2*1), &
      pbc(2*2), pad(3*1), pac(3*2), prim(3*2*2*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*1) = 0.0_dp
      kbc(1:2*2) = 0.0_dp
      kad(1:3*1) = 0.0_dp
      kac(1:3*2) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,2
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_2_2_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*2), kad(3*md_max), kac(3*2), &
      pbd(2*md_max), pbc(2*2), pad(3*md_max), pac(3*2), prim(3*2*2*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*2) = 0.0_dp
      kad(1:3*md_max) = 0.0_dp
      kac(1:3*2) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,2
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_2_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_2_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*1), kbc(2*3), kad(3*1), kac(3*3), pbd(2*1), &
      pbc(2*3), pad(3*1), pac(3*3), prim(3*2*3*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*1) = 0.0_dp
      kbc(1:2*3) = 0.0_dp
      kad(1:3*1) = 0.0_dp
      kac(1:3*3) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,3
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_2_3_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_2_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*3), kad(3*md_max), kac(3*3), &
      pbd(2*md_max), pbc(2*3), pad(3*md_max), pac(3*3), prim(3*2*3*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*3) = 0.0_dp
      kad(1:3*md_max) = 0.0_dp
      kac(1:3*3) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,3
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_2_3
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*mc_max), kad(3*md_max), &
      kac(3*mc_max), pbd(2*md_max), pbc(2*mc_max), pad(3*md_max), &
      pac(3*mc_max), prim(3*2*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*mc_max) = 0.0_dp
      kad(1:3*md_max) = 0.0_dp
      kac(1:3*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_3_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(3*1), kbc(3*1), kad(3*1), kac(3*1), pbd(3*1), &
      pbc(3*1), pad(3*1), pac(3*1), prim(3*3*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*1) = 0.0_dp
      kbc(1:3*1) = 0.0_dp
      kad(1:3*1) = 0.0_dp
      kac(1:3*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_3_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_3_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(3*2), kbc(3*1), kad(3*2), kac(3*1), pbd(3*2), &
      pbc(3*1), pad(3*2), pac(3*1), prim(3*3*1*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*2) = 0.0_dp
      kbc(1:3*1) = 0.0_dp
      kad(1:3*2) = 0.0_dp
      kac(1:3*1) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,1
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_3_1_2
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(3*md_max), kbc(3*1), kad(3*md_max), kac(3*1), &
      pbd(3*md_max), pbc(3*1), pad(3*md_max), pac(3*1), prim(3*3*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*md_max) = 0.0_dp
      kbc(1:3*1) = 0.0_dp
      kad(1:3*md_max) = 0.0_dp
      kac(1:3*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_3_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_3_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(3*1), kbc(3*2), kad(3*1), kac(3*2), pbd(3*1), &
      pbc(3*2), pad(3*1), pac(3*2), prim(3*3*2*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*1) = 0.0_dp
      kbc(1:3*2) = 0.0_dp
      kad(1:3*1) = 0.0_dp
      kac(1:3*2) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,2
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_3_2_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_3_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(3*md_max), kbc(3*2), kad(3*md_max), kac(3*2), &
      pbd(3*md_max), pbc(3*2), pad(3*md_max), pac(3*2), prim(3*3*2*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*md_max) = 0.0_dp
      kbc(1:3*2) = 0.0_dp
      kad(1:3*md_max) = 0.0_dp
      kac(1:3*2) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,2
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_3_2
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(3*md_max), kbc(3*mc_max), kad(3*md_max), &
      kac(3*mc_max), pbd(3*md_max), pbc(3*mc_max), pad(3*md_max), &
      pac(3*mc_max), prim(3*3*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*md_max) = 0.0_dp
      kbc(1:3*mc_max) = 0.0_dp
      kad(1:3*md_max) = 0.0_dp
      kac(1:3*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_3
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_4_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(4*1), kbc(4*1), kad(3*1), kac(3*1), pbd(4*1), &
      pbc(4*1), pad(3*1), pac(3*1), prim(3*4*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*1) = 0.0_dp
      kbc(1:4*1) = 0.0_dp
      kad(1:3*1) = 0.0_dp
      kac(1:3*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_4_1_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(4*md_max), kbc(4*1), kad(3*md_max), kac(3*1), &
      pbd(4*md_max), pbc(4*1), pad(3*md_max), pac(3*1), prim(3*4*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*md_max) = 0.0_dp
      kbc(1:4*1) = 0.0_dp
      kad(1:3*md_max) = 0.0_dp
      kac(1:3*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_4_1
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(4*md_max), kbc(4*mc_max), kad(3*md_max), &
      kac(3*mc_max), pbd(4*md_max), pbc(4*mc_max), pad(3*md_max), &
      pac(3*mc_max), prim(3*4*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*md_max) = 0.0_dp
      kbc(1:4*mc_max) = 0.0_dp
      kad(1:3*md_max) = 0.0_dp
      kac(1:3*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_4
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_5_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(5*1), kbc(5*1), kad(3*1), kac(3*1), pbd(5*1), &
      pbc(5*1), pad(3*1), pac(3*1), prim(3*5*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:5*1) = 0.0_dp
      kbc(1:5*1) = 0.0_dp
      kad(1:3*1) = 0.0_dp
      kac(1:3*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,5
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*5+mb)
            p_bc = pbc((mc-1)*5+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*5+mb) = kbd((md-1)*5+mb) - ks_bd
            kbc((mc-1)*5+mb) = kbc((mc-1)*5+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_5_1_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_5_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(5*md_max), kbc(5*1), kad(3*md_max), kac(3*1), &
      pbd(5*md_max), pbc(5*1), pad(3*md_max), pac(3*1), prim(3*5*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:5*md_max) = 0.0_dp
      kbc(1:5*1) = 0.0_dp
      kad(1:3*md_max) = 0.0_dp
      kac(1:3*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,5
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*5+mb)
            p_bc = pbc((mc-1)*5+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*5+mb) = kbd((md-1)*5+mb) - ks_bd
            kbc((mc-1)*5+mb) = kbc((mc-1)*5+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_5_1
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(5*md_max), kbc(5*mc_max), kad(3*md_max), &
      kac(3*mc_max), pbd(5*md_max), pbc(5*mc_max), pad(3*md_max), &
      pac(3*mc_max), prim(3*5*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:5*md_max) = 0.0_dp
      kbc(1:5*mc_max) = 0.0_dp
      kad(1:3*md_max) = 0.0_dp
      kac(1:3*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,5
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*5+mb)
            p_bc = pbc((mc-1)*5+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*5+mb) = kbd((md-1)*5+mb) - ks_bd
            kbc((mc-1)*5+mb) = kbc((mc-1)*5+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_5
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_6_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(6*1), kbc(6*1), kad(3*1), kac(3*1), pbd(6*1), &
      pbc(6*1), pad(3*1), pac(3*1), prim(3*6*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:6*1) = 0.0_dp
      kbc(1:6*1) = 0.0_dp
      kad(1:3*1) = 0.0_dp
      kac(1:3*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,6
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*6+mb)
            p_bc = pbc((mc-1)*6+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*6+mb) = kbd((md-1)*6+mb) - ks_bd
            kbc((mc-1)*6+mb) = kbc((mc-1)*6+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_6_1_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_6_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(6*md_max), kbc(6*1), kad(3*md_max), kac(3*1), &
      pbd(6*md_max), pbc(6*1), pad(3*md_max), pac(3*1), prim(3*6*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:6*md_max) = 0.0_dp
      kbc(1:6*1) = 0.0_dp
      kad(1:3*md_max) = 0.0_dp
      kac(1:3*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,6
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*6+mb)
            p_bc = pbc((mc-1)*6+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*6+mb) = kbd((md-1)*6+mb) - ks_bd
            kbc((mc-1)*6+mb) = kbc((mc-1)*6+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_6_1
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(6*md_max), kbc(6*mc_max), kad(3*md_max), &
      kac(3*mc_max), pbd(6*md_max), pbc(6*mc_max), pad(3*md_max), &
      pac(3*mc_max), prim(3*6*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:6*md_max) = 0.0_dp
      kbc(1:6*mc_max) = 0.0_dp
      kad(1:3*md_max) = 0.0_dp
      kac(1:3*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,6
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*6+mb)
            p_bc = pbc((mc-1)*6+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*6+mb) = kbd((md-1)*6+mb) - ks_bd
            kbc((mc-1)*6+mb) = kbc((mc-1)*6+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_6
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(7*md_max), kbc(7*mc_max), kad(3*md_max), &
      kac(3*mc_max), pbd(7*md_max), pbc(7*mc_max), pad(3*md_max), &
      pac(3*mc_max), prim(3*7*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:7*md_max) = 0.0_dp
      kbc(1:7*mc_max) = 0.0_dp
      kad(1:3*md_max) = 0.0_dp
      kac(1:3*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,7
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*7+mb)
            p_bc = pbc((mc-1)*7+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*7+mb) = kbd((md-1)*7+mb) - ks_bd
            kbc((mc-1)*7+mb) = kbc((mc-1)*7+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_7
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(9*md_max), kbc(9*mc_max), kad(3*md_max), &
      kac(3*mc_max), pbd(9*md_max), pbc(9*mc_max), pad(3*md_max), &
      pac(3*mc_max), prim(3*9*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:9*md_max) = 0.0_dp
      kbc(1:9*mc_max) = 0.0_dp
      kad(1:3*md_max) = 0.0_dp
      kac(1:3*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,9
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*9+mb)
            p_bc = pbc((mc-1)*9+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*9+mb) = kbd((md-1)*9+mb) - ks_bd
            kbc((mc-1)*9+mb) = kbc((mc-1)*9+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_9
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(10*md_max), kbc(10*mc_max), kad(3*md_max), &
      kac(3*mc_max), pbd(10*md_max), pbc(10*mc_max), pad(3*md_max), &
      pac(3*mc_max), prim(3*10*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:10*md_max) = 0.0_dp
      kbc(1:10*mc_max) = 0.0_dp
      kad(1:3*md_max) = 0.0_dp
      kac(1:3*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,10
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*10+mb)
            p_bc = pbc((mc-1)*10+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*10+mb) = kbd((md-1)*10+mb) - ks_bd
            kbc((mc-1)*10+mb) = kbc((mc-1)*10+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_10
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(11*md_max), kbc(11*mc_max), kad(3*md_max), &
      kac(3*mc_max), pbd(11*md_max), pbc(11*mc_max), pad(3*md_max), &
      pac(3*mc_max), prim(3*11*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:11*md_max) = 0.0_dp
      kbc(1:11*mc_max) = 0.0_dp
      kad(1:3*md_max) = 0.0_dp
      kac(1:3*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,11
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*11+mb)
            p_bc = pbc((mc-1)*11+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*11+mb) = kbd((md-1)*11+mb) - ks_bd
            kbc((mc-1)*11+mb) = kbc((mc-1)*11+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_11
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_3_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(15*md_max), kbc(15*mc_max), kad(3*md_max), &
      kac(3*mc_max), pbd(15*md_max), pbc(15*mc_max), pad(3*md_max), &
      pac(3*mc_max), prim(3*15*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:15*md_max) = 0.0_dp
      kbc(1:15*mc_max) = 0.0_dp
      kad(1:3*md_max) = 0.0_dp
      kac(1:3*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,15
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*15+mb)
            p_bc = pbc((mc-1)*15+mb)
            DO ma = 1,3
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*3+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*3+ma)
              kad((md-1)*3+ma) =  kad((md-1)*3+ma)-tmp*p_bc
              kac((mc-1)*3+ma) = kac((mc-1)*3+ma)-tmp*p_bd
            END DO
            kbd((md-1)*15+mb) = kbd((md-1)*15+mb) - ks_bd
            kbc((mc-1)*15+mb) = kbc((mc-1)*15+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_3_15
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_1_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*1), kad(4*1), kac(4*1), pbd(1*1), &
      pbc(1*1), pad(4*1), pac(4*1), prim(4*1*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:4*1) = 0.0_dp
      kac(1:4*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_1_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_1_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*2), kbc(1*1), kad(4*2), kac(4*1), pbd(1*2), &
      pbc(1*1), pad(4*2), pac(4*1), prim(4*1*1*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*2) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:4*2) = 0.0_dp
      kac(1:4*1) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_1_1_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_1_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*3), kbc(1*1), kad(4*3), kac(4*1), pbd(1*3), &
      pbc(1*1), pad(4*3), pac(4*1), prim(4*1*1*3), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*3) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:4*3) = 0.0_dp
      kac(1:4*1) = 0.0_dp
      p_index = 0
      DO md = 1,3
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_1_1_3
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_1_1_4(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*4), kbc(1*1), kad(4*4), kac(4*1), pbd(1*4), &
      pbc(1*1), pad(4*4), pac(4*1), prim(4*1*1*4), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*4) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:4*4) = 0.0_dp
      kac(1:4*1) = 0.0_dp
      p_index = 0
      DO md = 1,4
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_1_1_4
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*1), kad(4*md_max), kac(4*1), &
      pbd(1*md_max), pbc(1*1), pad(4*md_max), pac(4*1), prim(4*1*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:4*md_max) = 0.0_dp
      kac(1:4*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_1_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*2), kad(4*1), kac(4*2), pbd(1*1), &
      pbc(1*2), pad(4*1), pac(4*2), prim(4*1*2*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:4*1) = 0.0_dp
      kac(1:4*2) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_1_2_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_1_2_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*2), kbc(1*2), kad(4*2), kac(4*2), pbd(1*2), &
      pbc(1*2), pad(4*2), pac(4*2), prim(4*1*2*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*2) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:4*2) = 0.0_dp
      kac(1:4*2) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_1_2_2
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*2), kad(4*md_max), kac(4*2), &
      pbd(1*md_max), pbc(1*2), pad(4*md_max), pac(4*2), prim(4*1*2*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:4*md_max) = 0.0_dp
      kac(1:4*2) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_1_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_1_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*3), kad(4*1), kac(4*3), pbd(1*1), &
      pbc(1*3), pad(4*1), pac(4*3), prim(4*1*3*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*3) = 0.0_dp
      kad(1:4*1) = 0.0_dp
      kac(1:4*3) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,3
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_1_3_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*3), kad(4*md_max), kac(4*3), &
      pbd(1*md_max), pbc(1*3), pad(4*md_max), pac(4*3), prim(4*1*3*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*3) = 0.0_dp
      kad(1:4*md_max) = 0.0_dp
      kac(1:4*3) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,3
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_1_3
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_1_4_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*4), kad(4*1), kac(4*4), pbd(1*1), &
      pbc(1*4), pad(4*1), pac(4*4), prim(4*1*4*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*4) = 0.0_dp
      kad(1:4*1) = 0.0_dp
      kac(1:4*4) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,4
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_1_4_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_1_4(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*4), kad(4*md_max), kac(4*4), &
      pbd(1*md_max), pbc(1*4), pad(4*md_max), pac(4*4), prim(4*1*4*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*4) = 0.0_dp
      kad(1:4*md_max) = 0.0_dp
      kac(1:4*4) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,4
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_1_4
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*mc_max), kad(4*md_max), &
      kac(4*mc_max), pbd(1*md_max), pbc(1*mc_max), pad(4*md_max), &
      pac(4*mc_max), prim(4*1*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*mc_max) = 0.0_dp
      kad(1:4*md_max) = 0.0_dp
      kac(1:4*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_2_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*1), kbc(2*1), kad(4*1), kac(4*1), pbd(2*1), &
      pbc(2*1), pad(4*1), pac(4*1), prim(4*2*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*1) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:4*1) = 0.0_dp
      kac(1:4*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_2_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_2_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*2), kbc(2*1), kad(4*2), kac(4*1), pbd(2*2), &
      pbc(2*1), pad(4*2), pac(4*1), prim(4*2*1*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*2) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:4*2) = 0.0_dp
      kac(1:4*1) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_2_1_2
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*1), kad(4*md_max), kac(4*1), &
      pbd(2*md_max), pbc(2*1), pad(4*md_max), pac(4*1), prim(4*2*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:4*md_max) = 0.0_dp
      kac(1:4*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_2_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_2_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*1), kbc(2*2), kad(4*1), kac(4*2), pbd(2*1), &
      pbc(2*2), pad(4*1), pac(4*2), prim(4*2*2*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*1) = 0.0_dp
      kbc(1:2*2) = 0.0_dp
      kad(1:4*1) = 0.0_dp
      kac(1:4*2) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,2
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_2_2_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_2_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*2), kad(4*md_max), kac(4*2), &
      pbd(2*md_max), pbc(2*2), pad(4*md_max), pac(4*2), prim(4*2*2*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*2) = 0.0_dp
      kad(1:4*md_max) = 0.0_dp
      kac(1:4*2) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,2
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_2_2
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*mc_max), kad(4*md_max), &
      kac(4*mc_max), pbd(2*md_max), pbc(2*mc_max), pad(4*md_max), &
      pac(4*mc_max), prim(4*2*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*mc_max) = 0.0_dp
      kad(1:4*md_max) = 0.0_dp
      kac(1:4*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_3_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(3*1), kbc(3*1), kad(4*1), kac(4*1), pbd(3*1), &
      pbc(3*1), pad(4*1), pac(4*1), prim(4*3*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*1) = 0.0_dp
      kbc(1:3*1) = 0.0_dp
      kad(1:4*1) = 0.0_dp
      kac(1:4*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_3_1_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(3*md_max), kbc(3*1), kad(4*md_max), kac(4*1), &
      pbd(3*md_max), pbc(3*1), pad(4*md_max), pac(4*1), prim(4*3*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*md_max) = 0.0_dp
      kbc(1:3*1) = 0.0_dp
      kad(1:4*md_max) = 0.0_dp
      kac(1:4*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_3_1
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(3*md_max), kbc(3*mc_max), kad(4*md_max), &
      kac(4*mc_max), pbd(3*md_max), pbc(3*mc_max), pad(4*md_max), &
      pac(4*mc_max), prim(4*3*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*md_max) = 0.0_dp
      kbc(1:3*mc_max) = 0.0_dp
      kad(1:4*md_max) = 0.0_dp
      kac(1:4*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_3
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_4_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(4*1), kbc(4*1), kad(4*1), kac(4*1), pbd(4*1), &
      pbc(4*1), pad(4*1), pac(4*1), prim(4*4*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*1) = 0.0_dp
      kbc(1:4*1) = 0.0_dp
      kad(1:4*1) = 0.0_dp
      kac(1:4*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_4_1_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_4_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(4*md_max), kbc(4*1), kad(4*md_max), kac(4*1), &
      pbd(4*md_max), pbc(4*1), pad(4*md_max), pac(4*1), prim(4*4*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*md_max) = 0.0_dp
      kbc(1:4*1) = 0.0_dp
      kad(1:4*md_max) = 0.0_dp
      kac(1:4*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_4_1
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(4*md_max), kbc(4*mc_max), kad(4*md_max), &
      kac(4*mc_max), pbd(4*md_max), pbc(4*mc_max), pad(4*md_max), &
      pac(4*mc_max), prim(4*4*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*md_max) = 0.0_dp
      kbc(1:4*mc_max) = 0.0_dp
      kad(1:4*md_max) = 0.0_dp
      kac(1:4*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_4
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(5*md_max), kbc(5*mc_max), kad(4*md_max), &
      kac(4*mc_max), pbd(5*md_max), pbc(5*mc_max), pad(4*md_max), &
      pac(4*mc_max), prim(4*5*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:5*md_max) = 0.0_dp
      kbc(1:5*mc_max) = 0.0_dp
      kad(1:4*md_max) = 0.0_dp
      kac(1:4*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,5
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*5+mb)
            p_bc = pbc((mc-1)*5+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*5+mb) = kbd((md-1)*5+mb) - ks_bd
            kbc((mc-1)*5+mb) = kbc((mc-1)*5+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_5
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(6*md_max), kbc(6*mc_max), kad(4*md_max), &
      kac(4*mc_max), pbd(6*md_max), pbc(6*mc_max), pad(4*md_max), &
      pac(4*mc_max), prim(4*6*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:6*md_max) = 0.0_dp
      kbc(1:6*mc_max) = 0.0_dp
      kad(1:4*md_max) = 0.0_dp
      kac(1:4*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,6
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*6+mb)
            p_bc = pbc((mc-1)*6+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*6+mb) = kbd((md-1)*6+mb) - ks_bd
            kbc((mc-1)*6+mb) = kbc((mc-1)*6+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_6
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(7*md_max), kbc(7*mc_max), kad(4*md_max), &
      kac(4*mc_max), pbd(7*md_max), pbc(7*mc_max), pad(4*md_max), &
      pac(4*mc_max), prim(4*7*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:7*md_max) = 0.0_dp
      kbc(1:7*mc_max) = 0.0_dp
      kad(1:4*md_max) = 0.0_dp
      kac(1:4*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,7
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*7+mb)
            p_bc = pbc((mc-1)*7+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*7+mb) = kbd((md-1)*7+mb) - ks_bd
            kbc((mc-1)*7+mb) = kbc((mc-1)*7+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_7
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(9*md_max), kbc(9*mc_max), kad(4*md_max), &
      kac(4*mc_max), pbd(9*md_max), pbc(9*mc_max), pad(4*md_max), &
      pac(4*mc_max), prim(4*9*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:9*md_max) = 0.0_dp
      kbc(1:9*mc_max) = 0.0_dp
      kad(1:4*md_max) = 0.0_dp
      kac(1:4*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,9
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*9+mb)
            p_bc = pbc((mc-1)*9+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*9+mb) = kbd((md-1)*9+mb) - ks_bd
            kbc((mc-1)*9+mb) = kbc((mc-1)*9+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_9
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(10*md_max), kbc(10*mc_max), kad(4*md_max), &
      kac(4*mc_max), pbd(10*md_max), pbc(10*mc_max), pad(4*md_max), &
      pac(4*mc_max), prim(4*10*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:10*md_max) = 0.0_dp
      kbc(1:10*mc_max) = 0.0_dp
      kad(1:4*md_max) = 0.0_dp
      kac(1:4*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,10
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*10+mb)
            p_bc = pbc((mc-1)*10+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*10+mb) = kbd((md-1)*10+mb) - ks_bd
            kbc((mc-1)*10+mb) = kbc((mc-1)*10+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_10
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(11*md_max), kbc(11*mc_max), kad(4*md_max), &
      kac(4*mc_max), pbd(11*md_max), pbc(11*mc_max), pad(4*md_max), &
      pac(4*mc_max), prim(4*11*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:11*md_max) = 0.0_dp
      kbc(1:11*mc_max) = 0.0_dp
      kad(1:4*md_max) = 0.0_dp
      kac(1:4*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,11
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*11+mb)
            p_bc = pbc((mc-1)*11+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*11+mb) = kbd((md-1)*11+mb) - ks_bd
            kbc((mc-1)*11+mb) = kbc((mc-1)*11+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_11
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_4_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(15*md_max), kbc(15*mc_max), kad(4*md_max), &
      kac(4*mc_max), pbd(15*md_max), pbc(15*mc_max), pad(4*md_max), &
      pac(4*mc_max), prim(4*15*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:15*md_max) = 0.0_dp
      kbc(1:15*mc_max) = 0.0_dp
      kad(1:4*md_max) = 0.0_dp
      kac(1:4*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,15
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*15+mb)
            p_bc = pbc((mc-1)*15+mb)
            DO ma = 1,4
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*4+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*4+ma)
              kad((md-1)*4+ma) =  kad((md-1)*4+ma)-tmp*p_bc
              kac((mc-1)*4+ma) = kac((mc-1)*4+ma)-tmp*p_bd
            END DO
            kbd((md-1)*15+mb) = kbd((md-1)*15+mb) - ks_bd
            kbc((mc-1)*15+mb) = kbc((mc-1)*15+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_4_15
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_5_1_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*1), kad(5*1), kac(5*1), pbd(1*1), &
      pbc(1*1), pad(5*1), pac(5*1), prim(5*1*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:5*1) = 0.0_dp
      kac(1:5*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,5
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*5+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*5+ma)
              kad((md-1)*5+ma) =  kad((md-1)*5+ma)-tmp*p_bc
              kac((mc-1)*5+ma) = kac((mc-1)*5+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_5_1_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_5_1_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*2), kbc(1*1), kad(5*2), kac(5*1), pbd(1*2), &
      pbc(1*1), pad(5*2), pac(5*1), prim(5*1*1*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*2) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:5*2) = 0.0_dp
      kac(1:5*1) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,5
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*5+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*5+ma)
              kad((md-1)*5+ma) =  kad((md-1)*5+ma)-tmp*p_bc
              kac((mc-1)*5+ma) = kac((mc-1)*5+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_5_1_1_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_5_1_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*3), kbc(1*1), kad(5*3), kac(5*1), pbd(1*3), &
      pbc(1*1), pad(5*3), pac(5*1), prim(5*1*1*3), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*3) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:5*3) = 0.0_dp
      kac(1:5*1) = 0.0_dp
      p_index = 0
      DO md = 1,3
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,5
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*5+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*5+ma)
              kad((md-1)*5+ma) =  kad((md-1)*5+ma)-tmp*p_bc
              kac((mc-1)*5+ma) = kac((mc-1)*5+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_5_1_1_3
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_5_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*1), kad(5*md_max), kac(5*1), &
      pbd(1*md_max), pbc(1*1), pad(5*md_max), pac(5*1), prim(5*1*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:5*md_max) = 0.0_dp
      kac(1:5*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,5
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*5+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*5+ma)
              kad((md-1)*5+ma) =  kad((md-1)*5+ma)-tmp*p_bc
              kac((mc-1)*5+ma) = kac((mc-1)*5+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_5_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_5_1_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*2), kad(5*1), kac(5*2), pbd(1*1), &
      pbc(1*2), pad(5*1), pac(5*2), prim(5*1*2*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:5*1) = 0.0_dp
      kac(1:5*2) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,5
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*5+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*5+ma)
              kad((md-1)*5+ma) =  kad((md-1)*5+ma)-tmp*p_bc
              kac((mc-1)*5+ma) = kac((mc-1)*5+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_5_1_2_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_5_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*2), kad(5*md_max), kac(5*2), &
      pbd(1*md_max), pbc(1*2), pad(5*md_max), pac(5*2), prim(5*1*2*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:5*md_max) = 0.0_dp
      kac(1:5*2) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,5
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*5+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*5+ma)
              kad((md-1)*5+ma) =  kad((md-1)*5+ma)-tmp*p_bc
              kac((mc-1)*5+ma) = kac((mc-1)*5+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_5_1_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_5_1_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*3), kad(5*1), kac(5*3), pbd(1*1), &
      pbc(1*3), pad(5*1), pac(5*3), prim(5*1*3*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*3) = 0.0_dp
      kad(1:5*1) = 0.0_dp
      kac(1:5*3) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,3
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,5
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*5+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*5+ma)
              kad((md-1)*5+ma) =  kad((md-1)*5+ma)-tmp*p_bc
              kac((mc-1)*5+ma) = kac((mc-1)*5+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_5_1_3_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_5_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*3), kad(5*md_max), kac(5*3), &
      pbd(1*md_max), pbc(1*3), pad(5*md_max), pac(5*3), prim(5*1*3*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*3) = 0.0_dp
      kad(1:5*md_max) = 0.0_dp
      kac(1:5*3) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,3
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,5
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*5+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*5+ma)
              kad((md-1)*5+ma) =  kad((md-1)*5+ma)-tmp*p_bc
              kac((mc-1)*5+ma) = kac((mc-1)*5+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_5_1_3
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_5_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*mc_max), kad(5*md_max), &
      kac(5*mc_max), pbd(1*md_max), pbc(1*mc_max), pad(5*md_max), &
      pac(5*mc_max), prim(5*1*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*mc_max) = 0.0_dp
      kad(1:5*md_max) = 0.0_dp
      kac(1:5*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,5
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*5+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*5+ma)
              kad((md-1)*5+ma) =  kad((md-1)*5+ma)-tmp*p_bc
              kac((mc-1)*5+ma) = kac((mc-1)*5+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_5_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_5_2_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*1), kbc(2*1), kad(5*1), kac(5*1), pbd(2*1), &
      pbc(2*1), pad(5*1), pac(5*1), prim(5*2*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*1) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:5*1) = 0.0_dp
      kac(1:5*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,5
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*5+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*5+ma)
              kad((md-1)*5+ma) =  kad((md-1)*5+ma)-tmp*p_bc
              kac((mc-1)*5+ma) = kac((mc-1)*5+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_5_2_1_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_5_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*1), kad(5*md_max), kac(5*1), &
      pbd(2*md_max), pbc(2*1), pad(5*md_max), pac(5*1), prim(5*2*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:5*md_max) = 0.0_dp
      kac(1:5*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,5
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*5+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*5+ma)
              kad((md-1)*5+ma) =  kad((md-1)*5+ma)-tmp*p_bc
              kac((mc-1)*5+ma) = kac((mc-1)*5+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_5_2_1
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_5_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*mc_max), kad(5*md_max), &
      kac(5*mc_max), pbd(2*md_max), pbc(2*mc_max), pad(5*md_max), &
      pac(5*mc_max), prim(5*2*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*mc_max) = 0.0_dp
      kad(1:5*md_max) = 0.0_dp
      kac(1:5*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,5
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*5+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*5+ma)
              kad((md-1)*5+ma) =  kad((md-1)*5+ma)-tmp*p_bc
              kac((mc-1)*5+ma) = kac((mc-1)*5+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_5_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_5_3_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(3*1), kbc(3*1), kad(5*1), kac(5*1), pbd(3*1), &
      pbc(3*1), pad(5*1), pac(5*1), prim(5*3*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*1) = 0.0_dp
      kbc(1:3*1) = 0.0_dp
      kad(1:5*1) = 0.0_dp
      kac(1:5*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,5
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*5+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*5+ma)
              kad((md-1)*5+ma) =  kad((md-1)*5+ma)-tmp*p_bc
              kac((mc-1)*5+ma) = kac((mc-1)*5+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_5_3_1_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_5_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(3*md_max), kbc(3*1), kad(5*md_max), kac(5*1), &
      pbd(3*md_max), pbc(3*1), pad(5*md_max), pac(5*1), prim(5*3*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*md_max) = 0.0_dp
      kbc(1:3*1) = 0.0_dp
      kad(1:5*md_max) = 0.0_dp
      kac(1:5*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,5
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*5+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*5+ma)
              kad((md-1)*5+ma) =  kad((md-1)*5+ma)-tmp*p_bc
              kac((mc-1)*5+ma) = kac((mc-1)*5+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_5_3_1
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_5_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(3*md_max), kbc(3*mc_max), kad(5*md_max), &
      kac(5*mc_max), pbd(3*md_max), pbc(3*mc_max), pad(5*md_max), &
      pac(5*mc_max), prim(5*3*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*md_max) = 0.0_dp
      kbc(1:3*mc_max) = 0.0_dp
      kad(1:5*md_max) = 0.0_dp
      kac(1:5*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,5
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*5+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*5+ma)
              kad((md-1)*5+ma) =  kad((md-1)*5+ma)-tmp*p_bc
              kac((mc-1)*5+ma) = kac((mc-1)*5+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_5_3
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_5_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(4*md_max), kbc(4*mc_max), kad(5*md_max), &
      kac(5*mc_max), pbd(4*md_max), pbc(4*mc_max), pad(5*md_max), &
      pac(5*mc_max), prim(5*4*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*md_max) = 0.0_dp
      kbc(1:4*mc_max) = 0.0_dp
      kad(1:5*md_max) = 0.0_dp
      kac(1:5*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,5
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*5+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*5+ma)
              kad((md-1)*5+ma) =  kad((md-1)*5+ma)-tmp*p_bc
              kac((mc-1)*5+ma) = kac((mc-1)*5+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_5_4
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_5_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(5*md_max), kbc(5*mc_max), kad(5*md_max), &
      kac(5*mc_max), pbd(5*md_max), pbc(5*mc_max), pad(5*md_max), &
      pac(5*mc_max), prim(5*5*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:5*md_max) = 0.0_dp
      kbc(1:5*mc_max) = 0.0_dp
      kad(1:5*md_max) = 0.0_dp
      kac(1:5*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,5
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*5+mb)
            p_bc = pbc((mc-1)*5+mb)
            DO ma = 1,5
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*5+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*5+ma)
              kad((md-1)*5+ma) =  kad((md-1)*5+ma)-tmp*p_bc
              kac((mc-1)*5+ma) = kac((mc-1)*5+ma)-tmp*p_bd
            END DO
            kbd((md-1)*5+mb) = kbd((md-1)*5+mb) - ks_bd
            kbc((mc-1)*5+mb) = kbc((mc-1)*5+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_5_5
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_5_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(6*md_max), kbc(6*mc_max), kad(5*md_max), &
      kac(5*mc_max), pbd(6*md_max), pbc(6*mc_max), pad(5*md_max), &
      pac(5*mc_max), prim(5*6*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:6*md_max) = 0.0_dp
      kbc(1:6*mc_max) = 0.0_dp
      kad(1:5*md_max) = 0.0_dp
      kac(1:5*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,6
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*6+mb)
            p_bc = pbc((mc-1)*6+mb)
            DO ma = 1,5
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*5+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*5+ma)
              kad((md-1)*5+ma) =  kad((md-1)*5+ma)-tmp*p_bc
              kac((mc-1)*5+ma) = kac((mc-1)*5+ma)-tmp*p_bd
            END DO
            kbd((md-1)*6+mb) = kbd((md-1)*6+mb) - ks_bd
            kbc((mc-1)*6+mb) = kbc((mc-1)*6+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_5_6
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_5_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(7*md_max), kbc(7*mc_max), kad(5*md_max), &
      kac(5*mc_max), pbd(7*md_max), pbc(7*mc_max), pad(5*md_max), &
      pac(5*mc_max), prim(5*7*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:7*md_max) = 0.0_dp
      kbc(1:7*mc_max) = 0.0_dp
      kad(1:5*md_max) = 0.0_dp
      kac(1:5*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,7
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*7+mb)
            p_bc = pbc((mc-1)*7+mb)
            DO ma = 1,5
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*5+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*5+ma)
              kad((md-1)*5+ma) =  kad((md-1)*5+ma)-tmp*p_bc
              kac((mc-1)*5+ma) = kac((mc-1)*5+ma)-tmp*p_bd
            END DO
            kbd((md-1)*7+mb) = kbd((md-1)*7+mb) - ks_bd
            kbc((mc-1)*7+mb) = kbc((mc-1)*7+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_5_7
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_5_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(9*md_max), kbc(9*mc_max), kad(5*md_max), &
      kac(5*mc_max), pbd(9*md_max), pbc(9*mc_max), pad(5*md_max), &
      pac(5*mc_max), prim(5*9*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:9*md_max) = 0.0_dp
      kbc(1:9*mc_max) = 0.0_dp
      kad(1:5*md_max) = 0.0_dp
      kac(1:5*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,9
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*9+mb)
            p_bc = pbc((mc-1)*9+mb)
            DO ma = 1,5
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*5+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*5+ma)
              kad((md-1)*5+ma) =  kad((md-1)*5+ma)-tmp*p_bc
              kac((mc-1)*5+ma) = kac((mc-1)*5+ma)-tmp*p_bd
            END DO
            kbd((md-1)*9+mb) = kbd((md-1)*9+mb) - ks_bd
            kbc((mc-1)*9+mb) = kbc((mc-1)*9+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_5_9
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_5_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(10*md_max), kbc(10*mc_max), kad(5*md_max), &
      kac(5*mc_max), pbd(10*md_max), pbc(10*mc_max), pad(5*md_max), &
      pac(5*mc_max), prim(5*10*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:10*md_max) = 0.0_dp
      kbc(1:10*mc_max) = 0.0_dp
      kad(1:5*md_max) = 0.0_dp
      kac(1:5*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,10
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*10+mb)
            p_bc = pbc((mc-1)*10+mb)
            DO ma = 1,5
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*5+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*5+ma)
              kad((md-1)*5+ma) =  kad((md-1)*5+ma)-tmp*p_bc
              kac((mc-1)*5+ma) = kac((mc-1)*5+ma)-tmp*p_bd
            END DO
            kbd((md-1)*10+mb) = kbd((md-1)*10+mb) - ks_bd
            kbc((mc-1)*10+mb) = kbc((mc-1)*10+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_5_10
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_5_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(11*md_max), kbc(11*mc_max), kad(5*md_max), &
      kac(5*mc_max), pbd(11*md_max), pbc(11*mc_max), pad(5*md_max), &
      pac(5*mc_max), prim(5*11*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:11*md_max) = 0.0_dp
      kbc(1:11*mc_max) = 0.0_dp
      kad(1:5*md_max) = 0.0_dp
      kac(1:5*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,11
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*11+mb)
            p_bc = pbc((mc-1)*11+mb)
            DO ma = 1,5
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*5+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*5+ma)
              kad((md-1)*5+ma) =  kad((md-1)*5+ma)-tmp*p_bc
              kac((mc-1)*5+ma) = kac((mc-1)*5+ma)-tmp*p_bd
            END DO
            kbd((md-1)*11+mb) = kbd((md-1)*11+mb) - ks_bd
            kbc((mc-1)*11+mb) = kbc((mc-1)*11+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_5_11
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_5_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(15*md_max), kbc(15*mc_max), kad(5*md_max), &
      kac(5*mc_max), pbd(15*md_max), pbc(15*mc_max), pad(5*md_max), &
      pac(5*mc_max), prim(5*15*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:15*md_max) = 0.0_dp
      kbc(1:15*mc_max) = 0.0_dp
      kad(1:5*md_max) = 0.0_dp
      kac(1:5*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,15
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*15+mb)
            p_bc = pbc((mc-1)*15+mb)
            DO ma = 1,5
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*5+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*5+ma)
              kad((md-1)*5+ma) =  kad((md-1)*5+ma)-tmp*p_bc
              kac((mc-1)*5+ma) = kac((mc-1)*5+ma)-tmp*p_bd
            END DO
            kbd((md-1)*15+mb) = kbd((md-1)*15+mb) - ks_bd
            kbc((mc-1)*15+mb) = kbc((mc-1)*15+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_5_15
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_6_1_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*1), kad(6*1), kac(6*1), pbd(1*1), &
      pbc(1*1), pad(6*1), pac(6*1), prim(6*1*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:6*1) = 0.0_dp
      kac(1:6*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,6
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*6+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*6+ma)
              kad((md-1)*6+ma) =  kad((md-1)*6+ma)-tmp*p_bc
              kac((mc-1)*6+ma) = kac((mc-1)*6+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_6_1_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_6_1_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*2), kbc(1*1), kad(6*2), kac(6*1), pbd(1*2), &
      pbc(1*1), pad(6*2), pac(6*1), prim(6*1*1*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*2) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:6*2) = 0.0_dp
      kac(1:6*1) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,6
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*6+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*6+ma)
              kad((md-1)*6+ma) =  kad((md-1)*6+ma)-tmp*p_bc
              kac((mc-1)*6+ma) = kac((mc-1)*6+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_6_1_1_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_6_1_1_3(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*3), kbc(1*1), kad(6*3), kac(6*1), pbd(1*3), &
      pbc(1*1), pad(6*3), pac(6*1), prim(6*1*1*3), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*3) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:6*3) = 0.0_dp
      kac(1:6*1) = 0.0_dp
      p_index = 0
      DO md = 1,3
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,6
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*6+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*6+ma)
              kad((md-1)*6+ma) =  kad((md-1)*6+ma)-tmp*p_bc
              kac((mc-1)*6+ma) = kac((mc-1)*6+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_6_1_1_3
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_6_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*1), kad(6*md_max), kac(6*1), &
      pbd(1*md_max), pbc(1*1), pad(6*md_max), pac(6*1), prim(6*1*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:6*md_max) = 0.0_dp
      kac(1:6*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,6
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*6+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*6+ma)
              kad((md-1)*6+ma) =  kad((md-1)*6+ma)-tmp*p_bc
              kac((mc-1)*6+ma) = kac((mc-1)*6+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_6_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_6_1_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*2), kad(6*1), kac(6*2), pbd(1*1), &
      pbc(1*2), pad(6*1), pac(6*2), prim(6*1*2*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:6*1) = 0.0_dp
      kac(1:6*2) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,6
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*6+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*6+ma)
              kad((md-1)*6+ma) =  kad((md-1)*6+ma)-tmp*p_bc
              kac((mc-1)*6+ma) = kac((mc-1)*6+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_6_1_2_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_6_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*2), kad(6*md_max), kac(6*2), &
      pbd(1*md_max), pbc(1*2), pad(6*md_max), pac(6*2), prim(6*1*2*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:6*md_max) = 0.0_dp
      kac(1:6*2) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,6
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*6+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*6+ma)
              kad((md-1)*6+ma) =  kad((md-1)*6+ma)-tmp*p_bc
              kac((mc-1)*6+ma) = kac((mc-1)*6+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_6_1_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_6_1_3_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*3), kad(6*1), kac(6*3), pbd(1*1), &
      pbc(1*3), pad(6*1), pac(6*3), prim(6*1*3*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*3) = 0.0_dp
      kad(1:6*1) = 0.0_dp
      kac(1:6*3) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,3
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,6
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*6+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*6+ma)
              kad((md-1)*6+ma) =  kad((md-1)*6+ma)-tmp*p_bc
              kac((mc-1)*6+ma) = kac((mc-1)*6+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_6_1_3_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_6_1_3(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*3), kad(6*md_max), kac(6*3), &
      pbd(1*md_max), pbc(1*3), pad(6*md_max), pac(6*3), prim(6*1*3*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*3) = 0.0_dp
      kad(1:6*md_max) = 0.0_dp
      kac(1:6*3) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,3
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,6
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*6+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*6+ma)
              kad((md-1)*6+ma) =  kad((md-1)*6+ma)-tmp*p_bc
              kac((mc-1)*6+ma) = kac((mc-1)*6+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_6_1_3
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_6_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*mc_max), kad(6*md_max), &
      kac(6*mc_max), pbd(1*md_max), pbc(1*mc_max), pad(6*md_max), &
      pac(6*mc_max), prim(6*1*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*mc_max) = 0.0_dp
      kad(1:6*md_max) = 0.0_dp
      kac(1:6*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,6
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*6+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*6+ma)
              kad((md-1)*6+ma) =  kad((md-1)*6+ma)-tmp*p_bc
              kac((mc-1)*6+ma) = kac((mc-1)*6+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_6_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_6_2_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*1), kbc(2*1), kad(6*1), kac(6*1), pbd(2*1), &
      pbc(2*1), pad(6*1), pac(6*1), prim(6*2*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*1) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:6*1) = 0.0_dp
      kac(1:6*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,6
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*6+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*6+ma)
              kad((md-1)*6+ma) =  kad((md-1)*6+ma)-tmp*p_bc
              kac((mc-1)*6+ma) = kac((mc-1)*6+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_6_2_1_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_6_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*1), kad(6*md_max), kac(6*1), &
      pbd(2*md_max), pbc(2*1), pad(6*md_max), pac(6*1), prim(6*2*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:6*md_max) = 0.0_dp
      kac(1:6*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,6
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*6+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*6+ma)
              kad((md-1)*6+ma) =  kad((md-1)*6+ma)-tmp*p_bc
              kac((mc-1)*6+ma) = kac((mc-1)*6+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_6_2_1
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_6_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*mc_max), kad(6*md_max), &
      kac(6*mc_max), pbd(2*md_max), pbc(2*mc_max), pad(6*md_max), &
      pac(6*mc_max), prim(6*2*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*mc_max) = 0.0_dp
      kad(1:6*md_max) = 0.0_dp
      kac(1:6*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,6
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*6+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*6+ma)
              kad((md-1)*6+ma) =  kad((md-1)*6+ma)-tmp*p_bc
              kac((mc-1)*6+ma) = kac((mc-1)*6+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_6_2
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_6_3_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(3*1), kbc(3*1), kad(6*1), kac(6*1), pbd(3*1), &
      pbc(3*1), pad(6*1), pac(6*1), prim(6*3*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*1) = 0.0_dp
      kbc(1:3*1) = 0.0_dp
      kad(1:6*1) = 0.0_dp
      kac(1:6*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,6
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*6+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*6+ma)
              kad((md-1)*6+ma) =  kad((md-1)*6+ma)-tmp*p_bc
              kac((mc-1)*6+ma) = kac((mc-1)*6+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_6_3_1_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_6_3_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(3*md_max), kbc(3*1), kad(6*md_max), kac(6*1), &
      pbd(3*md_max), pbc(3*1), pad(6*md_max), pac(6*1), prim(6*3*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*md_max) = 0.0_dp
      kbc(1:3*1) = 0.0_dp
      kad(1:6*md_max) = 0.0_dp
      kac(1:6*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,6
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*6+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*6+ma)
              kad((md-1)*6+ma) =  kad((md-1)*6+ma)-tmp*p_bc
              kac((mc-1)*6+ma) = kac((mc-1)*6+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_6_3_1
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_6_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(3*md_max), kbc(3*mc_max), kad(6*md_max), &
      kac(6*mc_max), pbd(3*md_max), pbc(3*mc_max), pad(6*md_max), &
      pac(6*mc_max), prim(6*3*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*md_max) = 0.0_dp
      kbc(1:3*mc_max) = 0.0_dp
      kad(1:6*md_max) = 0.0_dp
      kac(1:6*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,6
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*6+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*6+ma)
              kad((md-1)*6+ma) =  kad((md-1)*6+ma)-tmp*p_bc
              kac((mc-1)*6+ma) = kac((mc-1)*6+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_6_3
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_6_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(4*md_max), kbc(4*mc_max), kad(6*md_max), &
      kac(6*mc_max), pbd(4*md_max), pbc(4*mc_max), pad(6*md_max), &
      pac(6*mc_max), prim(6*4*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*md_max) = 0.0_dp
      kbc(1:4*mc_max) = 0.0_dp
      kad(1:6*md_max) = 0.0_dp
      kac(1:6*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,6
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*6+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*6+ma)
              kad((md-1)*6+ma) =  kad((md-1)*6+ma)-tmp*p_bc
              kac((mc-1)*6+ma) = kac((mc-1)*6+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_6_4
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_6_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(5*md_max), kbc(5*mc_max), kad(6*md_max), &
      kac(6*mc_max), pbd(5*md_max), pbc(5*mc_max), pad(6*md_max), &
      pac(6*mc_max), prim(6*5*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:5*md_max) = 0.0_dp
      kbc(1:5*mc_max) = 0.0_dp
      kad(1:6*md_max) = 0.0_dp
      kac(1:6*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,5
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*5+mb)
            p_bc = pbc((mc-1)*5+mb)
            DO ma = 1,6
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*6+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*6+ma)
              kad((md-1)*6+ma) =  kad((md-1)*6+ma)-tmp*p_bc
              kac((mc-1)*6+ma) = kac((mc-1)*6+ma)-tmp*p_bd
            END DO
            kbd((md-1)*5+mb) = kbd((md-1)*5+mb) - ks_bd
            kbc((mc-1)*5+mb) = kbc((mc-1)*5+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_6_5
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_6_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(6*md_max), kbc(6*mc_max), kad(6*md_max), &
      kac(6*mc_max), pbd(6*md_max), pbc(6*mc_max), pad(6*md_max), &
      pac(6*mc_max), prim(6*6*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:6*md_max) = 0.0_dp
      kbc(1:6*mc_max) = 0.0_dp
      kad(1:6*md_max) = 0.0_dp
      kac(1:6*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,6
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*6+mb)
            p_bc = pbc((mc-1)*6+mb)
            DO ma = 1,6
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*6+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*6+ma)
              kad((md-1)*6+ma) =  kad((md-1)*6+ma)-tmp*p_bc
              kac((mc-1)*6+ma) = kac((mc-1)*6+ma)-tmp*p_bd
            END DO
            kbd((md-1)*6+mb) = kbd((md-1)*6+mb) - ks_bd
            kbc((mc-1)*6+mb) = kbc((mc-1)*6+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_6_6
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_6_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(7*md_max), kbc(7*mc_max), kad(6*md_max), &
      kac(6*mc_max), pbd(7*md_max), pbc(7*mc_max), pad(6*md_max), &
      pac(6*mc_max), prim(6*7*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:7*md_max) = 0.0_dp
      kbc(1:7*mc_max) = 0.0_dp
      kad(1:6*md_max) = 0.0_dp
      kac(1:6*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,7
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*7+mb)
            p_bc = pbc((mc-1)*7+mb)
            DO ma = 1,6
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*6+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*6+ma)
              kad((md-1)*6+ma) =  kad((md-1)*6+ma)-tmp*p_bc
              kac((mc-1)*6+ma) = kac((mc-1)*6+ma)-tmp*p_bd
            END DO
            kbd((md-1)*7+mb) = kbd((md-1)*7+mb) - ks_bd
            kbc((mc-1)*7+mb) = kbc((mc-1)*7+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_6_7
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_6_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(9*md_max), kbc(9*mc_max), kad(6*md_max), &
      kac(6*mc_max), pbd(9*md_max), pbc(9*mc_max), pad(6*md_max), &
      pac(6*mc_max), prim(6*9*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:9*md_max) = 0.0_dp
      kbc(1:9*mc_max) = 0.0_dp
      kad(1:6*md_max) = 0.0_dp
      kac(1:6*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,9
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*9+mb)
            p_bc = pbc((mc-1)*9+mb)
            DO ma = 1,6
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*6+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*6+ma)
              kad((md-1)*6+ma) =  kad((md-1)*6+ma)-tmp*p_bc
              kac((mc-1)*6+ma) = kac((mc-1)*6+ma)-tmp*p_bd
            END DO
            kbd((md-1)*9+mb) = kbd((md-1)*9+mb) - ks_bd
            kbc((mc-1)*9+mb) = kbc((mc-1)*9+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_6_9
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_6_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(10*md_max), kbc(10*mc_max), kad(6*md_max), &
      kac(6*mc_max), pbd(10*md_max), pbc(10*mc_max), pad(6*md_max), &
      pac(6*mc_max), prim(6*10*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:10*md_max) = 0.0_dp
      kbc(1:10*mc_max) = 0.0_dp
      kad(1:6*md_max) = 0.0_dp
      kac(1:6*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,10
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*10+mb)
            p_bc = pbc((mc-1)*10+mb)
            DO ma = 1,6
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*6+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*6+ma)
              kad((md-1)*6+ma) =  kad((md-1)*6+ma)-tmp*p_bc
              kac((mc-1)*6+ma) = kac((mc-1)*6+ma)-tmp*p_bd
            END DO
            kbd((md-1)*10+mb) = kbd((md-1)*10+mb) - ks_bd
            kbc((mc-1)*10+mb) = kbc((mc-1)*10+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_6_10
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_6_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(11*md_max), kbc(11*mc_max), kad(6*md_max), &
      kac(6*mc_max), pbd(11*md_max), pbc(11*mc_max), pad(6*md_max), &
      pac(6*mc_max), prim(6*11*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:11*md_max) = 0.0_dp
      kbc(1:11*mc_max) = 0.0_dp
      kad(1:6*md_max) = 0.0_dp
      kac(1:6*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,11
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*11+mb)
            p_bc = pbc((mc-1)*11+mb)
            DO ma = 1,6
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*6+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*6+ma)
              kad((md-1)*6+ma) =  kad((md-1)*6+ma)-tmp*p_bc
              kac((mc-1)*6+ma) = kac((mc-1)*6+ma)-tmp*p_bd
            END DO
            kbd((md-1)*11+mb) = kbd((md-1)*11+mb) - ks_bd
            kbc((mc-1)*11+mb) = kbc((mc-1)*11+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_6_11
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_6_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(15*md_max), kbc(15*mc_max), kad(6*md_max), &
      kac(6*mc_max), pbd(15*md_max), pbc(15*mc_max), pad(6*md_max), &
      pac(6*mc_max), prim(6*15*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:15*md_max) = 0.0_dp
      kbc(1:15*mc_max) = 0.0_dp
      kad(1:6*md_max) = 0.0_dp
      kac(1:6*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,15
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*15+mb)
            p_bc = pbc((mc-1)*15+mb)
            DO ma = 1,6
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*6+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*6+ma)
              kad((md-1)*6+ma) =  kad((md-1)*6+ma)-tmp*p_bc
              kac((mc-1)*6+ma) = kac((mc-1)*6+ma)-tmp*p_bd
            END DO
            kbd((md-1)*15+mb) = kbd((md-1)*15+mb) - ks_bd
            kbc((mc-1)*15+mb) = kbc((mc-1)*15+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_6_15
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_7_1_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*1), kad(7*1), kac(7*1), pbd(1*1), &
      pbc(1*1), pad(7*1), pac(7*1), prim(7*1*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:7*1) = 0.0_dp
      kac(1:7*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,7
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*7+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*7+ma)
              kad((md-1)*7+ma) =  kad((md-1)*7+ma)-tmp*p_bc
              kac((mc-1)*7+ma) = kac((mc-1)*7+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_7_1_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_7_1_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*2), kbc(1*1), kad(7*2), kac(7*1), pbd(1*2), &
      pbc(1*1), pad(7*2), pac(7*1), prim(7*1*1*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*2) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:7*2) = 0.0_dp
      kac(1:7*1) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,7
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*7+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*7+ma)
              kad((md-1)*7+ma) =  kad((md-1)*7+ma)-tmp*p_bc
              kac((mc-1)*7+ma) = kac((mc-1)*7+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_7_1_1_2
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_7_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*1), kad(7*md_max), kac(7*1), &
      pbd(1*md_max), pbc(1*1), pad(7*md_max), pac(7*1), prim(7*1*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:7*md_max) = 0.0_dp
      kac(1:7*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,7
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*7+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*7+ma)
              kad((md-1)*7+ma) =  kad((md-1)*7+ma)-tmp*p_bc
              kac((mc-1)*7+ma) = kac((mc-1)*7+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_7_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_7_1_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*2), kad(7*1), kac(7*2), pbd(1*1), &
      pbc(1*2), pad(7*1), pac(7*2), prim(7*1*2*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:7*1) = 0.0_dp
      kac(1:7*2) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,7
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*7+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*7+ma)
              kad((md-1)*7+ma) =  kad((md-1)*7+ma)-tmp*p_bc
              kac((mc-1)*7+ma) = kac((mc-1)*7+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_7_1_2_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_7_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*2), kad(7*md_max), kac(7*2), &
      pbd(1*md_max), pbc(1*2), pad(7*md_max), pac(7*2), prim(7*1*2*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:7*md_max) = 0.0_dp
      kac(1:7*2) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,7
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*7+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*7+ma)
              kad((md-1)*7+ma) =  kad((md-1)*7+ma)-tmp*p_bc
              kac((mc-1)*7+ma) = kac((mc-1)*7+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_7_1_2
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_7_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*mc_max), kad(7*md_max), &
      kac(7*mc_max), pbd(1*md_max), pbc(1*mc_max), pad(7*md_max), &
      pac(7*mc_max), prim(7*1*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*mc_max) = 0.0_dp
      kad(1:7*md_max) = 0.0_dp
      kac(1:7*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,7
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*7+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*7+ma)
              kad((md-1)*7+ma) =  kad((md-1)*7+ma)-tmp*p_bc
              kac((mc-1)*7+ma) = kac((mc-1)*7+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_7_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_7_2_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*1), kbc(2*1), kad(7*1), kac(7*1), pbd(2*1), &
      pbc(2*1), pad(7*1), pac(7*1), prim(7*2*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*1) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:7*1) = 0.0_dp
      kac(1:7*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,7
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*7+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*7+ma)
              kad((md-1)*7+ma) =  kad((md-1)*7+ma)-tmp*p_bc
              kac((mc-1)*7+ma) = kac((mc-1)*7+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_7_2_1_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_7_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*1), kad(7*md_max), kac(7*1), &
      pbd(2*md_max), pbc(2*1), pad(7*md_max), pac(7*1), prim(7*2*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:7*md_max) = 0.0_dp
      kac(1:7*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,7
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*7+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*7+ma)
              kad((md-1)*7+ma) =  kad((md-1)*7+ma)-tmp*p_bc
              kac((mc-1)*7+ma) = kac((mc-1)*7+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_7_2_1
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_7_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*mc_max), kad(7*md_max), &
      kac(7*mc_max), pbd(2*md_max), pbc(2*mc_max), pad(7*md_max), &
      pac(7*mc_max), prim(7*2*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*mc_max) = 0.0_dp
      kad(1:7*md_max) = 0.0_dp
      kac(1:7*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,7
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*7+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*7+ma)
              kad((md-1)*7+ma) =  kad((md-1)*7+ma)-tmp*p_bc
              kac((mc-1)*7+ma) = kac((mc-1)*7+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_7_2
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_7_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(3*md_max), kbc(3*mc_max), kad(7*md_max), &
      kac(7*mc_max), pbd(3*md_max), pbc(3*mc_max), pad(7*md_max), &
      pac(7*mc_max), prim(7*3*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*md_max) = 0.0_dp
      kbc(1:3*mc_max) = 0.0_dp
      kad(1:7*md_max) = 0.0_dp
      kac(1:7*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,7
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*7+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*7+ma)
              kad((md-1)*7+ma) =  kad((md-1)*7+ma)-tmp*p_bc
              kac((mc-1)*7+ma) = kac((mc-1)*7+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_7_3
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_7_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(4*md_max), kbc(4*mc_max), kad(7*md_max), &
      kac(7*mc_max), pbd(4*md_max), pbc(4*mc_max), pad(7*md_max), &
      pac(7*mc_max), prim(7*4*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*md_max) = 0.0_dp
      kbc(1:4*mc_max) = 0.0_dp
      kad(1:7*md_max) = 0.0_dp
      kac(1:7*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,7
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*7+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*7+ma)
              kad((md-1)*7+ma) =  kad((md-1)*7+ma)-tmp*p_bc
              kac((mc-1)*7+ma) = kac((mc-1)*7+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_7_4
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_7_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(5*md_max), kbc(5*mc_max), kad(7*md_max), &
      kac(7*mc_max), pbd(5*md_max), pbc(5*mc_max), pad(7*md_max), &
      pac(7*mc_max), prim(7*5*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:5*md_max) = 0.0_dp
      kbc(1:5*mc_max) = 0.0_dp
      kad(1:7*md_max) = 0.0_dp
      kac(1:7*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,5
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*5+mb)
            p_bc = pbc((mc-1)*5+mb)
            DO ma = 1,7
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*7+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*7+ma)
              kad((md-1)*7+ma) =  kad((md-1)*7+ma)-tmp*p_bc
              kac((mc-1)*7+ma) = kac((mc-1)*7+ma)-tmp*p_bd
            END DO
            kbd((md-1)*5+mb) = kbd((md-1)*5+mb) - ks_bd
            kbc((mc-1)*5+mb) = kbc((mc-1)*5+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_7_5
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_7_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(6*md_max), kbc(6*mc_max), kad(7*md_max), &
      kac(7*mc_max), pbd(6*md_max), pbc(6*mc_max), pad(7*md_max), &
      pac(7*mc_max), prim(7*6*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:6*md_max) = 0.0_dp
      kbc(1:6*mc_max) = 0.0_dp
      kad(1:7*md_max) = 0.0_dp
      kac(1:7*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,6
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*6+mb)
            p_bc = pbc((mc-1)*6+mb)
            DO ma = 1,7
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*7+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*7+ma)
              kad((md-1)*7+ma) =  kad((md-1)*7+ma)-tmp*p_bc
              kac((mc-1)*7+ma) = kac((mc-1)*7+ma)-tmp*p_bd
            END DO
            kbd((md-1)*6+mb) = kbd((md-1)*6+mb) - ks_bd
            kbc((mc-1)*6+mb) = kbc((mc-1)*6+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_7_6
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_7_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(7*md_max), kbc(7*mc_max), kad(7*md_max), &
      kac(7*mc_max), pbd(7*md_max), pbc(7*mc_max), pad(7*md_max), &
      pac(7*mc_max), prim(7*7*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:7*md_max) = 0.0_dp
      kbc(1:7*mc_max) = 0.0_dp
      kad(1:7*md_max) = 0.0_dp
      kac(1:7*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,7
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*7+mb)
            p_bc = pbc((mc-1)*7+mb)
            DO ma = 1,7
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*7+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*7+ma)
              kad((md-1)*7+ma) =  kad((md-1)*7+ma)-tmp*p_bc
              kac((mc-1)*7+ma) = kac((mc-1)*7+ma)-tmp*p_bd
            END DO
            kbd((md-1)*7+mb) = kbd((md-1)*7+mb) - ks_bd
            kbc((mc-1)*7+mb) = kbc((mc-1)*7+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_7_7
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_7_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(9*md_max), kbc(9*mc_max), kad(7*md_max), &
      kac(7*mc_max), pbd(9*md_max), pbc(9*mc_max), pad(7*md_max), &
      pac(7*mc_max), prim(7*9*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:9*md_max) = 0.0_dp
      kbc(1:9*mc_max) = 0.0_dp
      kad(1:7*md_max) = 0.0_dp
      kac(1:7*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,9
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*9+mb)
            p_bc = pbc((mc-1)*9+mb)
            DO ma = 1,7
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*7+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*7+ma)
              kad((md-1)*7+ma) =  kad((md-1)*7+ma)-tmp*p_bc
              kac((mc-1)*7+ma) = kac((mc-1)*7+ma)-tmp*p_bd
            END DO
            kbd((md-1)*9+mb) = kbd((md-1)*9+mb) - ks_bd
            kbc((mc-1)*9+mb) = kbc((mc-1)*9+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_7_9
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_7_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(10*md_max), kbc(10*mc_max), kad(7*md_max), &
      kac(7*mc_max), pbd(10*md_max), pbc(10*mc_max), pad(7*md_max), &
      pac(7*mc_max), prim(7*10*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:10*md_max) = 0.0_dp
      kbc(1:10*mc_max) = 0.0_dp
      kad(1:7*md_max) = 0.0_dp
      kac(1:7*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,10
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*10+mb)
            p_bc = pbc((mc-1)*10+mb)
            DO ma = 1,7
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*7+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*7+ma)
              kad((md-1)*7+ma) =  kad((md-1)*7+ma)-tmp*p_bc
              kac((mc-1)*7+ma) = kac((mc-1)*7+ma)-tmp*p_bd
            END DO
            kbd((md-1)*10+mb) = kbd((md-1)*10+mb) - ks_bd
            kbc((mc-1)*10+mb) = kbc((mc-1)*10+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_7_10
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_7_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(11*md_max), kbc(11*mc_max), kad(7*md_max), &
      kac(7*mc_max), pbd(11*md_max), pbc(11*mc_max), pad(7*md_max), &
      pac(7*mc_max), prim(7*11*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:11*md_max) = 0.0_dp
      kbc(1:11*mc_max) = 0.0_dp
      kad(1:7*md_max) = 0.0_dp
      kac(1:7*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,11
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*11+mb)
            p_bc = pbc((mc-1)*11+mb)
            DO ma = 1,7
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*7+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*7+ma)
              kad((md-1)*7+ma) =  kad((md-1)*7+ma)-tmp*p_bc
              kac((mc-1)*7+ma) = kac((mc-1)*7+ma)-tmp*p_bd
            END DO
            kbd((md-1)*11+mb) = kbd((md-1)*11+mb) - ks_bd
            kbc((mc-1)*11+mb) = kbc((mc-1)*11+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_7_11
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_7_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(15*md_max), kbc(15*mc_max), kad(7*md_max), &
      kac(7*mc_max), pbd(15*md_max), pbc(15*mc_max), pad(7*md_max), &
      pac(7*mc_max), prim(7*15*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:15*md_max) = 0.0_dp
      kbc(1:15*mc_max) = 0.0_dp
      kad(1:7*md_max) = 0.0_dp
      kac(1:7*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,15
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*15+mb)
            p_bc = pbc((mc-1)*15+mb)
            DO ma = 1,7
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*7+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*7+ma)
              kad((md-1)*7+ma) =  kad((md-1)*7+ma)-tmp*p_bc
              kac((mc-1)*7+ma) = kac((mc-1)*7+ma)-tmp*p_bd
            END DO
            kbd((md-1)*15+mb) = kbd((md-1)*15+mb) - ks_bd
            kbc((mc-1)*15+mb) = kbc((mc-1)*15+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_7_15
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_9_1_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*1), kad(9*1), kac(9*1), pbd(1*1), &
      pbc(1*1), pad(9*1), pac(9*1), prim(9*1*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:9*1) = 0.0_dp
      kac(1:9*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,9
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*9+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*9+ma)
              kad((md-1)*9+ma) =  kad((md-1)*9+ma)-tmp*p_bc
              kac((mc-1)*9+ma) = kac((mc-1)*9+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_9_1_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_9_1_1_2(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*2), kbc(1*1), kad(9*2), kac(9*1), pbd(1*2), &
      pbc(1*1), pad(9*2), pac(9*1), prim(9*1*1*2), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*2) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:9*2) = 0.0_dp
      kac(1:9*1) = 0.0_dp
      p_index = 0
      DO md = 1,2
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,9
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*9+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*9+ma)
              kad((md-1)*9+ma) =  kad((md-1)*9+ma)-tmp*p_bc
              kac((mc-1)*9+ma) = kac((mc-1)*9+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_9_1_1_2
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_9_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*1), kad(9*md_max), kac(9*1), &
      pbd(1*md_max), pbc(1*1), pad(9*md_max), pac(9*1), prim(9*1*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:9*md_max) = 0.0_dp
      kac(1:9*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,9
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*9+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*9+ma)
              kad((md-1)*9+ma) =  kad((md-1)*9+ma)-tmp*p_bc
              kac((mc-1)*9+ma) = kac((mc-1)*9+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_9_1_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_9_1_2_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*2), kad(9*1), kac(9*2), pbd(1*1), &
      pbc(1*2), pad(9*1), pac(9*2), prim(9*1*2*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:9*1) = 0.0_dp
      kac(1:9*2) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,9
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*9+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*9+ma)
              kad((md-1)*9+ma) =  kad((md-1)*9+ma)-tmp*p_bc
              kac((mc-1)*9+ma) = kac((mc-1)*9+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_9_1_2_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_9_1_2(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*2), kad(9*md_max), kac(9*2), &
      pbd(1*md_max), pbc(1*2), pad(9*md_max), pac(9*2), prim(9*1*2*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*2) = 0.0_dp
      kad(1:9*md_max) = 0.0_dp
      kac(1:9*2) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,2
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,9
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*9+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*9+ma)
              kad((md-1)*9+ma) =  kad((md-1)*9+ma)-tmp*p_bc
              kac((mc-1)*9+ma) = kac((mc-1)*9+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_9_1_2
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_9_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*mc_max), kad(9*md_max), &
      kac(9*mc_max), pbd(1*md_max), pbc(1*mc_max), pad(9*md_max), &
      pac(9*mc_max), prim(9*1*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*mc_max) = 0.0_dp
      kad(1:9*md_max) = 0.0_dp
      kac(1:9*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,9
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*9+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*9+ma)
              kad((md-1)*9+ma) =  kad((md-1)*9+ma)-tmp*p_bc
              kac((mc-1)*9+ma) = kac((mc-1)*9+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_9_1
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_9_2_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(2*1), kbc(2*1), kad(9*1), kac(9*1), pbd(2*1), &
      pbc(2*1), pad(9*1), pac(9*1), prim(9*2*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*1) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:9*1) = 0.0_dp
      kac(1:9*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,9
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*9+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*9+ma)
              kad((md-1)*9+ma) =  kad((md-1)*9+ma)-tmp*p_bc
              kac((mc-1)*9+ma) = kac((mc-1)*9+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_9_2_1_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_9_2_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*1), kad(9*md_max), kac(9*1), &
      pbd(2*md_max), pbc(2*1), pad(9*md_max), pac(9*1), prim(9*2*1*md_max), &
      scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*1) = 0.0_dp
      kad(1:9*md_max) = 0.0_dp
      kac(1:9*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,9
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*9+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*9+ma)
              kad((md-1)*9+ma) =  kad((md-1)*9+ma)-tmp*p_bc
              kac((mc-1)*9+ma) = kac((mc-1)*9+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_9_2_1
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_9_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*mc_max), kad(9*md_max), &
      kac(9*mc_max), pbd(2*md_max), pbc(2*mc_max), pad(9*md_max), &
      pac(9*mc_max), prim(9*2*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*mc_max) = 0.0_dp
      kad(1:9*md_max) = 0.0_dp
      kac(1:9*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,9
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*9+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*9+ma)
              kad((md-1)*9+ma) =  kad((md-1)*9+ma)-tmp*p_bc
              kac((mc-1)*9+ma) = kac((mc-1)*9+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_9_2
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_9_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(3*md_max), kbc(3*mc_max), kad(9*md_max), &
      kac(9*mc_max), pbd(3*md_max), pbc(3*mc_max), pad(9*md_max), &
      pac(9*mc_max), prim(9*3*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*md_max) = 0.0_dp
      kbc(1:3*mc_max) = 0.0_dp
      kad(1:9*md_max) = 0.0_dp
      kac(1:9*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,9
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*9+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*9+ma)
              kad((md-1)*9+ma) =  kad((md-1)*9+ma)-tmp*p_bc
              kac((mc-1)*9+ma) = kac((mc-1)*9+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_9_3
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_9_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(4*md_max), kbc(4*mc_max), kad(9*md_max), &
      kac(9*mc_max), pbd(4*md_max), pbc(4*mc_max), pad(9*md_max), &
      pac(9*mc_max), prim(9*4*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*md_max) = 0.0_dp
      kbc(1:4*mc_max) = 0.0_dp
      kad(1:9*md_max) = 0.0_dp
      kac(1:9*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,9
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*9+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*9+ma)
              kad((md-1)*9+ma) =  kad((md-1)*9+ma)-tmp*p_bc
              kac((mc-1)*9+ma) = kac((mc-1)*9+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_9_4
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_9_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(5*md_max), kbc(5*mc_max), kad(9*md_max), &
      kac(9*mc_max), pbd(5*md_max), pbc(5*mc_max), pad(9*md_max), &
      pac(9*mc_max), prim(9*5*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:5*md_max) = 0.0_dp
      kbc(1:5*mc_max) = 0.0_dp
      kad(1:9*md_max) = 0.0_dp
      kac(1:9*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,5
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*5+mb)
            p_bc = pbc((mc-1)*5+mb)
            DO ma = 1,9
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*9+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*9+ma)
              kad((md-1)*9+ma) =  kad((md-1)*9+ma)-tmp*p_bc
              kac((mc-1)*9+ma) = kac((mc-1)*9+ma)-tmp*p_bd
            END DO
            kbd((md-1)*5+mb) = kbd((md-1)*5+mb) - ks_bd
            kbc((mc-1)*5+mb) = kbc((mc-1)*5+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_9_5
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_9_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(6*md_max), kbc(6*mc_max), kad(9*md_max), &
      kac(9*mc_max), pbd(6*md_max), pbc(6*mc_max), pad(9*md_max), &
      pac(9*mc_max), prim(9*6*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:6*md_max) = 0.0_dp
      kbc(1:6*mc_max) = 0.0_dp
      kad(1:9*md_max) = 0.0_dp
      kac(1:9*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,6
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*6+mb)
            p_bc = pbc((mc-1)*6+mb)
            DO ma = 1,9
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*9+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*9+ma)
              kad((md-1)*9+ma) =  kad((md-1)*9+ma)-tmp*p_bc
              kac((mc-1)*9+ma) = kac((mc-1)*9+ma)-tmp*p_bd
            END DO
            kbd((md-1)*6+mb) = kbd((md-1)*6+mb) - ks_bd
            kbc((mc-1)*6+mb) = kbc((mc-1)*6+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_9_6
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_9_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(7*md_max), kbc(7*mc_max), kad(9*md_max), &
      kac(9*mc_max), pbd(7*md_max), pbc(7*mc_max), pad(9*md_max), &
      pac(9*mc_max), prim(9*7*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:7*md_max) = 0.0_dp
      kbc(1:7*mc_max) = 0.0_dp
      kad(1:9*md_max) = 0.0_dp
      kac(1:9*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,7
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*7+mb)
            p_bc = pbc((mc-1)*7+mb)
            DO ma = 1,9
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*9+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*9+ma)
              kad((md-1)*9+ma) =  kad((md-1)*9+ma)-tmp*p_bc
              kac((mc-1)*9+ma) = kac((mc-1)*9+ma)-tmp*p_bd
            END DO
            kbd((md-1)*7+mb) = kbd((md-1)*7+mb) - ks_bd
            kbc((mc-1)*7+mb) = kbc((mc-1)*7+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_9_7
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_9_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(9*md_max), kbc(9*mc_max), kad(9*md_max), &
      kac(9*mc_max), pbd(9*md_max), pbc(9*mc_max), pad(9*md_max), &
      pac(9*mc_max), prim(9*9*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:9*md_max) = 0.0_dp
      kbc(1:9*mc_max) = 0.0_dp
      kad(1:9*md_max) = 0.0_dp
      kac(1:9*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,9
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*9+mb)
            p_bc = pbc((mc-1)*9+mb)
            DO ma = 1,9
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*9+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*9+ma)
              kad((md-1)*9+ma) =  kad((md-1)*9+ma)-tmp*p_bc
              kac((mc-1)*9+ma) = kac((mc-1)*9+ma)-tmp*p_bd
            END DO
            kbd((md-1)*9+mb) = kbd((md-1)*9+mb) - ks_bd
            kbc((mc-1)*9+mb) = kbc((mc-1)*9+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_9_9
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_9_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(10*md_max), kbc(10*mc_max), kad(9*md_max), &
      kac(9*mc_max), pbd(10*md_max), pbc(10*mc_max), pad(9*md_max), &
      pac(9*mc_max), prim(9*10*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:10*md_max) = 0.0_dp
      kbc(1:10*mc_max) = 0.0_dp
      kad(1:9*md_max) = 0.0_dp
      kac(1:9*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,10
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*10+mb)
            p_bc = pbc((mc-1)*10+mb)
            DO ma = 1,9
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*9+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*9+ma)
              kad((md-1)*9+ma) =  kad((md-1)*9+ma)-tmp*p_bc
              kac((mc-1)*9+ma) = kac((mc-1)*9+ma)-tmp*p_bd
            END DO
            kbd((md-1)*10+mb) = kbd((md-1)*10+mb) - ks_bd
            kbc((mc-1)*10+mb) = kbc((mc-1)*10+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_9_10
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_9_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(11*md_max), kbc(11*mc_max), kad(9*md_max), &
      kac(9*mc_max), pbd(11*md_max), pbc(11*mc_max), pad(9*md_max), &
      pac(9*mc_max), prim(9*11*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:11*md_max) = 0.0_dp
      kbc(1:11*mc_max) = 0.0_dp
      kad(1:9*md_max) = 0.0_dp
      kac(1:9*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,11
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*11+mb)
            p_bc = pbc((mc-1)*11+mb)
            DO ma = 1,9
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*9+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*9+ma)
              kad((md-1)*9+ma) =  kad((md-1)*9+ma)-tmp*p_bc
              kac((mc-1)*9+ma) = kac((mc-1)*9+ma)-tmp*p_bd
            END DO
            kbd((md-1)*11+mb) = kbd((md-1)*11+mb) - ks_bd
            kbc((mc-1)*11+mb) = kbc((mc-1)*11+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_9_11
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_9_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(15*md_max), kbc(15*mc_max), kad(9*md_max), &
      kac(9*mc_max), pbd(15*md_max), pbc(15*mc_max), pad(9*md_max), &
      pac(9*mc_max), prim(9*15*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:15*md_max) = 0.0_dp
      kbc(1:15*mc_max) = 0.0_dp
      kad(1:9*md_max) = 0.0_dp
      kac(1:9*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,15
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*15+mb)
            p_bc = pbc((mc-1)*15+mb)
            DO ma = 1,9
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*9+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*9+ma)
              kad((md-1)*9+ma) =  kad((md-1)*9+ma)-tmp*p_bc
              kac((mc-1)*9+ma) = kac((mc-1)*9+ma)-tmp*p_bd
            END DO
            kbd((md-1)*15+mb) = kbd((md-1)*15+mb) - ks_bd
            kbc((mc-1)*15+mb) = kbc((mc-1)*15+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_9_15
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_10_1_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*1), kad(10*1), kac(10*1), pbd(1*1), &
      pbc(1*1), pad(10*1), pac(10*1), prim(10*1*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:10*1) = 0.0_dp
      kac(1:10*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,10
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*10+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*10+ma)
              kad((md-1)*10+ma) =  kad((md-1)*10+ma)-tmp*p_bc
              kac((mc-1)*10+ma) = kac((mc-1)*10+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_10_1_1_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_10_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*1), kad(10*md_max), kac(10*1), &
      pbd(1*md_max), pbc(1*1), pad(10*md_max), pac(10*1), &
      prim(10*1*1*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:10*md_max) = 0.0_dp
      kac(1:10*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,10
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*10+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*10+ma)
              kad((md-1)*10+ma) =  kad((md-1)*10+ma)-tmp*p_bc
              kac((mc-1)*10+ma) = kac((mc-1)*10+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_10_1_1
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_10_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*mc_max), kad(10*md_max), &
      kac(10*mc_max), pbd(1*md_max), pbc(1*mc_max), pad(10*md_max), &
      pac(10*mc_max), prim(10*1*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*mc_max) = 0.0_dp
      kad(1:10*md_max) = 0.0_dp
      kac(1:10*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,10
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*10+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*10+ma)
              kad((md-1)*10+ma) =  kad((md-1)*10+ma)-tmp*p_bc
              kac((mc-1)*10+ma) = kac((mc-1)*10+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_10_1
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_10_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*mc_max), kad(10*md_max), &
      kac(10*mc_max), pbd(2*md_max), pbc(2*mc_max), pad(10*md_max), &
      pac(10*mc_max), prim(10*2*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*mc_max) = 0.0_dp
      kad(1:10*md_max) = 0.0_dp
      kac(1:10*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,10
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*10+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*10+ma)
              kad((md-1)*10+ma) =  kad((md-1)*10+ma)-tmp*p_bc
              kac((mc-1)*10+ma) = kac((mc-1)*10+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_10_2
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_10_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(3*md_max), kbc(3*mc_max), kad(10*md_max), &
      kac(10*mc_max), pbd(3*md_max), pbc(3*mc_max), pad(10*md_max), &
      pac(10*mc_max), prim(10*3*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*md_max) = 0.0_dp
      kbc(1:3*mc_max) = 0.0_dp
      kad(1:10*md_max) = 0.0_dp
      kac(1:10*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,10
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*10+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*10+ma)
              kad((md-1)*10+ma) =  kad((md-1)*10+ma)-tmp*p_bc
              kac((mc-1)*10+ma) = kac((mc-1)*10+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_10_3
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_10_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(4*md_max), kbc(4*mc_max), kad(10*md_max), &
      kac(10*mc_max), pbd(4*md_max), pbc(4*mc_max), pad(10*md_max), &
      pac(10*mc_max), prim(10*4*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*md_max) = 0.0_dp
      kbc(1:4*mc_max) = 0.0_dp
      kad(1:10*md_max) = 0.0_dp
      kac(1:10*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,10
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*10+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*10+ma)
              kad((md-1)*10+ma) =  kad((md-1)*10+ma)-tmp*p_bc
              kac((mc-1)*10+ma) = kac((mc-1)*10+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_10_4
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_10_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(5*md_max), kbc(5*mc_max), kad(10*md_max), &
      kac(10*mc_max), pbd(5*md_max), pbc(5*mc_max), pad(10*md_max), &
      pac(10*mc_max), prim(10*5*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:5*md_max) = 0.0_dp
      kbc(1:5*mc_max) = 0.0_dp
      kad(1:10*md_max) = 0.0_dp
      kac(1:10*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,5
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*5+mb)
            p_bc = pbc((mc-1)*5+mb)
            DO ma = 1,10
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*10+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*10+ma)
              kad((md-1)*10+ma) =  kad((md-1)*10+ma)-tmp*p_bc
              kac((mc-1)*10+ma) = kac((mc-1)*10+ma)-tmp*p_bd
            END DO
            kbd((md-1)*5+mb) = kbd((md-1)*5+mb) - ks_bd
            kbc((mc-1)*5+mb) = kbc((mc-1)*5+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_10_5
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_10_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(6*md_max), kbc(6*mc_max), kad(10*md_max), &
      kac(10*mc_max), pbd(6*md_max), pbc(6*mc_max), pad(10*md_max), &
      pac(10*mc_max), prim(10*6*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:6*md_max) = 0.0_dp
      kbc(1:6*mc_max) = 0.0_dp
      kad(1:10*md_max) = 0.0_dp
      kac(1:10*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,6
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*6+mb)
            p_bc = pbc((mc-1)*6+mb)
            DO ma = 1,10
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*10+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*10+ma)
              kad((md-1)*10+ma) =  kad((md-1)*10+ma)-tmp*p_bc
              kac((mc-1)*10+ma) = kac((mc-1)*10+ma)-tmp*p_bd
            END DO
            kbd((md-1)*6+mb) = kbd((md-1)*6+mb) - ks_bd
            kbc((mc-1)*6+mb) = kbc((mc-1)*6+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_10_6
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_10_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(7*md_max), kbc(7*mc_max), kad(10*md_max), &
      kac(10*mc_max), pbd(7*md_max), pbc(7*mc_max), pad(10*md_max), &
      pac(10*mc_max), prim(10*7*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:7*md_max) = 0.0_dp
      kbc(1:7*mc_max) = 0.0_dp
      kad(1:10*md_max) = 0.0_dp
      kac(1:10*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,7
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*7+mb)
            p_bc = pbc((mc-1)*7+mb)
            DO ma = 1,10
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*10+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*10+ma)
              kad((md-1)*10+ma) =  kad((md-1)*10+ma)-tmp*p_bc
              kac((mc-1)*10+ma) = kac((mc-1)*10+ma)-tmp*p_bd
            END DO
            kbd((md-1)*7+mb) = kbd((md-1)*7+mb) - ks_bd
            kbc((mc-1)*7+mb) = kbc((mc-1)*7+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_10_7
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_10_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(9*md_max), kbc(9*mc_max), kad(10*md_max), &
      kac(10*mc_max), pbd(9*md_max), pbc(9*mc_max), pad(10*md_max), &
      pac(10*mc_max), prim(10*9*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:9*md_max) = 0.0_dp
      kbc(1:9*mc_max) = 0.0_dp
      kad(1:10*md_max) = 0.0_dp
      kac(1:10*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,9
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*9+mb)
            p_bc = pbc((mc-1)*9+mb)
            DO ma = 1,10
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*10+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*10+ma)
              kad((md-1)*10+ma) =  kad((md-1)*10+ma)-tmp*p_bc
              kac((mc-1)*10+ma) = kac((mc-1)*10+ma)-tmp*p_bd
            END DO
            kbd((md-1)*9+mb) = kbd((md-1)*9+mb) - ks_bd
            kbc((mc-1)*9+mb) = kbc((mc-1)*9+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_10_9
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_10_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(10*md_max), kbc(10*mc_max), kad(10*md_max), &
      kac(10*mc_max), pbd(10*md_max), pbc(10*mc_max), pad(10*md_max), &
      pac(10*mc_max), prim(10*10*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:10*md_max) = 0.0_dp
      kbc(1:10*mc_max) = 0.0_dp
      kad(1:10*md_max) = 0.0_dp
      kac(1:10*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,10
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*10+mb)
            p_bc = pbc((mc-1)*10+mb)
            DO ma = 1,10
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*10+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*10+ma)
              kad((md-1)*10+ma) =  kad((md-1)*10+ma)-tmp*p_bc
              kac((mc-1)*10+ma) = kac((mc-1)*10+ma)-tmp*p_bd
            END DO
            kbd((md-1)*10+mb) = kbd((md-1)*10+mb) - ks_bd
            kbc((mc-1)*10+mb) = kbc((mc-1)*10+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_10_10
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_10_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(11*md_max), kbc(11*mc_max), kad(10*md_max), &
      kac(10*mc_max), pbd(11*md_max), pbc(11*mc_max), pad(10*md_max), &
      pac(10*mc_max), prim(10*11*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:11*md_max) = 0.0_dp
      kbc(1:11*mc_max) = 0.0_dp
      kad(1:10*md_max) = 0.0_dp
      kac(1:10*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,11
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*11+mb)
            p_bc = pbc((mc-1)*11+mb)
            DO ma = 1,10
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*10+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*10+ma)
              kad((md-1)*10+ma) =  kad((md-1)*10+ma)-tmp*p_bc
              kac((mc-1)*10+ma) = kac((mc-1)*10+ma)-tmp*p_bd
            END DO
            kbd((md-1)*11+mb) = kbd((md-1)*11+mb) - ks_bd
            kbc((mc-1)*11+mb) = kbc((mc-1)*11+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_10_11
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_10_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(15*md_max), kbc(15*mc_max), kad(10*md_max), &
      kac(10*mc_max), pbd(15*md_max), pbc(15*mc_max), pad(10*md_max), &
      pac(10*mc_max), prim(10*15*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:15*md_max) = 0.0_dp
      kbc(1:15*mc_max) = 0.0_dp
      kad(1:10*md_max) = 0.0_dp
      kac(1:10*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,15
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*15+mb)
            p_bc = pbc((mc-1)*15+mb)
            DO ma = 1,10
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*10+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*10+ma)
              kad((md-1)*10+ma) =  kad((md-1)*10+ma)-tmp*p_bc
              kac((mc-1)*10+ma) = kac((mc-1)*10+ma)-tmp*p_bd
            END DO
            kbd((md-1)*15+mb) = kbd((md-1)*15+mb) - ks_bd
            kbc((mc-1)*15+mb) = kbc((mc-1)*15+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_10_15
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_11_1_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*1), kad(11*1), kac(11*1), pbd(1*1), &
      pbc(1*1), pad(11*1), pac(11*1), prim(11*1*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:11*1) = 0.0_dp
      kac(1:11*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,11
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*11+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*11+ma)
              kad((md-1)*11+ma) =  kad((md-1)*11+ma)-tmp*p_bc
              kac((mc-1)*11+ma) = kac((mc-1)*11+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_11_1_1_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_11_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*1), kad(11*md_max), kac(11*1), &
      pbd(1*md_max), pbc(1*1), pad(11*md_max), pac(11*1), &
      prim(11*1*1*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:11*md_max) = 0.0_dp
      kac(1:11*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,11
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*11+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*11+ma)
              kad((md-1)*11+ma) =  kad((md-1)*11+ma)-tmp*p_bc
              kac((mc-1)*11+ma) = kac((mc-1)*11+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_11_1_1
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_11_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*mc_max), kad(11*md_max), &
      kac(11*mc_max), pbd(1*md_max), pbc(1*mc_max), pad(11*md_max), &
      pac(11*mc_max), prim(11*1*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*mc_max) = 0.0_dp
      kad(1:11*md_max) = 0.0_dp
      kac(1:11*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,11
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*11+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*11+ma)
              kad((md-1)*11+ma) =  kad((md-1)*11+ma)-tmp*p_bc
              kac((mc-1)*11+ma) = kac((mc-1)*11+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_11_1
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_11_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*mc_max), kad(11*md_max), &
      kac(11*mc_max), pbd(2*md_max), pbc(2*mc_max), pad(11*md_max), &
      pac(11*mc_max), prim(11*2*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*mc_max) = 0.0_dp
      kad(1:11*md_max) = 0.0_dp
      kac(1:11*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,11
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*11+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*11+ma)
              kad((md-1)*11+ma) =  kad((md-1)*11+ma)-tmp*p_bc
              kac((mc-1)*11+ma) = kac((mc-1)*11+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_11_2
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_11_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(3*md_max), kbc(3*mc_max), kad(11*md_max), &
      kac(11*mc_max), pbd(3*md_max), pbc(3*mc_max), pad(11*md_max), &
      pac(11*mc_max), prim(11*3*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*md_max) = 0.0_dp
      kbc(1:3*mc_max) = 0.0_dp
      kad(1:11*md_max) = 0.0_dp
      kac(1:11*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,11
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*11+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*11+ma)
              kad((md-1)*11+ma) =  kad((md-1)*11+ma)-tmp*p_bc
              kac((mc-1)*11+ma) = kac((mc-1)*11+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_11_3
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_11_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(4*md_max), kbc(4*mc_max), kad(11*md_max), &
      kac(11*mc_max), pbd(4*md_max), pbc(4*mc_max), pad(11*md_max), &
      pac(11*mc_max), prim(11*4*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*md_max) = 0.0_dp
      kbc(1:4*mc_max) = 0.0_dp
      kad(1:11*md_max) = 0.0_dp
      kac(1:11*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,11
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*11+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*11+ma)
              kad((md-1)*11+ma) =  kad((md-1)*11+ma)-tmp*p_bc
              kac((mc-1)*11+ma) = kac((mc-1)*11+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_11_4
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_11_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(5*md_max), kbc(5*mc_max), kad(11*md_max), &
      kac(11*mc_max), pbd(5*md_max), pbc(5*mc_max), pad(11*md_max), &
      pac(11*mc_max), prim(11*5*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:5*md_max) = 0.0_dp
      kbc(1:5*mc_max) = 0.0_dp
      kad(1:11*md_max) = 0.0_dp
      kac(1:11*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,5
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*5+mb)
            p_bc = pbc((mc-1)*5+mb)
            DO ma = 1,11
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*11+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*11+ma)
              kad((md-1)*11+ma) =  kad((md-1)*11+ma)-tmp*p_bc
              kac((mc-1)*11+ma) = kac((mc-1)*11+ma)-tmp*p_bd
            END DO
            kbd((md-1)*5+mb) = kbd((md-1)*5+mb) - ks_bd
            kbc((mc-1)*5+mb) = kbc((mc-1)*5+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_11_5
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_11_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(6*md_max), kbc(6*mc_max), kad(11*md_max), &
      kac(11*mc_max), pbd(6*md_max), pbc(6*mc_max), pad(11*md_max), &
      pac(11*mc_max), prim(11*6*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:6*md_max) = 0.0_dp
      kbc(1:6*mc_max) = 0.0_dp
      kad(1:11*md_max) = 0.0_dp
      kac(1:11*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,6
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*6+mb)
            p_bc = pbc((mc-1)*6+mb)
            DO ma = 1,11
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*11+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*11+ma)
              kad((md-1)*11+ma) =  kad((md-1)*11+ma)-tmp*p_bc
              kac((mc-1)*11+ma) = kac((mc-1)*11+ma)-tmp*p_bd
            END DO
            kbd((md-1)*6+mb) = kbd((md-1)*6+mb) - ks_bd
            kbc((mc-1)*6+mb) = kbc((mc-1)*6+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_11_6
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_11_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(7*md_max), kbc(7*mc_max), kad(11*md_max), &
      kac(11*mc_max), pbd(7*md_max), pbc(7*mc_max), pad(11*md_max), &
      pac(11*mc_max), prim(11*7*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:7*md_max) = 0.0_dp
      kbc(1:7*mc_max) = 0.0_dp
      kad(1:11*md_max) = 0.0_dp
      kac(1:11*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,7
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*7+mb)
            p_bc = pbc((mc-1)*7+mb)
            DO ma = 1,11
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*11+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*11+ma)
              kad((md-1)*11+ma) =  kad((md-1)*11+ma)-tmp*p_bc
              kac((mc-1)*11+ma) = kac((mc-1)*11+ma)-tmp*p_bd
            END DO
            kbd((md-1)*7+mb) = kbd((md-1)*7+mb) - ks_bd
            kbc((mc-1)*7+mb) = kbc((mc-1)*7+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_11_7
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_11_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(9*md_max), kbc(9*mc_max), kad(11*md_max), &
      kac(11*mc_max), pbd(9*md_max), pbc(9*mc_max), pad(11*md_max), &
      pac(11*mc_max), prim(11*9*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:9*md_max) = 0.0_dp
      kbc(1:9*mc_max) = 0.0_dp
      kad(1:11*md_max) = 0.0_dp
      kac(1:11*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,9
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*9+mb)
            p_bc = pbc((mc-1)*9+mb)
            DO ma = 1,11
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*11+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*11+ma)
              kad((md-1)*11+ma) =  kad((md-1)*11+ma)-tmp*p_bc
              kac((mc-1)*11+ma) = kac((mc-1)*11+ma)-tmp*p_bd
            END DO
            kbd((md-1)*9+mb) = kbd((md-1)*9+mb) - ks_bd
            kbc((mc-1)*9+mb) = kbc((mc-1)*9+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_11_9
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_11_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(10*md_max), kbc(10*mc_max), kad(11*md_max), &
      kac(11*mc_max), pbd(10*md_max), pbc(10*mc_max), pad(11*md_max), &
      pac(11*mc_max), prim(11*10*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:10*md_max) = 0.0_dp
      kbc(1:10*mc_max) = 0.0_dp
      kad(1:11*md_max) = 0.0_dp
      kac(1:11*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,10
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*10+mb)
            p_bc = pbc((mc-1)*10+mb)
            DO ma = 1,11
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*11+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*11+ma)
              kad((md-1)*11+ma) =  kad((md-1)*11+ma)-tmp*p_bc
              kac((mc-1)*11+ma) = kac((mc-1)*11+ma)-tmp*p_bd
            END DO
            kbd((md-1)*10+mb) = kbd((md-1)*10+mb) - ks_bd
            kbc((mc-1)*10+mb) = kbc((mc-1)*10+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_11_10
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_11_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(11*md_max), kbc(11*mc_max), kad(11*md_max), &
      kac(11*mc_max), pbd(11*md_max), pbc(11*mc_max), pad(11*md_max), &
      pac(11*mc_max), prim(11*11*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:11*md_max) = 0.0_dp
      kbc(1:11*mc_max) = 0.0_dp
      kad(1:11*md_max) = 0.0_dp
      kac(1:11*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,11
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*11+mb)
            p_bc = pbc((mc-1)*11+mb)
            DO ma = 1,11
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*11+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*11+ma)
              kad((md-1)*11+ma) =  kad((md-1)*11+ma)-tmp*p_bc
              kac((mc-1)*11+ma) = kac((mc-1)*11+ma)-tmp*p_bd
            END DO
            kbd((md-1)*11+mb) = kbd((md-1)*11+mb) - ks_bd
            kbc((mc-1)*11+mb) = kbc((mc-1)*11+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_11_11
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_11_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(15*md_max), kbc(15*mc_max), kad(11*md_max), &
      kac(11*mc_max), pbd(15*md_max), pbc(15*mc_max), pad(11*md_max), &
      pac(11*mc_max), prim(11*15*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:15*md_max) = 0.0_dp
      kbc(1:15*mc_max) = 0.0_dp
      kad(1:11*md_max) = 0.0_dp
      kac(1:11*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,15
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*15+mb)
            p_bc = pbc((mc-1)*15+mb)
            DO ma = 1,11
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*11+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*11+ma)
              kad((md-1)*11+ma) =  kad((md-1)*11+ma)-tmp*p_bc
              kac((mc-1)*11+ma) = kac((mc-1)*11+ma)-tmp*p_bd
            END DO
            kbd((md-1)*15+mb) = kbd((md-1)*15+mb) - ks_bd
            kbc((mc-1)*15+mb) = kbc((mc-1)*15+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_11_15
! *****************************************************************************
!> \brief ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_15_1_1_1(kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    REAL(KIND=dp) :: kbd(1*1), kbc(1*1), kad(15*1), kac(15*1), pbd(1*1), &
      pbc(1*1), pad(15*1), pac(15*1), prim(15*1*1*1), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*1) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:15*1) = 0.0_dp
      kac(1:15*1) = 0.0_dp
      p_index = 0
      DO md = 1,1
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,15
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*15+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*15+ma)
              kad((md-1)*15+ma) =  kad((md-1)*15+ma)-tmp*p_bc
              kac((mc-1)*15+ma) = kac((mc-1)*15+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_15_1_1_1
! *****************************************************************************
!> \brief ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_15_1_1(md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*1), kad(15*md_max), kac(15*1), &
      pbd(1*md_max), pbc(1*1), pad(15*md_max), pac(15*1), &
      prim(15*1*1*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*1) = 0.0_dp
      kad(1:15*md_max) = 0.0_dp
      kac(1:15*1) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,1
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,15
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*15+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*15+ma)
              kad((md-1)*15+ma) =  kad((md-1)*15+ma)-tmp*p_bc
              kac((mc-1)*15+ma) = kac((mc-1)*15+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_15_1_1
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_15_1(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(1*md_max), kbc(1*mc_max), kad(15*md_max), &
      kac(15*mc_max), pbd(1*md_max), pbc(1*mc_max), pad(15*md_max), &
      pac(15*mc_max), prim(15*1*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:1*md_max) = 0.0_dp
      kbc(1:1*mc_max) = 0.0_dp
      kad(1:15*md_max) = 0.0_dp
      kac(1:15*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,1
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*1+mb)
            p_bc = pbc((mc-1)*1+mb)
            DO ma = 1,15
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*15+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*15+ma)
              kad((md-1)*15+ma) =  kad((md-1)*15+ma)-tmp*p_bc
              kac((mc-1)*15+ma) = kac((mc-1)*15+ma)-tmp*p_bd
            END DO
            kbd((md-1)*1+mb) = kbd((md-1)*1+mb) - ks_bd
            kbc((mc-1)*1+mb) = kbc((mc-1)*1+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_15_1
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_15_2(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(2*md_max), kbc(2*mc_max), kad(15*md_max), &
      kac(15*mc_max), pbd(2*md_max), pbc(2*mc_max), pad(15*md_max), &
      pac(15*mc_max), prim(15*2*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:2*md_max) = 0.0_dp
      kbc(1:2*mc_max) = 0.0_dp
      kad(1:15*md_max) = 0.0_dp
      kac(1:15*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,2
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*2+mb)
            p_bc = pbc((mc-1)*2+mb)
            DO ma = 1,15
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*15+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*15+ma)
              kad((md-1)*15+ma) =  kad((md-1)*15+ma)-tmp*p_bc
              kac((mc-1)*15+ma) = kac((mc-1)*15+ma)-tmp*p_bd
            END DO
            kbd((md-1)*2+mb) = kbd((md-1)*2+mb) - ks_bd
            kbc((mc-1)*2+mb) = kbc((mc-1)*2+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_15_2
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_15_3(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(3*md_max), kbc(3*mc_max), kad(15*md_max), &
      kac(15*mc_max), pbd(3*md_max), pbc(3*mc_max), pad(15*md_max), &
      pac(15*mc_max), prim(15*3*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:3*md_max) = 0.0_dp
      kbc(1:3*mc_max) = 0.0_dp
      kad(1:15*md_max) = 0.0_dp
      kac(1:15*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,3
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*3+mb)
            p_bc = pbc((mc-1)*3+mb)
            DO ma = 1,15
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*15+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*15+ma)
              kad((md-1)*15+ma) =  kad((md-1)*15+ma)-tmp*p_bc
              kac((mc-1)*15+ma) = kac((mc-1)*15+ma)-tmp*p_bd
            END DO
            kbd((md-1)*3+mb) = kbd((md-1)*3+mb) - ks_bd
            kbc((mc-1)*3+mb) = kbc((mc-1)*3+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_15_3
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_15_4(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(4*md_max), kbc(4*mc_max), kad(15*md_max), &
      kac(15*mc_max), pbd(4*md_max), pbc(4*mc_max), pad(15*md_max), &
      pac(15*mc_max), prim(15*4*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:4*md_max) = 0.0_dp
      kbc(1:4*mc_max) = 0.0_dp
      kad(1:15*md_max) = 0.0_dp
      kac(1:15*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,4
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*4+mb)
            p_bc = pbc((mc-1)*4+mb)
            DO ma = 1,15
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*15+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*15+ma)
              kad((md-1)*15+ma) =  kad((md-1)*15+ma)-tmp*p_bc
              kac((mc-1)*15+ma) = kac((mc-1)*15+ma)-tmp*p_bd
            END DO
            kbd((md-1)*4+mb) = kbd((md-1)*4+mb) - ks_bd
            kbc((mc-1)*4+mb) = kbc((mc-1)*4+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_15_4
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_15_5(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(5*md_max), kbc(5*mc_max), kad(15*md_max), &
      kac(15*mc_max), pbd(5*md_max), pbc(5*mc_max), pad(15*md_max), &
      pac(15*mc_max), prim(15*5*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:5*md_max) = 0.0_dp
      kbc(1:5*mc_max) = 0.0_dp
      kad(1:15*md_max) = 0.0_dp
      kac(1:15*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,5
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*5+mb)
            p_bc = pbc((mc-1)*5+mb)
            DO ma = 1,15
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*15+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*15+ma)
              kad((md-1)*15+ma) =  kad((md-1)*15+ma)-tmp*p_bc
              kac((mc-1)*15+ma) = kac((mc-1)*15+ma)-tmp*p_bd
            END DO
            kbd((md-1)*5+mb) = kbd((md-1)*5+mb) - ks_bd
            kbc((mc-1)*5+mb) = kbc((mc-1)*5+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_15_5
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_15_6(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(6*md_max), kbc(6*mc_max), kad(15*md_max), &
      kac(15*mc_max), pbd(6*md_max), pbc(6*mc_max), pad(15*md_max), &
      pac(15*mc_max), prim(15*6*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:6*md_max) = 0.0_dp
      kbc(1:6*mc_max) = 0.0_dp
      kad(1:15*md_max) = 0.0_dp
      kac(1:15*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,6
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*6+mb)
            p_bc = pbc((mc-1)*6+mb)
            DO ma = 1,15
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*15+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*15+ma)
              kad((md-1)*15+ma) =  kad((md-1)*15+ma)-tmp*p_bc
              kac((mc-1)*15+ma) = kac((mc-1)*15+ma)-tmp*p_bd
            END DO
            kbd((md-1)*6+mb) = kbd((md-1)*6+mb) - ks_bd
            kbc((mc-1)*6+mb) = kbc((mc-1)*6+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_15_6
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_15_7(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(7*md_max), kbc(7*mc_max), kad(15*md_max), &
      kac(15*mc_max), pbd(7*md_max), pbc(7*mc_max), pad(15*md_max), &
      pac(15*mc_max), prim(15*7*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:7*md_max) = 0.0_dp
      kbc(1:7*mc_max) = 0.0_dp
      kad(1:15*md_max) = 0.0_dp
      kac(1:15*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,7
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*7+mb)
            p_bc = pbc((mc-1)*7+mb)
            DO ma = 1,15
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*15+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*15+ma)
              kad((md-1)*15+ma) =  kad((md-1)*15+ma)-tmp*p_bc
              kac((mc-1)*15+ma) = kac((mc-1)*15+ma)-tmp*p_bd
            END DO
            kbd((md-1)*7+mb) = kbd((md-1)*7+mb) - ks_bd
            kbc((mc-1)*7+mb) = kbc((mc-1)*7+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_15_7
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_15_9(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(9*md_max), kbc(9*mc_max), kad(15*md_max), &
      kac(15*mc_max), pbd(9*md_max), pbc(9*mc_max), pad(15*md_max), &
      pac(15*mc_max), prim(15*9*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:9*md_max) = 0.0_dp
      kbc(1:9*mc_max) = 0.0_dp
      kad(1:15*md_max) = 0.0_dp
      kac(1:15*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,9
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*9+mb)
            p_bc = pbc((mc-1)*9+mb)
            DO ma = 1,15
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*15+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*15+ma)
              kad((md-1)*15+ma) =  kad((md-1)*15+ma)-tmp*p_bc
              kac((mc-1)*15+ma) = kac((mc-1)*15+ma)-tmp*p_bd
            END DO
            kbd((md-1)*9+mb) = kbd((md-1)*9+mb) - ks_bd
            kbc((mc-1)*9+mb) = kbc((mc-1)*9+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_15_9
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_15_10(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(10*md_max), kbc(10*mc_max), kad(15*md_max), &
      kac(15*mc_max), pbd(10*md_max), pbc(10*mc_max), pad(15*md_max), &
      pac(15*mc_max), prim(15*10*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:10*md_max) = 0.0_dp
      kbc(1:10*mc_max) = 0.0_dp
      kad(1:15*md_max) = 0.0_dp
      kac(1:15*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,10
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*10+mb)
            p_bc = pbc((mc-1)*10+mb)
            DO ma = 1,15
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*15+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*15+ma)
              kad((md-1)*15+ma) =  kad((md-1)*15+ma)-tmp*p_bc
              kac((mc-1)*15+ma) = kac((mc-1)*15+ma)-tmp*p_bd
            END DO
            kbd((md-1)*10+mb) = kbd((md-1)*10+mb) - ks_bd
            kbc((mc-1)*10+mb) = kbc((mc-1)*10+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_15_10
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_15_11(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(11*md_max), kbc(11*mc_max), kad(15*md_max), &
      kac(15*mc_max), pbd(11*md_max), pbc(11*mc_max), pad(15*md_max), &
      pac(15*mc_max), prim(15*11*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:11*md_max) = 0.0_dp
      kbc(1:11*mc_max) = 0.0_dp
      kad(1:15*md_max) = 0.0_dp
      kac(1:15*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,11
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*11+mb)
            p_bc = pbc((mc-1)*11+mb)
            DO ma = 1,15
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*15+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*15+ma)
              kad((md-1)*15+ma) =  kad((md-1)*15+ma)-tmp*p_bc
              kac((mc-1)*15+ma) = kac((mc-1)*15+ma)-tmp*p_bd
            END DO
            kbd((md-1)*11+mb) = kbd((md-1)*11+mb) - ks_bd
            kbc((mc-1)*11+mb) = kbc((mc-1)*11+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_15_11
! *****************************************************************************
!> \brief ...
!> \param mc_max ...
!> \param md_max ...
!> \param kbd ...
!> \param kbc ...
!> \param kad ...
!> \param kac ...
!> \param pbd ...
!> \param pbc ...
!> \param pad ...
!> \param pac ...
!> \param prim ...
!> \param scale ...
! *****************************************************************************
  SUBROUTINE block_15_15(mc_max,md_max,kbd,kbc,kad,kac,pbd,pbc,pad,pac,prim,scale)
    INTEGER                                  :: mc_max, md_max
    REAL(KIND=dp) :: kbd(15*md_max), kbc(15*mc_max), kad(15*md_max), &
      kac(15*mc_max), pbd(15*md_max), pbc(15*mc_max), pad(15*md_max), &
      pac(15*mc_max), prim(15*15*mc_max*md_max), scale

    INTEGER                                  :: ma, mb, mc, md, p_index
    REAL(KIND=dp)                            :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:15*md_max) = 0.0_dp
      kbc(1:15*mc_max) = 0.0_dp
      kad(1:15*md_max) = 0.0_dp
      kac(1:15*mc_max) = 0.0_dp
      p_index = 0
      DO md = 1,md_max
        DO mc = 1,mc_max
          DO mb = 1,15
            ks_bd = 0.0_dp
            ks_bc = 0.0_dp
            p_bd = pbd((md-1)*15+mb)
            p_bc = pbc((mc-1)*15+mb)
            DO ma = 1,15
              p_index=p_index+1
              tmp = scale*prim(p_index)
              ks_bc = ks_bc + tmp* pad((md-1)*15+ma)
              ks_bd = ks_bd + tmp* pac((mc-1)*15+ma)
              kad((md-1)*15+ma) =  kad((md-1)*15+ma)-tmp*p_bc
              kac((mc-1)*15+ma) = kac((mc-1)*15+ma)-tmp*p_bd
            END DO
            kbd((md-1)*15+mb) = kbd((md-1)*15+mb) - ks_bd
            kbc((mc-1)*15+mb) = kbc((mc-1)*15+mb) - ks_bc
          END DO
        END DO
      END DO
  END SUBROUTINE block_15_15
#endif
END MODULE hfx_contract_block
