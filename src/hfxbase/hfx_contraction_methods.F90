!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Contains routines for contraction without dgemms.
!> \notes Contains specific routines for contraction. The compiler flag
!>        -D__MAX_CONTR defines the maximum angular momentum up to which
!>        specialized code will be used. Default setting is d-functions.
!>        Increasing -D__MAX_CONTR produces faster code but might overburden
!>        the optimization capabilities of some poor compilers.
!>        Templates are available up to l=10. The Fypp parameters lmax_specific
!>        and lmax_generic control the generated specialized routines and sparse
!>        cases in the generic routine; both default to 4 (g-functions).
!>        Fypp preserves compile-time dimensions and explicit sparse sums.
!> \par History
!>      07.2009 created [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************

MODULE hfx_contraction_methods

!** This defines the default behaviour
#ifndef __MAX_CONTR
#define __MAX_CONTR 2
#endif

   USE kinds, ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: contract

CONTAINS

   #! Generation limits, independent of the compiler's __MAX_CONTR setting.
   #:set lmax_specific = 4
   #:set lmax_generic = 4

   #! Shell names and dimensions from the legacy generator, for l=0 through 10.
   #:set shells = 'spdfghijklm'
   #:set cartesian = [(l + 1)*(l + 2)//2 for l in range(len(shells))]
   #:set spherical = [2*l + 1 for l in range(len(shells))]
   #:assert 0 <= lmax_specific < len(shells)
   #:assert 0 <= lmax_generic < len(shells)
   #! Nonzero (Cartesian, spherical) indices, in the original accumulation order.
   #:set terms = [ &
      [(1, 1)], &
      [(1, 3), (2, 1), (3, 2)], &
      [(1, 3), (1, 5), (2, 1), (3, 4), (4, 3), (4, 5), (5, 2), (6, 3)], &
      [(1, 5), (1, 7), (2, 1), (2, 3), (3, 4), (3, 6), (4, 5), (4, 7), &
       (5, 2), (6, 5), (7, 1), (7, 3), (8, 4), (8, 6), (9, 3), (10, 4)], &
      [(1, 5), (1, 7), (1, 9), (2, 1), (2, 3), (3, 6), (3, 8), (4, 5), &
       (4, 9), (5, 2), (5, 4), (6, 5), (6, 7), (7, 1), (7, 3), (8, 6), &
       (8, 8), (9, 3), (10, 6), (11, 5), (11, 7), (11, 9), (12, 2), &
       (12, 4), (13, 5), (13, 7), (14, 4), (15, 5)], &
      [(1, 7), (1, 9), (1, 11), (2, 1), (2, 3), (2, 5), (3, 6), (3, 8), &
       (3, 10), (4, 7), (4, 9), (4, 11), (5, 2), (5, 4), (6, 7), (6, 9), &
       (7, 1), (7, 3), (7, 5), (8, 6), (8, 10), (9, 3), (9, 5), (10, 6), &
       (10, 8), (11, 7), (11, 9), (11, 11), (12, 2), (12, 4), (13, 7), (13, 9), &
       (14, 4), (15, 7), (16, 1), (16, 3), (16, 5), (17, 6), (17, 8), (17, 10), &
       (18, 3), (18, 5), (19, 6), (19, 8), (20, 5), (21, 6)], &
      [(1, 7), (1, 9), (1, 11), (1, 13), (2, 1), (2, 3), (2, 5), (3, 8), &
       (3, 10), (3, 12), (4, 7), (4, 9), (4, 11), (4, 13), (5, 2), (5, 4), &
       (5, 6), (6, 7), (6, 9), (6, 11), (7, 1), (7, 5), (8, 8), (8, 10), &
       (8, 12), (9, 3), (9, 5), (10, 8), (10, 10), (11, 7), (11, 9), (11, 11), &
       (11, 13), (12, 2), (12, 4), (12, 6), (13, 7), (13, 11), (14, 4), (14, 6), &
       (15, 7), (15, 9), (16, 1), (16, 3), (16, 5), (17, 8), (17, 10), (17, 12), &
       (18, 3), (18, 5), (19, 8), (19, 10), (20, 5), (21, 8), (22, 7), (22, 9), &
       (22, 11), (22, 13), (23, 2), (23, 4), (23, 6), (24, 7), (24, 9), (24, 11), &
       (25, 4), (25, 6), (26, 7), (26, 9), (27, 6), (28, 7)], &
      [(1, 9), (1, 11), (1, 13), (1, 15), (2, 1), (2, 3), (2, 5), (2, 7), &
       (3, 8), (3, 10), (3, 12), (3, 14), (4, 9), (4, 11), (4, 13), (4, 15), &
       (5, 2), (5, 4), (5, 6), (6, 9), (6, 11), (6, 13), (7, 1), (7, 3), &
       (7, 5), (7, 7), (8, 8), (8, 10), (8, 12), (8, 14), (9, 3), (9, 5), &
       (9, 7), (10, 8), (10, 10), (10, 12), (11, 9), (11, 11), (11, 13), (11, 15), &
       (12, 2), (12, 6), (13, 9), (13, 11), (13, 13), (14, 4), (14, 6), (15, 9), &
       (15, 11), (16, 1), (16, 3), (16, 5), (16, 7), (17, 8), (17, 10), (17, 12), &
       (17, 14), (18, 3), (18, 5), (18, 7), (19, 8), (19, 12), (20, 5), (20, 7), &
       (21, 8), (21, 10), (22, 9), (22, 11), (22, 13), (22, 15), (23, 2), (23, 4), &
       (23, 6), (24, 9), (24, 11), (24, 13), (25, 4), (25, 6), (26, 9), (26, 11), &
       (27, 6), (28, 9), (29, 1), (29, 3), (29, 5), (29, 7), (30, 8), (30, 10), &
       (30, 12), (30, 14), (31, 3), (31, 5), (31, 7), (32, 8), (32, 10), (32, 12), &
       (33, 5), (33, 7), (34, 8), (34, 10), (35, 7), (36, 8)], &
      [(1, 9), (1, 11), (1, 13), (1, 15), (1, 17), (2, 1), (2, 3), (2, 5), &
       (2, 7), (3, 10), (3, 12), (3, 14), (3, 16), (4, 9), (4, 11), (4, 13), &
       (4, 15), (4, 17), (5, 2), (5, 4), (5, 6), (5, 8), (6, 9), (6, 11), &
       (6, 13), (6, 15), (7, 1), (7, 3), (7, 5), (7, 7), (8, 10), (8, 12), &
       (8, 14), (8, 16), (9, 3), (9, 5), (9, 7), (10, 10), (10, 12), (10, 14), &
       (11, 9), (11, 13), (11, 17), (12, 2), (12, 4), (12, 6), (12, 8), (13, 9), &
       (13, 11), (13, 13), (13, 15), (14, 4), (14, 6), (14, 8), (15, 9), (15, 11), &
       (15, 13), (16, 1), (16, 3), (16, 5), (16, 7), (17, 10), (17, 12), (17, 14), &
       (17, 16), (18, 3), (18, 7), (19, 10), (19, 12), (19, 14), (20, 5), (20, 7), &
       (21, 10), (21, 12), (22, 9), (22, 11), (22, 13), (22, 15), (22, 17), (23, 2), &
       (23, 4), (23, 6), (23, 8), (24, 9), (24, 11), (24, 13), (24, 15), (25, 4), &
       (25, 6), (25, 8), (26, 9), (26, 13), (27, 6), (27, 8), (28, 9), (28, 11), &
       (29, 1), (29, 3), (29, 5), (29, 7), (30, 10), (30, 12), (30, 14), (30, 16), &
       (31, 3), (31, 5), (31, 7), (32, 10), (32, 12), (32, 14), (33, 5), (33, 7), &
       (34, 10), (34, 12), (35, 7), (36, 10), (37, 9), (37, 11), (37, 13), (37, 15), &
       (37, 17), (38, 2), (38, 4), (38, 6), (38, 8), (39, 9), (39, 11), (39, 13), &
       (39, 15), (40, 4), (40, 6), (40, 8), (41, 9), (41, 11), (41, 13), (42, 6), &
       (42, 8), (43, 9), (43, 11), (44, 8), (45, 9)], &
      [(1, 11), (1, 13), (1, 15), (1, 17), (1, 19), (2, 1), (2, 3), (2, 5), &
       (2, 7), (2, 9), (3, 10), (3, 12), (3, 14), (3, 16), (3, 18), (4, 11), &
       (4, 15), (4, 17), (4, 19), (5, 2), (5, 4), (5, 6), (5, 8), (6, 11), &
       (6, 13), (6, 15), (6, 17), (7, 1), (7, 3), (7, 7), (7, 9), (8, 10), &
       (8, 12), (8, 14), (8, 16), (8, 18), (9, 3), (9, 5), (9, 7), (9, 9), &
       (10, 10), (10, 12), (10, 14), (10, 16), (11, 11), (11, 13), (11, 15), (11, 17), &
       (11, 19), (12, 2), (12, 4), (12, 6), (12, 8), (13, 11), (13, 13), (13, 15), &
       (13, 17), (14, 4), (14, 6), (14, 8), (15, 11), (15, 13), (15, 15), (16, 1), &
       (16, 3), (16, 5), (16, 7), (16, 9), (17, 10), (17, 14), (17, 18), (18, 3), &
       (18, 5), (18, 7), (18, 9), (19, 10), (19, 12), (19, 14), (19, 16), (20, 5), &
       (20, 7), (20, 9), (21, 10), (21, 12), (21, 14), (22, 11), (22, 13), (22, 17), &
       (22, 19), (23, 2), (23, 4), (23, 6), (23, 8), (24, 11), (24, 13), (24, 15), &
       (24, 17), (25, 4), (25, 8), (26, 11), (26, 13), (26, 15), (27, 6), (27, 8), &
       (28, 11), (28, 13), (29, 1), (29, 3), (29, 5), (29, 9), (30, 10), (30, 12), &
       (30, 14), (30, 16), (30, 18), (31, 3), (31, 5), (31, 7), (31, 9), (32, 10), &
       (32, 12), (32, 14), (32, 16), (33, 5), (33, 7), (33, 9), (34, 10), (34, 14), &
       (35, 7), (35, 9), (36, 10), (36, 12), (37, 11), (37, 13), (37, 15), (37, 17), &
       (37, 19), (38, 2), (38, 4), (38, 6), (38, 8), (39, 11), (39, 13), (39, 15), &
       (39, 17), (40, 4), (40, 6), (40, 8), (41, 11), (41, 13), (41, 15), (42, 6), &
       (42, 8), (43, 11), (43, 13), (44, 8), (45, 11), (46, 1), (46, 3), (46, 5), &
       (46, 7), (46, 9), (47, 10), (47, 12), (47, 14), (47, 16), (47, 18), (48, 3), &
       (48, 5), (48, 7), (48, 9), (49, 10), (49, 12), (49, 14), (49, 16), (50, 5), &
       (50, 7), (50, 9), (51, 10), (51, 12), (51, 14), (52, 7), (52, 9), (53, 10), &
       (53, 12), (54, 9), (55, 10)], &
      [(1, 11), (1, 13), (1, 15), (1, 17), (1, 19), (1, 21), (2, 1), (2, 3), &
       (2, 5), (2, 7), (2, 9), (3, 12), (3, 14), (3, 16), (3, 18), (3, 20), &
       (4, 11), (4, 13), (4, 15), (4, 17), (4, 19), (4, 21), (5, 2), (5, 4), &
       (5, 6), (5, 8), (5, 10), (6, 11), (6, 13), (6, 15), (6, 17), (6, 19), &
       (7, 1), (7, 3), (7, 5), (7, 7), (7, 9), (8, 12), (8, 16), (8, 18), &
       (8, 20), (9, 3), (9, 5), (9, 7), (9, 9), (10, 12), (10, 14), (10, 16), &
       (10, 18), (11, 11), (11, 13), (11, 15), (11, 17), (11, 19), (11, 21), (12, 2), &
       (12, 4), (12, 8), (12, 10), (13, 11), (13, 13), (13, 15), (13, 17), (13, 19), &
       (14, 4), (14, 6), (14, 8), (14, 10), (15, 11), (15, 13), (15, 15), (15, 17), &
       (16, 1), (16, 5), (16, 9), (17, 12), (17, 14), (17, 16), (17, 18), (17, 20), &
       (18, 3), (18, 5), (18, 7), (18, 9), (19, 12), (19, 14), (19, 16), (19, 18), &
       (20, 5), (20, 7), (20, 9), (21, 12), (21, 14), (21, 16), (22, 11), (22, 13), &
       (22, 15), (22, 17), (22, 19), (22, 21), (23, 2), (23, 4), (23, 6), (23, 8), &
       (23, 10), (24, 11), (24, 15), (24, 19), (25, 4), (25, 6), (25, 8), (25, 10), &
       (26, 11), (26, 13), (26, 15), (26, 17), (27, 6), (27, 8), (27, 10), (28, 11), &
       (28, 13), (28, 15), (29, 1), (29, 3), (29, 5), (29, 7), (29, 9), (30, 12), &
       (30, 14), (30, 18), (30, 20), (31, 3), (31, 5), (31, 7), (31, 9), (32, 12), &
       (32, 14), (32, 16), (32, 18), (33, 5), (33, 9), (34, 12), (34, 14), (34, 16), &
       (35, 7), (35, 9), (36, 12), (36, 14), (37, 11), (37, 13), (37, 15), (37, 17), &
       (37, 19), (37, 21), (38, 2), (38, 4), (38, 6), (38, 10), (39, 11), (39, 13), &
       (39, 15), (39, 17), (39, 19), (40, 4), (40, 6), (40, 8), (40, 10), (41, 11), &
       (41, 13), (41, 15), (41, 17), (42, 6), (42, 8), (42, 10), (43, 11), (43, 15), &
       (44, 8), (44, 10), (45, 11), (45, 13), (46, 1), (46, 3), (46, 5), (46, 7), &
       (46, 9), (47, 12), (47, 14), (47, 16), (47, 18), (47, 20), (48, 3), (48, 5), &
       (48, 7), (48, 9), (49, 12), (49, 14), (49, 16), (49, 18), (50, 5), (50, 7), &
       (50, 9), (51, 12), (51, 14), (51, 16), (52, 7), (52, 9), (53, 12), (53, 14), &
       (54, 9), (55, 12), (56, 11), (56, 13), (56, 15), (56, 17), (56, 19), (56, 21), &
       (57, 2), (57, 4), (57, 6), (57, 8), (57, 10), (58, 11), (58, 13), (58, 15), &
       (58, 17), (58, 19), (59, 4), (59, 6), (59, 8), (59, 10), (60, 11), (60, 13), &
       (60, 15), (60, 17), (61, 6), (61, 8), (61, 10), (62, 11), (62, 13), (62, 15), &
       (63, 8), (63, 10), (64, 11), (64, 13), (65, 10), (66, 11)]]

   #:def call_generic()
      CALL contract_generic(ncoa, ncob, ncoc, ncod, nsoa, nsob, nsoc, nsod, &
                            work, nl_a, nl_b, nl_c, nl_d, &
                            sphi_a, sphi_b, sphi_c, sphi_d, primitives, buffer1, buffer2)
   #:enddef

   #! The first three stages transpose into a work buffer; the last accumulates
   #! into the four-dimensional output. All sparse terms are expanded by Fypp.
   #:def transform(angular, stage, ns)
      #:set axis = 'abcd'[stage]
      #:set source = ['work', 'buffer1', 'buffer2', 'buffer1'][stage]
      #:if stage < 3
         #:set target = ['buffer1', 'buffer2', 'buffer1'][stage]
         DO i = 1, imax
            #:for k, j in terms[angular]
               ${target}$ (i + imax*(${j}$-1)) = ${target}$ (i + imax*(${j}$-1)) &
                                                 + ${source}$ (${k}$+(i - 1)*kmax)*sphi_${axis}$ (${k}$, ${j}$+s_offset_${axis}$1)
            #:endfor
         END DO
      #:else
         DO i1 = 1, ${ns[2]}$
         DO i2 = 1, ${ns[1]}$
         DO i3 = 1, ${ns[0]}$
            i = i + 1
            #:for k, j in terms[angular]
               primitives(s_offset_a1 + i3, s_offset_b1 + i2, s_offset_c1 + i1, s_offset_d1 + ${j}$) = &
                  primitives(s_offset_a1 + i3, s_offset_b1 + i2, s_offset_c1 + i1, s_offset_d1 + ${j}$) &
                  + buffer1(${k}$+(i - 1)*kmax)*sphi_d(${k}$, ${j}$+s_offset_d1)
            #:endfor
         END DO
         END DO
         END DO
      #:endif
   #:enddef

   #! Share the loop structure between all specialized kernels and the fallback.
   #! Runtime dimensions and dense fallback loops are emitted only for generic.
   #:def contraction(angular=None)
      #:set generic = angular is None
      #:if generic
         #:set name = 'generic'
         #:set nc = ['nco' + axis for axis in 'abcd']
         #:set ns = ['nso' + axis for axis in 'abcd']
      #:else
         #:set name = ''.join(shells[l] for l in angular)
         #:set nc = [str(cartesian[l]) for l in angular]
         #:set ns = [str(spherical[l]) for l in angular]
         #:set max_contr = max(angular)
#if __MAX_CONTR >= ${max_contr}$
      #:endif
! **************************************************************************************************
!> \brief Contract Cartesian integrals into spherical basis functions.
      #:if generic
         #:for axis in 'abcd'
!> \param nco${axis}$ Number of Cartesian functions on center ${axis}$
         #:endfor
         #:for axis in 'abcd'
!> \param nso${axis}$ Number of spherical functions on center ${axis}$
         #:endfor
      #:endif
!> \param work Cartesian integrals
      #:for axis in 'abcd'
!> \param nl_${axis}$ Number of contractions on center ${axis}$
      #:endfor
      #:for axis in 'abcd'
!> \param sphi_${axis}$ Contraction and spherical transformation coefficients
      #:endfor
!> \param primitives Accumulated spherical integrals
!> \param buffer1 Work buffer
!> \param buffer2 Work buffer
! **************************************************************************************************
      #:if generic
         SUBROUTINE contract_generic(ncoa, ncob, ncoc, ncod, nsoa, nsob, nsoc, nsod, &
                                     work, nl_a, nl_b, nl_c, nl_d, &
                                     sphi_a, sphi_b, sphi_c, sphi_d, primitives, buffer1, buffer2)

            INTEGER, INTENT(IN) :: ncoa, ncob, ncoc, ncod, nsoa, nsob, nsoc, nsod
            #:else
            SUBROUTINE contract_${name}$ (work, nl_a, nl_b, nl_c, nl_d, &
                                          sphi_a, sphi_b, sphi_c, sphi_d, primitives, buffer1, buffer2)

               #:endif
               REAL(dp), DIMENSION(${ '*'.join(nc) }$), INTENT(IN) :: work
               INTEGER :: nl_a, nl_b, nl_c, nl_d
               #:for axis, ncart, nsph in zip('abcd', nc, ns)
                  REAL(dp), DIMENSION(${ncart}$, ${nsph}$*nl_${axis}$), INTENT(IN) :: sphi_${axis}$
               #:endfor
               REAL(dp), DIMENSION(${ns[0]}$*nl_a, ${ns[1]}$*nl_b, ${ns[2]}$*nl_c, ${ns[3]}$*nl_d) :: primitives
               REAL(dp), DIMENSION(${ '*'.join(nc) }$) :: buffer1, buffer2

               INTEGER :: i, i1, i2, i3, ia, ib, ic, id, imax, &
   #:if generic
                          j, jmax, k, &
   #:endif
                          kmax, s_offset_a1, s_offset_b1, s_offset_c1, s_offset_d1

               s_offset_a1 = 0
               DO ia = 1, nl_a
                  s_offset_b1 = 0
                  DO ib = 1, nl_b
                     s_offset_c1 = 0
                     DO ic = 1, nl_c
                        s_offset_d1 = 0
                        DO id = 1, nl_d
                           #:for stage, axis in enumerate('abcd')
                              #:set source = ['work', 'buffer1', 'buffer2', 'buffer1'][stage]
                              #:if stage < 3
                                 #:set target = ['buffer1', 'buffer2', 'buffer1'][stage]
                                 ${target}$ = 0.0_dp
                              #:endif
                              imax = ${ '*'.join(ns[:stage] + nc[stage + 1:]) }$
                              #:if generic
                                 jmax = ${ns[stage]}$
                              #:endif
                              kmax = ${nc[stage]}$
                              #:if stage == 3
                                 i = 0
                              #:endif
                              #:if generic
                                 SELECT CASE (kmax)
                                    #:for l, ncart in enumerate(cartesian[:lmax_generic + 1])
                                       CASE (${ncart}$)
                                       $:transform(l, stage, ns)
                                    #:endfor
                                 CASE DEFAULT
                                    DO j = 1, jmax
                                       #:if stage < 3
                                          DO i = 1, imax
                                             DO k = 1, kmax
                                                ${target}$ (i + imax*(j - 1)) = ${target}$ (i + imax*(j - 1)) &
                                                           + ${source}$ (k + (i - 1)*kmax)*sphi_${axis}$ (k, j + s_offset_${axis}$1)
                                             END DO
                                          END DO
                                       #:else
                                          i = 0
                                          DO i1 = 1, nsoc
                                          DO i2 = 1, nsob
                                          DO i3 = 1, nsoa
                                             i = i + 1
                                             DO k = 1, kmax
                                               primitives(s_offset_a1 + i3, s_offset_b1 + i2, s_offset_c1 + i1, s_offset_d1 + j) = &
                                                 primitives(s_offset_a1 + i3, s_offset_b1 + i2, s_offset_c1 + i1, s_offset_d1 + j) &
                                                   + buffer1(k + (i - 1)*kmax)*sphi_d(k, j + s_offset_d1)
                                             END DO
                                          END DO
                                          END DO
                                          END DO
                                       #:endif
                                    END DO
                                 END SELECT
                              #:else
                                 $:transform(angular[stage], stage, ns)
                              #:endif
                           #:endfor
                           s_offset_d1 = s_offset_d1 + ${ns[3]}$
                        END DO
                        s_offset_c1 = s_offset_c1 + ${ns[2]}$
                     END DO
                     s_offset_b1 = s_offset_b1 + ${ns[1]}$
                  END DO
                  s_offset_a1 = s_offset_a1 + ${ns[0]}$
               END DO
            END SUBROUTINE contract_${name}$
            #:if not generic
#endif
            #:endif
            #:enddef

! **************************************************************************************************
!> \brief ...
!> \param ncoa ...
!> \param ncob ...
!> \param ncoc ...
!> \param ncod ...
!> \param nsoa ...
!> \param nsob ...
!> \param nsoc ...
!> \param nsod ...
!> \param n_a ...
!> \param n_b ...
!> \param n_c ...
!> \param n_d ...
!> \param nl_a ...
!> \param nl_b ...
!> \param nl_c ...
!> \param nl_d ...
!> \param work ...
!> \param sphi_a ...
!> \param sphi_b ...
!> \param sphi_c ...
!> \param sphi_d ...
!> \param primitives ...
!> \param buffer1 ...
!> \param buffer2 ...
! **************************************************************************************************
            SUBROUTINE contract(ncoa, ncob, ncoc, ncod, nsoa, nsob, nsoc, nsod, &
                                n_a, n_b, n_c, n_d, nl_a, nl_b, nl_c, nl_d, work, &
                                sphi_a, sphi_b, sphi_c, sphi_d, &
                                primitives, &
                                buffer1, buffer2)

               INTEGER, INTENT(IN)         :: ncoa, ncob, ncoc, ncod, nsoa, nsob, nsoc, nsod, &
                                              n_a, n_b, n_c, n_d, nl_a, nl_b, nl_c, nl_d
               REAL(dp), DIMENSION(ncoa*ncob*ncoc*ncod), INTENT(IN) :: work
               REAL(dp), DIMENSION(ncoa, nsoa*nl_a), INTENT(IN)   :: sphi_a
               REAL(dp), DIMENSION(ncob, nsob*nl_b), INTENT(IN)   :: sphi_b
               REAL(dp), DIMENSION(ncoc, nsoc*nl_c), INTENT(IN)   :: sphi_c
               REAL(dp), DIMENSION(ncod, nsod*nl_d), INTENT(IN)   :: sphi_d

               REAL(dp), DIMENSION(nsoa*nl_a, nsob*nl_b, nsoc*nl_c, nsod*nl_d) :: primitives
               REAL(dp), DIMENSION(ncoa*ncob*ncoc*ncod)  :: buffer1, buffer2

#if !defined(__LIBINT)
               MARK_USED(ncoa)
               MARK_USED(ncob)
               MARK_USED(ncoc)
               MARK_USED(ncod)
               MARK_USED(nsoa)
               MARK_USED(nsob)
               MARK_USED(nsoc)
               MARK_USED(nsod)
               MARK_USED(n_a)
               MARK_USED(n_b)
               MARK_USED(n_c)
               MARK_USED(n_d)
               MARK_USED(nl_a)
               MARK_USED(nl_b)
               MARK_USED(nl_c)
               MARK_USED(nl_d)
               MARK_USED(sphi_a)
               MARK_USED(sphi_b)
               MARK_USED(sphi_c)
               MARK_USED(sphi_d)
               MARK_USED(work)
               MARK_USED(primitives)
               MARK_USED(buffer1)
               MARK_USED(buffer2)
               CPABORT("libint not compiled in")
#else
               SELECT CASE (n_a)
                  #:for a in range(lmax_specific + 1)
                     CASE (${a}$)
                     SELECT CASE (n_b)
                        #:for b in range(lmax_specific + 1)
                           CASE (${b}$)
                           SELECT CASE (n_c)
                              #:for c in range(lmax_specific + 1)
                                 CASE (${c}$)
                                 SELECT CASE (n_d)
                                    #:for d in range(lmax_specific + 1)
                                       CASE (${d}$)
                                       #:set max_contr = max(a, b, c, d)
#if __MAX_CONTR >= ${max_contr}$
                                       CALL contract_${''.join(shells[l] for l in (a, b, c, d))}$ (work, &
                                                                                                   nl_a, nl_b, nl_c, nl_d, &
                                                                       sphi_a, sphi_b, sphi_c, sphi_d, primitives, buffer1, buffer2)
#else
                                       $:call_generic()
#endif
                                    #:endfor
                                 CASE DEFAULT
                                    $:call_generic()
                                 END SELECT
                              #:endfor
                           CASE DEFAULT
                              $:call_generic()
                           END SELECT
                        #:endfor
                     CASE DEFAULT
                        $:call_generic()
                     END SELECT
                  #:endfor
               CASE DEFAULT
                  $:call_generic()
               END SELECT
#endif
            END SUBROUTINE contract

#if defined(__LIBINT)
            #:for a in range(lmax_specific + 1)
               #:for b in range(lmax_specific + 1)
                  #:for c in range(lmax_specific + 1)
                     #:for d in range(lmax_specific + 1)
                        $:contraction((a, b, c, d))
                     #:endfor
                  #:endfor
               #:endfor
            #:endfor

            $:contraction()
#endif
            END MODULE hfx_contraction_methods
