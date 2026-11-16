!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief routines to contract density matrix blocks with the four-center
!>        integrals to yield the Kohn-Sham matrix. The specialized routines
!>        are about 1.2-2.0 as fast as the default one.
!> \par History
!>      10.2009 created [Joost VandeVondele]
!> \author Joost VandeVondele
! **************************************************************************************************
MODULE hfx_contract_block
   USE kinds, ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_contract_block'
   PUBLIC :: contract_block
CONTAINS

   #! Keep the original dimensions and limit for further specialization.
   #:set block_sizes = [1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 15]
   #:set max_product = 20
   #:set block_args = ['kbd', 'kbc', 'kad', 'kac', 'pbd', 'pbc', 'pad', 'pac', 'prim', 'scale']

   #:def call_block(fixed=())
      #:set name = 'block_' + ('_'.join(map(str, fixed)) if fixed else 'default')
      #:set runtime_dims = ['m' + axis + '_max' for axis in 'abcd'[len(fixed):]]
      CALL ${name}$ (${', '.join(runtime_dims + block_args)}$)
   #:enddef

   #! The same loop body generates the default kernel and every specialization.
   #! Fixed dimensions remain literal constants in the generated Fortran code.
   #:def block(fixed=())
      #:set name = 'block_' + ('_'.join(map(str, fixed)) if fixed else 'default')
      #:set runtime_dims = ['m' + axis + '_max' for axis in 'abcd'[len(fixed):]]
      #:set a, b, c, d = list(map(str, fixed)) + runtime_dims
! **************************************************************************************************
!> \brief Contract density matrix blocks with four-center integrals.
      #:for dim in runtime_dims
!> \param ${dim}$ Block dimension
      #:endfor
      #:for arg in block_args
!> \param ${arg}$ ...
      #:endfor
! **************************************************************************************************
      SUBROUTINE ${name}$ (${', '.join(runtime_dims + block_args)}$)
      #:if runtime_dims
         INTEGER :: ${', '.join(runtime_dims)}$
      #:endif
      REAL(KIND=dp) :: kbd(${b}$*${d}$), kbc(${b}$*${c}$), kad(${a}$*${d}$), &
                       kac(${a}$*${c}$), pbd(${b}$*${d}$), pbc(${b}$*${c}$), pad(${a}$*${d}$), &
                       pac(${a}$*${c}$), prim(${a}$*${b}$*${c}$*${d}$), scale

      INTEGER :: ma, mb, mc, md, p_index
      REAL(KIND=dp) :: ks_bc, ks_bd, p_bc, p_bd, tmp

      kbd(1:${b}$*${d}$) = 0.0_dp
      kbc(1:${b}$*${c}$) = 0.0_dp
      kad(1:${a}$*${d}$) = 0.0_dp
      kac(1:${a}$*${c}$) = 0.0_dp
      p_index = 0
      DO md = 1, ${d}$
         DO mc = 1, ${c}$
            DO mb = 1, ${b}$
               ks_bd = 0.0_dp
               ks_bc = 0.0_dp
               p_bd = pbd((md - 1)*${b}$+mb)
               p_bc = pbc((mc - 1)*${b}$+mb)
               DO ma = 1, ${a}$
                  p_index = p_index + 1
                  tmp = scale*prim(p_index)
                  ks_bc = ks_bc + tmp*pad((md - 1)*${a}$+ma)
                  ks_bd = ks_bd + tmp*pac((mc - 1)*${a}$+ma)
                  kad((md - 1)*${a}$+ma) = kad((md - 1)*${a}$+ma) - tmp*p_bc
                  kac((mc - 1)*${a}$+ma) = kac((mc - 1)*${a}$+ma) - tmp*p_bd
               END DO
               kbd((md - 1)*${b}$+mb) = kbd((md - 1)*${b}$+mb) - ks_bd
               kbc((mc - 1)*${b}$+mb) = kbc((mc - 1)*${b}$+mb) - ks_bc
            END DO
         END DO
      END DO
      END SUBROUTINE ${name}$
      #:enddef

! **************************************************************************************************
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
! **************************************************************************************************
      SUBROUTINE contract_block(ma_max, mb_max, mc_max, md_max, kbd, kbc, kad, kac, pbd, pbc, pad, pac, prim, scale)
         INTEGER                                  :: ma_max, mb_max, mc_max, md_max
         REAL(KIND=dp) :: kbd(mb_max*md_max), kbc(mb_max*mc_max), &
                          kad(ma_max*md_max), kac(ma_max*mc_max), pbd(mb_max*md_max), &
                          pbc(mb_max*mc_max), pad(ma_max*md_max), pac(ma_max*mc_max), &
                          prim(ma_max*mb_max*mc_max*md_max), scale

#if !defined (__LIBINT)
         MARK_USED(ma_max)
         MARK_USED(mb_max)
         MARK_USED(mc_max)
         MARK_USED(md_max)
         MARK_USED(kbd)
         MARK_USED(kbc)
         MARK_USED(kad)
         MARK_USED(kac)
         MARK_USED(pbd)
         MARK_USED(pbc)
         MARK_USED(pad)
         MARK_USED(pac)
         MARK_USED(prim)
         MARK_USED(scale)
         CPABORT("libint not compiled in")
#else
         SELECT CASE (ma_max)
            #:for a in block_sizes
               CASE (${a}$)
               SELECT CASE (mb_max)
                  #:for b in block_sizes
                     CASE (${b}$)
                     #:if a*b < max_product
                        SELECT CASE (mc_max)
                           #:for c in block_sizes
                              CASE (${c}$)
                              #:if a*b*c < max_product
                                 SELECT CASE (md_max)
                                    #:for d in block_sizes
                                       CASE (${d}$)
                                       #:if a*b*c*d < max_product
                                          $:call_block((a, b, c, d))
                                       #:else
                                          $:call_block((a, b, c))
                                       #:endif
                                    #:endfor
                                 CASE DEFAULT
                                    $:call_block()
                                 END SELECT
                              #:else
                                 $:call_block((a, b))
                              #:endif
                           #:endfor
                        CASE DEFAULT
                           $:call_block()
                        END SELECT
                     #:else
                        $:call_block((a, b))
                     #:endif
                  #:endfor
               CASE DEFAULT
                  $:call_block()
               END SELECT
            #:endfor
         CASE DEFAULT
            $:call_block()
         END SELECT
#endif
      END SUBROUTINE contract_block

#if defined(__LIBINT)
      $:block()

      #:for a in block_sizes
         #:for b in block_sizes
            #:if a*b < max_product
               #:for c in block_sizes
                  #:if a*b*c < max_product
                     #:for d in block_sizes
                        #:if a*b*c*d < max_product
                           $:block((a, b, c, d))
                        #:endif
                     #:endfor
                     #:if a*b*c*block_sizes[-1] >= max_product
                        $:block((a, b, c))
                     #:endif
                  #:endif
               #:endfor
            #:endif
            #:if a*b*block_sizes[-1] >= max_product
               $:block((a, b))
            #:endif
         #:endfor
      #:endfor
#endif
   END MODULE hfx_contract_block
