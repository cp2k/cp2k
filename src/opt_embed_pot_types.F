!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE opt_embed_pot_types
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE kinds,                           ONLY: dp
   USE pw_types,                        ONLY: pw_p_type

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: opt_embed_pot_type

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'opt_embed_pot_types'

   TYPE opt_embed_pot_type
      TYPE(cp_fm_type), POINTER                 :: embed_pot_coef, embed_pot_grad, &
                                                   step, embed_pot_hess, &
                                                   prev_embed_pot_coef, prev_embed_pot_grad, &
                                                   prev_embed_pot_hess, kinetic_mat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)  :: w_func, max_diff, int_diff
      INTEGER                                   :: n_iter, i_iter, dimen_aux, last_accepted
      REAL(KIND=dp)                             :: lambda, allowed_decrease, max_trad, min_trad, &
                                                   grad_norm
      TYPE(pw_p_type), POINTER                  :: const_pot
      REAL(KIND=dp)                             :: reg_term, trust_rad, conv_max, conv_int, step_len
      LOGICAL                                   :: accept_step, newton_step, level_shift, &
                                                   add_const_pot, converged

   END TYPE opt_embed_pot_type

END MODULE opt_embed_pot_types

