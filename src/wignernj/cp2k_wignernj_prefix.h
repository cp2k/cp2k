/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: BSD-3-Clause                                     */
/*----------------------------------------------------------------------------*/

/* Namespace the bundled libwignernj.
 *
 * libwignernj exports its C API under short, unqualified names (gaunt,
 * clebsch_gordan, wigner3j, xmalloc, ...). That is fine for a shared library,
 * but CP2K compiles this bundled copy straight into libcp2k, next to every
 * other library CP2K links, where such names collide: libpace defines its own
 * clebsch_gordan, and xmalloc is a name half the C world uses. Rename every
 * symbol the bundled objects export, so that the bundled copy contributes
 * nothing to the global namespace but cp2k_wignernj_*.
 *
 * Every bundled .c file includes this header before anything else, so the
 * renaming covers both the definitions and the declarations in wignernj.h.
 * src/common/wignernj_interface.F binds to the prefixed names when
 * __LIBWIGNERNJ_BUNDLED is set and to the plain ones when CP2K links against
 * an installed libwignernj.
 *
 * When resynchronising with a newer libwignernj, regenerate this list from
 *   nm -g --defined-only *.o
 * over the bundled objects built without this header.
 */

#ifndef CP2K_WIGNERNJ_PREFIX_H
#define CP2K_WIGNERNJ_PREFIX_H

#define bigint_add cp2k_wignernj_bigint_add
#define bigint_bit_length cp2k_wignernj_bigint_bit_length
#define bigint_cmp cp2k_wignernj_bigint_cmp
#define bigint_copy cp2k_wignernj_bigint_copy
#define bigint_div_prime_pow cp2k_wignernj_bigint_div_prime_pow
#define bigint_div_u64 cp2k_wignernj_bigint_div_u64
#define bigint_div_u64_exact cp2k_wignernj_bigint_div_u64_exact
#define bigint_free cp2k_wignernj_bigint_free
#define bigint_frexp cp2k_wignernj_bigint_frexp
#define bigint_frexpf cp2k_wignernj_bigint_frexpf
#define bigint_frexpl cp2k_wignernj_bigint_frexpl
#define bigint_init cp2k_wignernj_bigint_init
#define bigint_is_zero cp2k_wignernj_bigint_is_zero
#define bigint_mul cp2k_wignernj_bigint_mul
#define bigint_mul_prime_pow cp2k_wignernj_bigint_mul_prime_pow
#define bigint_mul_prime_pow_ws cp2k_wignernj_bigint_mul_prime_pow_ws
#define bigint_mul_u64 cp2k_wignernj_bigint_mul_u64
#define bigint_mul_ws cp2k_wignernj_bigint_mul_ws
#define bigint_reserve cp2k_wignernj_bigint_reserve
#define bigint_set_u64 cp2k_wignernj_bigint_set_u64
#define bigint_set_zero cp2k_wignernj_bigint_set_zero
#define bigint_sub_signed cp2k_wignernj_bigint_sub_signed
#define bigint_to_double cp2k_wignernj_bigint_to_double
#define bigint_to_float cp2k_wignernj_bigint_to_float
#define bigint_to_long_double cp2k_wignernj_bigint_to_long_double
#define bigint_words_for_factorial cp2k_wignernj_bigint_words_for_factorial
#define bigint_ws_free cp2k_wignernj_bigint_ws_free
#define bigint_ws_init cp2k_wignernj_bigint_ws_init
#define bigint_ws_reserve cp2k_wignernj_bigint_ws_reserve
#define clebsch_gordan cp2k_wignernj_clebsch_gordan
#define clebsch_gordan_f cp2k_wignernj_clebsch_gordan_f
#define clebsch_gordan_l cp2k_wignernj_clebsch_gordan_l
#define clebsch_gordan_max_factorial cp2k_wignernj_clebsch_gordan_max_factorial
#define fano_x cp2k_wignernj_fano_x
#define fano_x_f cp2k_wignernj_fano_x_f
#define fano_x_l cp2k_wignernj_fano_x_l
#define fano_x_max_factorial cp2k_wignernj_fano_x_max_factorial
#define g_nprimes cp2k_wignernj_g_nprimes
#define g_pi_table cp2k_wignernj_g_pi_table
#define g_primes cp2k_wignernj_g_primes
#define gaunt cp2k_wignernj_gaunt
#define gaunt_f cp2k_wignernj_gaunt_f
#define gaunt_l cp2k_wignernj_gaunt_l
#define gaunt_max_factorial cp2k_wignernj_gaunt_max_factorial
#define gaunt_real cp2k_wignernj_gaunt_real
#define gaunt_real_f cp2k_wignernj_gaunt_real_f
#define gaunt_real_l cp2k_wignernj_gaunt_real_l
#define gaunt_real_max_factorial cp2k_wignernj_gaunt_real_max_factorial
#define legendre_valuation cp2k_wignernj_legendre_valuation
#define pfrac_bigint_mul_prime_pow_array                                       \
  cp2k_wignernj_pfrac_bigint_mul_prime_pow_array
#define pfrac_copy cp2k_wignernj_pfrac_copy
#define pfrac_div_factorial cp2k_wignernj_pfrac_div_factorial
#define pfrac_free cp2k_wignernj_pfrac_free
#define pfrac_init cp2k_wignernj_pfrac_init
#define pfrac_lcm_scaled_product cp2k_wignernj_pfrac_lcm_scaled_product
#define pfrac_mul_factorial cp2k_wignernj_pfrac_mul_factorial
#define pfrac_mul_int cp2k_wignernj_pfrac_mul_int
#define pfrac_mul_pow_into_acc cp2k_wignernj_pfrac_mul_pow_into_acc
#define pfrac_to_sqrt_rational cp2k_wignernj_pfrac_to_sqrt_rational
#define pfrac_to_sqrt_rational_ws cp2k_wignernj_pfrac_to_sqrt_rational_ws
#define pfrac_zero cp2k_wignernj_pfrac_zero
#define racah_w cp2k_wignernj_racah_w
#define racah_w_f cp2k_wignernj_racah_w_f
#define racah_w_l cp2k_wignernj_racah_w_l
#define racah_w_max_factorial cp2k_wignernj_racah_w_max_factorial
#define wigner3j cp2k_wignernj_wigner3j
#define wigner3j_exact cp2k_wignernj_wigner3j_exact
#define wigner3j_f cp2k_wignernj_wigner3j_f
#define wigner3j_l cp2k_wignernj_wigner3j_l
#define wigner3j_max_factorial cp2k_wignernj_wigner3j_max_factorial
#define wigner6j cp2k_wignernj_wigner6j
#define wigner6j_exact cp2k_wignernj_wigner6j_exact
#define wigner6j_f cp2k_wignernj_wigner6j_f
#define wigner6j_l cp2k_wignernj_wigner6j_l
#define wigner6j_max_factorial cp2k_wignernj_wigner6j_max_factorial
#define wigner9j cp2k_wignernj_wigner9j
#define wigner9j_exact cp2k_wignernj_wigner9j_exact
#define wigner9j_f cp2k_wignernj_wigner9j_f
#define wigner9j_l cp2k_wignernj_wigner9j_l
#define wigner9j_max_factorial cp2k_wignernj_wigner9j_max_factorial
#define wignernj_exact_free cp2k_wignernj_wignernj_exact_free
#define wignernj_exact_init cp2k_wignernj_wignernj_exact_init
#define wignernj_exact_mul_sqrt_int cp2k_wignernj_wignernj_exact_mul_sqrt_int
#define wignernj_exact_reset cp2k_wignernj_wignernj_exact_reset
#define wignernj_exact_to_double cp2k_wignernj_wignernj_exact_to_double
#define wignernj_exact_to_float cp2k_wignernj_wignernj_exact_to_float
#define wignernj_exact_to_long_double                                          \
  cp2k_wignernj_wignernj_exact_to_long_double
#define wignernj_factorial_cache_fill                                          \
  cp2k_wignernj_wignernj_factorial_cache_fill
#define wignernj_factorial_cache_release                                       \
  cp2k_wignernj_wignernj_factorial_cache_release
#define wignernj_max_factorial_arg cp2k_wignernj_wignernj_max_factorial_arg
#define wignernj_real_ylm_in_complex_ylm                                       \
  cp2k_wignernj_wignernj_real_ylm_in_complex_ylm
#define wignernj_real_ylm_in_complex_ylm_f                                     \
  cp2k_wignernj_wignernj_real_ylm_in_complex_ylm_f
#define wignernj_real_ylm_in_complex_ylm_l                                     \
  cp2k_wignernj_wignernj_real_ylm_in_complex_ylm_l
#define wignernj_scratch_acquire cp2k_wignernj_wignernj_scratch_acquire
#define wignernj_scratch_lcm_clear cp2k_wignernj_wignernj_scratch_lcm_clear
#define wignernj_scratch_lcm_dirty cp2k_wignernj_wignernj_scratch_lcm_dirty
#define wignernj_scratch_release cp2k_wignernj_wignernj_scratch_release
#define wignernj_scratch_relinquish cp2k_wignernj_wignernj_scratch_relinquish
#define wignernj_scratch_terms_reserve                                         \
  cp2k_wignernj_wignernj_scratch_terms_reserve
#define wignernj_thread_cleanup cp2k_wignernj_wignernj_thread_cleanup
#define wignernj_thread_local_scratch_available                                \
  cp2k_wignernj_wignernj_thread_local_scratch_available
#define wignernj_warmup_to cp2k_wignernj_wignernj_warmup_to
#define xalloc_set_test_failure_countdown                                      \
  cp2k_wignernj_xalloc_set_test_failure_countdown
#define xcalloc cp2k_wignernj_xcalloc
#define xmalloc cp2k_wignernj_xmalloc
#define xrealloc cp2k_wignernj_xrealloc

#endif /* CP2K_WIGNERNJ_PREFIX_H */
