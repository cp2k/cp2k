/*-----------------------------------------------------------------------------!
 *   CP2K: A general program to perform molecular dynamics simulations         !
 *   Copyright (C) 2002,2003,2004  CP2K developers group                       !
 *-----------------------------------------------------------------------------!
 *
 * c bridge of the f77 interface of cp2k
 *
 * 11.2004 created [fawzi]
 * 
 * NOTES mpi bridging not yet done
 */


/* compiler dependent declarations, now hardcoded for nag-gcc on mac */
typedef double f_real;
typedef int f_integer;
/* end of compiler dependent declarations (.c files contains some more) */

#if defined(__cplusplus) && ! defined(cp2k_direct_connection)
extern "C" {
#endif
  int ccp_create_fenv(int *env_id,const char * input_path, const char * output_path);
  int ccp_run_input(const char *input_path, const char *output_path);
  int ccp_init_cp2k(int init_mpi);
  int ccp_finalize_cp2k(int init_mpi);
  int ccp_destroy_fenv(int env_id);
  int ccp_set_pos(int env_id, f_real *new_pos, int n_el);
  int ccp_get_pos(int env_id, f_real *pos, int n_el);
  int ccp_get_natom(int env_id, int *natom);
  int ccp_get_force(int env_id, f_real *force, int n_el);
  int ccp_get_energy(int env_id, f_real *e_pot);
  int ccp_calc_energy_force(int env_id, int calc_force);
  int ccp_calc_force(int env_id, f_real *new_pos, int n_el_pos, 
		     f_real *e_pot, f_real *force, int n_el_force);
  int ccp_calc_energy(int env_id, f_real *new_pos, int n_el, f_real *e_pot);
  /* mpi tranfer not yet done */ 
#if defined(__cplusplus) && ! defined(cp2k_direct_connection)
}
#endif
