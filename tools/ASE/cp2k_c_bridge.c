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

#include "cp2k_c_bridge.h"
#include <string.h>

/* compiler dependent declarations, now hardcoded for nag-gcc on mac */
#define F77_NAME(name) name ## _
#define TRAIL(var) ,int _ ## var ## _len
#define DECLARE_CHARACTER_DYN(var) int _ ## var ## _len
#define F77_CREATE_CHARACTER_DYN(var,length) _ ## var ## _len = length
#define CHARACTER_ARG(var) var
#define TRAIL_ARG(var) ,_ ## var ## _len
#define F77_FREE_CHARACTER_DYN(var)
/* end of compiler dependent declarations */

#if defined(__cplusplus)
extern "C" {
#endif
  void F77_NAME(cp_create_fenv)(f_integer *env_id, const char *input_path,
			     const char *ouput_path,f_integer *ierr
			     TRAIL(input_path) TRAIL(ouput_path));
  void F77_NAME(cp_run_input)(const char *input_path,
			   const char *ouput_path,f_integer *ierr
			   TRAIL(input_path) TRAIL(ouput_path));
  void F77_NAME(cp_init_cp2k)(f_integer *init_mpi, f_integer *ierr);
  void F77_NAME(cp_finalize_cp2k)(f_integer *finalize_mpi, f_integer *ierr);
  void F77_NAME(cp_destroy_fenv)(f_integer *env_id, f_integer *ierr);
  void F77_NAME(cp_set_pos)(f_integer *env_id, f_real *pos,
			 f_integer *n_el,f_integer *ierr);
  void F77_NAME(cp_get_pos)(f_integer *env_id, f_real *pos,
			 f_integer *n_el,f_integer *ierr);
  void F77_NAME(cp_get_natom)(f_integer *env_id, f_integer *natom,
			      f_integer *ierr);
  void F77_NAME(cp_get_force)(f_integer *env_id, f_real *force,
			 f_integer *n_el,f_integer *ierr);
  void F77_NAME(cp_get_energy)(f_integer *env_id, f_real *e_pot,
			    f_integer *ierr);
  void F77_NAME(cp_calc_energy_force)(f_integer *env_id, 
				   f_integer * calc_force, f_integer *ierr);
  void F77_NAME(cp_calc_force)(f_integer *env_id, f_real *pos, 
			    f_integer *n_el_pos,f_real *e_pot,
			    f_real *force, f_integer *n_el_force,
			    f_integer *ierr);
  void F77_NAME(cp_calc_energy)(f_integer *env_id, f_real *pos, 
				f_integer *n_el_pos,f_real *e_pot,
				f_integer *ierr);
#if defined(__cplusplus)
}
#endif

/****** character len* functions *******/
int ccp_create_fenv(int *env_id,const char * input_path, const char * output_path)
{
  f_integer e_id,ierr;
  DECLARE_CHARACTER_DYN(input_path);
  DECLARE_CHARACTER_DYN(output_path);

  e_id=(f_integer)env_id;
  F77_CREATE_CHARACTER_DYN(input_path,strlen(input_path));
  F77_CREATE_CHARACTER_DYN(output_path,strlen(output_path));

  F77_NAME(cp_create_fenv)(&e_id,CHARACTER_ARG(input_path),
                           CHARACTER_ARG(output_path), 
                           &ierr TRAIL_ARG(input_path)
		           TRAIL_ARG(output_path));

  F77_FREE_CHARACTER_DYN(input_path);
  F77_FREE_CHARACTER_DYN(output_path);
  *env_id=(int)e_id;
  return (int)ierr;
}

int ccp_run_input(const char *input_path, const char *output_path)
{
  f_integer ierr;
  DECLARE_CHARACTER_DYN(input_path);
  DECLARE_CHARACTER_DYN(output_path);

  F77_CREATE_CHARACTER_DYN(input_path,strlen(input_path));
  F77_CREATE_CHARACTER_DYN(output_path,strlen(output_path));

  F77_NAME(cp_run_input)(CHARACTER_ARG(input_path),
                         CHARACTER_ARG(output_path), 
                         &ierr TRAIL_ARG(input_path)
		         TRAIL_ARG(output_path));

  F77_FREE_CHARACTER_DYN(input_path);
  F77_FREE_CHARACTER_DYN(output_path);
  return (int)ierr;
}

/****** end of character len* functions *******/

int ccp_init_cp2k(int init_mpi)
{
  f_integer i_mpi,ierr;
  i_mpi=init_mpi;
  
  F77_NAME(cp_init_cp2k)(&i_mpi,&ierr);
  return (int)ierr;
}

int ccp_finalize_cp2k(int init_mpi)
{
  f_integer i_mpi,ierr;
  i_mpi=(f_integer)init_mpi;
  
  F77_NAME(cp_finalize_cp2k)(&i_mpi,&ierr);
  return (int)ierr;
}

int ccp_destroy_fenv(int env_id)
{
  f_integer e_id,ierr;
  
  e_id=(f_integer)env_id;
  F77_NAME(cp_destroy_fenv)(&e_id,&ierr);
  return (int)ierr;
}

int ccp_set_pos(int env_id, f_real *new_pos,int n_el)
{
  f_integer e_id, ierr, nel;
  
  e_id=(f_integer)env_id;
  nel=(f_integer)n_el;
  F77_NAME(cp_set_pos)(&e_id,new_pos,&nel,&ierr);
  return (int)ierr;
}

int ccp_get_pos(int env_id, f_real *pos,int n_el)
{
  f_integer e_id, ierr, nel;
  
  e_id=(f_integer)env_id;
  nel=(f_integer)n_el;
  F77_NAME(cp_get_pos)(&e_id,pos,&nel,&ierr);
  return (int)ierr;
}

int ccp_get_natom(int env_id, int *natom)
{
  f_integer e_id, ierr, nat;
  
  e_id=(f_integer)env_id;
  F77_NAME(cp_get_natom)(&e_id,&nat,&ierr);
  *natom=(int)nat;
  return (int)ierr;
}

int ccp_get_force(int env_id, f_real *force,int n_el)
{
  f_integer e_id, ierr, nel;
  
  e_id=(f_integer)env_id;
  nel=(f_integer)n_el;
  F77_NAME(cp_get_force)(&e_id,force,&nel,&ierr);
  return (int)ierr;
}

int ccp_get_energy(int env_id, f_real *e_pot)
{
  f_integer e_id, ierr;
  
  e_id=(f_integer)env_id;
  F77_NAME(cp_get_energy)(&e_id,e_pot,&ierr);
  return (int)ierr;
}

int ccp_calc_energy_force(int env_id, int calc_force)
{
  f_integer e_id, ierr, calc_f;
  
  e_id=(f_integer)env_id;
  calc_f=(f_integer)(calc_f==0?0:1);
  F77_NAME(cp_calc_energy_force)(&e_id,&calc_f,&ierr);
  return (int)ierr;
}

int ccp_calc_force(int env_id, f_real *new_pos, int n_el_pos, 
		   f_real *e_pot, f_real *force, int n_el_force)
{
  f_integer e_id, nel_pos,nel_force,ierr;
  
  e_id=(f_integer)env_id;
  nel_pos=(f_integer)n_el_pos;
  nel_force=(f_integer)n_el_force;
  F77_NAME(cp_calc_force)(&e_id,new_pos,&nel_pos,e_pot,force,&nel_force,&ierr);
  return (int)ierr;
}

int ccp_calc_energy(int env_id, f_real *new_pos, int n_el,f_real *e_pot)
{
  f_integer e_id,nel,ierr;
  
  e_id=(f_integer)env_id;
  nel=(f_integer)n_el;
  F77_NAME(cp_calc_energy)(&e_id,new_pos,&nel,e_pot,&ierr);
  return (int)ierr;
}
