/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2016  CP2K developers group                         *
 *****************************************************************************/

#include <stdbool.h>

#ifndef CP2K_C_INTERFACE_H
#define CP2K_C_INTERFACE_H

#ifdef __cplusplus
 extern "C" {
#endif

/* transport parameters read from a CP2K input file. */
/* DO NOT change the ORDERING or the NAMES           */
    typedef struct{
           int     n_occ;
           int     n_atoms;
           double  energy_diff;
           double  evoltfactor;
           int     method;
           int     injection_method;
           int     rlaxis_integration_method;
           int     linear_solver;
           int     n_abscissae;
           int     n_kpoint;
           int     num_interval;
           int     num_contacts;
           int     tasks_per_point;
           int     gpus_per_point;
           int     n_blocks;
           int     n_points_beyn;
           int     ncrc_beyn;
           int     cutout[2];
           double  colzero_threshold;
           double  eps_limit;
           double  eps_limit_cc;
           double  eps_decay;
           double  eps_singularity_curvatures;
           double  eps_mu;
           double  eps_eigval_degen;
           double  eps_fermi;
           double  energy_interval;
           double  min_interval;
           double  temperature;
           double  n_rand_beyn;
           double  n_rand_cc_beyn;
           double  svd_cutoff;
           int*    contacts_data;
           int*    nsgf;
           double* zeff;
           int*    tridiag_blocks;
           bool    extra_scf;
           } cp2k_transport_parameters;

/* CP2K's C-interoperable CSR matrix                 */
/* DO NOT change the ORDERING or the NAMES           */
   typedef struct{
           int     nrows_total;
           int     ncols_total;
           int     nze_total;
           int     nze_local;
           int     nrows_local;
           int     data_type;
           int     first_row;
           int*    rowptr_local;
           int*    colind_local;
           int*    nzerow_local;
           double* nzvals_local;
           } cp2k_csr_interop_type;

/* funtion pointer type pointing to a C routine that takes the S and H matrices as input and outputs a P matrix */
   typedef void(*ptr2function) (cp2k_transport_parameters, cp2k_csr_interop_type, cp2k_csr_interop_type, cp2k_csr_interop_type *, cp2k_csr_interop_type *);
   void c_scf_method(cp2k_transport_parameters cp2k_transport_params, cp2k_csr_interop_type S, cp2k_csr_interop_type KS, cp2k_csr_interop_type * P, cp2k_csr_interop_type * PImag);

/* routines for initializing and finalizing a CP2K run and setting up and destroying CP2K environments */
   void cp_c_init_cp2k(int *init_mpi, int *ierr);
   void cp_c_get_version(char* version_str, int str_length);
   void cp_c_finalize_cp2k(int *finalize_mpi, int *ierr);
   void cp_c_create_fenv(int *new_env_id, char *input_file_path, char *output_file_path, int *ierr);
   void cp_c_create_fenv_comm(int *new_env_id, char *input_file_path, char *output_file_path, int *mpi_comm, int *ierr);
   void cp_c_destroy_fenv(int *env_id, int *ierr);

   void cp_c_get_natom(int *env_id, int *natom, int *ierr);
   void cp_c_get_nparticle(int *env_id, int *nparticle, int *ierr);
   void cp_c_get_energy(int *env_id, double *e_pot, int *ierr);
   void cp_c_get_force(int *env_id, double force[], int *n_el, int *ierr);
   void cp_c_get_pos(int *env_id, double pos[], int *n_el, int *ierr);

   void cp_c_calc_energy_force(int *env_id, int *calc_force, int *ierr);

   void cp_c_run_input(char *input_file_path, char *output_file_path, int *ierr);
   void cp_c_run_input_comm(char *input_file_path, char *output_file_path, int *mpi_comm, int *ierr);

/* routine that gets a C function pointer and passes it to CP2K */
   void cp_c_ext_method_set_ptr(int *f_env_id, ptr2function, int *ierr);

#ifdef __cplusplus
 }
#endif

#endif
