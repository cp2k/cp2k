/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2016  CP2K developers group                         *
 *****************************************************************************/

#include <stdbool.h>

/** \file cp2k_c_interface.h
 * \brief Definitions for the functions exported in cp2k_c_interface.F
 * \sa cp2k_c_interface.F
 */

#ifndef CP2K_C_INTERFACE_H
#define CP2K_C_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

/** \brief Get the CP2K version string
 * \param version_str The buffer to write the version string into
 * \param str_length The size of the buffer (must be large enough)
 */
void cp_c_get_version(char* version_str, int str_length);

/** \brief Initialize CP2K
 * \param init_mpi Set to 1 if CP2K should initialize MPI.
 *                 If set to 0, you have to provide the MPI_COMM when creating environments.
 * \param ierr     Non-zero if an error occurred
 * \sa cp_init_cp2k
 * \warning You are supposed to call cp_c_finalize_cp2k() before terminating the program.
 */
void cp_c_init_cp2k(int* init_mpi, int* ierr);

/** \brief Finalize CP2K
 * \param init_mpi Set to 1 if CP2K should also finalize MPI.
 *                 If set to 0, your code has to call MPI_Finalize (or equivalent).
 * \param ierr Non-zero if an error occurred
 * \sa cp_finalize_cp2k
 */
void cp_c_finalize_cp2k(int* finalize_mpi, int* ierr);

/** \brief Create a new force environment
 * \param new_env_id Tag for the created force environment
 * \param input_file_path Path to a CP2K input file
 * \param output_file_path Path to a file where CP2K is going to append its output (created if non-existent)
 * \param ierr Non-zero if an error occurred
 * \sa cp_create_fenv
 * \warning You are supposed to call cp_c_destroy_fenv() to cleanup, before cp_c_finalize_cp2k().
 */
void cp_c_create_fenv(int* new_env_id, char* input_file_path, char* output_file_path, int* ierr);

/** \brief Create a new force environment (custom managed MPI)
 * \param new_env_id Tag for the created force environment
 * \param input_file_path Path to a CP2K input file
 * \param output_file_path Path to a file where CP2K is going to write its output.
 *                         Will be created if not existent, otherwise appended.
 * \param mpi_comm MPI_COMM if MPI is not managed by CP2K
 * \param ierr Non-zero if an error occurred
 * \sa cp_create_fenv_comm
 * \warning You are supposed to call cp_c_destroy_fenv() to cleanup, before cp_c_finalize_cp2k().
 */
void cp_c_create_fenv_comm(int * new_env_id, char* input_file_path, char* output_file_path, int* mpi_comm, int* ierr);

/** \brief Destroy/cleanup a force environment
 * \param env_id Tag of the force environment as returned by either cp_c_create_fenv() or cp_c_create_fenv_comm()
 * \param ierr Non-zero if an error occurred
 * \sa cp_destroy_fenv
 */
void cp_c_destroy_fenv(int* env_id, int* ierr);

/** \brief Set positions of the particles
 * \param env_id Tag of the force environment as returned by either cp_c_create_fenv() or cp_c_create_fenv_comm()
 * \param new_pos Array containing the new positions of the particles
 * \param n_el Size of the new_pos array
 * \param ierr Non-zero if an error occurred
 * \sa cp_set_pos
 */
void cp_c_set_pos(int* env_id, double* new_pos, int* n_el, int* ierr);

/** \brief Set velocity of the particles
 * \param env_id Tag of the force environment as returned by either cp_c_create_fenv() or cp_c_create_fenv_comm()
 * \param new_vel Array containing the new velocities of the particles
 * \param n_el Size of the new_vel array
 * \param ierr Non-zero if an error occurred
 * \sa cp_set_vel
 */
void cp_c_set_vel(int* env_id, double* new_vel, int* n_el, int* ierr);

/** \brief Get an arbitrary result as 1D array from CP2K
 * \param env_id Tag of the force environment as returned by either cp_c_create_fenv() or cp_c_create_fenv_comm()
 * \param description The string tag of the result
 * \param N Size of the RESULT array
 * \param RESULT Pre-allocated array
 * \param ierr Non-zero if an error occurred
 * \sa cp_get_result_r1
 */
void cp_c_get_result_r1(int* env_id, char* description, int* N, double* RESULT, int* ierr);

/** \brief Get the number of atoms
 * \param env_id Tag of the force environment as returned by either cp_c_create_fenv() or cp_c_create_fenv_comm()
 * \param natom The number of atoms
 * \param ierr Non-zero if an error occurred
 * \sa cp_get_natom
 */
void cp_c_get_natom(int* env_id, int* natom, int* ierr);

/** \brief Get the number of particles
 * \param env_id Tag of the force environment as returned by either cp_c_create_fenv() or cp_c_create_fenv_comm()
 * \param nparticle The number of particles
 * \param ierr Non-zero if an error occurred
 * \sa cp_get_nparticle
 */
void cp_c_get_nparticle(int* env_id, int* nparticle, int* ierr);

/** \brief Get the positions of the particles
 * \param env_id Tag of the force environment as returned by either cp_c_create_fenv() or cp_c_create_fenv_comm()
 * \param pos Pre-allocated array of at least 3*nparticle elements. Use cp_c_get_nparticle() to get the number of particles.
 * \param n_el Size of the force array
 * \param ierr Non-zero if an error occurred
 * \sa cp_get_pos
 */
void cp_c_get_pos(int* env_id, double* pos, int* n_el, int* ierr);

/** \brief Get the forces for the particles
 * \param env_id Tag of the force environment as returned by either cp_c_create_fenv() or cp_c_create_fenv_comm()
 * \param force Pre-allocated array of at least 3*nparticle elements. Use cp_c_get_nparticle() to get the number of particles.
 * \param n_el Size of the force array
 * \param ierr Non-zero if an error occurred
 * \sa cp_get_force
 */
void cp_c_get_force(int* env_id, double* force, int* n_el, int* ierr);

/** \brief Get the potential energy of the system
 * \param env_id Tag of the force environment as returned by either cp_c_create_fenv() or cp_c_create_fenv_comm()
 * \param e_pot The potential energy
 * \param ierr Non-zero if an error occurred
 * \sa cp_get_energy
 */
void cp_c_get_energy(int* env_id, double* e_pot, int* ierr);

/** \brief Calculate energy and forces (optional) of the system
 * \param env_id Tag of the force environment as returned by either cp_c_create_fenv() or cp_c_create_fenv_comm()
 * \param calc_force Set to a non-zero value if forces should be calculated (otherwise only energies)
 * \param ierr Non-zero if an error occurred
 * \sa cp_calc_energy_force
 */
void cp_c_calc_energy_force(int* env_id, int* calc_force, int* ierr);

/** \brief Calculate and return energy of the configuration given by the specified positions
 * \param env_id Tag of the force environment as returned by either cp_c_create_fenv() or cp_c_create_fenv_comm()
 * \param pos Array containing the positions of the particles
 * \param n_el Size of the pos array
 * \param e_pot The potential energy
 * \param ierr Non-zero if an error occurred
 * \sa cp_calc_energy, cp_c_get_nparticle
 */
void cp_c_calc_energy(int* env_id, double* pos, int* n_el, double* e_pot, int* ierr);

/** \brief Calculate and return energy and forces of the configuration given by the provided positions
 * \param env_id Tag of the force environment as returned by either cp_c_create_fenv() or cp_c_create_fenv_comm()
 * \param pos Array containing the positions of the particles
 * \param n_el_pos Size of the pos array
 * \param e_pot The potential energy
 * \param force Pre-allocated array of at least 3*nparticle elements
 * \param n_el_force Size of the force array
 * \param ierr Non-zero if an error occurred
 * \sa cp_calc_energy, cp_c_get_nparticle, cp_c_calc_energy
 */
void cp_c_calc_force(int* env_id, double* pos, int* n_el_pos, double* e_pot, double* force, int* n_el_force, int* ierr);

/** \brief Make a CP2K run with the given input file
 * \param input_file_path Path to a CP2K input file
 * \param output_file_path Path to a file where CP2K is going to append its output (created if non-existent)
 * \param ierr Non-zero if an error occurred
 * \sa cp_run_input
 */
void cp_c_run_input(char* input_file_path, char* output_file_path, int* ierr);

/** \brief Make a CP2K run with the given input file (custom managed MPI)
 * \param input_file_path Path to a CP2K input file
 * \param output_file_path Path to a file where CP2K is going to append its output (created if non-existent)
 * \param mpi_comm MPI_COMM if MPI is not managed by CP2K
 * \param ierr Non-zero if an error occurred
 * \sa cp_run_input_comm
 */
void cp_c_run_input_comm(char* input_file_path, char* output_file_path, int* mpi_comm, int* ierr);

/** \brief Perform the shake procedure (enforce constraints)
 * \param env_id Tag of the force environment as returned by either cp_c_create_fenv() or cp_c_create_fenv_comm()
 * \param dt Timestep
 * \param shake_tol Tolerance
 * \param ierr Non-zero if an error occurred
 * \sa cp_do_shake
 */
void cp_c_do_shake(int* f_env_id, double* dt, double* shake_tol, int* ierr);

/** \brief Transport parameters read from a CP2K input file
 *
 * This definition matches the respective type definition in the transport_env_types module
 */
typedef struct {
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

/** \brief CP2K's C-interoperable CSR matrix
 *
 * This definition matches the respective type definition in the transport_env_types module
 */
typedef struct {
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

/** \brief Function pointer type for the externally evaluated density matrix
 *
 * Function pointer type pointing to a C routine that takes the S and H matrices as input and outputs a P and PImag matrix.
 *
 * Function definition example:
 * \code{.c}
 * void c_scf_method(
 *     cp2k_transport_parameters cp2k_transport_params,
 *     cp2k_csr_interop_type S,
 *     cp2k_csr_interop_type KS,
 *     cp2k_csr_interop_type* P,
 *     cp2k_csr_interop_type* PImag
 *     );
 * \endcode
 * \sa cp2k_transport_parameters, cp2k_csr_interop_type
 */
typedef void (*ext_method_callback_f_ptr) (
    cp2k_transport_parameters, // Transport parameters
    cp2k_csr_interop_type,  // S-Matrix
    cp2k_csr_interop_type,  // H-Matrix
    cp2k_csr_interop_type*, // P-Matrix
    cp2k_csr_interop_type*  // PImag-Matrix
    );

/** \brief Set the function callback for the externally evaluated density matrix
 */
void cp_c_ext_method_set_ptr(int* f_env_id, ext_method_callback_f_ptr func, int* ierr);

#ifdef __cplusplus
}
#endif

#endif
/* vim: set ts=4 sw=4 tw=0 :*/
