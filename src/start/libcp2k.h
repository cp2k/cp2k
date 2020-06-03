/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2020  CP2K developers group                         *
 *****************************************************************************/

#include <stdbool.h>

/**
 * \brief Definitions for the functions exported in libcp2k.F
 * \author Mohammad Hossein Bani-Hashemian
 */

#ifndef LIBCP2K_H
#define LIBCP2K_H

#ifdef __cplusplus
extern "C" {
#endif

typedef int force_env_t;

/** \brief Get the CP2K version string
 * \param version_str The buffer to write the version string into
 * \param str_length The size of the buffer (must be large enough)
 */
void cp2k_get_version(char* version_str, int str_length);

/** \brief Initialize CP2K and MPI
 * \warning You are supposed to call cp2k_finalize() before terminating the program.
 */
void cp2k_init();

/** \brief Initialize CP2K without initializing MPI
 * \warning You are supposed to call cp2k_finalize() before terminating the program.
 */
void cp2k_init_without_mpi();

/** \brief Finalize CP2K and MPI
 */
void cp2k_finalize();

/** \brief Finalize CP2K and without finalizing MPI
 */
void cp2k_finalize_without_mpi();

/** \brief Create a new force environment
 * \param new_force_env the created force environment
 * \param input_file_path Path to a CP2K input file
 * \param output_file_path Path to a file where CP2K is going to append its output (created if non-existent)
 * \warning You are supposed to call cp2k_destroy_force_env() to cleanup, before cp2k_finalize().
 */
void cp2k_create_force_env(force_env_t* new_force_env, const char* input_file_path, const char* output_file_path);

/** \brief Create a new force environment (custom managed MPI)
 * \param new_force_env the created force environment
 * \param input_file_path Path to a CP2K input file
 * \param output_file_path Path to a file where CP2K is going to write its output.
 *                         Will be created if not existent, otherwise appended.
 * \param mpi_comm Fortran MPI communicator if MPI is not managed by CP2K
 * \warning You are supposed to call cp2k_destroy_force_env() to cleanup, before cp2k_finalize().
 */
void cp2k_create_force_env_comm(force_env_t* new_force_env, const char* input_file_path, const char* output_file_path, int mpi_comm);

/** \brief Destroy/cleanup a force environment
 * \param force_env the force environment
 */
void cp2k_destroy_force_env(force_env_t force_env);

/** \brief Set positions of the particles
 * \param force_env the force environment
 * \param new_pos Array containing the new positions of the particles
 * \param n_el Size of the new_pos array
 */
void cp2k_set_positions(force_env_t force_env, const double* new_pos, int n_el);

/** \brief Set velocity of the particles
 * \param force_env the force environment
 * \param new_vel Array containing the new velocities of the particles
 * \param n_el Size of the new_vel array
 */
void cp2k_set_velocities(force_env_t force_env, const double* new_vel, int n_el);

/** \brief Set the size of the cell
 * \param force_env the force environment
 * \param new_cell Array containing the new cell
 */
void cp2k_set_cell(force_env_t force_env, const double* new_cell);

/** \brief Get an arbitrary result as 1D array from CP2K
 * \param force_env the force environment
 * \param description The string tag of the result
 * \param results Pre-allocated array
 * \param n_el size of the results array
 */
void cp2k_get_result(force_env_t force_env, const char* description, double* result, int n_el);

/** \brief Get the number of atoms
 * \param force_env the force environment
 * \param natom The number of atoms
 */
void cp2k_get_natom(force_env_t force_env, int* natom);

/** \brief Get the number of particles
 * \param force_env the force environment
 * \param nparticle The number of particles
 */
void cp2k_get_nparticle(force_env_t force_env, int* nparticle);

/** \brief Get the positions of the particles
 * \param force_env the force environment
 * \param pos Pre-allocated array of at least 3*nparticle elements. Use cp2k_get_nparticle() to get the number of particles.
 * \param n_el Size of the force array
 */
void cp2k_get_positions(force_env_t force_env, double* pos, int n_el);

/** \brief Get the forces for the particles
 * \param force_env the force environment
 * \param force Pre-allocated array of at least 3*nparticle elements. Use cp2k_get_nparticle() to get the number of particles.
 * \param n_el Size of the force array
 */
void cp2k_get_forces(force_env_t force_env, double* force, int n_el);

/** \brief Get the potential energy of the system
 * \param force_env the force environment
 * \param e_pot The potential energy
 */
void cp2k_get_potential_energy(force_env_t force_env, double* e_pot);

/** \brief Calculate energy and forces of the system
 * \param force_env the force environment
 */
void cp2k_calc_energy_force(force_env_t force_env);

/** \brief Calculate only the energy of the system
 * \param force_env the force environment
 */
void cp2k_calc_energy(force_env_t force_env);

/** \brief Make a CP2K run with the given input file
 * \param input_file_path Path to a CP2K input file
 * \param output_file_path Path to a file where CP2K is going to append its output (created if non-existent)
 */
void cp2k_run_input(const char* input_file_path, const char* output_file_path);

/** \brief Make a CP2K run with the given input file (custom managed MPI)
 * \param input_file_path Path to a CP2K input file
 * \param output_file_path Path to a file where CP2K is going to append its output (created if non-existent)
 * \param mpi_comm Fortran MPI communicator if MPI is not managed by CP2K
 */
void cp2k_run_input_comm(const char* input_file_path, const char* output_file_path, int mpi_comm);

/** \brief Transport parameters read from a CP2K input file
 *
 * This definition matches the respective type definition in the transport_env_types module
 */
typedef struct {
    int     n_occ;
    int     n_atoms;
    double  energy_diff;
    double  evoltfactor;
    double  e_charge;
    double  boltzmann;
    double  h_bar;
    int     iscf;
    int     method;
    int     qt_formalism;
    int     injection_method;
    int     rlaxis_integration_method;
    int     linear_solver;
    int     matrixinv_method;
    int     transport_neutral;
    int     num_pole;
    int     ordering;
    int     row_ordering;
    int     verbosity;
    int     pexsi_np_symb_fact;
    int     n_kpoint;
    int     num_interval;
    int     num_contacts;
    int     stride_contacts;
    int     tasks_per_energy_point;
    int     tasks_per_pole;
    int     gpus_per_point;
    int     n_points_beyn;
    int     ncrc_beyn;
    int     tasks_per_integration_point;
    int     n_points_inv;
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
    double  dens_mixing;
    double  n_rand_beyn;
    double  n_rand_cc_beyn;
    double  svd_cutoff;
    int*    contacts_data;
    int*    nsgf;
    double* zeff;
    bool    obc_equilibrium;
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
 * Function pointer type pointing to a C routine that takes the S and H matrices as input and outputs a P matrix.
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
void cp2k_transport_set_callback(force_env_t force_env, ext_method_callback_f_ptr func);

/** \brief Get the number of molecular orbitals in the active space
 * \param force_env the force environment
 * \returns The number of elements or -1 if unavailable
 */
int cp2k_active_space_get_mo_count(force_env_t force_env);

/** \brief Get the Fock submatrix for the active space
 * \param force_env the force environment
 * \param buf Pre-allocated array of at least mo_count^2 elements.
 *            Use `cp2k_active_space_get_mo_count()` to get the number of molecular orbitals.
 * \param buf_len Size of the buf array
 * \returns The number of elements written to buf or -1 if unavailable
 */
long int cp2k_active_space_get_fock_sub(force_env_t force_env, double* buf, long int buf_len);

/** \brief Get the number of non-zero elements in the ERI matrix
 * \param force_env the force environment
 * \returns The number of elements or -1 if unavailable
 */
long int cp2k_active_space_get_eri_nze_count(force_env_t force_env);

/** \brief Get the non-zero elements of the ERI matrix
 *
 * The buf_coords will contain the coordinates in the format `[i1, j1, k1, l1, i2, j2, k2, l2, ... ]`.
 * \param force_env the force environment
 * \param buf_coords Pre-allocated array of at least 4*nze_count elements.
 *                   Use `cp2k_active_space_get_eri_nze_count()` to get the number of non-zero elements.
 * \param buf_coords_len Size of the buf_coords array
 * \param buf_values Pre-allocated array of at least nze_count elements.
 *                   Use `cp2k_active_space_get_eri_nze_count()` to get the number of non-zero elements.
 * \param buf_values_len Size of the buf_values array
 * \returns The number of elements written to buf_values or -1 if unavailable
 */
int cp2k_active_space_get_eri(force_env_t force_env,
                              int* buf_coords, long int buf_coords_len,
                              double* buf_values, long int buf_values_len);

#ifdef __cplusplus
}
#endif

#endif
/* vim: set ts=4 sw=4 tw=0 :*/
