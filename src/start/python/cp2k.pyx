# cython: language_level=2
#  vim: set ts=4 sw=4 tw=0 :

import numpy as np

from cpython.mem cimport PyMem_Malloc, PyMem_Free

cdef extern from "../libcp2k.h":
    ctypedef int force_env_t

    void cp2k_get_version(char* version_str, int str_length)
    void cp2k_init()
    void cp2k_init_without_mpi()
    void cp2k_finalize()
    void cp2k_finalize_without_mpi()
    void cp2k_create_force_env(force_env_t* new_force_env, const char* input_file_path, const char* output_file_path)
    void cp2k_create_force_env_comm(force_env_t* new_force_env, const char* input_file_path, const char* output_file_path, int mpi_comm)
    void cp2k_destroy_force_env(force_env_t force_env)
    void cp2k_set_positions(force_env_t force_env, const double* new_pos, int n_el)
    void cp2k_set_velocities(force_env_t force_env, const double* new_vel, int n_el)
    void cp2k_get_result(force_env_t force_env, const char* description, double* result, int n_el)
    void cp2k_get_natom(force_env_t force_env, int* natom)
    void cp2k_get_nparticle(force_env_t force_env, int* nparticle)
    void cp2k_get_positions(force_env_t force_env, double* pos, int n_el)
    void cp2k_get_forces(force_env_t force_env, double* force, int n_el)
    void cp2k_get_potential_energy(force_env_t force_env, double* e_pot)
    void cp2k_calc_energy_force(force_env_t force_env)
    void cp2k_calc_energy(force_env_t force_env)
    void cp2k_run_input(const char* input_file_path, const char* output_file_path)
    void cp2k_run_input_comm(const char* input_file_path, const char* output_file_path, int mpi_comm)

def get_version_string():
    n = 255 * sizeof(char)

    data = <char *>PyMem_Malloc(n)
    if not data:
        raise MemoryError()

    versionstr = ''
    try:
        cp2k_get_version(data, n)
        versionstr = data.decode('UTF-8')
    finally:
        PyMem_Free(data)

    return versionstr

def init(manage_mpi = True):
    if manage_mpi:
        cp2k_init()
    else:
        cp2k_init_without_mpi()

def finalize(manage_mpi = True):
    if manage_mpi:
        cp2k_finalize()
    else:
        cp2k_finalize_without_mpi()

def run_input(input_file_path, output_file_path = None, mpi_comm = None):
    input_file_path = input_file_path.encode('UTF-8')

    if output_file_path is None:
        output_file_path = u'__STD_OUT__'.encode('UTF-8')
    else:
        output_file_path = output_file_path.encode('UTF-8')

    if mpi_comm:
        cp2k_run_input_comm(input_file_path, output_file_path, mpi_comm)
    else:
        cp2k_run_input(input_file_path, output_file_path)

def create_force_env(input_file_path, output_file_path, mpi_comm = None):
    cdef force_env_t fenv
    if mpi_comm:
        cp2k_create_force_env_comm(&fenv, input_file_path, output_file_path, mpi_comm)
    else:
        cp2k_create_force_env(&fenv, input_file_path, output_file_path)

cdef class ForceEnvironment(object):
    cdef force_env_t _force_env

    def __init__(self, input_file_path not None, output_file_path = None, mpi_comm = None):
        input_file_path = input_file_path.encode('UTF-8')

        if output_file_path is None:
            output_file_path = u'__STD_OUT__'.encode('UTF-8')
        else:
            output_file_path = output_file_path.encode('UTF-8')

        if mpi_comm:
            cp2k_create_force_env_comm(&self._force_env, input_file_path, output_file_path, mpi_comm)
        else:
            cp2k_create_force_env(&self._force_env, input_file_path, output_file_path)

    def destroy(self):
        cp2k_destroy_force_env(self._force_env)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.destroy()

    
    property positions:

        def __get__(self):
            positions = np.zeros([3*self.nparticle], dtype=np.double)
            cdef double [::1] positions_view = positions
            cp2k_get_positions(self._force_env, &positions_view[0], positions_view.shape[0])
            return positions

        def __set__(self, double[::1] positions not None):
            if  positions.shape[0] != 3*self.nparticle:
                raise ValueError('the positions array must have exactly {} (3*nparticle) elements'.format(3*self.nparticle))

            cp2k_set_positions(self._force_env, &positions[0], positions.shape[0])

    def set_velocities(self, double[::1] velocities not None):
        if  velocities.shape[0] != 3*self.nparticle:
            raise ValueError('the velocities array must have exactly {} (3*nparticle) elements'.format(3*self.nparticle))

        cp2k_set_velocities(self._force_env, &velocities[0], velocities.shape[0])

    property natom:
        def __get__(self):
            cdef int natom
            cp2k_get_natom(self._force_env, &natom)
            return natom

    property nparticle:
        def __get__(self):
            cdef int nparticle
            cp2k_get_nparticle(self._force_env, &nparticle)
            return nparticle

    def get_result(self, description):
        results = np.zeros([3*self.nparticle], dtype=np.double)
        cdef double [::1] results_view = results
        cp2k_get_result(self._force_env, description.encode('UTF-8'), &results_view[0], results_view.shape[0])
        return results

    property forces:
        def __get__(self):
            forces = np.zeros([3*self.nparticle], dtype=np.double)
            cdef double [::1] forces_view = forces
            cp2k_get_forces(self._force_env, &forces_view[0], forces_view.shape[0])
            return forces

    def calc_energy_force(self):
        cp2k_calc_energy_force(self._force_env)

    def calc_energy(self):
        cp2k_calc_energy(self._force_env)

    property potential_energy:
        def __get__(self):
            cdef double e_pot
            cp2k_get_potential_energy(self._force_env, &e_pot)
            return e_pot
