#ifndef CP2K_OPENPMD_AUXILIARY_H
#define CP2K_OPENPMD_AUXILIARY_H

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif // defined(__cplusplus)

/* Wrapper for MPI_Comm_f2c to ensure it's a compiled symbol. */
void CP2K_MPI_Comm_f2c(
    // in
    MPI_Fint comm_f,
    // out
    MPI_Comm *comm_c);

size_t CP2K_strlen(char const *);

void CP2K_free(void *);

void CP2K_log(char const *);

#ifdef __cplusplus
} // extern "C"
#endif // defined(__cplusplus)

#endif // !defined(CP2K_OPENPMD_AUXILIARY_H)
