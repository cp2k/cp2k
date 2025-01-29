#ifndef CP2K_OPENPMD_H
#define CP2K_OPENPMD_H

#include <stdint.h>

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif // defined(__cplusplus)

typedef enum { openPMD_Access_create, openPMD_Access_read_only } openPMD_Access;

typedef struct openPMD_Attributable_opaque *openPMD_Attributable;

typedef struct openPMD_Series_opaque *openPMD_Series;

typedef struct openPMD_Iteration_opaque *openPMD_Iteration;

// Actually uint64_t, but ISO C bindings for Fortran only have that as a GNU
// extension, so we use signed integers and convert
typedef int64_t openPMD_Iteration_Index_t;

/*******************
 * Series members. *
 *******************/

int openPMD_Series_create(
    // in
    char const *filename, openPMD_Access access,
    // out
    openPMD_Series *series);

int openPMD_Series_create_mpi(
    // in
    char const *filename, openPMD_Access access, MPI_Comm comm,
    // out
    openPMD_Series *series);

int openPMD_Series_close(
    // in
    openPMD_Series series);

int openPMD_Series_write_Iteration(
    // in
    openPMD_Series series, openPMD_Iteration_Index_t index,
    // out
    openPMD_Iteration *iteration);

int openPMD_Series_upcast_to_attributable(
    // in
    openPMD_Series series,
    // out
    openPMD_Attributable *attr);

char const *openPMD_get_default_extension();

/*************************
 * Attributable members. *
 *************************/

int openPMD_attributable_set_attribute_vec_int(openPMD_Attributable attr,
                                               char const *attr_name,
                                               int const *begin, int length);

/**********************
 * Iteration members. *
 **********************/

int openPMD_Iteration_upcast_to_attributable(
    // in
    openPMD_Iteration iteration,
    // out
    openPMD_Attributable *attr);

#ifdef __cplusplus
} // extern "C"
#endif // defined(__cplusplus)

#endif // !defined(CP2K_OPENPMD_H)
