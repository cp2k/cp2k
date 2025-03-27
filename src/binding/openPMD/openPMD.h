#ifndef CP2K_OPENPMD_H
#define CP2K_OPENPMD_H

#include <stdint.h>

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif // defined(__cplusplus)

typedef enum { openPMD_Access_create, openPMD_Access_read_only } openPMD_Access;

typedef enum {
  openPMD_Type_CHAR,
  openPMD_Type_UCHAR,
  openPMD_Type_SCHAR,
  openPMD_Type_SHORT,
  openPMD_Type_INT,
  openPMD_Type_LONG,
  openPMD_Type_LONGLONG,
  openPMD_Type_USHORT,
  openPMD_Type_UINT,
  openPMD_Type_ULONG,
  openPMD_Type_ULONGLONG,
  openPMD_Type_FLOAT,
  openPMD_Type_DOUBLE,
  openPMD_Type_LONG_DOUBLE,
  openPMD_Type_CFLOAT,
  openPMD_Type_CDOUBLE,
  openPMD_Type_CLONG_DOUBLE,
  openPMD_Type_STRING,
  openPMD_Type_VEC_CHAR,
  openPMD_Type_VEC_SHORT,
  openPMD_Type_VEC_INT,
  openPMD_Type_VEC_LONG,
  openPMD_Type_VEC_LONGLONG,
  openPMD_Type_VEC_UCHAR,
  openPMD_Type_VEC_USHORT,
  openPMD_Type_VEC_UINT,
  openPMD_Type_VEC_ULONG,
  openPMD_Type_VEC_ULONGLONG,
  openPMD_Type_VEC_FLOAT,
  openPMD_Type_VEC_DOUBLE,
  openPMD_Type_VEC_LONG_DOUBLE,
  openPMD_Type_VEC_CFLOAT,
  openPMD_Type_VEC_CDOUBLE,
  openPMD_Type_VEC_CLONG_DOUBLE,
  openPMD_Type_VEC_SCHAR,
  openPMD_Type_VEC_STRING,
  openPMD_Type_ARR_DBL_7,
  openPMD_Type_BOOL
} openPMD_Datatype;

typedef struct openPMD_Attributable_opaque *openPMD_Attributable;

typedef struct openPMD_Series_opaque *openPMD_Series;

typedef struct openPMD_Iteration_opaque *openPMD_Iteration;

typedef struct openPMD_Mesh_opaque *openPMD_Mesh;

typedef struct openPMD_ParticleSpecies_opaque *openPMD_ParticleSpecies;

typedef struct openPMD_RecordComponent_opaque *openPMD_MeshRecordComponent;

typedef struct openPMD_RecordComponent_opaque *openPMD_RecordComponent;

typedef struct openPMD_Record_opaque *openPMD_Record;

typedef struct openPMD_DynamicMemoryView_opaque *openPMD_DynamicMemoryView;

// Actually uint64_t, but ISO C bindings for Fortran only have that as a GNU
// extension, so we use signed integers and convert
typedef uint64_t openPMD_Iteration_Index_t;

/*******************
 * Series members. *
 *******************/

int openPMD_Series_create(
    // in
    char const *filename, openPMD_Access access, char const *config,
    // out
    openPMD_Series *series);

int openPMD_Series_create_mpi(
    // in
    char const *filename, openPMD_Access access, MPI_Comm comm,
    char const *config,
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

int openPMD_Series_upcast_to_Attributable(
    // in
    openPMD_Series series,
    // out
    openPMD_Attributable *attr);

int openPMD_Series_present(openPMD_Series series);

char const *openPMD_get_default_extension();

/*************************
 * Attributable members. *
 *************************/

int openPMD_attributable_set_attribute_vec_int(openPMD_Attributable attr,
                                               char const *attr_name,
                                               int const *begin, int length);

int openPMD_attributable_set_attribute_vec_string(openPMD_Attributable attr,
                                                  char const *attr_name,
                                                  char const **begin,
                                                  int length);

int openPMD_Attributable_series_flush(openPMD_Attributable attr);

/**********************
 * Iteration members. *
 **********************/

int openPMD_Iteration_upcast_to_Attributable(
    // in
    openPMD_Iteration iteration,
    // out
    openPMD_Attributable *attr);

int openPMD_Iteration_get_mesh(
    // in
    openPMD_Iteration iteration, char const *name,
    // out
    openPMD_Mesh *mesh);

int openPMD_Iteration_get_particle_species(
    // in
    openPMD_Iteration iteration, char const *name,
    // out
    openPMD_ParticleSpecies *particle_species);

/****************
 * Mesh members *
 ****************/

int openPMD_Mesh_upcast_to_RecordComponent(
    // in
    openPMD_Mesh mesh,
    // out
    openPMD_RecordComponent *rc);

int openPMD_Mesh_upcast_to_MeshRecordComponent(
    // in
    openPMD_Mesh mesh,
    // out
    openPMD_MeshRecordComponent *mrc);

int openPMD_Mesh_set_axis_labels(
    // in
    openPMD_Mesh mesh, char const **labels, int len_labels, int invert);

int openPMD_Mesh_setGridGlobalOffset(
    // in
    openPMD_Mesh mesh, double const *labels, int len_labels, int invert);

int openPMD_Mesh_setGridSpacing(
    // in
    openPMD_Mesh mesh, double const *labels, int len_labels, int invert);

int openPMD_Mesh_setPosition(
    // in
    openPMD_Mesh mesh, double const *labels, int len_labels, int invert);

/***************************
 * RecordComponent members *
 ***************************/

int openPMD_RecordComponent_upcast_to_Attributable(
    // in
    openPMD_RecordComponent rc,
    // out
    openPMD_Attributable *attr);

int openPMD_RecordComponent_resetDataset(
    // in
    openPMD_RecordComponent rc, openPMD_Datatype, int dimensions,
    int const *extent, int invert, char const *cfg);

int openPMD_RecordComponent_makeEmpty(
    // in
    openPMD_RecordComponent rc, openPMD_Datatype, int dimensions);

int openPMD_RecordComponent_makeConstant(
    // in
    openPMD_RecordComponent rc, openPMD_Datatype, int dimensions,
    int const *extent, int invert, void const *value);

int openPMD_RecordComponent_storeChunk(
    // in
    openPMD_RecordComponent rc, openPMD_Datatype, int dimensions,
    int const *offset, int const *extent, int invert, void *data);

int openPMD_RecordComponent_storeChunkSpan(
    // in
    openPMD_RecordComponent rc, openPMD_Datatype, int dimensions,
    int const *offset, int const *extent, int invert,
    // out
    openPMD_DynamicMemoryView *);

int openPMD_DynamicMemoryView_resolve(
    // in
    openPMD_DynamicMemoryView, int deallocate,
    // out
    void **write_buffer);

/*******************************
 * MeshRecordComponent members *
 *******************************/

int openPMD_MeshRecordComponent_upcast_to_RecordComponent(
    // in
    openPMD_MeshRecordComponent mrc,
    // out
    openPMD_RecordComponent *rc);

/***************************
 * ParticleSpecies members *
 ***************************/

int openPMD_ParticleSpecies_get_Record(
    // in
    openPMD_ParticleSpecies, char const *name,
    // out
    openPMD_Record *);

/******************
 * Record members *
 ******************/

int openPMD_Record_upcast_to_RecordComponent(
    // in
    openPMD_Record record,
    // out
    openPMD_RecordComponent *rc);

int openPMD_Record_get_Component(
    // in
    openPMD_Record record, char const *name,
    // out
    openPMD_RecordComponent *rc);

#ifdef __cplusplus
} // extern "C"
#endif // defined(__cplusplus)

#endif // !defined(CP2K_OPENPMD_H)
