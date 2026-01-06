/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

// Single-file library, logically containing 4 files:
// openPMD.h openPMD.cpp auxiliary.h auxiliary.c
// Implements C bindings for a subselection of the openPMD-api.

#ifdef __OPENPMD

///////////////
// openPMD.h //
///////////////

// #ifndef CP2K_OPENPMD_H
// #define CP2K_OPENPMD_H

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

int openPMD_Series_get_Iteration(
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

int openPMD_Attributable_series_flush(openPMD_Attributable attr, char const *);

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

int openPMD_Iteration_close(
    // in
    openPMD_Iteration iteration);

int openPMD_Iteration_closed(
    // in
    openPMD_Iteration iteration);

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

/***********
 * Helpers *
 ***********/

char *openPMD_json_merge(char const *into, char const *from,
                         MPI_Comm maybe_comm);

#ifdef __cplusplus
} // extern "C"
#endif // defined(__cplusplus)

// #endif // !defined(CP2K_OPENPMD_H)

/////////////////
// openPMD.cpp //
/////////////////

#include <openPMD/openPMD.hpp>
#if !OPENPMDAPI_VERSION_GE(0, 17, 0)
#include <openPMD/auxiliary/Filesystem.hpp>
#endif

#include <any>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <variant>

namespace implementation {
namespace {
struct DynamicMemoryView {
  openPMD_Datatype dtype;
  std::any memory_view;
};

constexpr auto datatype_c_to_cxx(openPMD_Datatype dt) -> openPMD::Datatype {
  switch (dt) {
  case openPMD_Type_CHAR:
    return openPMD::Datatype::CHAR;
  case openPMD_Type_UCHAR:
    return openPMD::Datatype::UCHAR;
  case openPMD_Type_SCHAR:
    return openPMD::Datatype::SCHAR;
  case openPMD_Type_SHORT:
    return openPMD::Datatype::SHORT;
  case openPMD_Type_INT:
    return openPMD::Datatype::INT;
  case openPMD_Type_LONG:
    return openPMD::Datatype::LONG;
  case openPMD_Type_LONGLONG:
    return openPMD::Datatype::LONGLONG;
  case openPMD_Type_USHORT:
    return openPMD::Datatype::USHORT;
  case openPMD_Type_UINT:
    return openPMD::Datatype::UINT;
  case openPMD_Type_ULONG:
    return openPMD::Datatype::ULONG;
  case openPMD_Type_ULONGLONG:
    return openPMD::Datatype::ULONGLONG;
  case openPMD_Type_FLOAT:
    return openPMD::Datatype::FLOAT;
  case openPMD_Type_DOUBLE:
    return openPMD::Datatype::DOUBLE;
  case openPMD_Type_LONG_DOUBLE:
    return openPMD::Datatype::LONG_DOUBLE;
  case openPMD_Type_CFLOAT:
    return openPMD::Datatype::CFLOAT;
  case openPMD_Type_CDOUBLE:
    return openPMD::Datatype::CDOUBLE;
  case openPMD_Type_CLONG_DOUBLE:
    return openPMD::Datatype::CLONG_DOUBLE;
  case openPMD_Type_STRING:
    return openPMD::Datatype::STRING;
  case openPMD_Type_VEC_CHAR:
    return openPMD::Datatype::VEC_CHAR;
  case openPMD_Type_VEC_SHORT:
    return openPMD::Datatype::VEC_SHORT;
  case openPMD_Type_VEC_INT:
    return openPMD::Datatype::VEC_INT;
  case openPMD_Type_VEC_LONG:
    return openPMD::Datatype::VEC_LONG;
  case openPMD_Type_VEC_LONGLONG:
    return openPMD::Datatype::VEC_LONGLONG;
  case openPMD_Type_VEC_UCHAR:
    return openPMD::Datatype::VEC_UCHAR;
  case openPMD_Type_VEC_USHORT:
    return openPMD::Datatype::VEC_USHORT;
  case openPMD_Type_VEC_UINT:
    return openPMD::Datatype::VEC_UINT;
  case openPMD_Type_VEC_ULONG:
    return openPMD::Datatype::VEC_ULONG;
  case openPMD_Type_VEC_ULONGLONG:
    return openPMD::Datatype::VEC_ULONGLONG;
  case openPMD_Type_VEC_FLOAT:
    return openPMD::Datatype::VEC_FLOAT;
  case openPMD_Type_VEC_DOUBLE:
    return openPMD::Datatype::VEC_DOUBLE;
  case openPMD_Type_VEC_LONG_DOUBLE:
    return openPMD::Datatype::VEC_LONG_DOUBLE;
  case openPMD_Type_VEC_CFLOAT:
    return openPMD::Datatype::VEC_CFLOAT;
  case openPMD_Type_VEC_CDOUBLE:
    return openPMD::Datatype::VEC_CDOUBLE;
  case openPMD_Type_VEC_CLONG_DOUBLE:
    return openPMD::Datatype::VEC_CLONG_DOUBLE;
  case openPMD_Type_VEC_SCHAR:
    return openPMD::Datatype::VEC_SCHAR;
  case openPMD_Type_VEC_STRING:
    return openPMD::Datatype::VEC_STRING;
  case openPMD_Type_ARR_DBL_7:
    return openPMD::Datatype::ARR_DBL_7;
  case openPMD_Type_BOOL:
    return openPMD::Datatype::BOOL;
  }
  return openPMD::Datatype::UNDEFINED;
}

constexpr auto access_c_to_cxx(openPMD_Access access) -> openPMD::Access {
  switch (access) {
  case openPMD_Access_create:
    return openPMD::Access::CREATE;
  case openPMD_Access_read_only:
    return openPMD::Access::READ_ONLY;
  }
  // unreachable
  return static_cast<openPMD::Access>(0);
}

template <typename... Args>
auto Series_create(openPMD_Series *series_param,
                   Args &&...constructor_args) -> int {
  auto series = reinterpret_cast<openPMD::Series **>(series_param);
  try {
    *series = new openPMD::Series(std::forward<Args>(constructor_args)...);
  } catch (std::exception const &e) {
    std::cout << "[Series_create] Caught error: '" << e.what() << "'\n";
    delete *series;
    return 1;
  } catch (...) {
    std::cout << "[Series_create] Caught unknown error.\n";
    delete *series;
    return 1;
  }
  return 0;
}

/*
 * Use to static_cast<> from a C++ subclass From to its parent To.
 * FromOpaque and ToOpaque are the corresponding opaque types from
 * the C header, pointers need to be reinterpret_cast<>ed to/from the
 * actual underlying C++ types.
 */
template <typename From, typename To, typename FromOpaque, typename ToOpaque>
int do_upcast(FromOpaque from_param, ToOpaque *to_param) {
  auto from = reinterpret_cast<From *>(from_param);
  auto to = reinterpret_cast<To **>(to_param);
  auto from_upcasted = static_cast<To *>(from);
  *to = from_upcasted;
  return 0;
}

template <typename res_t = void, typename T>
auto pointer_to_vector(T *ptr, size_t len, bool invert) {
  using resolved_type = std::conditional_t<std::is_same_v<res_t, void>,
                                           std::remove_cv_t<T>, res_t>;
  if (!invert) {
    return std::vector<resolved_type>{ptr, ptr + len};
  } else {
    return std::vector<resolved_type>{std::reverse_iterator(ptr + len),
                                      std::reverse_iterator(ptr)};
  }
}

template <typename res_t = void, typename T>
auto pointer_to_vector(T *ptr, size_t len, int invert) {
  return pointer_to_vector<res_t, T>(ptr, len, static_cast<bool>(invert));
}

#if !OPENPMDAPI_VERSION_GE(0, 17, 0)
auto resolveFilename(char const *unparsed) -> std::optional<std::string> {
  char const *current_char = unparsed;
  // First, trim front
  while (std::isspace(*current_char) && *current_char != '\0' &&
         *current_char != '@') {
    ++current_char;
  }
  if (*current_char == '@') {
    // Now, trim back
    auto begin_filename = current_char + 1;
    auto end_filename = begin_filename;
    size_t i = 0;
    while (*end_filename != '\0') {
      ++end_filename;
      ++i;
      if (i > FILENAME_MAX) {
        throw std::runtime_error(
            "resolveFilename: Longer than maximum allowed filename.");
      }
    }
    // *end_filename is now equal to '\0'
    --end_filename;
    while (std::isspace(*end_filename)) {
      --end_filename;
    }
    // need past-the-end iterator
    ++end_filename;
    return std::string(begin_filename, end_filename);
  } else {
    return std::nullopt;
  }
}

auto resolveFileContent(char const *maybeStringMaybeFilename,
                        std::optional<MPI_Comm> maybe_comm)
    -> std::variant<char const *, std::string> {
  auto filename = resolveFilename(maybeStringMaybeFilename);
  if (!filename.has_value()) {
    return maybeStringMaybeFilename;
  }
  auto serialImplementation = [&]() {
    std::fstream handle;
    handle.open(*filename, std::ios_base::in);
    std::stringstream stream;
    stream << handle.rdbuf();
    if (!handle.good()) {
      throw std::runtime_error("Failed acessing file '" + *filename + "'.");
    }
    handle.close();
    return stream.str();
  };
#if openPMD_HAVE_MPI
  auto parallelImplementation = [&](MPI_Comm comm) {
    return openPMD::auxiliary::collective_file_read(*filename, comm);
  };
  if (maybe_comm.has_value()) {
    return parallelImplementation(*maybe_comm);
  } else {
    return serialImplementation();
  }
#else
  return serialImplementation();
#endif
}
#endif

struct RecordComponent_makeConstant {
  template <typename Type>
  static int call(openPMD::RecordComponent &rc, void const *value_in) {
    if (value_in) {
      auto value = static_cast<Type const *>(value_in);
      rc.makeConstant(*value);
    } else {
      rc.makeConstant<Type>(Type{});
    }
    return 0;
  }
  static constexpr char const *errorMsg = "RecordComponent_makeConstant";
};

struct RecordComponent_storeChunk {
  template <typename Type>
  static int call(openPMD::RecordComponent &rc, openPMD::Offset const &o,
                  openPMD::Extent const &e, void *data) {
    rc.storeChunkRaw<Type>(static_cast<Type *>(data), o, e);
    return 0;
  }
  static constexpr char const *errorMsg = "RecordComponent_storeChunk";
};

struct RecordComponent_storeChunkSpan {
  template <typename Type>
  static int call(openPMD::RecordComponent &rc, openPMD::Offset const &o,
                  openPMD::Extent const &e, std::any &memory_view) {
    auto buffer = rc.storeChunk<Type>(o, e);
    memory_view =
        std::make_any<openPMD::DynamicMemoryView<Type>>(std::move(buffer));
    return 0;
  }
  static constexpr char const *errorMsg = "RecordComponent_storeChunkSpan";
};

struct DynamicMemoryView_resolve {
  template <typename Type> static void *call(DynamicMemoryView &memory_view) {
    auto &get_buffer = std::any_cast<openPMD::DynamicMemoryView<Type> &>(
        memory_view.memory_view);
    auto span = get_buffer.currentBuffer();
    return span.data();
  }
  static constexpr char const *errorMsg = "DynamicMemoryView_resolve";
};
} // namespace
} // namespace implementation

extern "C" {
int openPMD_Series_create(char const *filename, openPMD_Access access,
                          char const *config, openPMD_Series *series_param) {
  if (!config) {
    config = "";
  }
  return implementation::Series_create(
      series_param, filename, implementation::access_c_to_cxx(access), config);
}

int openPMD_Series_create_mpi(char const *filename, openPMD_Access access,
                              MPI_Comm comm, char const *config,
                              openPMD_Series *series_param) {
  if (!config) {
    config = "";
  }
  return implementation::Series_create(series_param, filename,
                                       implementation::access_c_to_cxx(access),
                                       comm, config);
}

int openPMD_Series_close(openPMD_Series series_param) {
  auto series = reinterpret_cast<openPMD::Series *>(series_param);
  series->close();
  delete series;
  return 0;
}

int openPMD_Series_write_Iteration(openPMD_Series series_param,
                                   openPMD_Iteration_Index_t index,
                                   openPMD_Iteration *iteration) {
  auto series = reinterpret_cast<openPMD::Series *>(series_param);
  auto &res_iteration = series->writeIterations()[index];
  *reinterpret_cast<openPMD::Iteration **>(iteration) = &res_iteration;
  return 0;
}

int openPMD_Series_get_Iteration(openPMD_Series series_param,
                                 openPMD_Iteration_Index_t index,
                                 openPMD_Iteration *iteration) {
  auto series = reinterpret_cast<openPMD::Series *>(series_param);
  auto &res_iteration = series->iterations[index];
  *reinterpret_cast<openPMD::Iteration **>(iteration) = &res_iteration;
  return 0;
}

int openPMD_Series_upcast_to_Attributable(
    // in
    openPMD_Series series,
    // out
    openPMD_Attributable *attr) {
  return implementation::do_upcast<openPMD::Series, openPMD::Attributable>(
      series, attr);
}

int openPMD_Series_present(openPMD_Series series_param) {
  auto series = reinterpret_cast<openPMD::Series *>(series_param);
  bool res = series && *series;
  return res ? 1 : 0;
}

char const *openPMD_get_default_extension() {
  constexpr static char const *preferred_order[] = {"bp5", "bp4", "bp", "h5",
                                                    "json"};
  auto const &available_extensions = openPMD::getFileExtensions();
  for (auto const &ext : preferred_order) {
    if (std::find(available_extensions.begin(), available_extensions.end(),
                  ext) != available_extensions.end()) {
      return ext;
    }
  }
  return nullptr;
}

int openPMD_attributable_set_attribute_vec_int(openPMD_Attributable attr,
                                               char const *attr_name,
                                               const int *begin, int length) {
  auto attributable = reinterpret_cast<openPMD::Attributable *>(attr);
  attributable->setAttribute(attr_name,
                             std::vector<int>{begin, begin + size_t(length)});
  return 0;
}

int openPMD_attributable_set_attribute_vec_string(openPMD_Attributable attr,
                                                  char const *attr_name,
                                                  char const **begin,
                                                  int length) {
  auto attributable = reinterpret_cast<openPMD::Attributable *>(attr);
  attributable->setAttribute(
      attr_name, std::vector<std::string>{begin, begin + size_t(length)});
  return 0;
}

int openPMD_Attributable_series_flush(openPMD_Attributable attr,
                                      char const *config) {
  auto attributable = reinterpret_cast<openPMD::Attributable *>(attr);
  attributable->seriesFlush(config ? config : "");
  return 0;
}

int openPMD_Iteration_upcast_to_Attributable(
    // in
    openPMD_Iteration iteration,
    // out
    openPMD_Attributable *attr) {
  return implementation::do_upcast<openPMD::Iteration, openPMD::Attributable>(
      iteration, attr);
}

int openPMD_Iteration_get_mesh(
    // in
    openPMD_Iteration iteration_param, char const *name,
    // out
    openPMD_Mesh *mesh) {
  auto iteration = reinterpret_cast<openPMD::Iteration *>(iteration_param);
  auto &res = iteration->meshes[name];
  *reinterpret_cast<openPMD::Mesh **>(mesh) = &res;
  return 0;
}

int openPMD_Iteration_get_particle_species(
    // in
    openPMD_Iteration iteration_param, char const *name,
    // out
    openPMD_ParticleSpecies *particle_species) {
  auto iteration = reinterpret_cast<openPMD::Iteration *>(iteration_param);
  auto &res = iteration->particles[name];
  *reinterpret_cast<openPMD::ParticleSpecies **>(particle_species) = &res;
  return 0;
}

int openPMD_Iteration_close(openPMD_Iteration iteration_param) {
  auto iteration = reinterpret_cast<openPMD::Iteration *>(iteration_param);
  iteration->close();
  return 0;
}

int openPMD_Iteration_closed(openPMD_Iteration iteration_param) {
  auto iteration = reinterpret_cast<openPMD::Iteration *>(iteration_param);
  return iteration->closed();
}

int openPMD_Mesh_upcast_to_RecordComponent(
    // in
    openPMD_Mesh mesh_param,
    // out
    openPMD_RecordComponent *rc) {
  return implementation::do_upcast<openPMD::Mesh, openPMD::RecordComponent>(
      mesh_param, rc);
}

int openPMD_Mesh_upcast_to_MeshRecordComponent(
    // in
    openPMD_Mesh mesh,
    // out
    openPMD_MeshRecordComponent *mrc) {
  return implementation::do_upcast<openPMD::Mesh, openPMD::MeshRecordComponent>(
      mesh, mrc);
}

int openPMD_Mesh_set_axis_labels(
    // in
    openPMD_Mesh mesh_param, char const **labels, int len_labels, int invert) {
  auto mesh = reinterpret_cast<openPMD::Mesh *>(mesh_param);
  mesh->setAxisLabels(implementation::pointer_to_vector<std::string>(
      labels, len_labels, invert));
  return 0;
}

int openPMD_Mesh_setGridGlobalOffset(
    // in
    openPMD_Mesh mesh_param, double const *labels, int len_labels, int invert) {
  auto mesh = reinterpret_cast<openPMD::Mesh *>(mesh_param);
  mesh->setGridGlobalOffset(
      implementation::pointer_to_vector(labels, len_labels, invert));
  return 0;
}

int openPMD_Mesh_setGridSpacing(
    // in
    openPMD_Mesh mesh_param, double const *labels, int len_labels, int invert) {
  auto mesh = reinterpret_cast<openPMD::Mesh *>(mesh_param);
  mesh->setGridSpacing(
      implementation::pointer_to_vector(labels, len_labels, invert));
  return 0;
}

int openPMD_Mesh_setPosition(
    // in
    openPMD_Mesh mesh_param, double const *labels, int len_labels, int invert) {
  auto mesh = reinterpret_cast<openPMD::Mesh *>(mesh_param);
  mesh->setPosition(
      implementation::pointer_to_vector(labels, len_labels, invert));
  return 0;
}

int openPMD_MeshRecordComponent_upcast_to_RecordComponent(
    // in
    openPMD_MeshRecordComponent mrc,
    // out
    openPMD_RecordComponent *rc) {
  return implementation::do_upcast<openPMD::MeshRecordComponent,
                                   openPMD::RecordComponent>(mrc, rc);
}

int openPMD_RecordComponent_upcast_to_Attributable(
    // in
    openPMD_RecordComponent rc,
    // out
    openPMD_Attributable *attr) {
  return implementation::do_upcast<openPMD::RecordComponent,
                                   openPMD::Attributable>(rc, attr);
}

int openPMD_RecordComponent_resetDataset(
    // in
    openPMD_RecordComponent rc_param, openPMD_Datatype dtype, int dimensions,
    int const *extent, int invert, char const *cfg) {
  auto rc = reinterpret_cast<openPMD::RecordComponent *>(rc_param);
  openPMD::Dataset ds(
      implementation::datatype_c_to_cxx(dtype),
      implementation::pointer_to_vector<openPMD::Extent::value_type>(
          extent, dimensions, invert),
      cfg ? cfg : "{}");
  rc->resetDataset(std::move(ds));
  return 0;
}

int openPMD_RecordComponent_makeEmpty(
    // in
    openPMD_RecordComponent rc_param, openPMD_Datatype dtype, int dimensions) {
  auto rc = reinterpret_cast<openPMD::RecordComponent *>(rc_param);
  rc->makeEmpty(implementation::datatype_c_to_cxx(dtype), dimensions);
  return 0;
}

int openPMD_RecordComponent_makeConstant(
    // in
    openPMD_RecordComponent rc_param, openPMD_Datatype dt, int dimensions,
    int const *extent, int invert, void const *value) {
  openPMD_RecordComponent_resetDataset(rc_param, dt, dimensions, extent, invert,
                                       NULL);
  auto rc = reinterpret_cast<openPMD::RecordComponent *>(rc_param);

  openPMD::switchDatasetType<implementation::RecordComponent_makeConstant>(
      implementation::datatype_c_to_cxx(dt), *rc, value);

  return 0;
}

int openPMD_RecordComponent_storeChunk(
    // in
    openPMD_RecordComponent rc_param, openPMD_Datatype dt, int dimensions,
    int const *offset, int const *extent, int invert, void *data) {
  auto rc = reinterpret_cast<openPMD::RecordComponent *>(rc_param);
  return openPMD::switchDatasetType<implementation::RecordComponent_storeChunk>(
      implementation::datatype_c_to_cxx(dt), *rc,
      implementation::pointer_to_vector<openPMD::Offset::value_type>(
          offset, dimensions, invert),
      implementation::pointer_to_vector<openPMD::Extent::value_type>(
          extent, dimensions, invert),
      data);
}

int openPMD_RecordComponent_storeChunkSpan(
    // in
    openPMD_RecordComponent rc_param, openPMD_Datatype dt, int dimensions,
    int const *offset, int const *extent, int invert,
    // out
    openPMD_DynamicMemoryView *memory_view_param) {
  auto rc = reinterpret_cast<openPMD::RecordComponent *>(rc_param);
  auto memory_view =
      reinterpret_cast<implementation::DynamicMemoryView **>(memory_view_param);
  *memory_view = new implementation::DynamicMemoryView();
  (*memory_view)->dtype = dt;
  return openPMD::switchDatasetType<
      implementation::RecordComponent_storeChunkSpan>(
      implementation::datatype_c_to_cxx(dt), *rc,
      implementation::pointer_to_vector<openPMD::Offset::value_type>(
          offset, dimensions, invert),
      implementation::pointer_to_vector<openPMD::Extent::value_type>(
          extent, dimensions, invert),
      (*memory_view)->memory_view);
}

int openPMD_DynamicMemoryView_resolve(
    // in
    openPMD_DynamicMemoryView memory_view_param, int deallocate,
    // out
    void **write_buffer) {
  auto memory_view =
      reinterpret_cast<implementation::DynamicMemoryView *>(memory_view_param);
  *write_buffer =
      openPMD::switchDatasetType<implementation::DynamicMemoryView_resolve>(
          implementation::datatype_c_to_cxx(memory_view->dtype), *memory_view);
  if (deallocate) {
    delete memory_view;
  }
  return 0;
}

int openPMD_ParticleSpecies_get_Record(
    // in
    openPMD_ParticleSpecies ps_param, char const *name,
    // out
    openPMD_Record *record) {
  auto ps = reinterpret_cast<openPMD::ParticleSpecies *>(ps_param);
  auto &res = ps->operator[](name);
  *reinterpret_cast<openPMD::Record **>(record) = &res;
  return 0;
}

int openPMD_Record_upcast_to_RecordComponent(
    // in
    openPMD_Record record,
    // out
    openPMD_RecordComponent *rc) {
  return implementation::do_upcast<openPMD::Record, openPMD::RecordComponent>(
      record, rc);
}

int openPMD_Record_get_Component(
    // in
    openPMD_Record record_param, char const *name,
    // out
    openPMD_RecordComponent *rc) {
  auto record = reinterpret_cast<openPMD::Record *>(record_param);
  auto &res = record->operator[](name);
  *reinterpret_cast<openPMD::RecordComponent **>(rc) = &res;
  return 0;
}

char *openPMD_json_merge(char const *into, char const *from,
                         MPI_Comm maybe_comm) {
#if OPENPMDAPI_VERSION_GE(0, 17, 0)
  auto const res = maybe_comm ? openPMD::json::merge(into, from, maybe_comm)
                              : openPMD::json::merge(into, from);
#else
  // openPMD::json::merge in release 0.16.1 has an API bug in which it does not
  // resolve filenames, so we do it manually.
  auto comm =
      maybe_comm ? std::make_optional<MPI_Comm>(maybe_comm) : std::nullopt;
  auto const into_resolved = implementation::resolveFileContent(into, comm);
  auto const from_resolved = implementation::resolveFileContent(from, comm);
  auto const res = std::visit(
      [&](auto &&i) {
        return std::visit(
            [&](auto &&f) {
              return openPMD::json::merge(
                  // poor man's std::forward()
                  static_cast<decltype(i)>(i), static_cast<decltype(f)>(f));
            },
            from_resolved);
      },
      into_resolved);
#endif
  char *c_res = static_cast<char *>(malloc(sizeof(char) * (res.size() + 1)));
  std::copy_n(res.c_str(), res.size() + 1, c_res);
  return c_res;
}
}

/////////////////
// auxiliary.h //
/////////////////

// #ifndef CP2K_OPENPMD_AUXILIARY_H
// #define CP2K_OPENPMD_AUXILIARY_H

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

#ifdef __cplusplus
} // extern "C"
#endif // defined(__cplusplus)

// #endif // !defined(CP2K_OPENPMD_AUXILIARY_H)

/////////////////
// auxiliary.c //
/////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif // defined(__cplusplus)

void CP2K_MPI_Comm_f2c(MPI_Fint comm_f, MPI_Comm *comm_c) {
  *comm_c = MPI_Comm_f2c(comm_f);
}

size_t CP2K_strlen(char const *s) { return strlen(s); }

void CP2K_free(void *ptr) { free(ptr); }

#ifdef __cplusplus
}
#endif // defined(__cplusplus)
#endif // defined(__OPENPMD)
