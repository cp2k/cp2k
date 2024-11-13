#ifndef CP2K_OPENPMD_H
#define CP2K_OPENPMD_H

#ifdef __cplusplus
extern "C" {
#endif // defined(__cplusplus)

typedef enum { openPMD_Access_create, openPMD_Access_read_only } openPMD_Access;

typedef struct openPMD_Series_opaque *openPMD_Series;

int openPMD_Series_create(
    // in
    char const *filename, openPMD_Access access,
    // out
    openPMD_Series *series);

int openPMD_Series_close(
    // in
    openPMD_Series series);

#ifdef __cplusplus
} // extern "C"
#endif // defined(__cplusplus)

#endif // !defined(CP2K_OPENPMD_H)
