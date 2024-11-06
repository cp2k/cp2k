#ifdef __cplusplus
extern "C" {
#endif

typedef enum { openPMD_Access_create, openPMD_Access_read_only } openPMD_Access;

struct openPMD_Series_s; // opaque
typedef struct openPMD_Series_s *openPMD_Series;

int openPMD_Series_create(
    // in
    char const *filename, openPMD_Access access,
    // out
    openPMD_Series *series);

int openPMD_Series_close(
    // in
    openPMD_Series series);

#ifdef __cplusplus
}
#endif
