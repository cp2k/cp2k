#ifndef GRID_COLLOCATE_GPU__
#define GRID_COLLOCATE_GPU__

void grid_collocate_task_list_gpu(void *const ptr,
    const bool orthorhombic, const int func, const int nlevels,
    const int npts_global[nlevels][3], const int npts_local[nlevels][3],
    const int shift_local[nlevels][3], const int border_width[nlevels][3],
    const double dh[nlevels][3][3], const double dh_inv[nlevels][3][3],
    double *grid[nlevels]);
#endif
