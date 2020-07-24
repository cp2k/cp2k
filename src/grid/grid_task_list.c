/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2020  CP2K developers group                         *
 *****************************************************************************/

#include "grid_task_list.h"
#include "grid_common.h"
#include "grid_collocate_cpu.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <omp.h>

// *****************************************************************************
typedef struct {
    int nset;
    int nsgf;
    int maxco;
    int maxpgf;
    int* lmin;
    int* lmax;
    int* npgf;
    int* nsgf_set;
    int* first_sgf;
    double* sphi;
    double* zet;
} intl_basis_set_t;


// *****************************************************************************
void grid_create_basis_set(const int nset,
                           const int nsgf,
                           const int maxco,
                           const int maxpgf,
                           const int lmin[nset],
                           const int lmax[nset],
                           const int npgf[nset],
                           const int nsgf_set[nset],
                           const int first_sgf[nset],
                           const double sphi[nsgf][maxco],
                           const double zet[nset][maxpgf],
                           grid_basis_set_t* basis_set_ext) {

    intl_basis_set_t* basis_set = malloc(sizeof(intl_basis_set_t));

    basis_set->nset = nset;
    basis_set->nsgf = nsgf;
    basis_set->maxco = maxco;
    basis_set->maxpgf = maxpgf;

    size_t size = nset * sizeof(int);
    basis_set->lmin = malloc(size);
    memcpy(basis_set->lmin, lmin, size);
    basis_set->lmax = malloc(size);
    memcpy(basis_set->lmax, lmax, size);
    basis_set->npgf = malloc(size);
    memcpy(basis_set->npgf, npgf, size);
    basis_set->nsgf_set = malloc(size);
    memcpy(basis_set->nsgf_set, nsgf_set, size);
    basis_set->first_sgf = malloc(size);
    memcpy(basis_set->first_sgf, first_sgf, size);
    size = nsgf * maxco * sizeof(double);
    basis_set->sphi = malloc(size);
    memcpy(basis_set->sphi, sphi, size);
    size = nset * maxpgf * sizeof(double);
    basis_set->zet = malloc(size);
    memcpy(basis_set->zet, zet, size);

    basis_set_ext->internal = basis_set;
}


// *****************************************************************************
void grid_free_basis_set(grid_basis_set_t basis_set_ext) {
    intl_basis_set_t* basis_set = (intl_basis_set_t*) basis_set_ext.internal;

    free(basis_set->lmin);
    free(basis_set->lmax);
    free(basis_set->npgf);
    free(basis_set->nsgf_set);
    free(basis_set->first_sgf);
    free(basis_set->sphi);
    free(basis_set->zet);
    free(basis_set);
}


// *****************************************************************************
typedef struct {
    int level;
    int iatom;
    int jatom;
    int iset;
    int jset;
    int ipgf;
    int jpgf;
    int subpatch;
    int dist_type;
    int block_num;
    double radius;
    double rab[3];
} intl_task_t;


// *****************************************************************************
typedef struct {
    int ntasks;
    int nlevels;
    int natoms;
    int nkinds;
    int nblocks;
    int buffer_size;
    double* blocks_buffer;
    int* block_offsets;
    double* atom_positions;
    int* atom_kinds;
    intl_basis_set_t** basis_sets;
    intl_task_t* tasks;
    int* tasks_per_level;
    int maxco;
} intl_task_list_t;


// *****************************************************************************
void grid_create_task_list(const int ntasks,
                           const int nlevels,
                           const int natoms,
                           const int nkinds,
                           const int nblocks,
                           const int buffer_size,
                           const int block_offsets[nblocks],
                           const double atom_positions[natoms][3],
                           const int atom_kinds[natoms],
                           const grid_basis_set_t basis_sets[nkinds],
                           const int level_list[ntasks],
                           const int iatom_list[ntasks],
                           const int jatom_list[ntasks],
                           const int iset_list[ntasks],
                           const int jset_list[ntasks],
                           const int ipgf_list[ntasks],
                           const int jpgf_list[ntasks],
                           const int subpatch_list[ntasks],
                           const int dist_type_list[ntasks],
                           const int block_num_list[ntasks],
                           const double radius_list[ntasks],
                           const double rab_list[ntasks][3],
                           double** blocks_buffer,
                           grid_task_list_t* task_list_ext) {

    if (task_list_ext->internal != NULL) {
        // TODO this is an opportunity to reuse some buffers
        grid_free_task_list(*task_list_ext);
    }

    intl_task_list_t* task_list = malloc(sizeof(intl_task_list_t));

    task_list->ntasks = ntasks;
    task_list->nlevels = nlevels;
    task_list->natoms = natoms;
    task_list->nkinds = nkinds;
    task_list->nblocks = nblocks;
    task_list->buffer_size = buffer_size;

    size_t size = buffer_size * sizeof(double);
    task_list->blocks_buffer = malloc(size);

    size = nblocks * sizeof(int);
    task_list->block_offsets = malloc(size);
    memcpy(task_list->block_offsets, block_offsets, size);

    size = 3 * natoms * sizeof(double);
    task_list->atom_positions = malloc(size);
    memcpy(task_list->atom_positions, atom_positions, size);

    size = natoms * sizeof(int);
    task_list->atom_kinds = malloc(size);
    memcpy(task_list->atom_kinds, atom_kinds, size);

    size = nkinds * sizeof(intl_basis_set_t*);
    task_list->basis_sets = malloc(size);
    memcpy(task_list->basis_sets, basis_sets, size);

    size = ntasks * sizeof(intl_task_t);
    task_list->tasks = malloc(size);
    for (int i=0; i < ntasks; i++) {
        task_list->tasks[i].level = level_list[i];
        task_list->tasks[i].iatom = iatom_list[i];
        task_list->tasks[i].jatom = jatom_list[i];
        task_list->tasks[i].iset = iset_list[i];
        task_list->tasks[i].jset = jset_list[i];
        task_list->tasks[i].ipgf = ipgf_list[i];
        task_list->tasks[i].jpgf = jpgf_list[i];
        task_list->tasks[i].subpatch = subpatch_list[i];
        task_list->tasks[i].dist_type = dist_type_list[i];
        task_list->tasks[i].block_num = block_num_list[i];
        task_list->tasks[i].radius = radius_list[i];
        task_list->tasks[i].rab[0] = rab_list[i][0];
        task_list->tasks[i].rab[1] = rab_list[i][1];
        task_list->tasks[i].rab[2] = rab_list[i][2];
    }

    // Count tasks per level.
    size = nlevels * sizeof(int);
    task_list->tasks_per_level = malloc(size);
    memset(task_list->tasks_per_level, 0, size);
    for (int i=0; i < ntasks; i++) {
        task_list->tasks_per_level[level_list[i]-1]++;
        assert(i==0 || level_list[i] >= level_list[i-1]);  // expect ordered list
    }

    // Find largest Cartesian subblock size.
    task_list->maxco = 0;
    for (int i=0; i < nkinds; i++) {
        task_list->maxco = max(task_list->maxco, task_list->basis_sets[i]->maxco);
    }

    *blocks_buffer = task_list->blocks_buffer;
    task_list_ext->internal = task_list;
}


// *****************************************************************************
void grid_free_task_list(grid_task_list_t task_list_ext) {
    intl_task_list_t* task_list = (intl_task_list_t*) task_list_ext.internal;

    free(task_list->blocks_buffer);
    free(task_list->block_offsets);
    free(task_list->atom_positions);
    free(task_list->atom_kinds);
    free(task_list->basis_sets);
    free(task_list->tasks);
    free(task_list->tasks_per_level);
    free(task_list);
}


// *****************************************************************************
static void collocate_one_grid_level(const intl_task_list_t* task_list,
                                     const int first_task,
                                     const int last_task,
                                     const bool orthorhombic,
                                     const int func,
                                     const int npts_global[3],
                                     const int npts_local[3],
                                     const int shift_local[3],
                                     const int border,
                                     const bool distributed,
                                     const double dh[3][3],
                                     const double dh_inv[3][3],
                                     double* grid) {

    // Using default(shared) because with GCC 9 the behavior around const changed:
    // https://www.gnu.org/software/gcc/gcc-9/porting_to.html
    #pragma omp parallel default(shared)
    {

    // Allocate thread local copy of the grid.
    const size_t npts_local_total = npts_local[0] * npts_local[1] * npts_local[2];
    const size_t grid_size = npts_local_total * sizeof(double);
    double* threadlocal_grid = malloc(grid_size);
    memset(threadlocal_grid, 0, grid_size);

    // Allocate pab matrix for re-use across tasks.
    const size_t max_pab_size = task_list->maxco * task_list->maxco * sizeof(double);
    double* pab = malloc(max_pab_size);

    // Initialize variables to detect when a new subblock has to be fetched.
    int prev_block_num= -1, prev_iset = -1, prev_jset = -1;

    #pragma omp for schedule(static)
    for (int itask=first_task; itask <= last_task; itask++) {
        // Define some convenient aliases.
        const intl_task_t* task = &task_list->tasks[itask];
        const int iatom = task->iatom - 1;
        const int jatom = task->jatom - 1;
        const int iset = task->iset - 1;
        const int jset = task->jset - 1;
        const int ipgf = task->ipgf -1;
        const int jpgf = task->jpgf -1;
        const int ikind = task_list->atom_kinds[iatom] - 1;
        const int jkind = task_list->atom_kinds[jatom] - 1;
        const intl_basis_set_t* ibasis = task_list->basis_sets[ikind];
        const intl_basis_set_t* jbasis = task_list->basis_sets[jkind];
        const int ncoseta = ncoset[ibasis->lmax[iset]];
        const int ncosetb = ncoset[jbasis->lmax[jset]];
        const int ncoa = ibasis->npgf[iset] * ncoseta; // size of carthesian set
        const int ncob = jbasis->npgf[jset] * ncosetb;
        const int block_num = task->block_num - 1;

        // Load subblock from buffer and decontract into Cartesian sublock pab.
        // The previous pab can be reused when only ipgf or jpgf has changed.
        if (block_num != prev_block_num || iset != prev_iset || jset != prev_jset) {
            prev_block_num = block_num; prev_iset = iset; prev_jset = jset;

            // Define some more convenient aliases.
            const int nsgf_seta = ibasis->nsgf_set[iset]; // size of spherical set
            const int nsgf_setb = jbasis->nsgf_set[jset];
            const int nsgfa = ibasis->nsgf; // size of entire spherical basis
            const int nsgfb = jbasis->nsgf;
            const int sgfa = ibasis->first_sgf[iset] - 1; // start of spherical set
            const int sgfb = jbasis->first_sgf[jset] - 1;
            const int maxcoa = ibasis->maxco;
            const int maxcob = jbasis->maxco;

            // Locate current matrix block within the buffer.
            const int block_offset = task_list->block_offsets[block_num];  // zero based
            const double* block = &task_list->blocks_buffer[block_offset];

            // Copy sub block for current sets and transpose it if needed.
            double subblock[nsgf_setb][nsgf_seta];
            if (iatom <= jatom) {
                for (int i=0; i < nsgf_setb; i++)
                    for (int j=0; j < nsgf_seta; j++)
                        subblock[i][j] = block[(i+sgfb)*nsgfa + j+sgfa];
            } else {
                for (int i=0; i < nsgf_setb; i++)
                    for (int j=0; j < nsgf_seta; j++)
                        subblock[i][j] = block[(j+sgfa)*nsgfb + i+sgfb];  // transposed
            }

            // work = MATMUL(ibasis->sphi, subblock)
            double work[nsgf_setb][ncoa];
            memset(work, 0, sizeof(double) * nsgf_setb * ncoa);
            for (int i=0; i < nsgf_setb; i++) {
              for (int j=0; j < ncoa ; j++) {
                for (int k=0; k < nsgf_seta; k++) {
                   work[i][j] += subblock[i][k] * ibasis->sphi[(k+sgfa)*maxcoa + j];
                }
              }
            }

            // double pab[ncob][ncoa]
            // pab = MATMUL(work, TRANSPOSE(jbasis->sphi))
            memset(pab, 0, sizeof(double) * ncob * ncoa);
            for (int i=0; i < ncob; i++) {
              for (int j=0; j < ncoa ; j++) {
                for (int k=0; k < nsgf_setb; k++) {
                   pab[i*ncoa + j] += work[k][j] * jbasis->sphi[(k+sgfb)*maxcob + i];
                }
              }
            }
        }  // end of block loading

        // TODO: task->dist_type should be encoded in task->subpatch.
        //       Also, do we really need distributed[level] ?
        const bool use_subpatch = distributed && task->dist_type==2;
        const double zeta = ibasis->zet[iset * ibasis->maxpgf + ipgf];
        const double zetb = jbasis->zet[jset * jbasis->maxpgf + jpgf];

        grid_collocate_pgf_product_cpu(/*orthorhombic=*/orthorhombic,
                                       /*use_subpatch=*/use_subpatch,
                                       /*subpatch=*/task->subpatch,
                                       /*border=*/border,
                                       /*func=*/func,
                                       /*la_max=*/ibasis->lmax[iset],
                                       /*la_min=*/ibasis->lmin[iset],
                                       /*lb_max=*/jbasis->lmax[jset],
                                       /*lb_min=*/jbasis->lmin[jset],
                                       /*zeta=*/zeta,
                                       /*zetb=*/zetb,
                                       /*rscale=*/(iatom==jatom)?1:2,
                                       /*dh=*/dh,
                                       /*dh_inv=*/dh_inv,
                                       /*ra=*/&task_list->atom_positions[3*iatom],
                                       /*rab=*/task->rab,
                                       /*npts_global=*/npts_global,
                                       /*npts_local=*/npts_local,
                                       /*shift_local=*/shift_local,
                                       /*radius=*/task->radius,
                                       /*o1=*/ipgf * ncoseta,
                                       /*o2=*/jpgf * ncosetb,
                                       /*n1=*/ncoa,
                                       /*n2=*/ncob,
                                       /*pab=*/(double (*)[ncoa])pab,
                                       /*grid=*/threadlocal_grid);
    }  // end of task loop

    // Merge thread local grids into shared grid.
    #pragma omp critical
    for (size_t i=0; i < npts_local_total; i++) {
        grid[i] += threadlocal_grid[i];
    }

    free(pab);
    free(threadlocal_grid);

    }  // end of omp parallel
}


// *****************************************************************************
void grid_collocate_task_list(const grid_task_list_t task_list_ext,
                              const bool orthorhombic,
                              const int func,
                              const int nlevels,
                              const int npts_global[nlevels][3],
                              const int npts_local[nlevels][3],
                              const int shift_local[nlevels][3],
                              const int border[nlevels],
                              const bool distributed[nlevels],
                              const double dh[nlevels][3][3],
                              const double dh_inv[nlevels][3][3],
                              double* grid[nlevels]) {

    const intl_task_list_t* task_list = (intl_task_list_t*) task_list_ext.internal;

    assert(task_list->nlevels == nlevels);

    int first_task = 0;
    for (int level=0; level < task_list->nlevels; level++) {
        const int last_task = first_task + task_list->tasks_per_level[level] - 1;

        collocate_one_grid_level(task_list,
                                 first_task,
                                 last_task,
                                 orthorhombic,
                                 func,
                                 npts_global[level],
                                 npts_local[level],
                                 shift_local[level],
                                 border[level],
                                 distributed[level],
                                 dh[level],
                                 dh_inv[level],
                                 grid[level]);

        first_task = last_task + 1;
    }
}

//EOF
