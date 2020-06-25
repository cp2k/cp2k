/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2020  CP2K developers group                         *
 *****************************************************************************/

#define _XOPEN_SOURCE 700   /* Enable POSIX 2008/13 */

#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <stdlib.h>

#include "grid_collocate_replay.h"
#include "grid_collocate_cpu.h"

// *****************************************************************************
void grid_collocate_record(const bool orthorhombic,
                           const bool use_subpatch,
                           const int subpatch,
                           const int border,
                           const int func,
                           const int la_max,
                           const int la_min,
                           const int lb_max,
                           const int lb_min,
                           const double zeta,
                           const double zetb,
                           const double rscale,
                           const double dh[3][3],
                           const double dh_inv[3][3],
                           const double ra[3],
                           const double rab[3],
                           const int npts_global[3],
                           const int npts_local[3],
                           const int shift_local[3],
                           const double radius,
                           const int o1,
                           const int o2,
                           const int n1,
                           const int n2,
                           const double pab[n2][n1],
                           const double* grid){

    static int counter = 0;
    counter++;
    char filename[100];
    snprintf(filename, sizeof(filename), "grid_collocate_%05i.task", counter);

    const int D = DECIMAL_DIG;  // In C11 we could use DBL_DECIMAL_DIG.
    FILE *fp = fopen(filename, "w+");
    fprintf(fp, "#Grid collocate task v8\n");
    fprintf(fp, "orthorhombic %i\n", orthorhombic);
    fprintf(fp, "use_subpatch %i\n", use_subpatch);
    fprintf(fp, "subpatch %i\n", subpatch);
    fprintf(fp, "border %i\n", border);
    fprintf(fp, "func %i\n", func);
    fprintf(fp, "la_max %i\n", la_max);
    fprintf(fp, "la_min %i\n", la_min);
    fprintf(fp, "lb_max %i\n", lb_max);
    fprintf(fp, "lb_min %i\n", lb_min);
    fprintf(fp, "zeta %.*e\n", D, zeta);
    fprintf(fp, "zetb %.*e\n", D, zetb);
    fprintf(fp, "rscale %.*e\n", D, rscale);
    for (int i=0; i<3; i++)
        fprintf(fp, "dh %i %.*e %.*e %.*e\n", i, D, dh[i][0], D, dh[i][1], D, dh[i][2]);
    for (int i=0; i<3; i++)
        fprintf(fp, "dh_inv %i %.*e %.*e %.*e\n", i, D, dh_inv[i][0], D, dh_inv[i][1], D, dh_inv[i][2]);
    fprintf(fp, "ra %.*e %.*e %.*e\n", D, ra[0], D, ra[1], D, ra[2]);
    fprintf(fp, "rab %.*e %.*e %.*e\n", D, rab[0], D, rab[1], D, rab[2]);
    fprintf(fp, "npts_global %i  %i %i\n", npts_global[0], npts_global[1], npts_global[2]);
    fprintf(fp, "npts_local %i  %i %i\n", npts_local[0], npts_local[1], npts_local[2]);
    fprintf(fp, "shift_local %i  %i %i\n", shift_local[0], shift_local[1], shift_local[2]);
    fprintf(fp, "radius %.*e\n", D, radius);
    fprintf(fp, "o1 %i\n", o1);
    fprintf(fp, "o2 %i\n", o2);
    fprintf(fp, "n1 %i\n", n1);
    fprintf(fp, "n2 %i\n", n2);

    for (int i=0; i < n2; i++) {
    for (int j=0; j < n1; j++) {
        fprintf(fp, "pab %i %i %.*e\n", i, j, D, pab[i][j]);
    }
    }

    const int npts_local_total = npts_local[0] * npts_local[1] * npts_local[2];

    int ngrid_nonzero = 0;
    for (int i=0; i<npts_local_total; i++) {
        if (grid[i] != 0.0) {
            ngrid_nonzero++;
        }
    }
    fprintf(fp, "ngrid_nonzero %i\n", ngrid_nonzero);

    for (int k=0; k<npts_local[2]; k++) {
    for (int j=0; j<npts_local[1]; j++) {
    for (int i=0; i<npts_local[0]; i++) {
        const double val =  grid[k*npts_local[1]*npts_local[0] + j*npts_local[0] + i];
        if (val != 0.0) {
            fprintf(fp, "grid %i %i %i %.*e\n", i, j, k, D, val);
        }
    }
    }
    }
    fprintf(fp, "#THE_END\n");
    fclose(fp);
    printf("Wrote %s\n", filename);

}

// *****************************************************************************
double grid_collocate_replay(const char* filename, const int cycles){
    FILE *fp = fopen(filename, "r");
    assert(fp != NULL && "Could not open task file.");

    char line[100], key[100];

    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(strcmp(line, "#Grid collocate task v8\n") == 0);

    int orthorhombic_i;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &orthorhombic_i) == 2);
    assert(strcmp(key, "orthorhombic") == 0);
    bool orthorhombic = orthorhombic_i;

    int use_subpatch_i;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &use_subpatch_i) == 2);
    assert(strcmp(key, "use_subpatch") == 0);
    bool use_subpatch = use_subpatch_i;

    int subpatch;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &subpatch) == 2);
    assert(strcmp(key, "subpatch") == 0);

    int border;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &border) == 2);
    assert(strcmp(key, "border") == 0);

    int func;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &func) == 2);
    assert(strcmp(key, "func") == 0);

    int la_max;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &la_max) == 2);
    assert(strcmp(key, "la_max") == 0);

    int la_min;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &la_min) == 2);
    assert(strcmp(key, "la_min") == 0);

    int lb_max;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &lb_max) == 2);
    assert(strcmp(key, "lb_max") == 0);

    int lb_min;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &lb_min) == 2);
    assert(strcmp(key, "lb_min") == 0);

    double zeta;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %le", key, &zeta) == 2);
    assert(strcmp(key, "zeta") == 0);

    double zetb;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %le", key, &zetb) == 2);
    assert(strcmp(key, "zetb") == 0);

    double rscale;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %le", key, &rscale) == 2);
    assert(strcmp(key, "rscale") == 0);

    double dh[3][3];
    for (int i=0; i<3; i++) {
        int j;
        assert(fgets(line, sizeof(line), fp) != NULL);
        assert(sscanf(line, "%99s %i %le %le %le", key, &j, &dh[i][0], &dh[i][1], &dh[i][2]) == 5);
        assert(strcmp(key, "dh") == 0 && i == j);
    }

    double dh_inv[3][3];
    for (int i=0; i<3; i++) {
        int j;
        assert(fgets(line, sizeof(line), fp) != NULL);
        assert(sscanf(line, "%99s %i %le %le %le", key, &j, &dh_inv[i][0], &dh_inv[i][1], &dh_inv[i][2]) == 5);
        assert(strcmp(key, "dh_inv") == 0 && i == j);
    }

    double ra[3];
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %le %le %le", key, &ra[0], &ra[1], &ra[2]) == 4);
    assert(strcmp(key, "ra") == 0);

    double rab[3];
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %le %le %le", key, &rab[0], &rab[1], &rab[2]) == 4);
    assert(strcmp(key, "rab") == 0);

    int npts_global[3];
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i %i %i", key, &npts_global[0], &npts_global[1], &npts_global[2]) == 4);
    assert(strcmp(key, "npts_global") == 0);

    int npts_local[3];
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i %i %i", key, &npts_local[0], &npts_local[1], &npts_local[2]) == 4);
    assert(strcmp(key, "npts_local") == 0);

    int shift_local[3];
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i %i %i", key, &shift_local[0], &shift_local[1], &shift_local[2]) == 4);
    assert(strcmp(key, "shift_local") == 0);

    double radius;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %le", key, &radius) == 2);
    assert(strcmp(key, "radius") == 0);

    int o1;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &o1) == 2);
    assert(strcmp(key, "o1") == 0);

    int o2;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &o2) == 2);
    assert(strcmp(key, "o2") == 0);

    int n1;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &n1) == 2);
    assert(strcmp(key, "n1") == 0);

    int n2;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &n2) == 2);
    assert(strcmp(key, "n2") == 0);

    double pab[n2][n1];
    for (int i=0; i<n2; i++) {
    for (int j=0; j<n1; j++) {
        int i2, j2;
        double value;
        assert(fgets(line, sizeof(line), fp) != NULL);
        assert(sscanf(line, "%99s %i %i %le", key, &i2, &j2, &value) == 4);
        assert(strcmp(key, "pab") == 0 && i == i2 && j==j2);
        pab[i][j] = value;
    }
    }

    int ngrid_nonzero;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &ngrid_nonzero) == 2);
    assert(strcmp(key, "ngrid_nonzero") == 0);

    const int npts_local_total = npts_local[0] * npts_local[1] * npts_local[2];
    const size_t sizeof_grid = sizeof(double) * npts_local_total;
    double* grid_ref = malloc(sizeof_grid);
    memset(grid_ref, 0, sizeof_grid);

    for (int n=0; n < ngrid_nonzero; n++) {
        int i, j, k;
        double value;
        assert(fgets(line, sizeof(line), fp) != NULL);
        assert(sscanf(line, "%99s %i %i %i %le", key, &i, &j, &k, &value) == 5);
        assert(strcmp(key, "grid") == 0);
        grid_ref[k*npts_local[1]*npts_local[0] + j*npts_local[0] + i] = value;
    }

    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(strcmp(line, "#THE_END\n") == 0);

    double* grid_test = malloc(sizeof_grid);
    memset(grid_test, 0, sizeof_grid);

    struct timespec start_time;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_time);

    for (int i=0; i < cycles ; i++) {
        grid_collocate_pgf_product_cpu(orthorhombic,
                                       use_subpatch,
                                       subpatch,
                                       border,
                                       func,
                                       la_max,
                                       la_min,
                                       lb_max,
                                       lb_min,
                                       zeta,
                                       zetb,
                                       rscale,
                                       dh,
                                       dh_inv,
                                       ra,
                                       rab,
                                       npts_global,
                                       npts_local,
                                       shift_local,
                                       radius,
                                       o1,
                                       o2,
                                       n1,
                                       n2,
                                       pab,
                                       grid_test);
    }

    struct timespec end_time;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_time);
    const double delta_sec = (end_time.tv_sec - start_time.tv_sec) + 1e-9 * (end_time.tv_nsec - start_time.tv_nsec);

    double max_value = 0.0;
    double max_diff = 0.0;
    for (int i=0; i < npts_local_total; i++) {
        const double ref_value = cycles * grid_ref[i];
        const double diff = fabs(grid_test[i] - ref_value);
        max_diff = fmax(max_diff, diff);
        max_value = fmax(max_value, fabs(grid_test[i]));
    }

    printf("Task: %-65s   Cycles: %e   Max value: %le   Max diff: %le   Time: %le sec\n",
           filename, (float)cycles, max_value, max_diff, delta_sec);

    free(grid_ref);
    free(grid_test);

    return max_diff;
}

//EOF
