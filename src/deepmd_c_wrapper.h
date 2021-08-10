/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2021 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

/*******************************************************************************
 * \brief A minimal wrapper for DeepMD-kit C++ interface.
 * \author Yongbin Zhuang and Yunpei Liu
 ******************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif
struct NNP;
typedef struct NNP nnp;

nnp *create_nnp(char *model);

void delete_nnp(nnp *n);

void compute_nnp(nnp *n, int *vecsize, double *dener, double *dforce,
                 double *dvirial, double *datom_ener, double *datom_virial,
                 double *dcoord_, int *datype_, double *dbox);

#ifdef __cplusplus
}
#endif
