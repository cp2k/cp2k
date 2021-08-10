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

#ifdef __DEEPMD
#include "deepmd_c_wrapper.h"
#include "DeepPot.h"
#include <iostream>

struct NNP {
  void *obj;
};

nnp *create_nnp(char *model) {
  nnp *n;
  deepmd::DeepPot *obj;
  n = (typeof(n))malloc(sizeof(*n));
  obj = new deepmd::DeepPot(model);
  n->obj = obj;
  //      cout << "this means wrap successfully" <<endl;
  return n;
}
void delete_nnp(nnp *n) { free(n); }
void compute_nnp(nnp *n, int *vecsize, double *dener, double *dforce,
                 double *dvirial, double *datom_ener, double *datom_virial,
                 double *dcoord_, int *datype_, double *dbox) {

  deepmd::DeepPot *obj;
  if (n == NULL)
    return;
  obj = static_cast<deepmd::DeepPot *>(n->obj);
  // covert array to vector
  double ener = 0.0;
  int vsize = *vecsize;
  std::vector<double> force_(vsize * 3, 0.0);
  std::vector<double> virial(9, 0.0);
  std::vector<double> atom_energy_(vsize, 0.0);
  std::vector<double> atom_virial_(vsize * 9, 0.0);
  std::vector<double> coord_(dcoord_, dcoord_ + vsize * 3);
  std::vector<double> box(dbox, dbox + 9);
  std::vector<int> atype_(datype_, datype_ + vsize);
  //	cout << "define ok" << endl;

  obj->compute(ener, force_, virial, atom_energy_, atom_virial_, coord_, atype_,
               box);
  //	cout << "input ok" << endl;
  *dener = ener;
  //	cout << "energy is " << *dener << endl;
  for (int i = 0; i < vsize * 3; i++) {
    dforce[i] = force_[i];
  }
  for (int i = 0; i < 9; i++) {
    dvirial[i] = virial[i];
  }
  for (int i = 0; i < vsize; i++) {
    datom_ener[i] = atom_energy_[i];
  }
  for (int i = 0; i < vsize * 9; i++) {
    datom_virial[i] = atom_virial_[i];
  }
  //	cout << "this means vector function wrap successfully" << endl;
}
#endif
