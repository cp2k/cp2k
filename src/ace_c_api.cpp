/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#if defined(__ACE)

// tested with lammps-user-pace-v.2023.11.25.fix2
#define COMPUTE_B_GRAD
#define EXTRA_C_PROJECTIONS

#if 0
#include <stdio.h>
#endif
#include <string>

#include "ace-evaluator/ace_c_basis.h"
#include "ace-evaluator/ace_evaluator.h"
#include "ace-evaluator/ace_recursive.h"
#include "ace-evaluator/ace_version.h"
#include "ace/ace_b_basis.h"

struct ACEData {
  ACECTildeBasisSet *basis_set;
  ACERecursiveEvaluator *ace;
  Array1D<DOUBLE_TYPE> *virial;
  Array2D<DOUBLE_TYPE> *forces;
};

bool hasEnding(std::string const &fullString, std::string const &ending) {
  if (fullString.length() >= ending.length()) {
    return (0 == fullString.compare(fullString.length() - ending.length(),
                                    ending.length(), ending));
  } else {
    return false;
  }
}

extern "C" void AcePotInitialize(int ntypec, const char *symbolsc, int nlen,
                                 const char *potential_file_name, double *rcutc,
                                 void **acedata_ptr) {

// avoid mixing C++ I/O with Fortran I/O, TODO: return this data so it can
// be printed on the Fortran side.
#if 0
  printf("---PACE initialization----\n");

  printf("ACE version: %d.%d.%d\n", VERSION_YEAR, VERSION_MONTH, VERSION_DAY);
#endif

  std::string potential_file_name_str(potential_file_name);
  // trim potential_file_name
  potential_file_name_str.erase(
      potential_file_name_str.find_last_not_of(" \n\r\t") + 1);

  std::string symbols_str(symbolsc, ntypec * 2);
  std::vector<std::string> elements;
  for (int i = 0; i < ntypec; i++) {
    auto el_str = symbols_str.substr(2 * i, 2);
    el_str.erase(el_str.find_last_not_of(" \n\r\t") + 1);
    elements.push_back(el_str);
  }
  if (ntypec != elements.size())
    throw std::runtime_error(
        "Number of elements and elements list are inconsistent");

#if 0
  printf("Number of atom types:                       %d\n", ntypec);
  printf("Element mapping:                            ");
  for (int i = 0; i < ntypec; i++)
    printf(" `%s`", elements[i].c_str());
  printf("\n");
#endif

    // Elements are contained in a string of length 2*ntypec
    // Each element has two chars in the string
    // The sequence of the elements in the string corresponds to their mapping,
    // i.e., the first element in the string is element 1, the second element in
    // the string is element 2, the third element in the string is element 3,
    // ...

#if 0
  printf("Filename:                                   '%s'\n",
         potential_file_name_str.c_str());
#endif

  ACEData *aceData;

  if (hasEnding(potential_file_name_str, ".yaml")) {
    ACEBBasisSet bBasisSet = ACEBBasisSet(potential_file_name_str);
    ACECTildeBasisSet cTildeBasisSet = bBasisSet.to_ACECTildeBasisSet();
    aceData = new ACEData;
    aceData->basis_set = new ACECTildeBasisSet(cTildeBasisSet);
  } else if (hasEnding(potential_file_name_str, ".yace")) {
    aceData = new ACEData;
    aceData->basis_set = new ACECTildeBasisSet(potential_file_name_str.c_str());
  } else {
    throw std::invalid_argument("Unrecognized file format: '" +
                                potential_file_name_str + "'");
  }
  aceData->ace = new ACERecursiveEvaluator();
  aceData->ace->set_recursive(true);

  aceData->ace->element_type_mapping.init(1 + ntypec);
  //    ace->element_type_mapping = {0,1,0}; // 0->0, 1(CP2k)-> 1(ACE),
  //    2(CP2k)-> 0(ACE)
  for (int i = 1; i <= ntypec; i++) {
    auto elemname = elements.at(i - 1);
    SPECIES_TYPE mu = aceData->basis_set->get_species_index_by_name(elemname);
    if (mu != -1) {
#if 0
      printf("Mapping CP2K atom type #%d(%s) -> ACE species type #%d\n", i,
             elemname.c_str(), mu);
#endif
      // set up CP2K atom type to ACE species mapping for ace evaluator
      aceData->ace->element_type_mapping(i) = mu;
    } else {
      throw std::runtime_error("Element " + elemname +
                               " is not supported by ACE-potential from file " +
                               potential_file_name_str);
    }
  }

  // the cutoffs of all pairs are stored in a ntypec*ntypec array
  // the index 0 corresponds to element 1
  // the index 1 corresponds to element 2
  // the index 2 corresponds to element 3, ...
  int k = 0;
  for (int i = 1; i <= ntypec; i++) {
    for (int j = 1; j <= ntypec; j++) {
      rcutc[k] = aceData->basis_set->radial_functions->cut(
          aceData->ace->element_type_mapping(i),
          aceData->ace->element_type_mapping(j));
      k++;
    }
  }

  aceData->ace->set_basis(*aceData->basis_set, 1);
  aceData->virial = new Array1D<DOUBLE_TYPE>(6, "virial");
  aceData->forces = new Array2D<DOUBLE_TYPE>(1, 3, "forces");
  *acedata_ptr = (void *)aceData;
#if 0
  printf("---Done PACE initialization----\n");
#endif
}

extern "C" void AcePotFinalize(void **acedata_ptr) {
  ACEData *aceData = (ACEData *)*acedata_ptr;
  delete aceData->basis_set;
  delete aceData->ace;
  delete aceData->virial;
  delete aceData->forces;
  delete aceData;
}

extern "C" void AcePotCompute(int natomc, int nghostc, int neic, int *neiatc,
                              int *originc, int *nlistc, int *attypec,
                              double *atposc, double *forcec, double *virialc,
                              double *energyc, void **acedata_ptr) {
  ACEData *aceData = (ACEData *)*acedata_ptr;
  // re-point double **x (LAMMPS/C style of 2D array) to atposc
  int tot_nat = natomc + nghostc;
  double **x = new double *[tot_nat];
  for (int i = 0; i < tot_nat; i++) {
    x[i] = &atposc[3 * i];
  }

  std::vector<int> numneigh(natomc, 0);
  for (int i = 0; i < natomc; i++) {
    numneigh[i] = neiatc[i + 1] - neiatc[i];
  }

  // determine the maximum number of neighbours
  int i, jnum;
  int max_jnum = 0;
  int nei = 0;
  for (i = 0; i < natomc; i++) {
    jnum = numneigh[i];
    nei = nei + jnum;
    if (jnum > max_jnum)
      max_jnum = jnum;
  }

  aceData->ace->resize_neighbours_cache(max_jnum);

  double dx, dy, dz, fx, fy, fz;

  // resize forces array
  if (aceData->forces->get_dim(0) < natomc)
    aceData->forces->resize(natomc, 3);

  aceData->forces->fill(0);
  aceData->virial->fill(0);

  // main loop over atoms
  for (i = 0; i < natomc; i++) {
    jnum = numneigh[i];

    const int *jlist = &nlistc[neiatc[i]];

    aceData->ace->compute_atom(i, x, attypec, jnum, jlist);

    energyc[i] = aceData->ace->e_atom;

    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];

      dx = x[j][0] - xtmp;
      dy = x[j][1] - ytmp;
      dz = x[j][2] - ztmp;

      fx = aceData->ace->neighbours_forces(jj, 0);
      fy = aceData->ace->neighbours_forces(jj, 1);
      fz = aceData->ace->neighbours_forces(jj, 2);

      (*aceData->forces)(i, 0) += fx;
      (*aceData->forces)(i, 1) += fy;
      (*aceData->forces)(i, 2) += fz;

      // virial f_dot_r, identical to LAMMPS virial_fdotr_compute
      (*aceData->virial)(0) += dx * fx;
      (*aceData->virial)(1) += dy * fy;
      (*aceData->virial)(2) += dz * fz;
      (*aceData->virial)(3) += dx * fy;
      (*aceData->virial)(4) += dx * fz;
      (*aceData->virial)(5) += dy * fz;

      // update forces only for real atoms
      if (j < natomc) {
        (*aceData->forces)(j, 0) -= fx;
        (*aceData->forces)(j, 1) -= fy;
        (*aceData->forces)(j, 2) -= fz;
      } else {
        // map ghost j into true_j within periodic cell
        int true_j = originc[j];
        (*aceData->forces)(true_j, 0) -= fx;
        (*aceData->forces)(true_j, 1) -= fy;
        (*aceData->forces)(true_j, 2) -= fz;
      }
    }
  }

  double ene = 0.0;
  for (int i = 0; i < natomc; i++) {
    ene += energyc[i];

    forcec[3 * i] = (*aceData->forces)(i, 0);
    forcec[3 * i + 1] = (*aceData->forces)(i, 1);
    forcec[3 * i + 2] = (*aceData->forces)(i, 2);
  }

  // copy virials
  for (int i = 0; i < 6; i++) {
    virialc[i] = (*aceData->virial)(i);
  }

  delete[] x;
}

#endif // defined(__ACE)

// EOF
