! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

!> @file tblite/basis/slater.f90
!> Provides coefficients for contracted Gaussian type orbitals to approximate Slater
!> type orbitals.

!> Implements expansion of Slater functions to contracted Gaussian functions.
!>
!> All the coefficients are taken from:
!>
!> Robert F. Stewart, Small Gaussian Expansions of Slater-Type Orbitals,
!> J. Chem. Phys. 52, 431-438 (1970). DOI: 10.1063/1.1672702
module tblite_basis_slater
   use mctc_env, only : wp
   use mctc_io_constants, only : pi
   use tblite_basis_type, only : cgto_type
   implicit none
   private

   public :: slater_to_gauss


   interface slater_to_gauss
      module procedure :: slater_to_gauss_cgto
      module procedure :: slater_to_gauss_array
   end interface slater_to_gauss


   !> Number of functions
   integer, parameter :: nf = 15

   !> Two over pi
   real(wp), parameter :: top = 2.0_wp / pi

   !> Double factorial, see OEIS A001147
   real(wp), parameter :: dfactorial(8) = &
      & [1.0_wp,1.0_wp,3.0_wp,15.0_wp,105.0_wp,945.0_wp,10395.0_wp,135135.0_wp]

   !> Exponents from first row Table I-V.
   real(wp), parameter :: pAlpha1(nf) = [ &
      & 2.709498091e-1_wp, & ! 1s
      & 1.012151084e-1_wp, & ! 2s
      & 5.296881757e-2_wp, & ! 3s
      & 3.264600274e-2_wp, & ! 4s
      & 2.216912938e-2_wp, & ! 5s
      & 1.759666885e-1_wp, & ! 2p
      & 9.113614253e-2_wp, & ! 3p
      & 5.578350235e-2_wp, & ! 4p
      & 3.769845216e-2_wp, & ! 5p
      & 1.302270363e-1_wp, & ! 3d
      & 7.941656339e-2_wp, & ! 4d
      & 5.352200793e-2_wp, & ! 5d
      & 1.033434062e-1_wp, & ! 4f
      & 6.952785407e-2_wp, & ! 5f
      & 8.565417784e-2_wp]   ! 5g


   !> Exponents from second row Table I-V.
   real(wp), parameter :: pAlpha2(2, nf) = reshape([ &
      & 8.518186635e-1_wp, 1.516232927e-1_wp, & ! 1s
      & 1.292278611e-1_wp, 4.908584205e-2_wp, & ! 2s
      & 6.694095822e-1_wp, 5.837135094e-2_wp, & ! 3s
      & 2.441785453e-1_wp, 4.051097664e-2_wp, & ! 4s
      & 1.213425654e-1_wp, 3.133152144e-2_wp, & ! 5s
      & 4.323908358e-1_wp, 1.069439065e-1_wp, & ! 2p
      & 1.458620964e-1_wp, 5.664210742e-2_wp, & ! 3p
      & 6.190052680e-2_wp, 2.648418407e-2_wp, & ! 4p
      & 2.691294191e-1_wp, 3.980805011e-2_wp, & ! 5p
      & 2.777427345e-1_wp, 8.336507714e-2_wp, & ! 3d
      & 1.330958892e-1_wp, 5.272119659e-2_wp, & ! 4d
      & 6.906014388e-2_wp, 3.399457777e-2_wp, & ! 5d
      & 2.006693538e-1_wp, 6.865384900e-2_wp, & ! 4f
      & 1.156094555e-1_wp, 4.778940916e-2_wp, & ! 5f
      & 1.554531559e-1_wp, 5.854079811e-2_wp],& ! 5g
      & shape(pAlpha2))

   !> Coefficients from second row Table I-V.
   real(wp), parameter :: pCoeff2(2, nf) = reshape([ &
      & 4.301284983e-1_wp, 6.789135305e-1_wp, & ! 1s
      & 7.470867124e-1_wp, 2.855980556e-1_wp, & ! 2s
      &-1.529645716e-1_wp, 1.051370110e+0_wp, & ! 3s
      &-3.046656896e-1_wp, 1.146877294e+0_wp, & ! 4s
      &-5.114756049e-1_wp, 1.307377277e+0_wp, & ! 5s
      & 4.522627513e-1_wp, 6.713122642e-1_wp, & ! 2p
      & 5.349653144e-1_wp, 5.299607212e-1_wp, & ! 3p
      & 8.743116767e-1_wp, 1.513640107e-1_wp, & ! 4p
      &-1.034227010e-1_wp, 1.033376378e+0_wp, & ! 5p
      & 4.666137923e-1_wp, 6.644706516e-1_wp, & ! 3d
      & 4.932764167e-1_wp, 5.918727866e-1_wp, & ! 4d
      & 6.539405185e-1_wp, 3.948945302e-1_wp, & ! 5d
      & 4.769346276e-1_wp, 6.587383976e-1_wp, & ! 4f
      & 4.856637346e-1_wp, 6.125980914e-1_wp, & ! 5f
      & 4.848298074e-1_wp, 6.539381621e-1_wp],& ! 5g
      & shape(pCoeff2))


   !> Exponents from third row Table I-V.
   real(wp), parameter :: pAlpha3(3, nf) = reshape([ &
      & 2.227660584e+0_wp, 4.057711562e-1_wp, 1.098175104e-1_wp, & ! 1s
      & 2.581578398e+0_wp, 1.567622104e-1_wp, 6.018332272e-2_wp, & ! 2s
      & 5.641487709e-1_wp, 6.924421391e-2_wp, 3.269529097e-2_wp, & ! 3s
      & 2.267938753e-1_wp, 4.448178019e-2_wp, 2.195294664e-2_wp, & ! 4s
      & 1.080198458e-1_wp, 4.408119382e-2_wp, 2.610811810e-2_wp, & ! 5s
      & 9.192379002e-1_wp, 2.359194503e-1_wp, 8.009805746e-2_wp, & ! 2p
      & 2.692880368e+0_wp, 1.489359592e-1_wp, 5.739585040e-2_wp, & ! 3p
      & 4.859692220e-1_wp, 7.430216918e-2_wp, 3.653340923e-2_wp, & ! 4p
      & 2.127482317e-1_wp, 4.729648620e-2_wp, 2.604865324e-2_wp, & ! 5p
      & 5.229112225e-1_wp, 1.639595876e-1_wp, 6.386630021e-2_wp, & ! 3d
      & 1.777717219e-1_wp, 8.040647350e-2_wp, 3.949855551e-2_wp, & ! 4d
      & 4.913352950e-1_wp, 7.329090601e-2_wp, 3.594209290e-2_wp, & ! 5d
      & 3.483826963e-1_wp, 1.249380537e-1_wp, 5.349995725e-2_wp, & ! 4f
      & 1.649233885e-1_wp, 7.487066646e-2_wp, 3.735787219e-2_wp, & ! 5f
      & 2.545432122e-1_wp, 1.006544376e-1_wp, 4.624463922e-2_wp],& ! 5g
      & shape(pAlpha3))

   !> Coefficients from third row Table I-V.
   real(wp), parameter :: pCoeff3(3, nf) = reshape([ &
      & 1.543289673e-1_wp, 5.353281423e-1_wp, 4.446345422e-1_wp, & ! 1s
      &-5.994474934e-2_wp, 5.960385398e-1_wp, 4.581786291e-1_wp, & ! 2s
      &-1.782577972e-1_wp, 8.612761663e-1_wp, 2.261841969e-1_wp, & ! 3s
      &-3.349048323e-1_wp, 1.056744667e+0_wp, 1.256661680e-1_wp, & ! 4s
      &-6.617401158e-1_wp, 7.467595004e-1_wp, 7.146490945e-1_wp, & ! 5s
      & 1.623948553e-1_wp, 5.661708862e-1_wp, 4.223071752e-1_wp, & ! 2p
      &-1.061945788e-2_wp, 5.218564264e-1_wp, 5.450015143e-1_wp, & ! 3p
      &-6.147823411e-2_wp, 6.604172234e-1_wp, 3.932639495e-1_wp, & ! 4p
      &-1.389529695e-1_wp, 8.076691064e-1_wp, 2.726029342e-1_wp, & ! 5p
      & 1.686596060e-1_wp, 5.847984817e-1_wp, 4.056779523e-1_wp, & ! 3d
      & 2.308552718e-1_wp, 6.042409177e-1_wp, 2.595768926e-1_wp, & ! 4d
      &-2.010175008e-2_wp, 5.899370608e-1_wp, 4.658445960e-1_wp, & ! 5d
      & 1.737856685e-1_wp, 5.973380628e-1_wp, 3.929395614e-1_wp, & ! 4f
      & 1.909729355e-1_wp, 6.146060459e-1_wp, 3.059611271e-1_wp, & ! 5f
      & 1.780980905e-1_wp, 6.063757846e-1_wp, 3.828552923e-1_wp],& ! 5g
      & shape(pCoeff3))


   !> Exponents from forth row Table I-V.
   real(wp), parameter :: pAlpha4(4, nf) = reshape([ &
      & 5.216844534e+0_wp, 9.546182760e-1_wp, & ! 1s
      & 2.652034102e-1_wp, 8.801862774e-2_wp, &
      & 1.161525551e+1_wp, 2.000243111e+0_wp, & ! 2s
      & 1.607280687e-1_wp, 6.125744532e-2_wp, &
      & 1.513265591e+0_wp, 4.262497508e-1_wp, & ! 3s
      & 7.643320863e-2_wp, 3.760545063e-2_wp, &
      & 3.242212833e-1_wp, 1.663217177e-1_wp, & ! 4s
      & 5.081097451e-2_wp, 2.829066600e-2_wp, &
      & 8.602284252e-1_wp, 1.189050200e-1_wp, & ! 5s
      & 3.446076176e-2_wp, 1.974798796e-2_wp, &
      & 1.798260992e+0_wp, 4.662622228e-1_wp, & ! 2p
      & 1.643718620e-1_wp, 6.543927065e-2_wp, &
      & 1.853180239e+0_wp, 1.915075719e-1_wp, & ! 3p
      & 8.655487938e-2_wp, 4.184253862e-2_wp, &
      & 1.492607880e+0_wp, 4.327619272e-1_wp, & ! 4p
      & 7.553156064e-2_wp, 3.706272183e-2_wp, &
      & 3.962838833e-1_wp, 1.838858552e-1_wp, & ! 5p
      & 4.943555157e-2_wp, 2.750222273e-2_wp, &
      & 9.185846715e-1_wp, 2.920461109e-1_wp, & ! 3d
      & 1.187568890e-1_wp, 5.286755896e-2_wp, &
      & 1.995825422e+0_wp, 1.823461280e-1_wp, & ! 4d
      & 8.197240896e-2_wp, 4.000634951e-2_wp, &
      & 4.230617826e-1_wp, 8.293863702e-2_wp, & ! 5d
      & 4.590326388e-2_wp, 2.628744797e-2_wp, &
      & 5.691670217e-1_wp, 2.074585819e-1_wp, & ! 4f
      & 9.298346885e-2_wp, 4.473508853e-2_wp, &
      & 2.017831152e-1_wp, 1.001952178e-1_wp, & ! 5f
      & 5.447006630e-2_wp, 3.037569283e-2_wp, &
      & 3.945205573e-1_wp, 1.588100623e-1_wp, & ! 5g
      & 7.646521729e-2_wp, 3.898703611e-2_wp],&
      & shape(pAlpha4))

   !> Coefficients from forth row Table I-V.
   real(wp), parameter :: pCoeff4(4, nf) = reshape([ &
      & 5.675242080e-2_wp, 2.601413550e-1_wp, & ! 1s
      & 5.328461143e-1_wp, 2.916254405e-1_wp, &
      &-1.198411747e-2_wp,-5.472052539e-2_wp, & ! 2s
      & 5.805587176e-1_wp, 4.770079976e-1_wp, &
      &-3.295496352e-2_wp,-1.724516959e-1_wp, & ! 3s
      & 7.518511194e-1_wp, 3.589627317e-1_wp, &
      &-1.120682822e-1_wp,-2.845426863e-1_wp, & ! 4s
      & 8.909873788e-1_wp, 3.517811205e-1_wp, &
      & 1.103657561e-2_wp,-5.606519023e-1_wp, & ! 5s
      & 1.179429987e+0_wp, 1.734974376e-1_wp, &
      & 5.713170255e-2_wp, 2.857455515e-1_wp, & ! 2p
      & 5.517873105e-1_wp, 2.632314924e-1_wp, &
      &-1.434249391e-2_wp, 2.755177589e-1_wp, & ! 3p
      & 5.846750879e-1_wp, 2.144986514e-1_wp, &
      &-6.035216774e-3_wp,-6.013310874e-2_wp, & ! 4p
      & 6.451518200e-1_wp, 4.117923820e-1_wp, &
      &-1.801459207e-2_wp,-1.360777372e-1_wp, & ! 5p
      & 7.533973719e-1_wp, 3.409304859e-1_wp, &
      & 5.799057705e-2_wp, 3.045581349e-1_wp, & ! 3d
      & 5.601358038e-1_wp, 2.432423313e-1_wp, &
      &-2.816702620e-3_wp, 2.177095871e-1_wp, & ! 4d
      & 6.058047348e-1_wp, 2.717811257e-1_wp, &
      &-2.421626009e-2_wp, 3.937644956e-1_wp, & ! 5d
      & 5.489520286e-1_wp, 1.190436963e-1_wp, &
      & 5.902730589e-2_wp, 3.191828952e-1_wp, & ! 4f
      & 5.639423893e-1_wp, 2.284796537e-1_wp, &
      & 9.174268830e-2_wp, 4.023496947e-1_wp, & ! 5f
      & 4.937432100e-1_wp, 1.254001522e-1_wp, &
      & 6.010484250e-2_wp, 3.309738329e-1_wp, & ! 5g
      & 5.655207585e-1_wp, 2.171122608e-1_wp],&
      & shape(pCoeff4))


   !> Exponents from fifth row Table I-V.
   real(wp), parameter :: pAlpha5(5, nf) = reshape([ &
      & 1.130563696e+1_wp, 2.071728178e+0_wp, 5.786484833e-1_wp, & ! 1s
      & 1.975724573e-1_wp, 7.445271746e-2_wp, &
      & 8.984956862e+0_wp, 1.673710636e+0_wp, 1.944726668e-1_wp, & ! 2s
      & 8.806345634e-2_wp, 4.249068522e-2_wp, &
      & 4.275877914e+0_wp, 1.132409433e+0_wp, 4.016256968e-1_wp, & ! 3s
      & 7.732370620e-2_wp, 3.800708627e-2_wp, &
      & 2.980263783e+0_wp, 3.792228833e-1_wp, 1.789717224e-1_wp, & ! 4s
      & 5.002110360e-2_wp, 2.789361681e-2_wp, &
      & 7.403763257e-1_wp, 1.367990863e-1_wp, 9.135301779e-2_wp, & ! 5s
      & 3.726907315e-2_wp, 2.241490836e-2_wp, &
      & 3.320386533e+0_wp, 8.643257633e-1_wp, 3.079819284e-1_wp, & ! 2p
      & 1.273309895e-1_wp, 5.606243164e-2_wp, &
      & 6.466803859e+0_wp, 1.555914802e+0_wp, 1.955925255e-1_wp, & ! 3p
      & 8.809647701e-2_wp, 4.234835707e-2_wp, &
      & 1.091977298e+0_wp, 3.719985051e-1_wp, 8.590019352e-2_wp, & ! 4p
      & 4.786503860e-2_wp, 2.730479990e-2_wp, &
      & 3.422168934e-1_wp, 1.665099900e-1_wp, 5.443732013e-2_wp, & ! 5p
      & 3.367775277e-2_wp, 2.091949042e-2_wp, &
      & 1.539033958e+0_wp, 4.922090297e-1_wp, 2.029756928e-1_wp, & ! 3d
      & 9.424112917e-2_wp, 4.569058269e-2_wp, &
      & 1.522122079e+0_wp, 2.173041823e-1_wp, 1.084876577e-1_wp, & ! 4d
      & 5.836797641e-2_wp, 3.206682246e-2_wp, &
      & 9.702946470e-1_wp, 3.603270196e-1_wp, 8.668717752e-2_wp, & ! 5d
      & 4.833708379e-2_wp, 2.751899341e-2_wp, &
      & 8.925960415e-1_wp, 3.277589120e-1_wp, 1.492869962e-1_wp, & ! 4f
      & 7.506099109e-2_wp, 3.892475795e-2_wp, &
      & 1.670735676e+0_wp, 2.072477219e-1_wp, 1.024709357e-1_wp, & ! 5f
      & 5.537913898e-2_wp, 3.072866652e-2_wp, &
      & 5.895429375e-1_wp, 2.393343780e-1_wp, 1.172646904e-1_wp, & ! 5g
      & 6.254074479e-2_wp, 3.411243214e-2_wp],&
      & shape(pAlpha5))

   !> Coefficients from fifth row Table I-V.
   real(wp), parameter :: pCoeff5(5, nf) = reshape([ &
      & 2.214055312e-2_wp, 1.135411520e-1_wp, 3.318161484e-1_wp, & ! 1s
      & 4.825700713e-1_wp, 1.935721966e-1_wp, &
      &-1.596349096e-2_wp,-5.685884883e-2_wp, 3.698265599e-1_wp, & ! 2s
      & 5.480512593e-1_wp, 1.472634893e-1_wp, &
      &-3.920358850e-3_wp,-4.168430506e-2_wp,-1.637440990e-1_wp, & ! 3s
      & 7.419373723e-1_wp, 3.724364929e-1_wp, &
      & 1.513948997e-3_wp,-7.316801518e-2_wp,-3.143703799e-1_wp, & ! 4s
      & 9.032615169e-1_wp, 3.294210848e-1_wp, &
      & 1.375523371e-2_wp,-3.097344179e-1_wp,-3.199192259e-1_wp, & ! 5s
      & 1.084547038e+0_wp, 3.345288361e-1_wp, &
      & 2.079051117e-2_wp, 1.235472099e-1_wp, 3.667738886e-1_wp, & ! 2p
      & 4.834930290e-1_wp, 1.653444074e-1_wp, &
      &-2.329023747e-3_wp,-1.357395221e-2_wp, 2.632185383e-1_wp, & ! 3p
      & 5.880427024e-1_wp, 2.242794445e-1_wp, &
      &-1.143929558e-2_wp,-6.322651538e-2_wp, 4.398907721e-1_wp, & ! 4p
      & 5.245859166e-1_wp, 1.017072253e-1_wp, &
      &-3.113958289e-2_wp,-1.374007017e-1_wp, 5.573881018e-1_wp, & ! 5p
      & 4.855428100e-1_wp, 6.605423564e-2_wp, &
      & 2.020869128e-2_wp, 1.321157923e-1_wp, 3.911240346e-1_wp, & ! 3d
      & 4.779609701e-1_wp, 1.463662294e-1_wp, &
      &-3.673711876e-3_wp, 1.167122499e-1_wp, 4.216476416e-1_wp, & ! 4d
      & 4.547673415e-1_wp, 1.037803318e-1_wp, &
      &-3.231527611e-3_wp,-2.434931372e-2_wp, 3.440817054e-1_wp, & ! 5d
      & 5.693674376e-1_wp, 1.511340183e-1_wp, &
      & 1.999839052e-2_wp, 1.395427440e-1_wp, 4.091508237e-1_wp, & ! 4f
      & 4.708252119e-1_wp, 1.328082566e-1_wp, &
      &-7.301193568e-4_wp, 8.414991343e-2_wp, 3.923683153e-1_wp, & ! 5f
      & 5.040033146e-1_wp, 1.328979300e-1_wp, &
      & 1.998085812e-2_wp, 1.460384050e-1_wp, 4.230565459e-1_wp, & ! 5g
      & 4.635699665e-1_wp, 1.226411691e-1_wp],&
      & shape(pCoeff5))


   !> Exponents from sixth row Table I-V.
   real(wp), parameter :: pAlpha6(6, nf) = reshape([ &
      & 2.310303149e+1_wp, 4.235915534e+0_wp, 1.185056519e+0_wp, & ! 1s
      & 4.070988982e-1_wp, 1.580884151e-1_wp, 6.510953954e-2_wp, &
      & 2.768496241e+1_wp, 5.077140627e+0_wp, 1.426786050e+0_wp, & ! 2s
      & 2.040335729e-1_wp, 9.260298399e-2_wp, 4.416183978e-2_wp, &
      & 3.273031938e+0_wp, 9.200611311e-1_wp, 3.593349765e-1_wp, & ! 3s
      & 8.636686991e-2_wp, 4.797373812e-2_wp, 2.724741144e-2_wp, &
      & 3.232838646e+0_wp, 3.605788802e-1_wp, 1.717905487e-1_wp, & ! 4s
      & 5.277666487e-2_wp, 3.163400284e-2_wp, 1.874093091e-2_wp, &
      !& 1.365346e+00_wp,   4.393213e-01_wp,   1.877069e-01_wp,   & ! 4s (old)
      !& 9.360270e-02_wp,   5.052263e-02_wp,   2.809354e-02_wp,   &
      & 1.410128298e+0_wp, 5.077878915e-1_wp, 1.847926858e-1_wp, & ! 5s
      & 1.061070594e-1_wp, 3.669584901e-2_wp, 2.213558430e-2_wp, &
      & 5.868285913e+0_wp, 1.530329631e+0_wp, 5.475665231e-1_wp, & ! 2p
      & 2.288932733e-1_wp, 1.046655969e-1_wp, 4.948220127e-2_wp, &
      & 5.077973607e+0_wp, 1.340786940e+0_wp, 2.248434849e-1_wp, & ! 3p
      & 1.131741848e-1_wp, 6.076408893e-2_wp, 3.315424265e-2_wp, &
      & 2.389722618e+0_wp, 7.960947826e-1_wp, 3.415541380e-1_wp, & ! 4p
      & 8.847434525e-2_wp, 4.958248334e-2_wp, 2.816929784e-2_wp, &
      !& 1.365346e+00_wp,   4.393213e-01_wp,   1.877069e-01_wp,   & ! 4p (old)
      !& 9.360270e-02_wp,   5.052263e-02_wp,   2.809354e-02_wp,   &
      & 3.778623374e+0_wp, 3.499121109e-1_wp, 1.683175469e-1_wp, & ! 5p
      & 5.404070736e-2_wp, 3.328911801e-2_wp, 2.063815019e-2_wp, &
      & 2.488296923e+0_wp, 7.981487853e-1_wp, 3.311327490e-1_wp, & ! 3d
      & 1.559114463e-1_wp, 7.877734732e-2_wp, 4.058484363e-2_wp, &
      & 4.634239420e+0_wp, 1.341648295e+0_wp, 2.209593028e-1_wp, & ! 4d
      & 1.101467943e-1_wp, 5.904190370e-2_wp, 3.232628887e-2_wp, &
      & 8.820520428e-1_wp, 3.410838409e-1_wp, 9.204308840e-2_wp, & ! 5d
      & 5.472831774e-2_wp, 3.391202830e-2_wp, 2.108227374e-2_wp, &
      & 1.357718039e+0_wp, 5.004907278e-1_wp, 2.296565064e-1_wp, & ! 4f
      & 1.173146814e-1_wp, 6.350097171e-2_wp, 3.474556673e-2_wp, &
      & 1.334096840e+0_wp, 2.372312347e-1_wp, 1.269485744e-1_wp, & ! 5f
      & 7.290318381e-2_wp, 4.351355997e-2_wp, 2.598071843e-2_wp, &
      & 8.574668996e-1_wp, 3.497184772e-1_wp, 1.727917060e-1_wp, & ! 5g
      & 9.373643151e-2_wp, 5.340032759e-2_wp, 3.057364464e-2_wp],&
      & shape(pAlpha6))

   !> Coefficients from sixth row Table I-V.
   real(wp), parameter :: pCoeff6(6, nf) = reshape([ &
      & 9.163596280e-3_wp, 4.936149294e-2_wp, 1.685383049e-1_wp, & ! 1s
      & 3.705627997e-1_wp, 4.164915298e-1_wp, 1.303340841e-1_wp, &
      &-4.151277819e-3_wp,-2.067024148e-2_wp,-5.150303337e-2_wp, & ! 2s
      & 3.346271174e-1_wp, 5.621061301e-1_wp, 1.712994697e-1_wp, &
      &-6.775596947e-3_wp,-5.639325779e-2_wp,-1.587856086e-1_wp, & ! 3s
      & 5.534527651e-1_wp, 5.015351020e-1_wp, 7.223633674e-2_wp, &
      & 1.374817488e-3_wp,-8.666390043e-2_wp,-3.130627309e-1_wp, & ! 4s
      & 7.812787397e-1_wp, 4.389247988e-1_wp, 2.487178756e-2_wp, &
      !& 3.775056e-03_wp,  -5.585965e-02_wp,  -3.192946e-01_wp,   & ! 4s (old)
      !&-2.764780e-02_wp,   9.049199e-01_wp,   3.406258e-01_wp,   &
      & 2.695439582e-3_wp, 1.850157487e-2_wp,-9.588628125e-2_wp, & ! 5s
      &-5.200673560e-1_wp, 1.087619490e+0_wp, 3.103964343e-1_wp, &
      & 7.924233646e-3_wp, 5.144104825e-2_wp, 1.898400060e-1_wp, & ! 2p
      & 4.049863191e-1_wp, 4.012362861e-1_wp, 1.051855189e-1_wp, &
      &-3.329929840e-3_wp,-1.419488340e-2_wp, 1.639395770e-1_wp, & ! 3p
      & 4.485358256e-1_wp, 3.908813050e-1_wp, 7.411456232e-2_wp, &
      &-1.665913575e-3_wp,-1.657464971e-2_wp,-5.958513378e-2_wp, & ! 4p
      & 4.053115554e-1_wp, 5.433958189e-1_wp, 1.204970491e-1_wp, &
      !&-7.052075e-03_wp,  -5.259505e-02_wp,  -3.773450e-02_wp,   & ! 4p (old)
      !& 3.874773e-01_wp,   5.791672e-01_wp,   1.221817e-01_wp,   &
      & 1.163246387e-4_wp,-2.920771322e-2_wp,-1.381051233e-1_wp, & ! 5p
      & 5.706134877e-1_wp, 4.768808140e-1_wp, 6.021665516e-2_wp, &
      & 7.283828112e-3_wp, 5.386799363e-2_wp, 2.072139149e-1_wp, & ! 3d
      & 4.266269092e-1_wp, 3.843100204e-1_wp, 8.902827546e-2_wp, &
      &-4.749842876e-4_wp,-3.566777891e-3_wp, 1.108670481e-1_wp, & ! 4d
      & 4.159646930e-1_wp, 4.621672517e-1_wp, 1.081250196e-1_wp, &
      &-4.097377019e-3_wp,-2.508271857e-2_wp, 2.648458555e-1_wp, & ! 5d
      & 5.097437054e-1_wp, 2.654483467e-1_wp, 2.623132212e-2_wp, &
      & 6.930234381e-3_wp, 5.634263745e-2_wp, 2.217065797e-1_wp, & ! 4f
      & 4.411388883e-1_wp, 3.688112625e-1_wp, 7.787514504e-2_wp, &
      &-9.486751531e-4_wp, 4.624275998e-2_wp, 2.373699784e-1_wp, & ! 5f
      & 4.589112231e-1_wp, 3.205010548e-1_wp, 5.077063693e-2_wp, &
      & 6.729778096e-3_wp, 5.874145170e-2_wp, 2.339955227e-1_wp, & ! 5g
      & 4.512983737e-1_wp, 3.552053926e-1_wp, 6.974153145e-2_wp],&
      & shape(pCoeff6))

   real(wp), parameter :: pAlpha6s(6) = [&
      & 5.800292686e-1_wp, 2.718262251e-1_wp, 7.938523262e-2_wp, & ! 6s
      & 4.975088254e-2_wp, 2.983643556e-2_wp, 1.886067216e-2_wp]
   real(wp), parameter :: pCoeff6s(6) = [&
      & 4.554359511e-3_wp, 5.286443143e-2_wp,-7.561016358e-1_wp, & ! 6s
      &-2.269803820e-1_wp, 1.332494651e+0_wp, 3.622518293e-1_wp]

   real(wp), parameter :: pAlpha6p(6) = [&
      & 6.696537714e-1_wp, 1.395089793e-1_wp, 8.163894960e-2_wp, & ! 6p
      & 4.586329272e-2_wp, 2.961305556e-2_wp, 1.882221321e-2_wp]
   real(wp), parameter :: pCoeff6p(6) = [&
      & 2.782723680e-3_wp,-1.282887780e-1_wp,-2.266255943e-1_wp, & ! 6p
      & 4.682259383e-1_wp, 6.752048848e-1_wp, 1.091534212e-1_wp]


contains


!> Expand Slater function in primitive gaussian functions
pure subroutine slater_to_gauss_cgto(ng, n, l, zeta, cgto, norm, info)
   !> Number of Gaussian functions for the expansion
   integer, intent(in) :: ng
   !> Principal quantum number of shell
   integer, intent(in) :: n
   !> Azimudal quantum number of shell
   integer, intent(in) :: l
   !> Exponent of Slater function to expand
   real(wp), intent(in) :: zeta
   !> Contracted Gaussian basis function
   type(cgto_type), intent(out) :: cgto
   !> Include normalization in contraction coefficients
   logical, intent(in) :: norm
   !> Status of the expansion, returns zero for success or the position of the
   !> faulty dummy argument
   integer, intent(out) :: info

   cgto%ang = l
   cgto%nprim = ng
   call slater_to_gauss(ng, n, l, zeta, cgto%alpha, cgto%coeff, &
      & norm, info)

end subroutine slater_to_gauss_cgto


!> Expand Slater function in primitive gaussian functions
pure subroutine slater_to_gauss_array(ng, n, l, zeta, alpha, coeff, norm, info)
   !> Number of Gaussian functions for the expansion
   integer, intent(in) :: ng
   !> Principal quantum number of shell
   integer, intent(in) :: n
   !> Azimudal quantum number of shell
   integer, intent(in) :: l
   !> Exponent of Slater function to expand
   real(wp), intent(in) :: zeta
   !> Exponents of primitive gaussian functions
   real(wp), intent(out) :: alpha(:)
   !> Contraction coefficients of primitive gaussians, can contain normalization
   real(wp), intent(out) :: coeff(:)
   !> Include normalization in contraction coefficients
   logical, intent(in) :: norm
   !> Status of the expansion, returns zero for success or the position of the
   !> faulty dummy argument
   integer, intent(out) :: info

   ! Storage location of the STO-NG coefficients and exponents
   integer :: ityp

   ! Basic checks first, we cannot support principal QN of six or higher
   ! also you should not violate l ∊ [n-1, n-2, ..., 1, 0]
   if ((n > 5).or.(n <= l)) then
      if (.not.(n == 6 .and. ng == 6)) then
         info = 2
         return
      end if
   endif

   ! Don't allow negative exponents in the first place
   if (zeta <= 0.0_wp) then
      info = 4
      return
   endif

   ! we have to use a little hack here,
   ! if you pass n and l correctly, everything is fine
   ! ityp: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
   !    n: 1 2 3 4 5 2 3 4 5  3  4  5  4  5  5
   !    l: 0 0 0 0 0 1 1 1 1  2  2  2  3  3  4
   select case(l) ! integer hack:
   case(0); ityp = n    ! s
   case(1); ityp = 4+n  ! p
   case(2); ityp = 7+n  ! d
   case(3); ityp = 9+n  ! f
   case(4); ityp = 10+n ! g
   case default ! I'm sorry, no h-functions for you today
      info = 3
      return
   end select

   select case(ng)
   case default ! currently we cannot go beyond 6 primitives
      info = 1
      return
   case(1)
      alpha(1) = pAlpha1(ityp) * zeta**2
      coeff(1) = 1.0_wp
   case(2)
      alpha(:ng) = pAlpha2(:, ityp) * zeta**2
      coeff(:ng) = pCoeff2(:, ityp)
   case(3)
      alpha(:ng) = pAlpha3(:, ityp) * zeta**2
      coeff(:ng) = pCoeff3(:, ityp)
   case(4)
      alpha(:ng) = pAlpha4(:, ityp) * zeta**2
      coeff(:ng) = pCoeff4(:, ityp)
   case(5)
      alpha(:ng) = pAlpha5(:, ityp) * zeta**2
      coeff(:ng) = pCoeff5(:, ityp)
   case(6)
      if (n == 6) then
         if (l == 0) then
            alpha(:ng) = pAlpha6s(:) * zeta**2
            coeff(:ng) = pCoeff6s(:)
         else if (l == 1) then
            alpha(:ng) = pAlpha6p(:) * zeta**2
            coeff(:ng) = pCoeff6p(:)
         else
            info = 2
            return
         end if
      else
         alpha(:ng) = pAlpha6(:, ityp) * zeta**2
         coeff(:ng) = pCoeff6(:, ityp)
      end if
   end select

   ! normalize the gaussian if requested
   ! <φ|φ> = (2i-1)!!(2j-1)!!(2k-1)!!/(4α)^(i+j+k) · sqrt(π/2α)³
   ! N² = (4α)^(i+j+k)/((2i-1)!!(2j-1)!!(2k-1)!!)  · sqrt(2α/π)³
   ! N = (4α)^((i+j+k)/2) / sqrt((2i-1)!!(2j-1)!!(2k-1)!!) · (2α/π)^(3/4)
   if (norm) then
      coeff(:ng) = coeff(:ng) * (top*alpha(:ng))**0.75_wp &
         & * sqrt(4*alpha(:ng))**l / sqrt(dfactorial(l+1))
   endif

   ! success
   info = 0

end subroutine slater_to_gauss_array


end module tblite_basis_slater
