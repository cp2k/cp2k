Orbital-energy-based atomistic properties
=========================================

The orbital-energy-based properties are computed based on the orbital energies of the system.
A localization of the orbitals is performed to obtain the atomic properties.
The following properties are available:

.. table:: Orbital-energy-based properties
   :widths: auto

   =============== ================================= ==================================================================================================== 
   Dictionary Key  Description                       Equation            
   =============== ================================= ====================================================================================================
   ``response``    response function                 :math:`\chi_A = \sum_{i=1}^{n_{\text{occ.}}} \sum_{a=1}^{n_{\text{virt.}}} \chi_A^{ia}`
                                                     :math:`\chi_A^{ia} = \frac{p_{i,A} p_{a,A}}{\Delta\varepsilon^2_{ai}+\Delta^2} \text{ with: } i \in \text{occ. orbital}, a \in \text{virt. orbital}` 
   ``gap``         atomistic HL-gap in eV            :math:`\epsilon_{A,\text{HL-gap}} = \frac{1}{\sum_{i=1}^{n_{\text{occ.}}} \sum_{a=1}^{n_{\text{virt.}}} \lambda_A^{ia}\frac{1}{\Delta\varepsilon_{ia}+\Delta}}-\Delta \text{ with: } \lambda_A^{ia}=\frac{\chi_A^{ia}}{\sum_{j,b}\chi_A^{jb}}`
   ``chem_pot``    chemical potential of an atom     :math:`\epsilon_{A,F} = \sum_{i=1}^{n_{\text{occ.}}} \sum_{a=1}^{n_{\text{virt.}}} \lambda_A \frac{\varepsilon_a+\varepsilon_i}{2}`
   ``HOAO``        highest occupied atomic orbital   :math:`\epsilon_{A,\text{HOAO}} = \epsilon_{A,F} - \frac{\epsilon_{A,\text{HL-gap}}}{2}`
   ``LUAO``        lowest unoccupied atomic orbital  :math:`\epsilon_{A,\text{LUAO}} = \epsilon_{A,F} + \frac{\epsilon_{A,\text{HL-gap}}}{2}`
   =============== ================================= ====================================================================================================

where :math:`\Delta = 0.5 \text{ eV}` is a damping factor to avoid division by zero.

.. note:: 
   For calculations with unpaired electrons, the properties are available for the alpha and beta spin channel separately.
   Without spin polarization, ``tblite`` does not provide two sets of density matrices and occupation vectors.
   In that case, we define occupation vectors that match the number of unpaired electrons using integer occupations.
   When using spin-polarized Hamiltonians, we use the occupation vectors and orbital energies from solving the UHF equations.
   In both cases the properties are computed for the alpha and beta spin channel separately.
   In this case all properties are available with the suffix ``_alpha`` and ``_beta`` for the alpha and beta spin channel, respectively.

Extended orbital-energy-based atomistic properties
==================================================

The extended orbital-energy-based properties are computed based on the atomic properties seen above.
The convolution kernel :math:`f_{\text{log}}(R_A,R_B)` has been defined in the :doc:`xtbml <index>` section.
The following properties are available:

.. table:: Extended orbital-energy-based properties
   :widths: auto

   ================ ===================================== =====================================================================================================
   Dictionary Key   Description                           Equation            
   ================ ===================================== =====================================================================================================
   ``ext_gap``      ext. atomistic HL-gap in eV            :math:`\epsilon_{A,\text{HL-gap, ext}} = \sum_B \epsilon_{B,\text{HL-gap}} \cdot \beta(A,B)` with: :math:`\beta(A,B) = \frac{1}{f_{\text{log}}(R_A,R_B) \cdot(\text{CN}_A+1)}`
   ``ext_chem_pot`` ext. chemical potential of an atom     :math:`\epsilon_{A,F,\text{ext}} = \sum_{B} \mu_{B,F} \cdot \beta(A,B)`
   ``ext_HOAO``     ext. highest occupied atomic orbital   :math:`\epsilon_{A,\text{HOAO}} = \epsilon_{A,F,\text{ext}} - \frac{\epsilon_{A,\text{HL-gap, ext}}}{2}`
   ``ext_LUAO``     ext. lowest unoccupied atomic orbital  :math:`\epsilon_{A,\text{LUAO}} = \epsilon_{A,F,\text{ext}} + \frac{\epsilon_{A,\text{HL-gap, ext}}}{2}`
   ================ ===================================== =====================================================================================================