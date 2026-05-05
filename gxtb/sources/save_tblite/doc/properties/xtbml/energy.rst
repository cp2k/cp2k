Energy-based atomistic properties
=================================

The energy-based properties are based on the energy terms of the xTB Hamiltonian.
An atomistic partitioning is used to compute the properties.
By exposing all terms of the energy term, a ML model could combine those terms to build a better fitting model.
The following energy terms are computed:

.. table:: Energy-based properties
   :widths: auto

   =============== =============================== ==================================================================================================== 
   Dictionary Key  Description                     Equation            
   =============== =============================== ====================================================================================================
   ``E_EHT``       Extended Hueckel term           :math:`E_{A,\text{EHT}}= \sum_{\kappa \in A} \sum_{\kappa} P_{\kappa\lambda} H_{\kappa\lambda}`
   ``E_rep``       Repulsion term                  :math:`E_{A,\text{rep}}=  \sum_{B\neq A}^{M} \frac{1}{2} \frac{Y_A^{\text{eff}} Y_B^{\text{eff}}}{R_{AB}}e^{(\alpha_A\alpha_B)^{0.5}(R_{AB})^{k_{\text{rep}}}}`, see also :ref:`repulsion`
   ``E_ies_ixc``   isotropic electrostatic and     :math:`E_{A,\text{IES+IXC}}=\frac{1}{2}\sum_{B\neq A}\sum_{l\in A}\sum_{l'\in B} q_{A,l} q_{B,l'} \gamma_{AB,ll'} + \frac{1}{3}\sum_{l\in A} \Gamma_{A,l}q_{A,l}^3`, see also :ref:`coulomb`
                   exchange correlation term
   ``E_axc``       anisotropic exchange            :math:`E_{A,\text{AXC}} = \left(f_{\text{XC}}^{\mu_A}|\boldsymbol{\mu}_A|^2+f^{\Theta_A}_{\text{XC}}||\boldsymbol{\Theta}_A||^2\right)`
                   correlation term
   ``E_aes``       anisotropic electrostatic term  see :math:`E_{AES}` definition in :ref:`coulomb`, while only summing up over atom B
   ``E_disp2``     two body dispersion term        :math:`E_{A,\text{disp}}^{(2)} = \sum_{A\neq B}\sum_{n=6,8} - s_n \frac{1}{2} \frac{C_n^{AB}(q_a,CN^A_{\text{cov}},q_b,CN^B_{\text{cov}})}{R^n_{AB}} f^{\text{log,BJ}}_n(R_{AB})`
   ``E_disp3``     three body dispersion term      :math:`E_{A,\text{disp}}^{(3)} = \sum_{B\neq A} \sum_{C\neq B\neq A}-s_9\frac{1}{3} \frac{(3 \cos(\theta_{ABC}) \cos(\theta_{BCA}) \cos(\theta_{CAB})+1)C_9^{ABC} (CN_{\text{cov}}^A,CN_{\text{cov}}^B,CN_{\text{cov}}^C)}{(R_{AB}R_{AC}R_{BC})^3} \cdot f^{\text{log,zero}}_9(R_{AB},R_{AC},R_{BC})`
   ``E_hx``        Halogen bond term               
   ``E_tot``       Sum of energy terms
   ``w_tot``       weight of atom in total energy   :math:`w_{A,\text{tot}} = \frac{E_{A,\text{tot}}}{E_{\text{tot}}}` 
   =============== =============================== ====================================================================================================

Extended energy-based properties
================================

There are **no** extended energy-based properties available.