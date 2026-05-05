Density-based atomistic properties
==================================

The density-based properties are computed based on the electron density of the system.
We mainly focus on the cumulative multipole moments and the Mulliken population analysis.
The following properties are available:

.. table:: Density-based properties
   :widths: auto

   =============== =============================== ==================================================================================================== 
   Dictionary Key  Description                     Equation            
   =============== =============================== ====================================================================================================
   ``p_s``         Mulliken population s-shell     :math:`p_s = \sum_{\mu \in s \in A}^{N_{\text{BF}}}\sum_{\nu}^{N_{\text{BF}}}P_{\mu \nu} S_{\mu\nu}`
   ``p_p``         Mulliken population p-shell     :math:`p_p = \sum_{\mu \in p \in A}^{N_{\text{BF}}}\sum_{\nu}^{N_{\text{BF}}}P_{\mu \nu} S_{\mu\nu}`
   ``p_d``         Mulliken population d-shell     :math:`p_d = \sum_{\mu \in d \in A}^{N_{\text{BF}}}\sum_{\nu}^{N_{\text{BF}}}P_{\mu \nu} S_{\mu\nu}`   
   ``dipm_s``      dipole moment s-shell           :math:`\mu_s = |\boldsymbol{\mu}_s| =\sqrt{\sum_{\alpha \in \:(x,y,z)}\left( \sum_{\kappa\in s \in A}\sum_{\lambda} P_{\kappa\lambda} \left( \alpha_A S_{\lambda\kappa} - D^\alpha_{\lambda\kappa}\right) \right)^2}`
                                                   |br|  with the dipole moment integral: :math:`D^\alpha_{\lambda\kappa} = \left\langle \psi_\lambda |\alpha| \psi_\kappa\right\rangle` 
   ``dipm_p``      dipole moment p-shell           :math:`\mu_p = |\boldsymbol{\mu}_p| =\sqrt{\sum_{\alpha \in \:(x,y,z)}\left( \sum_{\kappa\in p \in A}\sum_{\lambda} P_{\kappa\lambda} \left( \alpha_A S_{\lambda\kappa} - D^\alpha_{\lambda\kappa}\right) \right)^2}` 
   ``dipm_d``      dipole moment d-shell           :math:`\mu_d = |\boldsymbol{\mu}_d| =\sqrt{\sum_{\alpha \in \:(x,y,z)}\left( \sum_{\kappa\in d \in A}\sum_{\lambda} P_{\kappa\lambda} \left( \alpha_A S_{\lambda\kappa} - D^\alpha_{\lambda\kappa}\right) \right)^2}` 
   ``qm_s``        quadrupole moment s-shell       :math:`\Theta_s = ||\boldsymbol{\Theta}_s|| = \sqrt{ \sum_{\alpha\beta \:\in\:(xx,xy,xz,yy,yz,zz)}2\left(\Theta^{\alpha\beta}_{s}-\delta_{\alpha\beta}\Theta^{\alpha\beta}_{l}\right)^2 }` a_ b_
   ``qm_p``        quadrupole moment p-shell       :math:`\Theta_p = ||\boldsymbol{\Theta}_p|| = \sqrt{ \sum_{\alpha\beta \:\in\:(xx,xy,xz,yy,yz,zz)}2\left(\Theta^{\alpha\beta}_{p}-\delta_{\alpha\beta}\Theta^{\alpha\beta}_{l}\right)^2 }` a_ b_
   ``qm_d``        quadrupole moment d-shell       :math:`\Theta_d = ||\boldsymbol{\Theta}_d|| = \sqrt{ \sum_{\alpha\beta \:\in\:(xx,xy,xz,yy,yz,zz)}2\left(\Theta^{\alpha\beta}_{d}-\delta_{\alpha\beta}\Theta^{\alpha\beta}_{l}\right)^2 }` a_ b_
   ``q_A``         atomic partial charge           :math:`q_A = Z'_A - \sum_{\mu\in A}^{N_{\text{BF}}} \sum_{\nu}^{N_{\text{BF}}} P_{\mu \nu} S_{\mu\nu}`
   ``dipm_A``      atomic dipole moment            :math:`\mu_A = |\boldsymbol{\mu}_A| = \sqrt{\sum_{\alpha \in \:(x,y,z)}\left(\mu_A^\alpha\right)^2} \text{ with: } \mu_A^\alpha = \sum_{l\in A} \mu_l^\alpha` with: :math:`l \in {s,p,d}`
   ``qm_A``        atomic quadrupole moment        :math:`\Theta_A = ||\boldsymbol{\Theta}_A||=\sqrt{ \sum_{\alpha\beta \:\in\:(xx,xy,xz,yy,yz,zz)}2\left(\Theta^{\alpha\beta}_{A}-\delta_{\alpha\beta}\Theta^{\alpha\beta}_{A}\right)^2 }` with: :math:`\Theta_A^{\alpha\beta} = \sum_{l\in A} \Theta_l^{\alpha\beta}`                
   =============== =============================== ====================================================================================================
   
.. _a:

a. with :math:`\Theta_l^{\alpha\beta} = \frac{3}{2} \theta_l^{\alpha\beta} - \frac{\delta_{\alpha\beta}}{2}\left( \theta_l^{xx}+ \theta_l^{yy}+\theta_l^{zz}\right)`
|br| and :math:`\theta_l^{\alpha\beta} = \sum_{\kappa\in l\in A}\sum_{\lambda} P_{\kappa\lambda} \left( \alpha_A D^\beta_{\lambda\kappa}+\beta_A D^\alpha_{\lambda\kappa} -\alpha_A\beta_A S_{\lambda\kappa} - Q^{\alpha\beta}_{\lambda\kappa}\right)`
|br| and the quadrupole moment integral: :math:`Q^{\alpha\beta}_{\lambda\kappa} = \left\langle \psi_\lambda |\alpha\beta| \psi_\kappa\right\rangle`

.. _b:

b. Kronecker delta: :math:`\delta_{\alpha\beta} = 1` if :math:`\alpha = \beta` and 0 otherwise.

.. note:: 
   For unrestricted calculations, the properties are computed for each spin channel separately.
   This is currently only the case, if a spin-polarized calculation has been performed.
   In this case all properties are available with the suffix ``_alpha`` and ``_beta`` for the alpha and beta spin channel, respectively.


Extended density-based properties
=================================

.. |br| raw:: html

  <br/>

The extended density-based properties are computed based on the atomic properties seen above.
The convolution kernel :math:`f_{\text{log}}(R_A,R_B)` has been defined in the :doc:`xtbml <index>` section.
The following properties are available:

.. table:: Density-based extended properties
   :widths: auto

   =============== =============================== ==================================================================================================== 
   Dictionary Key  Description                     Equation            
   =============== =============================== ====================================================================================================
   ``ext_q_A``     ext. atomic partial charge      :math:`q_{A,\text{ext}} = \frac{q_A}{(CN_A+1)} + \sum_{B\neq A}^{N_{\text{atoms}}} q_B \cdot \frac{f_{\text{log}(R_A,R_B)}}{(CN_B+1)}`
   ``ext_dipm_A``  ext. atomic dipole moment       :math:`\mu_{A,\text{ext}} = | \boldsymbol{\mu}_{A,\text{ext}} | = \sqrt{\sum_{\alpha \in \:(x,y,z)} \left( \mu_{A,\text{ext}}^\alpha \right)^2}` 
                                                   with: :math:`\mu_{A,\text{ext}}^\alpha = \frac{\mu_A^\alpha}{(CN_A+1)} + \sum_{B\neq A} \left( \mu_B^\alpha - \alpha_{AB} q_B\right) \frac{f_{\text{log}(R_A,R_B)}}{(CN_B+1)}`
   ``ext_qm_A``    ext. atomic quadrupole moment   :math:`\Theta_{A,\text{ext}}= ||\boldsymbol{\Theta}_{A,\text{ext}}||= 	\sqrt{ \sum_{\alpha\beta \:\in\:(xx,xy,xz,yy,yz,zz)}2\left(\Theta^{\alpha\beta}_{A,\text{ext}}-\delta_{\alpha\beta}\Theta^{\alpha\beta}_{A,ext}\right)^2 }`             
                                                   |br| with: :math:`\Theta^{\alpha\beta}_{A,\text{ext}}= \frac{\Theta^{\alpha\beta}_{A}} {(CN_A+1)}+\sum_{B\neq A} \left(\Theta_B^{\alpha\beta}+\delta\Theta_B^{\alpha\beta}\right) \frac{f_{\text{log}}}{(CN_B+1)}`
                                                   |br| where: :math:`\delta\Theta_B^{\alpha\beta} = \frac{3}{2} \delta\theta_B^{\alpha\beta} - \frac{\delta_{\alpha\beta}}{2}\left( \delta\theta_B^{xx}+ \delta\theta_B^{yy}+\delta\theta_B^{zz}\right)`
                                                   |br| and: :math:`\delta\theta_B^{\alpha\beta} = -\beta_{AB}\mu^\alpha_B + \alpha_{AB}\mu^\beta_B+\alpha_{AB}\beta_{AB}q_B \text{ where: } \beta_{AB} = \beta_{A}-\beta_B`
   ``ext_dipm_e``  ext. atomic dipole              :math:`\mu_{A,\text{ext},e} =|\boldsymbol{\mu}_{A,\text{ext},e}| = \sqrt{\sum_{\alpha \in \:(x,y,z)}\left(\mu_{A,\text{ext},e}^\alpha\right)^2}`
                   (only electronic effects)       |br| with: :math:`\mu_{A,\text{ext},e}^\alpha = \frac{\mu_A^\alpha} {(CN_A+1)}+ \sum_{B\neq A} \left(\mu_B^\alpha + \alpha_{AB} p_B\right)  \frac{f_{\text{log}}}{(CN_B+1)}`
   ``ext_qm_e``    ext. atomic quadrupole          :math:`\Theta_{A,\text{ext},e}= ||\boldsymbol{\Theta}_{A,\text{ext},e}||= 	\sqrt{ \sum_{\alpha\beta \:\in\:(xx,xy,xz,yy,yz,zz)}2\left(\Theta^{\alpha\beta}_{A,\text{ext},e}\right)^2-\delta_{\alpha\beta}\left(\Theta^{\alpha\beta}_{A,\text{ext},e}\right)^2 }`
                   (only electronic effects)       |br| with: :math:`\Theta^{\alpha\beta}_{A,\text{ext},e} =  \sum_{B\neq A} \left(\Theta^{\alpha\beta}_{A} +\Theta^{\alpha\beta}_{B}+\delta\Theta_{B,e}^{\alpha\beta}\right) \frac{f_{\text{log}}}{(CN_B+1)}`
                                                   |br| with: :math:`\delta\Theta_{B,e}^{\alpha\beta} = \frac{3}{2} \delta\theta_{B,e}^{\alpha\beta} - \frac{\delta_{\alpha\beta}}{2}\left( \delta\theta_{B,e}^{xx}+ \delta\theta_{B,e}^{yy}+\delta\theta_{B,e}^{zz}\right)`
                                                   |br| with: :math:`\delta\theta_{B,e}^{\alpha\beta} = -\beta_{AB}\mu^\alpha_B + \alpha_{AB}\mu^\beta_B-\alpha_{AB}\beta_{AB}p_B`
   ``ext_dipm_Z``  ext. atomic dipole              :math:`\mu_{A,\text{ext},Z} =|\boldsymbol{\mu}_{A,\text{ext},Z}| = \sqrt{\sum_{\alpha \in \:(x,y,z)}\left(\mu_{A,\text{ext},Z}^\alpha\right)^2}`
                   (only nuclear effects)          |br| with: :math:`\mu_{A,\text{ext},Z}^\alpha = \sum_{B\neq A} \left(- \alpha_{AB} Z'_B \right) \frac{f_{\text{log}}}{(CN_B+1)}`
   ``ext_qm_Z``    ext. atomic quadrupole          :math:`\Theta_{A,\text{ext},Z}= ||\boldsymbol{\Theta}_{A,\text{ext},Z}||= \sqrt{ \sum_{\alpha\beta \:\in\:(xx,xy,xz,yy,yz,zz)}2\left(\Theta^{\alpha\beta}_{A,\text{ext},Z}\right)^2-\delta_{\alpha\beta}\left(\Theta^{\alpha\beta}_{A,\text{ext},Z}\right)^2 }`
                   (only nuclear effects)          |br| with: :math:`\Theta^{\alpha\beta}_{A,\text{ext},Z} = \sum_{B\neq A}  \left( \delta\Theta_{B,Z}^{\alpha\beta} \right)  \frac{f_{\text{log}}}{(CN_B+1)}`
                                                   |br| with: :math:`\delta\Theta_{B,Z}^{\alpha\beta} = \frac{3}{2} \delta\theta_{B,Z}^{\alpha\beta} - \frac{\delta_{\alpha\beta}}{2}\left( \delta\theta_{B,Z}^{xx}+ \delta\theta_{B,Z}^{yy}+\delta\theta_{B,Z}^{zz}\right)`
                                                   |br| with: :math:`\delta\theta_{B,Z}^{\alpha\beta} =\alpha_{AB}\beta_{AB}Z'_B`                     
   =============== =============================== ====================================================================================================

.. note:: 
   For unrestricted calculations, the properties are computed for each spin channel separately.
   This is currently only the case, if a spin-polarized calculation has been performed.
   In this case all properties are available with the suffix ``_alpha`` and ``_beta`` for the alpha and beta spin channel, respectively.