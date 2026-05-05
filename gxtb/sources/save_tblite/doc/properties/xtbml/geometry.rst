Geometry-based atomistic properties
===================================

The geometry-based properties are based solely on geometry information.
The following property is computed:

.. table:: Geometry-based properties
   :widths: auto

   =============== =============================== ==================================================================================================== 
   Dictionary Key  Description                     Equation            
   =============== =============================== ====================================================================================================
   ``CN_A``           D3-Coordination number         :math:`CN_A = \sum_{B\neq A}^{N_{\text{atoms}}}\frac{1}{1+\exp\left(\frac{-16a\cdot4(R_{A,\text{cov}}+R_{B,\text{cov}})}{3R_{AB}}-1\right)}`
                                                   where: :math:`R_{AB} = \sqrt{(x_a-x_b)^2+(y_a-y_b)^2+(z_a-z_b)^2}`
   =============== =============================== ====================================================================================================

Extended geometry-based properties
==================================

The extended geometry-based properties are computed based on the atomic properties seen above.
The convolution kernel :math:`f_{\text{log}}(R_A,R_B)` has been defined in the :doc:`xtbml <index>` section.
The following property is available:

.. table:: Geometry-based extended properties
   :widths: auto

   =============== =============================== ==================================================================================================== 
   Dictionary Key  Description                     Equation            
   =============== =============================== ====================================================================================================
   ``ext_CN_A``      ext. D3-Coordination number    :math:`CN_{A,\text{ext}} = \sum_{B\neq A}^{N_{\text{atoms}}} CN_B \cdot f_{\text{log}(R_A,R_B)}`
   =============== =============================== ====================================================================================================