Spectral Boundary Integral Algorithm
====================================

This section provides a detailed description of the spectral boundary integral algorithm as implemented in uguca. For more details about the derivation of the SBI method, please refer to :doc:`./references`.

Definitions
-----------

  - :math:`u_{n}` : nodal displacement vector at time :math:`n \Delta t`
  - :math:`v_{n}` : nodal velocity vector at time :math:`n \Delta t`
  - :math:`f_{n}` : nodal internal force vector at time :math:`n \Delta t` (due to deformation)
  - :math:`e_{n}` : nodal external load vector at time :math:`n \Delta t`
  - :math:`i_{n}` : nodal interface traction vector at time :math:`n \Delta t`
  - :math:`r_{n}` : nodal residual vector at time :math:`n \Delta t`
  - :math:`c_d` and :math:`c_s` : primary and secondary wave speeds
  - :math:`\mu` : shear modulus of material


Time integration
----------------

The time integration is computed with the following iterative procedure:

.. math::
   e_{n+1} &= \textrm{update external loading}\\
   u_{n+1}^{\pm} &= u_n^{\pm} + \Delta t \cdot v_{n}^{\pm} \\
   f_{n+1}^{\pm} &= \textrm{through convolution etc.}\\
   i_{n+1} &= \textrm{interface algorithm}\\
   r_{n+1}^{\pm} &= \pm (e_{n+1} + f_{n+1}^{\pm}- i_{n+1}) \\
   v_{n+1}^{\pm} &= A \cdot r_{n+1}^{\pm}

where the superscript :math:`\pm` indicates top and bottom interface, and :math:`A` is given for tangential direction by :math:`A^T = c_s / \mu` and normal direction by :math:`A^N = c_s^2 / \mu / c_d`

Comparing to the notation of [Barras et al., 2014], equation (10) provides the relation to compute the residual, with :math:`\tau_j(x,t)` being :math:`i`, :math:`\tau_j^0(x,t)` being :math:`e`, and :math:`f_j(x,t)` being :math:`f`. The velocity term is then used to compute the velocity based on the residual.

Convolution term
----------------

We use the independent spectral formulation for 2D and 3D as introduced by [Breitenfeld and Geubelle, 1998].

2D formulation
^^^^^^^^^^^^^^

In the time domain, the 2D elastodynamic relations between the traction components for the stress :math:`\tau_\alpha` acting on the fracture plane and the resulting displacements :math:`u_\alpha^\pm` take the form 

.. math::

   \tau_0^\pm(x_0,t)&=\tau_0^0(x_0,t)\mp\frac{\mu}{c_s}\frac{\partial{u_0^\pm}}{\partial{t}}+f_0^\pm(x_0,t),\\
   \tau_1^\pm(x_0,t)&=\tau_1^0(x_0,t)\mp\frac{(\lambda+2\mu)}{c_p}\frac{\partial{u_1^\pm}}{\partial{t}}+f_1^\pm(x_0,t),


where :math:`f` is a linear functional of the prior deformation history and can be expressed in terms of its spectral components as :math:`f_j(x_1,t)=F_j(t;q)e^{iqx_0}`. These components are given in [Breitenfeld and Geubelle, 1998] as

.. math::
   F_0^\pm(t;q)= &\mp\mu^\pm |q| \int_{0}^{t} H_{00}(|q|c_s^\pm t') U_0^\pm(t-t';q)|q|c_s^\pm \mathrm{d}t' + i(2-\eta^\pm)\mu^\pm q U_1^\pm(t;q) \\
   &+ i\mu^\pm q \int_0^t H_{01}(|q|c_s^\pm t') U_1^\pm (t-t';q)|q|c_s^\pm \mathrm{d}t',\\[5pt]
   F_1^\pm(t;q)= &\mp\mu^\pm |q|\int_{0}^{t}H_{11}(|q|c_s^\pm t') U_1^\pm(t-t';q)|q|c_s^\pm \mathrm{d}t' - i(2-\eta^\pm)\mu^\pm q U_0^\pm(t;q) \\
   &- i\mu^\pm q \int_0^t H_{01}(|q|c_s^\pm t') U_0^\pm (t-t';q)|q|c_s^\pm dt',\\


where :math:`\eta = c_p/c_s` and :math:`H_{\alpha\beta}` are the convolution kernels.

*Note that we corrected a sign error present in [Breitenfeld and Geubelle, 1998].*


3D formulation
^^^^^^^^^^^^^^

We replace the mode number :math:`q`  by a mode vector :math:`\mathbf{q}=(k,m)` which spans the fracture plane

.. math::
   [u_j^\pm(x_0,x_2,t),\tau_j(x_0,x_2,t)]=[U_j^\pm(t;k,m), T_j(t;k,m)] \exp\left(i(k x_0 + m x_2)\right).


The elastodynamic relations between the traction components for the stress :math:`\tau_j` and the resulting displacements :math:`u_j^\pm` becomes:

.. math::
  \tau_j^\pm(x_0,x_2,t)=\tau_j^0(x_0,x_2,t)\mp V_{jk}^\pm\frac{\partial{u_k^\pm}(x_0,x_1,t)}{\partial{t}}+f_j^\pm(x_0,x_1,t)

where the matrix :math:`[V_{ij}]` is diagonal

.. math::
   V_{00}=V_{22}=\frac{\mu}{c_s}, \quad V_{11}=\eta \frac{\mu}{c_s}


The Fourier coefficients :math:`F_j^\pm(t;k,m)` of the :math:`f_j^\pm(x_0,x_2,t)` involve a series of convolutions, as it was the case in 2D.

.. math::
   \begin{Bmatrix}
   F_0^\pm(t;k,m)\\[5pt]
   F_1^\pm(t;k,m)\\[5pt]
   F_2^\pm(t;k,m)
   \end{Bmatrix}\\[10pt]
   = & -i \mu^\pm(2-\eta^\pm)
   \begin{pmatrix}
   0 & -k & 0 \\[5pt]
   k &  0 & m \\[5pt]
   0 & -m & 0 
   \end{pmatrix}
   \begin{Bmatrix}
   U_0^\pm(t;k,m)\\[5pt]
   U_1^\pm(t;k,m)\\[5pt]
   U_2^\pm(t;k,m)
   \end{Bmatrix}\\[10pt]
   & - \mu^\pm q \int_0^t \left [
   i \frac{H_{01}(qc_s^\pm t')}{|q|}
   \begin{pmatrix}
   0 & -k & 0 \\[5pt]
   k &  0 & m \\[5pt]
   0 & -m & 0 
   \end{pmatrix}
   \pm H_{11}(qc_s^\pm t')
   \begin{pmatrix}
   0 &  0 & 0 \\[5pt]
   0 &  1 & 0 \\[5pt]
   0 &  0 & 0 
   \end{pmatrix} \right.\\[10pt]
   &\quad \left.{} \pm \frac{H_{00}(qc_s^\pm t')}{q^2}
   \begin{pmatrix}
   k^2 & 0 & km \\[5pt]
   0   & 0 & 0 \\[5pt]
   km  & 0 & m^2 
   \end{pmatrix}
   \pm \frac{H_{22}(qc_s^\pm t')}{q^2}
   \begin{pmatrix}
   m^2  & 0 & -km \\[5pt]
   0    & 0 & 0 \\[5pt]
   -km  & 0 & k^2 
   \end{pmatrix} \right ]\\[10pt]
   &\times
   \begin{Bmatrix}
   U_0^\pm(t-t';k,m)\\[5pt]
   U_1^\pm(t-t';k,m)\\[5pt]
   U_2^\pm(t-t';k,m)
   \end{Bmatrix}
   |q| c_s^\pm \mathrm{d}t',

where :math:`q = \sqrt{k^2+m^2}`.

The convolution kernels are defined by the following inverse Laplace transforms

.. math::
   H_{00}(t) =& \mathfrak{L}^{-1} \left[ \frac{s^2 \sqrt{s^2+\eta^2}}{\sqrt{s^2+\eta^2}\sqrt{s^2+1}-\eta}-s  \right],\\
   H_{01}(t) =& \mathfrak{L}^{-1} \left[ \frac{-\eta s^2 }{\sqrt{s^2+\eta^2} \sqrt{s^2+1}-\eta} +\eta \right],\\
   H_{11}(t) =& \mathfrak{L}^{-1} \left[ \frac{ \eta s^2 \sqrt{s^2+1}}{\sqrt{s^2 + \eta^2} \sqrt{s^2+1} -\eta } -\eta s \right],

where :math:`s= p / |q|c_s` is the non-dimensional Laplace transform variable. The kernels can be inverted numerically.
The inverse Laplace transform is the following:

.. math::
   H(t) = \mathfrak{L}^{-1}[h(s)] = \frac{1}{2\pi i} \int_{\gamma-i\infty}^{\gamma+i\infty} \exp(s t ) h(s) \mathrm{d}s  .


We use a numerical inversion based on [de Hood, 1982].

The mode III kernel is defined as follows

.. math::
   H_{22}(t) = J_1(t)/t,


where :math:`J_1(t)`  is the Bessel function.

Quasi-dynamic formulation
^^^^^^^^^^^^^^^^^^^^^^^^^

The quasi-dynamic formulation neglects the wave-mediated stress transfer and only considers the radiation damping contribution, i.e., only a static stress transfer is considered. This is equivalent to saying that the displacement history was constant and corresponds to the current displacement. Hence, the convolutional terms become a simple multiplication of a constant displacement with the constant full integral of the convolution kernel.

For instance, 

.. math::
   \mu^\pm |q| \int_{0}^{t} H_{00}(|q|c_s^\pm t') U_0^\pm(t-t';q)|q|c_s^\pm \mathrm{d}t' = \mu^\pm |q| \int_{0}^{\infty} H_{00}(|q|c_s^\pm t')|q|c_s^\pm \mathrm{d}t' U_0^\pm(t;q)

where

.. math::
   \int_{0}^{t} H_{00}(|q|c_s^\pm t') |q|c_s^\pm \mathrm{d}t'

can be precomputed during initialitation, and remains constant throught the entire simulation.

Since the computation of the convolution term is one of the most computationally intensive operations in fully-dynamic simulations, this approximation reduces considerably the overal computational cost of the quasi-dynamic formulation compared to the fully-dynamic one. It also allows for adaptive time stepping, where time steps are larger than the dynamic stable time step, which constitutes another computational gain. 


Computation of interface tractions
----------------------------------

Bimaterial set-up
^^^^^^^^^^^^^^^^^

**What is the force needed to close the normal gap?**

Independent on whether an opening or an interpenetration would occur in the next time step. The condition to satisfy is that there is no gap after the computation of the displacement of the next time step. Therefore, it has to be predicted.

.. math::
  \delta_{n+2} &= 0 \\
  u_{n+2}^+ - u_{n+2}^- &= 0 \\
  u_{n+1}^+ - u_{n+1}^- + \Delta t (v_{n+1}^+ - v_{n+1}^-) &= 0 \\
  u_{n+1}^+ - u_{n+1}^- + \Delta t \frac{c_s^+ c_s^+}{\mu^+ c_d^+}\left( e_{n+1} + f_{n+1}^+ - i_{n+1} \right) + \Delta t \frac{c_s^- c_s^-}{\mu^- c_d^-}\left( e_{n+1} + f_{n+1}^- - i_{n+1} \right)&= 0 \\
  u_{n+1}^+ - u_{n+1}^- + \Delta t \frac{c_s^+ c_s^+}{\mu^+ c_d^+}\left( e_{n+1} + f_{n+1}^+\right) + \Delta t \frac{c_s^- c_s^-}{\mu^- c_d^-}\left( e_{n+1} + f_{n+1}^- \right) &= i_{n+1} \left( \Delta t \frac{c_s^+ c_s^+}{\mu^+ c_d^+} + \Delta t \frac{c_s^- c_s^-}{\mu^- c_d^-} \right)\\
  u_{n+1}^+ - u_{n+1}^- + C^+ \left( e_{n+1} + f_{n+1}^+\right) + C^-\left( e_{n+1} + f_{n+1}^- \right) &= i_{n+1} \left( C^++ C^- \right)\\

  
with :math:`C^+ = \Delta t \frac{c_s^+ c_s^+}{\mu^+ c_d^+}` and :math:`C^- = \Delta t \frac{c_s^- c_s^-}{\mu^- c_d^-}`.

We therefore can compute the nodal force needed to close the normal gap

.. math::
   i_{n+1} =  \left[ u_{n+1}^+ - u_{n+1}^- + C^+ \left( e_{n+1} + f_{n+1}^+\right) + C^- \left( e_{n+1} + f_{n+1}^- \right) \right] / \left( C^+ + C^- \right)

If :math:`i_{n+1} > 0`, then it is an adhesive force avoiding the gap to open, whereas :math:`i_{n+1} < 0` is the contact pressure that avoids interpenetration.


**What is the force needed to maintain the current tangential gap?**

Here the objective is not to close the gap, but to maintain the current opening. For this we only need to predict the velocity that is computed at the end of the current step.

.. math::
  v_{n+1}^+ &= v_{n+1}^- \\
  \frac{c_s^+}{\mu^+} \left( e_{n+1} + f_{n+1}^+ - i_{n+1} \right) &= -\frac{c_s^-}{\mu^-}\left( e_{n+1} + f_{n+1}^- - i_{n+1} \right)\\
  \left(\frac{c_s^+}{\mu^+} + \frac{c_s^-}{\mu^-} \right) i_{n+1} &= \frac{c_s^+}{\mu^+} \left( e_{n+1} + f_{n+1}^+\right) + \frac{c_s^-}{\mu^-}\left( e_{n+1} + f_{n+1}^-\right)\\
  

with :math:`D^+ = \frac{c_s^+}{\mu^+}` and :math:`D^- = \frac{c_s^-}{\mu^-}`.

We therefore can compute the nodal force needed to maintain the tangential gap

.. math::
   i_{n+1} = e_{n+1} + \frac{D^+ f_{n+1}^+ + D^- f_{n+1}^-}{D^+ + D^-} 

   
**How to apply strength?**

First, you compute the force that you need to either maintain the gap or close the gap, depending on the phenomenological law that you wish to apply. Then, you compute the strength based on your interface properties and state, and then make sure that the interface tractions are not larger than the strength. If they are, decrease it to the size of the strength, while maintaining orientation of the interface traction.


Deformable-rigid interface
^^^^^^^^^^^^^^^^^^^^^^^^^^

**What is the force needed to close the normal gap?**

Independent on whether an opening or an interpenetration would occur in the next time step. The condition to satisfy is that there is no gap after the computation of the displacement of the next time step. Therefore, it has to be predicted. In the deformable-rigid approach, we set :math:`u_{n+2}^- = 0`.

.. math::
   \delta_{n+2} &= 0 \\
   u_{n+2}^+ &= 0 \\
   u_{n+1}^+ + \Delta t v_{n+1}^+ &= 0 \\
   u_{n+1}^+ + \Delta t \frac{c_s^+ c_s^+}{\mu^+ c_d^+}\left( e_{n+1} + f_{n+1}^+ - i_{n+1} \right) &= 0 \\
   u_{n+1}^+ + \Delta t \frac{c_s^+ c_s^+}{\mu^+ c_d^+}\left( e_{n+1} + f_{n+1}^+\right) &= i_{n+1} \Delta t \frac{c_s^+ c_s^+}{\mu^+ c_d^+}\\
   u_{n+1}^+ + C^+ \left( e_{n+1} + f_{n+1}^+\right) &= i_{n+1} C^+\\
   
with :math:`C^+ = \Delta t \frac{c_s^+ c_s^+}{\mu^+ c_d^+}`.

We therefore can compute the nodal force needed to close the normal gap

.. math::
  i_{n+1} =  u_{n+1}^+ / C^+  + e_{n+1} + f_{n+1}^+

If :math:`i_{n+1} > 0`, then it is an adhesive force avoiding the gap to open, whereas :math:`i_{n+1} < 0` is the contact pressure that avoids interpenetration.


**What is the force needed to maintain the current tangential gap?**

Here the objective is not to close the gap, but to maintain the current opening. For this we only need to predict the velocity that is computed at the end of the current step. In the deformable-rigid approach, we set :math:`v_{n+1}^- = 0`.

.. math::
  v_{n+1}^+ &= v_{n+1}^- \\
  v_{n+1}^+ &= 0 \\
  \frac{c_s^+}{\mu^+} \left( e_{n+1} + f_{n+1}^+ - i_{n+1} \right) &= 0\\
  \frac{c_s^+}{\mu^+} i_{n+1} &= \frac{c_s^+}{\mu^+} \left( e_{n+1} + f_{n+1}^+\right)\\

We therefore can compute the nodal force needed to maintain the tangential gap

.. math::
   i_{n+1} = e_{n+1} + f_{n+1}^+   
  

Uni-material shear interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can use the antisymemtry condition (assumes that the interface remains in contact:

.. math::
   u_0^+(x) &= -u_0^-(x)\\
   f_0^+(x) &= -f_0^-(x)\\
   \quad\\
   u_1^+(x) &= u_1^-(x) \quad \text{no opening allowed}\\
   f_1^+(x) &= f_1^-(x) \quad \text{convolution response is symmetric}\\
   \quad\\
   u_2^+(x) &= -u_2^-(x)\\
   f_2^+(x) &= -f_2^-(x)\\
   
**What is the force needed to close the normal gap?**

We can use the same expression for :math:`i_{n+1}` as for the bi-material case, with :math:`C^+ = C^- = C`, which leads to:

.. math::
   i_{n+1} =  (u_{n+1}^+ - u_{n+1}^-) / 2C  + e_{n+1} + \frac{1}{2}(f_{n+1}^+ + f_{n+1}^-)



If :math:`i_{n+1} > 0`, then it is an adhesive force avoiding the gap to open, whereas :math:`i_{n+1} < 0` is the contact pressure that avoids interpenetration.


**What is the force needed to maintain the current tangential gap?**

We can use the same expression for :math:`i_{n+1}` as for the bi-material case, with :math:`D^+ = D^- = D = \frac{c_s}{\mu}`, which leads to:

.. math::
  i_{n+1} = e_{n+1} + \frac{1}{2} (f_{n+1}^+ + f_{n+1}^-)







