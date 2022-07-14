Interface Laws
==============

Here, we provide a description of all existing interface laws. In addition, you may add your own interface law (see :doc:`developer_guide`).

.. used subsubsections for laws so that they could be group with subsections underlined with -----

General interface mechanics
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before we describe the interface laws individually, we summarize here the general mechanics of interfaces.

**Normal component: contact pressure and adhesion**:

First, we introduce the :math:`\sigma_n^{close\,gap}`. This is the normal traction at the interface, which would be required to close the gap. On one hand, if the two solids are about to penetrate each other, *i.e.* the gap becomes negative, :math:`\sigma_n^{close\,gap}` prevents the penetration and closes the gap. It corresponds to the contact pressure. On the other hand, if the interface is about to open, *i.e.* the gap becomes positive, :math:`\sigma_n^{close\,gap}` prevents opening and closes the gap. Hence, it corresponds to adhesion/cohesion. 

The interface would never open if the adhesive force was perfect. However, most interfaces are characterized by an adhesive strength :math:`\tilde \sigma_n`, which cannot be exceeded. Hence, the actual normal interface tractions are given by:

.. math::
   \sigma_n = \min(\sigma_n^{close\,gap},\tilde \sigma_n)

which implies that:

- if :math:`\sigma_n^{close\,gap}>\tilde \sigma_n`  the gap opens

- if :math:`\sigma_n^{close\,gap}<\tilde \sigma_n`  the gap remains closed and interpenetration is prevented.

The constitutive interface law for normal direction describes therefore :math:`\tilde \sigma_n` as function of relevant properties. This is what is implemented in the interface laws described below.


**Shear component: friction (and possibly adhesion)**:

Similar to the normal direction, we first introduce :math:`\tau_s^{cst \, slip}`. This is the shear traction at the interface, which would be required to maintain the shear gap (*i.e.* slip) constant. In this case, the interface would be stuck. 

The interface would never slip if the shear traction was perfect. However, most interfaces have a limited shear strength :math:`\tilde \tau_s`. Hence, the actual shear interface traction is given by:

.. math::
   \tau_s = \min\left( ||\tau_s^{cst \, slip}|| , \tilde \tau_s\right)  \frac{\tau_s^{cst \, slip}}{||\tau_s^{cst \, slip}||}
   
which implies that:

- if :math:`\tau_s^{cst \, slip}>\tilde \tau_s`  the slip increases

- if :math:`\tau_s^{cst \, slip}<\tilde \tau_s`  the slip is maintained (*i.e.* stick)

The constitutive interface law for shear direction describes therefore :math:`\tilde \tau_s` as function of relevant properties. This is what is implemented in the interface laws described below.


Linear Normal Cohesive Law
^^^^^^^^^^^^^^^^^^^^^^^^^^

The `Linear Normal Cohesive Law` describes the normal cohesive strength :math:`\tilde \sigma_n` as a linear function of the opening, as given by the following expression:

.. math::
   \tilde \sigma_n (\delta_n) = \max\left( \sigma_c - \delta_n  W,  0\right)

where :math:`W=\sigma_c^2 / 2 \Gamma_c` is the weakening rate, :math:`\delta_n` the opening, :math:`\Gamma_c` the fracture energy, and :math:`\sigma_c` the peak strengths.

The ``LinearNormalCohesiveLaw`` constructor uses :math:`\Gamma_c`, and :math:`\sigma_c` as arguments.

In the shear direction, this law maintains the current shear opening. It does not allow it to grow more.



Linear Shear Cohesive Law
^^^^^^^^^^^^^^^^^^^^^^^^^

The `Linear Shear Cohesive Law` assumes that the normal cohesive strength is :math:`\tilde \sigma_n=0` (analogous to a frictional interface). Hence, the interface is allowed to freely open without any resistance against it.

Further, this interface law assumes that the cohesive strength in the shear direction :math:`\tau_s` and is given by the following expression:

.. math::
   \tilde \tau_s (\delta_s) = \max\left( \tau_c - \delta_s  W,  \tau_r\right)

where :math:`W=0.5(\tau_c - \tau_r)^2/ \Gamma_c` is the weakening rate, :math:`\delta_s` is the slip, :math:`\Gamma_c` the fracture energy, and :math:`\tau_c` and :math:`\tau_r` are the peak and residual strengths, respectively.

The ``LinearShearCohesiveLaw`` constructor uses :math:`\Gamma_c`, :math:`\tau_c`, and :math:`\tau_r` as arguments.



Linear Coulomb Friction Law
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `Linear Coulomb Friction Law` applies the exact same procedure as the `Linear Shear Cohesive Law` (including :math:`\tilde \sigma_n=0`) with the only exception that the shear peak and residual shear strengths :math:`\tau_c=\mu_s \sigma_n` and :math:`\tau_r=\mu_k \sigma_n` are coupled to the local normal traction (as in Coulomb friction), where :math:`\mu_s` and :math:`\mu_k` are the static and kinetic friction coefficients, respectively.

.. math::
   \tilde \tau_s=\sigma_n \max\left( \mu_s -  \delta_s\frac{\mu_s-\mu_k}{d_c}, \mu_k \right)

where :math:`d_c` is the characteristic slip weakening distance.

Additionally, there is the possibility to introduce a characteristic regularization time :math:`t^{reg}`, which smooths the temporal evolution of the normal stress as follows:

.. math::
   \sigma_n^{reg}(t+\Delta t) = \frac{\sigma_n^{reg}(t) +\frac{\Delta t}{ t^{reg}} \sigma_n}{1+\frac{\Delta t }{t^{reg}}}
   
As a consequence also the temporal variation of the friction strength is regularized. This can be useful to improve numerical stability `(Kammer et al. 2014) <http://www.sciencedirect.com/science/article/pii/S0022509613002159>`_.

The ``LinearCoulombFrictionLaw`` constructor uses :math:`\mu_s`, :math:`\mu_k`, :math:`d_c`, and optionally :math:`t^{reg}` as arguments.

Barras Law (Mixed-Mode Fracture)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `Barras Law` (Mixed-Mode Fracture) considers a cohesive interface, which strength decreases linearly with increasing gap opening. This cohesive law is equivalent to the law used in `Barras et al. (2014) <http://link.springer.com/article/10.1007/s10704-014-9967-z>`_ but neglects friction which only comes into play after failure occurred.

.. math::
   \tilde \tau = \tau_c \max \left(0,1- ||\delta||/d_c \right)

where :math:`\tilde \sigma_n = \tilde \tau` as well as :math:`\tilde \tau_s = \tilde \tau`.

The ``BarrasLaw`` constructor uses :math:`\tau_c`, and :math:`d_c` as arguments.


Rate- and State- Friction Law
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `Rate- and State- Friction Law` considers an interface with frictional strength dependent on slip rate and state. The normal strength is zero, *i.e.* :math:`\tilde \sigma_n=0`. We implement the rate- and state- friction law following `Lapusta et al. (2000) <http://doi.wiley.com/10.1029/2000JB900250>`_. The interface shear strength is given by:

.. math::
   \tilde \tau_s=a \sigma_n \mathrm{arcsinh}\left[ \frac{V}{2V_0} +\exp\left(\frac{f_0 + b \ln \frac{V_0 \theta}{d_c}}{a}\right)\right]

where :math:`a` and :math:`b` are friction properties of the interface. The friction depends on the slip rate :math:`V` and the state variable :math:`\theta`. :math:`V_0` and :math:`f_0` are reference slip velocity and friction coefficient, respectively.

The temporal evolution of the state variable can be characterized by the aging law

.. math::
   \dot \theta = 1- \frac{V\theta}{d_c}

or by the slip law

.. math::
   \dot \theta= \frac{V\theta}{d_c}\ln(V\theta/d_c)

The ``RateAndStateLaw`` constructor uses all of these parameters as arguments.

   
