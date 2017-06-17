Notations
*********

We first recall the *glt* symbols for the 1d mass, stiffness and advection respectively,

.. math::

    M^p_n &= \left[\int_0^1 N_{i_1}^p(t) ~ N_{j_1}^p(t) ~dt\right]_{i_1, j_1=1}^{n+p}, 
    \\
    S^p_n &= \left[\int_0^1 \left(N_{i_1}^p(t)\right)^{\prime} ~ \left(N_{j_1}^p(t)\right)^{\prime} ~dt\right]_{i_1, j_1=1}^{n+p}.
    \\
    A^p_n &= \left[\int \left(N_{i_1}^p\right)^{\prime}(t) ~ N_{j_1}^p(t) ~dt\right]_{i_1, j_1=1}^{n+p},

Their corresponding **GLT** symbols are

.. math::

  \{n M^p_n\}_n           & \sim_{\rm GLT} \mathfrak{m}_p 
  \\
  \{\frac{1}{n} S^p_n\}_n & \sim_{\rm GLT} \mathfrak{s}_p
  \\
  \{  A^p_n\}_n           & \sim_{\rm GLT} \mathfrak{a}_p 

where,

.. math::

   \mathfrak{m}_p(x, \theta) &:= \mathfrak{m}_p(\theta) = \phi_{2p+1}(p+1) + 2 \sum_{k=1}^p \phi_{2p+1}(p+1-k) \cos(k \theta).
   \\
   \mathfrak{s}_p(x, \theta) &:= \mathfrak{s}_p(\theta) = - {\phi}''_{2p+1}(p+1) - 2 \sum_{k=1}^p {\phi}''_{2p+1}(p+1-k) \cos(k \theta).
   \\
   \mathfrak{a}_p(x, \theta) &:= \mathfrak{a}_p(\theta) = - 2 \sum_{k=1}^p {\phi}'_{2p+1}(p+1-k) \sin(k \theta).

.. and :math:`\phi_{2p+1}` is the *Cardinal Spline* of degree :math:`2p+1`.
