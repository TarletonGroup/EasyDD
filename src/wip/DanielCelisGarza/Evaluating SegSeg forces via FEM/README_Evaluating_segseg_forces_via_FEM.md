# Evaluating SegSeg forces via FEM

Find the volumetric stress by averaging the field point stresses inside the element volume $V$,

$\sigma_v = \dfrac{1}{V} \int\limits_V \bm{\sigma}(\vec{x})\, \mathrm{d}V\, .$

We can approximate this integral using Gauss quadrature and the field point stresses,

$\sigma_v \approx \dfrac{1}{V} \sum\limits_{i=1}^I w_i \sum\limits_{j=1}^J w_j \sum\limits_{k=1}^K w_k \, \bm{\sigma}(x_i,\,x_j,\,x_k)\, .$

The forces acting on node $n$ are obtained by integrating the volumetric stress, multiplied by the derivative of the shape function $N_n$,

$F_n = \sigma_v \int\limits_V \nabla N_n\, \mathrm{d}V\, .$

For $n$-th degree polynomial shape functions, Gauss quadrature with $2n-1$ points per corresponding dimension gives an exact numerical answer,

$F_n \approx \sigma_v \, \sum\limits_{i=1}^I w_i \sum\limits_{j=1}^J w_j \sum\limits_{k=1}^K w_k \, \nabla N_n(x_i,y_j,z_k)\, .$

Since we use linear shape functions in each dimension, we only need a single quadrature point in the middle of the element to exactly compute the integral.
