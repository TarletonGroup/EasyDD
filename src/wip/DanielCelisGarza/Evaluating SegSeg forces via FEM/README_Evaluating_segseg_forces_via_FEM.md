# Evaluating SegSeg forces via FEM

Find the volumetric stress by averaging the field point stresses inside the element volume $V$,

$\bm{\sigma}_v = \dfrac{1}{V} \int\limits_V \bm{\sigma}(\vec{x})\, \mathrm{d}V\, .$

We can approximate this integral using Gauss quadrature and the field point stresses,

$\bm{\sigma}_v \approx \dfrac{1}{V} \sum\limits_{i=1}^I w_i \sum\limits_{j=1}^J w_j \sum\limits_{k=1}^K w_k \, \bm{\sigma}(x_i,\,x_j,\,x_k)\, .$

The forces acting on node $n$ are obtained by integrating the volumetric stress, multiplied by the derivative of the shape function $\vec{N}_n$,

$\vec{F}_n = \int\limits_V (\bm{\sigma}_v \cdot \nabla \vec{N}_n)\, \mathrm{d}V\, .$

For $n$-th degree polynomial shape functions, Gauss quadrature with $2n-1$ points per corresponding dimension gives an exact numerical answer,

$\vec{F}_n \approx \sum\limits_{i=1}^I w_i \sum\limits_{j=1}^J w_j \sum\limits_{k=1}^K w_k \, [\bm{\sigma}_v \cdot \nabla \vec{N}_n(x_i,y_j,z_k)]\, .$

Since we use linear shape functions in each dimension, we only need a single quadrature point in the middle of the element to exactly compute the integral.
