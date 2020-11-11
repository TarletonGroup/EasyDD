# Evaluating SegSeg forces via FEM

Take an FE element. Evaluate the stress on the nodes.

$\sigma_v = \dfrac{1}{V} \int\limits_V \bm{\sigma}(\vec{x})\, \mathrm{d}V\, .$

Using gauss quadrature we can approximate this integral,

$\sigma_v \approx \dfrac{1}{V} \sum\limits_{i=1}^I w_i \sum\limits_{j=1}^J w_j \sum\limits_{k=1}^K w_k \, \bm{\sigma}(x_i,\,x_j,\,x_k)\, .$

Then the forces in node $n$ are obtained by integrating the volumetric stress multiplied by the derivative of the shape function $N_n$
$F_n = \int\limits_V (\sigma_v \cdot \nabla N_n)\, \mathrm{d}V\, .$

Again we can approximate this integral using gauss quadrature,
$F_n \approx \sum\limits_{i=1}^I w_i \sum\limits_{j=1}^J w_j \sum\limits_{k=1}^K w_k \, \sigma_v \cdot \nabla N_n(x_i,y_j,z_k)\, .$