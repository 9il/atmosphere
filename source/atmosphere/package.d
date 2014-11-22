/**
This module publicly imports all classes.

Usage:
--------
import atmosphere;
--------

Class_hierarchy:
$(D_CODE
$(DPREF mixture, MixtureOptimizer) $(BLUE abstract)
|-$(DPREF mixture, GradientDescent)
| |-$(DPREF mixture, GradientLikelihoodMaximization) : $(DPREF mixture, LikelihoodMaximization)
|-$(DPREF mixture, CoordinateDescent)
|-$(DPREF mixture, CoordinateDescentPartial)
| |-$(DPREF mixture, CoordinateLikelihoodMaximization) : $(DPREF mixture, LikelihoodMaximization)
|-$(DPREF parametrized.nvmm, NormalVarianceMeanMixture) : $(DPREF mixture, LikelihoodMaximization) $(BLUE abstract)
  |-$(DPREF parametrized.nvmm, NormalVarianceMeanMixtureEM) $(BLUE final)
  |-$(DPREF parametrized.nvmm, NormalVarianceMeanMixtureEMAndGradient) $(BLUE final)
  |-$(DPREF parametrized.nvmm, NormalVarianceMeanMixtureEMAndCoordinate) $(BLUE final)
)

<table class = "table table-condensed table-bordered">
<tr><th>Class</th><th>BLAS Level 2 % of all computations</th><th>Non-BLAS parallel</th><th>  Speed (one thread & big data) </th><th>Global maximum for convex function</th><th>Parameterized</th></tr>
<tr><td>GradientDescent, number of mixture components > ~16</td><td>~95-99%</td><td>No</td><td> normal </td><td>Yes</td><td>No</td></tr>
<tr><td>GradientDescent, number of mixture components < ~16</td><td>~95-99%</td><td>No</td><td> fast </td><td>Yes</td><td>No</td></tr>
<tr><td>CoordinateDescent</td><td>~10-50%</td><td>No</td><td> fast </td><td>Yes</td><td>No</td></tr>
<tr><td>CoordinateDescentPartial</td><td>~10-50%</td><td>No</td><td> fast, always faster then CoordinateDescent </td><td>Yes</td><td>No</td></tr>
<tr><td>NormalVarianceMeanMixtureEM</td><td>~10%</td><td>Yes</td><td> slow </td><td>No</td><td>Yes</td></tr>
<tr><td>NormalVarianceMeanMixtureEMAndGradient</td><td>~10%</td><td>Yes</td><td> slow </td><td>No</td><td>Yes</td></tr>
<tr><td>NormalVarianceMeanMixtureEMAndCoordinate</td><td>~5%</td><td>Partial</td><td> slow </td><td>No</td><td>Yes</td></tr>
</table>

Example: 
	$(HTTP github.com/9il/atmosphere_gm/blob/master/examples/normal_variance_mean_mixture/source/app.d, GitHub) 

If you want to enable non-BLAS parallelism use $(D atmosphere_gm_parallel) compilation version.
*/
module atmosphere;

public import atmosphere.mixture;
public import atmosphere.parametrized.nvmm;
