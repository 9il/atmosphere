/++
This module publicly imports all classes.

Import:
--------
import atmosphere;
--------

Class_hierarchy:
- [MixtureOptimizer](atmosphere/mixture/MixtureOptimizer.html)
	- [EMDescent](atmosphere/mixture/EMDescent.html)
		- [EMLikelihoodMaximization](atmosphere/mixture/EMLikelihoodMaximization.html) : [LikelihoodMaximization](atmosphere/mixture/LikelihoodMaximization.html)
	- [GradientDescent](atmosphere/mixture/GradientDescent.html)
	- [GradientDescentPartial](atmosphere/mixture/GradientDescentPartial.html)
		- [GradientLikelihoodMaximization](atmosphere/mixture/GradientLikelihoodMaximization.html) : [LikelihoodMaximization](atmosphere/mixture/LikelihoodMaximization.html)
	- [CoordinateDescent](atmosphere/mixture/CoordinateDescent.html)
	- [CoordinateDescentPartial](atmosphere/mixture/CoordinateDescentPartial.html)
		- [CoordinateLikelihoodMaximization](atmosphere/mixture/CoordinateLikelihoodMaximization.html) : [LikelihoodMaximization](atmosphere/mixture/LikelihoodMaximization.html)
	- [NormalVarianceMeanMixture](atmosphere/parametrized/nvmm/NormalVarianceMeanMixture.html) : [LikelihoodMaximization](atmosphere/mixture/LikelihoodMaximization.html)
		- [NormalVarianceMeanMixtureEM](atmosphere/parametrized/nvmm/NormalVarianceMeanMixtureEM.html)
		- [NormalVarianceMeanMixtureEMAndGradient](atmosphere/parametrized/nvmm/NormalVarianceMeanMixtureEMAndGradient.html)
		- [NormalVarianceMeanMixtureEMAndCoordinate](atmosphere/parametrized/nvmm/NormalVarianceMeanMixtureEMAndCoordinate.html)

Table of separating mixtures algorithms:
	$(TABLE
		$(TR  $(TH Class) $(TH BLAS Level 2 % of all computations) $(TH Non-BLAS parallel) $(TH   Speed (one thread & big data) ) $(TH Global maximum for convex function) $(TH Parameterized))
		$(TR $(TD [GradientDescent](atmosphere/mixture/GradientDescent.html), number of mixture components > ~16) $(TD ~95-99%) $(TD No) $(TD  normal ) $(TD Yes) $(TD No) )
		$(TR $(TD [GradientDescent](atmosphere/mixture/GradientDescent.html), number of mixture components < ~16) $(TD ~95-99%) $(TD No) $(TD  fast ) $(TD Yes) $(TD No) )
		$(TR $(TD [CoordinateDescent](atmosphere/mixture/CoordinateDescent.html)) $(TD ~10-50%) $(TD No) $(TD  fast ) $(TD Yes) $(TD No) )
		$(TR $(TD [CoordinateDescentPartial](atmosphere/mixture/CoordinateLikelihoodMaximization.html)) $(TD ~10-50%) $(TD No) $(TD  fast, always faster then CoordinateDescent ) $(TD Yes) $(TD No) )
		$(TR $(TD [NormalVarianceMeanMixtureEM](atmosphere/parametrized/nvmm/NormalVarianceMeanMixtureEM.html)) $(TD ~10%) $(TD Yes) $(TD  slow ) $(TD No) $(TD Yes) )
		$(TR $(TD [NormalVarianceMeanMixtureEMAndGradient](atmosphere/parametrized/nvmm/NormalVarianceMeanMixtureEMAndGradient.html)) $(TD ~10%) $(TD Yes) $(TD  slow ) $(TD No) $(TD Yes) )
		$(TR $(TD [NormalVarianceMeanMixtureEMAndCoordinate](atmosphere/parametrized/nvmm/NormalVarianceMeanMixtureEMAndCoordinate.html)) $(TD ~5%) $(TD Partial) $(TD  slow ) $(TD No) $(TD Yes) ) )

See_Also: [Atmosphere GM Test](https://github.com/9il/atmosphere_gm_test)

+/
module atmosphere;

public import atmosphere.mixture;
public import atmosphere.parametrized.nvmm;
