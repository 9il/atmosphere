/++
This module publicly imports all classes.

Import:
--------
import atmosphere;
--------

Class_hierarchy:
- [MixtureOptimizer](atmosphere/mixture/MixtureOptimizer.html)
	- [ExpectationMaximization](atmosphere/mixture/ExpectationMaximization.html)
		- [LikelihoodAscentEM](atmosphere/mixture/LikelihoodAscentEM.html) : [LikelihoodAscent](atmosphere/mixture/LikelihoodAscent.html)
	- [GradientDescent](atmosphere/mixture/GradientDescent.html)
	- [GradientDescentPartial](atmosphere/mixture/GradientDescentPartial.html)
		- [LikelihoodAscentGradient](atmosphere/mixture/LikelihoodAscentGradient.html) : [LikelihoodAscent](atmosphere/mixture/LikelihoodAscent.html)
	- [CoordinateDescent](atmosphere/mixture/CoordinateDescent.html)
	- [CoordinateDescentPartial](atmosphere/mixture/CoordinateDescentPartial.html)
		- [LikelihoodAscentCoordinate](atmosphere/mixture/LikelihoodAscentCoordinate.html) : [LikelihoodAscent](atmosphere/mixture/LikelihoodAscent.html)
	- [NvmmLikelihoodAscentEM](atmosphere/parametrized/nvmm/NvmmLikelihoodAscentEM.html) : [LikelihoodAscent](atmosphere/mixture/LikelihoodAscent.html)
		- [NvmmLikelihoodAscentEMEM](atmosphere/parametrized/nvmm/NvmmLikelihoodAscentEMEM.html)
		- [NvmmLikelihoodAscentEMGradient](atmosphere/parametrized/nvmm/NvmmLikelihoodAscentEMGradient.html)
		- [NvmmLikelihoodAscentEMCoordinate](atmosphere/parametrized/nvmm/NvmmLikelihoodAscentEMCoordinate.html)

See_Also: [Atmosphere GM Test](https://github.com/9il/atmosphere_gm_test)

+/

/+

Table of separating mixtures algorithms:
	$(TABLE
		$(TR  $(TH Class) $(TH BLAS Level 2 % of all computations) $(TH Non-BLAS parallel) $(TH   Speed (one thread & big data) ) $(TH Global maximum for convex function) $(TH Parameterized))
		$(TR $(TD [GradientDescent](atmosphere/mixture/GradientDescent.html), number of mixture components > ~16) $(TD ~95-99%) $(TD No) $(TD  normal ) $(TD Yes) $(TD No) )
		$(TR $(TD [GradientDescent](atmosphere/mixture/GradientDescent.html), number of mixture components < ~16) $(TD ~95-99%) $(TD No) $(TD  fast ) $(TD Yes) $(TD No) )
		$(TR $(TD [CoordinateDescent](atmosphere/mixture/CoordinateDescent.html)) $(TD ~10-50%) $(TD No) $(TD  fast ) $(TD Yes) $(TD No) )
		$(TR $(TD [CoordinateDescentPartial](atmosphere/mixture/LikelihoodAscentCoordinate.html)) $(TD ~10-50%) $(TD No) $(TD  fast, always faster then CoordinateDescent ) $(TD Yes) $(TD No) )
		$(TR $(TD [NvmmLikelihoodAscentEMEM](atmosphere/parametrized/nvmm/NvmmLikelihoodAscentEMEM.html)) $(TD ~10%) $(TD Yes) $(TD  slow ) $(TD No) $(TD Yes) )
		$(TR $(TD [NvmmLikelihoodAscentEMGradient](atmosphere/parametrized/nvmm/NvmmLikelihoodAscentEMGradient.html)) $(TD ~10%) $(TD Yes) $(TD  slow ) $(TD No) $(TD Yes) )
		$(TR $(TD [NvmmLikelihoodAscentEMCoordinate](atmosphere/parametrized/nvmm/NvmmLikelihoodAscentEMCoordinate.html)) $(TD ~5%) $(TD Partial) $(TD  slow ) $(TD No) $(TD Yes) ) )
+/
module atmosphere;

public import atmosphere.mixture;
public import atmosphere.parametrized.nvmm;
