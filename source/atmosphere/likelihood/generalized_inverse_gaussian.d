/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.likelihood.generalized_inverse_gaussian;

import std.traits;
import std.typecons;
import std.math;

import atmosphere.statistic: GeneralizedInverseGaussinStatistic;

/++
Normalized log-likelihood function of the generalized inverse Gaussian distribution.
See_Also: `distribution.params.GIGEtaOmega`
+/
T properGeneralizedInverseGaussianLikelihood(T)(in T lambda, in T eta, in T omega, in T[] sample)
	if(isFloatingPoint!T)
{
	return properGeneralizedInverseGaussianLikelihood(lambda, eta, omega, GeneralizedInverseGaussinStatistic!T(sample));
}

///ditto
T properGeneralizedInverseGaussianLikelihood(T)(in T lambda, in T eta, in T omega, in T[] sample, in T[] weights)
	if(isFloatingPoint!T)
{
	return properGeneralizedInverseGaussianLikelihood(lambda, eta, omega, GeneralizedInverseGaussinStatistic!T(sample, weights));
}

///
unittest {
	immutable l = properGeneralizedInverseGaussianLikelihood!double(1,2,3,[1,2,2]);
	immutable m = properGeneralizedInverseGaussianLikelihood!double(1,2,3,[1,2],[2,4]);
	assert(l == m);
}

///ditto
T properGeneralizedInverseGaussianLikelihood(T)(in T lambda, in T eta, in T omega, in GeneralizedInverseGaussinStatistic!T stat)
	if(isFloatingPoint!T)
in {
	assert(omega >= T.min_normal);
}
body {
	version(LDC)
	{
		import ldc.intrinsics: log = llvm_log;
	}
	import std.math : LN2;
	import atmosphere.math: logBesselK;
	with(stat) return
		- cast(T)LN2 - lambda * log(eta) - logBesselK(lambda, omega)
		+ (lambda - 1) * stat.meanl
		- omega * (mean / eta + meani * eta) / 2;
}
