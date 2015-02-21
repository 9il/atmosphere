/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.likelihood.generalized_inverse_gaussian;

import core.stdc.tgmath;

import std.traits;
import std.typecons;

import atmosphere.statistic: GeneralizedInverseGaussinStatistic;

/++
Normalized log-likelihood function of the generalized inverse Gaussian distribution.
See_Also: `distribution.params.GIGEtaOmega`
+/
T generalizedInverseGaussianLikelihood(T)(T lambda, T eta, T omega, in T[] sample)
	if(isFloatingPoint!T)
{
	return generalizedInverseGaussianLikelihood(lambda, eta, omega, GeneralizedInverseGaussinStatistic!T(sample));
}

///ditto
T generalizedInverseGaussianLikelihood(T)(T lambda, T eta, T omega, in T[] sample, in T[] weights)
	if(isFloatingPoint!T)
{
	return generalizedInverseGaussianLikelihood(lambda, eta, omega, GeneralizedInverseGaussinStatistic!T(sample, weights));
}

///
unittest {
	immutable l = generalizedInverseGaussianLikelihood!double(1,2,3,[1,2,2]);
	immutable m = generalizedInverseGaussianLikelihood!double(1,2,3,[1,2],[2,4]);
	assert(l == m);
}

///ditto
T generalizedInverseGaussianLikelihood(T)(T lambda, T eta, T omega, GeneralizedInverseGaussinStatistic!T stat)
	if(isFloatingPoint!T)
{
	import bessel;
	import std.typecons : Flag;
	with(stat) return
		- log(2 * eta * besselK(omega, lambda, Flag!"ExponentiallyScaled".yes))
		+ (lambda - 1) * (stat.meanl - log(eta))
		+ omega * (1 - (mean / eta + meani * eta) / 2);
}
