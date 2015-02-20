/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.likelihood.generalized_inverse_gaussian;

import core.stdc.tgmath;

import std.traits;
import std.typecons;
import std.math : LN2;
import atmosphere.utilities : sumOfLog2s, dotProduct;

/++
Normalized log-likelihood function of the generalized inverse Gaussian distribution.
Params:
	lambda = parameter lambda
	eta = scale parameter, eta
	omega = parameter omega
	sample = sample

See_Also: `distribution.params.GIGEtaOmega`
+/
T generalizedInverseGaussianLikelihood(T)(T lambda, T eta, T omega, in T[] sample)
	if(isFloatingPoint!T)
{
	import std.algorithm : sum, map;
	immutable n = sample.length;
	immutable one = sample.sum() / n;
	immutable mone = sample.map!"1/a".sum() / n;
	immutable dzero = T(LN2) * sample.sumOfLog2s() / n;
	return generalizedInverseGaussianLikelihood(lambda, eta, omega, one, mone, dzero);
}


/++
Normalized log-likelihood function of the generalized inverse Gaussian distribution.
Params:
	lambda = parameter lambda
	eta = scale parameter, eta
	omega = parameter omega
	sample = sample
	weights = weights for the sample

See_Also: `distribution.params.GIGEtaOmega`
+/
T generalizedInverseGaussianLikelihood(T)(T lambda, T eta, T omega, in T[] sample, in T[] weights)
	if(isFloatingPoint!T)
in {
	assert(weights.length == sample.length);
}
body {
	import std.algorithm : sum, map;
	immutable n = weights.sum;
	immutable one = sample.dotProduct(weights) / n;
	immutable mone = sample.map!"1/a".dotProduct(weights) / n;
	immutable dzero = T(LN2) * sample.map!log2.dotProduct(weights) / n;
	return generalizedInverseGaussianLikelihood(lambda, eta, omega, one, mone, dzero);
}

///
unittest {
	immutable l = generalizedInverseGaussianLikelihood!double(1,2,3,[1,2,2]);
	immutable m = generalizedInverseGaussianLikelihood!double(1,2,3,[1,2],[2,4]);
	assert(l == m);
}


/++
Normalized log-likelihood function of the generalized inverse Gaussian distribution.
Params:
	lambda = parameter lambda
	eta = scale parameter, eta
	omega = parameter omega
	one = `Σ weights[j] * sample[j] / Σ weights[j]`
	mone = `Σ weights[j] / sample[j] / Σ weights[j]`
	dzero = `Σ weights[j] * log(sample[j]) / Σ weights[j]`

See_Also: `distribution.params.GIGEtaOmega`
+/
T generalizedInverseGaussianLikelihood(T)(T lambda, T eta, T omega, T one, T mone, T dzero)
	if(isFloatingPoint!T)
{
	import bessel;
	import std.typecons : Flag;
	return 
		- log(2 * eta * besselK(omega, lambda, Flag!"ExponentiallyScaled".yes))
		+ (lambda - 1) * (dzero - log(eta))
		+ omega * (1 - (one / eta + mone * eta) / 2);
}
