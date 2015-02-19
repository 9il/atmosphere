/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.distribution.likelihood.inverse_gamma;

import core.stdc.tgmath;

import std.traits;
import std.typecons;
import std.mathspecial : LN2, gamma;


/++
Normalized log-likelihood function of the inverse-gamma distribution.
Params:
	shape = shape parameter
	scale = scale parameter
	sample = sample
+/
T inverseGammaLikelihood(T)(T shape, T scale, in T[] sample)
	if(isFloatingPoint!T)
in {
	import std.algorithm : all;
	assert(sample.all!"a > 0 && isNormal(a)");
}
body {
	T a = 0, b = 0;
	foreach(j; 0..sample.length)
	{
		a += 1 / sample[j];
		b += log2(sample[j]);
	}
	b *= T(LN2);
	immutable n = sample.length;
	return inverseGammaLikelihood!T(shape, scale, a/n, b/n);
}

///
unittest {
	import atmosphere.distribution.likelihood.generalized_gamma;
	immutable l = inverseGammaLikelihood!double(2,3,[1,2]);
	immutable m = generalizedGammaLikelihood!double(2,-1,3,[1,2]);
	assert(l == m);
}


/++
Normalized log-likelihood function of the inverse-gamma distribution.
Params:
	shape = shape parameter
	scale = scale parameter
	sample = sample
	weights = weights for the sample
+/
T inverseGammaLikelihood(T)(T shape, T scale, in T[] sample, in T[] weights)
	if(isFloatingPoint!T)
in {
	assert(weights.length == sample.length);
	import std.algorithm : all, any;
	assert(weights.length == sample.length);
	assert(sample.all!"a > 0 && isNormal(a)");
	assert(weights.all!"a >= 0 && isFinite(a)");
	assert(weights.any!"a > 0");
}
body {
	T a = 0, b = 0, n = 0;
	foreach(j; 0..sample.length)
	{
		immutable w = weights[j];
		n += w;
		a += w / sample[j];
		b += w * log2(sample[j]);
	}
	b *= T(LN2);
	return inverseGammaLikelihood!T(shape, scale, a/n, b/n);
}

///
unittest {
	import atmosphere.distribution.likelihood.generalized_gamma;
	immutable l = inverseGammaLikelihood!double(2,3,[1,2],[3,4]);
	immutable m = generalizedGammaLikelihood!double(2,-1,3,[1,2],[3,4]);
	assert(l == m);
}

///
unittest {
	immutable l = inverseGammaLikelihood!double(1,3,[1,2,2]);
	immutable m = inverseGammaLikelihood!double(1,3,[1,2],[2,4]);
	assert(l == m);
}


/++
Normalized log-likelihood function of the inverse-gamma distribution.
Params:
	shape = shape parameter
	scale = scale parameter
	a = `Σ weights[j] / sample[j] / Σ weights[j]`
	b = `Σ weights[j] * log(sample[j]) / Σ weights[j]`
+/
T inverseGammaLikelihood(T)(T shape, T scale, T a, T b)
	if(isFloatingPoint!T)
{
	return 
		- log(scale * gamma(shape))
		- (1 + shape) * (b - log(scale)) 
		- a * scale;
}

///
unittest {
	import atmosphere.distribution.likelihood.generalized_gamma;
	immutable l = inverseGammaLikelihood(4.0, 3.0, 2.0, 1.0);
	immutable m = generalizedGammaLikelihood(4.0, -1.0, 3.0, 2.0, 1.0);
	assert(l == m);
}
