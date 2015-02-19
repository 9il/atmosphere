/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.likelihood.generalized_gamma;

import core.stdc.tgmath;

import std.traits;
import std.typecons;
import std.mathspecial : LN2, gamma;


/++
Normalized log-likelihood function of the generalized gamma distribution.
Params:
	shape = shape parameter
	power = power parameter
	scale = scale parameter
	sample = sample
+/
T generalizedGammaLikelihood(T)(T shape, T power, T scale, in T[] sample)
	if(isFloatingPoint!T)
in {
	import std.algorithm : all;
	assert(sample.all!"a > 0 && isNormal(a)");
}
body {
	T a = 0, b = 0;
	foreach(j; 0..sample.length)
	{
		a += pow(sample[j], power);
		b += log2(sample[j]);
	}
	b *= T(LN2);
	immutable n = sample.length;
	return generalizedGammaLikelihood!T(shape, power, scale, a/n, b/n);
}


/++
Normalized log-likelihood function of the generalized gamma distribution.
Params:
	shape = shape parameter
	power = power parameter
	scale = scale parameter
	sample = sample
	weights = weights for the sample
+/
T generalizedGammaLikelihood(T)(T shape, T power, T scale, in T[] sample, in T[] weights)
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
		a += w * pow(sample[j], power);
		b += w * log2(sample[j]);
	}
	b *= T(LN2);
	return generalizedGammaLikelihood!T(shape, power, scale, a/n, b/n);
}

///
unittest {
	immutable l = generalizedGammaLikelihood!double(1,2,3,[1,2,2]);
	immutable m = generalizedGammaLikelihood!double(1,2,3,[1,2],[2,4]);
	assert(l == m);
}


/++
Normalized log-likelihood function of the generalized gamma distribution.
Params:
	shape = shape parameter
	power = power parameter
	scale = scale parameter
	a = `Σ weights[j] * sample[j] ^^ power`
	b = `Σ weights[j] * log(sample[j])`
+/
T generalizedGammaLikelihood(T)(T shape, T power, T scale, T a, T b)
	if(isFloatingPoint!T)
{
	return 
		- log((scale * gamma(shape)) / fabs(power)) 
		- (1 - shape * power) * (b - log(scale)) 
		- (power > 0 ? a / pow(scale, power) : a * pow(scale, -power)); //precise unification with inverse-gamma and gamma
}
