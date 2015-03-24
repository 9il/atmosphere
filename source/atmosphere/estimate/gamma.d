/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.estimate.gamma;

import core.stdc.tgmath;

import std.traits;
import std.typecons;

import atmosphere.statistic: GammaStatistic;

/++
Estimates parameters of the gamma distribution.
+/
Tuple!(T, "shape", T, "scale")
gammaEstimate(T)(in T[] sample)
	if(isFloatingPoint!T)
{
	return gammaEstimate(GammaStatistic!T(sample));
}

///
unittest {
	import atmosphere.estimate.generalized_gamma;
	immutable sample = [1.0, 0.5, 0.75];
	immutable p0 = gammaEstimate(sample);
	immutable p1 = generalizedGammaFixedPowerEstimate(1.0, sample);
	assert(p0.shape == p1.shape);
	assert(p0.scale == p1.scale);
}


///
unittest
{
	import std.range;
	import std.random;
	import atmosphere.random;
	import atmosphere.likelihood.gamma;
	auto length = 1000;
	auto shape = 2.0, scale = 3.0;
	auto rng = Random(1234);
	auto sample = GammaSRNG!double(rng, shape, scale).take(length).array;
	auto weights = iota(1.0, length + 1.0).array;
	auto params = gammaEstimate!double(sample, weights);
	auto lh0 = gammaLikelihood(shape, scale, sample, weights);
	auto lh1 = gammaLikelihood(params.shape, params.scale, sample, weights);
	assert(lh0 <= lh1);
}

///ditto
Tuple!(T, "shape", T, "scale")
gammaEstimate(T)(in T[] sample, in T[] weights)
	if(isFloatingPoint!T)
{
	return gammaEstimate(GammaStatistic!T(sample, weights));
}

///
unittest {
	import atmosphere.estimate.generalized_gamma;
	immutable sample = [1.0, 0.5, 0.75];
	immutable weights = [1.0, 4, 3];
	immutable p0 = gammaEstimate(sample, weights);
	immutable p1 = generalizedGammaFixedPowerEstimate(1.0, sample, weights);
	assert(p0.shape == p1.shape);
	assert(p0.scale == p1.scale);
}

///ditto
Tuple!(T, "shape", T, "scale")
gammaEstimate(T)(in GammaStatistic!T stat)
	if(isFloatingPoint!T)
{
	import atmosphere.math: logmdigammaInverse;
	with(stat)
	{
		immutable T shape = logmdigammaInverse(log(mean) - meanl);
		immutable T scale = mean / shape;
		return typeof(return)(shape, scale);
	}
}
