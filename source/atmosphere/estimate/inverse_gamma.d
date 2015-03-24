/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.estimate.inverse_gamma;

import core.stdc.tgmath;

import std.traits;
import std.typecons;

import atmosphere.statistic: InverseGammaStatistic;

/++
Estimates parameters of the inverse-gamma distribution.
+/
Tuple!(T, "shape", T, "scale")
inverseGammaEstimate(T)(in T[] sample)
	if(isFloatingPoint!T)
{
	return inverseGammaEstimate(InverseGammaStatistic!T(sample));
}

///
unittest {
	import atmosphere.estimate.generalized_gamma;
	immutable sample = [1.0, 0.5, 0.75];
	immutable p0 = inverseGammaEstimate(sample);
	immutable p1 = generalizedGammaFixedPowerEstimate(-1.0, sample);
	assert(p0.shape == p1.shape);
	assert(p0.scale == p1.scale);
}


///
unittest
{
	import std.range;
	import std.random;
	import atmosphere.random;
	import atmosphere.likelihood.inverse_gamma;
	auto length = 1000;
	auto shape = 2.0, scale = 3.0;
	auto rng = Random(1234);
	auto sample = InverseGammaSRNG!double(rng, shape, scale).take(length).array;
	auto weights = iota(1.0, length + 1.0).array;
	auto params = inverseGammaEstimate!double(sample, weights);
	auto lh0 = inverseGammaLikelihood(shape, scale, sample, weights);
	auto lh1 = inverseGammaLikelihood(params.shape, params.scale, sample, weights);
	assert(lh0 <= lh1);
}

///ditto
Tuple!(T, "shape", T, "scale")
inverseGammaEstimate(T)(in T[] sample, in T[] weights)
	if(isFloatingPoint!T)
{
	return inverseGammaEstimate(InverseGammaStatistic!T(sample, weights));
}

///
unittest {
	import atmosphere.estimate.generalized_gamma;
	immutable sample = [1.0, 0.5, 0.75];
	immutable weights = [1.0, 4, 3];
	immutable p0 = inverseGammaEstimate(sample, weights);
	immutable p1 = generalizedGammaFixedPowerEstimate(-1.0, sample, weights);
	assert(p0.shape == p1.shape);
	assert(p0.scale == p1.scale);
}

///ditto
Tuple!(T, "shape", T, "scale")
inverseGammaEstimate(T)(in InverseGammaStatistic!T stat)
	if(isFloatingPoint!T)
{
	import atmosphere.math: logmdigammaInverse;
	with(stat)
	{
		immutable T shape = logmdigammaInverse(log(meani) + meanl);
		immutable T scale = shape / meani;
		return typeof(return)(shape, scale);		
	}
}
