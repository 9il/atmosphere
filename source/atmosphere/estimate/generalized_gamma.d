/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.estimate.generalized_gamma;

import core.stdc.tgmath;

import std.traits;
import std.typecons;

import atmosphere.statistic: GeneralizedGammaFixedPowerStatistic;


/++
Estimates parameters of the generalized gamma distribution.
+/
Tuple!(T, "shape", T, "scale")
generalizedGammaFixedPowerEstimate(T)(T power, in T[] sample)
	if(isFloatingPoint!T)
{
	return generalizedGammaFixedPowerEstimate(power, GeneralizedGammaFixedPowerStatistic!T(power, sample));
}

///
unittest
{
	import std.range;
	import std.random;
	import atmosphere.random;
	import atmosphere.likelihood;
	size_t length = 1000;
	auto shape = 2.0, power = 1.4, scale = 2.3;
	auto rng = Random(1234);
	auto sample = GeneralizedGammaSRNG!double(rng, shape, power, scale).take(length).array;
	auto params = generalizedGammaFixedPowerEstimate!double(power, sample);
	auto lh0 = generalizedGammaLikelihood(shape, power, scale, sample);
	auto lh1 = generalizedGammaLikelihood(params.shape, power, params.scale, sample);
	assert(lh0 <= lh1);
}


///ditto
Tuple!(T, "shape", T, "scale")
generalizedGammaFixedPowerEstimate(T)(T power, in T[] sample, in T[] weights)
	if(isFloatingPoint!T)
{
	return generalizedGammaFixedPowerEstimate(power, GeneralizedGammaFixedPowerStatistic!T(power, sample, weights));
}

///
unittest
{
	import std.range;
	import std.random;
	import atmosphere.random;
	import atmosphere.likelihood;
	size_t length = 1000;
	auto shape = 2.0, power = 1.4, scale = 2.3;
	auto rng = Random(1234);
	auto sample = GeneralizedGammaSRNG!double(rng, shape, power, scale).take(length).array;
	auto weights = iota(1.0, length + 1.0).array;
	auto params = generalizedGammaFixedPowerEstimate!double(power, sample, weights);
	auto lh0 = generalizedGammaLikelihood(shape, power, scale, sample, weights);
	auto lh1 = generalizedGammaLikelihood(params.shape, power, params.scale, sample, weights);
	assert(lh0 <= lh1);
}


///ditto
Tuple!(T, "shape", T, "scale")
generalizedGammaFixedPowerEstimate(T)(T power, GeneralizedGammaFixedPowerStatistic!T stat)
	if(isFloatingPoint!T)
{
	import atmosphere.math: logmdigammaInverse;
	with(stat)
	{
		immutable T shape = logmdigammaInverse(log(meanp) - power * meanl);
		immutable T scale = pow(meanp / shape, 1 / power);
		return typeof(return)(shape, scale);		
	}
}
