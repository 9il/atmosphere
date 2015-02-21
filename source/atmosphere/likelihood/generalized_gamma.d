/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.likelihood.generalized_gamma;

import core.stdc.tgmath;

import std.traits;
import std.typecons;

import atmosphere.statistic: GeneralizedGammaFixedPowerStatistic;

/++
Normalized log-likelihood function of the generalized gamma distribution.
+/
T generalizedGammaLikelihood(T)(T shape, T power, T scale, in T[] sample)
	if(isFloatingPoint!T)
{
	return generalizedGammaLikelihood!T(shape, power, scale, GeneralizedGammaFixedPowerStatistic!T(power, sample));
}

///ditto
T generalizedGammaLikelihood(T)(T shape, T power, T scale, in T[] sample, in T[] weights)
	if(isFloatingPoint!T)
{
	return generalizedGammaLikelihood!T(shape, power, scale, GeneralizedGammaFixedPowerStatistic!T(power, sample, weights));
}

///
unittest {
	immutable l = generalizedGammaLikelihood!double(1,2,3,[1,2,2]);
	immutable m = generalizedGammaLikelihood!double(1,2,3,[1,2],[2,4]);
	assert(l == m);
}

///ditto
T generalizedGammaLikelihood(T)(T shape, T power, T scale, GeneralizedGammaFixedPowerStatistic!T stat)
	if(isFloatingPoint!T)
{
	with(stat) return 
		- log((scale * tgamma(shape)) / fabs(power)) 
		- (1 - shape * power) * (meanl - log(scale)) 
		- (power > 0 ? meanp / pow(scale, power) : meanp * pow(scale, -power)); //precise unification with inverse-gamma and gamma
}
