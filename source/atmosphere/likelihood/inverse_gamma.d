/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.likelihood.inverse_gamma;

import core.stdc.tgmath;

import std.traits;
import std.typecons;

import atmosphere.statistic: InverseGammaStatistic;

/++
Normalized log-likelihood function of the inverse-gamma distribution.
+/
T inverseGammaLikelihood(T)(T shape, T scale, in T[] sample)
	if(isFloatingPoint!T)
{
	return inverseGammaLikelihood!T(shape, scale, InverseGammaStatistic!T(sample));
}

///
unittest {
	import atmosphere.likelihood.generalized_gamma;
	immutable l = inverseGammaLikelihood!double(2,3,[1,2]);
	immutable m = generalizedGammaLikelihood!double(2,-1,3,[1,2]);
	assert(l == m);
}

///ditto
T inverseGammaLikelihood(T)(T shape, T scale, in T[] sample, in T[] weights)
	if(isFloatingPoint!T)
{
	return inverseGammaLikelihood!T(shape, scale, InverseGammaStatistic!T(sample, weights));
}

///
unittest {
	import atmosphere.likelihood.generalized_gamma;
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


///ditto
T inverseGammaLikelihood(T)(T shape, T scale, InverseGammaStatistic!T stat)
	if(isFloatingPoint!T)
{
	with(stat) return 
		- log(scale * tgamma(shape))
		- (1 + shape) * (meanl - log(scale)) 
		- meani * scale;
}
