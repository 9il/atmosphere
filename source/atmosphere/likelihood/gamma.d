/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.likelihood.gamma;

import core.stdc.tgmath;

import std.traits;
import std.typecons;

import atmosphere.statistic: GammaStatistic;


/++
Normalized log-likelihood function of the gamma distribution.
+/
T gammaLikelihood(T)(T shape, T scale, in T[] sample)
	if(isFloatingPoint!T)
{
	return gammaLikelihood!T(shape, scale, GammaStatistic!T(sample));
}

///
unittest {
	import atmosphere.likelihood.generalized_gamma;
	immutable l = gammaLikelihood!double(2,3,[1,2]);
	immutable m = generalizedGammaLikelihood!double(2,1,3,[1,2]);
	assert(l == m);
}

///ditto
T gammaLikelihood(T)(T shape, T scale, in T[] sample, in T[] weights)
	if(isFloatingPoint!T)
{
	return gammaLikelihood!T(shape, scale, GammaStatistic!T(sample, weights));
}

///
unittest {
	import atmosphere.likelihood.generalized_gamma;
	immutable l = gammaLikelihood!double(2,3,[1,2],[3,4]);
	immutable m = generalizedGammaLikelihood!double(2,1,3,[1,2],[3,4]);
	assert(l == m);
}

///
unittest {
	immutable l = gammaLikelihood!double(1,3,[1,2,2]);
	immutable m = gammaLikelihood!double(1,3,[1,2],[2,4]);
	assert(l == m);
}

///ditto
T gammaLikelihood(T)(T shape, T scale, GammaStatistic!T stat)
	if(isFloatingPoint!T)
{
	with(stat) return 
		- log(scale * tgamma(shape)) 
		- (1 - shape) * (meanl - log(scale)) 
		- mean / scale;
}
