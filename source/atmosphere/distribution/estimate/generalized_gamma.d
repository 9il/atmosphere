/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.distribution.estimate.generalized_gamma;

import core.stdc.tgmath;

import std.traits;
import std.typecons;
import std.math : LN2;


/++
Estimates parameters of the generalized gamma distribution.
Params:
	power = power parameter
	sample = sample
+/
Tuple!(T, "shape", T, "power", T, "scale")
generalizedGammaEstimate(T)(T power, in T[] sample)
	if(isFloatingPoint!T)
in {
	import std.algorithm : all;
	assert(sample.all!"a > 0 && isNormal(a)");
}
body {
	import std.algorithm : sum;
	Unqual!T a = 0, b = 0, c = 0;
	foreach(x; sample)
	{
		immutable l = log2(x);
		immutable y = x.pow(power);
		a += y;
		b += l * y;
		c += l;
	}
	b *= T(LN2);
	c *= T(LN2);
	immutable n = sample.length;
	return generalizedGammaEstimate(power, a/n, b/n, c/n);
}


/++
Estimates parameters of the generalized gamma distribution.
Params:
	power = power parameter
	sample = sample
	weights = weights
+/
Tuple!(T, "shape", T, "power", T, "scale")
generalizedGammaEstimate(T)(T power, in T[] sample, in T[] weights)
	if(isFloatingPoint!T)
in {
	import std.algorithm : all, any;
	assert(weights.length == sample.length);
	assert(sample.all!"a > 0 && isNormal(a)");
	assert(weights.all!"a >= 0 && isFinite(a)");
	assert(weights.any!"a > 0");
}
body {
	import std.algorithm : sum;
	Unqual!T a = 0, b = 0, c = 0, n = 0;
	foreach(i, x; sample)
	{
		immutable w = weights[i];
		immutable l = log2(x);
		immutable y = x.pow(power);
		n += w;
		a += w * y;
		b += w * (l * y);
		c += w * l;
	}
	b *= T(LN2);
	c *= T(LN2);
	return generalizedGammaEstimate(power, a/n, b/n, c/n);
}


/++
Estimates parameters of the generalized gamma distribution.
Params:
	power = power
	a = `Σ weights[j] * sample[j] ^^ power`
	b = `Σ weights[j] * log(sample[j]) * sample[j] ^^ power`
	c = `Σ weights[j] * log(sample[j])`
+/
Tuple!(T, "shape", T, "power", T, "scale")
generalizedGammaEstimate(T)(T power, T a, T b, T c)
	if(isFloatingPoint!T)
{
	immutable d = power * (b / a - c);
	return typeof(return)(1 / d, power, pow(a * d, 1 / power));
}
