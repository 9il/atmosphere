/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.distribution.estimate.gamma;

import core.stdc.tgmath;

import std.traits;
import std.typecons;
import std.math : LN2;


/++
Estimates parameters of the gamma distribution.
+/
Tuple!(T, "shape", T, "scale")
gammaEstimate(T)(in T[] sample)
	if(isFloatingPoint!T)
in {
	import std.algorithm : all;
	assert(sample.all!"a > 0 && isNormal(a)");
}
body {
	Unqual!T a = 0, b = 0, c = 0;
	foreach(x; sample)
	{
		immutable l = log2(x);
		a += x;
		b += l * x;
		c += l;
	}
	b *= T(LN2);
	c *= T(LN2);
	immutable n = sample.length;
	return gammaEstimate(a/n, b/n, c/n);
}

///
unittest {
	import atmosphere.distribution.estimate.generalized_gamma;
	immutable sample = [1.0, 0.5, 0.75];
	immutable p0 = gammaEstimate(sample);
	immutable p1 = generalizedGammaEstimate(1.0, sample);
	assert(p0.shape == p1.shape);
	assert(p0.scale == p1.scale);
}


///ditto
Tuple!(T, "shape", T, "scale")
gammaEstimate(T)(in T[] sample, in T[] weights)
	if(isFloatingPoint!T)
in {
	import std.algorithm : all, any;
	assert(weights.length == sample.length);
	assert(sample.all!"a > 0 && isNormal(a)");
	assert(weights.all!"a >= 0 && isFinite(a)");
	assert(weights.any!"a > 0");
}
body {
	Unqual!T a = 0, b = 0, c = 0, n = 0;
	foreach(i, x; sample)
	{
		immutable w = weights[i];
		immutable l = log2(x);
		n += w;
		a += w * x;
		b += w * (l * x);
		c += w * l;
	}
	b *= T(LN2);
	c *= T(LN2);
	return gammaEstimate(a/n, b/n, c/n);
}

///
unittest {
	import atmosphere.distribution.estimate.generalized_gamma;
	immutable sample = [1.0, 0.5, 0.75];
	immutable weights = [1.0, 4, 3];
	immutable p0 = gammaEstimate(sample, weights);
	immutable p1 = generalizedGammaEstimate(1.0, sample, weights);
	assert(p0.shape == p1.shape);
	assert(p0.scale == p1.scale);
}


/++
Estimates parameters of the inverse-gamma distribution.
Params:
	a = `Σ weights[j] * sample[j] / Σ weights[j]`
	b = `Σ weights[j] * log(sample[j]) * sample[j] / Σ weights[j]`
	c = `Σ weights[j] * log(sample[j]) / Σ weights[j]`
+/
Tuple!(T, "shape", T, "scale")
gammaEstimate(T)(T a, T b, T c)
	if(isFloatingPoint!T)
{
	immutable d = b / a - c;
	return typeof(return)(1 / d, a * d);
}

///
unittest {
	import atmosphere.distribution.estimate.generalized_gamma;
	immutable p0 = gammaEstimate(3.0, 2.0, 1.0);
	immutable p1 = generalizedGammaEstimate(1.0, 3.0, 2.0, 1.0);
	assert(p0.shape == p1.shape);
	assert(p0.scale == p1.scale);
}
