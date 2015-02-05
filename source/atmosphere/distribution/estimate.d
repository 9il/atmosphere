module atmosphere.distribution.estimate;

import atmosphere.distribution.likelihood;

import std.traits;
import std.typecons;
import std.mathspecial;

/++
Estimates parameters of the gamma distribution.
Params:
	a = `Σ weights[j] * sample[j]`
	b = `Σ weights[j] * log(sample[j]) * sample[j]`
	c = `Σ weights[j] * log(sample[j])`
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
	immutable p0 = gammaEstimate(3.0, 2.0, 1.0);
	immutable p1 = generalizedGammaEstimate(1.0, 3.0, 2.0, 1.0);
	assert(p0.shape == p1.shape);
	assert(p0.scale == p1.scale);
}

/++
Estimates parameters of the inverse-gamma distribution.
Params:
	a = `Σ weights[j] * sample[j]`
	b = `Σ weights[j] * log(sample[j]) / sample[j]`
	c = `Σ weights[j] * log(sample[j])`
+/
Tuple!(T, "shape", T, "scale")
inverseGammaEstimate(T)(T a, T b, T c)
	if(isFloatingPoint!T)
{
	immutable d = c - b / a;
	return typeof(return)(1 / d, 1 / (a * d));
}

///
unittest {
	immutable p0 = inverseGammaEstimate(3.0, 2.0, 1.0);
	immutable p1 = generalizedGammaEstimate(-1.0, 3.0, 2.0, 1.0);
	assert(p0.shape == p1.shape);
	assert(p0.scale == p1.scale);
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
