module atmosphere.distribution.likelihood;

import std.traits;
import std.mathspecial;


/++
Log-likelihood function of the gamma distribution.
Params:
	shape = shape
	scale = scale
	a = `Σ weights[j] * sample[j]`
	b = `Σ weights[j] * log(sample[j])`
+/
T gammaLikelihood(T)(T shape, T scale, T a, T b)
	if(isFloatingPoint!T)
{
	return 
		- log(scale * gamma(shape)) 
		- (1 - shape) * (b - log(scale)) 
		- a / scale;
}

///
unittest {
	immutable l = gammaLikelihood(4.0, 3.0, 2.0, 1.0);
	immutable m = generalizedGammaLikelihood(4.0, 1.0, 3.0, 2.0, 1.0);
	assert(l == m);
}


/++
Log-likelihood function of the gamma distribution.
Params:
	shape = shape
	scale = scale
	sample = sample
+/
T gammaLikelihood(T)(T shape, T scale, in T[] sample)
	if(isFloatingPoint!T)
{
	T a = 0, b = 0;
	foreach(j; 0..sample.length)
	{
		a += sample[j];
		b += log2(sample[j]);
	}
	return gammaLikelihood!T(shape, scale, a, T(LN2) * b);
}

///
unittest {
	immutable l = gammaLikelihood!double(1,3,[1,2]);
	immutable m = gammaLikelihood!double(1,3,[1,2],[1,1]);
	assert(l == m);
}


/++
Log-likelihood function of the gamma distribution.
Params:
	shape = shape
	scale = scale
	sample = sample
	weights = weights for sample
+/
T gammaLikelihood(T)(T shape, T scale, in T[] sample, in T[] weights)
	if(isFloatingPoint!T)
in {
	assert(weights.length == sample.length);
}
body {
	T a = 0, b = 0;
	foreach(j; 0..sample.length)
	{
		immutable w = weights[j];
		a += w * sample[j];
		b += w * log2(sample[j]);
	}
	return gammaLikelihood!T(shape, scale, a, T(LN2) * b);
}

///
unittest {
	immutable l = gammaLikelihood!double(2,3,[1,2],[3,4]);
	immutable m = generalizedGammaLikelihood!double(2,1,3,[1,2],[3,4]);
	assert(l == m);
}


/++
Log-likelihood function of the inverse-gamma distribution.
Params:
	shape = shape
	scale = scale
	a = `Σ weights[j] / sample[j]`
	b = `Σ weights[j] * log(sample[j])`
+/
T inverseGammaLikelihood(T)(T shape, T scale, T a, T b)
	if(isFloatingPoint!T)
{
	return 
		- log(scale * gamma(shape))
		- (1 + shape) * (b - log(scale)) 
		- a * scale;
}

///
unittest {
	immutable l = inverseGammaLikelihood(4.0, 3.0, 2.0, 1.0);
	immutable m = generalizedGammaLikelihood(4.0, -1.0, 3.0, 2.0, 1.0);
	assert(l == m);
}


/++
Log-likelihood function of the inverse-gamma distribution.
Params:
	shape = shape
	scale = scale
	sample = sample
+/
T inverseGammaLikelihood(T)(T shape, T scale, in T[] sample)
	if(isFloatingPoint!T)
{
	T a = 0, b = 0;
	foreach(j; 0..sample.length)
	{
		a += 1 / sample[j];
		b += log2(sample[j]);
	}
	return inverseGammaLikelihood!T(shape, scale, a, T(LN2) * b);
}

///
unittest {
	immutable l = inverseGammaLikelihood!double(1,3,[1,2]);
	immutable m = inverseGammaLikelihood!double(1,3,[1,2],[1,1]);
	assert(l == m);
}


/++
Log-likelihood function of the inverse-gamma distribution.
Params:
	shape = shape
	scale = scale
	sample = sample
	weights = weights for sample
+/
T inverseGammaLikelihood(T)(T shape, T scale, in T[] sample, in T[] weights)
	if(isFloatingPoint!T)
in {
	assert(weights.length == sample.length);
}
body {
	T a = 0, b = 0;
	foreach(j; 0..sample.length)
	{
		immutable w = weights[j];
		a += w / sample[j];
		b += w * log2(sample[j]);
	}
	return inverseGammaLikelihood!T(shape, scale, a, T(LN2) * b);
}

///
unittest {
	immutable l = inverseGammaLikelihood!double(2,3,[1,2],[3,4]);
	immutable m = generalizedGammaLikelihood!double(2,-1,3,[1,2],[3,4]);
	assert(l == m);
}


/++
Log-likelihood function of the generalized gamma distribution.
Params:
	shape = shape
	power = power
	scale = scale
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


/++
Log-likelihood function of the generalized gamma distribution.
Params:
	shape = shape
	power = power
	scale = scale
	sample = sample
+/
T generalizedGammaLikelihood(T)(T shape, T power, T scale, in T[] sample)
	if(isFloatingPoint!T)
{
	T a = 0, b = 0;
	foreach(j; 0..sample.length)
	{
		a += pow(sample[j], power);
		b += log2(sample[j]);
	}
	return generalizedGammaLikelihood!T(shape, power, scale, a, T(LN2) * b);
}

///
unittest {
	immutable l = generalizedGammaLikelihood!double(1,2,3,[1,2]);
	immutable m = generalizedGammaLikelihood!double(1,2,3,[1,2],[1,1]);
	assert(l == m);
}


/++
Log-likelihood function of the generalized gamma distribution.
Params:
	shape = shape
	power = power
	scale = scale
	sample = sample
	weights = weights for sample
+/
T generalizedGammaLikelihood(T)(T shape, T power, T scale, in T[] sample, in T[] weights)
	if(isFloatingPoint!T)
in {
	assert(weights.length == sample.length);
}
body {
	T a = 0, b = 0;
	foreach(j; 0..sample.length)
	{
		immutable w = weights[j];
		a += w * pow(sample[j], power);
		b += w * log2(sample[j]);
	}
	return generalizedGammaLikelihood!T(shape, power, scale, a, T(LN2) * b);
}


///++
//Log-likelihood function of the generalized inverse Gaussian distribution.
//Params:
//	lambda =
//	eta = scale
//	omega =
//	a = `Σ weights[j] * sample[j]`
//	b = `Σ weights[j] / sample[j]`
//	c = `Σ weights[j] * log(sample[j])`
//See_Also: distribution.params.GIGEtaOmega
//+/
//T generalizedInverseGaussianLikelihood(T)(T lambda, T eta, T omega, T a, T b, T c)
//	if(isFloatingPoint!T)
//{
//	import bessel;
//	import std.typecons : Flag;
//	return 
//		- log(2 * eta * besselK(omega, lambda, Flag!"ExponentiallyScaled".yes)) / omega
//		+ (lambda - 1) * (c - log(eta))
//		- 0.5f * omega * (b * eta + a / eta);
//}


/++
Log-likelihood function of the generalized inverse Gaussian distribution.
Params:
	lambda =
	eta = scale
	omega =
	sample = sample
See_Also: `distribution.params.GIGEtaOmega`
+/
T generalizedInverseGaussianLikelihood(T)(T lambda, T eta, T omega, in T[] sample)
	if(isFloatingPoint!T)
{
	import bessel;
	import std.typecons : Flag;
	import atmosphere.utilities;
	import std.algorithm : map, sum;
	auto normolizedSample = sample.map!(x => x / eta);
	immutable n = sample.length;
	T a = normolizedSample.map!"a + 1/a".sum / n;
	T c = T(LN2) * normolizedSample.sumOfLog2s() / n;
	return 
		- log(2 * besselK(omega, lambda, Flag!"ExponentiallyScaled".yes)) / omega
		+ (lambda - 1) * c
		- 0.5f * omega * a;
}

///
unittest {
	immutable l = generalizedInverseGaussianLikelihood!double(1,2,3,[1,2]);
	immutable m = generalizedInverseGaussianLikelihood!double(1,2,3,[1,2],[1,1]);
//	assert(l == m);
}


/++
Log-likelihood function of the generalized inverse Gaussian distribution.
Params:
	lambda =
	eta = scale
	omega =
	sample = sample
	weights = weights for sample
See_Also: `distribution.params.GIGEtaOmega`
+/
T generalizedInverseGaussianLikelihood(T)(T lambda, T eta, T omega, in T[] sample, in T[] weights)
	if(isFloatingPoint!T)
in {
	assert(weights.length == sample.length);
}
body {
	import bessel;
	import std.typecons : Flag;
	import atmosphere.utilities;
	import std.algorithm : map, sum;
	auto normolizedSample = sample.map!(x => x / eta);
	immutable n = weights.sum();
	T a = normolizedSample.map!"a + 1/a".dotProduct(weights) / n;
	T c = T(LN2) * normolizedSample.map!log2.dotProduct(weights) / n;
	return 
		- log(2 * eta * besselK(omega, lambda, Flag!"ExponentiallyScaled".yes)) / omega
		+ (lambda - 1) * c
		- 0.5f * omega * a;
}


private T logBesselKParamDerivative(T)(T x, T alpha)
{
	import bessel;
	import std.typecons;
	import scid.calculus;
	//T f(T t) {
	//	return exp(-x * cosh(t)) * sinh(alpha * t) * t;
	//}
	//auto r = integrate(&f, T(0), T.infinity);
	//return r.value;
	return diff((T alpha) { return log(besselK(x, alpha, Flag!"ExponentiallyScaled".yes)) - x; }, alpha, T(1), 16).value;

}

unittest {
	auto l = logBesselKParamDerivative!real(.60L, -2.0L);
	assert(isNormal(l));
}

private T logBesselKDerivative(T)(T x, T alpha)
{
	import bessel;
	import std.typecons;
	immutable a = besselK(x, alpha-1, Flag!"ExponentiallyScaled".yes);
	immutable b = besselK(x, alpha+1, Flag!"ExponentiallyScaled".yes);
	immutable c = besselK(x, alpha  , Flag!"ExponentiallyScaled".yes);
	return -0.5f*(a+b)/c;
}

unittest {
	auto l = logBesselKDerivative!real(.60L, -2.0L);
	assert(isNormal(l));
	import std.stdio;
	//writeln(l)
	import bessel;
	writeln(besselK(1, 1));
	writeln(besselK(2, 0));
}
