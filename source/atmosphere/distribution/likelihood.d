/++
Derivatives of probability density functions
+/
module distribution.likelihood;

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
T GammaLikelihood(T)(T shape, T scale, T a, T b)
	if(isFloatingPoint!T)
{
	return log(1 / (scale * gamma(shape))) 
		+ (shape - 1) * (b - log(scale)) 
		- a / scale;
}

/++
Log-likelihood function of the gamma distribution.
Params:
	shape = shape
	scale = scale
	sample = sample
+/
T GammaLikelihood(T)(T shape, T scale, in T[] sample)
	if(isFloatingPoint!T)
{
	T a = 0, b = 0;
	foreach(j; 0..sample.length)
	{
		a += sample[j];
		b += log2(sample[j]);
	}
	return GammaLikelihood!T(shape, scale, a, b);
}

unittest {
	assert(GammaLikelihood!double(1,3,[1,2]).isFinite);
}

/++
Log-likelihood function of the gamma distribution.
Params:
	shape = shape
	scale = scale
	sample = sample
	weights = weights for sample
+/
T GammaLikelihood(T)(T shape, T scale, in T[] sample, in T[] weights)
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
	return GammaLikelihood!T(shape, scale, a, T(LN2) * b);
}

unittest {
	assert(GammaLikelihood!double(1,3,[1,2],[2,3]).isFinite);
}


/++
Log-likelihood function of the inverse-gamma distribution.
Params:
	shape = shape
	scale = scale
	a = `Σ weights[j] / sample[j]`
	b = `Σ weights[j] * log(sample[j])`
+/
T InverseGammaLikelihood(T)(T shape, T scale, T a, T b)
	if(isFloatingPoint!T)
{
	return log(scale / gamma(shape))
		+ (shape - 1) * (b + log(scale)) 
		- a * scale;
}

/++
Log-likelihood function of the inverse-gamma distribution.
Params:
	shape = shape
	scale = scale
	sample = sample
+/
T InverseGammaLikelihood(T)(T shape, T scale, in T[] sample)
	if(isFloatingPoint!T)
{
	T a = 0, b = 0;
	foreach(j; 0..sample.length)
	{
		a += 1 / sample[j];
		b += log2(sample[j]);
	}
	return InverseGammaLikelihood!T(shape, scale, a, b);
}

unittest {
	assert(InverseGammaLikelihood!double(1,3,[1,2]).isFinite);
}

/++
Log-likelihood function of the inverse-gamma distribution.
Params:
	shape = shape
	scale = scale
	sample = sample
	weights = weights for sample
+/
T InverseGammaLikelihood(T)(T shape, T scale, in T[] sample, in T[] weights)
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
	return InverseGammaLikelihood!T(shape, scale, a, T(LN2) * b);
}

unittest {
	assert(InverseGammaLikelihood!double(1,3,[1,2],[2,3]).isFinite);
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
T GeneralizedGammaLikelihood(T)(T shape, T power, T scale, T a, T b)
	if(isFloatingPoint!T)
{
	return log(fabs(power) / (scale * gamma(shape))) 
		+ (shape * power - 1) * (b - log(scale)) 
		- a / pow(scale, power);
}

/++
Log-likelihood function of the generalized gamma distribution.
Params:
	shape = shape
	power = power
	scale = scale
	sample = sample
+/
T GeneralizedGammaLikelihood(T)(T shape, T power, T scale, in T[] sample)
	if(isFloatingPoint!T)
{
	T a = 0, b = 0;
	foreach(j; 0..sample.length)
	{
		a += pow(sample[j], power);
		b += log2(sample[j]);
	}
	return GeneralizedGammaLikelihood!T(shape, power, scale, a, b);
}

unittest {
	assert(GeneralizedGammaLikelihood!double(1,2,3,[1,2]).isFinite);
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
T GeneralizedGammaLikelihood(T)(T shape, T power, T scale, in T[] sample, in T[] weights)
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
	return GeneralizedGammaLikelihood!T(shape, power, scale, a, T(LN2) * b);
}

unittest {
	assert(GeneralizedGammaLikelihood!double(1,2,3,[1,2],[2,3]).isFinite);
}


/++
Log-likelihood function of the generalized inverse Gaussian distribution.
Params:
	lambda =
	eta = scale
	omega =
	a = `Σ weights[j] * sample[j]`
	b = `Σ weights[j] / sample[j]`
	c = `Σ weights[j] * log(sample[j])`
See_Also: distribution.params.GIGEtaOmega
+/
T GeneralizedInverseGaussianLikelihood(T)(T lambda, T eta, T omega, T a, T b, T c)
	if(isFloatingPoint!T)
{
	import bessel;
	import std.typecons : Flag;
	return 
		- log(2*eta*besselK(omega, lambda, Flag!"ExponentiallyScaled".yes)) / omega
		+ (lambda - 1) * (c - log(eta))
		- 0.5f * omega * (b * eta + a / eta);
}


/++
Log-likelihood function of the generalized inverse Gaussian distribution.
Params:
	lambda =
	eta = scale
	omega =
	sample = sample
See_Also: `distribution.params.GIGEtaOmega`
+/
T GeneralizedInverseGaussianLikelihood(T)(T lambda, T eta, T omega, in T[] sample)
	if(isFloatingPoint!T)
{
	T a = 0, b = 0, c = 0;
	foreach(j; 0..sample.length)
	{
		immutable s = sample[j];
		a += s;
		b += 1 / s;
		c += log2(s);
	}
	return GeneralizedInverseGaussianLikelihood!T(lambda, eta, omega, a, b, T(LN2) * c);
}

unittest {
	assert(GeneralizedInverseGaussianLikelihood!double(1,2,3,[1,2]).isFinite);
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
T GeneralizedInverseGaussianLikelihood(T)(T lambda, T eta, T omega, in T[] sample, in T[] weights)
	if(isFloatingPoint!T)
in {
	assert(weights.length == sample.length);
}
body {
	T a = 0, b = 0, c = 0;
	foreach(j; 0..sample.length)
	{
		immutable w = weights[j];
		immutable s = sample[j];
		a += w * s;
		b += w / s;
		c += w * log2(s);
	}
	return GeneralizedInverseGaussianLikelihood!T(lambda, eta, omega, a, b, T(LN2) * c);
}

unittest {
	assert(GeneralizedInverseGaussianLikelihood!double(1,2,3,[1,2],[2,3]).isFinite);
}
