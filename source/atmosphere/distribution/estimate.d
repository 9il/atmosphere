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

enum k01 = 0.421024438240708333335627379212609036136219748226660472298969L;

//struct GIGEstParamsFunction(T)
//{
//	private T mone, dzero, a, c, s;
//	this(T one, T mone, T dzero)
//	{

//		this.mone = mone;
//		this.dzero = dzero;
//		a = one * mone;
//		c = a / (a - 1);
//		s = sqrt(a) - c + 1; //Flag!"ExponentiallyScaled".yes

//	}
//	T opCall(T tau)
//	{
//		import bessel;
//		immutable e = sqrt(tau^^2 + a);
//		immutable omega = c / e;
//		immutable lambda = tau * omega;
//		immutable eta = (e - tau) / mone;
//		return 
//			besselK(omega, lambda, Flag!"ExponentiallyScaled".yes)
//			- k01 * exp(s + lambda*(dzero - log(eta)));
//	}
//}

struct GIGEstParamsFunction(T)
{
	private T one, mone, dzero, a, c, s;
	this(T one, T mone, T dzero)
	{

		this.one = one;
		this.mone = mone;
		this.dzero = dzero;
		a = one * mone;
		c = a / (a - 1);
		s = sqrt(a) - c + 1; //Flag!"ExponentiallyScaled".yes

	}
	T opCall(T lambda)
	{
		import bessel;
		immutable omega = sqrt((c^^2-lambda^^2)/a);
		immutable eta = sqrt((one*(c-lambda))/(mone*(c+lambda)));
		return 
			besselK(omega, lambda, Flag!"ExponentiallyScaled".yes)
			- k01 * exp(s + lambda*(dzero - log(eta)));
	}
}


/++
Estimates parameters of the generalized inverse Gaussian distribution
+/
//Tuple!(T, "lambda", T, "eta", T, "omega")
//generalizedInverseGaussian(T)(T one, T mone, T dzero, T tauLeft, T tauRight)
//	if(isFloatingPoint!T)
//{
//	import std.numeric;
//	import bessel;
//	immutable a = one * mone;
//	immutable c = a / (a - 1);
//	immutable s = sqrt(a) - c + 1; //Flag!"ExponentiallyScaled".yes
//	T f(T tau) {
//		immutable e = sqrt(tau^^2 + a);
//		immutable omega = c / e;
//		immutable lambda = tau * omega;
//		immutable eta = (e - tau) / mone;
//		return 
//			besselK(omega, lambda, Flag!"ExponentiallyScaled".yes)
//			- k01 * exp(s + lambda*(dzero - log(eta)));
//	}
//	import std.stdio;
//	writeln("1, ", f(tauLeft));
//	writeln("2, ", f(tauRight));
//	immutable tau = findRoot(&f, tauLeft, tauRight);
//	immutable e = sqrt(tau^^2 + a);
//	immutable omega = c / e;
//	immutable lambda = tau * omega;
//	immutable eta = (e - tau) / mone;
//	return typeof(return)(lambda, eta, omega);
//}

Tuple!(T, "lambda", T, "eta", T, "omega")
generalizedInverseGaussian(T)(T one, T mone, T dzero, T tauLeft, T tauRight)
	if(isFloatingPoint!T)
{
	import std.numeric;
	import bessel;
	immutable a = one * mone;
	immutable c = a / (a - 1);
	immutable s = sqrt(a) - c + 1; //Flag!"ExponentiallyScaled".yes
	T f(T lambda) {
		immutable omega = sqrt((c^^2-lambda^^2)/a);
		immutable eta = sqrt((one*(c-lambda))/(mone*(c+lambda)));
		return 
			besselK(omega, lambda, Flag!"ExponentiallyScaled".yes)
			- k01 * exp(s + lambda*(dzero - log(eta)));
	}
	import std.stdio;
	writeln("1, ", f(tauLeft));
	writeln("2, ", f(tauRight));
	immutable lambda = findRoot(&f, tauLeft, tauRight);
	immutable omega = sqrt((c^^2-lambda^^2)/a);
	immutable eta = sqrt((one*(c-lambda))/(mone*(c+lambda)));
	return typeof(return)(lambda, eta, omega);
}

unittest {
	alias f = generalizedInverseGaussian!double;
	import std.stdio;
	import std.range;
	auto fp = GIGEstParamsFunction!double(4, 2, -1);
	foreach(e; iota(-8.0/7, 8.0/7, 0.01))
	{
		fp(e).writeln(" ", e);
	}
	auto l = f(4, 2, -1, -8.0/7*0.99, 8.0/7*0.99);
	writeln(l);
}



unittest {
	import atmosphere;
	import std.random;
	import std.algorithm;
	import std.range;
	import std.stdio;
	double lambda = 2;
	double eta = 1;
	double omega = 3;
	auto rng = ProperGeneralizedInverseGaussianSRNG!double(rndGen, lambda, eta, omega);
	auto sample = rng.take(100000).array;
	double one = sample.sum() / sample.length;
	double mone = sample.map!"1/a".sum() / sample.length;
	double dzero = double(LN2) * sample.sumOfLog2s() / sample.length;
	immutable a = one * mone;
	immutable c = a / (a - 1);
	auto fp = GIGEstParamsFunction!double(one, mone, dzero);
	foreach(e; iota(-c, c, 0.01))
	{
		fp(e).writeln(" ", e);
	}
	writefln("one= %s  mone= %s  dzero= %s  a= %s  c= %s  ", one, mone, dzero, a, c);
	auto params = generalizedInverseGaussian!double(one, mone, dzero, 0, c*0.99);
	writeln(params);
	writeln(generalizedInverseGaussianLikelihood(lambda, eta, omega, sample));
	writeln(generalizedInverseGaussianLikelihood(params.lambda, params.eta, params.omega, sample));
}