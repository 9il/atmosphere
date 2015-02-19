/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.distribution.estimate.generalized_inverse_gaussian;

import core.stdc.tgmath;

import std.traits;
import std.typecons;
import std.math : LN2, signbit, isNormal;
import atmosphere.utilities : sumOfLog2s;

/++
Estimates parameters of the generalized inverse Gaussian distribution.
Params:
	sample = sample
	lambdaBounds = bounds of the parameter lambda
	omegaBounds = bounds of the parameter omega

See_Also: `distribution.params.GIGEtaOmega`
+/
Tuple!(T, "lambda", T, "eta", T, "omega") 
generalizedInverseGaussianEstimate(T)
	(in T[] sample, T[2] lambdaBounds = [-25, 25], T[2] omegaBounds = [1e-10, 1e10])
	if(isFloatingPoint!T)
in {
	import std.algorithm : all;
	assert(sample.all!"a > 0 && isNormal(a)");
}
body {
	return generalizedInverseGaussianEstimate!(T)(lambdaBounds, (T lambda, T one, T mone, T dzero){ return omegaBounds; }, sample);
}

unittest {
 	import atmosphere.distribution.likelihood;
	immutable lambda = -14.0;
	immutable eta = 1.0;
	immutable omega = 30.0;
	immutable one = 0x1.68ca419d651acp-1;
	immutable mone = 0x1.74ef782e761d4p+0;
	immutable dzero = -0x1.73d7daf81500dp-2;
	immutable params = generalizedInverseGaussianEstimate(lambda, [1e-10, 1e10], one, mone);
	immutable lhPrior = generalizedInverseGaussianLikelihood(lambda, eta, omega, one, mone, dzero);
	immutable lhCalc  = generalizedInverseGaussianLikelihood(params.lambda, params.eta, params.omega, one, mone, dzero);
	assert(lhPrior <= lhCalc);
}

/++
Estimates parameters of the generalized inverse Gaussian distribution.
Params:
	sample = sample
	weights = weights for the sample
	lambdaBounds = bounds of the parameter lambda
	omegaBounds = bounds of the parameter omega

See_Also: `distribution.params.GIGEtaOmega`
+/
Tuple!(T, "lambda", T, "eta", T, "omega") 
generalizedInverseGaussianEstimate(T)
	(in T[] sample, in T[] weights, T[2] lambdaBounds = [-25, 25], T[2] omegaBounds = [1e-10, 1e10])
	if(isFloatingPoint!T)
in {
	import std.algorithm : all, any;
	assert(weights.length == sample.length);
	assert(sample.all!"a > 0 && isNormal(a)");
	assert(weights.all!"a >= 0 && isFinite(a)");
	assert(weights.any!"a > 0");
}
body {
	return generalizedInverseGaussianEstimate!((T lambda, T one, T mone, T dzero){ return omegaBounds; }, T)(lambdaBounds, sample, weights);
}

/++
Estimates parameters of the generalized inverse Gaussian distribution.
Params:
	one = `Σ weights[j] * sample[j] / Σ weights[j]`
	mone = `Σ weights[j] / sample[j] / Σ weights[j]`
	dzero = `Σ weights[j] * log(sample[j]) / Σ weights[j]`
	lambdaBounds = bounds of the parameter lambda
	omegaBounds = bounds of the parameter omega

See_Also: `distribution.params.GIGEtaOmega`
+/
Tuple!(T, "lambda", T, "eta", T, "omega") 
generalizedInverseGaussianEstimate(T)
	(T one, T mone, T dzero, T[2] lambdaBounds = [-25, 25], T[2] omegaBounds = [1e-10, 1e10])
	if(isFloatingPoint!T)
{
	return generalizedInverseGaussianEstimate!((T lambda, T one, T mone, T dzero){ return omegaBounds; }, T)(lambdaBounds, one, mone, dzero);
}

/++
Estimates parameters of the generalized inverse Gaussian distribution.
Params:
	lambdaBounds = bounds of the parameter lambda
	omegaBoundsFun = function to calculate bounds of the parameter omega
	sample = sample

See_Also: `distribution.params.GIGEtaOmega`
+/
Tuple!(T, "lambda", T, "eta", T, "omega") 
generalizedInverseGaussianEstimate(T, OmegaBoundsFun)
	(T[2] lambdaBounds, scope OmegaBoundsFun omegaBoundsFun, in T sample[])
	if(isFloatingPoint!T)
in {
	import std.algorithm : all;
	assert(sample.all!"a > 0 && isNormal(a)");
}
body {

	import std.algorithm : sum, map;
	immutable n = sample.length;
	immutable one = sample.sum / n;
	immutable mone = sample.map!"1/a".sum / n;
	immutable dzero = T(LN2) * sample.sumOfLog2s() / n;
	return generalizedInverseGaussianEstimate!(T, OmegaBoundsFun)(lambdaBounds, omegaBoundsFun, one, mone, dzero);
}

/++
Estimates parameters of the generalized inverse Gaussian distribution.
Params:
	lambdaBounds = bounds of the parameter lambda
	omegaBoundsFun = function to calculate bounds of the parameter omega
	sample = sample
	weights = weights for the sample

See_Also: `distribution.params.GIGEtaOmega`
+/
Tuple!(T, "lambda", T, "eta", T, "omega") 
generalizedInverseGaussianEstimate(T, OmegaBoundsFun)
	(T[2] lambdaBounds, scope OmegaBoundsFun omegaBoundsFun, in T sample[], in T weights[])
	if(isFloatingPoint!T)
in {
	import std.algorithm : all, any;
	assert(weights.length == sample.length);
	assert(sample.all!"a > 0 && isNormal(a)");
	assert(weights.all!"a >= 0 && isFinite(a)");
	assert(weights.any!"a > 0");
}
body {

	import std.algorithm : sum, map;
	immutable n = weights.sum;
	immutable one = sample.dotProduct(weights) / n;
	immutable mone = sample.map!"1/a".dotProduct(weights) / n;
	immutable dzero = T(LN2) * sample.map!log2.dotProduct(weights) / n;
	return generalizedInverseGaussianEstimate!(T, OmegaBounds)(lambdaBounds, omegaBoundsFun, one, mone, dzero);
}

/++
Estimates parameters of the generalized inverse Gaussian distribution.
Params:
	lambdaBounds = bounds of the parameter lambda
	omegaBoundsFun = function to calculate bounds of the parameter omega
	one = `Σ weights[j] * sample[j] / Σ weights[j]`
	mone = `Σ weights[j] / sample[j] / Σ weights[j]`
	dzero = `Σ weights[j] * log(sample[j]) / Σ weights[j]`

See_Also: `distribution.params.GIGEtaOmega`
+/
Tuple!(T, "lambda", T, "eta", T, "omega") 
generalizedInverseGaussianEstimate(T, OmegaBoundsFun)
	(T[2] lambdaBounds, scope OmegaBoundsFun omegaBoundsFun, T one, T mone, T dzero)
	if(isFloatingPoint!T)
{
	import std.numeric : findRoot;
	auto f(T lambda)
	{
		immutable omegaBounds = omegaBoundsFun(lambda, one, mone, dzero);
		immutable params = generalizedInverseGaussianEstimate!T(lambda, omegaBounds, one, mone);
		immutable eta = params.eta;
		T dzeroScaled;
		if(eta <= T.min_normal*16)
		{
			immutable tau = lambda / params.omega;
			dzeroScaled = dzero - (log(one/(2*tau)) - one*mone/(4*tau^^2));
		}
		else
		{
			dzeroScaled = dzero - log(eta);
		}
		return GIGLikelihoodLambdaDerivative!T(lambda, params.omega, dzeroScaled);
	}
	immutable fax = f(lambdaBounds[0]);
	immutable fbx = f(lambdaBounds[1]);
	T lambda;
	if(signbit(fax) == signbit(fbx))
	{
		lambda = !(fabs(fax) > fabs(fbx)) ? lambdaBounds[0] : lambdaBounds[1];
	}
	else
	{
	    auto r = findRoot(&f, lambdaBounds[0], lambdaBounds[1], fax, fbx, (T a, T b) => false);
	    lambda = !(fabs(r[2]) > fabs(r[3])) ? r[0] : r[1];
	}
	immutable omegaBounds = omegaBoundsFun(lambda, one, mone, dzero);
	immutable params = generalizedInverseGaussianEstimate!T(lambda, omegaBounds, one, mone);
	return typeof(return)(lambda, params.eta, params.omega);
}

/++
Estimates parameters of the generalized inverse Gaussian distribution.
Params:
	lambda = parameter lambda
	omegaBounds = bounds of the parameter omega
	sample = sample

See_Also: `distribution.params.GIGEtaOmega`
+/
Tuple!(T, "lambda", T, "eta", T, "omega") 
generalizedInverseGaussianEstimate(T)
	(T lambda, T[2] omegaBounds, in T[] sample)
	if(isFloatingPoint!T)
in {
	import std.algorithm : all;
	assert(sample.all!"a > 0 && isNormal(a)");
}
body {
	import std.algorithm : sum, map;
	immutable n = sample.length;
	immutable one = sample.sum / n;
	immutable mone = sample.map!"1/a".sum / n;
	return generalizedInverseGaussianEstimate!T(lambda, omegaBounds, one, mone);
}

/++
Estimates parameters of the generalized inverse Gaussian distribution.
Params:
	lambda = parameter lambda
	omegaBounds = bounds of the parameter omega
	sample = sample
	weights = weights for the sample

See_Also: `distribution.params.GIGEtaOmega`
+/
Tuple!(T, "lambda", T, "eta", T, "omega") 
generalizedInverseGaussianEstimate(T)
	(T lambda, T[2] omegaBounds, in T[] sample, in T[] weights)
	if(isFloatingPoint!T)
in {
	import std.algorithm : all, any;
	assert(weights.length == sample.length);
	assert(sample.all!"a > 0 && isNormal(a)");
	assert(weights.all!"a >= 0 && isFinite(a)");
	assert(weights.any!"a > 0");
}
body {
	immutable n = weights.sum;
	immutable one = sample.dotProduct(weights) / n;
	immutable mone = sample.map!"1/a".dotProduct(weights) / n;
	return generalizedInverseGaussianEstimate!T(lambda, omegaBounds, one, mone);
}

/++
Estimates parameters of the generalized inverse Gaussian distribution.
Params:
	lambda = parameter lambda
	omegaBounds = bounds of the parameter omega
	one = `Σ weights[j] * sample[j] / Σ weights[j]`
	mone = `Σ weights[j] / sample[j] / Σ weights[j]`

See_Also: `distribution.params.GIGEtaOmega`
+/
Tuple!(T, "lambda", T, "eta", T, "omega") 
generalizedInverseGaussianEstimate(T)
	(T lambda, T[2] omegaBounds, T one, T mone)
	if(isFloatingPoint!T)
{
	import std.algorithm : sum, map;
	immutable prod = one * mone;
	immutable omega = GIGLikelihoodOmegaDerivativeZero!T(lambda, omegaBounds, prod);
	immutable tau = lambda / omega;
	immutable eta = (-tau + sqrt(tau^^2 + prod)) / mone;
	return typeof(return)(lambda, eta, omega);
}

private:

T GIGLikelihoodOmegaDerivative(T)(T lambda, T omega, T prod)
{
	import bessel;
	immutable alambda = fabs(lambda);
	immutable tau = alambda / omega;
	return 
		-besselK(omega, alambda-1, Flag!"ExponentiallyScaled".yes) / 
		besselK(omega, alambda, Flag!"ExponentiallyScaled".yes) 
		- (tau-sqrt(tau^^2 + prod));
}

T GIGLikelihoodLambdaDerivative(T)(T lambda, T omega, T dzeroScaled)
{
	return -logBesselKParamDerivative!T(omega, lambda) + dzeroScaled;
}

T GIGLikelihoodOmegaDerivativeZero(T)(T lambda, T[2] omegaBounds, T prod)
{
	import std.numeric : findRoot;
	T f(T omega)
	{
		return GIGLikelihoodOmegaDerivative!T(lambda, omega, prod);
	}
	immutable fax = f(omegaBounds[0]);
	immutable fbx = f(omegaBounds[1]);
	version(none)
	debug {
		import std.stdio;
		immutable step = (omegaBounds[1]-omegaBounds[0])/100;
		foreach(om; iota(omegaBounds[0], omegaBounds[1]+step/2, step)) {
			writeln(om, " ", f(om));
		}
	}
	if(signbit(fax) == signbit(fbx))
	{
		return !(fabs(fax) > fabs(fbx)) ? omegaBounds[0] : omegaBounds[1];
	}
    auto r = findRoot(&f, omegaBounds[0], omegaBounds[1], fax, fbx, (T a, T b) => false);
    return !(fabs(r[2]) > fabs(r[3])) ? r[0] : r[1];
}

T logBesselKParamDerivative(T)(T x, T alpha)
{
	import bessel;
	import std.typecons;
	import scid.calculus;
	immutable ret = diff((T alpha) => log(besselK(x, alpha, Flag!"ExponentiallyScaled".yes)), alpha, T(1), 16);
	return ret.value;
}

unittest {
	auto l = logBesselKParamDerivative!real(.60L, -2.0L);
	assert(isNormal(l));
}

version(none)
unittest {
	import atmosphere;
	import std.random;
	import std.algorithm;
	import std.range;
	import std.stdio;
	import std.conv;
	double lambda = -14;
	double eta = 1;
	double omega = 30;
	auto rng = ProperGeneralizedInverseGaussianSRNG!double(rndGen, lambda, eta, omega);
	auto sample = rng.take(100).array;
	auto params = generalizedInverseGaussianEstimate(lambda, [1e-10, 1e10], sample);
	writeln(params);
	writefln("\t%a", params.lambda);
	writefln("\t%a", params.eta);
	writefln("\t%a", params.omega);
//	writefln("\t%a", params.one);
//	writefln("\t%a", params.mone);
	auto lhPrior = generalizedInverseGaussianLikelihood(lambda, eta, omega, sample);
	auto lhCalc  = generalizedInverseGaussianLikelihood(params.lambda, params.eta, params.omega, sample);
	//assert(lhPrior <= lhCalc, text(lhPrior-lhCalc));
	writefln("ист = %s оцен = %s крит = %s", lhPrior, lhCalc, lhPrior <= lhCalc);
	auto params2 = generalizedInverseGaussianEstimate(sample);
	writeln(params2);
	writefln("\t%a", params2.lambda);
	writefln("\t%a", params2.eta);
	writefln("\t%a", params2.omega);
	//writefln("\t%a", params2.one);
	//writefln("\t%a", params2.mone);
	//writefln("\t%a", params2.dzero);

	auto lhCalc2  = generalizedInverseGaussianLikelihood(params2.lambda, params2.eta, params2.omega, sample);
	////assert(lhPrior <= lhCalc);
	writefln("ист = %s оцен = %s крит = %s", lhPrior, lhCalc2, lhPrior <= lhCalc2 && lhCalc <= lhCalc2);
}
