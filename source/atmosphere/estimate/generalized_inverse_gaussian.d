/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.estimate.generalized_inverse_gaussian;

import core.stdc.tgmath;

import std.traits;
import std.typecons;
import std.math : isNormal, isFinite;

import atmosphere.statistic: GeneralizedInverseGaussinStatistic, GeneralizedInverseGaussinFixedLambdaStatistic;


/++
Estimates parameters of the proper generalized inverse Gaussian distribution.

Params:
	stat = GIG statistica.
	sample = obeservation.
	weights = sample weights.
	relTolerance = Relative tolerance.
	absTolerance = Absolute tolerance.

Preconditions: $(D ax) and $(D bx) shall be finite reals. $(BR)
	$(D relTolerance) shall be normal positive real. $(BR)
	$(D absTolerance) shall be normal positive real no less then $(D T.epsilon*2).

References: "Algorithms for Minimization without Derivatives", Richard Brent, Prentice-Hall, Inc. (1973)

See_Also: `atmosphere.params`
+/
Tuple!(T, "lambda", T, "eta", T, "omega")
properGeneralizedInverseGaussianEstimate(T)
	(
		in T[] sample, 
		in T relTolerance = sqrt(T.epsilon),
		in T absTolerance = sqrt(T.epsilon),
	)
	if(isFloatingPoint!T)
{
	return properGeneralizedInverseGaussianEstimate(GeneralizedInverseGaussinStatistic!T(sample), relTolerance, absTolerance);
}

///
unittest
{
	import std.range;
	import std.random;
	import atmosphere.random;
	import atmosphere.likelihood;
	auto length = 1000;
	auto lambda = -2.0, eta = 1.4, omega = 2.3;
	auto rng = Random(1234);
	auto sample = ProperGeneralizedInverseGaussianSRNG!double(rng, lambda, eta, omega).take(length).array;
	auto params1 = properGeneralizedInverseGaussianEstimate(sample);
	auto params2 = properGeneralizedInverseGaussianFixedLambdaEstimate!double(lambda, sample);
	auto lh0 = properGeneralizedInverseGaussianLikelihood(lambda, eta, omega, sample);
	auto lh1 = properGeneralizedInverseGaussianLikelihood(params1.lambda, params1.eta, params1.omega, sample);
	auto lh2 = properGeneralizedInverseGaussianLikelihood(lambda, params2.eta, params2.omega, sample);
	assert(lh0 <= lh1);
	assert(lh0 <= lh2);
	assert(lh2 <= lh1);
}


///ditto
Tuple!(T, "lambda", T, "eta", T, "omega")
properGeneralizedInverseGaussianEstimate(T)
	(
		in T[] sample, 
		in T[] weights,
		in T relTolerance = sqrt(T.epsilon),
		in T absTolerance = sqrt(T.epsilon),
	)
	if(isFloatingPoint!T)
{
	return properGeneralizedInverseGaussianEstimate(GeneralizedInverseGaussinStatistic!T(sample, weights), relTolerance, absTolerance);
}


///ditto
Tuple!(T, "lambda", T, "eta", T, "omega")
properGeneralizedInverseGaussianEstimate(T)
	(
		in GeneralizedInverseGaussinStatistic!T stat, 		
		in T relTolerance = sqrt(T.epsilon),
		in T absTolerance = sqrt(T.epsilon),
	)
	if(isFloatingPoint!T)
body {
	import atmosphere.likelihood.generalized_inverse_gaussian;
	immutable prod = stat.mean * stat.meani;
	if(!isFinite(prod) || prod <= 1)
		return typeof(return)(T.nan, T.nan, T.nan);
	immutable u = prod / (prod - 1);
	assert(isFinite(u));
	immutable statFL = cast(GeneralizedInverseGaussinFixedLambdaStatistic!T) stat;
	T f(T lambda) {
		immutable params = properGeneralizedInverseGaussianFixedLambdaEstimate!T(lambda, statFL);
		with(params) return -properGeneralizedInverseGaussianLikelihood!T(lambda, eta, omega, stat);
	}
	import atmosphere.math: findLocalMin;
	immutable flm = findLocalMin!T(&f, -u, u, relTolerance, absTolerance);
	immutable lambda = flm.x;
	immutable params = properGeneralizedInverseGaussianFixedLambdaEstimate!T(lambda, statFL);
	with(params) return typeof(return)(lambda, eta, omega);
}

/++
Estimates parameters of the generalized inverse Gaussian distribution.

Params:
	stat = GIG statistica.
	sample = obeservation.
	weights = sample weights.
	relTolerance = Relative tolerance.
	absTolerance = Absolute tolerance.

Preconditions: $(D ax) and $(D bx) shall be finite reals. $(BR)
	$(D relTolerance) shall be normal positive real. $(BR)
	$(D absTolerance) shall be normal positive real no less then $(D T.epsilon*2).

References: "Algorithms for Minimization without Derivatives", Richard Brent, Prentice-Hall, Inc. (1973)

See_Also: `atmosphere.params`
+/
Tuple!(T, "lambda", T, "chi", T, "psi")
generalizedInverseGaussianEstimate(T)
	(
		in T[] sample,
		in T relTolerance = sqrt(T.epsilon),
		in T absTolerance = sqrt(T.epsilon),
	)
	if(isFloatingPoint!T)
{
	return generalizedInverseGaussianEstimate(GeneralizedInverseGaussinStatistic!T(sample), relTolerance, absTolerance);
}

///
unittest
{
	import std.range;
	import std.random;
	import atmosphere.random;
	import atmosphere.likelihood;
	import atmosphere.params;
	auto length = 1000;
	auto lambda = -2.0, chi = 1.4, psi = 2.3;
	auto rng = Random(1234);
	auto sample = new GeneralizedInverseGaussianRNG!double(rng, lambda, chi, psi).take(length).array;
	auto params1 = generalizedInverseGaussianEstimate(sample);
	auto params2 = generalizedInverseGaussianFixedLambdaEstimate!double(lambda, sample);
	auto p0 = GIGChiPsi!double(chi, psi);
	auto p1 = GIGChiPsi!double(params1.chi, params1.psi);
	auto p2 = GIGChiPsi!double(params2.chi, params2.psi);
	auto lh0 = properGeneralizedInverseGaussianLikelihood(lambda, p0.eta, p0.omega, sample);
	auto lh1 = properGeneralizedInverseGaussianLikelihood(params1.lambda, p1.eta, p1.omega, sample);
	auto lh2 = properGeneralizedInverseGaussianLikelihood(lambda, p2.eta, p2.omega, sample);
	assert(lh0 <= lh1);
	assert(lh0 <= lh2);
	assert(lh2 <= lh1);
}


///ditto
Tuple!(T, "lambda", T, "chi", T, "psi")
generalizedInverseGaussianEstimate(T)
	(
		in T[] sample,
		in T[] weights,
		in T relTolerance = sqrt(T.epsilon),
		in T absTolerance = sqrt(T.epsilon),
	)
	if(isFloatingPoint!T)
{
	return generalizedInverseGaussianEstimate(GeneralizedInverseGaussinStatistic!T(sample, weights), relTolerance, absTolerance);
}


///ditto
Tuple!(T, "lambda", T, "chi", T, "psi")
generalizedInverseGaussianEstimate(T)
	(
		in GeneralizedInverseGaussinStatistic!T stat, 		
		in T relTolerance = sqrt(T.epsilon),
		in T absTolerance = sqrt(T.epsilon),
	)
	if(isFloatingPoint!T)
body {
	import atmosphere.params : GIGEtaOmega;
	import atmosphere.estimate.gamma;
	import atmosphere.estimate.inverse_gamma;
	import atmosphere.likelihood.gamma;
	import atmosphere.likelihood.inverse_gamma;
	import atmosphere.likelihood.generalized_inverse_gaussian;
	import atmosphere.statistic: GammaStatistic, InverseGammaStatistic;
	
	immutable prod = stat.mean * stat.meani;
	if(!isFinite(prod) || prod <= 1)
		return typeof(return)(T.nan, T.nan, T.nan);
	immutable u = prod / (prod - 1);
	assert(isFinite(u));
	
	immutable gammaStat = cast(GammaStatistic!T)stat;
	immutable inverseGammaStat = cast(InverseGammaStatistic!T)stat;
	
	immutable gammaParams = gammaEstimate!T(gammaStat);
	immutable inverseGammaParams = inverseGammaEstimate!T(inverseGammaStat);
	immutable properParams = properGeneralizedInverseGaussianEstimate!T(stat, relTolerance, absTolerance);
	
	T gammaLikelihood = gammaLikelihood!T(gammaParams.shape, gammaParams.scale, gammaStat);
	T inverseGammaLikelihood = inverseGammaLikelihood!T(inverseGammaParams.shape, inverseGammaParams.scale, inverseGammaStat);
	T properLikelihood = properGeneralizedInverseGaussianLikelihood!T(properParams.lambda, properParams.eta, properParams.omega,  stat);
	
	if(!(gammaLikelihood > -T.infinity))
		gammaLikelihood = T.infinity;
	if(!(inverseGammaLikelihood > -T.infinity))
		inverseGammaLikelihood = T.infinity;
	if(!(properLikelihood > -T.infinity))
		properLikelihood = T.infinity;
	if(properLikelihood > gammaLikelihood && properLikelihood > inverseGammaLikelihood)
		with(GIGEtaOmega!T(properParams.eta, properParams.omega)) return typeof(return)(properParams.lambda, chi, psi);
	if(gammaLikelihood > inverseGammaLikelihood)
		with(       gammaParams) return typeof(return)(shape,         0, 2 / scale);
	if(gammaLikelihood < inverseGammaLikelihood)
		with(inverseGammaParams) return typeof(return)(shape, scale / 2,         0);
	return typeof(return).init;
}


/++
Estimates parameters of the proper generalized inverse Gaussian distribution for fixed `lambda`.

See_Also: `atmosphere.params`
+/
Tuple!(T, "eta", T, "omega")
properGeneralizedInverseGaussianFixedLambdaEstimate(T)(in T lambda, in T[] sample)
	if(isFloatingPoint!T)
{
	return properGeneralizedInverseGaussianFixedLambdaEstimate(lambda, GeneralizedInverseGaussinFixedLambdaStatistic!T(sample));
}

///
unittest
{
	import std.range;
	import std.random;
	import atmosphere.random;
	import atmosphere.likelihood;
	auto length = 1000;
	auto lambda = -2.0, eta = 1.4, omega = 2.3;
	auto rng = Random(1234);
	auto sample = ProperGeneralizedInverseGaussianSRNG!double(rng, lambda, eta, omega).take(length).array;
	auto params = properGeneralizedInverseGaussianFixedLambdaEstimate!double(lambda, sample);
	auto lh0 = properGeneralizedInverseGaussianLikelihood(lambda, eta, omega, sample);
	auto lh1 = properGeneralizedInverseGaussianLikelihood(lambda, params.eta, params.omega, sample);
	assert(lh0 <= lh1);
}

///ditto
Tuple!(T, "eta", T, "omega")
properGeneralizedInverseGaussianFixedLambdaEstimate(T)(in T lambda, in T[] sample, in T[] weights)
	if(isFloatingPoint!T)
{
	return properGeneralizedInverseGaussianFixedLambdaEstimate(lambda, GeneralizedInverseGaussinFixedLambdaStatistic!T(sample, weights));
}

///ditto
Tuple!(T, "eta", T, "omega")
properGeneralizedInverseGaussianFixedLambdaEstimate(T)
	(in T lambda, in GeneralizedInverseGaussinFixedLambdaStatistic!T stat)
	if(isFloatingPoint!T)
in {
	assert(isNormal(lambda));
}
body {
	import std.numeric: findRoot;
	import atmosphere.math;
	immutable prod = stat.mean * stat.meani;
	if(!isFinite(prod) || prod <= 1)
		return typeof(return)(T.nan, T.nan);
	immutable u = prod / (prod - 1);
	assert(isFinite(u));
	immutable alambda = fabs(lambda);
	if(alambda >= u)
		return typeof(return)(alambda > 0 ? 0 : T.infinity, 0);
	immutable omega = findRoot((T omega) => besselKD(alambda, omega) - prod, T.min_normal, T.max);
	immutable eta = sqrt(stat.mean / stat.meani) / besselKRM(lambda, omega);
	return typeof(return)(eta, omega);
}


/++
Estimates parameters of the generalized inverse Gaussian distribution for fixed `lambda`.

See_Also: `atmosphere.params`
+/
Tuple!(T, "chi", T, "psi")
generalizedInverseGaussianFixedLambdaEstimate(T)(in T lambda, in T[] sample)
	if(isFloatingPoint!T)
{
	return generalizedInverseGaussianFixedLambdaEstimate(lambda, GeneralizedInverseGaussinFixedLambdaStatistic!T(sample));
}

///
unittest
{
	import std.range;
	import std.random;
	import atmosphere.random;
	import atmosphere.likelihood;
	import atmosphere.params;
	auto length = 1000;
	auto lambda = -2.0, chi = 1.4, psi = 2.3;
	auto rng = Random(1234);
	auto sample = new GeneralizedInverseGaussianRNG!double(rng, lambda, chi, psi).take(length).array;
	auto params = generalizedInverseGaussianFixedLambdaEstimate!double(lambda, sample);
	auto p0 = GIGChiPsi!double(chi, psi);
	auto p1 = GIGChiPsi!double(params.chi, params.psi);
	auto lh0 = properGeneralizedInverseGaussianLikelihood(lambda, p0.chi, p0.psi, sample);
	auto lh1 = properGeneralizedInverseGaussianLikelihood(lambda, p1.chi, p1.psi, sample);
	assert(lh0 <= lh1);
}


///ditto
Tuple!(T, "chi", T, "psi")
generalizedInverseGaussianFixedLambdaEstimate(T)(in T lambda, in T[] sample, in T[] weights)
	if(isFloatingPoint!T)
{
	return generalizedInverseGaussianFixedLambdaEstimate(lambda, GeneralizedInverseGaussinFixedLambdaStatistic!T(sample, weights));
}


///ditto
Tuple!(T, "chi", T, "psi")
generalizedInverseGaussianFixedLambdaEstimate(T)
	(in T lambda, in GeneralizedInverseGaussinFixedLambdaStatistic!T stat)
	if(isFloatingPoint!T)
in {
	assert(isNormal(lambda));
}
body {
	import atmosphere.params : GIGEtaOmega;
	immutable prod = stat.mean * stat.meani;
	if(!isFinite(prod) || prod <= 1)
		return typeof(return)(T.nan, T.nan);
	immutable u = prod / (prod - 1);
	assert(isFinite(u));
	immutable alambda = fabs(lambda);
	if(alambda >= u)
		return lambda > 0 ? typeof(return)(0, 2 * alambda / stat.mean) : typeof(return)(2 * alambda / stat.meani, 0);
	immutable properParams = properGeneralizedInverseGaussianFixedLambdaEstimate!T(lambda, stat);
	with(GIGEtaOmega!T(properParams.eta, properParams.omega)) return typeof(return)(chi, psi);
}


/++
Estimates parameter `omega` of the proper generalized inverse Gaussian distribution for fixed `lambda` and `eta`.

See_Also: `atmosphere.params`
+/
T properGeneralizedInverseGaussianFixedLambdaEtaEstimate(T)(in T lambda, in T eta, in T[] sample)
	if(isFloatingPoint!T)
{
	return properGeneralizedInverseGaussianFixedLambdaEtaEstimate(lambda, eta, GeneralizedInverseGaussinFixedLambdaStatistic!T(sample));
}

///
unittest
{
	import std.range;
	import std.random;
	import atmosphere.random;
	import atmosphere.likelihood;
	auto length = 1000;
	auto lambda = -2.0, eta = 1.4, omega = 2.3;
	auto rng = Random(1234);
	auto sample = ProperGeneralizedInverseGaussianSRNG!double(rng, lambda, eta, omega).take(length).array;
	auto omega1 = properGeneralizedInverseGaussianFixedLambdaEtaEstimate!double(lambda, eta, sample);
	auto lh0 = properGeneralizedInverseGaussianLikelihood(lambda, eta, omega , sample);
	auto lh1 = properGeneralizedInverseGaussianLikelihood(lambda, eta, omega1, sample);
	assert(lh0 <= lh1);
}

///ditto
T properGeneralizedInverseGaussianFixedLambdaEtaEstimate(T)(in T lambda, in T eta, in T[] sample, in T[] weights)
	if(isFloatingPoint!T)
{
	return properGeneralizedInverseGaussianFixedLambdaEtaEstimate(lambda, eta, GeneralizedInverseGaussinFixedLambdaStatistic!T(sample, weights));
}

///ditto
T properGeneralizedInverseGaussianFixedLambdaEtaEstimate(T)
	(in T lambda, in T eta, in GeneralizedInverseGaussinFixedLambdaStatistic!T stat)
	if(isFloatingPoint!T)
in {
	assert(isNormal(lambda));
}
body {
	import std.numeric: findRoot;
	import atmosphere.math;
	immutable d = stat.mean / eta + stat.meani * eta;
	return findRoot((T omega) => besselKRS(lambda, omega) - d, double.min_normal, double.max);
}


version(none)
unittest
{
	import std.range;
	import std.random;
	import atmosphere.random;
	import atmosphere.likelihood;
	import atmosphere.params;
	auto length = 1000;
	auto lambda = -2.0, chi = 1.4, psi = 2.3;
	auto rng = Random(1234);
	auto sample = new GeneralizedInverseGaussianRNG!double(rng, lambda, chi, psi).take(length).array;
	auto chi1 = generalizedInverseGaussianFixedLambdaPsiEstimate!double(lambda, psi, sample);
	auto p0 = GIGChiPsi!double(chi , psi);
	auto p1 = GIGChiPsi!double(chi1, psi);
	auto lh0 = properGeneralizedInverseGaussianLikelihood(lambda, p0.eta, p0.omega , sample);
	auto lh1 = properGeneralizedInverseGaussianLikelihood(lambda, p1.eta, p1.omega, sample);
	assert(lh0 <= lh1);
}
