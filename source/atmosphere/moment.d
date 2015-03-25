/++
Moment of probability distributions.
+/
/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.moment;

import atmosphere.math: besselKR;

import core.stdc.tgmath;

import std.typecons;


/++
Moments of the generalized gamma distribution
+/
T generalizedGammaMean(T)(T shape, T scale, T power, uint order = 1)
{
	import std.math : pow;
	return scale.pow(order) * (tgamma(shape+order/power) / tgamma(shape));
}

///
unittest
{
	auto expectation     = generalizedGammaMean!double(2.0, 3.0, 4.0);
	auto secondRowMoment = generalizedGammaMean!double(2.0, 3.0, 4.0, 2);
}

/++
Moments of the generalized variance-gamma distribution
+/
T generalizedVarianceGammaMean(T)(T shape, T scale, T power, T beta)
{
	return beta * generalizedGammaMean(shape, scale, power);
}

///
unittest
{
	auto expectation = generalizedVarianceGammaMean!double(2.0, 3.0, 4.0, 5.0);
}


/++
Moments of the generalized inverse Gaussian distribution
+/
T properGeneralizedInverseGaussianMean(T)(T lambda, T eta, T omega, uint order = 1)
in {
	assert(order);
}
body {
	import std.math : pow;
	return eta.pow(order) * besselKR(lambda, omega, order-1);
}

///
unittest
{
	auto expectation     = properGeneralizedInverseGaussianMean!double(2.0, 3.0, 4.0);
	auto secondRowMoment = properGeneralizedInverseGaussianMean!double(2.0, 3.0, 4.0, 2);
}


/++
Variance of the generalized inverse Gaussian distribution
+/
//(besselK[x-1, y] * besselK[x+1, y] - besselK[x, y]^2) / besselK[x-1, y]^2
T properGeneralizedInverseGaussianVariance(T)(T lambda, T eta, T omega)
{
	immutable r = besselKR(lambda, omega);
	return eta^^2 * (r * (2 * (lambda + 1) / omega - r) + 1);
}

///
unittest
{
	import std.math : approxEqual;
	auto expectation     = properGeneralizedInverseGaussianMean!double(2.0, 3.0, 4.0);
	auto secondRowMoment = properGeneralizedInverseGaussianMean!double(2.0, 3.0, 4.0, 2);
	auto variance        = properGeneralizedInverseGaussianVariance!double(2.0, 3.0, 4.0);
	assert(approxEqual(secondRowMoment - expectation^^2, variance));
}


/++
Mean of the generalized Hyperbolic distribution
+/
T generalizedHyperbolicMean(T)(T lambda, T beta, T eta, T omega)
{
	return beta * properGeneralizedInverseGaussianMean(lambda, eta, omega);
}

///
unittest
{
	auto expectation = generalizedHyperbolicMean!double(2.0, 5.0, 3.0, 4.0);
}


/++
Variance of the generalized Hyperbolic distribution
+/
T generalizedHyperbolicVariance(T)(T lambda, T beta, T eta, T omega)
{
	immutable r = besselKR(lambda, omega);
	return eta * r + (beta * eta) ^^ 2 * (r * (2 * (lambda + 1) / omega - r) + 1);
}

///
unittest
{
	auto variance = generalizedHyperbolicVariance!double(2.0, 5.0, 3.0, 4.0);
}
