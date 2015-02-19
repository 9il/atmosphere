/++
Moment of probability distributions.
+/
/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.distribution.moment;

import bessel;
import std.typecons;
import std.mathspecial;


/++
Moments of the generalized gamma distribution
+/
T generalizedGammaMean(T)(T shape, T scale, T power, uint order = 1)
{
	return scale.pow(order) * (gamma(shape+order/power) / gamma(shape));
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
{
	immutable k0 = besselK(omega, lambda        , Flag!"ExponentiallyScaled".yes);
	immutable kr = besselK(omega, lambda + order, Flag!"ExponentiallyScaled".yes);
	return eta.pow(order) * (kr / k0);
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
T properGeneralizedInverseGaussianVariance(T)(T lambda, T eta, T omega)
{
	immutable k0 = besselK(omega, lambda    , Flag!"ExponentiallyScaled".yes);
	immutable k1 = besselK(omega, lambda + 1, Flag!"ExponentiallyScaled".yes);
	immutable k2 = besselK(omega, lambda + 2, Flag!"ExponentiallyScaled".yes);
	return eta^^2 * ((k0*k2 - k1^^2) / k0^^2);
}

///
unittest
{
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
	immutable k0 = besselK(omega, lambda    , Flag!"ExponentiallyScaled".yes);
	immutable k1 = besselK(omega, lambda + 1, Flag!"ExponentiallyScaled".yes);
	immutable k2 = besselK(omega, lambda + 2, Flag!"ExponentiallyScaled".yes);
	return eta * (k1 / k0) + (beta * eta)^^2 * ((k0*k2 - k1^^2) / k0^^2);
}

///
unittest
{
	auto variance = generalizedHyperbolicVariance!double(2.0, 5.0, 3.0, 4.0);
}
