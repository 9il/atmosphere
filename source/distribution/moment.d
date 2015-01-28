/++
Moment of probability distributions.
+/
module distribution.moment;

import bessel;
import std.typecons;
import std.mathspecial;


/++
Moments of the Generalized gamma distribution
+/
T GeneralizedGammaMean(T)(T shape, T scale, T power, uint order = 1)
{
	return scale^^order * (gamma(shape+order/power) / gamma(shape));
}

///
unittest
{
	auto expectation     = GeneralizedGammaMean!double(2.0, 3.0, 4.0);
	auto secondRowMoment = GeneralizedGammaMean!double(2.0, 3.0, 4.0, 2);
}

/++
Moments of the Generalized variance-gamma distribution
+/
T GeneralizedVarianceGammaMean(T)(T shape, T scale, T power, T beta)
{
	return beta * GeneralizedGammaMean(shape, scale, power);
}

///
unittest
{
	auto expectation = GeneralizedVarianceGammaMean!double(2.0, 3.0, 4.0, 5.0);
}


/++
Moments of the Generalized inverse Gaussian distribution
+/
T ProperGeneralizedInverseGaussianMean(T)(T lambda, T eta, T omega, uint order = 1)
{
	immutable k0 = besselK(omega, lambda        , Flag!"ExponentiallyScaled".yes);
	immutable kr = besselK(omega, lambda + order, Flag!"ExponentiallyScaled".yes);
	return eta.pow(order) * (kr / k0);
}

///
unittest
{
	auto expectation     = ProperGeneralizedInverseGaussianMean!double(2.0, 3.0, 4.0);
	auto secondRowMoment = ProperGeneralizedInverseGaussianMean!double(2.0, 3.0, 4.0, 2);
}


/++
Variance of the Generalized inverse Gaussian distribution
+/
T ProperGeneralizedInverseGaussianVariance(T)(T lambda, T eta, T omega)
{
	immutable k0 = besselK(omega, lambda    , Flag!"ExponentiallyScaled".yes);
	immutable k1 = besselK(omega, lambda + 1, Flag!"ExponentiallyScaled".yes);
	immutable k2 = besselK(omega, lambda + 2, Flag!"ExponentiallyScaled".yes);
	return eta^^2 * ((k0*k2 - k1^^2) / k0^^2);
}

///
unittest
{
	auto expectation     = ProperGeneralizedInverseGaussianMean!double(2.0, 3.0, 4.0);
	auto secondRowMoment = ProperGeneralizedInverseGaussianMean!double(2.0, 3.0, 4.0, 2);
	auto variance        = ProperGeneralizedInverseGaussianVariance!double(2.0, 3.0, 4.0);
	assert(approxEqual(secondRowMoment - expectation^^2, variance));
}


/++
Mean of the Generalized Hyperbolic distribution
+/
T GeneralizedHyperbolicMean(T)(T lambda, T beta, T eta, T omega)
{
	return beta * ProperGeneralizedInverseGaussianMean(lambda, eta, omega);
}

///
unittest
{
	auto expectation = GeneralizedHyperbolicMean!double(2.0, 5.0, 3.0, 4.0);
}


/++
Variance of the Generalized Hyperbolic distribution
+/
T GeneralizedHyperbolicVariance(T)(T lambda, T beta, T eta, T omega)
{
	immutable k0 = besselK(omega, lambda    , Flag!"ExponentiallyScaled".yes);
	immutable k1 = besselK(omega, lambda + 1, Flag!"ExponentiallyScaled".yes);
	immutable k2 = besselK(omega, lambda + 2, Flag!"ExponentiallyScaled".yes);
	return eta * (k1 / k0) + (beta * eta)^^2 * ((k0*k2 - k1^^2) / k0^^2);
}

///
unittest
{
	auto variance = GeneralizedHyperbolicVariance!double(2.0, 5.0, 3.0, 4.0);
}
