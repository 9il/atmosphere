///
module distribution.moment;

import bessel;
import std.typecons;
import std.math;


/++
Raw moment
+/
T E_GIG(T)(T lambda, T eta, T omega, uint order = 1)
{
	immutable k0 = besselK(omega, lambda        , Flag!"ExponentiallyScaled".yes);
	immutable kr = besselK(omega, lambda + order, Flag!"ExponentiallyScaled".yes);
	return eta.pow(order) * (kr / k0);
}

///
unittest
{
	auto expectation     = E_GIG!double(2.0, 3.0, 4.0);
	auto secondRowMoment = E_GIG!double(2.0, 3.0, 4.0, 2);
}


/++
Variance
+/
T V_GIG(T)(T lambda, T eta, T omega)
{
	immutable k0 = besselK(omega, lambda    , Flag!"ExponentiallyScaled".yes);
	immutable k1 = besselK(omega, lambda + 1, Flag!"ExponentiallyScaled".yes);
	immutable k2 = besselK(omega, lambda + 2, Flag!"ExponentiallyScaled".yes);
	return eta^^2 * ((k0*k2 - k1^^2) / k0^^2);
}

///
unittest
{
	auto expectation     = E_GIG!double(2.0, 3.0, 4.0);
	auto secondRowMoment = E_GIG!double(2.0, 3.0, 4.0, 2);
	auto variance        = V_GIG!double(2.0, 3.0, 4.0);
	assert(approxEqual(secondRowMoment - expectation^^2, variance));
}


/++
Expectation
+/
T E_GHyp(T)(T lambda, T beta, T eta, T omega)
{
	return beta * E_GIG(lambda, eta, omega);
}

///
unittest
{
	auto expectation = E_GHyp!double(2.0, 5.0, 3.0, 4.0);
}


/++
Variance
+/
T V_GHyp(T)(T lambda, T beta, T eta, T omega)
{
	immutable k0 = besselK(omega, lambda    , Flag!"ExponentiallyScaled".yes);
	immutable k1 = besselK(omega, lambda + 1, Flag!"ExponentiallyScaled".yes);
	immutable k2 = besselK(omega, lambda + 2, Flag!"ExponentiallyScaled".yes);
	return eta * (k1 / k0) + (beta * eta)^^2 * ((k0*k2 - k1^^2) / k0^^2);
}

///
unittest
{
	auto variance = V_GHyp!double(2.0, 5.0, 3.0, 4.0);
}
