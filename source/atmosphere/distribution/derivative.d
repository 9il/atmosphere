/++
Derivatives of probability density functions
+/
module atmosphere.distribution.derivative;

import bessel;
import std.traits;
import std.numeric;
import std.mathspecial;
import std.typecons;

import atmosphere.distribution.utilities;

/++
Derivative of probability density function interface
+/
interface Derivative(T)
{
	/**
	Call operator
	*/
	T opCall(T x);
}

///
alias toDerivative = convertTo!Derivative;
