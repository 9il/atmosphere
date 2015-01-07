/++
Derivatives of probability density functions
+/
module distribution.derivative;

import bessel;
import std.traits;
import std.numeric;
import std.mathspecial;
import std.typecons;


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
