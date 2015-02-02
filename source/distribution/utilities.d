module distribution.utilities;

import std.traits;

///
template convertTo(alias InterfaceTemp)
{
	InterfaceTemp!F convertTo(Fun, F = ReturnType!Fun)(Fun fun)
	{
		static assert(isFloatingPoint!F);
		return new class (fun) InterfaceTemp!F {
			Fun fun;
			this(Fun fun) { this.fun = fun; }
			F opCall(F x) { return fun(x); }
		};
	}
}

///
unittest
{
	import std.math;
	import distribution.pdf;

	real fun(real x)
	{
		// 1/sqrt(2 PI)
		enum c = 0.398942280401432677939946L;
		return c * exp(-0.5f * x * x);
	}

	PDF!real pdf = convertTo!PDF(&fun);
}
