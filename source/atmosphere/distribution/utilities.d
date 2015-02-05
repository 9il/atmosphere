module atmosphere.distribution.utilities;

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
	import atmosphere.distribution.pdf;

	real fun(real x)
	{
		// 1/sqrt(2 PI)
		enum c = 0.398942280401432677939946L;
		return c * exp(-0.5f * x * x);
	}

	PDF!real pdf = convertTo!PDF(&fun);
}

/++
Generates random permutation
+/
size_t[] randomPermutation(size_t length)
{
	import core.memory;
	import std.random : uniform;
	import std.algorithm : makeIndex;
	auto indexesR = new size_t[length];
	scope(exit) 
		GC.free(indexesR.ptr);
	auto indexesS = new size_t[length];
	foreach(j, ref index; indexesR)
	{
		index = uniform!"[]"(0, size_t.max);
		indexesS[j] = j;
	}
	makeIndex(indexesR, indexesS);
	return indexesS;
}
