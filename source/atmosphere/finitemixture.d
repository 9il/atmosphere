/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.finitemixture;

import std.range;
import std.traits;

/++
Returns: Callable structure
Params:
	funcs = input range of functions
	weights = mixture weights
+/
auto finiteMixture(FRange, T)(FRange funcs, const(T)[] weights)
	if(
		isForwardRange!FRange && 
		isCallable!(ElementType!FRange) &&
		isFloatingPoint!T)
in {
	assert(weights.length);
	static if(hasLength!FRange)
		assert(funcs.length == weights.length);
}
body {
	struct Result(T, FRange)
		if(
			isForwardRange!FRange && 
			isCallable!(ElementType!FRange) &&
			isFloatingPoint!T)
	{
		private FRange funcs;
		private const(T)[] weights;
		this(T)(FRange funcs, const(T)[] weights)
		{
			this.funcs = funcs;
			this.weights = weights;
		}
		T opCall(T x)
		{
			T s = 0;
			auto funcs = this.funcs.save;
			foreach(w; weights)
			{
				assert(!funcs.empty);
				s += funcs.front()(x);
				funcs.popFront;
			}
			return s;
		}
	}
	return Result!(T, FRange)(funcs, weights);
}

///
unittest
{
	import atmosphere.pdf;
	import atmosphere.utilities;
	auto pdfs = sequence!"n+1"().map!(shape => GammaSPDF!real(shape, 1));
	double[] weights = [0.25, 0.5, 0.125, 0.125];
	auto pdf0 = finiteMixture(pdfs, weights);           //struct
	auto pdf1 = finiteMixture(pdfs[3..4].chain(pdfs[0..3]), [0.125, 0.25, 0.5, 0.125]); //struct
	PDF!double pdf2 = toPDF(pdf0);                      //object
	foreach(x; iota(1, 4))
		assert(approxEqual(pdf1(x), pdf2(x)));
}
