module distribution.finitemixture;

import std.range;
import std.traits;

/++
+/
template finiteMixture(alias InterfaceTemp)
{
	class Result(T, FRange) : InterfaceTemp!T
		if(
			isForwardRange!FRange && 
			isCallable!(ElementType!FRange) &&
			isFloatingPoint!T)
	{
		private FRange funcs;
		private const(T)[] weights;

		///
		this(T)(FRange funcs, const(T)[] weights)
		{
			this.funcs = funcs;
			this.weights = weights;
		}

		///
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

	///
	auto finiteMixture(FRange, T)(FRange funcs, const(T)[] weights)
		if(
			isForwardRange!FRange && 
			isCallable!(ElementType!FRange) &&
			isFloatingPoint!T)
	{
		return new Result!(T, FRange)(funcs, weights);
	}
}



///
unittest
{
	import distribution.pdf;
	auto pdfs = sequence!"n+1"().map!(shape => new GammaPDF!real(shape, 1));
	double[] weights = [0.25, 0.5, 0.125, 0.125];
	PDF!double pdf = finiteMixture!PDF(pdfs, weights);
}
