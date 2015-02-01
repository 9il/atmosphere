module distribution.finitemixture;

import std.range;
import std.traits;

/++
+/
template finiteMixture(alias InterfaceTemp)
{
	class FiniteMixture(T, FRange) : InterfaceTemp!T
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
			auto r = funcs.save;
			foreach(w; weights)
			{
				assert(!r.empty);
				static if(isSomeFunction!(ElementType!FRange))
				{
					auto e = r.front;
					s += w * e(x);
				}
				else
					s += r.front.opCall(x);
				r.popFront;
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
		return new FiniteMixture!(T, FRange)(funcs, weights);
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
