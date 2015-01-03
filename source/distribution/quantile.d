/**
*/
module distribution.quantile;

import std.traits;
import std.mathspecial;

/**
*/
interface Quantile(T)
{
	/**
	Call operator
	*/
	abstract T opCall(T x);
}

///
unittest
{
	import std.traits, std.mathspecial;

	class NormalQuantile : Quantile!real
	{
		real opCall(real x)
		in {
			assert(x >= 0);
			assert(x <= 1);
		}
		body {
			return normalDistributionInverse(x);
		}
	}

	auto qf = new NormalQuantile;
	auto x = qf(0.1);
}



/**
*/
final class GammaQuantile(T) : Quantile!T
	if(isFloatingPoint!T)
{
	private T shape, scale;

	///
	this(T shape, T scale)
	in {
		assert(shape.isNormal);
		assert(shape > 0);
		assert(scale.isNormal);
		assert(scale > 0);
	}
	body {
		this.shape = shape;
		this.scale = scale;
	}

	T opCall(T x)
	{
		assert(x >= 0);
		assert(x <= 1);
		return scale * gammaIncompleteComplInverse(shape, 1-x);
	}
}

/**
*/
unittest 
{
	auto qf = new GammaQuantile!double(3, 2);
	auto x = qf(0.1);
}


/**
*/
final class InverseGammaQuantile(T) : Quantile!T
	if(isFloatingPoint!T)
{
	private T shape, scale;

	///
	this(T shape, T scale)
	in {
		assert(shape.isNormal);
		assert(shape > 0);
		assert(scale.isNormal);
		assert(scale > 0);
	}
	body {
		this.shape = shape;
		this.scale = scale;
	}

	T opCall(T x)
	{
		assert(x >= 0);
		assert(x <= 1);
		return scale / gammaIncompleteComplInverse(shape, 1-x);
	}
}

/**
*/
unittest 
{
	auto qf = new InverseGammaQuantile!double(3, 2);
	auto x = qf(0.1);
}
