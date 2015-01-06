/++
+/
module distribution.cdf;


import std.traits;
import std.mathspecial;

/++
+/
interface CDF(T)
{
	/++
	Call operator
	+/
	T opCall(T x);
}


///
unittest 
{
	import std.traits, std.mathspecial;

	class NormalCDF : CDF!real
	{
		real opCall(real x)
		{
			return normalDistribution(x);
		}
	}

	auto cdf = new NormalCDF;
	auto x = cdf(0.1);
	assert(isNormal(x));
}

/++
+/
final class GammaCDF(T) : CDF!T
	if(isFloatingPoint!T)
{
	private T shape, scale;

	///Constructor
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
		return x < 0 ? 0 : gammaIncomplete(shape, x / scale);
	}
}

/++
+/
unittest 
{
	auto cdf = new GammaCDF!double(3, 2);
	auto x = cdf(0.1);
	assert(isNormal(x));
}


/++
+/
final class InverseGammaCDF(T) : CDF!T
	if(isFloatingPoint!T)
{
	private T shape, scale;

	///Constructor
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
		return x < 0 ? 0 : gammaIncomplete(shape, scale / x);
	}
}

/++
+/
unittest 
{
	auto cdf = new InverseGammaCDF!double(3, 2);
	auto x = cdf(0.1);
	assert(isNormal(x));
}
