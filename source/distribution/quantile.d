/++
Quantile functions
+/
module distribution.quantile;

import std.traits;
import std.mathspecial;

/++
Quantile function interface
+/
interface Quantile(T)
{
	/++
	Call operator
	+/
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
	assert(isNormal(x));
}


/++
Class to compute quantile function as root of it's cumulative density function
+/
abstract class NumericQuantile(T) : Quantile!T
{
	import distribution.cdf;

	private CDF!T cdf;
	private T a, b;
	private scope bool delegate(T lo, T hi) tolerance;

	/++
    Params:
		cdf     = The CDF to _integrate.
		a	= (optional) The lower limit of integration.
		b	= (optional) The upper limit of integration.
	+/
	this(CDF!T cdf, T a = -T.max, T b = T.max)
	{
		this.cdf = cdf;
		this.a = a;
		this.b = b;
		this.tolerance = tolerance;
	}

	/++
	Call operator
	+/
	final T opCall(T x)
	in {
		assert(x >= 0);
		assert(x <= 1);
	}
	out(result) {
		assert(!result.isNaN);
	}
	body {
		import std.numeric : findRoot;
		T f(T y)
			in     { assert(y.isFinite); }
			out(r) { assert(!r.isNaN);   }
			body   { return cdf(y) - x;  }
		//return tolerance ? findRoot(&f, a, b) : findRoot(&f, a, b, tolerance);
		return findRoot(&f, a, b);
	}
}

///
unittest
{
	import std.traits, std.mathspecial;
	import distribution.pdf;
	import distribution.cdf;

	class NormalPDF : PDF!double
	{
		double opCall(double x)
		{
			// 1/sqrt(2 PI)
			enum c = 0.398942280401432677939946L;
			return c * exp(-0.5f * x * x);
		}
	}

	class NormalCDF : NumericCDF!double
	{
		this()
		{
			super(new NormalPDF);
		}
	}

	class NormalQuantile : NumericQuantile!double
	{
		this()
		{
			super(new NormalCDF, -20, 20);
		}
	}
	import scid.calculus;
	auto l = integrate(new NormalPDF, -double.infinity, 20);

	auto qf = new NormalQuantile;

	assert(approxEqual(qf(0.3), normalDistributionInverse(0.3)));
}


/++
Quantile function of the gamma distribution
+/
final class GammaQuantile(T) : Quantile!T
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
		assert(x >= 0);
		assert(x <= 1);
		return scale * gammaIncompleteComplInverse(shape, 1-x);
	}
}

///
unittest 
{
	auto qf = new GammaQuantile!double(3, 2);
	auto x = qf(0.1);
	assert(isNormal(x));
}


/++
Quantile function of the inverse-gamma distribution
+/
final class InverseGammaQuantile(T) : Quantile!T
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
		assert(x >= 0);
		assert(x <= 1);
		return scale / gammaIncompleteComplInverse(shape, 1-x);
	}
}

///
unittest 
{
	auto qf = new InverseGammaQuantile!double(3, 2);
	auto x = qf(0.1);
	assert(isNormal(x));
}
