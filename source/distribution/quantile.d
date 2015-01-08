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
	Constructor
    Params:
		cdf	= The CDF to to inverse.
		a	= (optional) The lower limit.
		b	= (optional) The upper limit.
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

/// Numeric quantile function of standard normal distribution
unittest
{
	import std.traits, std.mathspecial;
	import distribution.pdf;
	import distribution.cdf;

	class NormalPDF : PDF!real
	{
		real opCall(real x)
		{
			// 1/sqrt(2 PI)
			enum c = 0.398942280401432677939946L;
			return c * exp(-0.5f * x * x);
		}
	}

	class NormalCDF : NumericCDF!real
	{
		this()
		{
			super(new NormalPDF, [-3, -1, 0, 1, 3]);
		}
	}

	class NormalQuantile : NumericQuantile!real
	{
		this()
		{
			super(new NormalCDF);
		}
	}

	auto qf = new NormalQuantile;

	assert(approxEqual(qf(0.3), normalDistributionInverse(0.3)));
}

/// Numeric quantile function of Generalized Hyperbolic distribution
unittest
{
	import distribution.pdf;
	import distribution.cdf;
	import distribution.params;
	import distribution.moment;

	class GHypCDF: NumericCDF!real
	{
		this(real lambda, GHypChiPsi!real params)
		{
			immutable mu = 0;
			auto pdf = new GeneralizedHyperbolicPDF!real(lambda, params.alpha, params.beta, params.delta, mu);
			immutable expectation = E_GHyp!real(lambda, params.beta, params.chi, params.psi);
			super(pdf, [expectation]);	
		}
	}

	class GHypQuantile : NumericQuantile!real
	{
		this(real lambda, GHypChiPsi!real params)
		{
			super(new GHypCDF(lambda, params), -1000, 1000);	
		}
	}

	auto qf = new GHypQuantile(0.5, GHypChiPsi!real(5, 0.7, 0.6));
	import std.stdio, std.conv;
	assert(approxEqual(qf(0.95), 40.9263));
	assert(approxEqual(qf(0.99), 64.977));
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
