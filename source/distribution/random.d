/++
Random number generators
+/
module distribution.random;

import std.mathspecial;
import std.random;
import std.traits;
import distribution.params;

/++
Interface for infinity input range of random numbers.
+/
interface DistributionRNG(T)
{
	///always false
	enum empty = false;
	///do nothing
	static void popFront() @safe pure nothrow @nogc {}
	/++
	Returns: new random number.
	+/
	T front() @property;
}

///
unittest
{
	class NormalRNG : DistributionRNG!double
	{
		double front() @property
		{
			return rNormal();
		}
	}

	import std.range;
	auto rng = new NormalRNG;
	auto sigma = 4.0;
	auto mu = 2.0;
	auto sample = rng.map!(x => x * sigma + mu).take(9).array;
}


/++
Class to create normal variance-mean mixture random number generators.
Assume `U` has mixing probability density, `Y ~ N(0, 1)`.
Class constructs `RNG` for `Z = Y*U^(1/2)+beta*U`.
+/
class NormalVarianceMeanMixtureRNG(T, UniformRNG = Random) : DistributionRNG!T
	if (isFloatingPoint!T)
{
	private UniformRNG* rng;
	
	private DistributionRNG!T components;
	private T beta;

	/++
	Constructor
	Params:
		rng = uniform random number generator (`Y`)
		components = mixing random generator (`U`)
		beta = mixture scale parameter: `Y*U^(1/2)+beta*U`
	+/
	this(ref UniformRNG rng, DistributionRNG!T components, T beta)
	in {
		assert(beta.isFinite);
	}
	body {
		this.rng = &rng;
		this.components = components;
		this.beta = beta;
	}

	final T front() @property
	{
		immutable y = rng.rNormal!T;
		immutable u = components.front;
		return y * u.sqrt + beta * u;
	}
}

///
unittest
{
	import std.random;
	class MyVarianceGammaRNG : NormalVarianceMeanMixtureRNG!double
	{
		this(double beta, double shape, double scale)
		{
			auto components = new GammaRNG!double(rndGen, shape, scale);
			super(rndGen, components, beta);
		}
	}
}


/++
Class to generate random observations from a
gamma
distribution.
 +/
final class GammaRNG(T, UniformRNG = Random) : DistributionRNG!T
	if (isFloatingPoint!T)
{
	private UniformRNG* rng;
	
	private T shape, scale;

	/++
	Constructor
	Params:
		rng = uniform random number generator
		shape = shape parameter
		scale = scale parameter
	+/
	this(ref UniformRNG rng, T shape, T scale = 1)
	{
		this.rng = &rng;
		this.shape = shape;
		this.scale = scale;
	}

	T front() @property
	{
		return rng.rGamma(shape) * scale;
	}
}

///
unittest
{
	import std.range;
	auto rng = new GammaRNG!double(rndGen, 1.1, 1.1);
	auto sample = rng.map!(x => x + 4).take(9).array;
}


/++
Class to generate random observations from a
inverse-gamma
distribution.
 +/
final class InverseGammaRNG(T, UniformRNG = Random) : DistributionRNG!T
	if (isFloatingPoint!T)
{
	private UniformRNG* rng;
	
	private T shape, scale;

	/++
	Constructor
	Params:
		rng = uniform random number generator
		shape = shape parameter
		scale = scale parameter
	+/
	this(ref UniformRNG rng, T shape, T scale = 1)
	{
		this.rng = &rng;
		this.shape = shape;
		this.scale = scale;
	}

	T front() @property
	{
		return rng.rInverseGamma(shape) * scale;
	}
}

///
unittest
{
	import std.range;
	auto rng = new InverseGammaRNG!double(rndGen, 1.1, 1.1);
	auto sample = rng.map!(x => x + 4).take(9).array;
}


/++
Class to generate random observations from a
generalized gamma
distribution.
 +/
final class GeneralizedGammaRNG(T, UniformRNG = Random) : DistributionRNG!T
	if (isFloatingPoint!T)
{
	private UniformRNG* rng;
	
	private T shape, power, scale;

	/++
	Constructor
	Params:
		rng = uniform random number generator
		shape = shape parameter
		power = power parameter
		scale = scale parameter
	+/
	this(ref UniformRNG rng, T shape, T power, T scale = 1)
	{
		this.rng = &rng;
		this.shape = shape;
		this.power = power;
		this.scale = scale;
	}

	T front() @property
	{
		return rng.rGeneralizedGamma(shape, power) * scale;
	}
}

///
unittest
{
	import std.range;
	auto rng = new GeneralizedGammaRNG!double(rndGen, 1.1, 1.1, 1.1);
	auto sample = rng.map!(x => x + 4).take(9).array;
}


/++
Class to generate random observations from a
inverse Gaussian
distribution.
 +/
final class InverseGaussianRNG(T, UniformRNG = Random) : DistributionRNG!T
	if (isFloatingPoint!T)
{
	private UniformRNG* rng;
	
	private T mu, lambda;

	/++
	Constructor
	Params:
		rng = uniform random number generator
		mu = mu parameter
		lambda = lambda parameter
	+/
	this(ref UniformRNG rng, T mu, T lambda)
	{
		this.rng = &rng;
		this.mu = mu;
		this.lambda = lambda;
	}

	T front() @property
	{
		return rng.rInverseGaussian(mu, lambda);
	}
}

///
unittest
{
	import std.range;
	auto rng = new InverseGaussianRNG!double(rndGen, 1.1, 1.1);
	auto sample = rng.map!(x => x + 4).take(9).array;
}


/++
Class to generate random observations from a
proper (chi > 0, psi > 0) generalized inverse Gaussian distribution. 
The algorithm is based on that given by Dagpunar (1989).

References: $(LINK2 https://www.stat.auckland.ac.nz/~dscott/, Original R source code).
+/
final class ProperGeneralizedInverseGaussianRNG(T, UniformRNG = Random) : DistributionRNG!T
	if (isFloatingPoint!T)
{
	private UniformRNG* rng;
	
	private T lambdam1, eta, omega, a, b, c, m;

	/++
	Constructor
	Params:
		rng = uniform random number generator
		lambda = lambda parameter
		eta = sqrt(chi / psi)
		omega = sqrt(chi * psi)
	+/
	this(ref UniformRNG rng, T lambda, T eta, T omega)
	in {
		assert(eta.isNormal);
		assert(eta > 0);
		assert(omega.isNormal);
		assert(omega > 0);
	}
	body {
		import std.numeric : findRoot;

		this.rng = &rng;
		this.lambdam1 = lambda - 1;
		this.omega= omega;
		this.eta = eta;
		this.m = (lambdam1 + hypot(lambdam1, omega)) / omega;

		immutable a0 = 0.5f * omega;
		immutable a1 = -1 - lambda - m * a0;
		immutable a2 = lambdam1 * m - a0;
		immutable a3 = m * a0;
		//g = function(y) {
		//0.5f*omega*y^3 - y^2*(0.5f*omega*m + lambda + 1) +
		//y*(lambdam1*m - 0.5f*omega) + 0.5f*omega*m
		//}
		T g(T y) { return y * (y * (y * a0 + a1) + a2) + a3; }
		T upper = m;
		while (g(upper) <= 0)
			upper *= 2;
		immutable yM = findRoot(&g, 0, m);
		immutable yP = findRoot(&g, m, upper);
		immutable d = m + 1/m;
		immutable e = -0.25f * omega;

		this.a = (yP - m) * pow(yP/m, 0.5f * lambdam1) * exp(e * (yP + 1/yP - d));
		this.b = (yM - m) * pow(yM/m, 0.5f * lambdam1) * exp(e * (yM + 1/yM - d));
		this.c = e * d + 0.5f * lambdam1 * log(m);
	}


	T front() @property
	{
		immutable c0 = -0.5f*lambdam1;
		immutable c1 = 0.25f * omega;
		for (;;)
		{
			immutable r1 = rng.uniform01!T;
			immutable r2 = rng.uniform01!T;
			immutable x = m + a * r2 / r1 + b * (1 - r2) / r1;
			if (x > 0 && -log(r1) >= c0 * log(x) + c1 * (x + 1/x) + c)
				return x * eta;
		}
	}
}

///
unittest
{
	import std.range;
	auto rng = new ProperGeneralizedInverseGaussianRNG!double(rndGen, -2, 5.0, 2);
	auto sample = rng.map!(x => x + 4).take(9).array;
}


/++
Class to generate random observations from a
generalized inverse Gaussian 
distribution. 
+/
final class GeneralizedInverseGaussianRNG(T, UniformRNG = Random) : DistributionRNG!T
	if (isFloatingPoint!T)
{
	private DistributionRNG!T rng;
	
	/++
	Constructor
	Params:
		rng = uniform random number generator
		lambda = lambda parameter
		chi = chi parameter
		psi = psi parameter
	+/
	this(ref UniformRNG rng, T lambda, T chi, T psi)
	in {
		assert(lambda.isFinite);
		assert(chi.isFinite);
		assert(chi >= 0);
		assert(psi.isFinite);
		assert(psi >= 0);
	}
	body {
		immutable params = GIGChiPsi!T(chi, psi);
		if (chi <= T.min_normal)
			this.rng = new GammaRNG!(T, UniformRNG)(rng, lambda, 2 / psi);
		else if (psi <= T.min_normal)
			this.rng = new InverseGammaRNG!(T, UniformRNG)(rng, -lambda, chi / 2);
		else if (lambda == -0.5f)
			this.rng = new InverseGaussianRNG!(T, UniformRNG)(rng, params.eta, chi);
		else
			this.rng = new ProperGeneralizedInverseGaussianRNG!(T, UniformRNG)(rng, lambda, params.eta, params.omega);
	}

	T front() @property
	{
		return rng.front;
	}
}

///
unittest
{
	import std.range;
	auto rng = new GeneralizedInverseGaussianRNG!double(rndGen, -2, 5.0, 2);
	auto sample = rng.map!(x => x + 4).take(9).array;
}


/++
Class to generate random observations from a
variance-gamma 
distribution using normal variance-mean mixture of 
gamma
distribution.
+/
final class VarianceGammaRNG(T, UniformRNG = Random) : NormalVarianceMeanMixtureRNG!T
	if (isFloatingPoint!T)
{
	/++
	Constructor
	Params:
		rng = uniform random number generator
		beta = mixture scale parameter: `Y*U^(1/2)+beta*U`
		shape = gamma shape parameter
		scale = gamma scale parameter
	+/
	this(ref UniformRNG rng, T beta, T shape, T scale = 1)
	{
		super(rng, new GammaRNG!(T, UniformRNG)(rng, shape, scale), beta);
	}
}

///
unittest
{
	import std.range;
	auto rng = new VarianceGammaRNG!double(rndGen, 1.1, 1.1, 1.1);
	auto sample = rng.map!(x => x + 4).take(9).array;
}


/++
Class to generate random observations from a
hyperbolic asymmetric 
t-distribution using normal variance-mean mixture of 
inverse-gamma
distribution.
+/
final class HyperbolicAsymmetricTRNG(T, UniformRNG = Random) : NormalVarianceMeanMixtureRNG!T
	if (isFloatingPoint!T)
{
	/++
	Constructor
	Params:
		rng = uniform random number generator
		beta = mixture scale parameter: `Y*U^(1/2)+beta*U`
		shape = inverse-gamma shape parameter
		scale = inverse-gamma scale parameter
	+/
	this(ref UniformRNG rng, T beta, T shape, T scale = 1)
	{
		super(rng, new InverseGammaRNG!(T, UniformRNG)(rng, shape, scale), beta);
	}
}

///
unittest
{
	import std.range;
	auto rng = new HyperbolicAsymmetricTRNG!double(rndGen, 1.1, 1.1, 1.1);
	auto sample = rng.map!(x => x + 4).take(9).array;
}


/++
Class to generate random observations from a
generalized variance-gamma 
distribution using normal variance-mean mixture of 
generalized gamma
distribution.
+/
final class GeneralizedVarianceGammaRNG(T, UniformRNG = Random) : NormalVarianceMeanMixtureRNG!T
	if (isFloatingPoint!T)
{
	/++
	Constructor
	Params:
		rng = uniform random number generator
		beta = mixture scale parameter: `Y*U^(1/2)+beta*U`
		shape = generalized gamma shape parameter
		power = generalized gamma power parameter
		scale = generalized gamma scale parameter
	+/
	this(ref UniformRNG rng, T beta, T shape, T power, T scale = 1)
	{
		super(rng, new GeneralizedGammaRNG!(T, UniformRNG)(rng, shape, power, scale), beta);
	}
}

///
unittest
{
	import std.range;
	auto rng = new GeneralizedVarianceGammaRNG!double(rndGen, 1.1, 1.1, 1.1, 1.1);
	auto sample = rng.map!(x => x + 4).take(9).array;
}


/++
Class to generate random observations from a
normal inverse Gaussian 
distribution using normal variance-mean mixture of 
inverse Gaussian
distribution.
+/
final class NormalInverseGaussianRNG(T, UniformRNG = Random) : NormalVarianceMeanMixtureRNG!T
	if (isFloatingPoint!T)
{
	/++
	Constructor
	Params:
		rng = uniform random number generator
		beta = mixture scale parameter: `Y*U^(1/2) + beta*U`
		mu = inverse Gaussian mu parameter
		lambda = inverse Gaussian lambda parameter
	+/
	this(ref UniformRNG rng, T beta, T mu, T lambda)
	{
		super(rng, new InverseGaussianRNG!(T, UniformRNG)(rng, mu, lambda), beta);
	}
}

///
unittest
{
	import std.range;
	auto rng = new NormalInverseGaussianRNG!double(rndGen, 1.1, 1.1, 1.1);
	auto sample = rng.map!(x => x + 4).take(9).array;
}


/++
Class to generate random observations from a
proper generalized hyperbolic 
distribution using normal variance-mean mixture of 
proper generalized inverse Gaussian
distribution.
+/
final class ProperGeneralizedHyperbolicRNG(T, UniformRNG = Random) : NormalVarianceMeanMixtureRNG!T
	if (isFloatingPoint!T)
{
	/++
	Constructor
	Params:
		rng = uniform random number generator
		lambda = proper generalized inverse Gaussian lambda parameter
		beta = mixture scale parameter: `Y*U^(1/2)+beta*U`
		eta = proper generalized inverse Gaussian eta parameter
		omega = proper generalized inverse Gaussian omega parameter
	+/
	this(ref UniformRNG rng, T beta, T lambda, T eta, T omega)
	{
		super(rng, new ProperGeneralizedInverseGaussianRNG!(T, UniformRNG)(rng, lambda, eta, omega), beta);
	}
}

///
unittest
{
	import std.range;
	auto rng = new ProperGeneralizedHyperbolicRNG!double(rndGen, 1.1, 1.1, 1.1, 1.1);
	auto sample = rng.map!(x => x + 4).take(9).array;
}


/++
Class to generate random observations from a
generalized hyperbolic 
distribution using normal variance-mean mixture of 
generalized inverse Gaussian
distribution.
+/
final class GeneralizedHyperbolicRNG(T, UniformRNG = Random) : NormalVarianceMeanMixtureRNG!T
	if (isFloatingPoint!T)
{
	/++
	Constructor
	Params:
		rng = uniform random number generator
		lambda = generalized inverse Gaussian lambda parameter
		beta = mixture scale parameter: `Y*U^(1/2)+beta*U`
		chi = generalized inverse Gaussian chi parameter
		psi = generalized inverse Gaussian psi parameter
	+/
	this(ref UniformRNG rng, T lambda, T beta, T chi, T psi)
	{
		super(rng, new GeneralizedInverseGaussianRNG!(T, UniformRNG)(rng, lambda, chi, psi), beta);
	}
}

///
unittest
{
	import std.range;
	auto rng = new GeneralizedHyperbolicRNG!double(rndGen, 1.1, 1.1, 1.1, 1.1);
	auto sample = rng.map!(x => x + 4).take(9).array;
}



/++
Returns: random number from `[-1, +1] \ {0}`.
+/
T uniformM11E0(T = double)() 
	if (isFloatingPoint!T)
{
	return uniformM11E0!T(rndGen);
}

///ditto
T uniformM11E0(T = double, UniformRNG)(ref UniformRNG rng) 
	if (isFloatingPoint!T && isUniformRNG!UniformRNG)
{
	immutable x = rng.uniform01!T - 1;
	return uniform(0, 2, rng) ? -x : +x;
}

///
unittest
{
	auto x = uniformM11E0();
}


/++
Function to generate random observationb from
standard normal 
distribution.
+/
T rNormal(T = double)() 
	if (isFloatingPoint!T)
{
	return rndGen.rNormal!T;
}

///ditto
T rNormal(T = double, UniformRNG)(ref UniformRNG rng) 
	if (isFloatingPoint!T && isUniformRNG!UniformRNG)
{
	T p, p0, p1;
	do 
	{
		p0 = rng.uniformM11E0!T;
		p1 = rng.uniformM11E0!T;
		p = p0 * p0 + p1 * p1;
	} 
	while (p >= 1);
	return p0 * sqrt(-2 * log(p)/p);
}

///
unittest
{
	auto x = rNormal() * 5 + 1;
}


/++
Function to generate random observation from
standard exponential 
distribution.
+/
T rExponential(T = double)() 
	if (isFloatingPoint!T)
{
	return rndGen.rExponential!T;
}

///ditto
T rExponential(T = double, UniformRNG)(ref UniformRNG rng) 
	if (isFloatingPoint!T && isUniformRNG!UniformRNG)
{
	return -rng.uniform01!T.log;
}

///
unittest
{
	auto x = rExponential() * 5;
}


/++
Function to generate random observation from a
gamma 
distribution.

References
	"Computer Generation of Statistical Distributions" by Richard Saucier
+/
T rGamma(T = double)(T shape) 
	if (isFloatingPoint!T)
{
	return rndGen.rGamma!T(shape);
}

///ditto
T rGamma(T = double, UniformRNG)(ref UniformRNG rng, T shape) 
	if (isFloatingPoint!T && isUniformRNG!UniformRNG)
in {
	assert(shape.isNormal);
	assert(shape > 0);
}
body {
	enum T E = 4.5;
	enum T D = 2.5040773967762740733732583523868748412194809812852436493487L; //1 + log( E );

	if (shape == 1) 
		return rng.rExponential!T;

	immutable A = 1 / sqrt( 2. * shape - 1 );
	immutable B = shape - log( 4. );
	immutable Q = shape + 1 / A;
	immutable C = 1 + shape / E;
	if (shape < 1) 
		for (;;)
		{
			T p = C * rng.uniform01!T;
			T y = void, f = void;
			if ( p > 1 )
			{
				y = -log( ( C - p ) / shape );
				f = pow(y, shape - 1);
			}
			else 
			{
				y = pow( p, 1 / shape );
				f = exp(-y);
			}
			if ( rng.uniform01!T <=  f) 
				return y;
		}
	else 
		for (;;)
		{
			T p0 = rng.uniform01!T;
			T p1 = rng.uniform01!T;
			T v = A * log( p0 / ( 1. - p0 ) );
			T y = shape * exp( v );
			T z = p0 * p0 * p1;
			T w = B + Q * v - y;
			if ( w + D - E * z >= 0.0 || w >= log(z) ) 
				return y;
		}
}

///
unittest
{
	auto x = rGamma(2.0) * 5;
}


/++
Function to generate random observation from a
generalized gamma 
distribution.
+/
T rGeneralizedGamma(T = double)(T shape, T power) 
	if (isFloatingPoint!T)
{
	return rndGen.rGeneralizedGamma!T(shape, power);
}

///ditto
T rGeneralizedGamma(T = double, UniformRNG)(ref UniformRNG rng, T shape, T power) 
	if (isFloatingPoint!T && isUniformRNG!UniformRNG)
in {
	assert(power.isNormal);
	assert(shape.isNormal);
	assert(power > 0);
	assert(shape > 0);
}
body {
	return rng.rGamma!T(shape).pow(1/power);
}

///
unittest
{
	auto x = rGeneralizedGamma(2.0, 3.0) * 5;
}


/++
Function to generate random observation from a
inverse-gamma 
distribution.
+/
T rInverseGamma(T = double)(T shape) 
	if (isFloatingPoint!T)
{
	return rndGen.rInverseGamma!T(shape);
}

///ditto
T rInverseGamma(T = double, UniformRNG)(ref UniformRNG rng, T shape) 
	if (isFloatingPoint!T && isUniformRNG!UniformRNG)
in {
	assert(shape.isNormal);
	assert(shape > 0);
}
body {
	return 1 / rng.rGamma!T(shape);
}

///
unittest
{
	auto x = rInverseGamma(2.0) * 5;
}


/++
Function to generate random observation from a
inverse Gaussian
distribution.

References:
	Michael, John R.; Schucany, William R.; Haas, Roy W. (May 1976). "Generating Random Variates Using Transformations with Multiple Roots".
 +/
T rInverseGaussian(T = double)(T mu, T lambda)
	if (isFloatingPoint!T)
{
	return rndGen.rInverseGaussian!T(mu, lambda);
}

///ditto
T rInverseGaussian(T = double, UniformRNG)(ref UniformRNG rng, T mu, T lambda)
	if (isFloatingPoint!T && isUniformRNG!UniformRNG)
in {
	assert(mu.isNormal);
	assert(mu > 0);
	assert(lambda.isNormal);
	assert(lambda > 0);
}
body {
	immutable nu = mu * mu;
	immutable v = rng.rNormal!T;
	lambda *= 2;
	immutable y = v * v * nu / lambda;
	immutable x = mu + (y - sqrt(y * (mu + y)));
	return rng.uniform01!T <= mu / (mu + x) ? x : nu / x;
}

///
unittest
{
	auto x = rInverseGaussian(2.0, 3.0) * 5;
}


/++
Function to generate random observation from a
Chi-squared
distribution.
 +/
T rChiSquare(T = double)(T shape) 
	if (isFloatingPoint!T)
{
	return rndGen.rChiSquare!T(shape);
}

///ditto
T rChiSquare(T = double, UniformRNG)(ref UniformRNG rng, T shape) 
	if (isFloatingPoint!T && isUniformRNG!UniformRNG)
in {
	assert(shape.isNormal);
	assert(shape > 0);
}
body {
	return rng.rGamma(shape / 2) * 2;
}

///
unittest
{
	auto x = rChiSquare(2.0) * 5;
}


/++
Function to generate random observation from a
Student's t-distribution.
 +/
T rStudentT(T = double)(T shape) 
	if (isFloatingPoint!T)
{
	return rndGen.rStudentT!T(shape);
}

///ditto
T rStudentT(T = double, UniformRNG)(ref UniformRNG rng, T shape) 
	if (isFloatingPoint!T && isUniformRNG!UniformRNG)
in {
	assert(shape.isNormal);
	assert(shape > 0);
}
body {
	return rng.rNormal!T / sqrt(rng.rChiSquare!T(shape)/shape);
}

///
unittest
{
	auto x = rStudentT(2.0) * 5;
}


/++
Function to generate random observation from a
Weibull
distribution.
+/
T rWeibull(T = double)(T shape) 
	if (isFloatingPoint!T)
{
	return rndGen.rWeibull!T(shape);
}

///ditto
T rWeibull(T = double, UniformRNG)(ref UniformRNG rng, T power) 
	if (isFloatingPoint!T && isUniformRNG!UniformRNG)
in {
	assert(power.isNormal);
	assert(power > 0);
}
body {
	return  rng.rExponential!T.pow(1/power);
}

///
unittest
{
	auto x = rWeibull(2.0) * 5 + 1;
}
