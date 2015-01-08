/++
Probability density functions
+/
module distribution.pdf;

import bessel;
import std.traits;
import std.numeric;
import std.math;
import std.mathspecial;
import std.typecons;

import distribution.params;

/++
Probability density function interface
+/
interface PDF(T)
{
	/**
	Call operator
	*/
	T opCall(T x);
}

///
unittest
{
	import std.traits, std.mathspecial;

	class NormalPDF : PDF!real
	{
		real opCall(real x)
		{
			// 1/sqrt(2 PI)
			enum c = 0.398942280401432677939946L;
			return c * exp(-0.5f * x * x);
		}
	}

	auto pdf = new NormalPDF;
	auto x = pdf(0.1);
	assert(x.isNormal);

	import scid.calculus : integrate;
	auto result = pdf.integrate(-double.infinity, double.infinity);
	assert(abs(result.value - 1) < result.error);
}


/++
Gamma PDF
+/
final class GammaPDF(T) : PDF!T
	if(isFloatingPoint!T)
{
	private T shapem1, scale, c;

	/++
	Constructor
	Params:
		shape = gamma shape parameter
		scale = gamma scale parameter
	+/
	this(T shape, T scale)
	in {
		assert(shape.isNormal);
		assert(shape > 0);
		assert(scale.isNormal);
		assert(scale > 0);
	}
	body {
		this.c = 1 / (gamma(shape) * scale);
		assert(c.isNormal);
		this.shapem1 = shape - 1;
		this.scale = scale;
	}

	T opCall(T x)
	{
		if(x < 0)
			return 0;
		x /= scale;
		return c 
			* pow(x, shapem1) 
			* exp(-x);
	}
}

///
unittest 
{
	auto pdf = new GammaPDF!double(3, 2);
	auto x = pdf(0.1);
	assert(x.isNormal);

	import scid.calculus : integrate;
	auto result = pdf.integrate(-double.infinity, double.infinity);
	assert(abs(result.value - 1) < result.error);
}


/++
Inverse-gamma PDF
+/
final class InverseGammaPDF(T) : PDF!T
	if(isFloatingPoint!T)
{
	private T shapem1, scale, c;

	/++
	Constructor
	Params:
		shape = gamma shape parameter
		scale = gamma scale parameter
	+/
	this(T shape, T scale)
	in {
		assert(shape.isNormal);
		assert(shape > 0);
		assert(scale.isNormal);
		assert(scale > 0);
	}
	body {
		this.c = 1 / (gamma(shape) * scale);
		assert(c.isNormal);
		this.shapem1 = shape + 1;
		this.scale = scale;
	}

	T opCall(T x)
	{
		if(x < 0)
			return 0;
		x /= scale;
		return c 
			* pow(x, -shapem1) 
			* exp(-1/x);
	}
}

///
unittest 
{
	auto pdf = new InverseGammaPDF!double(3, 4);
	auto x = pdf(0.1);
	assert(x.isNormal);

	import scid.calculus : integrate;
	auto result = pdf.integrate(-double.infinity, double.infinity);
	assert(abs(result.value - 1) < result.error);
}


/++
Generalized gamma PDF
+/
final class GeneralizedGammaPDF(T) : PDF!T
	if(isFloatingPoint!T)
{
	private T shapem1, power, scale, c, e;

	/++
	Constructor
	Params:
		shape = shape parameter
		power = power parameter
		scale = scale parameter
	+/
	this(T shape, T power, T scale = 1)
	in {
		assert(shape.isNormal);
		assert(shape > 0);
		assert(power.isFinite);
		assert(scale.isNormal);
		assert(scale > 0);
	}
	body {
		this.power = power;
		this.scale = scale;
		this.c = abs(power) / gamma(shape) / scale;
		this.e = power * shape - 1;
	}

	T opCall(T x)
	{
		if(x < 0)
			return 0;
		x /= scale;
		return c
			* pow(x, e)
			* exp(-pow(x, power));
	}
}

///
unittest 
{
	auto pdf = new GeneralizedGammaPDF!double(3, 2, 0.5);
	auto x = pdf(0.1);
	assert(x.isNormal);

	import scid.calculus : integrate;
	auto result = pdf.integrate(-double.infinity, double.infinity);
	assert(abs(result.value - 1) < result.error);
}


/++
Inverse Gaussian PDF

See_Also: [distribution.params](distribution/params.html)
+/
final class InverseGaussianPDF(T) : PDF!T
	if(isFloatingPoint!T)
{
	private T omega, chi, psi;

	///Constructor
	this(T chi, T psi)
	in {
		assert(chi.isNormal);
		assert(chi > 0);
		assert(psi.isNormal);
		assert(psi > 0);
	}
	body {
		this.chi = chi;
		this.psi = psi;
		this.omega = GIGChiPsi!T(chi, psi).omega;
	}

	T opCall(T x)
	{
		return x < 0 ? 0 : sqrt(chi / (2*PI*x^^3)) * exp(omega - (chi / x + psi * x) / 2);
	}
}

///
unittest 
{
	auto pdf = new InverseGaussianPDF!double(3, 2);
	auto x = pdf(0.1);
	assert(x.isNormal);

	import scid.calculus : integrate;
	auto result = pdf.integrate(-double.infinity, double.infinity);
	assert(abs(result.value - 1) < result.error);
}


/++
Proper generalized inverse Gaussian PDF
+/
final class ProperGeneralizedInverseGaussianPDF(T) : PDF!T
	if(isFloatingPoint!T)
{
	private T omega, eta, c, lambdam1;

	///Constructor
	this(T lambda, T eta, T omega)
	in {
		assert(lambda.isFinite);
		assert(eta.isNormal);
		assert(eta > 0);
		assert(omega.isNormal);
		assert(omega > 0);

	}
	body {
		this.c = 2 * eta * besselK(omega, lambda, Flag!"ExponentiallyScaled".yes);
		this.eta = eta;
		this.omega = omega;
		this.lambdam1 = lambda - 1;
	}

	T opCall(T x)
	{
		if(x <= 0)
			return 0;
		x /= eta;
		return x < 0 ? 0 : pow(x, lambdam1) * exp(omega * (1 - (1/x + x) / 2)) / c;
	}
}

///
unittest 
{
	auto pdf = new ProperGeneralizedInverseGaussianPDF!double(4, 3, 2);
	auto x = pdf(0.1);
	assert(x.isNormal);

	import scid.calculus : integrate;
	auto result = pdf.integrate(-double.infinity, double.infinity);
	assert(abs(result.value - 1) < result.error);
}


/++
Generalized inverse Gaussian PDF

See_Also: [distribution.params](distribution/params.html)
+/
final class GeneralizedInverseGaussianPDF(T) : PDF!T
	if(isFloatingPoint!T)
{
	private PDF!T pdf;

	///Constructor
	this(T lambda, T chi, T psi)
	in {
		assert(lambda.isFinite);
		assert(chi.isNormal);
		assert(chi >= 0);
		assert(psi.isNormal);
		assert(psi >= 0);
	}
	body {
		immutable params = GIGChiPsi!T(chi, psi);
		if (chi <= T.min_normal)
			this.pdf = new GammaPDF!T(lambda, 2 / psi);
		else if (psi <= T.min_normal)
			this.pdf = new InverseGammaPDF!T(-lambda, chi / 2);
		else if (lambda == -0.5f)
			this.pdf = new InverseGaussianPDF!T(params.eta, chi);
		else
			this.pdf = new ProperGeneralizedInverseGaussianPDF!T(lambda, params.eta, params.omega);
	}

	T opCall(T x)
	{
		return pdf(x);
	}
}

///
unittest 
{
	auto pdf = new GeneralizedInverseGaussianPDF!double(3, 2, 1);
	auto x = pdf(0.1);
	assert(x.isNormal);

	import scid.calculus : integrate;
	auto result = pdf.integrate(-double.infinity, double.infinity);
	assert(abs(result.value - 1) < result.error);
}


/++
Variance gamma (gamma mixture of normals) PDF

See_Also: [distribution.params](distribution/params.html)
+/
final class VarianceGammaPDF(T) : PDF!T
	if(isFloatingPoint!T)
{
	private T lambdamh, alpha, beta, mu, c;

	/++
	Constructor
	Params:
		lambda = 
		alpha = 
		beta = 
		delta = 0
		mu = location
	+/
	this(T lambda, T alpha, T beta, T mu)
	in {
		assert(lambda.isNormal);
		assert(lambda > 0);
		assert(alpha.isNormal);
		assert(alpha > 0);
		assert(beta.isFinite);
		assert(beta > -alpha);
		assert(beta < +alpha);
		assert(mu.isFinite);
	}
	body {
		immutable params = GHypAlphaDelta!T(alpha, beta, 0);
		this.c = M_2_SQRTPI / SQRT2 * pow(params.psi / 2, lambda) / gamma(lambda);
		assert(c.isNormal);
		this.lambdamh = lambda - 0.5f;
		this.alpha = alpha;
		this.beta = beta;
		this.mu = mu;
	}

	T opCall(T x)
	{
		immutable y = x - mu;
		immutable z = abs(y);
		immutable a = alpha * z;
		return c 
			* pow(z / alpha, lambdamh)
			* besselK(a, lambdamh, Flag!"ExponentiallyScaled".yes)
			* exp(beta * y - a);
	}
}

///
unittest 
{
	auto pdf = new VarianceGammaPDF!double(1.1, 1.1, 1.0, 1.1);
	auto x = pdf(0.1);
	assert(x.isNormal);

	import scid.calculus : integrate;
	auto result = pdf.integrate(-double.infinity, double.infinity);
	assert(abs(result.value - 1) < result.error);
}


/++
Hyperbolic asymmetric T (inverse gamma mixture of normals) PDF

See_Also: [distribution.params](distribution/params.html)
+/
final class HyperbolicAsymmetricTPDF(T) : PDF!T
	if(isFloatingPoint!T)
{
	private T lambdamh, beta, delta, mu, c;

	/++
	Constructor
	Params:
		lambda = 
		alpha = `|beta|`
		beta = 
		delta = 
		mu = location
	+/
	this(T lambda, T beta, T delta, T mu)
	in {
		assert(lambda.isNormal);
		assert(lambda < 0);
		assert(beta.isFinite);
		assert(delta.isNormal);
		assert(delta > 0);
		assert(mu.isFinite);
	}
	body {
		immutable params = GHypAlphaDelta!T(abs(beta), beta, delta);
		this.c = M_2_SQRTPI / SQRT2 * pow(params.chi / 2, -lambda) / gamma(-lambda);
		assert(c.isNormal);
		this.lambdamh = lambda - 0.5f;
		this.beta = beta;
		this.delta = delta;
		this.mu = mu;
	}

	T opCall(T x)
	{
		immutable y = x - mu;
		immutable z = hypot(delta, y);
		immutable alpha = abs(beta);
		immutable a = z * alpha;
		import std.stdio;
		return c
			* pow(z / alpha, lambdamh)
			* besselK(a, lambdamh, Flag!"ExponentiallyScaled".yes)
			* exp(beta * y - a);
	}
}

///
unittest 
{
	auto pdf = new HyperbolicAsymmetricTPDF!double(-1.1, 1.1, 1.1, 1.1);
	auto x = pdf(0.1);
	import std.conv;
	assert(x.isNormal, text(x));
}


/++
Normal-inverse Gaussian (inverse Gaussian mixture of normals) PDF

See_Also: [distribution.params](distribution/params.html)
+/
final class NormalInverseGaussianPDF(T) : PDF!T
	if(isFloatingPoint!T)
{
	private T alpha, beta, delta, mu, c;

	/++
	Constructor
	Params:
		lambda = -1/2
		alpha = 
		beta = 
		delta = 
		mu = location
	+/
	this(T alpha, T beta, T delta, T mu)
	in {
		assert(alpha.isNormal);
		assert(alpha > 0);
		assert(beta.isFinite);
		assert(beta > -alpha);
		assert(beta < +alpha);
		assert(delta.isNormal);
		assert(delta > 0);
		assert(mu.isFinite);
	}
	body {
		immutable params = GHypAlphaDelta!T(alpha, beta, delta);
		this.c = exp(params.omega) * alpha * delta / PI;
		assert(c.isNormal);
		this.alpha = alpha;
		this.beta = beta;
		this.delta = delta;
		this.mu = mu;
	}

	T opCall(T x)
	{
		immutable y = x - mu;
		immutable z = hypot(delta, y);
		immutable a = z * alpha;
		return c
			/ z
			* besselK(a, 1, Flag!"ExponentiallyScaled".yes)
			* exp(beta * y - a);
	}
}

///
unittest 
{
	auto pdf = new NormalInverseGaussianPDF!double(1.1, 0.8, 1.1, 1.1);
	auto x = pdf(0.1);
	assert(x.isNormal);

	import scid.calculus : integrate;
	auto result = pdf.integrate(-double.infinity, double.infinity);
	assert(abs(result.value - 1) < result.error);
}


/++
Proper generalized hyperbolic (generalized inverse Gaussian mixture of normals) PDF

See_Also: [distribution.params](distribution/params.html)
+/
final class ProperGeneralizedHyperbolicPDF(T) : PDF!T
	if(isFloatingPoint!T)
{
	private T lambdamh, alpha, beta, delta, mu, omega, c;

	/++
	Constructor
	Params:
		lambda =
		alpha = 
		beta = 
		delta = 
		mu = location
	+/
	this(T lambda, T alpha, T beta, T delta, T mu)
	in {
		assert(lambda.isFinite);
		assert(alpha.isNormal);
		assert(alpha > 0);
		assert(beta.isFinite);
		assert(beta > -alpha);
		assert(beta < +alpha);
		assert(delta.isNormal);
		assert(delta > 0);
		assert(mu.isFinite);
	}
	body {
		immutable params = GHypAlphaDelta!T(alpha, beta, delta);
		this.omega = params.omega;
		this.c = M_2_SQRTPI / (SQRT2 * 2)
			* pow(params.eta, -lambda)
			/ besselK(omega, lambda, Flag!"ExponentiallyScaled".yes);
		assert(c.isNormal);
		this.lambdamh = lambda - 0.5f;
		this.alpha = alpha;
		this.beta = beta;
		this.delta = delta;
		this.mu = mu;
	}

	T opCall(T x)
	{
		immutable y = x - mu;
		immutable z = hypot(delta, y);
		immutable a = z * alpha;
		return c
			* pow(z / alpha, lambdamh)
			* besselK(a, lambdamh, Flag!"ExponentiallyScaled".yes)
			* exp(beta * y - (a - omega));
	}
}

///
unittest 
{
	auto pdf = new ProperGeneralizedHyperbolicPDF!double(1.1, 1.1, 0.8, 1.1, 1.1);
	auto x = pdf(0.1);
	assert(x.isNormal);

	import scid.calculus : integrate;
	auto result = pdf.integrate(-double.infinity, double.infinity);
	assert(abs(result.value - 1) < result.error);
}


/++
Generalized hyperbolic (generalized inverse Gaussian mixture of normals) PDF

See_Also: [distribution.params](distribution/params.html)
+/
final class GeneralizedHyperbolicPDF(T) : PDF!T
	if(isFloatingPoint!T)
{
	private PDF!T pdf;

	/++
	Constructor
	Params:
		lambda =
		alpha = 
		beta = 
		delta = 
		mu = location
	+/
	this(T lambda, T alpha, T beta, T delta, T mu)
	in {
		assert(lambda.isFinite);
		assert(alpha.isNormal || alpha == 0);
		assert(alpha >= 0);
		assert(beta.isFinite);
		assert(beta >= -alpha);
		assert(beta <= +alpha);
		assert(delta.isNormal || delta == 0);
		assert(delta >= 0);
		assert(mu.isFinite);
	}
	body {
		if (delta == 0)
			this.pdf = new VarianceGammaPDF!T(lambda, alpha, beta, mu);
		else if (alpha == abs(beta))
			this.pdf = new HyperbolicAsymmetricTPDF!T(lambda, beta, delta, mu);
		else if (lambda == -0.5f)
			this.pdf = new NormalInverseGaussianPDF!T(alpha, beta, delta, mu);
		else
			this.pdf = new ProperGeneralizedHyperbolicPDF!T(lambda, alpha, beta, delta, mu);
	}

	T opCall(T x)
	{
		return pdf(x);
	}
}

///
unittest 
{
	auto pdf = new GeneralizedHyperbolicPDF!double(1.1, 1.1, 0.9, 1.1, 1.1);
	auto x = pdf(0.1);
	assert(x.isNormal);

	import scid.calculus : integrate;
	auto result = pdf.integrate(-double.infinity, double.infinity);
	assert(abs(result.value - 1) < result.error);
}
