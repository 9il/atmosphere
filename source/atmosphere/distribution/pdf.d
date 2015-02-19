/++
Probability density functions
+/
/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.distribution.pdf;

//import core.stdc.tgmath;

import bessel;
import std.traits;
import std.numeric;
//import std.math : isNormal, isFinite, M_2_SQRTPI, SQRT2, PI, approxEqual;
import std.mathspecial;
import std.typecons;
import std.algorithm;
import std.array;

import atmosphere.distribution.params;
import atmosphere.distribution.utilities;


/++
Normal PDF
+/
struct NormalSPDF(T)
{
	private T c, mu, sigma2;

	/++
	Params:
		mu = location
		sigma2 = sigma^2
	+/
	this(T mu, T sigma2) pure
	in {
		assert(sigma2 > 0);
		assert(mu.isFinite);
	}
	body {
		c = 1 / sqrt(2*PI*sigma2);
		this.mu = mu;
		this.sigma2 = sigma2;
	}

	///
	T opCall(T x) const
	{
		import core.stdc.tgmath : exp;
		return c * exp((x-mu)^^2 / (-2*sigma2));
	}
}

///
unittest
{
	auto pdf = NormalSPDF!double(1.9, 2.3);
}

/++
Gamma PDF
+/
struct GammaSPDF(T)
	if(isFloatingPoint!T)
{
	import std.mathspecial : gamma;
	private T shapem1, scale, c;

	/++
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

	///
	T opCall(T x) const
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
	auto pdf = GammaSPDF!double(3, 2);
	auto x = pdf(0.1);
	assert(x.isNormal);

	import scid.calculus : integrate;
	auto result = pdf.integrate(-double.infinity, double.infinity);
	assert(fabs(result.value - 1) < result.error);
}


/++
Inverse-gamma PDF
+/
struct InverseGammaSPDF(T)
	if(isFloatingPoint!T)
{
	import std.mathspecial : gamma;
	private T shapem1, scale, c;

	/++
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

	///
	T opCall(T x) const
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
	auto pdf = InverseGammaSPDF!double(3, 4);
	auto x = pdf(0.1);
	assert(x.isNormal);

	import scid.calculus : integrate;
	auto result = pdf.integrate(-double.infinity, double.infinity);
	assert(fabs(result.value - 1) < result.error);
}


/++
Generalized gamma PDF
+/
struct GeneralizedGammaSPDF(T)
	if(isFloatingPoint!T)
{
	import std.mathspecial : gamma;
	private T shapem1, power, scale, c, e;

	/++
	Params:
		shape = shape parameter
		power = power parameter
		scale = scale parameter
	+/
	this(T shape, T power, T scale)
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
		this.c = fabs(power) / gamma(shape) / scale;
		this.e = power * shape - 1;
	}

	///
	T opCall(T x) const
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
	auto pdf = GeneralizedGammaSPDF!double(3, 2, 0.5);
	auto x = pdf(0.1);
	assert(x.isNormal);

	import scid.calculus : integrate;
	auto result = pdf.integrate(-double.infinity, double.infinity);
	assert(fabs(result.value - 1) < result.error);
}


/++
Inverse Gaussian PDF

See_Also: [distribution.params](distribution/params.html)
+/
struct InverseGaussianSPDF(T)
	if(isFloatingPoint!T)
{
	private T omega, chi, psi;

	///
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

	///
	T opCall(T x) const
	{
		return x < 0 ? 0 : sqrt(chi / (2*PI*x^^3)) * exp(omega - (chi / x + psi * x) / 2);
	}
}

///
unittest 
{
	auto pdf = InverseGaussianSPDF!double(3, 2);
	auto x = pdf(0.1);
	assert(x.isNormal);

	import scid.calculus : integrate;
	auto result = pdf.integrate(-double.infinity, double.infinity);
	assert(fabs(result.value - 1) < result.error);
}


/++
Proper generalized inverse Gaussian PDF

See_Also: [distribution.params](distribution/params.html)
+/
struct ProperGeneralizedInverseGaussianSPDF(T)
	if(isFloatingPoint!T)
{
	private T omega, eta, c, lambdam1;

	///
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

	///
	T opCall(T x) const
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
	auto pdf = ProperGeneralizedInverseGaussianSPDF!double(4, 3, 2);
	auto x = pdf(0.1);
	assert(x.isNormal);

	import scid.calculus : integrate;
	auto result = pdf.integrate(-double.infinity, double.infinity);
	assert(fabs(result.value - 1) < result.error);
}


/++
Variance gamma (gamma mixture of normals) PDF

See_Also: [distribution.params](distribution/params.html)
+/
struct VarianceGammaSPDF(T)
	if(isFloatingPoint!T)
{
	import std.mathspecial : gamma;
	private T lambdamh, alpha, beta, mu, c;

	/++
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

	///
	T opCall(T x) const
	{
		immutable y = x - mu;
		immutable z = fabs(y);
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
	auto pdf = VarianceGammaSPDF!double(1.1, 1.1, 1.0, 1.1);
	auto x = pdf(0.1);
	assert(x.isNormal);

	import scid.calculus : integrate;
	auto result = pdf.integrate(-double.infinity, double.infinity);
	assert(fabs(result.value - 1) < result.error);
}


/++
Hyperbolic asymmetric T (inverse gamma mixture of normals) PDF

See_Also: [distribution.params](distribution/params.html)
+/
struct HyperbolicAsymmetricTSPDF(T)
	if(isFloatingPoint!T)
{
	import std.mathspecial : gamma;
	private T lambdamh, beta, delta, mu, c;

	/++
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
		immutable params = GHypAlphaDelta!T(fabs(beta), beta, delta);
		this.c = M_2_SQRTPI / SQRT2 * pow(params.chi / 2, -lambda) / gamma(-lambda);
		assert(c.isNormal);
		this.lambdamh = lambda - 0.5f;
		this.beta = beta;
		this.delta = delta;
		this.mu = mu;
	}

	///
	T opCall(T x) const
	{
		immutable y = x - mu;
		immutable z = hypot(delta, y);
		immutable alpha = fabs(beta);
		immutable a = z * alpha;
		return c
			* pow(z / alpha, lambdamh)
			* besselK(a, lambdamh, Flag!"ExponentiallyScaled".yes)
			* exp(beta * y - a);
	}
}

///
unittest 
{
	auto pdf = HyperbolicAsymmetricTSPDF!double(-1.1, 1.1, 1.1, 1.1);
	auto x = pdf(0.1);
	import std.conv;
	assert(x.isNormal, text(x));
}


/++
Normal-inverse Gaussian (inverse Gaussian mixture of normals) PDF

See_Also: [distribution.params](distribution/params.html)
+/
struct NormalInverseGaussianSPDF(T)
	if(isFloatingPoint!T)
{
	private T alpha, beta, delta, mu, c;

	/++
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

	///
	T opCall(T x) const
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
	auto pdf = NormalInverseGaussianSPDF!double(1.1, 0.8, 1.1, 1.1);
	auto x = pdf(0.1);
	assert(x.isNormal);

	import scid.calculus : integrate;
	auto result = pdf.integrate(-double.infinity, double.infinity);
	assert(fabs(result.value - 1) < result.error);
}


/++
Proper generalized hyperbolic (generalized inverse Gaussian mixture of normals) PDF

See_Also: [distribution.params](distribution/params.html)
+/
struct ProperGeneralizedHyperbolicSPDF(T)
	if(isFloatingPoint!T)
{
	private T lambdamh, alpha, beta, delta, mu, omega, c;

	/++
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

	///
	T opCall(T x) const
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
	auto pdf = ProperGeneralizedHyperbolicSPDF!double(1.1, 1.1, 0.8, 1.1, 1.1);
	auto x = pdf(0.1);
	assert(x.isNormal);

	import scid.calculus : integrate;
	auto result = pdf.integrate(-double.infinity, double.infinity);
	assert(fabs(result.value - 1) < result.error);
}


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
	assert(fabs(result.value - 1) < result.error);
}

///
alias toPDF = convertTo!PDF;

///
unittest
{
	PDF!double pdf = GammaSPDF!double(1, 3).toPDF;
}

/++
Variance-mean mixture of normals
+/
abstract class NormalVarianceMeanMixturePDF(T) : PDF!T
	if(isFloatingPoint!T)
{
	private PDF!T pdf;
	private T beta, mu;
	private T epsRel, epsAbs;
	private const(T)[] subdivisions;

	/++
    Params:
		pdf     = The PDF to _integrate.
		beta = NVMM scale
		mu = NVMM location
		subdivisions     = TODO.
		epsRel  = (optional) The requested relative accuracy.
		epsAbs  = (optional) The requested fabsolute accuracy.
	See_also: [struct Result](https://github.com/kyllingstad/scid/blob/a9f3916526e4bf9a4da35d14a969e1abfa17a496/source/scid/types.d)
	+/
	this(PDF!T pdf, T beta, T mu, const(T)[] subdivisions = null, T epsRel = 1e-6, T epsAbs = 0)
	in {
		assert(subdivisions.all!isFinite);
		assert(subdivisions.all!(s => s > 0));
		assert(subdivisions.isSorted);
		assert(subdivisions.findAdjacent.empty);
	}
	body {
		this.pdf = pdf;
		this.beta = beta;
		this.mu = mu;
		this.epsRel = epsRel;
		this.epsAbs = epsAbs;
		this.subdivisions = subdivisions;
	}

	T opCall(T x)
	{
		import scid.calculus : integrate;
		T f(T z) {
			return exp(-0.5f * (x - (mu + beta * z)) ^^ 2 / z) / sqrt(2*PI*z) * pdf(z);
		}
		T sum = 0;
		T a = 0;
		foreach(s; subdivisions)
		{
			sum += integrate(&f, a, s, epsRel, epsAbs);
			a = s;
		}
		sum += integrate(&f, a, T.infinity);
		return sum;
	}
}

///
unittest
{
	class MyGHypPDF(T) : NormalVarianceMeanMixturePDF!T
	{
		import atmosphere.distribution.moment;
		this(T lambda, GHypEtaOmega!T params, T mu)
		{
			with(params)
			{
				auto pgig = ProperGeneralizedInverseGaussianSPDF!T(lambda, eta, omega);
				auto e = mu + properGeneralizedInverseGaussianMean(lambda, eta, omega);
				super(pgig.convertTo!PDF, params.beta, mu, [e]);				
			}
		}
	}

	import atmosphere.distribution.params;
	immutable double lambda = 2;
	immutable double mu = 0.3;
	immutable params = GHypEtaOmega!double(2, 3, 4);
	auto pghyp = new GeneralizedHyperbolicPDF!double(lambda, params.alpha, params.beta, params.delta, mu);
	auto pnvmm = new MyGHypPDF!double(lambda, params, mu);
	foreach(i; 0..5)
		assert(approxEqual(pghyp(i), pnvmm(i)));
}


/++
Generalized variance-gamma (generalized gamma mixture of normals) PDF
+/
final class GeneralizedVarianceGammaPDF(T) : NormalVarianceMeanMixturePDF!T
	if(isFloatingPoint!T)
{
	/++
	Params:
		shape = shape parameter (generalized gamma)
		power = power parameter (generalized gamma)
		scale = scale parameter (generalized gamma)
		beta = NVMM scale
		mu = NVMM location
	+/
	this(T shape, T power, T scale, T beta, T mu)
	in {
		assert(shape.isNormal);
		assert(shape > 0);
		assert(power.isFinite);
		assert(scale.isNormal);
		assert(scale > 0);
		assert(beta.isFinite);
		assert(mu.isFinite);
	}
	body {
		import atmosphere.distribution.moment : generalizedGammaMean;
		auto pdf  = GeneralizedGammaSPDF!double(shape, power, scale);
		auto e = generalizedGammaMean(shape, power, scale);
		assert(e > 0);
		assert(e.isFinite);
		super(pdf.convertTo!PDF, beta, mu, [e]);
	}
}

///
unittest 
{
	auto pdf = new GeneralizedVarianceGammaPDF!double(1.1, 1.1, 0.9, 1.1, 1.1);
	auto x = pdf(0.1);
	assert(x.isNormal);

	import scid.calculus : integrate;
	auto result = pdf.integrate(-double.infinity, double.infinity);
	assert(fabs(result.value - 1) < result.error);
}


/++
Generalized inverse Gaussian PDF

See_Also: [distribution.params](distribution/params.html)
+/
final class GeneralizedInverseGaussianPDF(T) : PDF!T
	if(isFloatingPoint!T)
{
	private PDF!T pdf;

	///
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
			this.pdf = GammaSPDF!T(lambda, 2 / psi).convertTo!PDF;
		else if (psi <= T.min_normal)
			this.pdf = InverseGammaSPDF!T(-lambda, chi / 2).convertTo!PDF;
		else if (lambda == -0.5f)
			this.pdf = InverseGaussianSPDF!T(params.eta, chi).convertTo!PDF;
		else
			this.pdf = ProperGeneralizedInverseGaussianSPDF!T(lambda, params.eta, params.omega).convertTo!PDF;
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
	assert(fabs(result.value - 1) < result.error);
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
			this.pdf = VarianceGammaSPDF!T(lambda, alpha, beta, mu).convertTo!PDF;
		else if (alpha == fabs(beta))
			this.pdf = HyperbolicAsymmetricTSPDF!T(lambda, beta, delta, mu).convertTo!PDF;
		else if (lambda == -0.5f)
			this.pdf = NormalInverseGaussianSPDF!T(alpha, beta, delta, mu).convertTo!PDF;
		else
			this.pdf = ProperGeneralizedHyperbolicSPDF!T(lambda, alpha, beta, delta, mu).convertTo!PDF;
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
	assert(fabs(result.value - 1) < result.error);
}
