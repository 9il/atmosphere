/**
*/
module distribution.pdf;


import bessel;
import std.traits;
import std.numeric;
import std.math;
import std.mathspecial;
import std.typecons;

/**
*/
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
}


/**
*/
final class GammaPDF(T) : PDF!T
	if(isFloatingPoint!T)
{
	private T shapeGamma, shapem1, scale;

	///
	this(T shape, T scale)
	in {
		assert(shape.isNormal);
		assert(shape > 0);
		assert(scale.isNormal);
		assert(scale > 0);
	}
	body {
		this.shapeGamma = shape.gamma;
		this.shapem1 = shape - 1;
		this.scale = scale;
	}

	T opCall(T x)
	{
		if(x < 0)
			return 0;
		x /= scale;
		return pow(x, shapem1) * exp(-x)  / shapeGamma / scale;
	}
}

/**
*/
unittest 
{
	auto pdf = new GammaPDF!double(3, 2);
	auto x = pdf(0.1);
}


/**
*/
final class InverseGammaPDF(T) : PDF!T
	if(isFloatingPoint!T)
{
	private T shapeGamma, shapem1, scale;

	///
	this(T shape, T scale)
	in {
		assert(shape.isNormal);
		assert(shape > 0);
		assert(scale.isNormal);
		assert(scale > 0);
	}
	body {
		this.shapeGamma = shape.gamma;
		this.shapem1 = shape - 1;
		this.scale = scale;
	}

	T opCall(T x)
	{
		if(x < 0)
			return 0;
		x *= scale;
		return pow(x, -shapem1) * exp(-1/x)  / shapeGamma / scale;
	}
}

/**
*/
unittest 
{
	auto pdf = new InverseGammaPDF!double(3, 2);
	auto x = pdf(0.1);
}

/**
*/
final class InverseGaussianPDF(T) : PDF!T
	if(isFloatingPoint!T)
{
	private T c, cchi, cpsi;

	///
	this(T chi, T psi)
	in {
		assert(chi.isNormal);
		assert(chi > 0);
		assert(psi.isNormal);
		assert(psi > 0);
	}
	body {
		this.c = sqrt(chi/(2 * PI)) * exp(sqrt(chi*psi));
		this.cchi = -0.5f * chi;
		this.cpsi = -0.5f * psi;
	}

	T opCall(T x)
	{
		return x < 0 ? 0 : c * exp(cchi / x + cpsi * x);
	}
}

/**
*/
unittest 
{
	auto pdf = new InverseGaussianPDF!double(3, 2);
	auto x = pdf(0.1);
}


/**
*/
final class ProperGeneralizedInverseGaussianPDF(T) : PDF!T
	if(isFloatingPoint!T)
{
	private T omega, eta, c, lambdam1;
	///
	this(T lambda, T omega, T eta = 1)
	in {
		assert(lambda.isFinite);
		assert(omega.isNormal);
		assert(omega > 0);
		assert(eta.isNormal);
		assert(eta > 0);
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

/**
*/
unittest 
{
	auto pdf = new ProperGeneralizedInverseGaussianPDF!double(3, 2);
	auto x = pdf(0.1);
}


/**
*/
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
		if (chi <= T.min_normal)
			this.pdf = new GammaPDF!T(lambda, 2 / psi);
		else if (psi <= T.min_normal)
			this.pdf = new InverseGammaPDF!T(-lambda, chi / 2);
		else if (lambda == -0.5f)
			this.pdf = new InverseGaussianPDF!T(sqrt(chi/psi), chi);
		else
			this.pdf = new ProperGeneralizedInverseGaussianPDF!T(lambda, sqrt(chi*psi), sqrt(chi/psi));
	}


	T opCall(T x)
	{
		return pdf(x);
	}
}


/**
*/
unittest 
{
	auto pdf = new GeneralizedInverseGaussianPDF!double(3, 2, 1);
	auto x = pdf(0.1);
}