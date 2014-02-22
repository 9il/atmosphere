module atmosphere.kernel.pdfkernel;

import core.stdc.tgmath;
import std.math : PI;

import atmosphere.kernel.math;

alias gamma = tgamma;

struct NormalKernel(F)
{
	F mean, scale;

	this(in F mean, in F scale)
	{
		this.mean = mean;
		this.scale = scale;
	}

	auto opCall(in F x) const 
	{ 
		return x.normalPDF!F(mean, scale);
	}
}
 

struct GammaKernel(F)
{
	F shape_m_1, scale;
	typeof(gamma(F.init)) gamma_shape;

	this(in F shape, in F scale)
	{
		this.shape_m_1 = shape - 1;
		this.scale = scale;
		this.gamma_shape = gamma(shape);
	}

	auto opCall(in F x) const
	{
		const y = x / scale;
		return pow(y, shape_m_1 / (gamma_shape * exp(y))) / scale; 
	}
}