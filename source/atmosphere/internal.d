/**
This module contains numeric methods.
Each method aims to do one iteration over the task.

------
problem: p' = argmin f(p), p_i >= 0, Σ_i p_i = 1.

p - mixture weights,
f = u(Wp),
u(ω) - convex function,
W - matrix of features(n rows, k columns),
k - length of mixture weights,
n - length of sample, n may vary (sliding window).
------
*/
/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.internal;

import std.math: isNaN, fabs;
import std.numeric: findRoot;
import std.algorithm;

import atmosphere.utilities;

public import simple_matrix;

/**
mixture = Wp
*/
void mix(T)(Matrix!(const T) WTransposed, in T[] p, T[] mixture)
in
{
	import std.string;
	assert(WTransposed.height);
	assert(WTransposed.width);
	assert(WTransposed.height == p.length, format("p.length(%s) != WTransposed.height(%s)", p.length, WTransposed.height));
	assert(WTransposed.width == mixture.length, format("mixture.length(%s) != WTransposed.height(%s)", mixture.length, WTransposed.width));
}
body
{
	gemv(WTransposed.transposed, p, mixture);
}

/**
One iteration of gradient descent optimization algorithm.
Params:
	grad = ∇u(ω)
	WTransposed = transposed version of W. k rows, n columns.
	p = mixture weight. length = k.
	mixture = Wp, length = n.
	pi = temporary array, length = n.
	xi = temporary array, length = n.
	gamma = temporary array, length = n.
	c = temporary array, length = k.
	tolerance = Defines an early termination condition. 
				Receives the current upper and lower bounds on the root. 
				The delegate must return true when these bounds are acceptable. 
See_Also: $(MREF mix),
	$(STDREF numeric, findRoot)
Preconditions:
	mixture = Wp
*/
void gradientDescentIteration(alias grad, T)
	(
		Matrix!(const T) WTransposed,
		scope T[] p,
		in T[] mixture,
		scope T[] pi,
		scope T[] xi,
		scope T[] gamma,
		scope T[] c,
		in bool delegate(T, T) @nogc nothrow tolerance = (a, b) => false,
	)
in
{
	assert(WTransposed.height);
	assert(WTransposed.width);
	assert(WTransposed.height == p.length);
	assert(WTransposed.height == c.length);
	assert(WTransposed.width == mixture.length);
	assert(WTransposed.width == pi.length);
	assert(WTransposed.width == xi.length);
	assert(WTransposed.width == gamma.length);
}
body
{
	//gemv(WTransposed.transposed, p, mixture);
	grad(mixture, xi); // xi = grad(mixture);
	gemv(WTransposed, xi, c);
	size_t i;
	i = c.length - c.minPos.length;
	const columni = WTransposed[i];
	foreach(j; 0..mixture.length)
		pi[j] = columni[j] - mixture[j];
	immutable theta = gRoot!grad(mixture, pi, columni, xi, gamma, tolerance);
	immutable onemtheta = 1 - theta;
	foreach(ref e; p)
		e *= onemtheta;
	p[i] += theta;
	p.normalize;
}

/**
One iteration of gradient descent optimization algorithm.

Params:
	PartialDerivative = Partial derivative `y` of objective convex function `u`: `du/dω_j = y(ω_j), 1 <= j <= n.
	WTransposed = transposed version of W. n columns, k rows.
	p = mixture weight. length = k.
	mixture = Wp, length = n.
	pi = temporary array, length = n.
	tolerance = Defines an early termination condition. 
				Receives the current upper and lower bounds on the root. 
				The delegate must return true when these bounds are acceptable. 
See_Also: $(MREF mix),
	$(MREF coordinateDescentIteration),
	$(STDREF numeric, findRoot)
Preconditions:
	mixture = Wp
*/
void gradientDescentIterationPartial(alias PartialDerivative, T)
	(
		Matrix!(const T) WTransposed,
		scope T[] p,
		in T[] mixture,
		scope T[] pi,
		scope T[] xi,
		scope T[] c,
		in bool delegate(T, T) @nogc nothrow tolerance = (a, b) => false,
	)
in
{
	assert(WTransposed.height);
	assert(WTransposed.width);
	assert(WTransposed.height == p.length);
	assert(WTransposed.width == mixture.length);
	assert(WTransposed.width == pi.length);
}
body
{
	import std.functional : unaryFun;
	//gemv(WTransposed.transposed, p, mixture);
	// xi = grad(mixture);
	foreach(j; 0..mixture.length)
		xi[j] = unaryFun!PartialDerivative(mixture[j]);
	gemv(WTransposed, xi, c);
	size_t i;
	i = c.length - c.minPos.length;
	const columni = WTransposed[i];
	//foreach(j; 0..mixture.length)
	//	pi[j] = columni[j] - mixture[j];
	immutable theta = gRoot!PartialDerivative(mixture, pi, columni, tolerance);
	immutable onemtheta = 1 - theta;
	foreach(ref e; p)
		e *= onemtheta;
	p[i] += theta;
	p.normalize;
}


/**
k iterations of coordinate descent optimization algorithm.

For better performance permute rows of WTransposed rows and corresponding elements of p.
Similar rows (in context of u) of WTransposed should be held far from each other.
Params:
	grad = ∇u(ω)
	WTransposed = transposed version of W. n columns, k rows.
	p = mixture weight. length = k.
	mixture = Wp, length = n.
	pi = temporary array, length = n.
	xi = temporary array, length = n.
	gamma = temporary array, length = n.
	tolerance = Defines an early termination condition. 
				Receives the current upper and lower bounds on the root. 
				The delegate must return true when these bounds are acceptable. 
See_Also:	$(MREF mix),
	$(MREF coordinateDescentIterationPartial),
	$(STDREF numeric, findRoot)
Preconditions:
	mixture = Wp
*/
void coordinateDescentIteration(alias grad, T)
	(
		Matrix!(const T) WTransposed,
		scope T[] p,
		scope T[] mixture,
		scope T[] pi,
		scope T[] xi,
		scope T[] gamma,
		in bool delegate(T, T) @nogc nothrow tolerance = (a, b) => false,
	)
in
{
	assert(WTransposed.height);
	assert(WTransposed.width);
	assert(WTransposed.height == p.length);
	assert(WTransposed.width == mixture.length);
	assert(WTransposed.width == pi.length);
	assert(WTransposed.width == xi.length);
	assert(WTransposed.width == gamma.length);
}
body
{
	//gemv(WTransposed.transposed, p, mixture);
	foreach(i; 0..p.length)
	{
		auto columni = WTransposed[i];
		foreach(j; 0..mixture.length)
			pi[j] = columni[j] - mixture[j];
		immutable theta = gRoot!grad(mixture, pi, columni, xi, gamma, tolerance);
		if(theta != 0 && !tolerance(0, theta))
		{
			immutable onemtheta = 1 - theta;
			foreach(ref e; p)
				e *= onemtheta;
			p[i] += theta;
			foreach(j; 0..mixture.length)
				mixture[j] = onemtheta * mixture[j] + theta * columni[j];
		}
	}
	p.normalize;
}


/**
Perform k iterations of coordinate descent optimization algorithm.

For better performance permute rows of $(D_PARAM WTransposed) rows and corresponding elements of p.
Similar rows (in context of u) of $(D_PARAM WTransposed) should be held far from each other.
Params:
	PartialDerivative = Partial derivative `y` of objective convex function `u`: `du/dω_j = y(ω_j), 1 <= j <= n.
	WTransposed = transposed version of W. n columns, k rows.
	p = mixture weight. length = k.
	mixture = Wp, length = n.
	pi = temporary array, length = n.
	tolerance = Defines an early termination condition. 
				Receives the current upper and lower bounds on the root. 
				The delegate must return true when these bounds are acceptable. 
See_Also: $(MREF mix),
	$(MREF coordinateDescentIteration),
	$(STDREF numeric, findRoot)
Preconditions:
	mixture = Wp
*/
void coordinateDescentIterationPartial(alias PartialDerivative, T)
	(
		Matrix!(const T) WTransposed,
		scope T[] p,
		scope T[] mixture,
		scope T[] pi,
		in bool delegate(T, T) @nogc nothrow tolerance = (a, b) => false,
	)
in
{
	assert(WTransposed.height);
	assert(WTransposed.width);
	assert(WTransposed.height == p.length);
	assert(WTransposed.width == mixture.length);
	assert(WTransposed.width == pi.length);
}
body
{
	//gemv(WTransposed.transposed, p, mixture);
	foreach(i; 0..p.length)
	{
		auto columni = WTransposed[i];
		//foreach(j; 0..mixture.length)
		//	pi[j] = columni[j] - mixture[j];
		immutable theta = gRoot!PartialDerivative(mixture, pi, columni, tolerance);
		if(theta != 0 && !tolerance(0, theta))
		{
			immutable onemtheta = 1 - theta;
			foreach(ref e; p)
				e *= onemtheta;
			p[i] += theta;
			foreach(j; 0..mixture.length)
				mixture[j] = onemtheta * mixture[j] + theta * columni[j];
		}
	}
	p.normalize;
}

/**
One iteration of Expectation Maximization algorithm.
Params:
	grad = ∇u(ω)
	WTransposed = transposed version of W. k rows, n columns. W[i, j] >= 0.
	p = mixture weight. length = k.
	mixture = Wp, length = n.
	pi = temporary array, length = n.
	c = temporary array, length = k.
See_Also: $(MREF mix)
Preconditions:
	mixture = Wp
*/
void EMIteration(alias grad, T)
	(
		Matrix!(const T) WTransposed,
		scope T[] p,
		in T[] mixture,
		scope T[] pi,
		scope T[] c,
	)
in
{
	assert(WTransposed.height);
	assert(WTransposed.width);
	assert(WTransposed.height == p.length);
	assert(WTransposed.height == c.length);
	assert(WTransposed.width == mixture.length);
	assert(WTransposed.width == pi.length);
	foreach(row; WTransposed)
		foreach(elem; row)
		{
			assert(elem >= 0, "Components must be non negative.");
		}
}
body
{
	//gemv(WTransposed.transposed, p, mixture);
	grad(mixture, pi); // pi = grad(mixture);
	gemv(WTransposed, pi, c);
	p[] *= c[];
	p.normalize;
}

T gRoot
	(
		alias grad,
		T, 
	)
	(
		in T[] mixture,
		scope T[] pi,
		in T[] columni,
		scope T[] xi,
		scope T[] gamma,
		in bool delegate(T, T) @nogc nothrow tolerance = (a, b) => false,
	)
{
	//return .gRoot!("-1/a", T)(mixture, pi, columni);
	T g(T theta)
	{
		foreach(j; 0..mixture.length)
			xi[j] = mixture[j] + theta * pi[j];
		grad(xi, gamma);
		auto ret = dotProduct(gamma, pi);
		return gCorrectioin(ret);
	}

	T g_01(in T[] vec)
	{
		grad(vec, gamma);
		auto ret = dotProduct(gamma, pi);
		return gCorrectioin(ret);
	}
	immutable g0 = g_01(mixture);
	if(g0 > -T.min_normal)
		return 0;
	immutable g1 = g_01(columni);
	if(g1 < +T.min_normal)
		return 1;
	assert(g0 <= 0);
	assert(g1 >= 0);
	auto r = findRoot(&g, cast(T)0.0, cast(T)1.0, g0, g1, tolerance);
	return !(fabs(r[2]) > fabs(r[3])) ? r[0] : r[1];
}

T gRoot
	(
		alias PartialDerivative,
		T, 
	)
	(
		in T[] mixture,
		scope T[] pi,
		in T[] columni,
		in bool delegate(T, T) @nogc nothrow tolerance = (a, b) => false,
	)
out(theta)
{
	assert(theta >= 0);
	assert(theta <= 1);
}
body
{
	immutable g0 = gCorrectioin(mixture.length - dotProductInverse(columni, mixture)); //g_01!PartialDerivative(pi, mixture);
	import std.stdio;
	if(g0 > -T.min_normal)
		return 0;
	immutable g1 = gCorrectioin(dotProductInverse2(columni, mixture, pi) - mixture.length);
	import atmosphere.utilities;
	if(g1 < +T.min_normal)
		return 1;
	assert(g0 <= 0);
	assert(g1 >= 0);
	foreach(j; 0..mixture.length)
		pi[j] = columni[j] - mixture[j];
	T g(T theta) {return .g!PartialDerivative(pi, mixture, theta);}
	auto r = findRoot(&g, cast(T)0.0, cast(T)1.0, g0, g1, tolerance);
	return !(fabs(r[2]) > fabs(r[3])) ? r[0] : r[1];
}

T g(alias PartialDerivative, T)(in T[] pi, in T[] mixture, T theta)
{
	T ret = 0;
	static if(is(typeof(PartialDerivative) == string) && (PartialDerivative == "-1/a"))
	{
		foreach(j; 0..pi.length)
			version(LDC)
				static if(is(Unqual!T == double))
					ret = inlineIR!(`
						%m = fmul fast double %1, %3
						%a = fadd fast double %2, %m
						%d = fdiv fast double %1, %a
						%r = fadd fast double %0, %d
						ret double %r`, double)(ret, pi[j], mixture[j], theta);
				else
				static if(is(Unqual!T == float))
					ret = inlineIR!(`
						%m = fmul fast float %1, %3
						%a = fadd fast float %2, %m
						%d = fdiv fast float %1, %a
						%r = fadd fast float %0, %d
						ret float %r`, float)(ret, pi[j], mixture[j], theta);
				else
					ret += pi[j] / (mixture[j] + theta * pi[j]);
			else
				ret += pi[j] / (mixture[j] + theta * pi[j]);
		ret = -ret;
	}
	else
		foreach(j; 0..pi.length)
			ret += pi[j] * unaryFun!PartialDerivative(mixture[j] + theta * pi[j]);
	return gCorrectioin(ret);
}

T gCorrectioin(T)(T x)
{
	if(x == +T.infinity)
		x = +T.max;
	if(x == -T.infinity)
		x = -T.max;
	assert(!isNaN(x));
	return x;
}
