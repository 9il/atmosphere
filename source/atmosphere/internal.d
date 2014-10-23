/**
This module contains numeric methods.
Each method aims to do one iteration over the task.

Optimization problem f(p) -> min, 
where p is discrete probability distribution with k elements.
f = u(Wp),
W - matrix(n rows, k columns),
u - convex function.
*/
module atmosphere.internal;

import core.stdc.tgmath : fabs;

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
	p = discrete probability distribution. length = k.
	mixture = Wp, length = n.
	pi = temporary array, length = n.
	xi = temporary array, length = n.
	gamma = temporary array, length = n.
	c = temporary array, length = k.
	tolerance = Defines an early termination condition. 
				Receives the current upper and lower bounds on the root. 
				The delegate must return true when these bounds are acceptable. 
See_Also:
	$(MREF mix),
	$(STDREF numeric, findRoot)
Preconditions:
	mixture = Wp
*/
void gradientDescentIteration(alias grad, T)
	(
		Matrix!(const T) WTransposed,
		T[] p,
		in T[] mixture,
		T[] pi,
		T[] xi,
		T[] gamma,
		T[] c,
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
	immutable i = c.length - c.minPos.length;
	const columni = WTransposed[i];
	foreach(j, w; columni)
		pi[j] = w - mixture[j];
	immutable theta = gRoot!grad(mixture, pi, columni, xi, gamma, tolerance);
	immutable onemtheta = 1 - theta;
	p.scal(onemtheta);
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
	p = discrete probability distribution. length = k.
	mixture = Wp, length = n.
	pi = temporary array, length = n.
	xi = temporary array, length = n.
	gamma = temporary array, length = n.
	tolerance = Defines an early termination condition. 
				Receives the current upper and lower bounds on the root. 
				The delegate must return true when these bounds are acceptable. 
See_Also:
	$(MREF mix),
	$(MREF coordinateDescentIterationPartial),
	$(STDREF numeric, findRoot)
Preconditions:
	mixture = Wp
*/
void coordinateDescentIteration(alias grad, T)
	(
		Matrix!(const T) WTransposed,
		T[] p,
		T[] mixture,
		T[] pi,
		T[] xi,
		T[] gamma,
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
		if(!tolerance(0, theta))
		{
			immutable onemtheta = 1 - theta;
			p.scal(onemtheta);
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
	PartialDerivative = Partial derivative $(D y) of objective convex function $(D u): $(D du/dω_j = y(ω_j), 1 <= j <= n).
	WTransposed = transposed version of W. n columns, k rows.
	p = discrete probability distribution. length = k.
	mixture = Wp, length = n.
	pi = temporary array, length = n.
	tolerance = Defines an early termination condition. 
				Receives the current upper and lower bounds on the root. 
				The delegate must return true when these bounds are acceptable. 
See_Also:
	$(MREF mix),
	$(MREF coordinateDescentIteration),
	$(STDREF numeric, findRoot)
Preconditions:
	mixture = Wp
*/
void coordinateDescentIterationPartial(alias PartialDerivative, T)
	(
		Matrix!(const T) WTransposed,
		T[] p,
		T[] mixture,
		T[] pi,
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
		foreach(j; 0..mixture.length)
			pi[j] = columni[j] - mixture[j];
		immutable theta = gRoot!PartialDerivative(mixture, pi, columni, tolerance);
		if(!tolerance(0, theta))
		{
			immutable onemtheta = 1 - theta;
			p.scal(onemtheta);
			p[i] += theta;
			foreach(j; 0..mixture.length)
				mixture[j] = onemtheta * mixture[j] + theta * columni[j];
		}
	}
	p.normalize;
}



//Returns: estimation of θ.
T gRoot
	(
		alias grad,
		T, 
	)
	(
		in T[] mixture,
		in T[] pi,
		in T[] columni,
		T[] xi,
		T[] gamma,
		in bool delegate(T, T) @nogc nothrow tolerance = (a, b) => false,
	)
{
	T g(T theta)
	{
		foreach(j; 0..mixture.length)
			xi[j] = mixture[j] + theta * pi[j];
		grad(xi, gamma);
		auto ret = dot(gamma, pi);
		return gCorrectioin(ret);
	}

	T g_01(in T[] vec)
	{
		grad(vec, gamma);
		auto ret = dot(gamma, pi);
		return gCorrectioin(ret);
	}

	immutable g0 = g_01(mixture);
	if(g0 > -T.min_normal)
		return 0;
	immutable g1 = g_01(columni);
	if(g1 < +T.min_normal)
		return 1;
	auto r = findRoot(&g, cast(T)0.0, cast(T)1.0, g0, g1, tolerance);
	return !(fabs(r[2]) > fabs(r[3])) ? r[0] : r[1];
}


/**
One iteration of Expectation Maximization algorithm.
Params:
	grad = ∇u(ω)
	WTransposed = transposed version of W. k rows, n columns. W[i, j] >= 0.
	p = discrete probability distribution. length = k.
	mixture = Wp, length = n.
	pi = temporary array, length = n.
	c = temporary array, length = k.
See_Also:
	$(MREF mix)
Preconditions:
	mixture = Wp
*/
void EMIteration(alias grad, T)
	(
		Matrix!(const T) WTransposed,
		T[] p,
		in T[] mixture,
		T[] pi,
		T[] c,
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


//ditto
T gRoot
	(
		alias PartialDerivative,
		T, 
	)
	(
		in T[] mixture,
		in T[] pi,
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
	T g(T theta)
	{
		T ret0 = 0;
		T ret1 = 0;
		T ret2 = 0;
		T ret3 = 0;

		immutable L1= pi.length & -4;
		size_t j;
	 
		for(; j < L1; j += 4)
		{
			immutable pi0 = pi[j+0];
			immutable pi1 = pi[j+1];
			immutable pi2 = pi[j+2];
			immutable pi3 = pi[j+3];
			ret0 += pi0 * PartialDerivative(mixture[j+0] + theta * pi0);
			ret1 += pi1 * PartialDerivative(mixture[j+1] + theta * pi1);
			ret2 += pi2 * PartialDerivative(mixture[j+2] + theta * pi2);
			ret3 += pi3 * PartialDerivative(mixture[j+3] + theta * pi3);
		}

		for(; j < pi.length; j++)
		{
			immutable pi0 = pi[j+0];
			ret0 += pi0 * PartialDerivative(mixture[j+0] + theta * pi0);
		}

		auto ret = (ret0+ret1)+(ret2+ret3);
		return gCorrectioin(ret);
	}

	T g_01(in T[] vec)
	{
		T ret0 = 0;
		T ret1 = 0;
		T ret2 = 0;
		T ret3 = 0;

		immutable L1= pi.length & -4;
		size_t j;
	 
		for(; j < L1; j += 4)
		{
			ret0 += pi[j+0] * PartialDerivative(vec[j+0]);
			ret1 += pi[j+1] * PartialDerivative(vec[j+1]);
			ret2 += pi[j+2] * PartialDerivative(vec[j+2]);
			ret3 += pi[j+3] * PartialDerivative(vec[j+3]);
		}

		for(; j < pi.length; j++)
		{
			ret0 += pi[j+0] * PartialDerivative(vec[j+0]);
		}

		auto ret = (ret0+ret1)+(ret2+ret3);
		return gCorrectioin(ret);
	}

	immutable g0 = g_01(mixture);
	if(g0 > -T.min_normal)
		return 0;
	immutable g1 = g_01(columni);
	if(g1 < +T.min_normal)
		return 1;
	auto r = findRoot(&g, cast(T)0.0, cast(T)1.0, g0, g1, tolerance);
	return !(fabs(r[2]) > fabs(r[3])) ? r[0] : r[1];
}


/*
Returns:
	±T.max for ±∞ and x otherwise.
*/
T gCorrectioin(T)(T x)
{
	if(x == +T.infinity)
		x = +T.max;
	if(x == -T.infinity)
		x = -T.max;
	return x;
}


