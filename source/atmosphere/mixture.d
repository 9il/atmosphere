/**
This module contains numeric methods.
Each method aims to do one iteration over the task.

Optimization problem f(p) -> min, 
where p is discrete probability distribution with k elements.
f = u(Wp),
W - matrix(n rows, k columns),
u - convex function.
*/
module atmosphere.mixture;

import core.stdc.stdlib : free;
import core.stdc.tgmath : fabs;

import atmosphere.utilities;


/**
One iteration of gradient descent optimization algorithm.
Params:
	grad = ∇u(ω)
	WT = transposed version of W. k rows, n columns.
	p = discrete probability distribution with, length = k.
	chi = temporary array, length = n.
	pi = temporary array, length = n.
	xi = temporary array, length = n.
	gamma = temporary array, length = n.
	c = temporary array, length = k.
	tolerance = Defines an early termination condition. 
				Receives the current upper and lower bounds on the root. 
				The delegate must return true when these bounds are acceptable. 
*/
void gradientDescentIteration(alias grad, T, Matrix)
	(
		Matrix WT,
		T[] p,
		T[] chi,
		T[] pi,
		T[] xi,
		T[] gamma,
		T[] c,
		in bool delegate(T, T) @nogc nothrow tolerance = (a, b) => false,
	)
in
{
	assert(WT.height);
	assert(WT.width);
	assert(WT.height == p.length);
	assert(WT.height == c.length);
	assert(WT.width == chi.length);
	assert(WT.width == pi.length);
	assert(WT.width == xi.length);
	assert(WT.width == gamma.length);
}
body
{
	gemv(WT.transposed, p, chi);
	grad(chi, xi); // xi = grad(chi);
	gemv(WT, xi, c);
	immutable i = c.length - c.minPos.length;
	foreach(j, w; WT[i])
		pi[j] = w - chi[j];
	immutable theta = gRoot!grad(chi, pi, xi, gamma, tolerance);
	immutable onemtheta = 1 - theta;
	p.scal(onemtheta);
	p[i] += theta;
	p.normalize;
}


/**
One iteration of gradient descent optimization algorithm.
Params:
	simpleGrad = du/dω_1, where du/dω_j = du/dω_1, 1 <= j <= n.
	WT = transposed version of W. n columns, k rows.
	p = discrete probability distribution with, length = k.
	chi = temporary array, length = n.
	pi = temporary array, length = n.
	xi = temporary array, length = n.
	c = temporary array, length = k.
	tolerance = Defines an early termination condition. 
				Receives the current upper and lower bounds on the root. 
				The delegate must return true when these bounds are acceptable. 
*/
void simpleGradientDescentIteration(alias simpleGrad, T, Matrix)
	(
		Matrix WT,
		T[] p,
		T[] chi,
		T[] pi,
		T[] xi,
		T[] c,
		in bool delegate(T, T) @nogc nothrow tolerance = (a, b) => false,
	)
in
{
	assert(WT.height);
	assert(WT.width);
	assert(WT.height == p.length);
	assert(WT.height == c.length);
	assert(WT.width == chi.length);
	assert(WT.width == pi.length);
	assert(WT.width == xi.length);
}
body
{
	gemv(WT.transposed, p, chi);
	foreach(j; 0..chi.length)
		xi[j] = simpleGrad(chi[j]);
	gemv(WT, xi, c);
	immutable i = c.length - c.minPos.length;
	foreach(j, w; WT[i])
		pi[j] = w - chi[j];
	immutable theta = gRoot!simpleGrad(chi, pi, tolerance);
	immutable onemtheta = 1 - theta;
	p.scal(onemtheta);
	p[i] += theta;
	p.normalize;
}

/**
k iterations of coordinate descent optimization algorithm.

For better performance permute rows of WT rows and corresponding elements of p.
Similar rows (in context of u) of WT should be held far from each other.
Params:
	grad = ∇u(ω)
	WT = transposed version of W. n columns, k rows.
	p = discrete probability distribution with, length = k.
	chi = temporary array, length = n.
	pi = temporary array, length = n.
	xi = temporary array, length = n.
	gamma = temporary array, length = n.
	tolerance = Defines an early termination condition. 
				Receives the current upper and lower bounds on the root. 
				The delegate must return true when these bounds are acceptable. 
*/
void coordinateDescentIteration(alias grad, T, Matrix)
	(
		Matrix WT,
		T[] p,
		T[] chi,
		T[] pi,
		T[] xi,
		T[] gamma,
		in bool delegate(T, T) @nogc nothrow tolerance = (a, b) => false,
	)
in
{
	assert(WT.height);
	assert(WT.width);
	assert(WT.height == p.length);
	assert(WT.width == chi.length);
	assert(WT.width == pi.length);
	assert(WT.width == xi.length);
	assert(WT.width == gamma.length);
}
body
{
	gemv(WT.transposed, p, chi);
	foreach(i; 0..p.length)
	{
		auto column = WT[i];
		foreach(j; 0..chi.length)
			pi[j] = column[j] - chi[j];
		immutable theta = gRoot!grad(chi, pi, xi, gamma, tolerance);
		immutable onemtheta = 1 - theta;
		if(theta)
		{
			p.scal(onemtheta);
			p[i] += theta;
			foreach(j; 0..chi.length)
				chi[j] = onemtheta * chi[j] - theta * column[j];
		}
	}
	p.normalize;
}


/**
k iterations of coordinate descent optimization algorithm.

For better performance permute rows of WT rows and corresponding elements of p.
Similar rows (in context of u) of WT should be held far from each other.
Params:
	simpleGrad = du/dω_1, where du/dω_j = du/dω_1, 1 <= j <= n.
	WT = transposed version of W. n columns, k rows.
	p = discrete probability distribution with, length = k.
	chi = temporary array, length = n.
	pi = temporary array, length = n.
	tolerance = Defines an early termination condition. 
				Receives the current upper and lower bounds on the root. 
				The delegate must return true when these bounds are acceptable. 
*/
void simpleCoordinateDescentIteration(alias simpleGrad, T, Matrix)
	(
		Matrix WT,
		T[] p,
		T[] chi,
		T[] pi,
		in bool delegate(T, T) @nogc nothrow tolerance = (a, b) => false,
	)
in
{
	assert(WT.height);
	assert(WT.width);
	assert(WT.height == p.length);
	assert(WT.width == chi.length);
	assert(WT.width == pi.length);
}
body
{
	gemv(WT.transposed, p, chi);
	foreach(i; 0..p.length)
	{
		auto column = WT[i];
		foreach(j; 0..chi.length)
			pi[j] = column[j] - chi[j];
		immutable theta = gRoot!simpleGrad(chi, pi, tolerance);
		immutable onemtheta = 1 - theta;
		if(theta)
		{
			p.scal(onemtheta);
			p[i] += theta;
			foreach(j; 0..chi.length)
				chi[j] = onemtheta * chi[j] - theta * column[j];
		}
	}
	p.normalize;
}


/**
EM like optimization algorithm.
Params:
	grad = ∇u(ω)
	WT = transposed version of W. n columns, k rows.
	p = discrete probability distribution with, length = k.
	chi = temporary array, length = n.
	xi = temporary array, length = n.
	c = temporary array, length = k.
*/
void EMLikeDescentIteration(alias grad, T, Matrix)
	(
		Matrix WT,
		T[] p,
		T[] chi,
		T[] xi,
		T[] c,
	)
	@nogc nothrow
in
{
	assert(WT.height);
	assert(WT.width);
	assert(WT.height == p.length);
	assert(WT.height == c.length);
	assert(WT.width == chi.length);
	assert(WT.width == xi.length);
}
body
{
	gemv(WT.transposed, p, chi);
	grad(chi, xi);
	gemv(WT, xi, c);
	foreach(i, ref elem; p)
		elem *= c[i];
	p.normalize;
}



/**
EM like descent optimization algorithm.
Params:
	simpleGrad = du/dω_1, where du/dω_j = du/dω_1, 1 <= j <= n.
	WT = transposed version of W. n columns, k rows.
	p = discrete probability distribution with, length = k.
	xi = temporary array, length = n.
	c = temporary array, length = k.
*/
void simpleEMLikeDescentIteration(alias simpleGrad, T, Matrix)  
	(
		Matrix WT,
		T[] p,
		T[] xi,
		T[] c,
	)
	@nogc nothrow
in
{
	assert(WT.height);
	assert(WT.width);
	assert(WT.height == p.length);
	assert(WT.height == c.length);
	assert(WT.width == xi.length);
}
body
{
	gemv(WT.transposed, p, xi);
	foreach(ref elem; xi)
		elem = simpleGrad(elem);
	gemv(WT, xi, c);
	foreach(i, ref elem; p)
		elem *= c[i];
	p.normalize;
}


private:

///
T gRoot
	(
		alias grad,
		T, 
	)
	(
		in T[] chi,
		in T[] pi,
		T[] xi,
		T[] gamma,
		in bool delegate(T, T) @nogc nothrow tolerance = (a, b) => false,
	)
{
	T g(T theta)
	{
		foreach(j; 0..chi.length)
			xi[j] = chi[j] + theta * pi[j];
		grad(xi, gamma);
		return dot(gamma, pi);
	}

	T g0 = void, g1 = void;
	g0 = g(0);
	if(g0 < T.min_normal || tolerance(g0, 0))
		return 0;
	g1 = g(1);
	if(g1 >= 1           || tolerance(g1, 1))
		return 1;
	auto r = findRoot(&g, cast(T)0, cast(T)1, g0, g1, tolerance);
	return !(fabs(r[2]) > fabs(r[3])) ? r[0] : r[1];
}


///
T gRoot
	(
		alias simpleGrad,
		T, 
	)
	(
		in T[] chi,
		in T[] pi,
		in bool delegate(T, T) @nogc nothrow tolerance = (a, b) => false,
	)
{
	T g(T theta)
	{
		T ret0 = 0;
		T ret1 = 0;
		T ret2 = 0;
		T ret3 = 0;

		immutable L1= pi.length & -4;
		size_t i;
	 
		for(; i < L1; i += 4)
		{
			immutable pi0 = pi[i+0];
			immutable pi1 = pi[i+1];
			immutable pi2 = pi[i+2];
			immutable pi3 = pi[i+3];
			ret0 += pi0 * simpleGrad(chi[i+0] + theta * pi0);
			ret1 += pi1 * simpleGrad(chi[i+1] + theta * pi1);
			ret2 += pi2 * simpleGrad(chi[i+2] + theta * pi2);
			ret3 += pi3 * simpleGrad(chi[i+3] + theta * pi3);
		}

		for(; i < pi.length; i++)
		{
			immutable pi0 = pi[i+0];
			ret0 += pi0 * simpleGrad(chi[i+0] + theta * pi0);
		}

		return (ret0+ret1)+(ret2+ret3);
	}

	T g0 = void, g1 = void;
	g0 = g(0);
	if(g0 < T.min_normal || tolerance(g0, 0))
		return 0;
	g1 = g(1);
	if(g1 >= 1           || tolerance(g1, 1))
		return 1;
	auto r = findRoot(&g, cast(T)0, cast(T)1, g0, g1, tolerance);
	return !(fabs(r[2]) > fabs(r[3])) ? r[0] : r[1];
}
