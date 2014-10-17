/**
Extern C bindings.
*/
module atmosphere.catmosphere;

import std.functional : toDelegate;

import atmosphere.mixture;

import simple_matrix;

extern(C):


/**
One iteration of gradient descent optimization algorithm.
Params:
	grad = ∇u(ω)
	WTptr = transposed version of W.
	k = number of rows.
	n = number of columns.
	WTshift = count of elements between adjacent rows.
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
void catmosphere_gradientDescentIteration
	(
		in void function(size_t, in double*, double*) @nogc nothrow grad,
		in double* WTptr,
		size_t k,
		size_t n,
		size_t WTshift,
		double* p,
		double* chi,
		double* pi,
		double* xi,
		double* gamma,
		double* c,
		in bool function(double, double) @nogc nothrow tolerance,
	)
{
	gradientDescentIteration!((in a, b) => grad(a.length, a.ptr, b.ptr))
	(
		Matrix!(double, false)(cast(double*)WTptr, k, n, WTshift),
		p[0..k],
		chi[0..n],
		pi[0..n],
		xi[0..n],
		gamma[0..n],
		c[0..k],
		cast(bool delegate(double, double) @nogc nothrow)(a, b) => tolerance(a, b),
	);
}


/**
One iteration of gradient descent optimization algorithm.
Params:
	simpleGrad = du/dω_1, where du/dω_j = du/dω_1, 1 <= j <= n.
	WTptr = transposed version of W.
	k = number of rows.
	n = number of columns.
	WTshift = count of elements between adjacent rows.
	p = discrete probability distribution with, length = k.
	chi = temporary array, length = n.
	pi = temporary array, length = n.
	xi = temporary array, length = n.
	c = temporary array, length = k.
	tolerance = Defines an early termination condition. 
				Receives the current upper and lower bounds on the root. 
				The delegate must return true when these bounds are acceptable. 
*/
void catmosphere_simpleGradientDescentIteration
	(
		in double function(double) @nogc nothrow simpleGrad,
		in double* WTptr,
		size_t k,
		size_t n,
		size_t WTshift,
		double* p,
		double* chi,
		double* pi,
		double* xi,
		double* c,
		in bool function(double, double) @nogc nothrow tolerance,
	)
{
	simpleGradientDescentIteration!simpleGrad
	(
		Matrix!(double, false)(cast(double*)WTptr, k, n, WTshift),
		p[0..k],
		chi[0..n],
		pi[0..n],
		xi[0..n],
		c[0..k],
		cast(bool delegate(double, double) @nogc nothrow)(a, b) => tolerance(a, b),
	);
}

/**
One iteration of gradient descent optimization algorithm.
-----------
u(ω_1, ..., ω_n) = -Σ_j log(ω_j),
du/dω_j, where du/dω_j = -1/ω_j, 1 <= j <= n.
-----------
Params:
	WTptr = transposed version of W.
	k = number of rows.
	n = number of columns.
	WTshift = count of elements between adjacent rows.
	p = discrete probability distribution with, length = k.
	chi = temporary array, length = n.
	pi = temporary array, length = n.
	xi = temporary array, length = n.
	c = temporary array, length = k.
	tolerance = Defines an early termination condition. 
				Receives the current upper and lower bounds on the root. 
				The delegate must return true when these bounds are acceptable. 
*/
void catmosphere_minusSumOfLogs_gradientDescentIteration
	(
		in double* WTptr,
		size_t k,
		size_t n,
		size_t WTshift,
		double* p,
		double* chi,
		double* pi,
		double* xi,
		double* c,
		in bool function(double, double) @nogc nothrow tolerance,
	)
{
	simpleGradientDescentIteration!(a => -1/a)
	(
		Matrix!(double, false)(cast(double*)WTptr, k, n, WTshift),
		p[0..k],
		chi[0..n],
		pi[0..n],
		xi[0..n],
		c[0..k],
		cast(bool delegate(double, double) @nogc nothrow)(a, b) => tolerance(a, b),
	);
}



/**
k iterations of coordinate descent optimization algorithm.

For better performance permute rows of WT rows and corresponding elements of p.
Similar rows (in context of u) of WT should be held far from each other.
Params:
	grad = ∇u(ω)
	WTptr = transposed version of W.
	k = number of rows.
	n = number of columns.
	WTshift = count of elements between adjacent rows.
	p = discrete probability distribution with, length = k.
	chi = temporary array, length = n.
	pi = temporary array, length = n.
	xi = temporary array, length = n.
	gamma = temporary array, length = n.
	tolerance = Defines an early termination condition. 
				Receives the current upper and lower bounds on the root. 
				The delegate must return true when these bounds are acceptable. 
*/
void catmosphere_coordinateDescentIteration
	(
		in void function(size_t, in double*, double*) @nogc nothrow grad,
		in double* WTptr,
		size_t k,
		size_t n,
		size_t WTshift,
		double* p,
		double* chi,
		double* pi,
		double* xi,
		double* gamma,
		in bool function(double, double) @nogc nothrow tolerance,
	)
{
	coordinateDescentIteration!((in a, b) => grad(a.length, a.ptr, b.ptr))
	(
		Matrix!(double, false)(cast(double*)WTptr, k, n, WTshift),
		p[0..k],
		chi[0..n],
		pi[0..n],
		xi[0..n],
		gamma[0..n],
		cast(bool delegate(double, double) @nogc nothrow)(a, b) => tolerance(a, b),
	);
}

/**
k iterations of coordinate descent optimization algorithm.

For better performance permute rows of WT rows and corresponding elements of p.
Similar rows (in context of u) of WT should be held far from each other.
Params:
	simpleGrad = du/dω_1, where du/dω_j = du/dω_1, 1 <= j <= n.
	WTptr = transposed version of W.
	k = number of rows.
	n = number of columns.
	WTshift = count of elements between adjacent rows.
	p = discrete probability distribution with, length = k.
	chi = temporary array, length = n.
	pi = temporary array, length = n.
	tolerance = Defines an early termination condition. 
				Receives the current upper and lower bounds on the root. 
				The delegate must return true when these bounds are acceptable. 
*/
void catmosphere_simpleCoordinateDescentIteration
	(
		in double function(double) @nogc nothrow simpleGrad,
		in double* WTptr,
		size_t k,
		size_t n,
		size_t WTshift,
		double* p,
		double* chi,
		double* pi,
		in bool function(double, double) @nogc nothrow tolerance,
	)
{
	simpleCoordinateDescentIteration!simpleGrad
	(
		Matrix!(double, false)(cast(double*)WTptr, k, n, WTshift),
		p[0..k],
		chi[0..n],
		pi[0..n],
		cast(bool delegate(double, double) @nogc nothrow)(a, b) => tolerance(a, b),
	);
}


/**
k iterations of coordinate descent optimization algorithm.
-----------
u(ω_1, ..., ω_n) = -Σ_j log(ω_j),
du/dω_j, where du/dω_j = -1/ω_j, 1 <= j <= n.
-----------
For better performance permute rows of WT rows and corresponding elements of p.
Similar rows (in context of u) of WT should be held far from each other.
Params:
	WTptr = transposed version of W.
	k = number of rows.
	n = number of columns.
	WTshift = count of elements between adjacent rows.
	p = discrete probability distribution with, length = k.
	chi = temporary array, length = n.
	pi = temporary array, length = n.
	tolerance = Defines an early termination condition. 
				Receives the current upper and lower bounds on the root. 
				The delegate must return true when these bounds are acceptable. 
*/
void catmosphere_minusSumOfLogs_coordinateDescentIteration
	(
		in double* WTptr,
		size_t k,
		size_t n,
		size_t WTshift,
		double* p,
		double* chi,
		double* pi,
		in bool function(double, double) @nogc nothrow tolerance,
	)
{
	simpleCoordinateDescentIteration!(a => -1/a)
	(
		Matrix!(double, false)(cast(double*)WTptr, k, n, WTshift),
		p[0..k],
		chi[0..n],
		pi[0..n],
		cast(bool delegate(double, double) @nogc nothrow)(a, b) => tolerance(a, b),
	);
}


/**
One iteration of gradient descent optimization algorithm.
Params:
	grad = ∇u(ω)
	WTptr = transposed version of W.
	k = number of rows.
	n = number of columns.
	WTshift = count of elements between adjacent rows.
	p = discrete probability distribution with, length = k.
	chi = temporary array, length = n.
	xi = temporary array, length = n.
	c = temporary array, length = k.
*/
void catmosphere_EMLikeDescentIteration
	(
		in void function(size_t, in double*, double*) @nogc nothrow grad,
		in double* WTptr,
		size_t k,
		size_t n,
		size_t WTshift,
		double* p,
		double* chi,
		double* xi,
		double* c,
	)
{
	EMLikeDescentIteration!((in a, b) => grad(a.length, a.ptr, b.ptr))
	(
		Matrix!(double, false)(cast(double*)WTptr, k, n, WTshift),
		p[0..k],
		chi[0..n],
		xi[0..n],
		c[0..k],
	);
}


/**
One iteration of EMLike descent optimization algorithm.
Params:
	simpleGrad = du/dω_1, where du/dω_j = du/dω_1, 1 <= j <= n.
	WTptr = transposed version of W.
	k = number of rows.
	n = number of columns.
	WTshift = count of elements between adjacent rows.
	p = discrete probability distribution with, length = k.
	xi = temporary array, length = n.
	c = temporary array, length = k.
*/
void catmosphere_simpleEMLikeDescentIteration
	(
		in double function(double) @nogc nothrow simpleGrad,
		in double* WTptr,
		size_t k,
		size_t n,
		size_t WTshift,
		double* p,
		double* xi,
		double* c,
	)
{
	simpleEMLikeDescentIteration!simpleGrad
	(
		Matrix!(double, false)(cast(double*)WTptr, k, n, WTshift),
		p[0..k],
		xi[0..n],
		c[0..k],
	);
}

/**
One iteration of EMLike descent optimization algorithm.
-----------
u(ω_1, ..., ω_n) = -Σ_j log(ω_j),
du/dω_j, where du/dω_j = -1/ω_j, 1 <= j <= n.
-----------
Params:
	WTptr = transposed version of W.
	k = number of rows.
	n = number of columns.
	WTshift = count of elements between adjacent rows.
	p = discrete probability distribution with, length = k.
	xi = temporary array, length = n.
	c = temporary array, length = k.
*/
void catmosphere_minusSumOfLogs_EMLikeDescentIteration
	(
		in double* WTptr,
		size_t k,
		size_t n,
		size_t WTshift,
		double* p,
		double* xi,
		double* c,
	)
{
	simpleEMLikeDescentIteration!(a => -1/a)
	(
		Matrix!(double, false)(cast(double*)WTptr, k, n, WTshift),
		p[0..k],
		xi[0..n],
		c[0..k],
	);
}
