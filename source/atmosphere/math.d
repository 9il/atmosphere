/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.math;

//import core.stdc.tgmath;
import std.range;
import std.traits;
import std.stdio;
import std.math : approxEqual;

import std.stdio;


//import std.math: nextUp, nextDown;
import std.math;


//version(LDC)
//{
//	import ldc.intrinsics: log = llvm_log, log2 = llvm_log2, fabs = llvm_fabs;
//}
//else
//{
//	import std.math: log, log2, fabs;
//}

//version (LDC)
//{
//	pragma(LDC_intrinsic, "llvm.fmuladd.f#")
//		T llvm_fmuladd(T)(T vala, T valb, T valc) @safe pure nothrow @nogc;
//}

/**
Computes accurate sum of binary logarithms of input range $(D r).
 */
ElementType!Range sumOfLog2s(Range)(Range r)
	if (isInputRange!Range && isFloatingPoint!(ElementType!Range))
{
	import std.math: frexp;
	long exp = 0;
	Unqual!(typeof(return)) x = 1;
	foreach (e; r)
	{
		if (e < 0)
			return typeof(return).nan;
		int lexp = void;
		x *= frexp(e, lexp);
		exp += lexp;
		if (x < 0.5)
		{
			x *= 2;
			exp--;
		}
	}
	return exp + log2(x);
}

//Issue 14231
private T findRoot(T, DF, DT)(scope DF f, in T a, in T b,
	scope DT tolerance) //= (T a, T b) => false)
	if(
		isFloatingPoint!T &&
		is(typeof(tolerance(T.init, T.init)) : bool) &&
		is(typeof(f(T.init)) == R, R) && isFloatingPoint!R
	)
{
	import std.numeric: findRoot;
	import std.functional;
	immutable fa = f(a);
	if(fa == 0)
		return a;
	immutable fb = f(b);
	if(fb == 0)
		return b;
	immutable r = findRoot!(Unqual!T, Unqual!T)(f.toDelegate, a, b, fa, fb, tolerance.toDelegate);
	// Return the first value if it is smaller or NaN
	return !(fabs(r[2]) > fabs(r[3])) ? r[0] : r[1];
}

//Issue 14231
private T findRoot(T, DF)(scope DF f, in T a, in T b)
{
	return findRoot(f, a, b, (T a, T b) => false);
}

/** Inverse of the Log Minus Digamma function
 * 
 *   Returns x such $(D log(x) - digamma(x) == y).
 *
 * References:
 *   1. Abramowitz, M., and Stegun, I. A. (1970).
 *	 Handbook of mathematical functions. Dover, New York,
 *	 pages 258-259, equation 6.3.18.
 * 
 * Authors: Ilya Yaroshenko
 */
T logmdigammaInverse(T)(const T y)
{
	//Issue 14231
	//import std.numeric: findRoot;
	enum maxY = 1 / T.min_normal;
	static assert(maxY > 0 && maxY <= T.max);
	if (y >= maxY)
		return 1 / y; //lim x->0 (log(x)-digamma(x))*x == 1
	if (y < 0)
		return T.nan;
	if (y < 2*T.min_normal)
		return 0.5 / y; //6.3.18
	if (y > 0)
		// x/2 <= logmdigamma(1 / x) <= x, x > 0
		// calls logmdigamma ~6 times
		return 1 / findRoot((T x) => logmdigamma(1 / x) - y, y,  2*y); // ~6 times
	return y; //NaN
}


unittest {
	alias T = double;
	import std.range : iota;
	import std.math : approxEqual, nextDown;
	//import std.mathspecial: logmdigamma;
	foreach (x; iota(1.3L, 2L, 0.1L))
	{
		assert(logmdigammaInverse(logmdigamma(x)).approxEqual(x));
		assert(logmdigammaInverse(logmdigamma(1/x)).approxEqual(1/x));
	}
	//WolframAlpha, 22.02.2015
	immutable T[2][5] testData = [
		[1.0L, 0.615556766479594378978099158335549201923L],
		[1.0L/8, 4.15937801516894947161054974029150730555L],
		[1.0L/1024, 512.166612384991507850643277924243523243L],
		[0.000500083333325000003968249801594877323784632117L, 1000.0L],
		[1017.644138623741168814449776695062817947092468536L, 1.0L/1024],
	];
	foreach(test; testData)
	{
		import std.stdio;
		writeln(logmdigammaInverse(test[0]));
		assert(approxEqual(logmdigammaInverse(test[0]), test[1], 1e-10, 0));
	}
	assert(approxEqual(logmdigamma(logmdigammaInverse(1.0)), 1, 1e-15, 0));
	assert(approxEqual(logmdigamma(logmdigammaInverse(T.min_normal)), T.min_normal, 1e-15, 0));
	assert(approxEqual(logmdigamma(logmdigammaInverse(T.max/2)), T.max/2, 1e-15, 0));
	assert(approxEqual(logmdigammaInverse(logmdigamma(1.0)), 1, 1e-10, 0));
	assert(approxEqual(logmdigammaInverse(logmdigamma(T.min_normal)), T.min_normal, 1e-15, 0));
	assert(approxEqual(logmdigammaInverse(logmdigamma(T.max/2)), T.max/2, 1e-15, 0));
}

T logmdigamma(T)(const T x)
{
	static immutable T [7] Bn_n  = [
	1.0L/(6*2), -1.0L/(30*4), 1.0L/(42*6), -1.0L/(30*8),
	5.0L/(66*10), -691.0L/(2730*12), 7.0L/(6*14) ];

	if (x <= 0.0)
	{
		if (x == 0.0)
		{
		   return T.infinity;
	}
		return T.nan;
	}

	T s = x;
	T w = 0f;
	while ( s < 10f ) {
		w += 1f/s;
		s += 1f;
	}
	T y;
	if ( s < 1.0e17 ) {
		immutable T z = 1.0/(s * s);
		version(none)
		{
			y = z * 
				llvm_fmuladd(z,
				llvm_fmuladd(z,
				llvm_fmuladd(z,
				llvm_fmuladd(z,
				llvm_fmuladd(z,
				llvm_fmuladd(z,
				Bn_n[6] ,
				Bn_n[5]),
				Bn_n[4]),
				Bn_n[3]),
				Bn_n[2]),
				Bn_n[1]),
				Bn_n[0]);
		}
		else
		{
			y =   z * (Bn_n[0]
				+ z * (Bn_n[1]
				+ z * (Bn_n[2]
				+ z * (Bn_n[3]
				+ z * (Bn_n[4]
				+ z * (Bn_n[5]
				+ z *  Bn_n[6]))))));
		}
	} else
		y = 0f;

	return x == s ? y + 0.5f/s : (log(x/s) + 0.5f/s + y + w);
}

unittest {
	import std.typetuple;
	import std.math: isNaN, isIdentical, NaN, approxEqual;
	import std.mathspecial: digamma;
	foreach(T; TypeTuple!(real, double, float))
	{
		alias lmd = logmdigamma!T;
		assert(lmd(-5.0).isNaN());
		assert(isIdentical(lmd(NaN(0xABC)), NaN(0xABC)));
		assert(lmd(0.0) == T.infinity);
		for(auto x = 0.01; x < 1.0; x += 0.1)
			assert(approxEqual(digamma(x), log(x) - lmd(x)));
		for(auto x = 1.0; x < 15.0; x += 1.0)
			assert(approxEqual(digamma(x), log(x) - lmd(x)));		
	}
}



T besselKR(T)(T x, T lambda)
{
	import bessel;
	immutable b = besselK(x, lambda+1, Flag!"ExponentiallyScaled".yes);
	immutable c = besselK(x, lambda  , Flag!"ExponentiallyScaled".yes);
	return b / c;
}


import std.math : sqrt, isNormal, isFinite, isNaN;
alias fmin = findLocalMin!(double, double delegate(double));
import std.typecons;
import std.stdio;

import std.math : frexp, ldexp;

/++
Find a real minimum of a real function $(D f(x)) via bracketing. 
Given a function $(D f) and a range $(D (ax..bx)), 
returns the value of $(D x) in the range which is closest to a minimum of $(D f(x)).
$(D f) is never evaluted at the endpoints of $(D ax) and $(D bx).
If $(D f(x)) has more than one minimum in the range, one will be chosen arbitrarily. 
If $(D f(x)) returns NaN or -Infinity, $(D (x, f(x), NaN)) will be returned;
otherwise, this algorithm is guaranteed to succeed.

Params:
	f = Function to be analyzed
	ax = Left bound of initial range of f known to contain the minimum.
	bx = Right bound of initial range of f known to contain the minimum.
	relTolerance = Relative tolerance.
	absTolerance = Absolute tolerance.

Preconditions:
	$(D ax) and $(D bx) shall be finite reals. $(BR)
	$(D relTolerance) shall be normal positive real. $(BR)
	$(D absTolerance) shall be normal positive real no less then $(T.epsilon*2).

Returns:
	A tuple consisting of $(D x), $(D y = f(x)) and $(D error = 3 * (absTolerance * fabs(x) + relTolerance)).

	The method used is a combination of golden section search and
successive parabolic interpolation. Convergence is never much slower
than that for a Fibonacci search.

References:
	"Algorithms for Minimization without Derivatives", Richard Brent, Prentice-Hall, Inc. (1973)

See_Also: $(LREF findRoot), $(XREF math, isNormal)

Authors: Ilya Yaroshenko (adaptive golden section search)
+/

///
unittest
{
	int i;
	double f(double x)
	{
		i++;
		return fabs(x-1);
	}

	//slow
	auto minD = findLocalMin(&f, -double.max/2, double.max/2, double.min_normal, 2*double.epsilon);
	writeln(minD);
	writeln(10 * double.epsilon);
	writeln(i);
	with(minD)
	{
		assert(approxEqual(x, 1));
		assert(approxEqual(y, 0));
		assert(error <= 10 * double.epsilon);
		//assert(minD.error == 0);
	}
}

unittest
{
	auto ret = findLocalMin((double x) => double.init, 0.0, 1.0, double.min_normal, 2*double.epsilon);
	assert(!ret.x.isNaN);
	assert(ret.y.isNaN);
	assert(ret.error.isNaN);
}

unittest
{
	auto ret = findLocalMin((double x) => log(x), 0.0, 1.0, double.min_normal, 2*double.epsilon);
	assert(ret.error < 3.00001 * ((2*double.epsilon)*fabs(ret.x)+ double.min_normal));
	assert(ret.x >= 0 && ret.x <= ret.error);
}


unittest
{
	size_t i;
	auto ret = findLocalMin((double x) {i++; return log(x);}, 0.0, double.max, double.min_normal, 2*double.epsilon);
	writeln(i);
	assert(ret.y < -18);
	assert(ret.error < 5e-08);
	assert(ret.x >= 0 && ret.x <= ret.error);
}

//unittest
//{
//	size_t i;
//	auto ret = findLocalMin((double x) {i++; return fabs(x);}, -double.max, double.max, double.min_normal, 2*double.epsilon);
//	writeln(i);
//	assert(ret.x.approxEqual(0));
//	assert(ret.y.approxEqual(0));
//	assert(ret.error.approxEqual(0));
//}

unittest
{
	size_t i;
	auto ret = findLocalMin((double x) {i++; return -fabs(x);}, -1.0, 1.0, double.min_normal, 2*double.epsilon);
	writeln(i);
	assert(ret.x.fabs.approxEqual(1));
	assert(ret.y.fabs.approxEqual(1));
	assert(ret.error.approxEqual(0));
}


//unittest
//{
//	size_t i;
//	auto ret = findLocalMin((double x) {i++; return fabs(x);}, -double.max, double.max, double.min_normal, 2*double.epsilon);
//	writeln(i);
//	assert(ret.x.approxEqual(0));
//	assert(ret.y.approxEqual(0));
//	assert(ret.error.approxEqual(0));
//}

unittest
{
	size_t i;
	auto ret = findLocalMin((double x) {i++; return -fabs(x);}, -1.0, 1.0, double.min_normal, 2*double.epsilon);
	writeln(i);
	assert(ret.x.fabs.approxEqual(1));
	assert(ret.y.fabs.approxEqual(1));
	assert(ret.error.approxEqual(0));
}


version(none)
// check speed regressions and infinity loops
unittest
{
	import std.typetuple;
	import std.math;
	size_t i = 0;
	R f00(T, R)(T x) { i++; return 0; }
	R f01(T, R)(T x) { i++; return 1; }
	R f02(T, R)(T x) { i++; return R.min_normal/2; } //subnormal
	R f03(T, R)(T x) { i++; return R(PI); }
	R f04(T, R)(T x) { i++; return log2(cast(R)x); }
	R f05(T, R)(T x) { i++; return log(cast(R)x); }
	R f06(T, R)(T x) { i++; return exp2(cast(R)x); }
	R f07(T, R)(T x) { i++; return exp(cast(R)x); }
	R f08(T, R)(T x) { i++; return sqrt(cast(R)x); }
	R f09(T, R)(T x) { i++; return cbrt(cast(R)x); }
	R f10(T, R)(T x) { i++; return x; }
	R f11(T, R)(T x) { i++; return cast(R)(x)^^2; }
	R f12(T, R)(T x) { i++; return cast(R)(x)^^PI; }
	R f13(T, R)(T x) { i++; return cast(R)(x)^^(1/PI); }
	R f14(T, R)(T x) { i++; return sin(cast(R)x); }
	R f15(T, R)(T x) { i++; return cos(cast(R)x); }
	R f16(T, R)(T x) { i++; return sqrt(abs(cast(R)(x)^^2 - 1)); }
	R f17(T, R)(T x) { i++; return sqrt(abs(cast(R)(x)^^2 - 1)); }
	R f18(T, R)(T x) { i++; return floor(exp(cast(R)x)); } //multiminimum
	R f19(T, R)(T x) { i++; return floor(log(cast(R)x)); } //multiminimum
	R f20(T, R)(T x) { i++; return floor(cast(R)x); }      //multiminimum
	// vars for global checks
	int s1, s2;
	int c1, c2, c;
	foreach(T; TypeTuple!(
		float,
		double,
		real,
		)) // 1
	foreach(R; TypeTuple!(
		float,
		double,
		real,
		)) // 2
	{
		immutable ar1 = [ 
			&f00!(T, R),
			&f01!(T, R),
			&f02!(T, R),
			&f03!(T, R),
			&f04!(T, R),
			&f05!(T, R),
			&f06!(T, R),
			&f07!(T, R),
			&f08!(T, R),
			&f09!(T, R),
			&f10!(T, R),
			&f11!(T, R),
			&f12!(T, R),
			&f13!(T, R),
			&f14!(T, R),
			&f15!(T, R),
			&f16!(T, R),
			&f17!(T, R),
			&f18!(T, R),
			&f19!(T, R),
			&f20!(T, R),
		];
		
		foreach(relTol; [
			T.min_normal, 
			T.epsilon, 
			sqrt(T.epsilon), 
			T.epsilon^^0.25,
			]) // 3
		foreach(absTol; [
			2*T.epsilon, 
			sqrt(T.epsilon),
			]) // 4
		{
			foreach(sign; [-1, 1]) // 5
			foreach(rBound; [1, T.max]) // 6
			foreach(shift; [0, -1, 1]) // 7
			foreach(j, f; ar1) // 8
			{
				auto m2 = findLocalMin!T((T x) => sign * f(x-shift), shift, rBound+shift, relTol, absTol);
			}
		}
	}
}

Tuple!(T, "x", Unqual!(ReturnType!DF), "y", T, "error") 
	findLocalMin(T, DF)(
		scope DF f,
		in T ax,
		in T bx,
		in T relTolerance = T.epsilon^^0.25,
		in T absTolerance = sqrt(T.epsilon),
		)
	if(isFloatingPoint!T)
in
{
	assert(isFinite(ax) && isFinite(bx), "ax and bx shall be finite reals");
	assert(isNormal(absTolerance) && absTolerance >= T.epsilon*2, 
		"absTolerance shall be normal positive real no less then $(T.epsilon*2)");
	assert(isNormal(relTolerance) && relTolerance > 0, 
		"relTolerance shall be normal positive real");
}
out(result)
{
	assert(isFinite(result.x));
}
body
{
	alias R = Unqual!(CommonType!(ReturnType!DF, T));

	// c is the squared inverse of the golden ratio
	// (3 - sqrt(5))/2
	// Value obtained from Wolfram Alpha.
	immutable T c = 0x0.61c8864680b583ea0c633f9fa31237p+0L;
	immutable T cm1 = 0x0.9e3779b97f4a7c15f39cc0605cedc8p+0L;

	R tolerance;

	T a = ax > bx ? bx : ax;
	T b = ax > bx ? ax : bx;

	//sequence of declarations suitable for SIMD instructions

	T  v = a * cm1 + b * c;
	assert(isFinite(v));

	R fv = f(v);

	//R fv = f(v);
	if(isNaN(fv) || fv == -T.infinity)
		return typeof(return)(v, fv, T.init);

	T  w = v;
	R fw = fv;

	T  x = v;
	R fx = fv;

	for(R d = 0, e = 0;;)
	{
		T m = (a + b) / 2;

		// This fix is not part of original algorithm
		if(!isFinite(m)) // fix infinity loop. Issue can be reproduced in R.
		{
			m = a / 2 + b / 2;
			if(!isFinite(m)) // fast-math compiler switch is enabled
			{
				//SIMD instructions can be used by compiler, do not reduce declarations
				int a_exp = void;
				int b_exp = void;
				immutable an = frexp(a, a_exp);
				immutable bn = frexp(b, b_exp);
				immutable am = ldexp(an, a_exp-1);
				immutable bm = ldexp(bn, b_exp-1);
				m = am + bm;
				if(!isFinite(m)) // wrong input: constraints is disabled in release mode
				{
					return typeof(return).init;
				}
			}
		}
		assert(isFinite(m));
		tolerance = absTolerance * fabs(x) + relTolerance;
		immutable t2 = tolerance * 2;
		//assert(x + tolerance < b);
		//assert(x - tolerance > a);


		// check stopping criterion
		if (!(fabs(x - m) > t2 - (b - a) / 2))
			break;

		R p = 0;
		R q = 0;
		R r = 0;

		// fit parabola
		if (fabs(e) > tolerance)
		{
			immutable  xw =  x -  w;
			immutable fxw = fx - fw;
			immutable  xv =  x -  v;
			immutable fxv = fx - fv;
			immutable xwfxv = xw * fxv;
			immutable xvfxw = xv * fxw;
			p = xv * xvfxw - xw * xwfxv;
			q = (xvfxw - xwfxv) * 2;
			if (q > 0)
				p = -p; 
			else
				q = -q;
			r = e;
			e = d;
		}

		T u;
		// a parabolic-interpolation step
		if (fabs(p) < fabs(q * r / 2) && p > q * (a - x) && p < q * (b - x))
		{
			d = p / q;
			u = x + d;
			// f must not be evaluated too close to a or b
			if (u - a < t2 || b - u < t2)
				d = x < m ? tolerance : -tolerance;
		}
		// a golden-section step
		else
		{
			e = (x < m ? b : a) - x;
			d = c * e;
		}

		// f must not be evaluated too close to x
		u = x + (fabs(d) >= tolerance ? d : d > 0 ? tolerance : -tolerance);
		immutable fu = f(u);
		if(isNaN(fu) || fu == -T.infinity)
			return typeof(return)(u, fu, T.init);

		//  update  a, b, v, w, and x
		if (fu <= fx)
		{
			u < x ? b : a = x;
			v = w; fv = fw;
			w = x; fw = fx;
			x = u; fx = fu;
		}
		else
		{
			u < x ? a : b = u;
			if (fu <= fw || w == x) 
			{
				v = w; fv = fw;
				w = u; fw = fu;
			} 
			else if (fu <= fv || v == x || v == w) 
			{ // do not remove this braces
				v = u; fv = fu;
			}
		}
	}
	return typeof(return)(x, fx, tolerance * 3);
}

double chebev(in double[] c, in double x) 
{


	double d = 0;
	double dd=0;
	for (size_t j=c.length-1;j>0;j--)
	{
		immutable sv = d;
		d = 2 * x * d - dd + c[j];
		dd = sv;
	}
    return x * d - dd + c[0] / 2;
}


immutable double[7] c1 = [
	-1.142022680371168e0,
	6.5165112670737e-3,
	3.087090173086e-4,
	-3.4706269649e-6,
	6.9437664e-9,
	3.67795e-11,
	-1.356e-13,
	];
immutable double[8] c2 = [
	1.843740587300905e0,
	-7.68528408447867e-2,
	1.2719271366546e-3,
	-4.9717367042e-6,
	-3.31261198e-8,
	2.423096e-10,
	-1.702e-13,
	-1.49e-15,
	];


unittest
{
	import std.math;
	writeln(besselIK(1, 2));
	assert(approxEqual(besselIK(1, 2)[2], 0.139865881816522427284));
	assert(approxEqual(besselIK(1, 20)[2], 5.88305796955703817765028e-10));
	assert(approxEqual(besselIK(1.2, 1)[2], 0.701080));
	assert(approxEqual(besselIK(10, 400)[2], 1.359308280937876101049800424530298700140591995078599604e-175));
	writeln(besselIK(1, 1e+100)[2]);
}

// Sets io, ko, ipo, and kpo respectively to the Bessel functions I .x/, K .x/ and their derivatives I0 .x/, K0 .x/, 
// for positive x and for xnu D  0. 
// The relative accuracy is within one or two significant digits of EPS. 
// FPMIN is a number close to the machine’s smallest floating-point number.
auto  besselIK(in double nu, in double x)
{
	//writeln("nu = ", nu, " x = ", x);
	import std.math : abs, sin, sinh, cosh, exp, floor;
	immutable size_t MAXIT=10000;
	immutable double EPS= double.epsilon;
	immutable double FPMIN= double.min_normal / EPS;
	immutable double XMIN=2.0, PI=3.141592653589793;

	double a,a1,b,c,d,del1,delh,dels,e,f,fact,fact2,ff,gam1,gam2,
	gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1,ripl,
	ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2,xx;

	int i,l,nl;
	if (x <= 0.0 || nu < 0.0) 
		throw new Exception("bad arguments in besselik");
	nl = cast(int)floor(nu+0.5);
	xmu=nu-nl;
	//nl is the number of downward re- xmu=nu-nl; currences of the I ’s and upward
	//recurrences of K’s. xmu lies be- tween 􏰀1=2 and 1/2.
	//besselfrac.h
	//Copyright 2007 Numerical Recipes Software
	//4
	//Bessel Function Implementations
	xmu2 = xmu * xmu;
	xi=1.0/x;
	xi2=2.0*xi;
	h=nu*xi;
	if (h < FPMIN)
	{
		h=FPMIN;
	}
	b=xi2*nu;
	d=0.0;
	c=h;
	for (i=0;i<MAXIT;i++)
	{
		b += xi2;
		d=1.0/(b+d);
		c=b+1.0/c;
		immutable del=c*d;
		h=del*h;
		if (abs(del-1.0) <= EPS)
			break; 
		//writeln(del);
	}

	//Evaluate CF1 by modified Lentz’s method (÷5.2).
	//Denominators cannot be zero here, so no need for special precau- tions.
	//if (i >= MAXIT)
	//	throw new Exception("x too large in besselik; try asymptotic expansion");
	//else
		//writeln("i = ", i);
	ril=FPMIN;
	ripl=h*ril;
	ril1=ril;
	rip1=ripl;
	fact=nu*xi;
	for (l=nl-1;l >= 0;l--)
	{
		ritemp=fact*ril+ripl;
		fact -= xi;
		ripl=fact*ritemp+ril;
		ril=ritemp;
	}
	f=ripl/ril;
	//Initialize I and I0 for downward recurrence.
	//Store values for later rescaling.
	//Now have unnormalized I and I0 . Use series.
	if (x < XMIN)
	{
		x2=0.5*x;
		pimu=PI*xmu;
		fact = (abs(pimu) < EPS ? 1.0 : pimu/sin(pimu));
		d = -log(x2);
		e=xmu*d;
		fact2 = (abs(e) < EPS ? 1.0 : sinh(e)/e);
		xx=8.0*xmu^^2-1.0;
		gam1=chebev(c1, xx);
		gam2=chebev(c2, xx);
		gampl= gam2-xmu*gam1;
		gammi= gam2+xmu*gam1;
		ff=fact*(gam1*cosh(e)+gam2*fact2*d);
		sum=ff;
		e=exp(e);
		p=0.5*e/gampl;
		q=0.5/(e*gammi);
		c=1.0;
		d=x2*x2;
		sum1=p;
		for (i=1;i<=MAXIT;i++) {
			ff=(i*ff+p+q)/(i*i-xmu2);
			c *= (d/i);
			p /= (i-xmu);
			q /= (i+xmu);
			immutable del=c*ff;
			sum += del;
			del1=c*(p-i*ff);
			sum1 += del1;
			if (abs(del) < abs(sum)*EPS)
				break;
		}
		if (i > MAXIT) 
			throw new Exception("besselK series failed to converge");

		rkmu=sum;
		rk1=sum1*xi2;
	}
	else
	{
		b=2.0*(1.0+x);
		d=1.0/b;
		h=delh=d;
		q1=0.0;
		q2=1.0;
		a1=0.25-xmu2;
		q=c=a1;
		a = -a1;
		s=1.0+q*delh;
		for (i=1;i<MAXIT;i++)
		{
			//Bessel Function Implementations 5
			//Evaluate CF2 by Steed’s algorithm (÷5.2), which is OK because there can be no zero denominators.
			//Initializations for recurrence (6.6.35). First term in equation (6.6.34).
			a -= 2*i;
			c = -a*c/(i+1.0);
			qnew=(q1-b*q2)/a;
			q1=q2;
			q2=qnew;
			q += c*qnew;
			b += 2.0;
			d=1.0/(b+a*d);
			delh=(b*d-1.0)*delh;
			h += delh;
			dels=q*delh;
			s += dels;
			if (abs(dels/s) <= EPS)
				break;
			//Need only test convergence of sum since CF2 itself converges more quickly.
		}
		if (i >= MAXIT) 
			throw new Exception("besselIK: failure to converge in cf2");
		//else
		//	writeln("i = ", i);
		h=a1*h;
		rkmu=sqrt(PI/(2.0*x))*exp(-x)/s;
		rk1=rkmu*(xmu+x+0.5-h)*xi;
	}
	rkmup=xmu*xi*rkmu-rk1;
	rimu=xi/(f*rkmu-rkmup);
	immutable io=(rimu*ril1)/ril;
	immutable ipo=(rimu*rip1)/ril;
	for (i=1;i <= nl;i++)
	{
		rktemp=(xmu+i)*xi2*rk1+rkmu;
		rkmu=rk1;
		rk1=rktemp;
	}
	immutable ko=rkmu;
	immutable kpo=nu*xi*rkmu-rk1;
	immutable xik = x;
	immutable nuik = nu;
	return tuple(io, ipo, ko, kpo, xik, nuik);
}


double modifiedBesselCF2(double nu, double x)
in {
	assert(fabs(nu) <= 0.5f);
	assert(x >= 2);
	assert(isFinite(x));
}
out(result)
{
	//writeln("modifiedBesselCF2 = ", result);
}
body {
	//writeln("nu = ", nu, " x = ", x);
	//immutable size_t MAXIT=10000;
	//immutable nl = cast(int)floor(nu+0.5);
	//immutable xmu=nu-nl;
	//nl is the number of downward re- xmu=nu-nl; currences of the I ’s and upward
	//recurrences of K’s. xmu lies be- tween 􏰀1=2 and 1/2.
	//besselfrac.h
	//Copyright 2007 Numerical Recipes Software
	//4
	//Bessel Function Implementations
	//immutable xmu2 = xmu * xmu;
	double b = 2 * (1 + x);
	if(isInfinity(b))
		return 1;
	double d = 1 / b;
	double h = d;
	double delh = d;
	double a1 = 0.25f - nu^^2;
	double a = -a1;
	uint i;
	//for (i=1;i<MAXIT;i++)
	do
	{
		//Bessel Function Implementations 5
		//Evaluate CF2 by Steed’s algorithm (÷5.2), which is OK because there can be no zero denominators.
		//Initializations for recurrence (6.6.35). First term in equation (6.6.34).
		//writeln("a = ", a);
		//writeln("b = ", b);
		a -= 2 * i;
		b += 2;
		d = 1 / (b + a * d);
		//writeln("d = ", d);
		delh = (b * d - 1.0) * delh;
		h += delh;
		//writeln("delh = ", delh);
		//writeln("h = ", h);
		//writeln(fabs(delh / h));
		assert(!isNaN(fabs(delh / h)));
	}
	while(fabs(delh / h) > double.epsilon);
	//if (i >= MAXIT) 
	//	throw new Exception("modifiedBesselCF2: failure to converge");
	//else
	//	writeln("modifiedBesselCF2: i = ", i);
	return (nu + 0.5 + x - a1 * h) / x;
}



double[2] modifiedBesselCF2Full(double nu, double x)
in {
	assert(fabs(nu) <= 0.5f);
	assert(x >= 2);
	assert(isFinite(x));
}
out(result)
{
	//writeln("modifiedBesselCF2 = ", result);
}
body {
	//writeln("nu = ", nu, " x = ", x);
	//immutable size_t MAXIT=10000;
	//immutable nl = cast(int)floor(nu+0.5);
	//immutable xmu=nu-nl;
	//nl is the number of downward re- xmu=nu-nl; currences of the I ’s and upward
	//recurrences of K’s. xmu lies be- tween 􏰀1=2 and 1/2.
	//besselfrac.h
	//Copyright 2007 Numerical Recipes Software
	//4
	//Bessel Function Implementations
	//immutable xmu2 = xmu * xmu;

	double b = 2 * (1 + x);
	if(isInfinity(b))
		return [nu < 0.5f ? 0 : 1.2533141373155003, 1];
	double d = 1 / b;
	double h = d;
	double delh = d;
	double a1 = 0.25f - nu^^2;
	double a = -a1;
	uint i;
	double q1 = 0;
	double q2 = 1;
	double q = a1;
	double c = a1;
	double s = 1 + q*delh;
	double dels = void;
	do
	{
		//Bessel Function Implementations 5
		//Evaluate CF2 by Steed’s algorithm (÷5.2), which is OK because there can be no zero denominators.
		//Initializations for recurrence (6.6.35). First term in equation (6.6.34).

		a -= 2 * i;
		c = -a * c / (i + 1);
		immutable qnew = (q1 - b * q2) / a;
		q1 = q2;
		q2 = qnew;
		q += c * qnew;
		b += 2;
		d = 1 / (b + a *  d);
		delh = (b * d - 1) * delh;
		dels = q * delh;
		h += delh;
		s += dels;
		//Need only test convergence of sum since CF2 itself converges more quickly.
	}
	while(fabs(dels / s) > double.epsilon);

	enum double SQRTPI_2 = sqrt(PI / 2);
	return [SQRTPI_2 / s, (nu + 0.5 + x - a1 * h) / x];
}


unittest
{
	import std.math;
	//writeln(modifiedBesselCF2(0.4, 1));
	//assert(approxEqual(modifiedBesselCF2(0.4, 1), 1.87491));
	//assert(approxEqual(modifiedBesselCF2(0.499, 1), 1.99872));
	//assert(approxEqual(modifiedBesselCF2(0, 1), 1.429625398260401758028108023450365630080017018192662983459034));

	//assert(0.5.nextDown+0.5 < 1);
	writeln(modifiedBesselCF2(0.5, 4.0.nextUp));
	writeln(modifiedBesselCF2(0.49999999, 4.00001));
	writeln(modifiedBesselCF2(0, 4.00001));
	//writeln(besselIK(1, 2.00001)[2] / besselIK(0, 2.00001)[2]);

	assert(approxEqual(modifiedBesselCF2(0.4, 2.00001), 1.44212));
	assert(approxEqual(modifiedBesselCF2(0.5, 2.00001), 1.5));
	assert(approxEqual(modifiedBesselCF2(0, 2.00001), 1.22804));

	//assert(approxEqual(modifiedBesselCF2(0.4, 1000), 1.00100));
	//assert(approxEqual(modifiedBesselCF2(0.499, 1000), 1.00090));
	//assert(approxEqual(modifiedBesselCF2(0, 1000), 1.000499875124805092705355767307116656100795175108335495396004));

	assert(approxEqual(modifiedBesselCF2(0.4, 1000), 1.00100));
	assert(approxEqual(modifiedBesselCF2(0.499, 1000), 1.00090));
	assert(approxEqual(modifiedBesselCF2(0, 1000), 1.000499875124805092705355767307116656100795175108335495396004));
}

/++
Returns:
	[K(x, nu), 2 / x * K(x, nu+1)]
+/
double[2] modifiedBesselTemmeSeriesImpl(double EPS = double.epsilon)(double nu, double x)
in {
	assert(fabs(nu) <= 0.5f);
	assert(x <= 2);
}
body {
	import std.math;
	immutable size_t MAXIT=10000;
	immutable double x2=0.5*x;
	immutable double pimu=PI*nu;
	immutable double nu2 = nu^^2;
	immutable double fact = (abs(pimu) < EPS ? 1.0 : pimu / sin(pimu));
	double d = -log(x2);
	double e = nu * d;
	immutable double fact2 = (abs(e) < EPS ? 1.0 : sinh(e) / e);
	immutable double xx=8.0*nu^^2-1.0;
	immutable double gam1=chebev(c1, xx);
	immutable double gam2=chebev(c2, xx);
	immutable double gampl= gam2-nu*gam1;
	immutable double gammi= gam2+nu*gam1;
	double ff=fact*(gam1*cosh(e)+gam2*fact2*d);
	double sum=ff;
	e = exp(e);
	double p=0.5*e/gampl;
	double q=0.5/(e*gammi);
	double c=1.0;
	d = x2*x2;
	double sum1=p;
	int i;
	for (i=1;i<=MAXIT;i++) {
		ff=(i*ff+p+q)/(i*i-nu2);
		c *= (d/i);
		p /= (i-nu);
		q /= (i+nu);
		immutable del=c*ff;
		sum += del;
		immutable del1=c*(p-i*ff);
		sum1 += del1;
		if (abs(del) < abs(sum)*EPS)
			break;
	}
	if (i > MAXIT) 
		throw new Exception("besselK series failed to converge");
	//else
	//	writeln("modifiedBesselTemmeSeriesImpl: i = ", i);
	return [sum, sum1];
}

/++
Returns:
	2 / x * K(x, nu+1) / K(x, nu)
+/
double modifiedBesselTemmeSeries(double EPS = double.epsilon)(double nu, double x)
in {
	assert(fabs(nu) <= 0.5f);
	assert(x <= 2);
}
body {
	immutable sums = modifiedBesselTemmeSeriesImpl!EPS(nu, x);
	return sums[1] / sums[0] / x * 2;
}

unittest
{
	import std.math;
	writeln(modifiedBesselTemmeSeries(0.4, 1));
	assert(approxEqual(modifiedBesselTemmeSeries(0.4, 1), 1.87491));
	assert(approxEqual(modifiedBesselTemmeSeries(0.499, 1), 1.99872));
	assert(approxEqual(modifiedBesselTemmeSeries(0, 1), 1.429625398260401758028108023450365630080017018192662983459034));

	//assert(approxEqual(modifiedBesselCF2(0.4, 2.00001), 1.44212));
	//assert(approxEqual(modifiedBesselCF2(0.5, 2.00001), 1.5));
	//assert(approxEqual(modifiedBesselCF2(0, 2.00001), 1.22804));

}


double besselKD(double nu, double x)
in {

}
body {
	immutable anu = fabs(nu);
	int nl = cast(int)floor(anu+0.5);
	double mu=anu-nl;
	//writeln("nl = ", nl);
	//auto l = mu;

	if(x >= 2)
	{
		double r = modifiedBesselCF2(mu, x);
		mu *= 2;
		if(nl)
		{
			double d;
			do
			{
				d = r;
				mu += 2;
				r = mu / x + 1 / r;
			}
			while(--nl);
			return r / d;
		}
		else
		{
			return r * (r - mu / x);
		}
	}
	else
	{
		immutable sums = modifiedBesselTemmeSeriesImpl(mu, x);
		double r = sums[1] / sums[0];
		if(nl)
		{
			immutable x2_4 = x * x / 4;
			double d;
			do
			{
				d = r;
				mu += 1;
				r = mu + x2_4 / r;
			}
			while(--nl);
			return r / d;
		}
		else
		{
			//writeln("SSSS");
			return r / x * (r - mu) / x * 4;
		}
	}
}

unittest {
	assert(approxEqual(besselKD(2, 3), 1.296650971433360224728795947829976056477719397636727742896));
	assert(approxEqual(besselKD(2.2, 0.8), 1.624595391806645609233071414803368115353081130470326304255669));
	assert(approxEqual(besselKD(0.4, 3), 1.334349384832101337822320699631685431560403480505548447388101));
	assert(approxEqual(besselKD(0.4, 0.8), 2.275575715445329347201695021084846610619371526255383391028726));

	double[3][] data = [
		[1, 1, 1.888245647341297344338366806469818100073176922119421727713106268163595832389436074069065151380128808],
		[10, 1, 1.099687904173359294358694759872440909233436542663077290833713959244496776538436835904856529604568407],
		[1, 10, 1.110348556190446111957216138204031855721452559868300859853180425692157659153049491662221677392988991],
		[1.49011e-08, 10, 1.11111],
		[1.49011e-08, 1, 36.2755],
		[2.22044e-16, 1, 72.3192],
		[1.49166e-154, 1, 708.628],
		[1.49166e-154, 10, 1.111111111111111],
		[double.min_normal, 1.1, 11],

	];
	foreach(test; data)
	{
		assert(approxEqual(besselKD(test[1], test[0]), test[2]));
	}
	assert(isFinite(besselKD(0.1, double.epsilon)));
	assert(isFinite(besselKD(0, double.epsilon)));
}

unittest {
	import std.numeric;
	import std.datetime;
	import std.algorithm;
	size_t i;
	auto lambdas = iota(-100, 101, 1).map!"a/2".array;
	auto ret = new double[lambdas.length];
	StopWatch sw;
	sw.start;
	foreach(k, lambda; lambdas)
	{
		double f(double x)
		{
			//i++;
			//writeln("lambda = ", lambda);
			//writeln("x = ", x);
			//writeln("l = ", besselKD(lambda, double.min_normal)-1.02);
			//writeln("r = ", besselKD(lambda, double.max)-1.02);
			return besselKD(lambda, x)-1.02;
		}
		//ret[k] = findRoot(&f, double.min_normal, double.max);
		ret[k] = besselKD(lambda, 2);
	}
	sw.stop;
	writeln(ret);
	writeln(cast(Duration)sw.peek);
	//writeln(findRoot(&f, double.min_normal, double.max));
	//writeln("IC = ", i);
}


double logBesselK(const double nu, const double x)
in {

}
body {
	immutable anu = fabs(nu);
	int nl = cast(int)floor(anu+0.5);
	double mu = anu - nl;
	if(x >= 2)
	{
		immutable r2 = modifiedBesselCF2Full(mu, x);
		double r = r2[1];
		double ret = log(r2[0]) + log(x) / 2 - x;
		if(nl)
		{
			mu *= 2;
			double l = 1;
			do
			{
				int ex = void;
				l *= frexp(r, ex);
				ret += ex;
				mu += 2;
				r = mu / x + 1 / r;
			}
			while(--nl);
			ret += log(l);
		}
		return ret;
	}
	else
	{
		immutable sums = modifiedBesselTemmeSeriesImpl(mu, x);
		double r = sums[1] / sums[0];
		double ret = log(sums[0]);
		if(nl)
		{
			ret += nl * log(2 / x);
			immutable x2_4 = x * x / 4;
			double l = 1;
			do
			{
				int ex = void;
				l *= frexp(r, ex);
				ret += ex;
				mu += 1;
				r = mu + x2_4 / r;
			}
			while(--nl);
			ret += log(l);
		}
		return ret;
	}
}
