/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.math;

//import core.stdc.tgmath;
import std.range;
import std.traits;
import std.typecons;
import std.math;


//version(LDC)
//{
//	import ldc.intrinsics: log = llvm_log, log2 = llvm_log2, fabs = llvm_fabs;
//}
//else
//{
//	import std.math: log, log2, fabs;
//}

version (LDC)
{
	pragma(LDC_intrinsic, "llvm.fmuladd.f#")
		T llvm_fmuladd(T)(T vala, T valb, T valc) @safe pure nothrow @nogc;
}

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
package T findRoot(T, DF, DT)(scope DF f, in T a, in T b,
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
package T findRoot(T, DF)(scope DF f, in T a, in T b)
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
		version(LDC)
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
	$(D absTolerance) shall be normal positive real no less then $(D T.epsilon*2).

Returns:
	A tuple consisting of $(D x), $(D y = f(x)) and $(D error = 3 * (absTolerance * fabs(x) + relTolerance)).

	The method used is a combination of golden section search and
successive parabolic interpolation. Convergence is never much slower
than that for a Fibonacci search.

References:
	"Algorithms for Minimization without Derivatives", Richard Brent, Prentice-Hall, Inc. (1973)

See_Also: $(LREF findRoot), $(XREF math, isNormal)
+/

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
		"absTolerance shall be normal positive real no less then $(D T.epsilon*2)");
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
	assert(ret.y < -18);
	assert(ret.error < 5e-08);
	assert(ret.x >= 0 && ret.x <= ret.error);
}


unittest
{
	size_t i;
	auto ret = findLocalMin((double x) {i++; return -fabs(x);}, -1.0, 1.0, double.min_normal, 2*double.epsilon);
	assert(ret.x.fabs.approxEqual(1));
	assert(ret.y.fabs.approxEqual(1));
	assert(ret.error.approxEqual(0));
}

unittest
{
	size_t i;
	auto ret = findLocalMin((double x) {i++; return -fabs(x);}, -1.0, 1.0, double.min_normal, 2*double.epsilon);
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


Unqual!(CommonType!(T1, T2)) 
chebevReversed(T1, T2)(in T1[] c, in T2 x) 
{
	typeof(return) d = 0;
	typeof(return) dd = 0;
	foreach(e; c[0..$-1])
	{
		immutable sv = d;
		d = 2 * x * d - dd + e;
		dd = sv;
	}
    return x * d - dd + c[$-1] / 2;
}



T modifiedBesselCF2(T)(in T nu, in T x)
	if(isFloatingPoint!T)
in {
	assert(fabs(nu) <= 0.5f);
	assert(x >= 2);
	assert(isFinite(x));
}
body {
	T b = 2 * (1 + x);
	if(isInfinity(b))
		return 1;
	T d = 1 / b;
	T h = d;
	T delh = d;
	T a1 = 0.25f - nu^^2;
	T a = -a1;
	uint i;
	do
	{
		i++;
		a -= 2 * i;
		b += 2;
		d = 1 / (b + a * d);
		delh = (b * d - 1) * delh;
		h += delh;
		assert(!isNaN(fabs(delh / h)));
	}
	while(fabs(delh / h) > T.epsilon);
	return (nu + 0.5 + x - a1 * h) / x;
}



T[2] modifiedBesselCF2Full(T)(in T nu, in T x)
	if(isFloatingPoint!T)
in {
	assert(fabs(nu) <= 0.5f);
	assert(x >= 2);
	assert(isFinite(x));
}
body {
	T b = 2 * (1 + x);
	if(isInfinity(b))
		return [nu < 0.5f ? 0 : 1.2533141373155003, 1];
	T d = 1 / b;
	T h = d;
	T delh = d;
	T a1 = 0.25f - nu^^2;
	T a = -a1;
	T q1 = 0;
	T q2 = 1;
	T q = a1;
	T c = a1;
	T s = 1 + q*delh;
	T dels = void;
	uint i;

	do
	{
		i++;
		a -= 2 * i;
		c = -a * c / (i + 1);
		immutable qnew = (q1 - b * q2) / a;
		q1 = q2;
		q2 = qnew;
		q += c * qnew;
		b += 2;
		d = 1 / (b + a * d);
		delh = (b * d - 1) * delh;
		dels = q * delh;
		h += delh;
		s += dels;
		//Need only test convergence of sum since CF2 itself converges more quickly.
	}
	while(fabs(dels / s) > T.epsilon);

	enum T SQRTPI_2 = sqrt(PI / 2);
	return [SQRTPI_2 / s, (nu + 0.5 + x - a1 * h) / x];
}


unittest
{
	assert(approxEqual(modifiedBesselCF2(0.4, 2.00001), 1.44212, 0.0, 1e-4));
	assert(approxEqual(modifiedBesselCF2(0.5, 2.00001), 1.5, 0.0, 1e-4));
	assert(approxEqual(modifiedBesselCF2(0.0, 2.00001), 1.22804, 0.0, 1e-4));
	assert(approxEqual(modifiedBesselCF2(0.4, 1000.0), 1.00100, 0.0, 1e-3));
	assert(approxEqual(modifiedBesselCF2(0.499, 1000.0), 1.00090, 0.0, 1e-4));
	assert(approxEqual(modifiedBesselCF2(0.0, 1000), 1.000499875124805092705355767307116656100795175108335495396004, 0.0, 1e-14));
}

/++
Returns:
	[K(x, nu), 2 / x * K(x, nu+1)]
+/
T[2] modifiedBesselTemmeSeriesImpl(T)(in T nu, in T x)
in {
	assert(fabs(nu) <= 0.5f);
	assert(x <= 2);
}
body {
	//TODO: recalculate c1 and c2 for different FP formats
	static if(is(T == float))
	{
		static immutable float[4] c1 = [
			-3.4706269649e-6,
			3.087090173086e-4,
			6.5165112670737e-3,
			-1.142022680371168e0,
			];
		static immutable float[8] c2 = [
			-4.9717367042e-6,
			1.2719271366546e-3,
			-7.68528408447867e-2,
			1.843740587300905e0,
			];
	}
	else
	{
		static immutable double[7] c1 = [
			-1.356e-13,
			3.67795e-11,
			6.9437664e-9,
			-3.4706269649e-6,
			3.087090173086e-4,
			6.5165112670737e-3,
			-1.142022680371168e0,
			];
		static immutable double[8] c2 = [
			-1.49e-15,
			-1.702e-13,
			2.423096e-10,
			-3.31261198e-8,
			-4.9717367042e-6,
			1.2719271366546e-3,
			-7.68528408447867e-2,
			1.843740587300905e0,
			];
	}

	immutable T x2 = x / 2;
	immutable T pimu = cast(T) PI * nu;
	immutable T nu2 = nu^^2;
	immutable T fact = (abs(pimu) < T.epsilon ? 1 : pimu / sin(pimu));
	T d = -log(x2);
	T e = nu * d;
	immutable T fact2 = (abs(e) < T.epsilon ? 1 : sinh(e) / e);
	immutable T xx = 8 * nu^^2 - 1;
	immutable T gam1 = chebevReversed(c1, xx);
	immutable T gam2 = chebevReversed(c2, xx);
	immutable T gampl = gam2-nu*gam1;
	immutable T gammi = gam2+nu*gam1;
	T ff = fact * (gam1 * cosh(e) + gam2 * fact2 * d);
	T sum = ff;
	e = exp(e);
	T p = 0.5f * e / gampl;
	T q = 0.5f / (e * gammi);
	T c = 1;
	d = x2^^2;
	T sum1 = p;
	int i;
	T del = void;
	do {
		i++;
		ff = (i * ff + p + q) / (i * i - nu2);
		c *= d / i;
		p /= i - nu;
		q /= i + nu;
		del = c * ff;
		sum += del;
		immutable del1 = c * (p - i * ff);
		sum1 += del1;
	}
	while(abs(del) >= abs(sum) * T.epsilon);
	return [sum, sum1];
}

/++
Returns:
	K(x, nu+1) / K(x, nu)
+/
T modifiedBesselTemmeSeries(T)(in T nu, in T x)
in {
	assert(fabs(nu) <= 0.5f);
	assert(x <= 2);
}
body {
	immutable sums = modifiedBesselTemmeSeriesImpl(nu, x);
	return sums[1] / sums[0] / x * 2;
}

unittest
{
	assert(approxEqual(modifiedBesselTemmeSeries(0.4, 1.0), 1.87491, 0.0, 1e-4));
	assert(approxEqual(modifiedBesselTemmeSeries(0.499, 1.0), 1.9987, 0.0, 1e-4));
	assert(approxEqual(modifiedBesselTemmeSeries(0.0, 1.0), 1.429625398260401758028108023450365630080017018192662983459034, 0.0, 1e-14));
}


T besselKD(T)(in T nu, in T x)
	if(isFloatingPoint!T)
in {

}
body {
	immutable anu = fabs(nu);
	int nl = cast(int)floor(anu+0.5);
	T mu=anu-nl;

	if(x >= 2)
	{
		T r = modifiedBesselCF2(mu, x);
		mu *= 2;
		if(nl)
		{
			T d;
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
		T r = sums[1] / sums[0];
		if(nl)
		{
			immutable x2_4 = x * x / 4;
			T d;
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
			return r / x * (r - mu) / x * 4;
		}
	}
}

unittest {

	assert(approxEqual(besselKD(0.0, 1.0),  2.043828779351212337989476332008573915738881135574432645409779, 0.0, 1e-14));
	assert(approxEqual(besselKD(0.0, 2.0), 1.508074700999049512999886217349628893603326842760291910336296, 0.0, 1e-14));
	assert(approxEqual(besselKD(0.0, 0.5), 3.210807086475177172355017867637452321623564307416114645664507, 0.0, 1e-14));
	assert(approxEqual(besselKD(1.0, 0.5), 2.543749846614157209231745521424755087195088446583780182048739, 0.0, 1e-14));
	assert(approxEqual(besselKD(1.0, 1.0), 1.888245647341297344338366806469818100073176922119421727713106, 0.0, 1e-14));
	assert(approxEqual(besselKD(1.0, 2.0), 1.477404884746695468820498957141015915645967793158828983150281, 0.0, 1e-14));
	assert(approxEqual(besselKD(0.4, 2.0), 1.502873827552600466534737445302356963995715395265386806359790, 0.0, 1e-14));
	assert(approxEqual(besselKD(0.4, 1.0), 2.015349430323277001160289680464116945013807398639704169522157, 0.0, 1e-14));
	assert(approxEqual(besselKD(0.4, 0.5), 3.071599175777718550825656123800662591905836277372864541044359, 0.0, 1e-13));
	assert(approxEqual(besselKD(11.1, 1.11e+10), 1.000000000090090090090090090045136443975807419926316735329998, 0.0, 1e-14));
	assert(approxEqual(besselKD(11.1, 1.11e-10), 1.099009900990099009900983462621096186432918126639008518272442, 0.0, 1e-14));
	assert(approxEqual(besselKD(2.0, 3.0), 1.296650971433360224728795947829976056477719397636727742896, 0.0, 1e-14));
	assert(approxEqual(besselKD(2.2, 0.8), 1.624595391806645609233071414803368115353081130470326304255669, 0.0, 1e-14));
	assert(approxEqual(besselKD(0.4, 3.0), 1.334349384832101337822320699631685431560403480505548447388101, 0.0, 1e-14));
	assert(approxEqual(besselKD(0.4, 0.8), 2.275575715445329347201695021084846610619371526255383391028726, 0.0, 1e-14));
	assert(approxEqual(besselKD(1.0, 1.0), 1.888245647341297344338366806469818100073176922119421727713106268163595832389436074069065151380128808, 0.0, 1e-14));
	assert(approxEqual(besselKD(1.0, 10.0), 1.099687904173359294358694759872440909233436542663077290833713959244496776538436835904856529604568407, 0.0, 1e-14));
	assert(approxEqual(besselKD(10.0, 1.0), 1.110348556190446111957216138204031855721452559868300859853180425692157659153049491662221677392988991, 0.0, 1e-14));
	assert(approxEqual(besselKD(10.0, 1.49011e-08), 1.11111, 0.0, 1e-4));
	assert(approxEqual(besselKD(1.0, 1.49011e-08), 36.2755, 0.0, 1e-4));
	assert(approxEqual(besselKD(1.0, 2.22044e-16), 72.3192, 0.0, 1e-4));
	assert(approxEqual(besselKD(1.0, 1.49166e-154), 708.628, 0.0, 1e-3));
	assert(approxEqual(besselKD(10.0, 1.49166e-154), 1.111111111111111, 0.0, 1e-4));
	assert(approxEqual(besselKD(1.1, double.min_normal), 11.0, 0.0, 1e-4));
	assert(isFinite(besselKD(0.1, double.epsilon)));
	assert(isFinite(besselKD(0, double.epsilon)));
}


T besselKRM(T)(in T nu, in T x)
	if(isFloatingPoint!T)
in {

}
body {
	if(nu < 0)
		return 1 / besselKRM(-nu, x);
	immutable anu = fabs(nu);
	int nl = cast(int)floor(anu+0.5f);
	T mu=anu-nl;

	if(x >= 2)
	{
		T r = modifiedBesselCF2(mu, x);
		mu *= 2;
		if(nl)
		{
			T d;
			do
			{
				d = r;
				mu += 2;
				r = mu / x + 1 / r;
			}
			while(--nl);
			return sqrt(r) * sqrt(d);
		}
		else
		{
			return sqrt(r / (r - mu / x));
		}
	}
	else
	{
		immutable sums = modifiedBesselTemmeSeriesImpl(mu, x);
		T r = sums[1] / sums[0];
		if(nl)
		{
			immutable x2_4 = x * x / 4;
			T d;
			do
			{
				d = r;
				mu += 1;
				r = mu + x2_4 / r;
			}
			while(--nl);
			return sqrt(r) * sqrt(d) * 2 / x;
		}
		else
		{
			return sqrt(r / (r - mu));
		}
	}
}

unittest {
	assert(approxEqual(besselKRM(0.0, 1.0),  1, 0.0, 1e-14));
	assert(approxEqual(besselKRM(0.0, 2.0), 1, 0.0, 1e-14));
	assert(approxEqual(besselKRM(0.0, 0.5), 1, 0.0, 1e-14));
	assert(approxEqual(besselKRM(1.0, 0.5), 2.857882088842869132027932428611809106473984497664375438969650, 0.0, 1e-14));
	assert(approxEqual(besselKRM(1.0, 1.0), 1.964497593920848638120029967522901757996791347548236796086390, 0.0, 1e-14));
	assert(approxEqual(besselKRM(1.0, 2.0), 1.492661023078886446121280337304042151124577136818374559084502, 0.0, 1e-14));
	assert(approxEqual(besselKRM(0.4, 2.0), 1.176363556186439775048077295541564293759652878052851708370510, 0.0, 1e-14));
	assert(approxEqual(besselKRM(0.4, 1.0), 1.320700844681798447461890099116285783758321865242462619205560, 0.0, 1e-14));
	assert(approxEqual(besselKRM(0.4, 0.5), 1.555719773577105697258998484956087203635813355802430411799783, 0.0, 1e-14));
	assert(approxEqual(besselKRM(11.1, 1.11e+10), 1.000000001000000000454954954912953493935881553614516153308770, 0.0, 1e-14));
	assert(approxEqual(besselKRM(11.1, 1.11e-10), 1.9077839604210010339594792869174419072661582323995475082e+11, 0.0, 1e-14));
}


T logBesselK(T)(in T nu, in T x)
	if(isFloatingPoint!T)
in {
	assert(x.isFinite);
	assert(x >= T.min_normal);
}
body {
	immutable anu = fabs(nu);
	int nl = cast(int)floor(anu+0.5f);
	T mu = anu - nl;
	if(x >= 2)
	{
		immutable r2 = modifiedBesselCF2Full(mu, x);
		T r = r2[1];
		T ret = log(r2[0]) - log(x) / 2 - x;
		if(nl)
		{
			mu *= 2;
			T l = 1;
			long exs;
			do
			{
				int ex = void;
				l *= frexp(r, ex);
				exs += ex;
				mu += 2;
				r = mu / x + 1 / r;
			}
			while(--nl);
			ret += cast(T)LN2 * (exs + log2(l));
		}
		return ret;
	}
	else
	{
		immutable sums = modifiedBesselTemmeSeriesImpl(mu, x);
		T r = sums[1] / sums[0];
		T ret = log(sums[0]);
		if(nl)
		{
			ret += nl * (cast(T)LN2 - log(x));
			immutable x2_4 = x * x / 4;
			T l = 1;
			long exs;
			do
			{
				int ex = void;
				l *= frexp(r, ex);
				exs += ex;
				mu += 1;
				r = mu + x2_4 / r;
			}
			while(--nl);
			ret += cast(T)LN2 * (exs + log2(l));
		}
		return ret;
	}
}

unittest {
	assert(approxEqual(logBesselK(0.0, 1.0),  -0.86506439890678809679875790803368568022489315035161867839839, 0.0, 1e-14));
	assert(approxEqual(logBesselK(0.0, 2.0), -2.17248820497570993473841333643717923143973636448830725037508, 0.0, 1e-14));
	assert(approxEqual(logBesselK(0.0, 0.5), -0.07858976986908141689523697453802973224104044591707028228710, 0.0, 1e-14));
	assert(approxEqual(logBesselK(1.0, 0.5), 0.504671397304651177308416839874408827783443947120148152746119, 0.0, 1e-14));
	assert(approxEqual(logBesselK(1.0, 1.0), -0.50765194821075233094791485120634298590875979858568077818450, 0.0, 1e-14));
	assert(approxEqual(logBesselK(1.0, 2.0), -1.96707130256051389147686464265533209478058062247688850653003, 0.0, 1e-14));
	assert(approxEqual(logBesselK(0.4, 2.0), -2.13936877477972262972321591624176676086808535490162854194273, 0.0, 1e-14));
	assert(approxEqual(logBesselK(0.4, 1.0), -0.80679541168661951923202239125477303379665660825595791745329, 0.0, 1e-14));
	assert(approxEqual(logBesselK(0.4, 0.5), 0.018456437585950072585426399714417863872587115447191398524780, 0.0, 1e-14));
	assert(approxEqual(logBesselK(11.1, 1.11e+10), -1.110000001133931411444888363310129265544563412394720435e10, 0.0, 1e-14));
	assert(approxEqual(logBesselK(11.1, 1.11e-10), 276.7693978383758755547170249282452519578463074349871817713218, 0.0, 1e-14));
}


T besselK(Flag!"ExponentiallyScaled" expFlag = Flag!"ExponentiallyScaled".no, T)(in T nu, in T x)
	if(isFloatingPoint!T)
in {
	assert(x.isFinite);
	assert(x >= T.min_normal);
}
body {
	immutable anu = fabs(nu);
	int nl = cast(int)floor(anu+0.5f);
	T mu = anu - nl;
	T r = void, ret = void;
	if(x >= 2)
	{
		immutable r2 = modifiedBesselCF2Full(mu, x);
		r = r2[1];
		static if(expFlag)
			ret = r2[0] / sqrt(x);
		else
			ret = r2[0] / sqrt(x) * exp(-x);
	}
	else
	{
		immutable sums = modifiedBesselTemmeSeriesImpl(mu, x);
		r = sums[1] / sums[0] / x * 2;
		static if(expFlag)
			ret = sums[0] * exp(x);
		else
			ret = sums[0];
	}
	if(nl)
	{
		mu *= 2;
		do
		{
			ret *= r;
			mu += 2;
			r = mu / x + 1 / r;
		}
		while(--nl);
	}
	return ret;

}

unittest {
	assert(approxEqual(besselK(0.0, 1.0),  0.421024438240708333335627379212609036136219748226660472298969, 0.0, 1e-14));
	assert(approxEqual(besselK(0.0, 2.0), 0.113893872749533435652719574932481832998326624388808882892529, 0.0, 1e-14));
	assert(approxEqual(besselK(0.0, 0.5), 0.924419071227665861781924167530216989538768311953529684815019, 0.0, 1e-14));
	assert(approxEqual(besselK(1.0, 0.5), 1.656441120003300893696445403174091511534100759464077446055427, 0.0, 1e-14));
	assert(approxEqual(besselK(1.0, 1.0), 0.601907230197234574737540001535617339261586889968106456017767, 0.0, 1e-14));
	assert(approxEqual(besselK(1.0, 2.0), 0.139865881816522427284598807035411023887234584841515530384442, 0.0, 1e-14));
	assert(approxEqual(besselK(0.4, 2.0), 0.117729133170423325690699234439848483526261145336868148793547, 0.0, 1e-14));
	assert(approxEqual(besselK(0.4, 1.0), 0.446285939834668179310023409179745077844437039396299727961180, 0.0, 1e-14));
	assert(approxEqual(besselK(0.4, 0.5), 1.018627810316608462220290727117637366662549451305627759999544, 0.0, 1e-14));
	assert(approxEqual(besselK(11.1, 100000.0), 0.0, 0.0, 1e-14));
	assert(approxEqual(besselK(11.1, 1/100000.0), 1.5940696949048155233993419471665743544335615887236248959e+65, 1e51, 1e-14));

	assert(approxEqual(besselK!(Flag!"ExponentiallyScaled".yes)(0.4, 2.0), 0.869907169474735192472996963053995791191543506954884340734219, 0.0, 1e-14));
	assert(approxEqual(besselK!(Flag!"ExponentiallyScaled".yes)(0.4, 1.0), 1.213130960549345270517674559709580335203213288805317162444177, 0.0, 1e-14));
	assert(approxEqual(besselK!(Flag!"ExponentiallyScaled".yes)(0.4, 0.5), 1.679433337795687807090050796759423855788542163267174953125487, 0.0, 1e-14));
	assert(approxEqual(besselK!(Flag!"ExponentiallyScaled".yes)(11.1, 100000.0), 0.003965764688216290045701642597190265433950245246337093558273, 0.0, 1e-14));

	assert(approxEqual(besselK(1.0, 2.0), 0.139865881816522427284, 0.0, 1e-14));
	assert(approxEqual(besselK(1.0, 20.0), 5.88305796955703817765028e-10, 0.0, 1e-14));
	assert(approxEqual(besselK(1.2, 1.0), 0.701080));
	assert(approxEqual(besselK(10.0, 400.0), 1.359308280937876101049800424530298700140591995078599604e-175, 0.0, 1e-14));
}


T besselKRS(T)(in T nu, in T x)
	if(isFloatingPoint!T)
in {
	assert(x.isFinite);
	assert(x >= T.min_normal);
}
body {
	immutable anu = fabs(nu);
	int nl = cast(int)floor(anu+0.5f);
	T mu=anu-nl;
	if(x >= 2)
	{
		T r = modifiedBesselCF2(mu, x);
		mu *= 2;
		if(nl)
		{
			T d = void;
			do
			{
				d = r;
				mu += 2;
				r = mu / x + 1 / r;
			}
			while(--nl);
			return r + 1 / d;
		}
		else
		{
			return r + (r - mu / x);
		}
	}
	else
	{
		immutable sums = modifiedBesselTemmeSeriesImpl(mu, x);
		T r = sums[1] / sums[0];
		immutable x_2 = x / 2;
		if(nl)
		{
			immutable x2_4 = x_2^^2;
			T d = void;
			do
			{
				d = r;
				mu += 1;
				r = mu + x2_4 / r;
			}
			while(--nl);
			return r / x_2 + x_2 / d;
		}
		else
		{
			return r / x_2 + (r - mu) / x_2;
		}
	}
}

unittest {
	assert(approxEqual(besselKRS(0.0, 1.0), 2.859250796520803516056216046900731260160034036385325966918068, 0.0, 1e-14));
	assert(approxEqual(besselKRS(0.0, 2.0), 2.456073859637815951485344904163437808478733905321719388934076, 0.0, 1e-14));
	assert(approxEqual(besselKRS(0.0, 0.5), 3.583745016864440467305949390266788044526362750672196202869812, 0.0, 1e-14));
	assert(approxEqual(besselKRS(1.0, 0.5), 5.116150836953170679741316342708831083392598534819729276449587, 0.0, 1e-14));
	assert(approxEqual(besselKRS(1.0, 1.0), 3.398967871187544687785354799542729509747041034806318935543498, 0.0, 1e-14));
	assert(approxEqual(besselKRS(1.0, 2.0), 2.628615517527578979897586948721927750995112315346619061302174, 0.0, 1e-14));
	assert(approxEqual(besselKRS(0.4, 2.0), 2.484249446052147531028619260526756105681146308342088119522102, 0.0, 1e-14));
	assert(approxEqual(besselKRS(0.4, 1.0), 2.949813167184170666766713209833464015196615131065026837642303, 0.0, 1e-14));
	assert(approxEqual(besselKRS(0.4, 0.5), 3.853102218097889263846807474862866521371874509444807669438352, 0.0, 1e-14));
	assert(approxEqual(besselKRS(11.1, 1.11e+10), 2.000000000090090091088061033917072660444297152286161364522101, 0.0, 1e-14));
	assert(approxEqual(besselKRS(11.1, 1.11e-10), 2.0000000000000000000001099009900990099009900953267052034e+11, 0.0, 1e-14));
}


T besselKRX(T)(in T nu, in T x)
	if(isFloatingPoint!T)
in {
	assert(nu >= 0);
	assert(x.isFinite);
	assert(x >= T.min_normal);
}
body {
	immutable anu = fabs(nu);
	int nl = cast(int)floor(anu+0.5f);
	T mu = anu - nl;
	if(x >= 2)
	{
		T r = modifiedBesselCF2(mu, x);
		if(nl)
		{
			mu *= 2;
			do
			{
				mu += 2;
				r = mu / x + 1 / r;
			}
			while(--nl);
		}
		return r * x;
	}
	else
	{
		immutable sums = modifiedBesselTemmeSeriesImpl(mu, x);
		T r = sums[1] / sums[0];
		if(nl)
		{
			immutable x2_4 = x * x / 4;
			do
			{
				mu += 1;
				r = mu + x2_4 / r;
			}
			while(--nl);
		}
		return r * 2;			
	}
}

unittest {
	assert(approxEqual(besselKRX(10.0, double.min_normal), 20.0, 0.0, 1e-14));
	assert(approxEqual(besselKRX(0.0, 1.0), 1.429625398260401758028108023450365630080017018192662983459034, 0.0, 1e-14));
	assert(approxEqual(besselKRX(0.0, 2.0), 2.456073859637815951485344904163437808478733905321719388934076, 0.0, 1e-14));
	assert(approxEqual(besselKRX(0.0, 0.5), 0.895936254216110116826487347566697011131590687668049050717453, 0.0, 1e-14));
	assert(approxEqual(besselKRX(1.0, 0.5), 2.279037709238292669935329085677207770848149633704932319112396, 0.0, 1e-14));
	assert(approxEqual(besselKRX(1.0, 1.0), 2.699483935593772343892677399771364754873520517403159467771749, 0.0, 1e-14));
	assert(approxEqual(besselKRX(1.0, 2.0), 3.628615517527578979897586948721927750995112315346619061302174, 0.0, 1e-14));
	assert(approxEqual(besselKRX(0.4, 2.0), 2.884249446052147531028619260526756105681146308342088119522102, 0.0, 1e-14));
	assert(approxEqual(besselKRX(0.4, 1.0), 1.874906583592085333383356604916732007598307565532513418821151, 0.0, 1e-14));
	assert(approxEqual(besselKRX(0.4, 0.5), 1.363275554524472315961701868715716630342968627361201917359588, 0.0, 1e-14));
	assert(approxEqual(besselKRX(11.1, 11.1), 27.05342065662651512004991741552637606461741839973809374665713, 0.0, 1e-14));
	assert(approxEqual(besselKRX(11, 2.0^^30), 1.000000010710209660444967811966826255158285934025032904639974*2.0^^30, 0.0, 1e-14));
	assert(approxEqual(besselKRX(11, 2.0^^30), 1.000000010710209660444967811966826255158285934025032904639974*2.0^^30, 0.0, 1e-14));
	assert(approxEqual(besselKRX(11, 2.0^^40), 1.000000000010459189070438615765505792843518644341506580233109*2.0^^40, 0.0, 1e-14));
	assert(approxEqual(besselKRX(11, 1e+10), 1.0000000011500000006037499999396249998267992188210915625e+10, 1, 1e-14));
	assert(approxEqual(besselKRX(11.1, 1.11e+10), 1.1100000011600000005538738738239753265465849195188195573e+10, 1, 1e-14));
	assert(approxEqual(besselKRX(11.1, 1.11e-10), 22.20000000000000000000060995049504950495049502906321387905301, 0.0, 1e-14));
}


T besselKR(T)(in T nu, in T x)
	if(isFloatingPoint!T)
in {
}
body {
	if(nu < 0)
		return 1 / besselKR(-nu-1, x);
	if(x.isNaN() || nu.isNaN() || x < 0)
		return T.init;
	if(x < T.min_normal)
		return T.infinity;
	immutable anu = fabs(nu);
	int nl = cast(int)floor(anu+0.5f);
	T mu = anu - nl;
	if(x >= 2)
	{
		T r = modifiedBesselCF2(mu, x);
		if(nl)
		{
			mu *= 2;
			do
			{
				mu += 2;
				r = mu / x + 1 / r;
			}
			while(--nl);
		}
		return r;
	}
	else
	{
		immutable sums = modifiedBesselTemmeSeriesImpl(mu, x);
		T r = sums[1] / sums[0];
		if(nl)
		{
			immutable x2_4 = x * x / 4;
			do
			{
				mu += 1;
				r = mu + x2_4 / r;
			}
			while(--nl);
		}
		return r * 2 / x;			
	}
}

unittest {
	assert(approxEqual(besselKR(0.0, 1.0)    , 1.429625398260401758028108023450365630080017018192662983459034, 0.0, 1e-14));
	assert(approxEqual(besselKR(0.0, 2.0) * 2, 2.456073859637815951485344904163437808478733905321719388934076, 0.0, 1e-14));
	assert(approxEqual(besselKR(0.0, 0.5) / 2, 0.895936254216110116826487347566697011131590687668049050717453, 0.0, 1e-14));
	assert(approxEqual(besselKR(1.0, 0.5) / 2, 2.279037709238292669935329085677207770848149633704932319112396, 0.0, 1e-14));
	assert(approxEqual(besselKR(1.0, 1.0)    , 2.699483935593772343892677399771364754873520517403159467771749, 0.0, 1e-14));
	assert(approxEqual(besselKR(1.0, 2.0) * 2, 3.628615517527578979897586948721927750995112315346619061302174, 0.0, 1e-14));
	assert(approxEqual(besselKR(0.4, 2.0) * 2, 2.884249446052147531028619260526756105681146308342088119522102, 0.0, 1e-14));
	assert(approxEqual(besselKR(0.4, 1.0)    , 1.874906583592085333383356604916732007598307565532513418821151, 0.0, 1e-14));
	assert(approxEqual(besselKR(0.4, 0.5) / 2, 1.363275554524472315961701868715716630342968627361201917359588, 0.0, 1e-14));
	assert(approxEqual(besselKR(11.1, 1.11e+10), 1.1100000011600000005538738738239753265465849195188195573e+10 / 1.11e+10, 1, 1e-14));
	assert(approxEqual(besselKR(11.1, 1.11e-10) * 1.11e-10, 22.20000000000000000000060995049504950495049502906321387905301, 0.0, 1e-14));

	assert(approxEqual(besselKR(-1.0, 1.0), 1 / 1.429625398260401758028108023450365630080017018192662983459034, 0.0, 1e-14));
}


T besselKR(T)(in T nu, in T x, uint k)
{
	T l = besselKR(nu, x);
	T r = l;
	T mu = nu * 2;
	while(k--)
	{
		mu += 2;
		r = mu / x + 1 / r;
		l *= r;
	}
	return l;
}

unittest {
	assert(approxEqual(besselKR(0.0, 1.0, 0)    , 1.429625398260401758028108023450365630080017018192662983459034, 0.0, 1e-14));
	assert(approxEqual(besselKR(0.0, 2.0, 0) * 2, 2.456073859637815951485344904163437808478733905321719388934076, 0.0, 1e-14));
	assert(approxEqual(besselKR(0.0, 1.0, 3), 105.0590223025824984495740493132204752844809530187891270737059, 1e-10, 1e-14));
	assert(approxEqual(besselKR(0.0, 2.0, 3), 19.28036929818907975742672452081718904239366952660859694467038, 1e-10, 1e-14));
	assert(approxEqual(besselKR(0.0, 0.5, 3), 813.7490033728880934611898780533576089052725501344392405739625, 1e-10, 1e-14));
	assert(approxEqual(besselKR(1.0, 0.5, 3), 7303.597652823473130198226747312888245046226857159388835630678, 1e-10, 1e-14));
	assert(approxEqual(besselKR(1.0, 1.0, 3), 599.6947228611295581541061895533584099941981855502445314254368, 1e-10, 1e-14));
	assert(approxEqual(besselKR(1.0, 2.0, 3), 67.42923276291368469846380423082891626492668473019928591953262, 1e-10, 1e-14));
	assert(approxEqual(besselKR(0.4, 2.0, 3), 32.55703150637502077170415944139304552928545885327101882556329, 1e-10, 1e-14));
	assert(approxEqual(besselKR(0.4, 1.0, 3), 222.9905656901318819890519502437505989113682776582595951935857, 1e-10, 1e-14));
	assert(approxEqual(besselKR(0.4, 0.5, 3), 2177.389452959348919338879066729351907090043415959389603727847, 1e-10, 1e-14));
}
