/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.math;

//import core.stdc.tgmath;
import std.range;
import std.traits;

version(LDC)
{
	import ldc.intrinsics: log = llvm_log, log2 = llvm_log2, fabs = llvm_fabs;
}
else
{
	import std.math: log, log2, fabs;
}

//version (LDC)
//{
//	pragma(LDC_intrinsic, "llvm.fmuladd.f#")
//	    T llvm_fmuladd(T)(T vala, T valb, T valc) @safe pure nothrow @nogc;
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
 *      Handbook of mathematical functions. Dover, New York,
 *      pages 258-259, equation 6.3.18.
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

version(none):


T besselKR(T)(T x, T lambda)
{
	import bessel;
	immutable b = besselK(x, lambda+1, Flag!"ExponentiallyScaled".yes);
	immutable c = besselK(x, lambda  , Flag!"ExponentiallyScaled".yes);
	return b / c;
}


double besselKD(double x, double lambda)
{
	alias T = double;
	import std.math : sqrt;
	//enum sqeps = sqrt(T.epsilon);
	// x -> inf
	if(1 / (x*x) <= T.epsilon)
	{
		return 1 + 1 / x;
	}
	else if(x*x <= T.min_normal)
	{
		immutable alambda = fabs(lambda);
		if(alambda == 0)
		{
			return (x * log(x))^^-2;
		}
		else if(alambda < 1)
		{
			immutable onemalambda = 1 - alambda;
			return 
			 	tgamma(onemalambda)
				* x_gamma(alambda)
			 	* pow(T(4), onemalambda)
			 	* pow(x, -2*onemalambda);
		}
		else if(alambda == 1)
		{
			return 2 * log(1 / x);
			//import bessel;
			//immutable a = besselK(x, lambda-1, Flag!"ExponentiallyScaled".yes);
			//immutable b = besselK(x, lambda+1, Flag!"ExponentiallyScaled".yes);
			//immutable c = besselK(x, lambda  , Flag!"ExponentiallyScaled".yes);
			//return (a / c) * (b / c);

		}
		else if(alambda < 2)
		{
			immutable alambdamone = alambda - 1;
			immutable twomalambda = 2 - alambda;
			return  alambdamone / alambdamone * (1
			 	- tgamma(twomalambda)
				/ tgamma(alambda)
			 	* pow(T(4), -alambdamone)
			 	* pow(x, 2*alambdamone)
				);
		}
		else if(alambda == 2)
		{
			//return 2 * (1 + x * x * (T(3)/8 + log(x / 2) / 2));
			return 2 * (1 + x * (x * (T(3)/8 + log(x / 2) / 2)));
		}
		immutable alambdamone = alambda - 1;
		immutable alambdamtwo = alambda - 2;
		return alambda / alambdamone * (1
			- x
			/ (2 * alambda * alambdamone * alambdamtwo)
			* x
			);
	}
	import bessel;
	immutable a = besselK(x, lambda-1, Flag!"ExponentiallyScaled".yes);
	immutable b = besselK(x, lambda+1, Flag!"ExponentiallyScaled".yes);
	immutable c = besselK(x, lambda  , Flag!"ExponentiallyScaled".yes);
	return (a / c) * (b / c);
}


unittest {
	import std.mathspecial;
	import std.stdio;
	import bessel;
	double ll = besselK((real.epsilon), 10, Flag!"ExponentiallyScaled".yes);
	//import std.conv;
	//bool b = cast(bool)ll.isFinite;
	//if(b)
	//	writeln("ssss");
	//else
	//	writeln("rrrr");

	//writeln(1000);
	//writefln("%s", ll.isFinite);
	//writeln(gamma(double.epsilon));
	//writeln(tgamma(double.epsilon));
	//writeln(sqrt(double.epsilon));
	//writeln(1/sqrt(double.epsilon));
	double[3][] data = [
		[1, 1, 1.888245647341297344338366806469818100073176922119421727713106268163595832389436074069065151380128808],
		[10, 1, 1.099687904173359294358694759872440909233436542663077290833713959244496776538436835904856529604568407],
		[1, 10, 1.110348556190446111957216138204031855721452559868300859853180425692157659153049491662221677392988991],
		[1.49011e-08, 10, 1.11111],
		[1.49011e-08, 1, 36.2755],
		[2.22044e-16, 1, 72.3192],
		[1.49166e-154, 1, 708.628],
		[1.49166e-154, 10, 1.111111111111111],
		[0, 10, 10.0/9],

	];
	foreach(test; data)
	{
		writeln(besselKD(test[0], test[1]), " ", test[2]);
	}
	writeln(double.epsilon);
	writeln(sqrt(double.epsilon));
	writeln(double.min_normal);
	writeln(sqrt(double.min_normal));
	assert(1.49166e-154*1.49166e-154 < double.min_normal);

	writeln(tgamma(-PI));
	writeln(tgamma(-PI/2));
	writeln(gamma(-PI));
	writeln(gamma(-PI/2));
}

T besselKDDerivative(T)(T x, T lambda)
{
	import bessel;
	immutable a = besselK(x, lambda-1, Flag!"ExponentiallyScaled".yes);
	immutable b = besselK(x, lambda+1, Flag!"ExponentiallyScaled".yes);
	immutable c = besselK(x, lambda, Flag!"ExponentiallyScaled".yes);
	return (a / c) * (b / c);
}

immutable real EULERGAMMA = 0.57721_56649_01532_86060_65120_90082_40243_10421_59335_93992L; /** Euler-Mascheroni constant 0.57721566.. */
private immutable real X_GAMMA_A = -0.65587807152025388107701951514539048127976638047858434729236L; //1/12 (6 gamma^2-pi^2)
private immutable real X_GAMMA_B = -0.04200263503409523552900393487542981871139450040110609352206L; //1/12 (2 eulergamma^3-eulergamma pi^2-2 polygamma(2, 1))

T x_gamma(T)(T x)
{
	if(x < T.epsilon)
	{
		return x^^2 * (1 + x * (cast(T) EULERGAMMA + x * (cast(T) X_GAMMA_A + x * cast(T) X_GAMMA_B)));
	}
	else
	{
		return x / tgamma(x);
	}
}