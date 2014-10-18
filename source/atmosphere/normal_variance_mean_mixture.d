/**
Algorithms for saparating normal variance mean mixture.
*/
module atmosphere.normal_variance_mean_mixture;

//phobos std.numeric.findRoot is not @nogc nothrow
//@nogc nothrow:


import core.stdc.stdlib;
import core.stdc.tgmath;

import simple_matrix;
import atmosphere.utilities;

import atmosphere.mixture;

/**
Algorithm for saparating normal variance mean mixture.
*/
enum SNVMMAlgorithm
{
	/**
	Expectation–maximization algorithm.
	*/
	ExpectationMaximization,
	/**
	Gradient descent optimization algorithm.
	*/	
	GradientDescent,
	/**
	Coordinate descent optimization algorithm.
	*/
	CoordinateDescent,
}

import std.range, std.traits, std.math : frexp;
/**
Computes accurate sum of binary logarithms of input range $(D r).
 */
ElementType!Range sumOfLog2s(Range)(Range r) 
    if (isInputRange!Range && isFloatingPoint!(ElementType!Range))
{
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

/**
Separates normal variance mean mixture.
*/
T separateNormalVarianceMeanMixture
	(
		SNVMMAlgorithm Algorithm, 
		T,
	)
	(
		in T[] sample,
		in T[] grid,
		T[] probability,
		in bool delegate(T alphaSave, T alpha, in T[] probabilitySave, in T[] probability) @nogc nothrow tolerance,
		in bool delegate(T, T) @nogc nothrow findRootTolerance = (a, b) => false,
	)
out(alpha)
{
	assert(isFinite(alpha));
	foreach(p; probability)
	{
		import std.conv;
		assert(p >= 0, p.to!string);
		assert(p <= 1);
	}
}
in
{
	assert(grid.length);
	assert(grid.length == probability.length);
	assert(sample.length >= grid.length);
	foreach(u; grid)
	{
		assert(u.isFinite);
		assert(u > 0);
	}
}
body
{
	T[] createArray(size_t length)
	{
		return (cast(T*)malloc(double.sizeof * length))[0..length];
	}
	const sampleAvg = sample.avg;

	T alpha = sampleAvg / dot(probability, grid);
	T alphaSave = void;

	auto ptr = createArray(grid.length * sample.length).ptr;
	auto WT = Matrix!(double)(ptr, grid.length, sample.length);
	auto probabilitySave = createArray(grid.length);
	scope(exit)	
	{
		WT.ptr.free;
		probabilitySave.ptr.free;
	}
	static if(Algorithm == SNVMMAlgorithm.GradientDescent || Algorithm == SNVMMAlgorithm.CoordinateDescent)
	{
		auto chi = createArray(sample.length);
		auto pi = createArray(sample.length);
		scope(exit)
		{
			chi.ptr.free;
			pi.ptr.free;
		}
	}
	static if(Algorithm == SNVMMAlgorithm.GradientDescent || Algorithm == SNVMMAlgorithm.ExpectationMaximization)
	{
		auto xi = createArray(sample.length);
		auto c = createArray(grid.length);
		scope(exit) 
		{
			xi.ptr.free;
			c.ptr.free;
		}
	}
	do
	{

		probabilitySave[] = probability[];
		alphaSave = alpha;
		foreach(i, u; grid)
		{
			const sqrtu = sqrt(u);
			auto column = WT[i];
			//up to a constant
			foreach(j, x; sample)
			{
				immutable y = (x - alpha * u) / sqrtu;
				column[j] = exp2(y*y/-2) / sqrtu;
				import std.conv;
				assert(column[j] > 0, column[j].to!string);
			}
		}

		import std.stdio;
		auto v = new T[sample.length];
		gemv(WT.transposed, probability, v);
		writefln("L = %s", v.sumOfLog2s);

		static if(Algorithm == SNVMMAlgorithm.ExpectationMaximization)
			simpleExpectationMaximizationIteration !(a => 1/a)(cast(Matrix!(const double))WT, probability, xi, c);
		static if(Algorithm == SNVMMAlgorithm.GradientDescent)
			simpleGradientDescentIteration  !(a => -1/a)(cast(Matrix!(const double))WT, probability, chi, pi, xi, c, findRootTolerance);
		static if(Algorithm == SNVMMAlgorithm.CoordinateDescent)
			simpleCoordinateDescentIteration!(a => -1/a)(cast(Matrix!(const double))WT, probability, chi, pi, findRootTolerance);
		alpha = sampleAvg / dot(probability, grid);


		writeln(alphaSave, ' ', alpha);
	}
	while(!tolerance(alphaSave, alpha, probabilitySave, probability));

	return alpha;
}


unittest
{
	alias SNVMMA_EM = separateNormalVarianceMeanMixture!(SNVMMAlgorithm.ExpectationMaximization, double);
	alias SNVMMA_GD = separateNormalVarianceMeanMixture!(SNVMMAlgorithm.GradientDescent, double);
	alias SNVMMA_CD = separateNormalVarianceMeanMixture!(SNVMMAlgorithm.CoordinateDescent, double);
}
