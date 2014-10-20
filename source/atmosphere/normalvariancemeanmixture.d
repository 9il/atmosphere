module atmosphere.normalvariancemeanmixture;

import atmosphere.stationary;
import atmosphere.internal;
import std.algorithm : sum;
import std.numeric : dotProduct;
import std.math : isFinite;
import core.stdc.tgmath : exp, sqrt;


/**
*/
abstract class NormalVarianceMeanMixtureSeparator(T) : StationaryOptimizer!T
{

private:

	T[] _sample;
	const T[] _grid;

	T _mean;
	T _alpha;
	T _log2Likelihood;

final:

	void updateAlpha()
	in
	{
		assert(_distribution.length == _grid.length);
	}
	body
	{
		_alpha =  _mean / dotProduct(_distribution, _grid);
	}

	void update()
	{
		updateMixture;
		_log2Likelihood = _mixture.sumOfLog2s;
		updateAlpha;
	}

	void updateComponents()
	{
		import std.algorithm : map;

		static struct Kernel
		{
			T alphau;
			T sqrtu;

			this(T alphau, T sqrtu)
			{
				this.alphau = alphau;
				this.sqrtu = sqrtu;
			}

			T opCall(T x) const
			{
				immutable y = (x - alphau) / sqrtu;
				return exp(y * y / -2) / sqrtu;
			}
		}
		auto kernels = _grid.map!(u => Kernel(alpha*u, sqrt(u)));
		components(kernels, _sample);
		update;
	}

public:

	/**
	
	*/
	this(in T[] _grid, in T[] _sample)
	in
	{
		assert(_grid.length);
		assert(_sample.length);
	}
	body
	{
		super(_grid.length, _sample.length);
		this._grid = _grid.dup;
		this._sample = new T[_sample.length];
		sample = _sample;
	}

	/**
	Returns:
		α
	*/
	T alpha() @property const
	{
		return _alpha;
	}

	/**
	Returns:
		sample mean
	*/
	T mean() @property const
	{
		return _mean;
	}

	/**
	Returns:
		LogLikelihood base 2.
	*/
	T log2Likelihood() @property const
	{
		return _log2Likelihood;
	}

	/**
	Sets sample and recalculate α.
	*/
	void sample(in T[] _sample) @property
	in
	{
		assert(_sample.length == this._sample.length);
		foreach(s; _sample)
		{
			assert(s.isFinite);
			assert(s > 0);
		}
	}
	body
	{
		this._sample[] = _sample[];
		_mean = sample.sum/sample.length;
		updateAlpha;
		updateComponents;
	}

	/**
	Returns:
		sample
	*/
	const(T)[] sample() @property const
	{
		return _sample;
	}

	/**
	Returns:
		grid
	*/
	const(T)[] grid() @property const
	{
		return _grid;
	}

	/**
	
	*/
	void optimize(
			scope bool delegate (
				T alphaPrev, 
				T alpha, 
				T log2LikelihoodValuePrev, 
				T log2LikelihoodValue, 
				in T[] distributionPrev, 
				in T[] distribution) 
			tolerance,
			scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null,
		)
	{
		T log2LikelihoodPrev, alphaPrev;
		scope T[] distributionPrev = new T[_distribution.length];
		do
		{
			log2LikelihoodPrev = _log2Likelihood;
			alphaPrev = _alpha;
			assert(distribution.length == distributionPrev.length);
			distributionPrev[] = _distribution[];
			eval(findRootTolerance);
		}
		while(!tolerance(alphaPrev, _alpha, log2LikelihoodPrev, _log2Likelihood, distributionPrev, _distribution));
	}
}


/**
Expectation–maximization algorithm
*/
final class NormalVarianceMeanMixtureEMSeparator(T) : NormalVarianceMeanMixtureSeparator!T
{
	private T[] pi;
	private T[] c;

	///
	this(in T[] _grid, in T[] _sample)
	{
		super(_grid, _sample);
		pi = new T[_sample.length];
		c = new T[_grid.length];
	}

	~this()
	{
		pi.destroy;
		c.destroy;
	}

	override void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null)
	{
		EMIteration!
			((a, b) {foreach(i, ai; a) b[i] = 1/ai;}, T)
			(componentsT, _distribution, _mixture, pi, c);
		updateComponents;
	}
}


/**
Expectation–maximization algorithm with inner gradient descend optimization.
*/
final class NormalVarianceMeanMixtureEMAndGradientSeparator(T) : NormalVarianceMeanMixtureSeparator!T
{
	private T[] pi;
	private T[] xi;
	private T[] gamma;
	private T[] c;

	///
	this(in T[] _grid, in T[] _sample)
	{
		super(_grid, _sample);
		pi = new T[_sample.length];
		xi = new T[_sample.length];
		gamma = new T[_sample.length];
		c = new T[_grid.length];
	}

	~this()
	{
		pi.destroy;
		xi.destroy;
		gamma.destroy;
		c.destroy;
	}

	override void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null)
	{
		gradientDescentIteration!
			((a, b) {foreach(i, ai; a) b[i] = -1/ai;}, T)
			(componentsT, _distribution, _mixture, pi, xi, gamma, c, findRootTolerance is null ? (a, b) => false : findRootTolerance);
		updateComponents;
	}
}


/**
Expectation–maximization algorithm with inner coordinate descend optimization.
Speed depends on permutation of elements of $(grid).
*/
final class NormalVarianceMeanMixtureEMAndCoordinateSeparator(T) : NormalVarianceMeanMixtureSeparator!T
{
	private T[] pi;

	///
	this(in T[] _grid, in T[] _sample)
	{
		super(_grid, _sample);
		pi = new T[_sample.length];
	}

	~this()
	{
		pi.destroy;
	}

	override void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null)
	{
		coordinateDescentIterationPartial!
			(a => -1/a, T)
			(componentsT, _distribution, _mixture, pi, findRootTolerance is null ? (a, b) => false : findRootTolerance);
		updateComponents;
	}
}


unittest {
	alias C0 = NormalVarianceMeanMixtureEMSeparator!(double);
	alias C1 = NormalVarianceMeanMixtureEMAndCoordinateSeparator!(double);
	alias C2 = NormalVarianceMeanMixtureEMAndGradientSeparator!(double);
}


private:

import std.traits : Unqual;
import std.math : log2, frexp;
/*
Computes accurate sum of binary logarithms of input range $(D r).
Will be avalible in std.numeric with with DMD 2.068.
 */
T sumOfLog2s(T)(T[] r) 
{
	import std.math : frexp; 

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
