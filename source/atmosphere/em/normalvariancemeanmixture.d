/**
*/
module atmosphere.em.normalvariancemeanmixture;

import atmosphere.mixtureoptimizer;
import atmosphere.internal;
import std.algorithm;
import std.range;
import std.numeric;
import std.traits;
import core.stdc.tgmath;

static import std.math;

/**
------
F(x) = ∫_0^∞ Φ((x-αu_i)√u) dG(u) ≈ Σ_i p_i*Φ((x-αu_i)/sqrt(u))
α - alpha
Φ - standard normal distribution
G - mixture distribution
p - approximation of G
------
*/
abstract class NormalVarianceMeanMixtureSeparator(T) : MixtureOptimizer!T
{
private:

	static struct Kernel
	{
		T alphau;
		T sqrtu;

		this(T alpha, T u)
		{
			assert(u > 0);
			this.alphau = alpha*u;
			this.sqrtu = sqrt(u);
			assert(sqrtu > 0);
		}

		T opCall(T x) const
		{
			immutable y = (x - alphau) / sqrtu;
			return exp(y * y / -2) / sqrtu;
		}
	}

	T[] _sample;
	const T[] _grid;
	Kernel[] kernels;

	T _mean;
	T _alpha;
	T _log2Likelihood;

final:

	void updateAlpha()
	in
	{
		assert(distribution.length == _grid.length);
	}
	body
	{
		_alpha =  _mean / dotProduct(distribution, _grid);
	}

	void updateMixture()
	{
		super.updateMixture;
		_log2Likelihood = mixture.sumOfLog2s;
		updateAlpha;
	}

	void updateComponents()
	{
		reset;
		_grid.map!(u => Kernel(alpha, u)).copy(kernels);
		put(_sample.map!(a => kernels.map!(k => k(a))));
		updateMixture;
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
		this.kernels = new Kernel[_grid.length];
		sample(_sample);
	}

	~this()
	{
		kernels.destroy;
	}

	/**
	Returns:
		$(D alpha = mean / dotProduct(distribution, grid))
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
	Sets sample and recalculate alpha and mixture.
	-------
	///use the same length
	double[] newSample = new double[optimizer.sample.length];
	///init newSample
	...

	optimizer.sample = newSample;
	auto newAlpha = optimizer.alpha;
	///use dup to save mixture
	auto newMixture = optimizer.mixture.dup;
	///use dup to save distribution, the distribution is not changed
	auto distrubution = optimizer.distribution.dup;

	///performe one iteration
	optimizer.eval;

	///check new data
	auto secondAlpha = optimizer.alpha;
	auto secondMixture = optimizer.mixture.dup;
	///The distribution has changed after $(D eval).
	auto newDistribution = optimizer.distribution.dup;
	-------
	Params:
		_sample = new sample with the same length
	*/
	void sample(in T[] _sample) @property
	in
	{
		assert(_sample.length == this._sample.length);
		foreach(s; _sample)
		{
			assert(std.math.isFinite(s));
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
		scope T[] distributionPrev = new T[distribution.length];
		do
		{
			log2LikelihoodPrev = _log2Likelihood;
			alphaPrev = _alpha;
			assert(distribution.length == distributionPrev.length);
			distributionPrev[] = distribution[];
			eval(findRootTolerance);
		}
		while(!tolerance(alphaPrev, _alpha, log2LikelihoodPrev, _log2Likelihood, distributionPrev, distribution));
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
			(cast(Matrix!(const T))_componentsT.matrix, _distribution, _mixture, pi, c);
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
			(cast(Matrix!(const T))_componentsT.matrix, _distribution, _mixture, pi, xi, gamma, c, findRootTolerance is null ? (a, b) => false : findRootTolerance);
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
			(cast(Matrix!(const T))_componentsT.matrix, _distribution, _mixture, pi, findRootTolerance is null ? (a, b) => false : findRootTolerance);
		updateComponents;
	}
}


unittest {
	alias C0 = NormalVarianceMeanMixtureEMSeparator!(double);
	alias C1 = NormalVarianceMeanMixtureEMAndCoordinateSeparator!(double);
	alias C2 = NormalVarianceMeanMixtureEMAndGradientSeparator!(double);
}


private:

/*
Computes accurate sum of binary logarithms of input range $(D r).
Will be avalible in std.numeric with with DMD 2.068.
 */
T sumOfLog2s(T)(T[] r) 
{
	import std.math : frexp; 
	import std.traits : Unqual;

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
