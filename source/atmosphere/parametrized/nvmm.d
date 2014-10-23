/**
*/
module atmosphere.parametrized.nvmm;

import atmosphere.mixture;
import atmosphere.internal;
import atmosphere.utilities : sumOfLog2s;
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

	override void update()
	{
		_log2Likelihood = mixture.sumOfLog2s;
		updateAlpha;
	}

	package T[] _sample;
	package const T[] _grid;

	package T _mean;
	package T _alpha;
	package T _log2Likelihood;

	/**
	
	*/
	this(in T[] _grid, size_t maxN)
	in
	{
		assert(_grid.length);
		assert(maxN);
	}
	body
	{
		super(_grid.length, maxN);
		this._grid = _grid.dup;
		this._sample = new T[maxN];
	}

final:

	/**
	
	*/
	void optimize(
			scope bool delegate (
				T alphaPrev, 
				T alpha, 
				T log2LikelihoodValuePrev, 
				T log2LikelihoodValue, 
				in T[] distributionPrev, 
				in T[] distribution,
			) 
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


	/**
	
	*/
	void optimize(
			scope bool delegate (
				T alphaPrev, 
				T alpha, 
				T log2LikelihoodValuePrev, 
				T log2LikelihoodValue, 
			) 
			tolerance,
			scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null,
		)
	{
		T log2LikelihoodPrev, alphaPrev;
		do
		{
			log2LikelihoodPrev = _log2Likelihood;
			alphaPrev = _alpha;
			eval(findRootTolerance);
		}
		while(!tolerance(alphaPrev, _alpha, log2LikelihoodPrev, _log2Likelihood));
	}

	/**
	
	*/
	void optimize(
			scope bool delegate (
				T alphaPrev, 
				T alpha, 
				in T[] distributionPrev, 
				in T[] distribution,
			) 
			tolerance,
			scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null,
		)
	{
		T alphaPrev;
		scope T[] distributionPrev = new T[distribution.length];
		do
		{
			alphaPrev = _alpha;
			assert(distribution.length == distributionPrev.length);
			distributionPrev[] = distribution[];
			eval(findRootTolerance);
		}
		while(!tolerance(alphaPrev, _alpha, distributionPrev, distribution));
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
		assert(_sample.length <= this._sample.length);
		foreach(s; _sample)
		{
			assert(std.math.isFinite(s));
		}
		assert(_componentsT.matrix.shift >= _sample.length);
	}
	body
	{
		reset;
		_componentsT.reserveBackN(_sample.length);
		this._sample[0.._sample.length] = _sample[];
		_mean = sample.sum/sample.length;
		updateAlpha;
		assert(_componentsT.matrix.width == _sample.length);
		updateComponents;
	}

	/**
	Returns:
		sample
	*/
	const(T)[] sample() @property const
	{
		return _sample[0..mixture.length];
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
	Returns:
		$(D alpha = mean / dotProduct(distribution, grid))
	*/
	T alpha() @property const
	{
		return _alpha;
	}


	/**
	Returns:
		grid
	*/
	const(T)[] grid() @property const
	{
		return _grid;
	}

	package void updateAlpha()
	in
	{
		assert(distribution.length == _grid.length);
	}
	body
	{
		_alpha =  _mean / dotProduct(distribution, _grid);
	}


	package void updateComponents()
	{
		auto m = _componentsT.matrix;
		assert(m.width == sample.length);
		version(atmosphere_gm_parallel)
		{
			import std.parallelism;
			//TODO: choice workUnitSize
			debug pragma(msg, "NormalVarianceMeanMixtureSeparator.updateComponents: parallel");
			auto pdfs = _grid.map!(u => PDF(alpha, u)).parallel;
			foreach(i, pdf; _grid.map!(u => PDF(alpha, u)).parallel)
			{
				auto r = m[i];
				foreach(x; sample)
				{
					r.front = pdf(x);
					r.popFront;
				}
			}
		}
		else
		{
			foreach(pdf; _grid.map!(u => PDF(alpha, u)))
			{
				auto r = m.front;
				m.popFront;
				foreach(x; sample)
				{
					r.front = pdf(x);
					r.popFront;
				}
			}
		}
		updateMixture;
	}

	private static struct PDF
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
}


/**
Expectation–maximization algorithm
*/
final class NormalVarianceMeanMixtureEMSeparator(T) : NormalVarianceMeanMixtureSeparator!T
{
	private T[] pi;
	private T[] c;

	///
	this(in T[] _grid, size_t maxN)
	in
	{
		assert(maxN);
	}
	body
	{
		super(_grid, maxN);
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
			(components, _distribution, mixture, pi[0..length], c);
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
	this(in T[] _grid, size_t maxN)
	{
		super(_grid, maxN);
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
			(components, _distribution, mixture, pi[0..length], xi[0..length], gamma[0..length], c, findRootTolerance is null ? (a, b) => false : findRootTolerance);
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
	this(in T[] _grid, size_t maxN)
	{
		super(_grid, maxN);
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
			(components, _distribution, _mixture[0..mixture.length], pi[0..length], findRootTolerance is null ? (a, b) => false : findRootTolerance);
		updateComponents;
	}
}


unittest {
	alias C0 = NormalVarianceMeanMixtureEMSeparator!(double);
	alias C1 = NormalVarianceMeanMixtureEMAndCoordinateSeparator!(double);
	alias C2 = NormalVarianceMeanMixtureEMAndGradientSeparator!(double);
}
