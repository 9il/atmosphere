/++
Module contains classes that perform optimization of mixture weights over sliding window.
------
problem: p' = argmin f(p), p_i >= 0, Σ_i p_i = 1.

p - mixture weights,
f = u(Wp),
u(ω) - convex function,
W - matrix of features(n rows, k columns),
k - length of mixture weights,
n - length of sample, n may vary (sliding window).
------

The same for likelihood maximization:
------
problem: p' = argmax L(p), p_i >= 0, Σ_i p_i = 1.

L(ω) = u(Wp)
ω = Wp,
u(ω) = Σ_j log(ω_j)  (LogLikelihood)

p - mixture weights,
f_i - mixture components
W - matrix of posterior probabilities(n rows, k columns),
W[j, i] = f_i(x_j), 1 <= i <= k, 1 <= j <= n.

k - length of mixture weights (count of ),
n - length of sample, n may vary (sliding window).
------

+/
/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.mixture;

import core.stdc.tgmath;

import atmosphere.internal;
import std.range;
import std.traits;
import std.numeric : dotProduct;
import std.algorithm;
import std.typecons;

import atmosphere.internal;

import std.math : isNormal, LN2;

/++
Exception thrown for MixtureOptimizer.
+/
class MixtureOptimizerException : Exception
{
	/++
	Constructor which takes an error message.
	Params:
		msg  = Message describing the error.
		file = The file where the error occurred.
		line = The line where the error occurred.
	+/
	this(string msg, 
		string file = __FILE__, 
		size_t line = __LINE__) 
		@safe pure
	{
		super(msg, file, line);
	}
}

///
class FeaturesException : MixtureOptimizerException
{
	/++
	Constructor which takes an error message.
	Params:
		msg  = Message describing the error.
		file = The file where the error occurred.
		line = The line where the error occurred.
	+/
	this(string msg = "There is value in the sample that incompatible with all PDFs or machine precision is insufficient.", 
		string file = __FILE__, 
		size_t line = __LINE__) 
		@safe pure
	{
		super(msg, file, line);
	}
}


/++
+/
abstract class MixtureOptimizer(T)
	if(isFloatingPoint!T)
{
	import std.datetime;

	package SlidingWindow!T _featuresT;
	package T[] _weights;
	package T[] _mixture;

	/++
	Constructor
	Params:
		k = number of components
		maxLength = maximal length of features. In terms of likelihood maximization maxLength is maximal length of a sample.
	+/
	this(size_t k, size_t maxLength)
	in
	{
		assert(k);
		assert(maxLength);
	}
	body
	{
		_featuresT = SlidingWindow!T(k, maxLength);
		_weights = new T[k];
		_weights[] = T(1)/k;
		_mixture = new T[maxLength];
	}

	/++
	Perform k (1) iterations of coordinate (gradient or EM) descent optimization algorithm.
	Params:
		findRootTolerance = Defines an early termination condition. 
			Receives the current upper and lower bounds on the root. 
			The delegate must return true when these bounds are acceptable.
	See_Also: $(STDREF numeric, findRoot)
		evaluate
		optimize
	+/
	abstract void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null);

	/++
	update method is called when mixture changes occur.
	+/
	abstract void update();

final:

	/++
	Perform `itersCount` iterations of `eval`.
	Params:
		itersCount = count of iterations
		findRootTolerance = Defines an early termination condition. 
			Receives the current upper and lower bounds on the root. 
			The delegate must return true when these bounds are acceptable.
	See_Also: $(STDREF numeric, findRoot)
		eval
		optimize
	+/
	void evaluate(size_t itersCount, scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null)
	{
		foreach(_; 0..itersCount)
			eval(findRootTolerance);
	}

	/++
	Perform `eval` during `time`.
	Params:
		time = time duration
		findRootTolerance = Defines an early termination condition. 
			Receives the current upper and lower bounds on the root. 
			The delegate must return true when these bounds are acceptable.
	See_Also: $(STDREF numeric, findRoot)
		eval
		optimize
	+/
	Tuple!(TickDuration, "duration", size_t, "itersCount") 
	evaluate(TickDuration time, scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null)
	{
		TickDuration dur;
		size_t iterCount;
		auto sw = StopWatch(AutoStart.yes);
		do
		{
			eval(findRootTolerance);	
			iterCount++;
			dur = sw.peek;	
		}
		while(dur < time);
		return typeof(return)(dur, iterCount);
	}

	/++
	Performs optimization.
	Params:
		objectiveFunction = accepts mixture.
		tolerance = Defines an early termination condition. 
			Receives the current and previous versions of `objectiveFunction(mixture))` and weights. 
			The delegate must return true when mixture and weights are acceptable. 
		findRootTolerance = Tolerance for inner optimization.
	See_Also: $(STDREF numeric, findRoot)
		eval
		evaluate
	+/
	void optimize(
			scope T delegate(in T[] mixture) objectiveFunction, 
			scope bool delegate (T objectiveFunctionValuePrev, T objectiveFunctionValue, in T[] weightsPrev, in T[] weights) tolerance,
			scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null,
		)
	{
		assert(length);
		T objectiveFunctionValuePrev;
		T objectiveFunctionValue = objectiveFunction(mixture);
		T[] weightsPrev = new T[weights.length];
		do
		{
			objectiveFunctionValuePrev = objectiveFunctionValue;
			weightsPrev[] = weights[];
			eval(findRootTolerance);
			objectiveFunctionValue = objectiveFunction(mixture);
		}
		while(!tolerance(objectiveFunctionValuePrev, objectiveFunctionValue, weightsPrev, weights));
	}

	/++
	Performs optimization.
	Params:
		objectiveFunction = accepts mixture.
		tolerance = Defines an early termination condition. 
			Receives the current and previous versions of `objectiveFunction(mixture))`. 
			The delegate must return true when mixture are acceptable. 
		findRootTolerance = Tolerance for inner optimization.
	See_Also: $(STDREF numeric, findRoot)
		eval
		evaluate
	+/
	void optimize
	(
		scope T delegate(in T[] mixture) objectiveFunction, 
		scope bool delegate (T objectiveFunctionValuePrev, T objectiveFunctionValue) tolerance,
		scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null,
	)
	{
		assert(length);
		T objectiveFunctionValuePrev;
		T objectiveFunctionValue = objectiveFunction(mixture);
		do
		{
			objectiveFunctionValuePrev = objectiveFunctionValue;
			eval(findRootTolerance);
			objectiveFunctionValue = objectiveFunction(mixture);
		}
		while(!tolerance(objectiveFunctionValuePrev, objectiveFunctionValue));
	}

	/++
	Performs optimization.
	Params:
		tolerance = Defines an early termination condition. 
			Receives the current and previous versions of weights. 
			The delegate must return true when mixture and weights are acceptable. 
		findRootTolerance = Tolerance for inner optimization.
	See_Also: $(STDREF numeric, findRoot)
		eval
		evaluate
	+/
	void optimize
	(
		scope bool delegate (in T[] weightsPrev, in T[] weights) tolerance,
		scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null,
	)
	{
		assert(length);
		T[] weightsPrev = new T[weights.length];
		do
		{
			weightsPrev[] = weights[];
			eval(findRootTolerance);
		}
		while(!tolerance(weightsPrev, weights));
	}

	/++
	Puts back new feature for each components.
	Params:
		features = One feature per component.
	+/
	void put(Range)(Range features)
	if (isInputRange!Range && hasLength!Range && isNumeric!(ElementType!Range))
	in
	{
		assert(_featuresT.matrix.length == featuresROR.length);
		assert(_featuresT.matrix.width < _featuresT.matrix.shift);
	}
	body
	{
		_featuresT.put(features);
		updateMixtureBack;
	}

	/++
	Puts back new features for each components.
	Params:
		featuresROR = Range of ranges of features per component.
	+/
	void put(RangeOfRanges)(RangeOfRanges featuresROR)
	if (isInputRange!RangeOfRanges && hasLength!RangeOfRanges && 
		isInputRange!(ElementType!RangeOfRanges) && hasLength!(ElementType!RangeOfRanges) && 
		isNumeric!(ElementType!(ElementType!RangeOfRanges)))
	in
	{
		assert(_featuresT.matrix.length == featuresROR.front.length);
		assert(_featuresT.matrix.width + featuresROR.length <= _featuresT.matrix.shift);
	}
	body
	{
		const n = featuresROR.length;
		.put(_featuresT, featuresROR);
		updateMixtureBackN(n);
	}

	/++
	Returns: current length of features for each component
	+/
	size_t length() @property const
	{
		return _featuresT.matrix.width;
	}

	/++
	Returns: maximal allowed length of features
	+/
	size_t maxLength() @property const
	{
		return _featuresT.matrix.shift;
	}


	/++
	Reset length of features to zero.
	+/
	void reset()
	{
		_featuresT.reset;
	}

	/++
	Remove one front feature for each component.
	+/
	void popFront()
	in
	{
		assert(length);
	}
	body
	{
		_featuresT.popFront;
		_mixture[0.._featuresT.length] = _mixture[1.._featuresT.length+1];
	}

	/++
	Remove n front features for each component.
	Params:
		n = features will be removed
	+/
	void popFrontN(size_t n)
	in
	{
		assert(length >= n);
	}
	body
	{
		_featuresT.popFrontN(n);
		_mixture[0.._featuresT.length] = _mixture[n.._featuresT.length+n];
	}

	/++
	Returns: Range of range of features for each component. A matrix with k rows and n columns.
		This is internal representation, and can be discarded after any methods calls.
	+/
	Matrix!(const(T)) features() const
	{
		return cast(typeof(return))_featuresT.matrix;
	}

	/++
	Returns: Const slice of the internal mixture representation.
	Example:
	---
	double objectiveFunction(in double[])
	{
		...
	}

	//save slice
	auto mixture = optimizer.mixture;

	auto value0 = objectiveFunction(mixture);
	optimizer.eval;
	auto value1 = objectiveFunction(mixture);


	//use `.dup` or copy to save current mixture

	//1: .dup
	auto mixtureSave1 = mixture.dup;

	//2: create array
	auto mixtureSave2 = new double[mixture.length];
	//2: copy
	mixtureSave2[] = mixture[];
	---
	+/
	const(T)[] mixture() @property const
	{
		return _mixture[0.._featuresT.length];
	}

	/++
	Returns: Const slice of the internal weights representation.
	Example:
	---
	//save slice
	auto weights = optimizer.weights;

	//use `.dup` or copy to save current weights

	//1: .dup
	auto weightsSave1 = weights.dup;

	//2: create array
	auto weightsSave2 = new double[weights.length];
	//2: copy
	weightsSave2[] = weights[];
	---
	+/
	const(T)[] weights() @property const
	{
		return _weights;
	}

	/++
	Set the mixture weights and calls [update](update).
	Params:
		_weights = new mixture weights
	+/
	void weights(Range)(Range _weights) @property
		if(isInputRange!Range) 
	{
		_weights.copy(this._weights);
		updateMixture;
	}

	package void updateMixture()
	{
		if (length)
		{
			mix(cast(Matrix!(const T))_featuresT.matrix, _weights, _mixture[0.._featuresT.matrix.width]);
			update();			
		}
	}

	package void updateMixtureBack()
	in
	{
		assert(length);
	}
	body
	{
		_mixture[_featuresT.matrix.width-1] = dotProduct(_weights, _featuresT.back);
		update();
	}

	package void updateMixtureBackN(size_t n)
	in
	{
		assert(length >= n);
	}
	body
	{
		mix(cast(Matrix!(const T))_featuresT[$-n..$].matrix, _weights, _mixture[0.._featuresT.matrix.width]);
		update();
	}

	package void updateMixturePopBack()
	in
	{
		assert(length);
	}
	body
	{
		updateMixturePopBackN(1);
	}

	package void updateMixturePopBackN(size_t n)
	in
	{
		assert(length >= n);
	}
	body
	{
		_mixture[0.._featuresT.matrix.width-n] = _mixture[n.._featuresT.matrix.width];
		update();
	}
}


/++
Params:
	Gradient = Gradient of the objective function. `Gradient(a, b)` should perform `b = grad_f(a)`.
+/
class ExpectationMaximization(alias Gradient, T) : MixtureOptimizer!T
	if(isFloatingPoint!T)
{
	private T[] pi;
	private T[] c;

	/++
	Constructor
	Params:
		k = number of components
		maxLength = maximal length of features. In terms of likelihood maximization maxLength is maximal length of a sample.
	+/
	this(size_t k, size_t maxLength)
	{
		super(k, maxLength);
		pi = new T[maxLength];
		c = new T[k];
	}

	final override void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null)
	{
		EMIteration!(Gradient, T)(cast(Matrix!(const T))_featuresT.matrix, _weights, mixture, pi[0.._featuresT.matrix.width], c);
		updateMixture;
	}

	override void update(){};
}


/++
Params:
	Gradient = Gradient of the objective function. `Gradient(a, b)` should perform `b = grad_f(a)`.
+/
class GradientDescent(alias Gradient, T) : MixtureOptimizer!T
	if(isFloatingPoint!T)
{
	private T[] pi;
	private T[] xi;
	private T[] gamma;
	private T[] c;

	/++
	Constructor
	Params:
		k = number of components
		maxLength = maximal length of features. In terms of likelihood maximization maxLength is maximal length of a sample.
	+/
	this(size_t k, size_t maxLength)
	{
		super(k, maxLength);
		pi = new T[maxLength];
		xi = new T[maxLength];
		gamma = new T[maxLength];
		c = new T[k];
	}

	final override void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null)
	{
		gradientDescentIteration!(Gradient, T)(cast(Matrix!(const T))_featuresT.matrix, _weights, mixture, pi[0.._featuresT.matrix.width], xi[0.._featuresT.matrix.width], gamma[0.._featuresT.matrix.width], c, findRootTolerance is null ? (a, b) => false : findRootTolerance);
		updateMixture;
	}

	override void update(){};
}


/++
Params:
	PartialDerivative = Partial derivative `y` of objective convex function `u`: `du/dω_j = y(ω_j), 1 <= j <= n`.
+/
class GradientDescentPartial(alias PartialDerivative, T) : MixtureOptimizer!T
	if(isFloatingPoint!T)
{
	private T[] pi;
	private T[] xi;
	private T[] c;

	/++
	Constructor
	Params:
		k = number of components
		maxLength = maximal length of features. In terms of likelihood maximization maxLength is maximal length of a sample.
	+/
	this(size_t k, size_t maxLength)
	{
		super(k, maxLength);
		pi = new T[maxLength];
		xi = new T[maxLength];
		c = new T[k];
	}

	final override void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null)
	{
		gradientDescentIterationPartial!(PartialDerivative, T)(cast(Matrix!(const T))_featuresT.matrix, _weights, _mixture[0.._featuresT.matrix.width], pi[0.._featuresT.matrix.width], xi[0.._featuresT.matrix.width], c, findRootTolerance is null ? (a, b) => false : findRootTolerance);
		updateMixture;
	}

	override void update(){};
}


/++
Params:
	Gradient = Gradient of the objective function. `Gradient(a, b)` should perform `b = grad_f(a)`.
+/
class CoordinateDescent(alias Gradient, T) : MixtureOptimizer!T
	if(isFloatingPoint!T)
{
	private T[] pi;
	private T[] xi;
	private T[] gamma;

	/++
	Constructor
	Params:
		k = number of components
		maxLength = maximal length of features. In terms of likelihood maximization maxLength is maximal length of a sample.
	+/
	this(size_t k, size_t maxLength)
	{
		super(k, maxLength);
		pi = new T[maxLength];
		xi = new T[maxLength];
		gamma = new T[maxLength];
	}

	final override void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null)
	{
		immutable n = _featuresT.matrix.width;
		coordinateDescentIteration!(Gradient, T)(cast(Matrix!(const T))_featuresT.matrix, _weights, _mixture[0..n], pi[0..n], xi[0..n], gamma[0..n], findRootTolerance is null ? (a, b) => false : findRootTolerance);
		updateMixture;
	}

	override void update(){};
}


/++
Params:
	PartialDerivative = Partial derivative `y` of objective convex function `u`: `du/dω_j = y(ω_j), 1 <= j <= n`.
+/
class CoordinateDescentPartial(alias PartialDerivative, T) : MixtureOptimizer!T
	if(isFloatingPoint!T)
{
	private T[] pi;

	/++
	Constructor
	Params:
		k = number of components
		maxLength = maximal length of features. In terms of likelihood maximization maxLength is maximal length of a sample.
	+/
	this(size_t k, size_t maxLength)
	{
		super(k, maxLength);
		pi = new T[maxLength];
	}

	final override void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null)
	{
		coordinateDescentIterationPartial!(PartialDerivative, T)(cast(Matrix!(const T))_featuresT.matrix, _weights, _mixture[0.._featuresT.matrix.width], pi[0.._featuresT.matrix.width], findRootTolerance is null ? (a, b) => false : findRootTolerance);
		updateMixture;
	}

	override void update(){};
}


/++
+/
interface LikelihoodAscent(T)
	if(isFloatingPoint!T)
{
	/++
	Params:
		pdfs = Probability distribution *functions*
		sample = Sample
	See_Also: $(STDREF traits, isCallable)
	+/
	void put(PDFRange, SampleRange)(PDFRange pdfs, SampleRange sample)
		if (isInputRange!PDFRange && hasLength!PDFRange && isCallable!(ElementType!PDFRange));

	/++
	Performs optimization.
	Params:
		tolerance = Defines an early termination condition. 
			Receives the current and previous versions of normalized log-likelihood and weights. 
			The delegate must return true when likelihood and weights are acceptable. 
		findRootTolerance = Tolerance for inner optimization.
	See_Also: $(STDREF numeric, findRoot)
	+/
	void optimize
	(
		scope bool delegate (T likelihoodValuePrev, T likelihoodValue, in T[] weightsPrev, in T[] weights) tolerance,
		scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null,
	);

	/++
	Performs optimization.
	Params:
		tolerance = Defines an early termination condition. 
			Receives the current and previous versions of normalized log-likelihood. 
			The delegate must return true when likelihood are acceptable. 
		findRootTolerance = Tolerance for inner optimization.
	Throws: [FeaturesException](atmosphere/mixture/FeaturesException.html) if [isFeaturesCorrect](atmosphere/mixture/LikelihoodAscent.isFeaturesCorrect.html) is false.
	See_Also: $(STDREF numeric, findRoot)
	+/
	void optimize
	(
		scope bool delegate (T likelihoodValuePrev, T likelihoodValue) tolerance,
		scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null,
	);

	/++
	Returns true if for all values in sample exists probability density function from mixture such that `isNormal(PDF(value) * (1.0 / sample.length))`
	See_Also: $(STDREF math, isNormal)
	+/
	bool isFeaturesCorrect() const;


	/++
	Returns: normalized log-likelihood.
	+/
	T likelihood() @property const;
}


/++
+/
class LikelihoodAscentCoordinate(T) : CoordinateDescentPartial!("-1/a", T), LikelihoodAscent!T
	if(isFloatingPoint!T)
{
	/++
	Constructor
	Params:
		k = number of components
		maxLength = maximal length of features. In terms of likelihood maximization maxLength is maximal length of a sample.
	+/
	this(size_t k, size_t maxLength)
	{
		super(k, maxLength);
	}

	mixin LikelihoodAscentTemplate!T;
}


/++
+/
class LikelihoodAscentGradient(T) : GradientDescentPartial!("-1/a", T), LikelihoodAscent!T
	if(isFloatingPoint!T)
{
	/++
	Constructor
	Params:
		k = number of components
		maxLength = maximal length of features. In terms of likelihood maximization maxLength is maximal length of a sample.
	+/
	this(size_t k, size_t maxLength)
	{
		super(k, maxLength);
	}

	mixin LikelihoodAscentTemplate!T;
}


/++
+/
class LikelihoodAscentEM(T) : ExpectationMaximization!((a, b){foreach(i, e; a) b[i] = 1/e;}, T), LikelihoodAscent!T
	if(isFloatingPoint!T)
{
	/++
	Constructor
	Params:
		k = number of components
		maxLength = maximal length of features. In terms of likelihood maximization maxLength is maximal length of a sample.
	+/
	this(size_t k, size_t maxLength)
	{
		super(k, maxLength);
	}

	mixin LikelihoodAscentTemplate!T;
}


package mixin template LikelihoodAscentTemplate(T)
{
	void put(PDFRange, SampleRange)(PDFRange pdfs, SampleRange sample)
		if (isInputRange!PDFRange && hasLength!PDFRange && isCallable!(ElementType!PDFRange))
	in
	{
		assert(pdfs.length == _featuresT.matrix.height);
	}
	body
	{
		import std.algorithm : map;
		super.put(sample.map!(x => pdfs.map!(pdf => pdf(x))));
		if (!isFeaturesCorrect)
			throw new FeaturesException;
	}

	void optimize
	(
		scope bool delegate (T likelihoodValuePrev, T likelihoodValue, in T[] weightsPrev, in T[] weights) tolerance,
		scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null,
	)
	{
		import std.functional : toDelegate;
		if (!isFeaturesCorrect)
			throw new FeaturesException;
		super.optimize(toDelegate(&_likelihood_!T), tolerance, findRootTolerance);
	}

	void optimize
	(
		scope bool delegate (T likelihoodValuePrev, T likelihoodValue) tolerance,
		scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null,
	)
	{
		import std.functional : toDelegate;
		if (!isFeaturesCorrect)
			throw new FeaturesException;
		super.optimize(toDelegate(&_likelihood_!T), tolerance, findRootTolerance);
	}

	bool isFeaturesCorrect() const
	{
		immutable c = T(1) / features.length;
		return features.transposed.all!(any!(x => isNormal(c * x)));
	}

	T likelihood() @property const
	{
		return _likelihood_!T(mixture);
	}
}


unittest {
	alias C0 = CoordinateDescent!((a, b){}, float);
	alias C3 = ExpectationMaximization!((a, b){}, float);
	alias C1 = LikelihoodAscent!float;
	alias C10 = LikelihoodAscentGradient!float;
	alias C11 = LikelihoodAscentCoordinate!float;
	alias C12 = LikelihoodAscentEM!float;
	alias C2 = GradientDescent!((a, b){}, float);
}


unittest {
	alias C0 = CoordinateDescent!((a, b){}, double);
	alias C3 = ExpectationMaximization!((a, b){}, double);
	alias C1 = LikelihoodAscent!double;
	alias C10 = LikelihoodAscentGradient!double;
	alias C11 = LikelihoodAscentCoordinate!double;
	alias C12 = LikelihoodAscentEM!double;
	alias C2 = GradientDescent!((a, b){}, double);
}


unittest {
	alias C0 = CoordinateDescent!((a, b){}, real);
	alias C3 = ExpectationMaximization!((a, b){}, real);
	alias C1 = LikelihoodAscent!real;
	alias C10 = LikelihoodAscentGradient!real;
	alias C11 = LikelihoodAscentCoordinate!real;
	alias C12 = LikelihoodAscentEM!real;
	alias C2 = GradientDescent!((a, b){}, real);
}

package T _likelihood_(T)(in T[] mixture)
{
	import atmosphere.summation : sumOfLog2s;
	return T(LN2) * mixture.sumOfLog2s / mixture.length;
}
