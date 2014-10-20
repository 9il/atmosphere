
/**
In most cases coordinate descent is much more faster then gradient descent.
*/
module atmosphere.stationary;

import std.range : isInputRange, hasLength;

import atmosphere.mixtureoptimizer;
import atmosphere.internal;
import simple_matrix;


/**
Params:
	T = floating point type
*/
abstract class StationaryOptimizer(T) : MixtureOptimizer!T
{

private:

	Matrix!T _componentsT;

package:

	T[] _distribution;
	T[] _mixture;

final:

	void updateMixture()
	{
		mix(cast(Matrix!(const T))_componentsT, _distribution, _mixture);
	}

	Matrix!(const T) componentsT()
	{
		return cast(typeof(return))_componentsT;
	}

public:

	/**
	Params:
		k = number of components
		n = length of each component. In terms of likelihood maximization n is length of sample.
	*/
	this(size_t k, size_t n)
	{
		_componentsT = Matrix!T(k, n);
		_distribution = new T[k];
		_distribution[] = T(1)/k;
		_mixture = new T[n];
	}

	/**
	This function accept components and store them in internal transposed form.
	Params:
		_components = components
	*/
	void components(in T[][] _components) @property
	in
	{
		assert(_components.length == _componentsT.width);
		foreach(c; _components)
			assert(c.length == _componentsT.height);
	}
	body
	{
		foreach(j, w; _components)
			foreach(i, e; w)
				_componentsT[i, j] = e;
		updateMixture;
	}

	/**
	Compute internal components using kernels and grid. 
	The pseudocode of this mathod: $(D components[i, j] = kernel[i](sample[j])).
	Params:
		kernels = Input range of kernels. kernels length equals k. 
		Element type of kernels can be a $(D function) pointer, $(D delegate), 
			$(D struct) or $(D class) with $(D opCall) method.
			Kernel compute feature for each point of grid.
			In terms of likelihood maximization each kernel is a probability density function.
		grid = pointers of components. 
		In terms of likelihood maximization grid is just a sample of length n.
	*/
	void components(Kernels)(Kernels kernels, in T[] grid) @property
		if(isInputRange!Kernels && hasLength!Kernels)
	in
	{
		assert(kernels.length == _componentsT.width);
		assert(grid.length == _componentsT.height);
	}
	body
	{
		auto m = _componentsT;
		foreach(kernel; kernels)
		{
			auto r = m.front;
			foreach(i, e; grid)
				r[i] = kernel(e);
			m.popFront;
		}
		updateMixture;
	}

override:

	T[][] components() @property const
	{
		return _componentsT.transpose.arrays;
	}

	const(T)[] mixture() @property const
	{
		return _mixture;
	}

	const(T)[] distribution() @property const
	{
		return _distribution;
	}

	void distribution(in T[] _distribution) @property
	{
		this._distribution[] = _distribution[];
		updateMixture;
	}
}


/**
Params:
	Gradient = Gradient of the objective function. $(D Gradient(a, b) should perfrom b = grad_f(a)).
	T = floating point type
*/
final class GradientDescent(alias Gradient, T) : StationaryOptimizer!T
{
	private T[] pi;
	private T[] xi;
	private T[] gamma;
	private T[] c;

	///
	this(size_t k, size_t n)
	{
		super(k, n);
		pi = new T[n];
		xi = new T[n];
		gamma = new T[n];
		c = new T[k];
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
		gradientDescentIteration!(Gradient, T)(cast(Matrix!(const T))_componentsT, _distribution, _mixture, pi, xi, gamma, c, findRootTolerance is null ? (a, b) => false : findRootTolerance);
		updateMixture;
	}
}


/**
Params:
	Gradient = Gradient of the objective function. $(D Gradient(a, b) should perfrom b = grad_f(a)).
	T = floating point type
*/
final class CoordinateDescent(alias Gradient, T) : StationaryOptimizer!T
{
	private T[] pi;
	private T[] xi;
	private T[] gamma;

	///
	this(size_t k, size_t n)
	{
		super(k, n);
		pi = new T[n];
		xi = new T[n];
		gamma = new T[n];
	}

	~this()
	{
		pi.destroy;
		xi.destroy;
		gamma.destroy;
	}

	override void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null)
	{
		coordinateDescentIteration!(Gradient, T)(cast(Matrix!(const T))_componentsT, _distribution, _mixture, pi, xi, gamma, findRootTolerance is null ? (a, b) => false : findRootTolerance);
		updateMixture;
	}
}


/**
Params:
	PartialDerivative = Partial derivative of the objective function. 
		There are optimization reasons for coordinate descent to use partial derivative instead of gradient.
		It is applicable if and only if all partial derivatives are equal.
		$(D b = PartialDerivative(a) equals b = (d f/ d x_i)(a)).
	T = floating point type
*/
final class CoordinateDescentPartial(alias PartialDerivative, T) : StationaryOptimizer!T
{
	private T[] pi;

	///
	this(size_t k, size_t n)
	{
		super(k, n);
		pi = new T[n];
	}

	~this()
	{
		pi.destroy;
	}

	override void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null)
	{
		coordinateDescentIterationPartial!(PartialDerivative, T)(cast(Matrix!(const T))_componentsT, _distribution, _mixture, pi, findRootTolerance is null ? (a, b) => false : findRootTolerance);
		updateMixture;
	}
}

/**
Likelihood gradient maximization.
----------
alias LikelihoodMaximizationCoordinate(T) = CoordinateDescentPartial!(a => -1/a, T);
----------
*/
alias LikelihoodMaximizationCoordinate(T) = CoordinateDescentPartial!(a => -1/a, T);


/**
Likelihood gradient maximization.
----------
alias LikelihoodMaximizationGradient(T) = GradientDescent!((a, b) {foreach(i, ai; a) b[i] = -1/ai;});
----------
*/
alias LikelihoodMaximizationGradient(T) = GradientDescent!((a, b) {foreach(i, ai; a) b[i] = -1/ai;}, T);


unittest {
	alias C0 = CoordinateDescent!((a, b){}, double);
	alias C1 = LikelihoodMaximizationGradient!(double);
	alias C2 = LikelihoodMaximizationCoordinate!(double);
}
