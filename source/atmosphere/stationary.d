
/**
In most cases coordinate descent is much more faster then gradient descent.
*/
module atmosphere.stationary;

import std.range : isInputRange, hasLength;

import atmosphere.mixtureoptimizer;
import atmosphere.internal;


/**
Params:
	T = floating point type
*/
abstract class StationaryOptimizer(T) : MixtureOptimizer!T
{
	private Matrix!T _componentsT_;
	private T[] _distribution_;
	private T[] _mixture_;

final:

	/**
	Params:
		k = number of components
		n = length of each component. In terms of likelihood maximization n is length of sample.
	*/
	this(size_t k, size_t n)
	{
		_componentsT_ = Matrix!T(k, n);
		_distribution_ = new T[k];
		_distribution_[] = T(1)/k;
		_mixture_ = new T[n];
	}

	//alias components = super.components;
	/**
	This function accept components and store them in internal transposed form.
	Params:
		_components = components
	*/
	void components(in T[][] _components)
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
	The pseudocode of this method: $(D components[i, j] = kernel[i](sample[j])).
	Params:
		kernels = Input range of kernels. kernels length equals k. 
		Element type of kernels can be a $(D function) pointer, $(D delegate), 
			$(D struct) or $(D class) with $(D opCall) method.
			Kernel compute feature for each point of grid.
			In terms of likelihood maximization each kernel is a probability density function.
		grid = pointers of components. 
		In terms of likelihood maximization grid is just a sample of length n.
	*/
	void components(Kernels)(Kernels kernels, in T[] grid)
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

	inout(Matrix!T) _componentsT() inout @property
	{
		return _componentsT_;
	}

	inout(T)[] _mixture() inout @property 
	{
		return _mixture_;
	}

	inout(T)[] _distribution() inout @property 
	{
		return _distribution_;
	}

	void _distribution(T[] _distribution_) @property
	{
		this._distribution_ = _distribution_;
	}
}


/**
Params:
	Gradient = Gradient of the objective function. $(D Gradient(a, b) should perform b = grad_f(a)).
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
		gradientDescentIteration!(Gradient, T)(cast(Matrix!(const T))_componentsT_, _distribution_, _mixture_, pi, xi, gamma, c, findRootTolerance is null ? (a, b) => false : findRootTolerance);
		updateMixture;
	}
}


/**
Params:
	Gradient = Gradient of the objective function. $(D Gradient(a, b) should perform b = grad_f(a)).
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
		coordinateDescentIteration!(Gradient, T)(cast(Matrix!(const T))_componentsT_, _distribution_, _mixture_, pi, xi, gamma, findRootTolerance is null ? (a, b) => false : findRootTolerance);
		updateMixture;
	}
}


/**
Params:
	PartialDerivative = Partial derivative of the objective function. 
		There are optimization reasons for coordinate descent to use partial derivative instead of gradient.
		It is applicable if and only if all partial derivatives are equal.
		$(D b = PartialDerivative(a) should perform b = (d f/ d x_i)(a)).
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
		coordinateDescentIterationPartial!(PartialDerivative, T)(cast(Matrix!(const T))_componentsT_, _distribution_, _mixture_, pi, findRootTolerance is null ? (a, b) => false : findRootTolerance);
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


/////Example:
//unittest {
//	import std.range;

//	//import atmosphere.stationary;
//	auto optimizer = new LikelihoodMaximizationCoordinate!double(10, 100);

//	auto components = optimizer.MixtureOptimizer.components;
//	alias ComponentsType = typeof(components);
//	alias ComponentType = ElementType!ComponentsType;
	
//	static assert(isRandomAccessRange!ComponentType);
//	static assert(isRandomAccessRange!ComponentsType);
//	static assert(hasAssignableElements!ComponentType == false);
//	static assert(hasAssignableElements!ComponentsType == false);
//}