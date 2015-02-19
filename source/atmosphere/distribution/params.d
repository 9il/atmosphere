/++
Parameters conversions
+/
/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.distribution.params;

import std.traits;
import std.math;
import std.typecons;

@safe pure nothrow @nogc:

/++
<img style="width: 600px;" src="images/GIGChiPsi.png">
+/
struct GIGChiPsi(T)
	if(isFloatingPoint!T)
{
	///
	T chi;
	///
	T psi;

const @property:
	/// `sqrt(chi / psi)`
	T eta()   { return sqrt(chi / psi); }
	/// `sqrt(chi * psi)`
	T omega() { return sqrt(chi * psi); }

	/// `cast` operator overloading 
	R opCast(R : GIGChiPsi!F, F)()      { return R(chi, psi); }
	/// ditto
	R opCast(R : GIGEtaOmega!F, F)()    { return R(eta, omega); }
}

///
unittest 
{
	double _chi  = 3;
	double _psi  = 4;

	auto params  = GIGChiPsi!double(_chi, _psi);

	auto params1 = cast(GIGChiPsi    !real) params;
	auto params2 = cast(GIGEtaOmega  !real) params;

	double chi   = params.chi;
	double psi   = params.psi;

	double eta   = params.eta;
	double omega = params.omega;
}


/++
<img style="width: 600px;" src="images/GIGEtaOmega.png">
+/
struct GIGEtaOmega(T)
	if(isFloatingPoint!T)
{
	///
	T eta;
	///
	T omega;

const @property:
	/// `omega * eta`
	T chi()   { return omega * eta; }
	/// `omega / eta`
	T psi()   { return omega / eta; }

	/// `cast` operator overloading 
	R opCast(R : GIGChiPsi!F, F)()      { return R(chi, psi); }
	/// ditto
	R opCast(R : GIGEtaOmega!F, F)()    { return R(eta, omega); }
}

///
unittest 
{
	double _eta  = 3;
	double _omega= 4;

	auto params  = GIGEtaOmega!double(_eta, _omega);

	auto params1 = cast(GIGChiPsi    !real) params;
	auto params2 = cast(GIGEtaOmega  !real) params;

	double chi   = params.chi;
	double psi   = params.psi;

	double eta   = params.eta;
	double omega = params.omega;
}


/++
<img style="width: 600px;" src="images/GHypAlphaDelta.png">
<img style="width: 150px;" src="images/GHypYx.png">
+/
struct GHypAlphaDelta(T)
	if(isFloatingPoint!T)
{
	///
	T alpha;
	///
	T beta;
	///
	T delta;

const @property:
	/// `delta^^2`
	T chi()   { return delta^^2; }
	/// `alpha^^2 - beta^^2`
	T psi()   { return alpha^^2 - beta^^2; }
	/// `sqrt(chi / psi)`
	T eta()   { return sqrt(chi / psi); }
	/// `sqrt(chi * psi)`
	T omega() { return sqrt(chi * psi); }

	/// `cast` operator overloading 
	R opCast(R : GHypAlphaDelta!F, F)() { return R(alpha, beta, delta); }
	/// ditto 
	R opCast(R : GHypChiPsi!F, F)()     { return R(beta, chi, psi); }
	/// ditto
	R opCast(R : GHypEtaOmega!F, F)()   { return R(beta, eta, omega); }
	/// ditto
	R opCast(R : GIGChiPsi!F, F)()      { return R(chi, psi); }
	/// ditto
	R opCast(R : GIGEtaOmega!F, F)()    { return R(eta, omega); }
}

///
unittest 
{
	double beta  = 2;
	double _alpha= 3;
	double _delta= 4;

	auto params  = GHypAlphaDelta!double(_alpha, beta, _delta);

	auto params0 = cast(GHypAlphaDelta!real) params;
	auto params1 = cast(GHypChiPsi    !real) params;
	auto params2 = cast(GHypEtaOmega  !real) params;
	auto params3 = cast(GIGChiPsi     !real) params;
	auto params4 = cast(GIGEtaOmega   !real) params;

	double alpha = params.alpha;
	double delta = params.delta;

	double chi   = params.chi;
	double psi   = params.psi;

	double eta   = params.eta;
	double omega = params.omega;
}


/++
<img style="width: 400px;" src="images/GHypNVMM.png">
+/
struct GHypChiPsi(T)
	if(isFloatingPoint!T)
{
	///
	T beta; 
	///
	T chi;
	///
	T psi;

const @property:
	/// `sqrt(chi / psi)`
	T eta()   { return sqrt(chi / psi); }
	/// `sqrt(chi * psi)`
	T omega() { return sqrt(chi * psi); }
	/// `sqrt(chi)`
	T delta() { return sqrt(chi); }
	/// `sqrt(psi + beta^^2)`
	T alpha() { return sqrt(psi + beta^^2); }

	/// `cast` operator overloading 
	R opCast(R : GHypAlphaDelta!F, F)() { return R(alpha, beta, delta); }
	/// ditto 
	R opCast(R : GHypChiPsi!F, F)()     { return R(beta, chi, psi); }
	/// ditto
	R opCast(R : GHypEtaOmega!F, F)()   { return R(beta, eta, omega); }
	/// ditto
	R opCast(R : GIGChiPsi!F, F)()      { return R(chi, psi); }
	/// ditto
	R opCast(R : GIGEtaOmega!F, F)()    { return R(eta, omega); }
}

///
unittest 
{
	double beta  = 2;
	double _chi  = 3;
	double _psi  = 4;

	auto params  = GHypChiPsi!double(beta, _chi, _psi);

	auto params0 = cast(GHypAlphaDelta!real) params;
	auto params1 = cast(GHypChiPsi    !real) params;
	auto params2 = cast(GHypEtaOmega  !real) params;
	auto params3 = cast(GIGChiPsi     !real) params;
	auto params4 = cast(GIGEtaOmega   !real) params;

	double alpha = params.alpha;
	double delta = params.delta;

	double chi   = params.chi;
	double psi   = params.psi;

	double eta   = params.eta;
	double omega = params.omega;
}


/// `beta, eta, omega` generalized hyperbolic distribution parameters
struct GHypEtaOmega(T)
	if(isFloatingPoint!T)
{
	///
	T beta;
	///
	T eta;
	///
	T omega;

const @property:
	/// `omega * eta`
	T chi()   { return omega * eta; }
	/// `omega / eta`
	T psi()   { return omega / eta; }
	/// `sqrt(chi)`
	T delta() { return sqrt(chi); }
	/// `sqrt(psi + beta^^2)`
	T alpha() { return sqrt(psi + beta^^2); }

	/// `cast` operator overloading 
	R opCast(R : GHypAlphaDelta!F, F)() { return R(alpha, beta, delta); }
	/// ditto 
	R opCast(R : GHypChiPsi!F, F)()     { return R(beta, chi, psi); }
	/// ditto
	R opCast(R : GHypEtaOmega!F, F)()   { return R(beta, eta, omega); }
	/// ditto
	R opCast(R : GIGChiPsi!F, F)()      { return R(chi, psi); }
	/// ditto
	R opCast(R : GIGEtaOmega!F, F)()    { return R(eta, omega); }
}

///
unittest 
{
	double beta  = 2;
	double _eta  = 3;
	double _omega= 4;

	auto params  = GHypEtaOmega!double(beta, _eta, _omega);

	auto params0 = cast(GHypAlphaDelta!real) params;
	auto params1 = cast(GHypChiPsi    !real) params;
	auto params2 = cast(GHypEtaOmega  !real) params;
	auto params3 = cast(GIGChiPsi     !real) params;
	auto params4 = cast(GIGEtaOmega   !real) params;

	double alpha = params.alpha;
	double delta = params.delta;

	double chi   = params.chi;
	double psi   = params.psi;

	double eta   = params.eta;
	double omega = params.omega;
}
