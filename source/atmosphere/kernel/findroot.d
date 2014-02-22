module atmosphere.kernel.findroot;

///
import std.numeric : findRoot;

import atmosphere.kernel.vector : HImpl;

///
F findHRoot(F)(F[] pi, F[] chi, scope bool delegate(F, F) tolerance)
{
	import std.math : fabs;
	F H0 = void, H1 = void;
	scope H = (F theta) => HImpl(pi, chi, theta);
	H0 = H(0);
	if(H0 <= 0 || tolerance(H0, 0))
		return 0;
	H1 = H(1);
	if(H1 >= 1 || tolerance(H1, 1))
		return 1;		
	auto r = H.findRoot(cast(F)0, cast(F)1, H0, H1, tolerance);
	return !(fabs(r[2]) > fabs(r[3])) ? r[0] : r[1];
}