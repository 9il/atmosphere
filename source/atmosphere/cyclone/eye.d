module atmosphere.cyclone.eye;


import core.stdc.string : memmove;
import std.traits;
import std.range;

import std.math : isNaN, fabs, isSubnormal;

import matrix;

///
import atmosphere.kernel.vector;
import atmosphere.kernel.findroot;
import atmosphere.kernel.memory;

debug import std.stdio;

///
struct Information(F)
if(isFloatingPoint!F)
{
	F theta = F.nan;
	Likelihood likelihood;
	sizediff_t m = sizediff_t.min;
	bool isInverted = false;
}


///
struct Eye(F, size_t D)
if(isFloatingPoint!F)
{
	///
	MatrixAllocator!F phiAll;

	///
	F[] p, psi, xi, pi, chi, gamma, omega;

	///
	F[] _xi, _pi, _chi, _gamma, _omega;

	///
	bool delegate(F lo, F hi) tolerance;

	///
	size_t[D] gridLengths;


	///
	auto ref phi() @property { return phiAll.matrix; }

	///
	this(
		size_t maxN,
		size_t k,
		size_t add4k,
		size_t[D] gridLengths,
		bool delegate(F lo, F hi) tolerance
		)
	{
		debug writeln("this begin... ");
		phiAll = typeof(phiAll)(k+add4k, maxN, k);
		p = creatAlignedArray!F(k);
		p[] = 1.0/k;
		psi = creatAlignedArray!F(k);
		_xi = creatAlignedArray!F(maxN);
		_pi = creatAlignedArray!F(maxN);
		_chi = creatAlignedArray!F(maxN);
		_gamma = creatAlignedArray!F(maxN);
		_omega = creatAlignedArray!F(maxN);
		this.tolerance = tolerance;
		this.gridLengths = gridLengths;
		tolerance = (a, b) => fabs(a-b) < 0x1p-20;
	}



	///
	void popFrontN(size_t N)
	in {
		assert(N <= phiAll.matrix.width);
	}
	body {
		phiAll.popFrontN(N);
		memmove(_chi.ptr, _chi.ptr+N, chi.length*F.sizeof);
		xi = _xi[0..phiAll.matrix.width];
		pi = _pi[0..phiAll.matrix.width];
		chi = _chi[0..phiAll.matrix.width];
		gamma = _gamma[0..phiAll.matrix.width];
		omega = _omega[0..phiAll.matrix.width];
	}



	///
	Information!F 
	putBack(ROR)(ROR ror, Information!F inf, bool preciseChi)
	if(isInputRange!ROR && isInputRange!(ElementType!ROR) && isFloatingPoint!(ElementType!(ElementType!ROR)))
	{

		debug writeln("putBack begin... ");
		size_t u;
		foreach(r; ror)
		{
			phiAll.putBack;
			auto ce = phi.column(phi.width-1);
			foreach(i; 0..phi.height)
			{
				assert(!r.empty);
				ce[i] = r.front;
				r.popFront;
			}
			assert(r.empty);			
			u++;
		}
		debug writeln("putBack 2... ");
		if(preciseChi)
		{
			u = phi.width;
		}
		if(u || preciseChi)
		{
			inf.theta = F.nan;
			xi    = _xi   [0..phi.width];
			pi    = _pi   [0..phi.width];
			chi   = _chi  [0..phi.width];
			gamma = _gamma[0..phi.width];
			omega = _omega[0..phi.width];
			debug writeln("putBack 3... ");
			auto phit = phi.transposed;
			const v = phi.width-u;
			phi.popFrontN(v);
			gemv(phit, p, chi[v..$]);
			debug writeln("putBack 33... ");
		}
		debug writeln("putBack 4... ");
		Information!F rinf;
		debug writeln("putBack 5... ");
		rinf.likelihood = likelihood(chi);
		debug writeln(rinf.likelihood);
		debug writeln("putBack 6... ");
		return rinf;
	}




	///
	Information!F opCall
	(
		Information!F inf,
		sizediff_t m, 
		bool preciseChi,
		bool inverse,
		bool clustered,
		bool normolizeP,
	)
	body
	{
		if(normolizeP)
		{
			p.scale(1/p.sum);
		}
		typeof(return) iret;
		
		iret.likelihood = inf.likelihood;
		
		if(preciseChi)
		{
			debug writeln("preciseChi");
			gemv(phi.transposed, p, chi);
			iret.likelihood = iret.likelihood.init;
		}

		if((iret.m = m) == -1)
		{
			xi.inverse(chi);
			gemv(phi, xi, psi);
			iret.m = psi.maxIndexOf;
		}

		substract(
			pi,
			phi[iret.m],
			chi
		);
		
		iret.theta = findHRoot(pi, chi, tolerance);

		if(iret.theta.isSubnormal)
			iret.theta = 0;
		if((1-iret.theta).isSubnormal)
			iret.theta = 1;

		if(iret.theta > 0)
		{
			iret.likelihood = iret.likelihood.init;

			if(iret.theta < 1)
			{
				p.scale(1-iret.theta);
				p[iret.m] += iret.theta;
				muladdmul(
					chi, 
					phi[iret.m],
					1-iret.theta, 
					iret.theta
				);
			}
			else
			{
				p[] = 0;
				p[iret.m] = 1;
				chi[] = phi[iret.m][];
			}

			if(clustered)
			{
				foreach(d; 0..D)
				{
					size_t[3] delta = void;
					delta[0] = iret.m;
					size_t[D] indexes = void;
					unpack(iret.m, gridLengths, indexes);
					indexes[d] = (indexes[d]+1) % gridLengths[d];
					delta[1] = pack(gridLengths, indexes);
					indexes[d] = (indexes[d]+gridLengths[d]-2) % gridLengths[d];
					delta[2] = pack(gridLengths, indexes);
					indexes[d] = (indexes[d]+1) % gridLengths[d];
					const cm = clusterization!(F, 3, Likelihood)(
						delta, 
						phi, 
						p, 
						chi, 
						gamma, 
						omega, 
						iret.likelihood
					);
					if(cm >= 0)
						iret.m = delta[cm];
				}
			}
		}
		else
		{
			if(inverse && p[iret.m] > 0 && p[iret.m] != 1)
			{

				gamma.muladdmul(
					chi, 
					phi[iret.m], 
					1/(1-p[iret.m]), 
					-p[iret.m]/(1-p[iret.m]),
				);

				auto lh = likelihood(gamma);

				if(iret.likelihood.isNaN)
					iret.likelihood = likelihood(chi);

				if(lh > iret.likelihood)
				{
					chi[] = gamma[];
					p.scale(1/(1-p[iret.m]));
					p[iret.m] = 0;
					iret.isInverted = true;
					iret.likelihood = lh;
				}
			}
		}
		return iret;
	}
}
