module atmosphere.kernel.vector;

import core.stdc.tgmath;

import std.typetuple;
import std.traits;
import std.range;
import std.math : isNormal, isNaN;
import std.algorithm;
//import std.numeric : dotProduct;


version(LDC)
{
	import ldc.intrinsics;
	import ldc.simd: shufflevector;
	import core.simd;
	pragma(LDC_inline_ir)
		R inlineIR(string s, R, P...)(P);
}


version(unittest)
{
	import std.traits;
}


//alias fdsfsdfsdef = substract!double;
//alias HImplDOUBLE = HImpl!double;	    	
//alias likelyhoodDOUBLE = likelihood!double;
//alias clusterizationDOUBLE = clusterization!(double,3,double, double[][]);

import core.memory;
import std.array;


void gemv(M, F)(in M m, in F[] a, F[] b)
in {
	assert (m.width == a.length);
	assert (m.height == b.length);
}
body {
	import matrix, cblas;

	static if(is(M : Matrix!F))
	{
		assert(m.ptr);
		assert(m.shift >= m.width);
		gemv(
			Order.RowMajor,
			Transpose.NoTrans,
			cast(blasint)b.length,
		 	cast(blasint)a.length,
			1.0,
			m.ptr,
			cast(blasint)m.shift,
			a.ptr,
			1,
			0.0,
			b.ptr,
			1);
	}
	else
	static if(is(M : Transposed!F))
	{
		assert(m.matrix.ptr);
		assert(m.matrix.shift >= m.matrix.width);
		gemv(
			Order.RowMajor,
			Transpose.Trans,
			cast(blasint)a.length,
		 	cast(blasint)b.length,
			1.0,
			m.matrix.ptr,
			cast(blasint)m.matrix.shift,
			a.ptr,
			1,
			0.0,
			b.ptr,
			1);
	}
	else static assert(0, format("gemv for %s not implimented", M.stringof));
	//foreach(i; 0..b.length)
	//{
	//	F s = 0;
	//	auto vector = matrix[i];
	//	static if(isArray!(typeof(vector)))
	//		enum size_t shift = 1; 
	//	else
	//		const shift = vector.shift;
	//	auto ptr = vector.ptr;
	//	foreach(e; a)
	//	{
	//		s += e**ptr;
	//		ptr += shift;
	//	}
	//	b[i] = s;
	//}
}

///
F HImpl(F)(in F[] pi, in F[] chi, in F theta)
in{
	assert((cast(size_t)pi.ptr & 15) == 0);
	assert((cast(size_t)chi.ptr & 15) == 0);
}
body {
	version(LDC)
	{
		pragma(msg, "fast: HImpl");
		import ldc.simd;
		import core.simd;

	    double4 ret0 = 0;
	    double4 ret1 = 0;
	    double4 ret2 = 0;
	    double4 ret3 = 0;
	    double4 theta4 = theta;

	    const L1= pi.length & -16;
	    const L2= pi.length & - 4;

		size_t i;
	 
	    for(; i < L1; i += 16)
	    {
			auto pi0 = *cast(double4*)(pi.ptr+i+0*4);
			auto pi1 = *cast(double4*)(pi.ptr+i+1*4);
			auto pi2 = *cast(double4*)(pi.ptr+i+2*4);
			auto pi3 = *cast(double4*)(pi.ptr+i+3*4);
	        ret0 += pi0 / (*cast(double4*)(chi.ptr+i+0*4) + theta4 * pi0);
	        ret1 += pi1 / (*cast(double4*)(chi.ptr+i+1*4) + theta4 * pi1);
	        ret2 += pi2 / (*cast(double4*)(chi.ptr+i+2*4) + theta4 * pi2);
	        ret3 += pi3 / (*cast(double4*)(chi.ptr+i+3*4) + theta4 * pi3);
	    }

	    for(; i < L2; i += 4)
	    {
			auto pi0 = *cast(double4*)(pi.ptr+i+0*4);
	        ret1 += pi0 / (*cast(double4*)(chi.ptr+i+0*4) + theta4 * pi0);
	    }

		ret0 = (ret0+ret1)+(ret2+ret3);
		ret0 += shufflevector!(double4, 2, 3, 0, 1)(ret0, ret0);
		F ret = ret0.array[0] + ret0.array[1];

	    for(; i < pi.length; i++)
	    {
			const pi0 = pi[i+0];
	        ret += pi0 / (chi[i+0] + theta * pi0);
	    }

	    return ret;
	}
	else
	{
	    F ret0 = 0;
	    F ret1 = 0;
	    F ret2 = 0;
	    F ret3 = 0;

	    const L1= pi.length & -4;

		size_t i;
	 
	    for(; i < L1; i += 4)
	    {
			const pi0 = pi[i+0];
			const pi1 = pi[i+1];
			const pi2 = pi[i+2];
			const pi3 = pi[i+3];
	        ret0 += pi0 / (chi[i+0] + theta * pi0);
	        ret1 += pi1 / (chi[i+1] + theta * pi1);
	        ret2 += pi2 / (chi[i+2] + theta * pi2);
	        ret3 += pi3 / (chi[i+3] + theta * pi3);
	    }

	    for(; i < pi.length; i++)
	    {
			const pi0 = pi[i+0];
	        ret0 += pi0 / (chi[i+0] + theta * pi0);
	    }

	    return (ret0+ret1)+(ret2+ret3);

	}
 }



/**
struct Likelihood
{ 
	long exp;
	ulong mant;
}
*/
struct LikelihoodImpl
{
	long exp;
	ulong mant;


	///
	this(long exp, ulong mant)
	{
		this.exp = exp;
		this.mant = mant;
	}


	///
	this(in double[] x)
	{
		enum ulong expMask
			= double.max_exp * 2L - 1 << double.mant_dig - 1;
		mant = long.min; 
		exp = -(double.max_exp-1) * x.length;
		long nan;
		foreach(ulong h; cast(ulong[])x)
		{
			nan |= h + (1UL << double.mant_dig - 1); //nan
			nan |= h; //sign
			if((h & expMask) == 0) //normalize
			{
				exp -= double.mant_dig + 1;
				double d = *cast(double*)&h;
				d *= (1UL << double.mant_dig + 1);
				h = *cast(ulong*)&d;
				nan |= expMask + h; //zero (log := -inf)
			}
			exp += h >>> double.mant_dig - 1;
			auto u  = h << 11;
			u |= long.min;
			auto ui64x2 = mul64x2(u, mant); //ui64x2 := u * mant
			ui64x2[0] += ui64x2[1] >>> 63;
			auto k = ui64x2[0] >>> 63;
			exp += k;
			mant = ui64x2[0] << (k ^ 1);
		}
		mant &= ~(nan >> 63); //set mant zero if nan or -inf
	}


	///
	F opCast(F)() const
	if(isFloatingPoint!F)
	{
		//complie time
		version(CTFE) 
			import std.math : log2;
		// fast real time
		else
			import core.stdc.tgmath : log2;
		if(isNaN)
			return F.nan;
		return exp + log2(cast(F)mant * cast(F)0.5 ^^ 63);
	}


	///
	bool isNaN() const @safe pure nothrow
	{
		return cast(long)mant >= 0;
	}


	///
	string toString() const
	{
		import std.string : format;
		if(isNaN)
			return "nan";
		return format("%s+log2(%.63f)", exp, mant * 0.5L ^^ 63);
	}


	///
	double opCmp()(const auto ref Likelihood rhs) const
	//in {
	//	assert(!isNaN);
	//	assert(!rhs.isNaN);
	//}
	body {
		import std.exception;
		if(isNaN || rhs.isNaN) return double.nan;
	         if(exp  > rhs.exp)  return  1;
		else if(exp  < rhs.exp)  return -1;
		else if(mant > rhs.mant) return  1;
		else if(mant < rhs.mant) return -1;
		else                     return  0;
	}
}


alias Likelihood = LikelihoodImpl;

F likelihood(F = Likelihood)(in double[] a)
{
	return cast(F)LikelihoodImpl(a);
}


private ulong[2] mul64x2(ulong a, ulong b)
{
	//LLVM D Compiler
	version(LDC)
	{
		return inlineIR!(`
		%x = zext i64 %0 to i128
		%y = zext i64 %1 to i128
		%z = mul i128 %x, %y
		%l = trunc i128 %z to i64  
		%t = lshr i128 %z, 64
		%h = trunc i128 %t to i64 
		%r0 = insertvalue [2 x i64] undef, i64 %h, 0
		%r1 = insertvalue [2 x i64] %r0, i64 %l, 1
		ret [2 x i64] %r1`, ulong[2])(a, b);
	}
	else
	//DigitalMars D Compiler 
	version(D_InlineAsm_X86_64)
	{
		ulong h = void, l = void;
		asm
		{
			mov	RAX, a;
			mov	RDX, b;
			mul RDX;
			mov	h, RDX;
			mov	l, RAX;
		}
		return [h, l];
	}
	else static assert(0, "mul64x2 not implimented");
}


/++
	version(GNU)
	{
		version(X86_64)
		{
			static assert(0,);
			ulong h = void, l = void;
			asm
			{
				"
				movq %[a], %%rax;
				movq %[b], %%rdx;
				mulq %%rdx;
				movq %%rdx, %[h];
				movq %%rax, %[l];
				"
				: 
				[h] "=r" h, 
				[l] "=r" l 
				: 
				[a] "r" a,
				[b] "r" b
				: 
				"%rax", 
				"%rdx";
			}
			return [h, l];
		}
		else static assert(0, "mul64x2 not implimented");
	}
	else
+/

///
//alias Likelihood = double;

//F likelihood(F)(F[] a)
//if(isFloatingPoint!F)
//{
//    F ret = 0;

//	version(likelihoodLikeProduct)
//	{
//	    const L1 = a.length & -0x10;

//	    size_t i;

//	    for(; i < L1; i += 0x10)
//	    {
//	    	version(LDC)
//	    	{
//	    		pragma(msg, "fast: likelihood");
//	    		import ldc.simd;
//	    		import core.simd;
//	    		auto m0 = *cast(double4*)(a.ptr+i+0x0);
//	    		auto m1 = *cast(double4*)(a.ptr+i+0x8);
//	    		m0     *= *cast(double4*)(a.ptr+i+0x4);
//	    		m1     *= *cast(double4*)(a.ptr+i+0xC);
//	    		m0 *= m1;
//	    		m0 *=     shufflevector!(double4, 2, 3, 0, 1)(m0, m0);
//	    		ret += log2(m0.array[0] * m0.array[1]);
//	    	}
//	    	else
//	    	{
//		        typeof(return) m0 = a[i+0x0];
//		        typeof(return) m1 = a[i+0x1];
//		        typeof(return) m2 = a[i+0x2];
//		        typeof(return) m3 = a[i+0x3];

//		        m0 *= a[i+0x4];
//		        m1 *= a[i+0x5];
//		        m2 *= a[i+0x6];
//		        m3 *= a[i+0x7];

//		        m0 *= a[i+0x8];
//		        m1 *= a[i+0x9];
//		        m2 *= a[i+0xA];
//		        m3 *= a[i+0xB];

//		        m0 *= a[i+0xC];
//		        m1 *= a[i+0xD];
//		        m2 *= a[i+0xE];
//		        m3 *= a[i+0xF];

//		        ret += log2((m0*m1)*(m2*m3));
//	    	}
//	    }

//	    if(L1 != a.length)
//	    {
//	    	typeof(return) m0 = a[i++];
	    	
//		    for(; i < a.length; i += 0x1)
//		    {
//		        m0 *= a[i+0x0];
//		    }

//		    ret += log2(m0);
//	    }
//    }
//    else
//    {
//        foreach(elem; a)
//            ret += log2(elem);
//    }
//    return ret;
//}



///
F sum(F)(in F[] a)
{
	version(LDC)
	{
		pragma(msg, "fast: sum");

	    double4 ret0 = 0;
	    double4 ret1 = 0;
	    double4 ret2 = 0;
	    double4 ret3 = 0;

	    const L1= a.length & -16;
	    const L2= a.length & - 4;

	    size_t i;

	    for(; i < L1; i += 0x10)
	    {
	        ret0 += *cast(double4*)(a.ptr+i+0*4);
	        ret1 += *cast(double4*)(a.ptr+i+1*4);
	        ret2 += *cast(double4*)(a.ptr+i+2*4);
	        ret3 += *cast(double4*)(a.ptr+i+3*4);
	    }

	    for(; i < L2; i += 0x4)
	    {
	        ret0 += *cast(double4*)(a.ptr+i+0*4);
	    }


		ret0 = (ret0+ret1)+(ret2+ret3);
		ret0 += shufflevector!(double4, 2, 3, 0, 1)(ret0, ret0);
		F ret = ret0.array[0] + ret0.array[1];

	    for(; i < a.length; i += 0x1)
	    {
	        ret += a[i];
	    }
	    
	    return ret;
	}
	else
	{
	    F ret0 = 0;
	    F ret1 = 0;
	    F ret2 = 0;
	    F ret3 = 0;

	    const L1 = a.length & -0x10;
	    const L2 = a.length & -0x4;

	    size_t i;

	    for(; i < L1; i += 0x10)
	    {
	        ret0 += a[i+0x0];
	        ret1 += a[i+0x1];
	        ret2 += a[i+0x2];
	        ret3 += a[i+0x3];
	        ret0 += a[i+0x4];
	        ret1 += a[i+0x5];
	        ret2 += a[i+0x6];
	        ret3 += a[i+0x7];
	        ret0 += a[i+0x8];
	        ret1 += a[i+0x9];
	        ret2 += a[i+0xA];
	        ret3 += a[i+0xB];
	        ret0 += a[i+0xC];
	        ret1 += a[i+0xD];
	        ret2 += a[i+0xE];
	        ret3 += a[i+0xF];
	    }

	    for(; i < L2; i += 0x4)
	    {
	        ret0 += a[i+0x0];
	        ret1 += a[i+0x1];
	        ret2 += a[i+0x2];
	        ret3 += a[i+0x3];
	    }

	    for(; i < a.length; i += 0x1)
	    {
	        ret0 += a[i+0x0];
	    }

	    return (ret0+ret1)+(ret2+ret3);
	}
}




///
void muladdmul(T)(T[] b, in T[] a, in T bc, in T ac)
in
{
    assert(b.length == a.length);
}
body
{
    foreach(i, ref be; b)
        be = be*bc+a[i]*ac;
}



///
void muladdmul(T)(T[] c, in T[] b, in T[] a, in T bc, in T ac)
in
{
    assert(b.length == a.length);
}
body
{
    foreach(i, ref ce; c)
        ce = b[i]*bc+a[i]*ac;
}



///
void addmul(T)(T[] c, in T[] a, in T ac)
in
{
    assert(c.length == a.length);
}
body
{
    foreach(i, ref ce; c)
        ce += a[i]*ac;
}


///
void addmul(T)(T[] c, in T[] b, in T[] a, in T ac)
in
{
    assert(b.length == a.length);
}
body
{
    foreach(i, ref ce; c)
        ce = b[i]+a[i]*ac;
}

void inverse(T)(T[] b, in T[] a)
in
{
    assert(b.length == a.length);
}
body
{
    foreach(i, ref be; b)
        be = 1/a[i];
}



void substract(T)(T[] c, in T[] a, in T[] b)
in
{
    assert(b.length == a.length);
    assert(c.length == a.length);
}
body
{
    foreach(i, ref ce; c)
        ce = a[i]-b[i];
}



///
size_t maxIndexOf(T)(in T[] a)
{
    size_t i;
    T ai = a[0];
    foreach(j; 1..a.length)
    {
        auto aj = a[j]; 
        if(aj > ai) 
        {
            i = j;
            ai = aj;
        }
    }
    return i;
}



///
void scale(T)(T[] a, T b)
{
    foreach(ref ae; a)
        ae *= b;
}







size_t pack(size_t N)(in auto ref size_t[N] l, in auto ref size_t[N] i) @safe pure nothrow
{
	size_t index = i[0];
	foreach(k; 1..N)
	{
		index = l[k] * index + i[k];
	}
	return index;
}


///
unittest {
	size_t[1] l1 = 10, i1 = 9;
	assert(pack!1(l1, i1) == 9);
	size_t[5] l5 = 10, i5 = [1,2,3,4,5];
	assert(pack!5(l5, i5) == 12345);
}



void unpack(size_t N)(size_t index, in auto ref size_t[N] l, ref size_t[N] i) @safe pure nothrow
{
	foreach_reverse(k; 1..N)
	{
		auto lk = l[k];
		i[k] = index % lk;
		index /= lk;
	}
	i[0] = index;
}


///
unittest {
	size_t[5] i;
	size_t[5] l = 10;
	unpack!5(12345, l, i);
	assert(i[0] == 1);
	assert(i[1] == 2);
	assert(i[2] == 3);
	assert(i[3] == 4);
	assert(i[4] == 5);

	i = 0;
	l[2] = 20;

	unpack!5(24345, l, i);
	assert(i[0] == 1);
	assert(i[1] == 2);
	assert(i[2] == 3);
	assert(i[3] == 4);
}

///
sizediff_t clusterization(F, size_t N, L, M)
(
	in size_t[N] delta, 
	M phi,
	F[] p, 
	F[] chi, 
	F[] gamma, 
	F[] omega, 
	ref L lh
)
in{
	assert(omega.length == gamma.length);
	assert(omega.length == chi.length);
	assert(omega.length);
}
body {
	F[N] pDelta = void;
	const(F)[][N] phiDelta = void;
	foreach(si; 0..N)
		pDelta[si] = p[delta[si]];
	F pgamma = pDelta[0];
	foreach(si; 1..N)
		pgamma += pDelta[si];
	if(pgamma == 0)
		return -1;
	foreach(si; 0..N)
		phiDelta[si] = phi[delta[si]];
	addmul(
		gamma,
		chi,
		phiDelta[0],
		-pDelta[0],
	);
	foreach(si; 1..N)
		if(pDelta[si] > 0)
			addmul(
				gamma,
				phiDelta[si],
				-pDelta[si],
			);
	sizediff_t maxj = -2;
	if(lh.isNaN)
		lh = likelihood!L(chi);
	foreach(si; pgamma == pDelta[0] .. N)
	{
		addmul(
			omega,
			gamma,
			phiDelta[si],
			pgamma,
		);
		auto lhc = likelihood!L(omega);
		if(lhc > lh)
		{
			maxj = si;
			chi[] = omega[];
			lh = lhc;
		}
	}
	if(maxj > 0)
	{
		foreach(si; 0..N)
			p[delta[si]] = 0;
		p[delta[maxj]] = pgamma;
	}
	return maxj;
}



//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
///UNTILS

///
auto avg(F)(F[] a)
{
	real s = 0;
	foreach(e; a) s += e;
	return s / a.length;
}


///
auto stdev(F)(F[] a, F avg = F.init)
{
	if(!avg.isNormal)
	{
		avg = a.avg();
	}
	return (a.map!(a => (a-avg)^^2).reduce!((a,b)=>a+b)/a.length).sqrt();
}