module atmosphere.utilities;

import std.traits;
import core.stdc.tgmath;

package:

public import std.numeric : findRoot;
public import std.algorithm : minPos;
public import std.math : isFinite;

import cblas;
import simple_matrix;


auto sum(Range)(Range range) if(!isArray!Range)
{
	Unqual!(ForeachType!Range) s = 0;
	foreach(elem; range)
		s += elem;
	return s;
}

F sum(F)(in F[] a)
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

auto avg(Range)(Range range)
{
	return range.sum / range.length;
}


void normalize(F)(F[] range)
{
	range[] /= range.sum;
}


void gemv(M, F)(in M m, in F[] a, F[] b)
in {
	assert (m.width == a.length);
	assert (m.height == b.length);
}
body {

	static if(is(M : Matrix!(F, GCAddRoot), bool GCAddRoot))
	{
		assert(m.ptr);
		assert(m.shift >= m.width);
		cblas.gemv(
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
	static if(is(M : Transposed!(F, GCAddRoot), bool GCAddRoot))
	{
		assert(m.matrix.ptr);
		assert(m.matrix.shift >= m.matrix.width);
		cblas.gemv(
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
	else
	{
		import std.string : format;
		static assert(0, format("gemv for %s not implimented", M.stringof));
	}
}


auto dot(Range1, Range2)(Range1 r1, Range2 r2)
{
	return cblas.dot(cast(blasint)r1.length, r1.ptr, cast(blasint)r1.shift, r2.ptr, cast(blasint)r2.shift);
}

ptrdiff_t shift(F)(F[])
{
	return 1;
}

void scal(Range, T)(Range r, T alpha)
{
	cblas.scal(cast(blasint)r.length, alpha, r.ptr, cast(blasint)r.shift);
}
