module atmosphere.cyclone.generator;

import std.numeric : gcd;
import std.string : format;

import atmosphere.fractal;

import atmosphere.kernel.vector : pack;


///
class GeneratorException : Exception
{

	this(string msg)
	{
		super("GeneratorException: " ~ msg);
	}

	this(size_t a, size_t b)
	{
		this(format("gcd(%s, %s) != 1", a, b));
	}
}

private void check(size_t a, size_t b)
{
	if(gcd(a, b) != 1)
		throw new GeneratorException(a, b);
}


///
struct Generator(size_t N)
{
	///
	size_t index;

	///
	size_t[N] lengths;

	///
	size_t[][N] frac;

	///
	this(size_t[N] lengths)
	{
		this.lengths = lengths;
		foreach(a; 0..N)
			foreach(b; 0..a)
				check(this.lengths[a], this.lengths[b]);
		foreach(i; 0..N)
			frac[i] = createFractale001(this.lengths[i], 1);
	}

	///
	size_t opCall()
	{		
		size_t[N] kindexes = void;
		foreach(i, ref kindex; kindexes)
			kindex = frac[i][index % lengths[i]];
		index++;
		return pack(lengths, kindexes);
	}
}
