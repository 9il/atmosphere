import std.file, std.path, std.stdio, std.conv, std.algorithm, std.range, std.math, std.functional;
import atmosphere;

const begin = 0.1;
const end = 15.0;
const count = 50;
const step = (end-begin)/count;
const eps = 1e-5;

void main() 
{
	writefln("ε = %s", eps);
	const grid = iota(begin, end+step/2, step).array;
	writeln(grid);
	foreach(file; "data".dirEntries("*.txt", SpanMode.shallow))

	{
		size_t counter;
		const alpha = file.baseName(".txt").to!double;
		auto fout = File("coordinate.txt", "w");

		bool tolerance 
			(
				double alphaPrev, 
				double alpha, 
				double sumOfLog2sValuePrev, 
				double sumOfLog2sValue, 
				in double[] distributionPrev, 
				in double[] distribution,
			)
		{
			counter++;
			return sumOfLog2sValue - sumOfLog2sValuePrev <= 10*eps;
		}

		writeln("===========================");
		writeln(file);
		writefln("α in sample = %s", alpha);
		writeln("===========================");

		const sample = file.readText.splitter.map!(to!double).array;
		auto kernels = grid.map!(u => Kernel(alpha*u, sqrt(u)));
		auto cheater = new LikelihoodMaximizationCoordinate!double(kernels.length, sample.length);
		cheater.components(kernels, sample);
		foreach(i; 0..1000)
			cheater.eval;
		writefln("cheater's sumOfLog2s = %s", cheater.mixture.sumOfLog2s);
		writefln("cheater's distribution = %s", cheater.distribution);
		writeln;

		//auto optimizer = new NormalVarianceMeanMixtureEMSeparator!double(grid, sample);
		//auto optimizer = new NormalVarianceMeanMixtureEMAndGradientSeparator!double(grid, sample);
		auto optimizer = new NormalVarianceMeanMixtureEMAndCoordinateSeparator!double(grid, sample);

		writefln("mean = %s", optimizer.mean);
		writefln("α = %s", optimizer.alpha);
		writefln("sumOfLog2s = %s", optimizer.sumOfLog2s);
		writeln;
		
		writeln("optimize...");
		optimizer.optimize(toDelegate(&tolerance));
		writefln("α = %s", optimizer.alpha);
		writefln("sumOfLog2s = %s", optimizer.sumOfLog2s);
		writefln("distribution = %s", optimizer.distribution);
		writefln("total iterations: %s", counter);
		writeln;
	}
}

/**
Computes accurate sum of binary logarithms of input range $(D r).
TODO: Delete this with DMD 2.068.
 */
ElementType!Range sumOfLog2s(Range)(Range r) 
    if (isInputRange!Range && isFloatingPoint!(ElementType!Range))
{
	import std.math : frexp; 

    long exp = 0;
    Unqual!(typeof(return)) x = 1; 
    foreach (e; r)
    {
        if (e < 0)
            return typeof(return).nan;
        int lexp = void;
        x *= frexp(e, lexp);
        exp += lexp;
        if (x < 0.5) 
        {
            x *= 2;
            exp--;
        }
    }
    return exp + log2(x); 
}

struct Kernel
{
	double alphau;
	double sqrtu;

	this(double alphau, double sqrtu)
	{
		this.alphau = alphau;
		this.sqrtu = sqrtu;
	}

	double opCall(double x) const
	{
		immutable y = (x - alphau) / sqrtu;
		return exp(y * y / -2) / sqrtu;
	}
}