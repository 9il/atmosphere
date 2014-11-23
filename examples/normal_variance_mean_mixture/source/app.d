import std.file, std.path, std.stdio, std.conv, std.algorithm, 
	std.range, std.math, std.functional, std.datetime;

import atmosphere;

const begin = 0.1;
const end = 15.0;
const count = 50;
const step = (end-begin)/count;
const eps = 5e-5;

void main()
{
	//uniform grid
	const double[] grid = iota(begin, end+step/2, step).array;
	writefln("left bound = %s", begin);
	writefln("right bound = %s", end);
	writefln("step = %s", step);
	writefln("grid = %s", grid);

	//for each *.txt file in data folder
	foreach(file; "data".dirEntries("*.txt", SpanMode.shallow))
	{
		///reads alpha from file name
		const alpha = file.baseName(".txt").to!double;

		writeln("DATA =========================================");
		writefln(file);
		writefln("α in sample = %s", alpha);
		writeln("==============================================\n");

		//reads sample from file
		const double[] sample = file.readText.splitter.map!(to!double).array;

		//iteration counter
		uint counter;
		auto sw = StopWatch();

///Common algorithm

		//finds mixture weights
		writeln("mixture optimization =========================");
		//compute range of PDFs
		auto pdfs = grid.map!(u => PDF(alpha, u)).array;
		///runtime parameters: probability density functions (up to a common constant), sample
		auto optimizer = new CoordinateLikelihoodMaximization!double(pdfs.length, sample.length);
		optimizer.put(pdfs, sample);
		sw.start;
		optimizer.optimize( ///optimization
			//tolerance
			(sumOfLog2sValuePrev, sumOfLog2sValue) 
			{
				counter++;
				return sumOfLog2sValue - sumOfLog2sValuePrev <= eps;
			});
		sw.stop;
		writefln("time: %s ms", sw.peek.msecs);
		writefln("total iterations: %s", counter);
		writefln("log2Likelihood = %s", optimizer.mixture.sumOfLog2s);
		writefln("-----------\nweights = %s", optimizer.weights);
		writeln("==============================================\n");
		counter = 0;
		sw.reset;

///Special α-parametrized EM algorithm:

		///finds good (possibly not the best) value of parameter alpha and mixture weights
		auto spacialEMOptimizer = new NormalVarianceMeanMixtureEMAndCoordinate!double(grid, sample.length);
		spacialEMOptimizer.sample = sample;
		writeln("α-parametrized EM mixture optimization =======");
		sw.start;
		spacialEMOptimizer.optimize( ///optimization
			//tolerance
			(alphaPrev, alpha, double log2LikelihoodPrev, double log2Likelihood)
			{
				counter++;
				return log2Likelihood - log2LikelihoodPrev <= eps;
			});
		sw.stop;
		writefln("time: %s ms", sw.peek.msecs);
		writefln("total iterations: %s", counter);
		writefln("α = %s", spacialEMOptimizer.alpha);
		writefln("log2Likelihood = %s", spacialEMOptimizer.log2Likelihood);
		writefln("-----------\nweights = %s", spacialEMOptimizer.weights);
		//writefln("-----------\ngrid = %s", optimizer.grid);
		writeln("==============================================\n");
	}
}

//probability density function
struct PDF
{
	double alphau;
	double sqrtu;

	this(double alpha, double u)
	{
		alphau = alpha * u;
		sqrtu = sqrt(u);
	}

	///call operator overloading
	double opCall(double x) const
	{
		immutable y = (x - alphau) / sqrtu;
		//up to a constant!
		return exp(y * y / -2) / sqrtu;
	}
}


///Computes accurate sum of binary logarithms of input range $(D r).
//Will be available in std.numeric with DMD 2.068.
double sumOfLog2s(in double[] r) 
{
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
