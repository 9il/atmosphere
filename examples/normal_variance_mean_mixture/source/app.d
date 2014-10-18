import std.file, std.path, std.stdio, std.conv, std.algorithm, std.range, std.math;
import atmosphere.normal_variance_mean_mixture;

void main() 
{
	const begin = 0.1;
	const end = 12.0;
	const count = 50;
	const step = (end-begin)/count;
	const grid = iota(begin, end+step/2, step).array;
	const eps = 1e-5;

	writeln(grid);
	
	bool tolerance(double alphaSave, double alpha, in double[] pSave, in double[] p) @nogc nothrow
	out(result)
	{
		//printf("result =%i\n", result);
	}
	body
	{
		double s = 0;
		foreach(i; 0..p.length)
			s += (p[i]-pSave[i])^^2;
		//printf("s =%f\n", sqrt(s));
		return s <= eps^^2;
	}

	bool findRootTolerance(double a, double b) @nogc nothrow
	{
		return abs(a-b) <= eps;
	}

	foreach(file; "data".dirEntries("*.txt", SpanMode.shallow).take(1))
	with(SNVMMAlgorithm)
	{
		const sample = file.readText.splitter.map!(to!double).array;
		auto p = new double[grid.length];
		double alpha;

		writeln(file);

		p[] = 1.0/p.length;
		alpha = separateNormalVarianceMeanMixture!(ExpectationMaximization)(sample, grid, p, &tolerance, &findRootTolerance);
		writeln(alpha, ' ', p);

		p[] = 1.0/p.length;
		alpha = separateNormalVarianceMeanMixture!(GradientDescent)(sample, grid, p, &tolerance, &findRootTolerance);
		writeln(alpha, ' ', p);

		p[] = 1.0/p.length;
		alpha = separateNormalVarianceMeanMixture!(CoordinateDescent)(sample, grid, p, &tolerance, &findRootTolerance);
		writeln(alpha, ' ', p);

		writeln("========================");
		writeln("========================");
		writeln("========================");

	}
}