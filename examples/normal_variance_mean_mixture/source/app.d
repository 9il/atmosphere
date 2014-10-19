import std.file, std.path, std.stdio, std.conv, std.algorithm, std.range, std.math;
import atmosphere;

const begin = 0.1;
const end = 15.0;
const count = 50;
const step = (end-begin)/count;
const eps = 1e-5;

size_t[] index = [1, 2, 46, 48, 20, 44, 38, 17, 25, 41, 22, 10, 19, 15, 4, 28, 37, 23, 26, 43, 21, 24, 42, 30, 13, 27, 8, 40, 32, 14, 31, 6, 0, 35, 33, 12, 47, 7, 34, 11, 18, 29, 49, 36, 39, 3, 45, 5, 16, 9, ];
void main() 
{
	writeln("ε = %s", eps);
	const grid = iota(begin, end+step/2, step).indexed(index).array;
	writeln(grid);
	foreach(file; ["data/0.8.txt"])//"data".dirEntries("*.txt", SpanMode.shallow).take(1))
	{
		size_t counter;

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

			double s = 0;
			foreach(i; 0..distributionPrev.length)
				s += (distributionPrev[i]-distribution[i])^^2;

			fout.writeln(sumOfLog2sValue);
			//return s <= eps^^2;
			return sumOfLog2sValue >= -1381.84;
		}


		writeln("===========================");
		writeln(file);
		writeln("===========================");
		const sample = file.readText.splitter.map!(to!double).array;

		auto optimizer = new NormalVarianceMeanMixtureEMSeparator!double(grid, sample);
		writefln("mean = %s", optimizer.mean);
		writefln("α = %s", optimizer.alpha);
		writefln("sumOfLog2s = %s", optimizer.sumOfLog2s);
		writeln;
		
		//writeln("one iteration...");
		//optimizer.eval;
		//writefln("α = %s", optimizer.alpha);
		//writefln("sumOfLog2s = %s", optimizer.sumOfLog2s);
		//writeln;
		
		writeln("optimize...");
		optimizer.optimize(&tolerance);
		writefln("α = %s", optimizer.alpha);
		writefln("sumOfLog2s = %s", optimizer.sumOfLog2s);
		writefln("distribution = %s", optimizer.distribution);
		writefln("total iterations: %s", counter);
		writeln;
	}
}
