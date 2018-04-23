#include <cmath>
#include "gnuplot.h"

int main()
{
	using namespace std;
	using namespace gnuplot;

	double x[1200], y[1200], z[1200];
	for (int i = 0; i<1200; i++) {
		if (i < 1000) {
			double
				t = 0.1*i,
				s = 1 - 0.001*i;
			x[i] = s * cos(t);
			y[i] = s * sin(t);
			z[i] = t;
		}
		else {
			x[i] = 0;
			y[i] = 0;
			z[i] = 100 - (i - 1000);
		}
	}

	CGnuplot gp;
	gp.SetYRange(-2, 2);
	gp.SetXRange(-2, 2);
	gp.Plot(x, y, z);
	gp.DumpToFile("kinoko");

	__KEYWAIT__;
	return 0;
}