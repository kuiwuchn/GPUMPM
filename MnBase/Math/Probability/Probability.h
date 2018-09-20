#ifndef __PROBABILITY_H_
#define __PROBABILITY_H_

#include <cmath>
#include <cstdlib>

namespace mn {

    double PDF(int lambda, int k);
    int rand_p(double lambda);
	double anti_normal_PDF(double u, double o, int x);

	double PDF(double u, double o, int x);
    int rand_normal(double u, double o);
	int rand_anti_normal(double u, double o);
}

#endif