#include "Probability.h"
#include <cmath>

#ifdef _WIN32
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#endif

namespace mn {

    double PDF(int lambda, int k){
	    double pdf=1;
	    int i;
	    for (i=1;i<=k;++i)
		    pdf*=(double)lambda/i;
	    return pdf*exp(-1.0*lambda);
    }

	double PDF(double u, double o, int x) {
		static const double co = 1. / sqrt(2 * M_PI);
	    double index = -(x - u) * (x - u) / 2 / o / o;
		return co / o * exp(index);
	}
	double anti_normal_PDF(double u, double o, int x) {
		static const double co = 1. / sqrt(2 * M_PI);
	    double index = -(x - u) * (x - u) / 2 / o / o;
		return 1 - co / o * exp(index);
	}

    int rand_p(double lambda) {
	    double u=(double)rand()/RAND_MAX;
	    int x=0;
	    double cdf=exp(-1.0*lambda);
		while (u >= cdf) {
			x++;
			cdf += PDF(lambda, x);
		}
		return x;
    }
	int rand_normal(double u, double o) {
	    double val=(double)rand()/RAND_MAX;
	    int x=0;
	    double cdf=0;	// PDF(u, o, x)
	    while (val>=cdf) {
		    x++;
		    cdf+=PDF(u,o,x);
	    }
	    return x;
    }
	int rand_anti_normal(double u, double o) {
	    double val=(double)rand()/RAND_MAX;
	    int x=0;
	    double cdf=0;	// PDF(u, o, x)
	    while (val>=cdf) {
		    x++;
		    cdf+=anti_normal_PDF(u,o,x);
	    }
	    return x;
    }
}