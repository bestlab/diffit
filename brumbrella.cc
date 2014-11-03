/*
 * Brownian dynamics on 1D quartic potential V0 = eps*(x**2-1)**2
 * with added umbrella bias V1 = 1/2*k1*(x-x1)**2
 *
 * x(t+dt) = x(t) + beta*D*f(t)*dt + g*(2*D*dt)**0.5
 *
 * where g is chosen from a normal distribution with unit variance
 * D is diff coeff.
 * dt is time step
 * f is force, f=-d_x(V0+V1)
 */
#include <cstdio> 
#include <cstdlib> 
#include <vector>
#include <string>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_vector.h>

using namespace std;

const char *usage = "\n\n"
"        Usage\n"
"             brumbrella -o x.dat -D diff -d dt -e eps -x x1 -k k1 -n nsteps -p nprint\n"
"\n\n";

int main(int argc, char **argv)
{
	int c;
	double D, dt, eps, x1, k1;
	double t, V0, V1, f0, f1, x, tmp, g, sig, R, seed;
	int nsteps, nprint;
	FILE *outp;
	gsl_rng *twister;
	string outp_name;
	nsteps = 1000;
	eps = 10; 	// kBT
	D = 1.0;	// D
	x1 = -1.0;	// L
	k1 = 0.0;	// kBT/L
	dt = 0.01;	// L**2/D
	nprint = 100;
	outp_name = "junk.dat";
	seed = 27041994;

	while (1) {
		c=getopt(argc,argv,"ho:D:d:e:x:k:n:p:s:");
		if (c == -1)	// no more options
			break;
		switch (c) {
			case 'h':
				fprintf(stdout,"%s\n",usage);
				exit(0);
				break;
			case 'o':
				outp_name = optarg;
				break;
			case 'D':
				D = atof(optarg);
				break;
			case 'd':
				dt = atof(optarg);
				break;
			case 'e':
				eps = atof(optarg);
				break;
			case 'x':
				x1 = atof(optarg);
				break;
			case 'k':
				k1 = atof(optarg);
				break;
			case 'n':
				nsteps = atoi(optarg);
				break;
			case 'p':
				nprint = atoi(optarg);
				break;
			case 's':
				seed = atoi(optarg);
				break;
			default:
				fprintf(stderr,"?? getopt returned character %c ??\n", c);
				fprintf(stderr,"%s\n",usage);
				exit(1);
		}
	}
	// open output
	outp = fopen(outp_name.c_str(),"w");
	// set up RNG
	twister = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(twister, seed);

	x = x1; 	// may as well start at bottom of umbrella...
	sig  = sqrt(2.0*D*dt);
	for (int s=0; s<nsteps; s++) {
		t = double(s)*dt;
		tmp = x*x-1.0;
		V0 = eps*tmp*tmp;
		f0 = -4.0*eps*tmp*x;
		tmp = x-x1;
		V1 = 0.5*k1*tmp*tmp;
		f1 = -k1*tmp;
		R = gsl_ran_gaussian(twister, sig);
		x += D*(f0+f1)*dt+R;
		if (s%nprint == 0) {
			//fprintf(stdout,"%12.6f %8.3f\n",t,x);
			fprintf(outp,"%12.6f %8.3f\n",t,x);
		}
	}
}
