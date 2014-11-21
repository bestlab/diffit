/*
 * code for 1D diffusion on a surface defined by discrete F(Q_i), D(Q_i)
 */

#include <cstdio>
#include <cstdlib>
#include <string>
#include <getopt.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "ap.h"
#include "spline3.h"

using namespace std;

const char *Usage = "\n\nUsage: "
" to be written...\n\n";

void parse_1d(const char *f, ap::real_1d_array &x,
	ap::real_1d_array &y, int &n)
{
	const int buf_len = 1024;
	char buf[buf_len];
	double xi, yi;
	int i;
	//
	FILE *inp = fopen(f,"r");
	if (inp == NULL) {
		fprintf(stderr,"Couldn't open file %s\n",f);
		exit(1);
	}

	i = 0;
	while (fgets(buf,buf_len, inp) != NULL) {
		if (buf[0] == '#' || strlen(buf) <= 1)
			continue;
		i++;
	}
	n = i; //-1;
	//fprintf(stdout,"n = %i\n",n);
	x.setbounds(0,n-1);
	y.setbounds(0,n-1);
	rewind(inp);
	i = 0;
	while (fgets(buf,buf_len, inp) != NULL) {
		if (buf[0] == '#' || strlen(buf) <= 1)
			continue;
		xi = atof(strtok(buf," \t"));
		yi = atof(strtok(NULL," \t"));
		x(i) = xi;
		y(i) = yi;
		//fprintf(stdout,"%12.6e %12.6e\n",xi,yi);
		i++;
	}
	fclose(inp);
}

int main(int argc, char **argv)
{
	ap::real_1d_array xf, f, xd, d, f_coeffs, d_coeffs; 
	string f_file, d_file;
	int c, nf, nd;
	double dt =1.;
	int nsteps = 1000;
	double x, xlo, xhi, x0, dR;
	int seed=27041994;
	int nprint=1000;
	x0 = 0.5;
	gsl_rng *twister;

	if (argc == 1) {
		fprintf(stdout,"%s\n",Usage);
	}
	
	while (1) {
		c = getopt(argc,argv,"f:d:t:n:x:p:");
		if (c == -1) 
			break;
		switch(c) {
			case 'f':
				f_file = string(optarg);
				break;
			case 'd':
				d_file = string(optarg);
				break;
			case 't':
				dt = atof(optarg);
				break;
			case 'n':
				nsteps = atoi(optarg);
				break;
			case 'p':
				nprint = atoi(optarg);
				break;
			case 'x':
				x0 = atof(optarg);
				break;
			default:
				fprintf(stderr,"?? getopt returned character code 0%o ??\n", c);
				fprintf(stderr,"%s\n",Usage);
				exit(1);


		}
	}

	parse_1d(f_file.c_str(), xf, f, nf);
	parse_1d(d_file.c_str(), xd, d, nd);
	xlo = max(xf(0),xd(0));
	xhi = min(xf(nf-1),xd(nd-1));

	//fprintf(stdout,"%8.3f %8.3f %8.3f\n",xf(0),xd(0),xlo);
	//fprintf(stdout,"%8.3f %8.3f %8.3f\n",xf(nf-1),xd(nd-1),xhi);

	buildcubicspline(xf,f,nf,0,0.,0,0.,f_coeffs);
	//buildcubicspline(xd,d,nd,0,0.,0,0.,d_coeffs);
	buildcubicspline(xd,d,nd,2,0.,2,0.,d_coeffs);
	buildcubicspline(xf,f,nf,2,0.,2,0.,f_coeffs);

	int nx = 1000;
	double dx = (xhi-xlo)/double(nx);
	double F, dF, d2F, D, dD, d2D;
	//for (int k=0; k<nx; k++) {
	//	double xk = xlo+double(k)*dx;
	//	splinedifferentiation(f_coeffs,xk,F,dF,d2F);
	//	splinedifferentiation(d_coeffs,xk,D,dD,d2D);
	//	fprintf(stdout,"%8.3f %12.6e %12.6e %12.6e\n",xk,D,F,dF);
	//}
	
	// create rng:
	twister = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(twister, seed);
	//
	x = x0;
	for (int step=0; step<nsteps; step++) {
		if (step%nprint ==0) {
			fprintf(stdout,"%8.3f %8.3f\n",double(step)*dt,x);
		}
		splinedifferentiation(f_coeffs,x,F,dF,d2F);
		splinedifferentiation(d_coeffs,x,D,dD,d2D);
		dR = sqrt(2.*D*dt);
		x += -dF*D*dt + gsl_ran_gaussian(twister,dR);
	}
}
