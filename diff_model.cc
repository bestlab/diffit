

#include "diff_model.h"

tmat::tmat()
{
	dim = 0; // not initialized
}

tmat::tmat(const string & mat_file)
{
	// do the readin'
	read_mat(mat_file);
}

void tmat::read_mat(const string & mat_file)
{
	// do the readin'
	const int buflen=1024;
	char buf[buflen];
	double tmp;
	char *rtn;
	int minij,maxij;
	FILE *inp = fopen(mat_file.c_str(),"r");
	if (inp == NULL) {
		fprintf(stderr,"Could not open file %s\n",mat_file.c_str());
		exit(1);
	}
	// first line
	rtn=fgets(buf,buflen,inp);
	dim = atoi(strtok(buf," "));
	fprintf(stdout,"dim = %i\n",dim);
	qlo = atof(strtok(NULL," "));
	fprintf(stdout,"qlo = %8.3f\n",qlo);
	qhi = atof(strtok(NULL," "));
	fprintf(stdout,"qhi = %8.3f\n",qhi);
	dq = (qhi-qlo)/float(dim);
	lag = atof(strtok(NULL," "));
	fprintf(stdout,"lag = %8.3f\n",lag);
	// second line
	rtn=fgets(buf,buflen,inp);
	k1 = atof(strtok(buf," "));
	fprintf(stdout,"k1 = %8.3f\n",k1);
	q1 = atof(strtok(NULL," "));
	fprintf(stdout,"q1 = %8.3f\n",q1);
	//mat = gsl_matrix_alloc(dim,dim);
	mat.resize(dim);
	bin_lo = dim; bin_hi = 0;
	for (int i=0; i<dim; i++) {
		mat[i].resize(dim);
		rtn=fgets(buf,buflen,inp);
		for (int j=0; j<dim; j++) {
			if (j==0) {
				tmp = atoi(strtok(buf," "));
			} else {
				tmp = atoi(strtok(NULL," "));
			}
			//fprintf(stdout,"%5i %8.3f\n",j,tmp);
			//gsl_matrix_set(mat,i,j,tmp);
			mat[i][j] = tmp;
			if (tmp>0) {
				minij = ( i<j ? i : j );
				maxij = ( i<j ? j : i );
				if (minij<bin_lo) {
					bin_lo = minij;
				} 
				if (maxij>bin_hi) {
					bin_hi = maxij;
				} 
			}
		}
	}
	bin_hi+=1; // so we can use it like conventional C bound
}

tmat::~tmat()
{
	// deallocate gsl data
}

