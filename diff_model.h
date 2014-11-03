
#ifndef __DIFF_MODEL
#define __DIFF_MODEL

#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

using namespace std;

class tmat {
	public:
		//gsl_matrix *mat;
		vector<vector<int> > mat;
		double qlo, qhi, q1, k1, lag, dq;
		int dim, bin_lo, bin_hi;
		// bin_lo and bin_hi are the lowest/highest bins which actually have data!
		tmat();
		tmat(const string &mat_file);
		~tmat();
		void read_mat(const string &mat_file);
};


#endif	//ifndef __DIFF_MODEL
