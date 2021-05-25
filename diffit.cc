/*
 * maximum likelihood fit of diffusive model to simulation trajectories
 */

#define VERBOSE
#define GLOBALMOVES
// this does not seem helfpul: 
//#define OFFSET

#include <vector>
#include <string>
#include <algorithm>
#include <cstdio> 
#include <cstdlib> 
#include <cmath>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>

#include "diff_model.h"

using namespace std;

void parse_cmd(const int largc, char **largv, double& D0,
		double &dD, double &dF, double &stiffness, double &relative_stiffness, 
		double &min_relative_stiffness, int &nsteps, int &nprint, 
		string &outp_name, vector<string> &mat_files, string &restart_file,
		string &restart_save_file, string &propagator_file,
		bool &pbc, bool &t0_offset,
		int &seed, double &T1, double &T0, const char *usage)
{
	// defaults
	string NONE = "none";
	D0 = 1.0;			// initial guess at D
	dD = 0.01;			// scale for random moves in ln(D)
	dF = 0.01;			// scale for random moves in F
	nsteps = 10000;			// number of mc steps
	nprint = 100;			// frq for writing output
	stiffness = -1.0;		// i.e. no stiffness
	relative_stiffness = -1.0;		// i.e. no relative stiffness
	min_relative_stiffness = -1.0;		// i.e. no min relative stiffness
	outp_name = "default.dat";
	restart_file = "none";
	propagator_file = "none";
	restart_save_file = "none";
	pbc = false;
	seed = 27041994;
	T0 = 1.0;
	T1 = 1.0;
	t0_offset = false;
	int nmat = 0;
	int c;
	//
	if (largc == 1) {
		fprintf(stdout,"%s\n",usage);
		exit(0);
	}
	while (1) {
		c=getopt(largc,largv,"ho:D:d:f:n:p:q:Q:s:S:PT:e:r:R:A:g:O");
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
				D0 = atof(optarg);
				break;
			case 'd':
				dD = atof(optarg);
				break;
			case 'f':
				dF = atof(optarg);
				break;
			case 'n':
				nsteps = atoi(optarg);
				break;
			case 'p':
				nprint = atoi(optarg);
				break;
			case 'S':
				stiffness = atof(optarg);
				break;
			case 'q':
				relative_stiffness = atof(optarg);
				break;
			case 'Q':
				min_relative_stiffness = atof(optarg);
				break;
			case 's':
				seed = atoi(optarg);
				break;
			case 'r':
				restart_file = optarg; // read for this run...
				break;
			case 'R':
				restart_save_file = optarg; // write for the next run...
				break;
			case 'g':
				propagator_file = optarg; // write for the next run...
				break;
			case 'P':
				pbc = true;
				break;
			case 'A':
				T0 = atof(optarg);
				break;
			case 'T':
				T1 = atof(optarg);
				break;
#ifdef OFFSET
			case 'O':
				t0_offset = true;
				break;
#endif // OFFSET
			default:
				fprintf(stderr,"?? getopt returned character code 0%o ??\n", c);
				fprintf(stderr,"%s\n",usage);
				exit(1);
		}
	}
	nmat = largc-optind;
	mat_files.resize(nmat);
	for (int i=0; i<nmat; i++) {
		mat_files[i] = largv[optind+i];
	}
	if (stiffness>0 && relative_stiffness >0) {
		fprintf(stderr,"Cannot set both a stiffness (-S) and a relative stiffness (-q) - you need to choose!\n");
		exit(1);
	}
	if (stiffness>0 && min_relative_stiffness >0) {
		fprintf(stderr,"Cannot set both a stiffness (-S) and a relative stiffness (-Q) - you need to choose!\n");
		exit(1);
	}
	if (relative_stiffness>0 && min_relative_stiffness >0) {
		fprintf(stderr,"Cannot set both a rel. stiffness (-q) and a min relative stiffness (-Q) - you need to choose!\n");
		exit(1);
	}
#ifdef VERBOSE
	fprintf(stdout,"==================================================\n");
	fprintf(stdout,"Input read:\n");
	fprintf(stdout,"   Output file name: %s\n",outp_name.c_str());
	fprintf(stdout," Initial guess at D: %12.6e\n",D0);
	fprintf(stdout,"     Moves in ln(D): %12.6e\n",dD);
	fprintf(stdout,"      Moves in F(Q): %12.6e\n",dF);
	fprintf(stdout,"         # MC steps: %8i\n",nsteps);
	fprintf(stdout,"    Print frequency: %8i\n",nprint);
	fprintf(stdout,"Initial Temperature: %8.3f\n",T0);
	fprintf(stdout,"  Final Temperature: %8.3f\n",T1);
	fprintf(stdout,"Transition matrix files:\n");
	for (int i=0; i<nmat; i++) {
		fprintf(stdout,"\t%s\n",mat_files[i].c_str());
	}
	if (restart_file != NONE) {
		fprintf(stdout," Reading restart from: %s\n",restart_file.c_str());
	}
	if (restart_save_file != NONE) {
		fprintf(stdout," Saving restart to: %s\n",restart_save_file.c_str());
	}
	if (propagator_file != NONE) {
		fprintf(stdout," Saving final propagators to: %s\n",propagator_file.c_str());
	}
	if (pbc) {
		fprintf(stdout,"Using periodic boundaries\n");
	} else {
		fprintf(stdout,"Not using periodic boundaries\n");
	}

	fprintf(stdout,"==================================================\n");
#endif // VERBOSE
	return;
}

void clean_matrices(vector<tmat *> &TMAT)
{
	int nmat = TMAT.size();
	for (int i=0; i<nmat; i++)
		delete [] TMAT[i];
	return;
}

void read_matrices(vector<string> &mat_files, vector<tmat *> &TMAT)
{
	TMAT.resize(mat_files.size());
	for (int i=0; i<mat_files.size(); i++)
		TMAT[i] = new tmat(mat_files[i].c_str());
	return;
}

void initialize_data(const double &D0, vector<double> &DQ,
		vector<double> &FQ, vector<double> &DQ_trial, vector<double> &FQ_trial, 
		const int &nbin, bool pbc)
{
	int D_nbin;
	//double dQ;
	//dQ = (qhi-qlo)/double(nbin);
	D_nbin = ( pbc ? nbin : nbin - 1 );
	DQ.resize(D_nbin);
	FQ.resize(nbin);
	DQ_trial.resize(D_nbin);
	FQ_trial.resize(nbin);
	for (int i=0; i<D_nbin; i++) {
		DQ[i] = DQ_trial[i] = D0; 
	}
	for (int i=0; i<nbin; i++) {
		FQ[i] = FQ_trial[i] = 0.0;
	}
	return;
}

void initialize_data_restart(const double &D0, vector<double> &DQ,
		vector<double> &FQ, vector<double> &DQ_trial, vector<double> &FQ_trial, 
		const int &nbin, double &t0, bool pbc, bool t0_offset, const string &restart_file)
{
	int D_nbin;
	FILE *rst;
	const int buflen=1024;
	char buf[buflen];
	//char *rtn;
	int rtn;
	int nbin_read, D_nbin_read;
	//double dQ;
	//dQ = (qhi-qlo)/double(nbin);
	D_nbin = ( pbc ? nbin : nbin - 1 );
	DQ.resize(D_nbin);
	FQ.resize(nbin);
	DQ_trial.resize(D_nbin);
	FQ_trial.resize(nbin);
	rst = fopen(restart_file.c_str(),"r");
	if (rst == NULL) {
		fprintf(stderr,"Could not open file %s\n",restart_file.c_str());
		exit(1);
	}
	//rtn = fread(buf,buflen,rst);
	//nbin_read = atoi(strtok(buf," "));
	rtn = fread(&nbin_read,sizeof(int),1,rst);
	//rtn = fgets(buf,buflen,rst);
	//D_nbin_read = atoi(strtok(buf," "));
	rtn = fread(&D_nbin_read,sizeof(int),1,rst);
	if (nbin_read != nbin || D_nbin_read != D_nbin) {
		fprintf(stderr,"X dimensions of restart data do not agree with transition matrices\n");
		fprintf(stderr,"X dimensions of restart data (F) = %i\n",nbin_read);
		fprintf(stderr,"        X dimensions of tmat (F) = %i\n",nbin);
		fprintf(stderr,"X dimensions of restart data (D) = %i\n",D_nbin_read);
		fprintf(stderr,"        X dimensions of tmat (D) = %i\n",D_nbin);
		fprintf(stderr,"QUITTING!\n");
		exit(1);
	}
	for (int i=0; i<nbin; i++) {
		//rtn = fgets(buf,buflen,rst);
		//FQ[i] = atof(strtok(buf," "));
		rtn = fread(&FQ[i],sizeof(double),1,rst);
	}
	for (int i=0; i<D_nbin; i++) {
		//rtn = fgets(buf,buflen,rst);
		//DQ[i] = atof(strtok(buf," "));
		rtn = fread(&DQ[i],sizeof(double),1,rst);
	}
	if (t0_offset) {
		rtn = fread(&t0, sizeof(double),1,rst);
	} else {
		t0 = 0.;
	}
	return;
}

/* create rate matrix K from scratch */
void build_K_1D(gsl_matrix *K, vector<double> &PQ, vector<double> &DQ, int n,bool pbc,double dQ)
{
	int iplus,iminus;
	double sumK;
	gsl_matrix_set_zero(K);

	for (int i=0; i<n-1; i++) {
		// off diagonal elements
		gsl_matrix_set(K,i+1,i, DQ[i]*sqrt(PQ[i+1]/PQ[i])/(dQ*dQ));
		gsl_matrix_set(K,i,i+1, DQ[i]*sqrt(PQ[i]/PQ[i+1])/(dQ*dQ));
	}
	if (pbc) { // take care of boundary elements
		gsl_matrix_set(K,0,n-1, DQ[n-1]*sqrt(PQ[0]/PQ[n-1])/(dQ*dQ));
		gsl_matrix_set(K,n-1,0, DQ[n-1]*sqrt(PQ[n-1]/PQ[0])/(dQ*dQ));
	}
	for (int i=0; i<n; i++) {
		iplus = ( i<n-1 ? i+1 : 0 );
		iminus = ( i>0 ? i-1 : n-1 );
		sumK = gsl_matrix_get(K,iminus,i)+gsl_matrix_get(K,iplus,i);
		gsl_matrix_set(K,i,i, -sumK);
	}
	//for (int i=0; i<n;i++) {
	//	for (int j=0; j<n;j++) {
	//		fprintf(stdout,"%5.3f ",gsl_matrix_get(K,i,j));
	//	}
	//	fprintf(stdout,"\n");
	//}

	return;
}

void calc_propagators(gsl_matrix *K, gsl_matrix *expKt, vector<double> &Peq,
		int n, double tlag,
		gsl_matrix *Phalf, gsl_matrix *Pminushalf, gsl_matrix *Ksymm,
		gsl_matrix *tmp_a, gsl_matrix *tmp_b, gsl_matrix *evecs,
		gsl_vector *evals, gsl_eigen_symmv_workspace *w)
{
	double sqtP;
	gsl_matrix_set_zero(Phalf);
	gsl_matrix_set_zero(Pminushalf);
	for (int i=0; i<n; i++) {
		//fprintf(stdout,"peq[%i] = %12.6f\n",i,Peq[i]);
		sqtP = sqrt(Peq[i]);
		gsl_matrix_set(Phalf,i,i,sqtP);
		gsl_matrix_set(Pminushalf,i,i,1./sqtP);
	}
	// this is how we have to multiply two freaking matrices with
	// gsl, apparently
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.,
			Pminushalf,K,0.,expKt);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1., expKt,Phalf,0.,Ksymm);
	gsl_eigen_symmv( Ksymm,evals,evecs,w);
	gsl_eigen_symmv_sort(evals,evecs,GSL_EIGEN_SORT_ABS_ASC);
	//gsl_vector_scale(evals,tlag);
	gsl_matrix_set_zero(tmp_a);
	for (int i=0;i<n;i++) {
		gsl_matrix_set(tmp_a,i,i,exp(gsl_vector_get(evals,i)*tlag));
	}
	// exp(tKsymm)*U^T
	gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.,tmp_a,evecs,0.,tmp_b);
	// exp(tKsymm) = U*exp(tKsymm)*U^T
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.,evecs,tmp_b,0.,tmp_a);
	// exp(tKsymm)*Pminushalf
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.,tmp_a,Pminushalf,0.,tmp_b);
	// exp(tK) = Phalf*exp(tKsymm)*Pminushalf
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.,Phalf,tmp_b,0.,expKt);

}

void update_1D(vector<vector<double> > &PQ, const vector<double> &FQ, vector<double> &DQ,
		const vector<vector<double> > &WQ, vector<gsl_matrix *> &K, int nbin, bool pbc, double dQ)
{
	double sum_P;
	int nmat = PQ.size();
	const double pmin=1.e-20;
	for (int u=0; u<nmat; u++) {
		sum_P = 0.0;
		for (int i=0; i<nbin; i++) {
			PQ[u][i] = exp(-FQ[i])*WQ[u][i];
			if (PQ[u][i] < pmin)
		//		fprintf(stdout,"pmin\n");
				PQ[u][i] = pmin;
			sum_P+=PQ[u][i];
		}
		for (int i=0; i<nbin; i++) {
			PQ[u][i] /= sum_P;
		}
		build_K_1D(K[u],PQ[u],DQ,nbin,pbc,dQ);
	}
	return;
}

double stiffness_prior(vector<double>& DQ, bool pbc, double stiffness)
{
	double tmp,log_stiff;
	int D_nbin;
	D_nbin = DQ.size();
	for (int i=0; i<D_nbin-1; i++) {
		tmp = DQ[i]-DQ[i+1];
		log_stiff += tmp*tmp/(2.*stiffness*stiffness);
	}
	if (pbc) {
		tmp = DQ[D_nbin-1]-DQ[0];
		log_stiff += tmp*tmp/(2.*stiffness*stiffness);
	}
	return -log_stiff;
}

double relative_stiffness_prior(vector<double>& DQ, bool pbc, double relative_stiffness)
{
	double tmp,tmp2,log_stiff;
	int D_nbin;
	D_nbin = DQ.size();
	for (int i=0; i<D_nbin-1; i++) {
		tmp = DQ[i]-DQ[i+1];
		tmp2 = relative_stiffness*(DQ[i]+DQ[i+1])/2.;
		log_stiff += tmp*tmp/(2.*tmp2*tmp2);
	}
	if (pbc) {
		tmp = DQ[D_nbin-1]-DQ[0];
		tmp2 = relative_stiffness*(DQ[D_nbin-1]+DQ[0])/2.;
		log_stiff += tmp*tmp/(2.*tmp2*tmp2);
	}
	return -log_stiff;
}

double min_relative_stiffness_prior(vector<double>& DQ, bool pbc, double min_relative_stiffness)
{
	double tmp,tmp2,log_stiff;
	int D_nbin;
	D_nbin = DQ.size();
	for (int i=0; i<D_nbin-1; i++) {
		tmp = DQ[i]-DQ[i+1];
		tmp2 = min_relative_stiffness*min(DQ[i],DQ[i+1]);
		log_stiff += tmp*tmp/(2.*tmp2*tmp2);
	}
	if (pbc) {
		tmp = DQ[D_nbin-1]-DQ[0];
		tmp2 = min_relative_stiffness*min(DQ[D_nbin-1],DQ[0]);
		log_stiff += tmp*tmp/(2.*tmp2*tmp2);
	}
	return -log_stiff;
}

double log_likelihood(vector<gsl_matrix *> K, vector<tmat *> TMAT, vector<vector<double> > &PQ,
		gsl_matrix *expKt, int n,
		gsl_matrix *Phalf, gsl_matrix *Pminushalf, gsl_matrix *Ksymm,
		gsl_matrix *tmp_a, gsl_matrix *tmp_b, gsl_matrix *evecs,
		gsl_vector *evals, gsl_eigen_symmv_workspace *w, double t0)
//		vector<double> &DQ, double stiffness, bool pbc)
{
	int D_nbin;
	int nmat = TMAT.size();
	double log_like = 0.;
	double tmp;
	double log_stiff = 0.;
	for (int u=0; u<nmat; u++) {
		calc_propagators(K[u], expKt, PQ[u], n, TMAT[u]->lag + t0,
		Phalf, Pminushalf, Ksymm, tmp_a, tmp_b, evecs, evals, w);
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; j++) {

				// abs() included here in case expKt goes < 0 for numerical reasons...
				double tmp = double(TMAT[u]->mat[i][j])*log(abs(gsl_matrix_get(expKt,i,j)));
				log_like += tmp;
				//fprintf(stdout,"%5i %5i %12i %12.6e %12.6e %12.6e\n",i,j,TMAT[u]->mat[i][j],gsl_matrix_get(expKt,i,j),tmp,log_like);
			}
		}
	}
//	if (stiffness>0.0) {
//		D_nbin = DQ.size();
//		for (int i=0; i<D_nbin-1; i++) {
//			tmp = DQ[i]-DQ[i+1];
//			log_stiff += tmp*tmp/(stiffness*stiffness);
//		}
//		if (pbc) {
//			tmp = DQ[D_nbin-1]-DQ[0];
//			log_stiff += tmp*tmp/(stiffness*stiffness);
//		}
//	}
//	log_like -= log_stiff;
	return log_like;
}


const char *usage = "\n\n"
"        Usage\n"
"             diffit -o log.dat -D D0 -d dD -f dF -n nsteps -p nprint \n"
"                   [-S stiffness | -q relative_stiffness | -Q min_relative_stiffness ] \n"
"                   [-r restart_read_file] [-R restart_save_file] \n"
"                   [-s random_seed] [-g propagator_output ] \n"
"		    [-A initial_temperature] [-T final_temperature] \n"
"                   mat1 mat2 ... matN\n"
"        format of mat files:\n"
"        line 1: nbin Qlo Qhi lag\n"
"              nbin: number of bins along coordinate Q\n"
"              Qlo,Qhi: range of Q\n"
"              lag: lag time\n"
"        line : k1 Q1\n"
"              k1: spring constant k1 (units: kBT/L**2) \n"
"                  L is the unit of the reaction coordinate Q\n"
"              Q1: umbrella V1 = 0.5*k1*(Q-Q1)**2 \n"
"        remaining lines: square matrix of transitions (nbin x nbin)\n"
"\n\n";


int main(int argc, char **argv)
{
	string outp_name,restart_file,propagator_file,restart_save_file;
	string NONE = "none";
	vector<string> mat_files;
	//vector<gsl_matrix *> tmat;	// to read the data from mat_files into
	vector<tmat *> TMAT;	// to read the data from mat_files into
	vector<double> FQ, FQ_trial, DQ, DQ_trial, Peq;	// global F(Q), D(Q)
	vector<vector <double> > PQ, PQ_trial, WQ; 			// PQ for each umbrella
	double D0,dD,dF,dQ,stiffness,relative_stiffness,min_relative_stiffness;
	double Qlo,Qhi,kQ,Q0,L;
	double log_like, log_like_trial, prior, prior_trial, sum_P, T1, T, T0, crit,rate, global_scale;
	double E, E_trial;
	double min_lag, dt0;
	int nsteps, nprint,nbin,nbin_D, nmat;
	vector <gsl_matrix *> K, K_trial;
	gsl_matrix *tmp, *Phalf, *Pminushalf, *Ksymm, *tmp_a, *tmp_b, *evecs,*expKt, *Keq;
	gsl_vector *evals;
	gsl_eigen_symmv_workspace *w;
	bool pbc, t0_offset;
	double t0 = 0.0;
	double old_t0, new_t0; 
	FILE *outp;
	int seed;
	const double pmin = 1.e-20;
	gsl_rng *twister;

	parse_cmd(argc, argv, D0, dD, dF, stiffness, relative_stiffness, min_relative_stiffness,
			nsteps, nprint, 
			outp_name, mat_files, restart_file, restart_save_file,
			propagator_file,
			pbc, t0_offset, seed, T1, T0, usage);
	read_matrices(mat_files,TMAT);

	// TODO: include consistency check between input matrices (same nbin, dQ, etc.).
	nbin = TMAT[0]->dim;
	nbin_D = ( pbc ? nbin : nbin-1 );
	Qlo = TMAT[0]->qlo;
	Qhi = TMAT[0]->qhi;
	L = Qhi - Qlo;
	dQ = TMAT[0]->dq;
	nmat = TMAT.size();
	outp = fopen(outp_name.c_str(),"w");
	fprintf(outp,"# free energies discretized at these values of X\n#");
	for (int i=0; i<nbin; i++) {
		double QQ = Qlo + (double(i)+0.5)*dQ;
		fprintf(outp,"%12.6f ",QQ);
	}
	fprintf(outp,"\n");
	fprintf(outp,"# diffusion coefficients discretized at these values of X\n#");
	for (int i=0; i<nbin_D; i++) {
		double QQ = Qlo + double(i+1)*dQ;
		fprintf(outp,"%12.6f ",QQ);
	}
	fprintf(outp,"\n");

	if (restart_file == NONE) {
		initialize_data(D0,DQ,FQ,DQ_trial,FQ_trial,nbin,pbc);
	} else {
		initialize_data_restart(D0,DQ,FQ,DQ_trial,FQ_trial,nbin,t0,pbc,t0_offset,restart_file);
	}

	// set up initial rate matrices - one for each umbrella
	PQ.resize(nmat);
	Peq.resize(nbin);
	PQ_trial.resize(nmat);
	WQ.resize(nmat);
	K.resize(nmat);
	K_trial.resize(nmat);
	tmp = gsl_matrix_alloc(nbin,nbin);
	Phalf = gsl_matrix_alloc(nbin,nbin);
	Pminushalf = gsl_matrix_alloc(nbin,nbin);
	expKt = gsl_matrix_alloc(nbin,nbin);
	Ksymm = gsl_matrix_alloc(nbin,nbin);
	Keq = gsl_matrix_alloc(nbin,nbin);
	tmp_a = gsl_matrix_alloc(nbin,nbin);
	tmp_b = gsl_matrix_alloc(nbin,nbin);
	//gsl_matrix *expKt = gsl_matrix_alloc(n,n);
	evecs = gsl_matrix_alloc(nbin,nbin);
	evals = gsl_vector_alloc(nbin);
	w = gsl_eigen_symmv_alloc(nbin);
	for (int u=0; u<TMAT.size(); u++) {
		// TODO: put this in separate function
		PQ[u].resize(nbin);
		PQ_trial[u].resize(nbin);
		WQ[u].resize(nbin);
		kQ = TMAT[u]->k1;
		Q0 = TMAT[u]->q1;
		sum_P = 0.0;
		for (int i=0; i<nbin; i++) {
			double QQ = Qlo + (float(i)+0.5)*dQ;
			double OFFQ = QQ-Q0;
			// pbc
			if (pbc) {
				double oldoffq = OFFQ;
				OFFQ -=  L*round(OFFQ/L);
				fprintf(stdout,"%8.3f %8.3f %8.3f %8.3f\n",Q0,QQ,oldoffq,OFFQ);
			}
			WQ[u][i] = dQ*exp(-0.5*kQ*OFFQ*OFFQ);
			PQ[u][i] = exp(-FQ[i])*WQ[u][i];
			if (PQ[u][i] < pmin)
				//fprintf(stdout,"pmin\n");
				PQ[u][i] = pmin;
			sum_P+=PQ[u][i];
		}
		for (int i=0; i<nbin; i++) {
			PQ[u][i] /= sum_P;
			fprintf(stdout,"u= %i, i = %i, p(i) = %12.6e\n",u,i,PQ[u][i]);
		}
		K[u] = gsl_matrix_alloc(nbin,nbin);
		K_trial[u] = gsl_matrix_alloc(nbin,nbin);
		build_K_1D(K[u],PQ[u],DQ,nbin,pbc,dQ);
	}
	
	// calculate initial log likelihood
	prior = prior_trial = 0;
	log_like = log_likelihood(K, TMAT, PQ, expKt, nbin,
		Phalf, Pminushalf, Ksymm, tmp_a, tmp_b, evecs, evals, w, t0);
	if (stiffness>0) {
		prior = stiffness_prior(DQ, pbc, stiffness);
	} else if (relative_stiffness>0) {
		prior = relative_stiffness_prior(DQ, pbc, relative_stiffness);
	} else if (min_relative_stiffness>0) {
		prior = min_relative_stiffness_prior(DQ, pbc, min_relative_stiffness);
	}
	E = - (log_like + prior);
	fprintf(stdout,"Initial log-likelihood = %12.6e\n",log_like);
	 
	twister = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(twister, seed);

	if (t0_offset) {
		min_lag = -1.;
		for (int k=0; k<TMAT.size(); k++) {
			double llag = TMAT[k] -> lag;
			if (min_lag <0) {
				min_lag= llag;
			} else if (llag < min_lag) {
				min_lag = llag;
			}
		}
		dt0 = 0.1*min_lag;
	}

	int picked_bin;
	double new_f, old_f, new_d, old_d;
	// loop over nsteps
	double tgrad = (T1-T0)/double(nsteps);

	int F_try, F_success, D_try, D_success, D_glob_try, D_glob_success, t0_try, t0_success;
	F_try = F_success = D_try = D_success = D_glob_try = D_glob_success = t0_try = t0_success =  0;
	for (int step=0; step<nsteps; step++) {
	// (i) attempt moves in either F or D or global D (sequentially)
		T = T0 + float(step)*tgrad;
		// first F:
		F_try++;
		picked_bin = gsl_rng_uniform_int(twister,nbin);
		old_f = FQ[picked_bin];
		new_f = old_f + gsl_ran_gaussian(twister,dF);
		//new_f = old_f + gsl_ran_flat(twister,-dF,dF);
		FQ[picked_bin] = new_f;
		update_1D(PQ,FQ,DQ,WQ,K,nbin,pbc,dQ);
		log_like_trial = log_likelihood(K, TMAT, PQ, expKt, nbin,
				Phalf, Pminushalf, Ksymm, tmp_a, tmp_b, evecs,
				evals, w, t0);
		if (stiffness>0) {
			prior_trial = stiffness_prior(DQ, pbc, stiffness);
		} else if (relative_stiffness>0) {
			prior_trial = relative_stiffness_prior(DQ, pbc, relative_stiffness);
		} else if (min_relative_stiffness>0) {
			prior_trial = min_relative_stiffness_prior(DQ, pbc, min_relative_stiffness);
		}
		E_trial = -( log_like_trial + prior_trial);
		if (E_trial<E) {
			//accept
			E = E_trial; 
			log_like = log_like_trial;
			prior = prior_trial;
			F_success++;
		} else {
			// use MC crit.
			crit = exp( -(E_trial-E)/T );
			if (gsl_ran_flat(twister,0,1)<crit) {
				E = E_trial; 
				log_like = log_like_trial;
				prior = prior_trial;
				F_success++;
			} else {
				// go back...
				FQ[picked_bin] = old_f;
				//	no point updating for rejections
			}
		}

		// then D:
		D_try++;
		picked_bin = gsl_rng_uniform_int(twister,nbin_D);
		old_d = DQ[picked_bin];
		//new_d = abs(old_d*(1.+gsl_ran_gaussian(twister,dD)));
		new_d = old_d*exp(gsl_ran_flat(twister,-dD,dD));
		DQ[picked_bin] = new_d;
		update_1D(PQ,FQ,DQ,WQ,K,nbin,pbc,dQ);
		log_like_trial = log_likelihood(K, TMAT, PQ, expKt, nbin,
				Phalf, Pminushalf, Ksymm, tmp_a, tmp_b, evecs,
				evals, w, t0);
		if (stiffness>0) {
			prior_trial = stiffness_prior(DQ, pbc, stiffness);
		} else if (relative_stiffness>0) {
			prior_trial = relative_stiffness_prior(DQ, pbc, relative_stiffness);
		} else if (min_relative_stiffness>0) {
			prior_trial = min_relative_stiffness_prior(DQ, pbc, min_relative_stiffness);
		}
		E_trial = -(log_like_trial + prior_trial);
		if (E_trial<E) {
			//accept
			E = E_trial;
			log_like = log_like_trial;
			prior = prior_trial;
			D_success++;
		} else {
			crit = exp( -(E_trial-E)/T );
			if (gsl_ran_flat(twister,0,1)<crit) {
				E = E_trial;
				log_like = log_like_trial;
				prior = prior_trial;
				D_success++;
			} else {
				// go back...
				DQ[picked_bin] = old_d;
			}
		}

#ifdef GLOBALMOVES
		// then D (global scaling of all D)
		D_glob_try++;
		double global_fudge = 0.1/sqrt(float(nbin));
		global_scale = exp(gsl_ran_flat(twister,-dD*global_fudge,dD*global_fudge));
		for (int q=0; q<nbin_D;q++) {
			DQ[q] *= global_scale; 
		}
		update_1D(PQ,FQ,DQ,WQ,K,nbin,pbc,dQ);
		log_like_trial = log_likelihood(K, TMAT, PQ, expKt, nbin,
				Phalf, Pminushalf, Ksymm, tmp_a, tmp_b, evecs,
				evals, w, t0);
		if (stiffness>0) {
			prior_trial = stiffness_prior(DQ, pbc, stiffness);
		} else if (relative_stiffness>0) {
			prior_trial = relative_stiffness_prior(DQ, pbc, relative_stiffness);
		} else if (min_relative_stiffness>0) {
			prior_trial = min_relative_stiffness_prior(DQ, pbc, min_relative_stiffness);
		}
		E_trial = -(log_like_trial + prior_trial);
		if (E_trial<E) {
			//accept
			E = E_trial;
			log_like = log_like_trial;
			prior = prior_trial;
			D_glob_success++;
		} else {
			crit = exp( -(E_trial-E)/T );
			if (gsl_ran_flat(twister,0,1)<crit) {
				E = E_trial;
				log_like = log_like_trial;
				prior = prior_trial;
				D_glob_success++;
			} else {
				// go back...
				for (int q=0; q<nbin_D;q++) {
					DQ[q] /= global_scale;
				}
			}
		}
#endif //GLOBALMOVES

#ifdef OFFSET
		if (t0_offset) {
			t0_try++;
			old_t0 = t0;
			new_t0 = -2.0*min_lag;
			while (new_t0 < -min_lag) {
				new_t0 = old_t0 + gsl_ran_gaussian(twister,dt0);
			}

			log_like_trial = log_likelihood(K, TMAT, PQ, expKt, nbin,
				Phalf, Pminushalf, Ksymm, tmp_a, tmp_b, evecs,
				evals, w, new_t0);
			E_trial = -(log_like_trial + prior); // prior unchanged since D unchanged
			if (E_trial<E) {
				//accept
				E = E_trial;
				log_like = log_like_trial;
				t0_success++;
				t0 = new_t0;
			} else {
				crit = exp( -(E_trial-E)/T );
				if (gsl_ran_flat(twister,0,1)<crit) {
					E = E_trial;
					log_like = log_like_trial;
					t0_success++;
					t0 = new_t0;
				}
			}

		}

#endif // OFFSET
		if (step%nprint==0) {
			sum_P = 0.0;
			for (int i=0; i<nbin; i++) {
				Peq[i] = exp(-FQ[i]);
				sum_P+=Peq[i];
			}
			for (int i=0; i<nbin; i++) {
				Peq[i] /= sum_P;
			}
			build_K_1D(Keq,Peq,DQ,nbin,pbc,dQ);
			// do this to get eigenvalues
			calc_propagators(Keq, expKt, Peq, nbin, TMAT[0]->lag,
					Phalf, Pminushalf, Ksymm, tmp_a, tmp_b, evecs, evals, w);
			rate = gsl_vector_get(evals,1);
			fprintf(stdout,"step %i : T=%5.3f; -loglike = %12.6e; prior = %12.6e; k_slowest = %12.6e\n",
					step,T,-log_like,prior,rate);
			fprintf(outp,"%8i %12.6e %12.6e ",step,-log_like,rate); // not changing this yet as would break python scripts
			for (int i=0; i<nbin; i++) {
				fprintf(outp,"%12.6e ",FQ[i]);
			}
			for (int i=0; i<nbin_D; i++) {
				fprintf(outp,"%12.6e ",DQ[i]);
			}
			//
			// always add new stuff at the end!
			//
			if (t0_offset) {
				fprintf(outp,"%12.6e ",t0);
			}
			fprintf(outp,"\n");
		}
	}
	fclose(outp);
	fprintf(stdout,"================================================================================\n" );
	fprintf(stdout,"       Fraction D moves accepted = %12.6f\n", double(D_success)/double(D_try));
	fprintf(stdout,"       Fraction F moves accepted = %12.6f\n", double(F_success)/double(F_try));
	fprintf(stdout,"Fraction global D moves accepted = %12.6f\n", double(D_glob_success)/double(D_glob_try));
	if (t0_offset) {
		fprintf(stdout,"Fraction t0 moves accepted = %12.6f\n", double(t0_success)/double(t0_try));
	}
	fprintf(stdout,"================================================================================\n" );

	if (propagator_file != NONE) {
		outp = fopen(propagator_file.c_str(),"w");
		for (int bini=0;bini<nbin; bini++) { // row bin
			double Qbini = Qlo + (float(bini)+0.5)*dQ;
			fprintf(outp,"%8.3f ",Qbini);
			for (int binj=0; binj<nbin; binj++) { // column bin (origin)
				double p_it_j0 = gsl_matrix_get(expKt,bini,binj);
				fprintf(outp,"%12.6f ",p_it_j0);
			}
			fprintf(outp,"\n");
		}
		fclose(outp);
	}

	if (restart_save_file != NONE) {
		outp = fopen(restart_save_file.c_str(),"w");
		//fprintf(outp,"%i\n",nbin);
		fwrite(&nbin,sizeof(int),1,outp);
		//fprintf(outp,"%i\n",nbin_D);
		fwrite(&nbin_D,sizeof(int),1,outp);
		for (int i=0; i<nbin; i++) 
			fwrite(&FQ[i],sizeof(double),1,outp);
			//fprintf(outp,"%12.6e\n",FQ[i]);
		for (int i=0; i<nbin_D; i++) 
			fwrite(&DQ[i],sizeof(double),1,outp);
			//fprintf(outp,"%12.6e\n",DQ[i]);
		//
		// add new stuff to the end of file to retain backwards compatibility...
		//
		fwrite(&t0,sizeof(double),1,outp);
		fclose(outp);
	}

	/*
	FILE *TMPF = fopen("tmpf.out","w");
	for (int i=0; i<nbin; i++) {
		double QQ = Qlo+(float(i)+0.5)*dQ;
		fprintf(TMPF,"%12.6f %12.6f\n",QQ,FQ[i]);
	}
	fclose(TMPF);
	FILE *TMPD = fopen("tmpd.out","w");
	for (int i=0; i<nbin_D; i++) {
		double QQ = Qlo+float(i+1)*dQ;
		fprintf(TMPD,"%12.6f %12.6f\n",QQ,DQ[i]);
	}
	fclose(TMPD);
	*/
	// FREE SPACE
	for (int u=0; u<nmat; u++) {
		gsl_matrix_free(K[u]);
		gsl_matrix_free(K_trial[u]);
	}
	gsl_matrix_free(tmp);
	gsl_matrix_free(Phalf);
	gsl_matrix_free(Pminushalf);
	gsl_matrix_free(expKt);
	gsl_matrix_free(Ksymm);
	gsl_matrix_free(Keq);
	gsl_matrix_free(tmp_a);
	gsl_matrix_free(tmp_b);
	gsl_matrix_free(evecs);
	gsl_vector_free(evals);
	gsl_eigen_symmv_free(w);
	//clean_matrices(TMAT);
}
