/*
 * Bayesian construction of diffusion/free-energy profiles
 * from order parameters calculated over MD simulation
 *
 * thist.cc - just runs over the simulation and
 * constructs the transition matrix p(j,t|i,0) for a 1-D 
 * coordinate 'Q'
 */

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <cmath>


using namespace std;

const string Usage = "\n\
	Usage:\n\
		thist -h\n\
			(this message)\n\
		thist -c col -L lo -H hi -n nbin -o outp.dat -t dt [-C] rx1.dat ... rxN.dat\n\
			rx1...rxN.dat = files with reaction coordinate data\n\
			-c col = column in rx1...rxN.dat to take rc data from\n\
			-L lo  = lower histogram bound\n\
			-H hi = upper histogram bound\n\
			-n nbin = number of bins in histogram\n\
			-t dt = time step for constructing Markov matrix (an integer!)\n\
			-o outp.dat = output matrix\n\
			-C [optional] : center data in bins\n\
			";

void bin_transitions(const string &ifile, vector<vector<int> > &hist, int col, double lo,
		double hi, int nbin, int dt, bool centre_data)
{
	double dQ, Qi, Qj, Qic, Qjc;
	vector<double> Qbuf;
	ifstream inp;
	const int buf_len = 1024;
	char buf[buf_len];
	string tmps;
	double tmpf;
	int dim, Qptr, bin_i, bin_j;
	Qbuf.resize(dt);
	if (hist.size() != nbin) {
		hist.resize(nbin);
		for (int i=0;i<nbin;i++) {
			hist[i].resize(nbin);
			for (int j=0; j<nbin;j++) {
				hist[i][j] = 0;
			}
		}
	}
	dQ = (hi-lo)/double(nbin);
	inp.open(ifile.c_str());
	if (!inp.good()) {
		fprintf(stderr,"Could not open file %s to read data\n",ifile.c_str());
		exit(1);
	}
	dim = 0;
	while (!inp.eof()) {
		inp.getline(buf,buf_len,'\n');
		dim++;
	}
	dim--;
	inp.close();
	inp.clear();
	inp.open(ifile.c_str());
	for (int i=0; i<dt; i++) {	// fill Qbuf
		for (int t=1; t<col; t++)
			inp >> tmps;
		inp >> Qbuf[i];
		inp.getline(buf,buf_len,'\n');
	}
	Qptr = 0;
	for (int i=0; i<dim-dt; i++) {
		for (int t=1; t<col; t++)
			inp >> tmps;
		inp >> Qi;
		Qj = Qbuf[Qptr];
		Qbuf[Qptr] = Qi;
		Qptr++; 
		if (Qptr >= dt)
			Qptr = 0;
		if (centre_data) {
			bin_j = int(floor((Qj-lo)/dQ));
			Qjc = lo + dQ*(float(bin_j)+0.5);
			Qic = Qi + (Qjc-Qj);
			bin_i = int(floor((Qic-lo)/dQ));
		} else {
			bin_i = int(floor((Qi-lo)/dQ));
			bin_j = int(floor((Qj-lo)/dQ));
		}
		if (bin_i>=0 && bin_i<nbin && bin_j>=0 && bin_j<nbin) {
			hist[bin_i][bin_j]++;
		}
		//cout << "xx: " << tmpf << endl;
		inp.getline(buf,buf_len,'\n');
	}
	inp.close();
	return;
}

/* write matrix in indexed format */
void write_mat(const string &ofile, double lo, double hi, int nbin, int dt, 
		vector<vector<int> > &hist)
{
	FILE *outp = fopen(ofile.c_str(), "w");
	fprintf(outp,"%12.5f %12.5f %i %i\n",lo,hi,nbin,dt);
	for (int i=0; i<nbin; i++) {
		for (int j=0; j<nbin; j++) {
			fprintf(outp,"%5i %5i %12i\n",i,j,hist[i][j]);
		}
	}
	fclose(outp);
	return;
}

/* write matrix in square format with extra info */
void write_sqmat_extra(const string &ofile,vector<vector<int> > &hist,
		double lo, double hi, int nbin, double dt,
		double k1, double q1)
{
	FILE *outp = fopen(ofile.c_str(), "w");
	// write binning info 
	fprintf(outp,"%5i %12.5e %12.5e %12.5e\n",nbin,lo,hi,dt);
	// write umbrella info 
	fprintf(outp,"%12.5e %12.5e\n",k1,q1);
	//
	for (int i=0; i<nbin; i++) {
		for (int j=0; j<nbin; j++) {
			fprintf(outp," %i",hist[j][i]);
		}
		fprintf(outp,"\n");
	}
	fclose(outp);
	return;
}

/* write matrix in square format */
void write_sqmat(const string &ofile,vector<vector<int> > &hist)
{
	FILE *outp = fopen(ofile.c_str(), "w");
	int nbin = hist.size();
	for (int i=0; i<nbin; i++) {
		for (int j=0; j<nbin; j++) {
			fprintf(outp," %i",hist[j][i]);
		}
		fprintf(outp,"\n");
	}
	fclose(outp);
	return;
}

int main(int argc, char *argv[])
{
	vector<string> files;
	vector<vector<int> > hist;
	int nbin, nfile, col, dt, c;
	double lo, hi,tscale,k1,q1;
	bool centre_data,gh;
	string ofile;
	dt = 1;
	nbin = 50;
	lo = 0.0; hi = 1.0;
	col = 1;
	ofile = "junk.dat";
	centre_data = false;
	k1 = q1 =0;	// default: no umbrella
	tscale = 1.0;	// time between saved values of coordinate
	gh = false;
	while (1) {
		c=getopt(argc,argv,"hc:L:H:n:o:t:Cs:k:q:g");
		if (c == -1)	// no more options
			break;
		switch (c) {
			case 'h':
				fprintf(stdout,"%s\n",Usage.c_str());
				exit(0);
				break;
			case 'L':
				lo = atof(optarg);
				break;
			case 'H':
				hi = atof(optarg);
				break;
			case 't':
				dt = atoi(optarg);
				break;
			case 's':
				tscale = atof(optarg);
				break;
			case 'c':
				col = atoi(optarg);
				break;
			case 'n':
				nbin = atoi(optarg);
				break;
			case 'o':
				ofile = optarg;
				break;
			case 'k':
				k1 = atof(optarg);
				break;
			case 'q':
				q1 = atof(optarg);
				break;
			case 'C':
				centre_data = true;
				break;
			case 'g':
				gh = true;
				break;
			default:
				fprintf(stderr,"?? getopt returned character %c ??\n", c);
				fprintf(stderr,"%s\n",Usage.c_str());
				exit(1);
		}
	}
	nfile = argc-optind;
	files.resize(nfile);
	for (int i=0; i<argc-optind; i++) 
		files[i] = argv[optind+i];
	fprintf(stdout,"==========================================================\n");
	fprintf(stdout,"Reading data from column %i of the following files:\n",col);
	for (int i=0; i<nfile; i++) 
		fprintf(stdout,"%s\n",files[i].c_str());
	fprintf(stdout,"Binning data from %8.3f to %8.3f in %i bins\n",lo,hi,nbin);
	fprintf(stdout,"Using Markov time step = %i\n",dt);
	fprintf(stdout,"Writing output to file %s\n",ofile.c_str());
	if (centre_data) {
		fprintf(stdout,"Will centre data in bins\n");
	} else {
		fprintf(stdout,"Will not centre data in bins\n");
	}
	fprintf(stdout,"==========================================================\n");

	for (int i=0; i<nfile; i++) 
		bin_transitions(files[i], hist, col, lo, hi, nbin, dt, centre_data);

	/*write_mat(ofile, lo, hi, nbin, dt, hist);*/
	if (gh) {
		write_sqmat(ofile,hist);
	} else {
		write_sqmat_extra(ofile,hist,lo,hi,nbin,float(dt)*tscale,k1,q1);
	}

	return 0;
}

