#include "fftw3.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>
//-------------------------------
using namespace std;
//-------------------------------
double readData(double* v, const int Nx, const char* const fname);
void writeData(const fftw_complex* const f, const int N, const double T,const char* const fname);
//-------------------------------
int main(int argc, char** argv){

	if(argc != 3){
		cout << "Usage: " << argv[0] << " input_file \t output_file" << endl;
		exit(1);
	}

	char *in_file  = argv[1];
	char *out_file = argv[2];

	const int N = 16384;
	double T;

	// Allocate memory
	fftw_complex* f = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2+1));
	double* inR  = (double*) malloc(sizeof(double)*N);

  // Create plan
	fftw_plan FW  = fftw_plan_dft_r2c_1d(N, inR, f, FFTW_ESTIMATE);

	// Read input data
  // - Read data from file, store second column in input field for forward
	//   FFT. At the same time determine duration T of the signal.

  T = readData(inR, N, in_file);

  // Calculate FFT
  fftw_execute(FW);

  // write output file
  writeData(f, N,  T, out_file);

  // Clean up
	fftw_destroy_plan(FW);
  fftw_free(f);
  free(inR);

	return 0;
}
//-------------------------------
void writeData(const fftw_complex* const f, const int N, const double T,const char* const fname){
	ofstream out(fname);
	const double dOm = 2*M_PI/T;
	double pk;

	for(int i=0; i<=N/2; i++){
		pk = sqrt(f[i][0]*f[i][0] + f[i][1]*f[i][1])/N;
		out << i*dOm << "\t" << pk << "\t" << f[i][0] << "\t" << f[i][1] << endl;
	}
	out.close();
}
//-------------------------------
double readData(double* v, const int Nx,  const char* const fname){
	ifstream in(fname);
	double x,value;

	for(int i=1; i<=Nx; i++){
		in >> x;
		in >> value;
		v[i-1] = value;
	}

	in.close();
	return x;
}
