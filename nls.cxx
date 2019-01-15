#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <complex>

//---------------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//---------------------------------------
void initialize(cmplx* const psi0, const double eta, const double dx,
                const double dt,const double xmin, const int N);



void writeToFile(cmplx* const psi, const string s, const double dx,
                 const int Nx, const double xmin);
//---------------------------------------
int main(){

        const int N = 1000;
        const double L = 200;
        const double xmin = -L/2;
        const double Tend = 100;
        const double dx = L / (N - 1);
        const double dt = dx  / 10;
        const int Na = 50;
        int Nk = int(Tend / Na / dt + 0.5);
        double t=0;

        const double eta = 0.2;

        cmplx* psi = new cmplx[N];

        stringstream strm;

        initialize(psi, eta, dx , dt,xmin, N);

        writeToFile(psi, "psi_0", dx, N, xmin);

        cout << "Nk = " << Nk << endl;

        for(int i=1; i<=Na; i++)
        {


            for(int j=1; j<=Nk-2; j++){


                t +=dt;
            }



            strm.str("");
            strm << "psi_" << i;
            writeToFile(psi, strm.str(), dx, N, xmin);
        }

        cout << "t = " << t << endl;

        delete[] psi;
        return 0;
}
//-----------------------------------------------
void initialize(cmplx* const psi0, const double eta, const double dx,
                const double dt,const double xmin, const int N)
{
  const double x0 = 0;
	const double psi = sqrt(2) * eta;

	for(int i=0; i<N; i++){
		double x = i*dx + xmin;
		psi0[i] = psi/cosh(0.5*eta * x);
	}
}
//-----------------------------------------------



//-----------------------------------------------


//-----------------------------------
void writeToFile(cmplx* const psi, const string s, const double dx,
                 const int Nx, const double xmin)
{
  ofstream out(s.c_str());

	for(int i=0; i<Nx; i++){
		double x = xmin + i * dx;
    out << x << "\t" << norm(psi[i])
             << "\t" << psi[i].real()
             << "\t" << psi[i].imag() <<  endl;
	}
	out.close();
}
