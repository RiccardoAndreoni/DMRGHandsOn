#ifndef DMRG_HPP
#define DMRG_HPP

#include<iostream>
#include<cmath>
#include<fstream>
#include<iomanip>
#include<string>
#include<vector>
#include<complex.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_matrix_complex_double.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_complex_math.h>
#define chimax = 100

using namespace std;

// Classes

class model
{
	private:
		gsl_complex Jx, Jy, Jz, h;
		int dim;

		gsl_matrix_complex** O;	//Single site operators
		void InitOps();			//Initialize ssOps depending on spin 

		gsl_matrix_complex* Hs;	//Single site Hamiltonian
		gsl_matrix_complex* Hint;	//(NN) interaction Hamiltonian
		void InitH();

	public:
		model(double Jx_, double Jy_, double Jz_, double h_, int dim_);
		model(string parfile);

		int D(); 	// get dim
		gsl_complex JX();	// get Jx
		gsl_complex JY();	// get Jy
		gsl_complex JZ();	// get Jz
};

class block
{
	private:
  		
  		int l;										//Size of the block
  		int chi; 									//Number of the eigenvectors: min{chimax, d^l}
  		complex<double>* H;							//Hamiltonian
		void InitHamiltonian();						//Initialize H from operators
  		void SetHamiltonian(complex<double>* M);  	//Import the renormalized H (thilde) from system

  	public: 

  		double* GetHamiltonian();
  		void EnlargeBlock(double* M);
  		block(model &M);							//Default constructor as single site block

};

class sys
{
	private: 

		model M;
		block L, R;
		double* Hss;		//Site-site interaction Hamiltonian
		double* GS;
		void InitHss();

	public:

		sys();

		// GS computation
		void compute_GS();
		// auto mv_mul = [&](const vector<double>& in, vector<double>& out);
		void Ham_mul(gsl_matrix_complex * in, gsl_matrix_complex * out);

		void compute_Rmat();
		// double* GetHtot();

};

class DMRG
{

	private: 

		sys S;
		double** RL;
		double** RR;

	public: 

		void Infinite(); //Compute gs via Lanczos, reduced density matrix, chi highest eigenvectors, R, Hnew->Hthildenew			
		void Sweeps();
		DMRG();
};


#endif