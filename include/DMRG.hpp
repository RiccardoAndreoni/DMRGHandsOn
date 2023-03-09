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
#include <lambda_lanczos/lambda_lanczos.hpp>
#define chimax = 100

using namespace std;

// Classes

class model
{
	private:
		gsl_complex* J;
		gsl_complex h;
		int dim;

		gsl_matrix_complex** O;		//Single site operators
		void InitOps();				//Initialize ssOps depending on spin 

		gsl_matrix_complex* Hs;		//Single site Hamiltonian
		gsl_matrix_complex* Hint;	//(NN) interaction Hamiltonian
		void InitH();

	public:
		model(double Jx_, double Jy_, double Jz_, double h_, int dim_);
		model(string parfile);

		/* Get parameters */
		int getDim(){ return dim; }							//Get dim
		gsl_complex getJ(size_t i){ return J[i]; }			//Get J[i]
		gsl_matrix_complex * getO(size_t i){ return O[i]; }	//Get O[i]
		gsl_matrix_complex * getHs(){ return Hs; }			//Get H single site
};

class block
{
	private:
  		
  		int l;										//Size of the block
  		int chi; 									//Number of the eigenvectors: min{chimax, d^l}

  		gsl_matrix_complex * H;						//Hamiltonian
		void InitHamiltonian();						//Initialize H from operators
  		void SetHamiltonian(complex<double>* M);  	//Import the renormalized H (thilde) from system

		gsl_matrix_complex*** S;					//Renormalized single site operators acting on the whole block space

  	public: 

  		block(model &M);							//Default constructor as single site block

  		double* GetHamiltonian();
  		void EnlargeBlock(double* M);				//Substitute H with the new renormalized Hamiltonian

		// Get parameters
		int getChi() { return chi; }
		int getl() { return l; }
		gsl_matrix_complex * getS(int site, int a){ return S[site][a]; }
		gsl_matrix_complex * getH(){ return H; }	//Get O[i]

};

class sys
{
	private: 

		model* M;
		block* L;
		block* R;
		double* Hss;				//Site-site interaction Hamiltonian
		gsl_matrix_complex * GS;	// Ground state of the system
		void InitHss();

	public:

		sys(double Jx_, double Jy_, double Jz_, double h_, int dim_);

		// GS computation
		void compute_GS();
		// auto mv_mul = [&](const std::vector<double>& in, std::vector<double>& out);

		// Add sites
		void add_site_L(gsl_matrix_complex * H);	// Add site to the left block to build HLo from HL
		void add_site_R(gsl_matrix_complex * H);	// Add site to the left block to build HoR from HR

		void compute_Rmat();
		// double* GetHtot();

};

class DMRG
{

	private: 

		sys S;
		std::vector<gsl_matrix_complex*> RL;
		std::vector<gsl_matrix_complex*> RR;

	public: 

		void Infinite(); //Compute gs via Lanczos, reduced density matrix, chi highest eigenvectors, R, Hnew->Hthildenew			
		void Sweeps();
		DMRG();
};


#endif