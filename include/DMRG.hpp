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
#include<gsl/gsl_linalg.h>
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

		gsl_matrix_complex** O;		//Single site spin matrices
		gsl_matrix_complex* Id;		//Identity
		void InitOps();				//Initialize ssOps depending on spin 

		gsl_matrix_complex* Hs;		//Single site Hamiltonian
		//gsl_matrix_complex* Hint;	//(NN) interaction Hamiltonian
		void InitH();

	public:
		model(double Jx_, double Jy_, double Jz_, double h_, int dim_);
		// model(string parfile);

		/* Get parameters */
		int getDim(){ return dim; }								//Get dim
		gsl_complex getJ(size_t i){ return J[i]; }				//Get J[i]
		gsl_complex geth(){ return h; }				//Get J[i]
		gsl_matrix_complex * getO(size_t i){ return O[i]; }		//Get O[i]
		gsl_matrix_complex * getId(){ return Id; }	
		gsl_matrix_complex * getHs(){ return Hs; }				//Get H single site
		//gsl_matrix_complex * getHint(){ return Hint; }			//Get H single site
};

class block
{
	private:
  		
  		int l;										//Size of the block
  		int chi; 									//Number of the eigenvectors: min{chimax, d^l}

  		gsl_matrix_complex * H;						//Hamiltonian
  		void SetHamiltonian(complex<double>* M);  	//Import the renormalized H (thilde) from system

		gsl_matrix_complex*** S;					//Renormalized single site operators acting on the whole block space

  	public: 

  		block(model *M);							//Default constructor as single site block

  		double* GetHamiltonian();
  		void EnlargeBlock(double* M);				//Substitute H with the new renormalized Hamiltonian

		// Get parameters
		int getChi() { return chi; }
		int getl() { return l; }
		gsl_matrix_complex * getS(size_t site, size_t a){ return S[site][a]; } 
		gsl_matrix_complex * getH(){ return H; }	//Get O[i]
};

class sys
{
	private: 

		model* M;
		block* L;
		block* R;
		gsl_matrix_complex * GS;		// Ground state of the system

	public:

		sys(double Jx_, double Jy_, double Jz_, double h_, int dim_);

		// Add sites
		void computeHLo(gsl_matrix_complex * H);	// Add site to the left block to build HLo from HL
		void computeHoR(gsl_matrix_complex * H);	// Add site to the left block to build HoR from HR

		// GS computation
		double compute_GS();

		// Miao
		void compute_Rmat();
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