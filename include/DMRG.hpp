#ifndef DMRG_HPP
#define DMRG_HPP

#include<iostream>
#include<cmath>
#include<fstream>
#include<iomanip>
#include<string>
#include<vector>
#include <utility>
#include<complex.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_matrix_complex_double.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_complex_math.h>
#include <lambda_lanczos/lambda_lanczos.hpp>

#define chimax 6

using namespace std;

/******************* Classes *******************/

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
  		
  		int l;		//Size of the block
  		int chi; 	//Number of the eigenvectors: min{chimax, d^l}
		char pos;	//Position of the block: 'l' left, 'r' right
		model* M;

  		gsl_matrix_complex * H;							//Hamiltonian
  		void SetHamiltonian(gsl_matrix_complex* m);  	//Subsitute H with m

		std::vector<gsl_matrix_complex**> S;	//Renormalized single site operators acting on the whole block space
				//For block L -> S[i] refers to the ith site
				//For block R -> S is mirrored so S[0]:lth site, S[l-1]:0th site and so on

  	public: 

  		block(model *M, char p);	//Default constructor as single site block

		/* Renormalize block */
  		void Renormalize(gsl_matrix_complex* R);

		/* Add sites */
		void computeHLo(gsl_matrix_complex* H);	// Build HLo from HL
		void computeHoR(gsl_matrix_complex* H);	// Build HoR from HR
		void AddSite();

		/* Get stuff */
		int getChi() { return chi; }
		int getl() { return l; }
		gsl_matrix_complex * getS(size_t site, size_t a){ return S[site][a]; } 
		gsl_matrix_complex * getH(){ return H; }	
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

		// GS computation
		double compute_GS();

		// Computation of renormalization matrices
		std::pair<gsl_matrix_complex*, gsl_matrix_complex*> compute_Rmat();

		// Get stuff
		gsl_matrix_complex * getGS(){ return GS; }
		block* getL(){ return L; }
		block* getR(){ return R; }
};


class DMRG
{

	private: 

		sys *S;
		std::vector<gsl_matrix_complex*> RL;
		std::vector<gsl_matrix_complex*> RR;
		double Egs;

	public: 

		DMRG(double Jx_, double Jy_, double Jz_, double h_, int dim_);

		/* Infinite */
		void Infinite(); //Compute gs via Lanczos, reduced density matrix, chi highest eigenvectors, R, Hnew->Hthildenew	
		void LoadRmats(pair<gsl_matrix_complex*, gsl_matrix_complex*>);	

		/* Finite */
		void Sweeps();

		/*Return Egs*/
		double getEgs(){return Egs;}
		sys* getSYS(){ return S; }
};


#endif