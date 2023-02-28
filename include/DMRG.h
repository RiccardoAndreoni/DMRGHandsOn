#include<iostream>
#include<cmath>
#include<fstream>
#include<iomanip>
#include<complex>
#define dim 3
#define h 2.
#define Jx 1.
#define Jy 1.
#define Jz 1.
#define chimax = 100

using namespace std;

class block
{
	private:
  		
  		int l;										//Size of the block
  		int chi; 									//Number of the eigenvectors: min{chimax, d^l}
  		complex<double>* H;							//Hamiltonian
  		complex<double>** O;						//Single site operators {}	
		void InitHamiltonian();						//Initialize H from operators
  		void SetHamiltonian(complex<double>* M);  	//Import the renormalized H (thilde) from system

  	public: 

  		double* GetHamiltonian();
  		void EnlargeBlock(double* M);
  		block();							//Default constructor as single site block

};

class sys
{
	private: 

		block S,E,s1,s2;
		double* Hss;						//Site-site interaction Hamiltonian
		// double* Htot;
		void InitHss();

	public:

		sys();
		// void SetHtot();						//GetHamiltonian -> builds H_L+2 (superblock)
		double* GetHtot();

};

class DMRG
{

	private: 

		sys S;
		double** R;

	public: 

		void Infinite(); //Compute gs via Lanczos, reduced density matrix, chi highest eigenvectors, R, Hnew->Hthildenew			
		void Sweeps();
		DMRG();
};


