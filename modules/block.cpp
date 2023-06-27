#define BLOCK_CPP

#include "DMRG.hpp"
#include "service.hpp"

block::block(model *M, char p)
{
	l = 1;
	chi = M->getDim();
	pos = p;

	/* Same model in sys */
	this->M = M;

	/* Initialize single site Hamiltonian */
	H = gsl_matrix_complex_calloc(M->getDim(),M->getDim());
	gsl_matrix_complex_memcpy(H, M->getHs());

	/* Initialize single site spin operators */
	/////////////////////////////////////////////////
	S.emplace_back(new gsl_matrix_complex* [3]);
	for(size_t i=0; i<3; i++)
	{ 
		S[0][i] = gsl_matrix_complex_calloc(M->getDim(),M->getDim());
		gsl_matrix_complex_memcpy(S[0][i], M->getO(i)); 
	}
}

/* Replace block->H with a new Hamiltonian m */
void block::SetHamiltonian(gsl_matrix_complex * m, bool freemem)
{
		cout << "\tIn SetHamiltonian: " << H << endl; // TEST
		cout << "\t\t -> dim " << H->size1 << "x" << H->size2 << endl; // TEST
	if(freemem) gsl_matrix_complex_free(H); // false to keep memory about the discarded Hamiltonian
	H = m;
}

/********* Routines to add site to the block *********/

/* Build HLo from HL when a site is added on the RIGHT of the block */
void block::computeHLo(gsl_matrix_complex * m)
{
	gsl_matrix_complex * Id = gsl_matrix_complex_calloc(chi, chi);
	gsl_matrix_complex_set_identity(Id);

	gsl_blas_zgetp(gsl_complex_rect(1, 0), H, M->getId(), m);
	gsl_blas_zgetp(gsl_complex_rect(1, 0), Id , M->getHs(), m);

	for(size_t i=0; i<3; i++)
	{ gsl_blas_zgetp(M->getJ(i), S.back()[i], M->getO(i), m); }
}

/* Build HoR from HR when a site is added on the LEFT of the block */
void block::computeHoR(gsl_matrix_complex * m)
{
	gsl_matrix_complex * Id = gsl_matrix_complex_calloc(chi, chi);
	gsl_matrix_complex_set_identity(Id);
	gsl_blas_zgetp(gsl_complex_rect(1, 0), M->getId(), H, m);
	gsl_blas_zgetp(gsl_complex_rect(1, 0), M->getHs(), Id , m);
	
	for(size_t i=0; i<3; i++)
	{ gsl_blas_zgetp(M->getJ(i), M->getO(i), S.back()[i], m); }
}

/* Add site to a block */
gsl_matrix_complex * block::AddSite()
{
	int dim = M->getDim();

	// cout << "\t 011 " << H << endl; // TEST
	
	cout << "AddSite.1" << endl; // TEST

	/* Compute Hamiltonian block+site */
	gsl_matrix_complex* Htemp = gsl_matrix_complex_calloc(chi*dim, chi*dim);
	// cout << "\t 01a " << H << " " << Htemp << endl; // TEST

	switch(pos) 
	{
		case 'l':
			computeHLo(Htemp);
		break;
		case 'r':
			computeHoR(Htemp);
		break;
		default:
			error_message("Block position not allowed in AddSite");
	}
	SetHamiltonian(Htemp, true);   // we don't need info about the renormalized Hamiltonians

	// cout << "\t 012 " << H << endl; // TEST

	cout << "AddSite.2" << endl; // TEST

	/* Redefine block S -> S*Id */
	gsl_matrix_complex* temp;
	for(int i=0; i<l; i++)
	{
		for(size_t j=0; j<3; j++)
		{
			// cout << j << ":" << S[i][j]->size1 << "'" << S[i][j]->size2 << endl; // TEST
			temp = S[i][j];
			// cout << j << ":" << temp->size1 << "'" << temp->size2 << endl; // TEST
			S[i][j] = gsl_matrix_complex_calloc(chi*dim, chi*dim);
			switch(pos) 
			{
				case 'l':
					gsl_blas_zgetp(gsl_complex_rect(1,0), temp, M->getId(), S[i][j]);
				break;
				case 'r':
					gsl_blas_zgetp(gsl_complex_rect(1,0), M->getId(), temp, S[i][j]);
				break;
				default: error_message("Block position not allowed in AddSite");
			}
			gsl_matrix_complex_free(temp);
		}
	}

	// cout << "\t 013 " << H << endl; // TEST

	cout << "AddSite.3" << endl; // TEST

	/* Add new site's S */
	gsl_matrix_complex* Id = gsl_matrix_complex_alloc(chi, chi);
	gsl_matrix_complex_set_identity(Id);
	S.emplace_back(new gsl_matrix_complex*[3]);
	switch(pos) 
	{
		case 'l':
			for(size_t i=0; i<3; i++)
			{
				S.back()[i] = gsl_matrix_complex_calloc(dim*chi, dim*chi);
				gsl_blas_zgetp(gsl_complex_rect(1,0), Id, M->getO(i), S.back()[i]);
			} break;
		case 'r':
			for(size_t i=0; i<3; i++)
			{
				S.back()[i] = gsl_matrix_complex_calloc(dim*chi, dim*chi);
				gsl_blas_zgetp(gsl_complex_rect(1,0), M->getO(i), Id, S.back()[i]);
			} break;
		default: error_message("Block position not allowed in AddSite");
	}

	// /* Increase l */
	// l++;
	
	// cout << "\t 014 " << H << endl; // TEST

	return Htemp; 
}

/* Break the block for the Finite */ //TO TEST
void block::BreakBlock(gsl_matrix_complex * old_H)
{
	// Change l
	l -= 1;

	// Change back chi
	int dim = M->getDim();

		cout << "\t l=" << l << ", dim=" << dim << endl;
		
	if((int)pow(dim,l) > chimax) chi = chimax; 
	else chi = (int)pow(dim,l);

		cout << "\t chi = " << chi << endl; // TEST

	// Change Hamiltonian
		cout << "in BreakBlock: " << old_H << endl; // TEST
	SetHamiltonian(old_H, true);
}

/* Transform all block(l)+site operators in block(l+1) operators */
void block::Renormalize(gsl_matrix_complex* R)
{
	/* Update chi */
	int dim = M->getDim();
	int old_chi = chi;
	if(dim*old_chi > chimax) chi = chimax;
	else chi = old_chi*dim;

	/* Set H to renormalized ham (new hamiltonian for block+site) */
	gsl_matrix_complex * Hnew = gsl_matrix_complex_calloc(chi, chi);
	gsl_matrix_complex * temp = gsl_matrix_complex_calloc(chi, old_chi*dim);
	gsl_blas_zgemm(CblasTrans, CblasNoTrans, gsl_complex_rect(1,0), R, H, gsl_complex_rect(0,0), temp);
	gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), temp, R, gsl_complex_rect(0,0), Hnew);
	gsl_matrix_complex_free(temp);
	SetHamiltonian(Hnew, false); 	// DON'T FREE MEMORY ABT THE HAMILTONIAN: 
									// we need to store the block+site Ham for the finite part
	
	/* Renormalize S single site */
	gsl_matrix_complex * temp2 = gsl_matrix_complex_calloc(chi, old_chi*dim);
	for(int i=0; i<l+1; i++){
		for(int j=0; j<3; j++)
		{
			gsl_blas_zgemm(CblasTrans, CblasNoTrans, gsl_complex_rect(1,0), R, S[i][j], gsl_complex_rect(0,0), temp2);
			gsl_matrix_complex_free(S[i][j]);
			S[i][j] = gsl_matrix_complex_calloc(chi, chi);
			gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), temp2, R, gsl_complex_rect(0,0), S[i][j]);
		}
	}
	// for(int j=0; j<3; j++)
	// {
	// 	gsl_blas_zgemm(CblasTrans, CblasNoTrans, gsl_complex_rect(1,0), R, S.back()[j], gsl_complex_rect(0,0), temp2);
	// 	gsl_matrix_complex_free(S.back()[j]);
	// 	S.back()[j] = gsl_matrix_complex_calloc(chi, chi);
	// 	gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), temp2, R, gsl_complex_rect(0,0), S.back()[j]);
	// }
	gsl_matrix_complex_free(temp2);

	/* Increase l */
	l++;
}

void block::delS()
{
	for(size_t i=0; i<3; i++)
	{
		gsl_matrix_complex_free(S.back()[i]);
	}
	free(S.back());
	S.pop_back();
}