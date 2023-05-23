#define BLOCK_CPP

#include "DMRG.hpp"
#include "service.hpp"

block::block(model *M, char p){
	l = 1;
    chi = M->getDim();
    pos = p;

    // Same model in sys
    this->M = M;

	// Initialize single site Hamiltonian
	H = gsl_matrix_complex_calloc(M->getDim(),M->getDim());
    gsl_matrix_complex_memcpy(H, M->getHs());

    // Initialize single site spin operators
    /////////////////////////////////////////////////
    S.emplace_back(new gsl_matrix_complex* [3]);
    for(size_t i=0; i<3; i++)
    { 
        S[0][i] = gsl_matrix_complex_calloc(M->getDim(),M->getDim());
        gsl_matrix_complex_memcpy(S[0][i], M->getO(i)); 
    }
}

void block::SetHamiltonian(gsl_matrix_complex * m, bool freemem)
{
    if(freemem) gsl_matrix_complex_free(H);
    H = m;
}

/********* Routines to add site to the block *********/

// Add site to the left block to build HLo from HL
void block::computeHLo(gsl_matrix_complex * m)
{
    gsl_matrix_complex * Id = gsl_matrix_complex_calloc(chi, chi);
    gsl_matrix_complex_set_identity(Id);

    gsl_blas_zgetp(gsl_complex_rect(1, 0), H, M->getId(), m);
    gsl_blas_zgetp(gsl_complex_rect(1, 0), Id , M->getHs(), m);

    for(size_t i=0; i<3; i++)
    { gsl_blas_zgetp(M->getJ(i), S[l-1][i], M->getO(i), m); }
}

// Add site to the left block to build HoR from HR 
 void block::computeHoR(gsl_matrix_complex * m)
 {
     gsl_matrix_complex * Id = gsl_matrix_complex_calloc(chi, chi);
     gsl_matrix_complex_set_identity(Id);
     gsl_blas_zgetp(gsl_complex_rect(1, 0), M->getId(), H, m);
     gsl_blas_zgetp(gsl_complex_rect(1, 0), M->getHs(), Id , m);
     
     for(size_t i=0; i<3; i++)
     { gsl_blas_zgetp(M->getJ(i), M->getO(i), S[l-1][i], m); }
 }

void block::AddSite()
{
    int dim = M->getDim();
    
   // cout << "1" << endl;
    // Compute Hamiltonian block+site
    gsl_matrix_complex* Htemp = gsl_matrix_complex_calloc(chi*dim, chi*dim);
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
    SetHamiltonian(Htemp, false);   // DON'T FREE MEMORY ABT THE HAMILTONIAN: 
                                    // we need to store the renormalized one for the finite part

    //cout << "2" << endl;
    // Redefine block S -> S*Id
    // gsl_matrix_complex_print(S[0][1]);
    // cout << S[0][1]->size1 << "'" << S[0][1]->size2 << endl;

    gsl_matrix_complex* temp;
    for(int i=0; i<l; i++)
    {
        for(size_t j=0; j<3; j++)
        {
            // cout << j << ":" << S[i][j]->size1 << "'" << S[i][j]->size2 << endl; TEST
            temp = S[i][j];
            // cout << j << ":" << temp->size1 << "'" << temp->size2 << endl; TEST
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

    //cout << "3" << endl;
    // Add new site's S
    gsl_matrix_complex* Id = gsl_matrix_complex_alloc(chi, chi);
    gsl_matrix_complex_set_identity(Id);
    S.emplace_back(new gsl_matrix_complex*[3]);
    switch(pos) 
    {
        case 'l':
            for(size_t i=0; i<3; i++)
            {
                S[l][i] = gsl_matrix_complex_calloc(dim*chi, dim*chi);
                gsl_blas_zgetp(gsl_complex_rect(1,0), Id, M->getO(i), S[l][i]);
            } break;
        case 'r':
            for(size_t i=0; i<3; i++)
            {
                S[l][i] = gsl_matrix_complex_calloc(dim*chi, dim*chi);
                gsl_blas_zgetp(gsl_complex_rect(1,0), M->getO(i), Id, S[l][i]);
            } break;
        default: error_message("Block position not allowed in AddSite");
    }

    // Increase l
    l++;
}

void block::Renormalize(gsl_matrix_complex* R)
{
    // Update chi
    int dim = M->getDim();
    int old_chi = chi;
    if(dim*old_chi > chimax) chi = chimax;
    else chi = old_chi*dim;

    // Set H to renormalized ham (new hamiltonian for block+site)
    gsl_matrix_complex * Hnew = gsl_matrix_complex_calloc(chi, chi);
    gsl_matrix_complex * temp = gsl_matrix_complex_calloc(chi, old_chi*dim);
    gsl_blas_zgemm(CblasTrans, CblasNoTrans, gsl_complex_rect(1,0), R, H, gsl_complex_rect(0,0), temp);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), temp, R, gsl_complex_rect(0,0), Hnew);
    gsl_matrix_complex_free(temp);
    SetHamiltonian(Hnew, true); // we don't need info about the nonrenormalized Ham
    
    // Renormalize S single site
    gsl_matrix_complex * temp2 = gsl_matrix_complex_calloc(chi, old_chi*dim);
    for(int i=0; i<l; i++){
        for(int j=0; j<3; j++)
        {
            gsl_blas_zgemm(CblasTrans, CblasNoTrans, gsl_complex_rect(1,0), R, S[i][j], gsl_complex_rect(0,0), temp2);
            gsl_matrix_complex_free(S[i][j]);
            S[i][j] = gsl_matrix_complex_calloc(chi, chi);
            gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), temp2, R, gsl_complex_rect(0,0), S[i][j]);
        }
    }
    gsl_matrix_complex_free(temp2);
}