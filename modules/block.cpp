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

void block::SetHamiltonian(gsl_matrix_complex * m)
{
    gsl_matrix_complex_free(H);
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
    SetHamiltonian(Htemp);
    gsl_matrix_complex_free(Htemp);

    // Redefine block S -> S*Id
    gsl_matrix_complex* temp;
    for(int i=0; i<l; i++)
    {
        for(size_t j=0; j<3; j++)
        {
            temp = S[i][j];
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

    // Add new site's S
    gsl_matrix_complex* Id = gsl_matrix_complex_alloc(chi, chi);
    gsl_matrix_complex_set_identity(Id);
    switch(pos) 
    {
        case 'l':
            for(size_t i=0; i<3; i++)
            {
                S.emplace_back(new gsl_matrix_complex*[3]);
                S[l][i] = gsl_matrix_complex_calloc(dim*chi, dim*chi);
                gsl_blas_zgetp(gsl_complex_rect(1,0), Id, M->getO(i), S[l][i]);
            } break;
        case 'r':
            for(size_t i=0; i<3; i++)
            {
                S.emplace_back(new gsl_matrix_complex*[3]);
                S[l][i] = gsl_matrix_complex_calloc(dim*chi, dim*chi);
                gsl_blas_zgetp(gsl_complex_rect(1,0), M->getO(i), Id, S[l][i]);
            } break;
        default: error_message("Block position not allowed in AddSite");
    }

    // Increase l
    l++;
}

// void block::EnlargeBlock(model* M, gsl_matrix_complex* ham, gsl_matrix_complex* r)
// {

//     // Update chi
//     if(M->getDim()*chi > chimax) chi = chimax;
//     else chi = chi*M->getDim();

//     // set H to renormalized ham (new hamiltonian for block+site)


//     // renormalize old S
//     // add new S for the new site
// }