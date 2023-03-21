#define BLOCK_CPP

#include "DMRG.hpp"
#include "service.hpp"

block::block(model *M){
	l = 1;
    chi = M->getDim();

	// Initialize single site Hamiltonian
	H = gsl_matrix_complex_calloc(M->getDim(),M->getDim());
    gsl_matrix_complex_memcpy(H, M->getHs());

    // Initialize single site spin operators
    S = new gsl_matrix_complex** [l]; 
    S[0] = new gsl_matrix_complex* [3];
    for(size_t i=0; i<3; i++)
    { 
        S[0][i] = gsl_matrix_complex_calloc(M->getDim(),M->getDim());
        gsl_matrix_complex_memcpy(S[0][i], M->getO(i)); 
    }
}