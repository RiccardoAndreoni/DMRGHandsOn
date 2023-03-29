#define DMRG_CPP

#include "DMRG.hpp"
#include "service.hpp"

// DMRG::DRMG(){
//     // construct model
//     // construct sys
//     // allocate RL, RR
// }

// void DMRG::run_DMRG(){
//     // Infinite() until L
//     // Finite until convergence (dE<eps)
//     // Measure
// }

void DMRG::Infinite(){
    // Add sites
    
    // Compute GS

    // Compute renormalization matrices
    std::pair<gsl_matrix_complex*, gsl_matrix_complex*> R;
    R = S->compute_Rmat();

    // Store R mats in memory vectors
    RL.emplace_back(new gsl_matrix_complex*);
    RR.emplace_back(new gsl_matrix_complex*);
    RL[RL.size() -1] = R.first;
    RR[RR.size() -1] = R.second;

    // Renormalize

}

// DMRG::Infinite()
// {

// }