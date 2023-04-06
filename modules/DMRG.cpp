#define DMRG_CPP

#include "DMRG.hpp"
#include "service.hpp"

 DMRG::DMRG(double Jx_, double Jy_, double Jz_, double h_, int dim_){
     // construct sys
    S = new sys(Jx_, Jy_, Jz_, h_, dim_);
     // allocate RL, RR (?)

 }

// void DMRG::run_DMRG(){
//     // Infinite() until L
//     // Finite until convergence (dE<eps)
//     // Measure
// }

void DMRG::Infinite(){



    // Add sites
    S->getL()->AddSite();
    S->getR()->AddSite();
    

    // Compute GS
    Egs = S->compute_GS();

    // Compute renormalization matrices
    std::pair<gsl_matrix_complex*, gsl_matrix_complex*> R;
    R = S->compute_Rmat();

    // Store R mats in memory vectors
    RL.emplace_back(new gsl_matrix_complex);
    RR.emplace_back(new gsl_matrix_complex);
    RL[RL.size() -1] = R.first;
    RR[RR.size() -1] = R.second;

    // Renormalize
    S->getL()->Renormalize(R.first);
    S->getR()->Renormalize(R.second);

}

// DMRG::Finite()
// {

// }