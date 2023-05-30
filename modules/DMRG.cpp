#define DMRG_CPP

#include "DMRG.hpp"
#include "service.hpp"

DMRG::DMRG(double Jx_, double Jy_, double Jz_, double h_, int dim_){
    // construct sys
    S = new sys(Jx_, Jy_, Jz_, h_, dim_);
}

// void DMRG::run_DMRG(){
        // Infinite() until L
        
        // Remove last HmatsR from memory
        // HR.popback();

        // Finite until convergence (dE<eps)
        // Measure
// }

/****** Perform one step of the infinite algorithm ******/
void DMRG::Infinite()
{

    /* Add sites */
    gsl_matrix_complex * temp_HLo = S->getL()->AddSite();
    gsl_matrix_complex * temp_HoR = S->getR()->AddSite();
    
    /* Compute GS */
    Egs = S->compute_GS();

    /* Compute renormalization matrices */
    std::pair<gsl_matrix_complex*, gsl_matrix_complex*> R;
    R = S->compute_Rmat();

    /* Store R mats in memory vectors */
    RL.emplace_back(new gsl_matrix_complex);
    RR.emplace_back(new gsl_matrix_complex);
    RL[RL.size() -1] = R.first;
    RR[RR.size() -1] = R.second;

    /* Renormalize */
    S->getL()->Renormalize(R.first);
    S->getR()->Renormalize(R.second);

    /* Store H mats in memory vectors */
    HL.emplace_back(new gsl_matrix_complex);
    HR.emplace_back(new gsl_matrix_complex);
    HL[HL.size() -1] = temp_HLo;
    HR[HR.size() -1] = temp_HoR;

}

// void DMRG::Finite()
// {
//     // DIRECTION SWITCH


//     // Enlarge on left
//     S->getL()->AddSite();

//     // Shrink on right
//     S->getR()->BreakBlock(HR.back());
//     HR.pop_back();

//     // Reshape GS
//     reshape phi

//     // Compute new GS

//     // SVD

// }