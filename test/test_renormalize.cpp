#define MAIN

#include "DMRG.hpp"
#include "service.hpp"

// Move the blocks L and R to public in sys to run this code

int main(){

    std::pair<gsl_matrix_complex*, gsl_matrix_complex*> Rs;
    sys * SYS = new sys(1,2,3,4,2);

    SYS->L->AddSite();
    SYS->R->AddSite();

    cout << "########################" << endl;
    cout << "HLo =" << endl;
    gsl_matrix_complex_print(SYS->L->getH());
    cout << endl;
    cout << "HoR =" << endl;
    gsl_matrix_complex_print(SYS->R->getH());
    cout << endl;

    ////////////////////////////////////////////

    double En;
    En = SYS->compute_GS();

    cout << "########################" << endl;
    cout << "Ground state = " << endl;
    gsl_matrix_complex_print(SYS->getGS());
    cout << endl;
    cout << "Ground state energy = " << En << endl;

    ////////////////////////////////////////////

    cout << "########################" << endl;
    cout << "Renormalization matrices:" << endl;
    cout << endl;
    
    Rs = SYS->compute_Rmat();
    
    cout << "RL = " << endl;
    gsl_matrix_complex_print(Rs.first);
    cout << "RR = " << endl;
    gsl_matrix_complex_print(Rs.second);
    cout << endl;

    ////////////////////////////////////////////

    cout << "########################" << endl;
    cout << "Renormalized Hamiltonians:" << endl;
    cout << endl;

    SYS->L->Renormalize(Rs.first);
    SYS->R->Renormalize(Rs.second);
    
    cout << "HL = " << endl;
    gsl_matrix_complex_print(SYS->L->getH());
    cout << "HR = " << endl;
    gsl_matrix_complex_print(SYS->R->getH());

    return 0;
}