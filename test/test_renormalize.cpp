#define MAIN

#include "DMRG.hpp"
#include "service.hpp"

int main(){

    std::pair<gsl_matrix_complex*, gsl_matrix_complex*> Rs;
    sys * SYS = new sys(1,2,3,4,2);

    SYS->getL()->AddSite();
    SYS->getR()->AddSite();

    cout << "########################" << endl;
    cout << "HLo =" << endl;
    gsl_matrix_complex_print(SYS->getL()->getH());
    cout << endl;
    cout << "HoR =" << endl;
    gsl_matrix_complex_print(SYS->getR()->getH());
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


    Rs = SYS->compute_Rmat();


    cout << "########################" << endl;
    cout << "Renormalization matrices:" << endl;
    cout << endl;
    
    
    cout << "RL = " << endl;
    gsl_matrix_complex_print(Rs.first);
    cout << "RR = " << endl;
    gsl_matrix_complex_print(Rs.second);
    cout << endl;

    ////////////////////////////////////////////

    cout << "########################" << endl;
    cout << "Enlarged spin operators:" << endl;
    cout << "Sx:" << endl;
    gsl_matrix_complex_print(SYS->getL()->getS(0, 0));
    cout << "Sy:" << endl;
    gsl_matrix_complex_print(SYS->getL()->getS(0, 1));
    cout << "Sz:" << endl;
    gsl_matrix_complex_print(SYS->getL()->getS(0, 2));

    ////////////////////////////////////////////

     cout << "########################" << endl;
     cout << "Renormalized Hamiltonians:" << endl;
     cout << endl;
     
     SYS->getL()->Renormalize(Rs.first);
     SYS->getR()->Renormalize(Rs.second);
     
     cout << "HL = " << endl;
     gsl_matrix_complex_print(SYS->getL()->getH());
     cout << "HR = " << endl;
     gsl_matrix_complex_print(SYS->getR()->getH());

    ////////////////////////////////////////////

    cout << "########################" << endl;
    cout << "Renormalized spin operators:" << endl;
    cout << "Sx:" << endl;
    gsl_matrix_complex_print(SYS->getL()->getS(0, 0));
    cout << "Sy:" << endl;
    gsl_matrix_complex_print(SYS->getL()->getS(0, 1));
    cout << "Sz:" << endl;
    gsl_matrix_complex_print(SYS->getL()->getS(0, 2));




    return 0;
}