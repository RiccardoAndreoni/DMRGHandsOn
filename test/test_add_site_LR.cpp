#define MAIN

#include "DMRG.hpp"
#include "service.hpp"

// Move the blocks L and R to public in sys to run this code

int main(){
    sys * SYS = new sys(1,2,3,4,2);
    size_t dir = 1;

    SYS->L->AddSite();
    SYS->R->AddSite();
    
    cout << "L : S[0][" << dir << "] = " << endl;
    gsl_matrix_complex_print(SYS->L->getS(0, dir));

    cout << "R : S[0][" << dir << "] = " << endl;
    gsl_matrix_complex_print(SYS->R->getS(0, dir));

    cout << "L : S[1][" << dir << "] = " << endl;
    gsl_matrix_complex_print(SYS->L->getS(1, dir));

    cout << "R : S[1][" << dir << "] = " << endl;
    gsl_matrix_complex_print(SYS->R->getS(1, dir));
    
    return 0;
}