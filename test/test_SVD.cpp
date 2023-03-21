#define MAIN

#include "DMRG.hpp"
#include "service.hpp"

int main(){

    sys * SYS = new sys(1,2,3,4,2);

    // cout << "h = "  << GSL_REAL(SYS->M->geth()) << endl;    
    // cout << "Jx = " << GSL_REAL(SYS->M->getJ(0)) << endl;
    // cout << "Jy = " << GSL_REAL(SYS->M->getJ(1)) << endl;
    // cout << "Jz = " << GSL_REAL(SYS->M->getJ(2)) << endl;    

    SYS->compute_GS();
    SYS->compute_Rmat();
    
    return 0;
}