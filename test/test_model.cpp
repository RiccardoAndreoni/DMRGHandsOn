#define MAIN

#include "DMRG.hpp"
#include "service.hpp"

int main(){

    model * M = new model(1,2,3,4,2);

    cout << "h = "  << GSL_REAL(M->geth()) << endl;    
    cout << "Jx = " << GSL_REAL(M->getJ(0)) << endl;
    cout << "Jy = " << GSL_REAL(M->getJ(1)) << endl;
    cout << "Jz = " << GSL_REAL(M->getJ(2)) << endl;    

    cout << "Id = " << endl;
    gsl_matrix_complex_print(M->getId());
    cout << "Sx = " << endl;
    gsl_matrix_complex_print(M->getO(0));
    cout << "Sy = " << endl;
    gsl_matrix_complex_print(M->getO(1));
    cout << "Sz = " << endl;
    gsl_matrix_complex_print(M->getO(2));
    cout << endl; 

    cout << "Hs = " << endl;
    gsl_matrix_complex_print(M->getHs());
    cout << "Hint = " << endl;
    gsl_matrix_complex_print(M->getHint());


    return 0;
}