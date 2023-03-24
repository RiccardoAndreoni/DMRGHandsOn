#define MAIN

#include "DMRG.hpp"
#include "service.hpp"

int main(){

    // BEFORE
    gsl_matrix * A = gsl_matrix_calloc(4, 4);
    gsl_matrix_set_identity(A);
    gsl_matrix_set(A, 1, 2, 0.5);

    cout << "### BEFORE ###" << endl;
    cout << endl;

    cout << "A =" << endl;
    gsl_matrix_print(A);


    // AFTER
    gsl_matrix* temp = A;

    A = gsl_matrix_calloc(2,3);

    cout << "### BEFORE ###" << endl;
    cout << endl;
    
    cout << "A =" << endl;
    gsl_matrix_print(A);
    cout << "temp =" << endl;
    gsl_matrix_print(temp);

    return 0;
}