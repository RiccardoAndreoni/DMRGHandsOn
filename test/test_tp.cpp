#define MAIN

#include "DMRG.hpp"
#include "service.hpp"

int main(){

    gsl_matrix_complex * A = gsl_matrix_complex_calloc(2,2);
    gsl_matrix_complex_set(A, 0, 0, gsl_complex_rect(1, 0));
    gsl_matrix_complex_set(A, 0, 1, gsl_complex_rect(0, 1));
    gsl_matrix_complex_set(A, 1, 0, gsl_complex_rect(1, 2));
    gsl_matrix_complex_set(A, 1, 1, gsl_complex_rect(2, 0));
    gsl_matrix_complex * B = gsl_matrix_complex_calloc(2,2);
    gsl_matrix_complex_set_all(B, gsl_complex_rect(1, 0));
    gsl_matrix_complex * C = gsl_matrix_complex_calloc(4,4);
    gsl_matrix_complex_set_identity(C);

    gsl_matrix_complex_print(A);
    cout << endl;
    gsl_matrix_complex_print(B);
    cout << endl;
    gsl_matrix_complex_print(C);
    cout << endl;

    gsl_blas_zgetp(gsl_complex_rect(1, 0), A, B, C);

    gsl_matrix_complex_print(C);


    return 0;
}