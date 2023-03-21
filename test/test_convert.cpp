#define MAIN

#include "DMRG.hpp"
#include "service.hpp"

int main(){

    gsl_matrix * A = gsl_matrix_calloc(4, 4);
    gsl_matrix_set_identity(A);
    gsl_matrix_set(A, 1, 2, 0.5);

    gsl_matrix_print(A);

    gsl_matrix_complex * B = gsl_matrix_complex_calloc(4, 4);
    conv_real_comp(B, A);

    gsl_matrix_complex_print(B);

    gsl_matrix * C = gsl_matrix_calloc(4, 4);
    conv_comp_real(C, B);

    gsl_matrix_print(C);

    return 0;
}