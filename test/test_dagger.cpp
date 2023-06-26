#define MAIN

#include "DMRG.hpp"
#include "service.hpp"

int main(){

    gsl_matrix_complex * H = gsl_matrix_complex_calloc(2,3);
    gsl_matrix_complex_set(H, 0, 0, gsl_complex_rect(1, 1));
    gsl_matrix_complex_set(H, 0, 1, gsl_complex_rect(2, 0));
    gsl_matrix_complex_set(H, 0, 2, gsl_complex_rect(2, 3));
    gsl_matrix_complex_set(H, 1, 0, gsl_complex_rect(0, 2));
    gsl_matrix_complex_set(H, 1, 1, gsl_complex_rect(6, 0.5));
    gsl_matrix_complex_set(H, 1, 1, gsl_complex_rect(111, 0.98));

    gsl_matrix_complex_print(H);
    cout << endl;
    gsl_matrix_complex_print(dagger(H));

    return 0;
}