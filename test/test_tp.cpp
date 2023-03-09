#define MAIN

#include "DMRG.hpp"
#include "service.hpp"

int main(){

    // gsl_matrix_complex * A = gsl_matrix_complex_calloc(2,2);
    // gsl_matrix_complex_set(A, 0, 0, gsl_complex_rect(1, 0));
    // gsl_matrix_complex_set(A, 0, 1, gsl_complex_rect(0, 1));
    // gsl_matrix_complex_set(A, 1, 0, gsl_complex_rect(1, 2));
    // gsl_matrix_complex_set(A, 1, 1, gsl_complex_rect(2, 0));
    // gsl_matrix_complex * B = gsl_matrix_complex_calloc(2,2);
    // gsl_matrix_complex_set_all(B, gsl_complex_rect(1, 0));
    // gsl_matrix_complex * C = gsl_matrix_complex_calloc(4,4);
    // gsl_matrix_complex_set_identity(C);

    gsl_matrix_complex * Id = gsl_matrix_complex_calloc(2,2);
    gsl_matrix_complex_set_identity(Id);

    gsl_matrix_complex * H = gsl_matrix_complex_calloc(2,2);
    gsl_matrix_complex_set(H, 0, 0, gsl_complex_rect(1, 0));
    gsl_matrix_complex_set(H, 0, 1, gsl_complex_rect(2, 0));
    gsl_matrix_complex_set(H, 1, 0, gsl_complex_rect(2, 0));
    gsl_matrix_complex_set(H, 1, 1, gsl_complex_rect(1, 0));

    gsl_matrix_complex * Result = gsl_matrix_complex_calloc(4, 4);

    gsl_blas_zgetp(gsl_complex_rect(1, 0), H, Id, Result);
    gsl_blas_zgetp(gsl_complex_rect(1, 0), Id, H, Result);

    gsl_matrix_complex_print(Id);
    cout << endl;
    gsl_matrix_complex_print(H);
    cout << endl;
    gsl_matrix_complex_print(Result);
    cout << endl;

    // gsl_blas_zgetp(gsl_complex_rect(1, 0), A, B, C);

    // gsl_matrix_complex_print(C);


    return 0;
}