#define MAIN

#include "DMRG.hpp"
#include "service.hpp"

int main(){

    vector<double> A = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    gsl_matrix * B = gsl_matrix_calloc(3, 3);

    res_vec_mat(A, B);

    gsl_matrix_print(B);
    cout << endl;
    cout << "\t|" << endl;
    cout << "\tV" << endl;
    cout << endl;
    
    gsl_matrix_complex * C = gsl_matrix_complex_calloc(2, 2);
    gsl_matrix_view D = gsl_matrix_submatrix(B, 0, 0, 2, 2);
    conv_real_comp(C, &D.matrix);

    gsl_matrix_complex_print(C);

    return 0;
}