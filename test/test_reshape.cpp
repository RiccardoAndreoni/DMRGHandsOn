#define MAIN

#include "DMRG.hpp"
#include "service.hpp"

int main(){

    vector<double> A = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    gsl_matrix_complex * B = gsl_matrix_complex_calloc(3, 3);
    vector<double> C = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    cout << "\tA : " << endl;
    for(size_t i=0; i<size(A); i++){
        cout << A[i] << " ";
    }
    cout << endl;
    cout << "\t|" << endl;
    cout << "\tV" << endl;
    cout << endl;

    res_vec_mat(A, B);

    gsl_matrix_complex_print(B);
    cout << endl;
    cout << "\t|" << endl;
    cout << "\tV" << endl;
    cout << endl;
    
    res_mat_vec(B, C);

    for(size_t i=0; i<size(C); i++){
        cout << C[i] << ", ";
    }
    cout << endl;



    return 0;
}