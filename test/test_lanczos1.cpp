#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <lambda_lanczos/lambda_lanczos.hpp>
#include "DMRG.hpp"
#include "service.hpp"

using std::cout;
using std::endl;
using std::setprecision;
using lambda_lanczos::LambdaLanczos;
template<typename T>
using vector = std::vector<T>;


int main() {
    const int n = 4;
    // gsl_matrix_complex * matrix = gsl_matrix_complex_calloc(3, 3);
    
    // gsl_matrix_complex_set(matrix, 0, 0, gsl_complex_rect(2, 0));
    // gsl_matrix_complex_set(matrix, 0, 1, gsl_complex_rect(1, 0));
    // gsl_matrix_complex_set(matrix, 0, 2, gsl_complex_rect(1, 0));
    // gsl_matrix_complex_set(matrix, 1, 0, gsl_complex_rect(1, 0));
    // gsl_matrix_complex_set(matrix, 1, 1, gsl_complex_rect(2, 0));
    // gsl_matrix_complex_set(matrix, 1, 2, gsl_complex_rect(1, 0));
    // gsl_matrix_complex_set(matrix, 2, 0, gsl_complex_rect(1, 0));
    // gsl_matrix_complex_set(matrix, 2, 1, gsl_complex_rect(1, 0));
    // gsl_matrix_complex_set(matrix, 2, 2, gsl_complex_rect(3, 0));

    // double matrix[n][n] = { {2.0, 1.0, 2.0},
    //                         {1.0, 2.0, 1.0},
    //                         {2.0, 1.0, 2.0} };
    /* Its eigenvalues are {4, 1, 1} */

    // the matrix-vector multiplication routine

    gsl_matrix_complex * Id = gsl_matrix_complex_calloc(2,2);
    gsl_matrix_complex_set_identity(Id);

    gsl_matrix_complex * H = gsl_matrix_complex_calloc(2,2);
    gsl_matrix_complex_set(H, 0, 0, gsl_complex_rect(1, 0));
    gsl_matrix_complex_set(H, 0, 1, gsl_complex_rect(2, 0));
    gsl_matrix_complex_set(H, 1, 0, gsl_complex_rect(2, 0));
    gsl_matrix_complex_set(H, 1, 1, gsl_complex_rect(1, 0));

    gsl_matrix_complex * S = gsl_matrix_complex_calloc(2,2);
    gsl_matrix_complex_set(S, 0, 1, gsl_complex_rect(1, 0));
    gsl_matrix_complex_set(S, 1, 0, gsl_complex_rect(1, 0));

    auto mv_mul = [&](const std::vector<double>& in, std::vector<double>& out) {
        gsl_matrix_complex * in_mat = gsl_matrix_complex_calloc(2, 2);
        gsl_matrix_complex * out_mat = gsl_matrix_complex_calloc(2, 2);

        res_vec_mat(in, in_mat);
        
        gsl_matrix_complex * temp = gsl_matrix_complex_calloc(2,2);

        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), S, in_mat, gsl_complex_rect(0,0), temp);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), temp, S, gsl_complex_rect(0,0), out_mat);

        gsl_matrix_complex_free(temp);

        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), H, in_mat, gsl_complex_rect(1,0), out_mat);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), in_mat, H, gsl_complex_rect(1,0), out_mat);

        res_mat_vec(out_mat, out);
    };

    LambdaLanczos<double> engine(mv_mul, n, false, 1); // true means to calculate the largest eigenvalue.
    std::vector<double> eigenvalues;
    std::vector<std::vector<double>> eigenvectors;
    engine.run(eigenvalues, eigenvectors);

    cout << "Eigenvalue: " << setprecision(16) << eigenvalues[0] << endl;
    cout << "Eigenvector: ";
    for(int i = 0; i < n; ++i) {
    cout << eigenvectors[0][i] << " ";
    }
    cout << endl;

    return EXIT_SUCCESS;
}
