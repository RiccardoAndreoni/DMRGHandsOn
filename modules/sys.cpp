#define SYS_CPP

#include "DMRG.hpp"
#include "service.hpp"

// sys::sys(){
//     L  = block();
//     R  = block();
// }


/* GS computation */

void sys::compute_GS(){
    // define mv_mul
    // lanczos
}

// auto sys::mv_mul = [&](const vector<double>& in, vector<double>& out){

//     int chi_l = L->chi;
//     int chi_r = R->chi;
//     int dim = M->dim;

//     gsl_matrix_complex * in_mat = gsl_matrix_complex_calloc(chi_l*dim, dim*chi_r);
//     gsl_matrix_complex * out_mat = gsl_matrix_complex_calloc(chi_l*dim, dim*chi_r);

//     // Reshape in: vector -> gsl_matrix 
    
    
//     // Define Htot action
//     // Reshape out: gsl_matrix -> vector 
// }

void sys::Ham_mul(gsl_matrix_complex * in_mat, gsl_matrix_complex * out_mat){
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1, 0), )
}
