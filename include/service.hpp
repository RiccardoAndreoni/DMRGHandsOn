#ifndef SERVICE_HPP
#define SERVICE_HPP

#include "DMRG.hpp"

// Print matrices

void gsl_matrix_complex_print(const gsl_matrix_complex *m);
void gsl_matrix_print(const gsl_matrix *m);

// Tensor product

int  gsl_blas_zgetp (const gsl_complex alpha,
                     const gsl_matrix_complex * A,
                     const gsl_matrix_complex * B,
                     gsl_matrix_complex * C);

// Error messages

inline void error_message(const string err)
{
	cout << "ERROR: " << err << endl;
	exit(1);
}

// Reshape
void res_vec_mat(const vector<double>& in, gsl_matrix_complex * out);
void res_vec_mat(const vector<double>& in, gsl_matrix * out);
void res_mat_vec(gsl_matrix_complex * in, vector<double>& out);

// Convert
void conv_real_comp(gsl_matrix_complex * out, gsl_matrix * in);
void conv_comp_real(gsl_matrix * out, gsl_matrix_complex * in);


#endif