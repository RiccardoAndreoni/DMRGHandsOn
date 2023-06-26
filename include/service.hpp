#ifndef SERVICE_HPP
#define SERVICE_HPP

#include "DMRG.hpp"

/**************** Print matrices ***************/

void gsl_matrix_complex_print(const gsl_matrix_complex *m);
void gsl_matrix_print(const gsl_matrix *m);
void gsl_vector_print(const gsl_vector *v);



/**************** Tensor product ***************/

int  gsl_blas_zgetp (const gsl_complex alpha,
					const gsl_matrix_complex * A,
					const gsl_matrix_complex * B,
					gsl_matrix_complex * C);


/******************** Dagger *******************/

gsl_matrix_complex* dagger(gsl_matrix_complex* m);

/********************* Norm ********************/

double gsl_matrix_complex_norm(gsl_matrix_complex* m);
void gsl_matrix_complex_normalize(gsl_matrix_complex* m);

/**************** Error messages ***************/

inline void error_message(const string err)
{
	cout << "ERROR: " << err << endl;
	exit(1);
}

/******************* Reshape *******************/
void res_vec_mat(const vector<double>& in, gsl_matrix_complex * out);
void res_vec_mat(const vector<double>& in, gsl_matrix * out);
void res_mat_vec(gsl_matrix_complex * in, vector<double>& out);

/******************* Convert *******************/
void conv_real_comp(gsl_matrix_complex * out, gsl_matrix * in);
void conv_comp_real(gsl_matrix * out, gsl_matrix_complex * in);



/************** Lanczos utilities **************/

/**
 * @brief Template class to implement random vector initializer.
 *
 * "Partially specialization of function" is not allowed,
 * so here it is mimicked by wrapping the "init" function with a class template.
*/
template <typename T>
struct InitialGuess {
public:
/**
 * @brief Initialize given vector randomly in the range of [-1, 1].
 *
 * For complex type, the real and imaginary part of each element will be initialized in
 * the range of [-1, 1].
*/
static void init(std::vector<T>& v) {

	size_t n = v.size();
	for(size_t i = 0; i < n; ++i) {
	v[i] = 0;
	}
}
};

template <typename T>
struct InitialGuess<std::complex<T>> {
public:
static void init(std::vector<std::complex<T>>& v) {
	
	size_t n = v.size();
	for(size_t i = 0; i < n; ++i) {
	v[i] = 0;
	}
}
};


#endif