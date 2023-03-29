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

// Lanczos initial guess

template<typename T>
using vector = std::vector<T>;
/**
 * @brief Template class to implement random vector initializer.
 *
 * "Partially specialization of function" is not allowed,
 * so here it is mimicked by wrapping the "init" function with a class template.
 */
template <typename T>
struct InitialGuessZero {
	public:
	/**
	 * @brief Initialize given vector randomly in the range of [-1, 1].
	 *
	 * For complex type, the real and imaginary part of each element will be initialized in
	 * the range of [-1, 1].
	 */
	static void init(std::vector<T>& v) {

		// size_t n = v.size();
		// for(size_t i = 0; i < n; ++i) {
		// v[i] = 1;
		// }
		v[0] = 0;
		v[1] = 0;
		v[2] = 1;
		v[3] = 1;
	}
};

template <typename T>
struct InitialGuessZero<std::complex<T>> {
	public:
	static void init(std::vector<std::complex<T>>& v) {

		// size_t n = v.size();
		// for(size_t i = 0; i < n; ++i) {
		//   v[i] = 1;
		// }
		v[0] = 0;
		v[1] = 0;
		v[2] = 1;
		v[3] = 1;
	}
};

///////////////////////////////////////////////////////////////



int main() {
	const int n = 4;

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
	engine.init_vector = InitialGuessZero<double>::init;
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
