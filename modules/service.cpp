#define SERVICE_CPP

#include "DMRG.hpp"
#include "service.hpp"

// Print matrices

void gsl_matrix_complex_print(const gsl_matrix_complex *m)
{
	size_t i, j;
	const size_t N = m->size1;
	const size_t M = m->size2;

	for (i = 0; i < N; ++i){
		for (j = 0; j < M; ++j){
			gsl_complex z = gsl_matrix_complex_get(m, i, j);
			if(GSL_IMAG(z)<0) 	cout << GSL_REAL(z) << "-i" << -GSL_IMAG(z) << "\t";
			else				cout << GSL_REAL(z) << "+i" << GSL_IMAG(z) << "\t";
		}
		cout << endl;
	}
}

// Tensor product

int  gsl_blas_zgetp (const gsl_complex alpha,
                     const gsl_matrix_complex * A,
                     const gsl_matrix_complex * B,
                     gsl_matrix_complex * C){
	
	// sizes
	const size_t rA = A->size1;
	const size_t cA = A->size2;
	const size_t rB = B->size1;
	const size_t cB = B->size2;
	const size_t rC = C->size1;
	const size_t cC = C->size2;
	if(rC != rA*rB) if(cC != cA*cB) 
	{ error_message("C is "+to_string(rC)+"x"+to_string(cC)+", should be "+to_string(rA*rB)+"x"+to_string(cA*cB)); }

	// Set elements
	for (size_t i=0; i<rC ; ++i){
		for (size_t j=0; j<cC; ++j){
			gsl_matrix_complex_set(C, i, j, 
									gsl_complex_add(gsl_matrix_complex_get(C, i, j),
													gsl_complex_mul(alpha, 	
																	gsl_complex_mul(gsl_matrix_complex_get(A, i/rA, j/cA), 
																					gsl_matrix_complex_get(B, i%rB, j%cB))))
									);
		}
	}
	return 0;
}

// Reshaping

void res_vec_mat(const vector<double>& in, gsl_matrix_complex * out){
    size_t rows = out->size1;
    size_t cols = out->size2;
    for(size_t i=0; i<rows; i++){
        for(size_t j=0; j<cols; j++){
            gsl_matrix_complex_set(out, i, j, gsl_complex_rect(in[j + i*rows], 0));
        }
    }
}

void res_mat_vec(gsl_matrix_complex * in, vector<double>& out){
    size_t rows = in->size1;
    size_t cols = in->size2;
    for(size_t i=0; i<rows; i++){
        for(size_t j=0; j<cols; j++){
            out[j + i*rows] = GSL_REAL(gsl_matrix_complex_get(in, i, j));
        }
    }
}