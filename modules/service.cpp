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

void gsl_matrix_print(const gsl_matrix *m)
{
	size_t i, j;
	const size_t N = m->size1;
	const size_t M = m->size2;

	for (i = 0; i < N; ++i){
		for (j = 0; j < M; ++j){
			cout << gsl_matrix_get(m, i, j) << "\t";
		}
		cout << endl;
	}
}

void gsl_vector_print(const gsl_vector *v)
{
	size_t i;
	const size_t N = v->size;

	for (i = 0; i < N; ++i){
		cout << gsl_vector_get(v, i) << endl;
	}
}

/***** Tensor product *****/

int  gsl_blas_zgetp (const gsl_complex alpha,
                     const gsl_matrix_complex * A,
                     const gsl_matrix_complex * B,
                     gsl_matrix_complex * C)
{
	/* sizes */
	const size_t rA = A->size1;
	const size_t cA = A->size2;
	const size_t rB = B->size1;
	const size_t cB = B->size2;
	const size_t rC = C->size1;
	const size_t cC = C->size2;
	if(rC != rA*rB) if(cC != cA*cB) 
	{ error_message("C is "+to_string(rC)+"x"+to_string(cC)+", should be "+to_string(rA)+"*"+to_string(rB)+"x"+to_string(cA)+"*"+to_string(cB)); }

	/* Set elements */
	for (size_t i=0; i<rC ; ++i){
		for (size_t j=0; j<cC; ++j){
			gsl_matrix_complex_set(C, i, j, 
									gsl_complex_add(gsl_matrix_complex_get(C, i, j),
													gsl_complex_mul(alpha, 	
																	gsl_complex_mul(gsl_matrix_complex_get(A, i/rB, j/cB), 
																					gsl_matrix_complex_get(B, i%rB, j%cB))))
									);
		}
	}
	return 0;
}

/***** Dagger *****/

gsl_matrix_complex* dagger(gsl_matrix_complex* m) 
{
	size_t rows = m->size1;
	size_t cols = m->size2;

    gsl_matrix_complex* result = gsl_matrix_complex_alloc(cols, rows);

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            gsl_matrix_complex_set(result, j, i, gsl_complex_conjugate(gsl_matrix_complex_get(m, i, j)));
        }
    }

    return result;
}

/***** Norm *****/

double gsl_matrix_complex_norm(gsl_matrix_complex* m)
{
	gsl_matrix_complex* temp = gsl_matrix_complex_alloc(m->size1, m->size2);
	gsl_matrix* temp2		 = gsl_matrix_calloc(m->size1, m->size2);

	gsl_matrix_complex_memcpy(temp, m);
	gsl_matrix_complex_mul_elements(temp, m);

	conv_comp_real(temp2, temp);

	double norm = sqrt(gsl_matrix_norm1(temp2));

	gsl_matrix_complex_free(temp);
	gsl_matrix_free(temp2);

	return norm;
}

void gsl_matrix_complex_normalize(gsl_matrix_complex* m)
{
	double norm = gsl_matrix_complex_norm(m);
	cout << norm << endl;
	gsl_matrix_complex_scale(m, gsl_complex_rect(1./norm, 0));
}

/***** Reshaping *****/

void res_vec_mat(const vector<double>& in, gsl_matrix_complex * out)
{
    size_t rows = out->size1;
    size_t cols = out->size2;
    for(size_t i=0; i<rows; i++){
        for(size_t j=0; j<cols; j++){
            gsl_matrix_complex_set(out, i, j, gsl_complex_rect(in[j + i*rows], 0));
        }
    }
}

void res_vec_mat(const vector<double>& in, gsl_matrix * out)
{
    size_t rows = out->size1;
    size_t cols = out->size2;
    for(size_t i=0; i<rows; i++){
        for(size_t j=0; j<cols; j++){
            gsl_matrix_set(out, i, j, in[j + i*rows]);
        }
    }
}

void res_mat_vec(gsl_matrix_complex * in, vector<double>& out)
{
    size_t rows = in->size1;
    size_t cols = in->size2;
    for(size_t i=0; i<rows; i++){
        for(size_t j=0; j<cols; j++){
            out[j + i*rows] = GSL_REAL(gsl_matrix_complex_get(in, i, j));
        }
    }
}

/***** Convert *****/

void conv_real_comp(gsl_matrix_complex * out, gsl_matrix * in)
{
	if(in->size1 != out->size1) error_message("in and out must be the same size");
	if(in->size2 != out->size2) error_message("in and out must be the same size");

	size_t i, j;
	const size_t N = in->size1;
	const size_t M = in->size2;

	for (i = 0; i < N; ++i)
	{
		for (j = 0; j < M; ++j)
		{
			gsl_matrix_complex_set(out, i, j, gsl_complex_rect(gsl_matrix_get(in, i, j), 0));
		}
	}
}

void conv_comp_real(gsl_matrix * out, gsl_matrix_complex * in)
{
	if(in->size1 != out->size1) error_message("in and out must be the same size");
	if(in->size2 != out->size2) error_message("in and out must be the same size");

	size_t i, j;
	const size_t N = in->size1;
	const size_t M = in->size2;

	for (i = 0; i < N; ++i)
	{
		for (j = 0; j < M; ++j)
		{
			gsl_matrix_set(out, i, j, GSL_REAL(gsl_matrix_complex_get(in, i, j)));
		}
	}
}