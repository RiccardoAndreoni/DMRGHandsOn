#define SYS_CPP

#include "DMRG.hpp"
#include "service.hpp"

sys::sys(double Jx_, double Jy_, double Jz_, double h_, int dim_)
{
	M = new model(Jx_, Jy_, Jz_, h_, dim_);
	// cout << "model built" << endl; // TEST
	L = new block(M, 'l');
	// cout << "L built" << endl; // TEST
	R = new block(M, 'r');
	// cout << "R built" << endl; // TEST
}


/* Lanczos initial guess */

// template<typename T>
// using vector = std::vector<T>;
// /**
//  * @brief Template class to implement random vector initializer.
//  *
//  * "Partially specialization of function" is not allowed,
//  * so here it is mimicked by wrapping the "init" function with a class template.
//  */

// template <typename T>
// struct InitialGuessLanczos{
// 	public:
// 	static void init(std::vector<T>& v) {
		
// 		// We need access to RR and RL vectors of renormalization matrices which belongs to the DMRG class
// 			//Options:  -Take init as an input somehow and do the GS_tilde computation in dmrg? This implies 
// 			//			   creating a diff compute GS for the finite algorithm
// 			// 			-Call RR and RL somehow ?
//  		GS_tilde=transpose(U(most recent)).GS.(Id(2) x transpose(V(old infinite iterations)))
// 		GS_tilde=reshape(GS_tilde)
// 		RR 
// 		v=GS_tilde

// 		// res_mat_vec(GS)
		
// 		// int RandSeed = 13;
//     	// std::mt19937 mt(RandSeed);
//     	// std::uniform_real_distribution<T> rand((T)(-1.0), (T)(1.0));

//    		// size_t n = v.size();
//     	// for(size_t i = 0; i < n; ++i) {
//       	// 	v[i] = rand(mt); 
//     	// }
// 	}
// };

// template <typename T>
// struct InitialGuessLanczos<std::complex<T>> {
// 	public:
// 	static void init(std::vector<std::complex<T>>& v) {

// 		// res_mat_vec(GS)


// 		// int RandSeed = 13;
// 	    // std::mt19937 mt(RandSeed);
// 	    // std::uniform_real_distribution<T> rand((T)(-1.0), (T)(1.0));

// 	    // size_t n = v.size();
// 	    // for(size_t i = 0; i < n; ++i) {
// 	    // 	v[i] = std::complex<T>(rand(mt), rand(mt));

// 		// }
// 	}
// };

/****** GS computation ******/
double sys::compute_GS()
{
	/* Set dimensions */
	int chi_l = L->getChi();
	int chi_r = R->getChi();
	int dim = M->getDim();
	int n = dim*dim*chi_l*chi_r;

	/* Define mv_mul */
	auto mv_mul = [&](const std::vector<double>& in, std::vector<double>& out) 
	{
		gsl_matrix_complex * in_mat = gsl_matrix_complex_calloc(chi_l*dim, dim*chi_r);
		gsl_matrix_complex * out_mat = gsl_matrix_complex_calloc(chi_l*dim, dim*chi_r);
		
		/* Reshape in: vector -> gsl_matrix */
		res_vec_mat(in, in_mat);
		
		/* Sl phi Sr */
		gsl_matrix_complex * S_l = gsl_matrix_complex_calloc(chi_l*dim, dim*chi_l);
		gsl_matrix_complex * S_r = gsl_matrix_complex_calloc(chi_r*dim, dim*chi_r);

		gsl_matrix_complex * Id_l = gsl_matrix_complex_calloc(chi_l, chi_l);
		gsl_matrix_complex_set_identity(Id_l);
		gsl_matrix_complex * Id_r = gsl_matrix_complex_calloc(chi_r, chi_r);
		gsl_matrix_complex_set_identity(Id_r);

		gsl_matrix_complex * temp = gsl_matrix_complex_calloc(chi_l*dim, dim*chi_r);
		for(size_t i=0; i<3; i++){ 
			gsl_blas_zgetp(gsl_complex_rect(1,0), Id_l, M->getO(i), S_l);
			gsl_blas_zgetp(gsl_complex_rect(1,0), M->getO(i), Id_r, S_r);
			gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, M->getJ(i), S_l, in_mat, gsl_complex_rect(0,0), temp);
			gsl_blas_zgemm(CblasNoTrans, CblasTrans, gsl_complex_rect(1,0), temp, S_r, gsl_complex_rect(1,0), out_mat);
			gsl_matrix_complex_set_zero(temp);
			gsl_matrix_complex_set_zero(S_l);
			gsl_matrix_complex_set_zero(S_r);
		}

		gsl_matrix_complex_free(temp);
		gsl_matrix_complex_free(Id_l);
		gsl_matrix_complex_free(Id_r);
		gsl_matrix_complex_free(S_l);
		gsl_matrix_complex_free(S_r);

		/* Hl phi 1r */
		gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), L->getH(), in_mat, gsl_complex_rect(1,0), out_mat);
		/* 1l phi Hr */
		gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), in_mat, R->getH(), gsl_complex_rect(1,0), out_mat);  
		
		/* Reshape out: gsl_matrix -> vector */
		res_mat_vec(out_mat, out);

	};

	// std::vector<double> test_in{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};  // TEST
	// std::vector<double> test_in{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};            // TEST
	// std::vector<double> test_in{0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};            // TEST
	// std::vector<double> test_out{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};           // TEST
	// cout << "BEFORE:" << endl;                                                              // TEST
	// for(size_t i=0; i<size(test_out); i++) {cout << test_out[i] << endl;}                   // TEST
	// cout << endl;                                                                           // TEST
	// mv_mul(test_in, test_out);                                                              // TEST
	// cout << "AFTER: " << endl;                                                              // TEST
	// for(size_t i=0; i<size(test_out); i++) {cout << test_out[i] << endl;}                   // TEST
	// cout << endl;                                                                           // TEST

	/* Allocate GS */
	gsl_matrix_complex_free(GS);
	GS = gsl_matrix_complex_calloc(chi_l*dim, dim*chi_r);

	/* Lanczos */
	lambda_lanczos::LambdaLanczos<double> engine(mv_mul, n, false, 1); 
	// if(guess) engine.init_vector = InitialGuessLanczos<double>::init; // Initial guess taken from the stored GS

	std::vector<double> eigenvalues;
	std::vector<std::vector<double>> eigenvectors;

	// cout << "Run Lanczos" << endl;  // TEST
	engine.run(eigenvalues, eigenvectors);
	// cout << "End Lanczos" << endl;  // TEST
	
	res_vec_mat(eigenvectors[0], GS); // NB: GS is real

	// cout << endl;                                                                           // TEST
	// for(size_t i=0; i<size(eigenvectors[0]); i++) {cout << eigenvectors[0][i] << endl;}     // TEST
	// cout << "E0 = " << eigenvalues[0] << endl;                                              // TEST

	return eigenvalues[0];
}

/****** SVD and resize -> define RL and RR ******/
std::pair<gsl_matrix_complex*, gsl_matrix_complex*> sys::compute_Rmat()
{
	gsl_matrix_complex* tempRL;
	gsl_matrix_complex* tempRR;

	/* Set dimensions */
	int chi_l = L->getChi();
	int chi_r = R->getChi();
	int dim = M->getDim();
	// int n = dim*dim*chi_l*chi_r;

	/* Definition of temporary U and V, to be truncated after SVD */
	gsl_matrix * tempU = gsl_matrix_alloc(chi_l*dim, dim*chi_r);
	conv_comp_real(tempU, GS);
	gsl_matrix * tempV = gsl_matrix_alloc(dim*chi_r, dim*chi_r);
	gsl_vector * tempS = gsl_vector_alloc(dim*chi_r);
	gsl_vector * tempW = gsl_vector_alloc(dim*chi_r);

	gsl_linalg_SV_decomp(tempU, tempV, tempS, tempW);
	
	cout << "########################" << endl; //TEST
	cout << "Singular values = " << endl; // TEST
	gsl_vector_print(tempS);            // TEST
	cout << endl;                        //TEST

	gsl_vector_free(tempS);
	gsl_vector_free(tempW);

	// cout << "U = " << endl;     // TEST
	// gsl_matrix_print(tempU);    // TEST
	// cout << endl;               // TEST
	// cout << "V = " << endl;     // TEST
	// gsl_matrix_print(tempV);    // TEST
	// cout << endl;               // TEST
	// cout << "S = " << endl;     // TEST
	// gsl_vector_print(tempS);    // TEST

	/* Truncation */
	if(tempU->size2 > chimax)
	{
		tempRL = gsl_matrix_complex_alloc(tempU->size1, chimax);
		gsl_matrix_view truncU = gsl_matrix_submatrix(tempU, 0, 0, tempU->size1, chimax);
		conv_real_comp(tempRL, &truncU.matrix); 
	}
	else 
	{
		tempRL = gsl_matrix_complex_alloc(tempU->size1, tempU->size2);
		conv_real_comp(tempRL, tempU); 
	}
	if(tempV->size2 > chimax)
	{ 
		tempRR = gsl_matrix_complex_alloc(tempV->size1, chimax);
		gsl_matrix_view truncV = gsl_matrix_submatrix(tempV, 0, 0, tempV->size1, chimax);
		conv_real_comp(tempRR, &truncV.matrix); 
	}
	else 
	{
		tempRR = gsl_matrix_complex_alloc(tempV->size1, tempV->size2);
		conv_real_comp(tempRR, tempV);
	}

	return make_pair(tempRL, tempRR);
}

/***** Save *****/

void sys::saveH(char pos, gsl_matrix_complex* H)
{
	switch(pos) 
	{
		case 'l':
			HL.emplace_back(new gsl_matrix_complex);
			HL.back() = H;
		break;
		case 'r':
			HR.emplace_back(new gsl_matrix_complex);
			HR.back() = H;
		break;
		default:
			error_message("Block position not allowed in saveH");
	}
}

void sys::saveR(char pos, gsl_matrix_complex* R)
{
	switch(pos) 
	{
		case 'l':
			RL.emplace_back(new gsl_matrix_complex);
			RL.back() = R;
		break;
		case 'r':
			RR.emplace_back(new gsl_matrix_complex);
			RR.back() = R;
		break;
		default:
			error_message("Block position not allowed in saveR");
	}
}

/***** Load *****/

gsl_matrix_complex* sys::loadH(char pos)
{
	switch(pos) 
	{
		case 'l':
			return HL.back();
			HL.pop_back();
		break;
		case 'r':
			return HR.back();
			HR.pop_back();
		break;
		default:
			error_message("Block position not allowed in loadH");
			exit(1);
	}
}

gsl_matrix_complex* sys::loadR(char pos)
{
	switch(pos) 
	{
		case 'l':
			return RL.back();
			RL.pop_back();
		break;
		case 'r':
			return RR.back();
			RR.pop_back();
		break;
		default:
			error_message("Block position not allowed in loadR");
			exit(1);
	}
}
