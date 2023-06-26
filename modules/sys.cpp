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

/****** compute guess from old GS ******/

void sys::GStoGuess(char dir)
{
	int dim = M->getDim();
	int chi_l = L->getChi();
	int chi_r = R->getChi();

	gsl_matrix_complex * newGS;

	gsl_matrix_complex * U; // Renormalization matrices
	gsl_matrix_complex * V;

	gsl_matrix_complex * Lprod; // Submatrices for bases change product
	gsl_matrix_complex * Rprod;
	gsl_matrix_complex * temp;
	gsl_matrix_complex * id = gsl_matrix_complex_calloc(dim,dim);
	gsl_matrix_complex_set_identity(id);


    cout << "GStoGUESS.1" << endl; // TEST

	switch(dir) 
	{
		case 'r':	// U*GS*(1xV)
			U = loadR('l', false);
			V = loadR('r', false);
			cout << "U(" << U->size1 << "x" << U->size2 << ")" << endl;
			cout << "V(" << V->size1 << "x" << V->size2 << ")" << endl;
			Lprod = dagger(U);

			Rprod = gsl_matrix_complex_calloc(V->size2*dim, V->size1*dim); // exchange dimensions wrt V
			gsl_matrix_complex * temp_dag;
			temp_dag = dagger(V);
			gsl_blas_zgetp(gsl_complex_rect(1,0), id, temp_dag, Rprod);

			gsl_matrix_complex_free(temp_dag);
			// gsl_matrix_complex_free(V); // Free memory about the right renormalization matrix (popping out of RR is not enough)

			newGS = gsl_matrix_complex_calloc(chi_l, dim*dim*chi_r);
		break;
		case 'l':	// (1xU)*GS*V
			U = loadR('l', false);
			V = loadR('r', false);

			Lprod = gsl_matrix_complex_calloc(U->size1*dim, U->size2*dim); 
			gsl_blas_zgetp(gsl_complex_rect(1,0), U, id, Lprod);

			Rprod = V;

			// gsl_matrix_complex_free(U); 

			newGS = gsl_matrix_complex_calloc(chi_l*dim*dim, chi_r);			
		break;
		default:
			error_message("Direction not allowed in GStoGuess");
	}

	cout << "GStoGUESS.2" << endl; // TEST

	temp = gsl_matrix_complex_calloc(Lprod->size1, GS->size2);

	cout << "\tchi_l: " << chi_l << " chi_r: " << chi_r << endl; // TEST

	cout << "\tLprod: " << Lprod->size1 << "x" << Lprod->size2 << endl; // TEST
	cout << "\tGS: "    << GS->size1    << "x" << GS->size2    << endl; // TEST
	cout << "\tRprod: " << Rprod->size1 << "x" << Rprod->size2 << endl; // TEST
	cout << "\tnewGS: " << newGS->size1 << "x" << newGS->size2 << endl; // TEST
	cout << "GStoGUESS.3" << endl; // TEST
	gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), Lprod, GS, gsl_complex_rect(0,0), temp);
	cout << "GStoGUESS.4" << endl; // TEST
			// cout << "\t\ttemp(" << temp->size1 << "x" << temp->size2 << ")  -  Rprod(" 
			// << Rprod->size1 << "x" << Rprod->size2 << ")" << endl; // TEST
	gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), temp, Rprod, gsl_complex_rect(0,0), newGS);

	cout << "GStoGUESS.5" << endl; // TEST

	/* Free memory */
	gsl_matrix_complex_free(temp);
	gsl_matrix_complex_free(Lprod);
	gsl_matrix_complex_free(Rprod);

	cout << "GStoGUESS.6" << endl; // TEST

	gsl_matrix_complex_normalize(newGS);
	GUESS = newGS;
}

/****** GS computation ******/

double sys::compute_GS(bool guess)
{
	/* Set dimensions */
	int chi_l = L->getChi();
	int chi_r = R->getChi();
	int dim = M->getDim();
	int n = dim*dim*chi_l*chi_r;

	gsl_matrix_complex * tmp_HL = gsl_matrix_complex_alloc(L->getH()->size1,L->getH()->size2);
	gsl_matrix_complex_memcpy(tmp_HL, L->getH());
	gsl_matrix_complex * tmp_HR = gsl_matrix_complex_alloc(R->getH()->size1,R->getH()->size2);
	gsl_matrix_complex_memcpy(tmp_HR, R->getH());

	cout << "compute_GS.1" << endl; // TEST
	cout << "\tchi_l = " << chi_l << "; chi_r = " << chi_r << endl; // TEST

	/* Define mv_mul */
	auto mv_mul = [&](const std::vector<double>& in, std::vector<double>& out) 
	{
		// cout << "  Entering mv_mul..." << endl;	// TEST

		gsl_matrix_complex * in_mat = gsl_matrix_complex_calloc(chi_l*dim, dim*chi_r);
		gsl_matrix_complex * out_mat = gsl_matrix_complex_calloc(chi_l*dim, dim*chi_r);
		
		/* Reshape in: vector -> gsl_matrix */
		// cout << "\t reshape vec->mat" << endl; // TEST
		res_vec_mat(in, in_mat);


		/* Sl phi Sr */
		// cout << "\t Spin definitions" << endl; // TEST
		gsl_matrix_complex * S_l = gsl_matrix_complex_calloc(chi_l*dim, dim*chi_l);
		gsl_matrix_complex * S_r = gsl_matrix_complex_calloc(chi_r*dim, dim*chi_r);

		gsl_matrix_complex * Id_l = gsl_matrix_complex_calloc(chi_l, chi_l);
		gsl_matrix_complex_set_identity(Id_l);
		gsl_matrix_complex * Id_r = gsl_matrix_complex_calloc(chi_r, chi_r);
		gsl_matrix_complex_set_identity(Id_r);

		// cout << "\t interaction part" << endl; // TEST
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

		// cout << "\t free" << endl; // TEST
		gsl_matrix_complex_free(temp);
		gsl_matrix_complex_free(Id_l);
		gsl_matrix_complex_free(Id_r);
		gsl_matrix_complex_free(S_l);
		gsl_matrix_complex_free(S_r);

		// cout << "\t tensor product" << endl; // TEST
		// cout << "\t\tHl(" << tmp_HL->size1 << "x" << tmp_HL->size2 << ")  -  in_mat(";	// TEST
		// cout << in_mat->size1 << "x" << in_mat->size2 << ")   ##   " ;							// TEST
		// cout << "\tin_mat(" << in_mat->size1 << "x" << in_mat->size2 << ")  -  Hr("; 			// TEST
		// cout << tmp_HR->size1 << "x" << tmp_HR->size2 << ")" << endl; 					// TEST

		/* Hl phi 1r */
		// cout << "\t\t1." << endl; // TEST
		gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), tmp_HL, in_mat, gsl_complex_rect(1,0), out_mat);
		/* 1l phi Hr */
		// cout << "\t\t2." << endl; // TEST
		gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), in_mat, tmp_HR, gsl_complex_rect(1,0), out_mat);  
		

		/* Reshape out: gsl_matrix -> vector */
		// cout << "\t  reshape mat->vec" << endl; // TEST
		res_mat_vec(out_mat, out);

		// cout << "\t -- DONE --" << endl; // TEST
	};

	cout << "compute_GS.2" << endl; // TEST

	/* Allocate GS */
	gsl_matrix_complex_free(GS);
	GS = gsl_matrix_complex_calloc(chi_l*dim, dim*chi_r);

	cout << "\tGS_size = " << chi_l*dim << "x" << dim*chi_r << endl; // TEST
	cout << "compute_GS.3" << endl; // TEST

	/* Lanczos */
	lambda_lanczos::LambdaLanczos<double> engine(mv_mul, n, false, 1); 

	/* Guess for Lanczos */
	std::function<void(std::vector<double>& vec)> fn = [this](std::vector<double>& vec) { res_mat_vec(this->GUESS, vec); };
	if(guess) engine.init_vector = fn;
	cout << "\t" << guess << endl; // TEST

	/* Run Lanczos */
	std::vector<double> eigenvalues;
	std::vector<std::vector<double>> eigenvectors;
	engine.run(eigenvalues, eigenvectors);
	cout << "\teigval: " << eigenvalues[0] << endl; // TEST

	cout << "compute_GS.4" << endl; // TEST

	res_vec_mat(eigenvectors[0], GS); // NB: GS is real

	cout << "compute_GS.5" << endl; // TEST

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
		break;
		case 'r':
			return HR.back();
		break;
		default:
			error_message("Block position not allowed in loadH");
			exit(1);
	}
}

gsl_matrix_complex* sys::loadR(char pos, bool del)
{
	gsl_matrix_complex* res;
	switch(pos) 
	{
		case 'l':
			res = RL.back();
			if(del) RL.pop_back();
		break;
		case 'r':
			res = RR.back();
			if(del) RR.pop_back();
		break;
		default:
			error_message("Block position not allowed in loadR");
			exit(1);
	}
	return res;
}

/***** Delete memory inside intermediate step *****/

void sys::DeleteMem(char dir)
{
	switch(dir) 
    {
        case 'r':
            HR.pop_back();

			gsl_matrix_complex_free(RR.back());
            RR.pop_back();

            // R->reducel();
        break;
        case 'l':	
            HL.pop_back();

			gsl_matrix_complex_free(RL.back());
            RL.pop_back();
            
			// L->reducel();
        break;
        default:
            error_message("Direction not allowed in DeleteMem");
    }

}














// double sys::compute_GS()
// {
// 	/* Set dimensions */
// 	int chi_l = L->getChi();
// 	int chi_r = R->getChi();
// 	int dim = M->getDim();
// 	int n = dim*dim*chi_l*chi_r;

// 	/* Define mv_mul */
// 	auto mv_mul = [&](const std::vector<double>& in, std::vector<double>& out) 
// 	{
// 		gsl_matrix_complex * in_mat = gsl_matrix_complex_calloc(chi_l*dim, dim*chi_r);
// 		gsl_matrix_complex * out_mat = gsl_matrix_complex_calloc(chi_l*dim, dim*chi_r);
		
// 		/* Reshape in: vector -> gsl_matrix */
// 		res_vec_mat(in, in_mat);
		
// 		/* Sl phi Sr */
// 		gsl_matrix_complex * S_l = gsl_matrix_complex_calloc(chi_l*dim, dim*chi_l);
// 		gsl_matrix_complex * S_r = gsl_matrix_complex_calloc(chi_r*dim, dim*chi_r);

// 		gsl_matrix_complex * Id_l = gsl_matrix_complex_calloc(chi_l, chi_l);
// 		gsl_matrix_complex_set_identity(Id_l);
// 		gsl_matrix_complex * Id_r = gsl_matrix_complex_calloc(chi_r, chi_r);
// 		gsl_matrix_complex_set_identity(Id_r);

// 		gsl_matrix_complex * temp = gsl_matrix_complex_calloc(chi_l*dim, dim*chi_r);
// 		for(size_t i=0; i<3; i++){ 
// 			gsl_blas_zgetp(gsl_complex_rect(1,0), Id_l, M->getO(i), S_l);
// 			gsl_blas_zgetp(gsl_complex_rect(1,0), M->getO(i), Id_r, S_r);
// 			gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, M->getJ(i), S_l, in_mat, gsl_complex_rect(0,0), temp);
// 			gsl_blas_zgemm(CblasNoTrans, CblasTrans, gsl_complex_rect(1,0), temp, S_r, gsl_complex_rect(1,0), out_mat);
// 			gsl_matrix_complex_set_zero(temp);
// 			gsl_matrix_complex_set_zero(S_l);
// 			gsl_matrix_complex_set_zero(S_r);
// 		}

// 		gsl_matrix_complex_free(temp);
// 		gsl_matrix_complex_free(Id_l);
// 		gsl_matrix_complex_free(Id_r);
// 		gsl_matrix_complex_free(S_l);
// 		gsl_matrix_complex_free(S_r);

// 		/* Hl phi 1r */
// 		gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), L->getH(), in_mat, gsl_complex_rect(1,0), out_mat);
// 		/* 1l phi Hr */
// 		gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), in_mat, R->getH(), gsl_complex_rect(1,0), out_mat);  
		
// 		/* Reshape out: gsl_matrix -> vector */
// 		res_mat_vec(out_mat, out);

// 	};

// 	// std::vector<double> test_in{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};  // TEST
// 	// std::vector<double> test_in{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};            // TEST
// 	// std::vector<double> test_in{0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};            // TEST
// 	// std::vector<double> test_out{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};           // TEST
// 	// cout << "BEFORE:" << endl;                                                              // TEST
// 	// for(size_t i=0; i<size(test_out); i++) {cout << test_out[i] << endl;}                   // TEST
// 	// cout << endl;                                                                           // TEST
// 	// mv_mul(test_in, test_out);                                                              // TEST
// 	// cout << "AFTER: " << endl;                                                              // TEST
// 	// for(size_t i=0; i<size(test_out); i++) {cout << test_out[i] << endl;}                   // TEST
// 	// cout << endl;                                                                           // TEST

// 	/* Allocate GS */
// 	gsl_matrix_complex_free(GS);
// 	GS = gsl_matrix_complex_calloc(chi_l*dim, dim*chi_r);

// 	/* Lanczos */
// 	lambda_lanczos::LambdaLanczos<double> engine(mv_mul, n, false, 1); 

// 	std::vector<double> eigenvalues;
// 	std::vector<std::vector<double>> eigenvectors;

// 	// cout << "Run Lanczos" << endl;  // TEST
// 	engine.run(eigenvalues, eigenvectors);
// 	// cout << "End Lanczos" << endl;  // TEST
	
// 	res_vec_mat(eigenvectors[0], GS); // NB: GS is real

// 	// cout << endl;                                                                           // TEST
// 	// for(size_t i=0; i<size(eigenvectors[0]); i++) {cout << eigenvectors[0][i] << endl;}     // TEST
// 	// cout << "E0 = " << eigenvalues[0] << endl;                                              // TEST

// 	return eigenvalues[0];
// }
