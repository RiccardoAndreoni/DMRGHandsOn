#define SYS_CPP

#include "DMRG.hpp"
#include "service.hpp"

sys::sys(double Jx_, double Jy_, double Jz_, double h_, int dim_)
{
    M = new model(Jx_, Jy_, Jz_, h_, dim_);
    L = new block(M);
    R = new block(M);

    GS = gsl_matrix_complex_calloc(2,2);
}

/* Add site to the left block to build HLo from HL */

void sys::computeHLo(gsl_matrix_complex * H)
{
    int chi_l = L->getChi();
    gsl_matrix_complex * Id = gsl_matrix_complex_calloc(chi_l, chi_l);
    gsl_matrix_complex_set_identity(Id);

    gsl_blas_zgetp(gsl_complex_rect(1, 0), L->getH(), M->getId(), H);
    gsl_blas_zgetp(gsl_complex_rect(1, 0), Id , M->getHs(), H);

    for(size_t i=0; i<3; i++)
    { gsl_blas_zgetp(M->getJ(i), L->getS(L->getl()-1, i), M->getO(i), H); }
}

/* Add site to the left block to build HoR from HR */

void sys::computeHoR(gsl_matrix_complex * H)
{
    int chi_r = R->getChi();
    gsl_matrix_complex * Id = gsl_matrix_complex_calloc(chi_r, chi_r);
    gsl_matrix_complex_set_identity(Id);

    gsl_blas_zgetp(gsl_complex_rect(1, 0), M->getId(), R->getH(), H);
    gsl_blas_zgetp(gsl_complex_rect(1, 0), M->getHs(), Id , H);

    for(size_t i=0; i<3; i++)
    { gsl_blas_zgetp(M->getJ(i), M->getO(i), R->getS(0, i), H); }
}

/* GS computation */

double sys::compute_GS()
{
    // Set dimensions
    int chi_l = L->getChi();
    int chi_r = R->getChi();
    int dim = M->getDim();
    int n = dim*dim*chi_l*chi_r;

    // Allocate GS
    gsl_matrix_complex_free(GS);
    GS = gsl_matrix_complex_calloc(chi_l*dim, dim*chi_r);

    // Define Hlo and Hor starting from HL and HR
    gsl_matrix_complex * HLo = gsl_matrix_complex_calloc(chi_l*dim, chi_l*dim);
    computeHLo(HLo);
    gsl_matrix_complex * HoR = gsl_matrix_complex_calloc(dim*chi_r, dim*chi_r); 
    computeHoR(HoR);

    // cout << "ComputeGS" << endl;    // TEST
    // cout << "HLo = " << endl;       // TEST
    // gsl_matrix_complex_print(HLo);  // TEST
    // cout << "HoR = " << endl;       // TEST
    // gsl_matrix_complex_print(HoR);  // TEST

    // Define mv_mul
    auto mv_mul = [&](const std::vector<double>& in, std::vector<double>& out) 
    {
        gsl_matrix_complex * in_mat = gsl_matrix_complex_calloc(chi_l*dim, dim*chi_r);
        gsl_matrix_complex * out_mat = gsl_matrix_complex_calloc(chi_l*dim, dim*chi_r);
        
        // Reshape in: vector -> gsl_matrix 
        res_vec_mat(in, in_mat);
        
        // Sl phi Sr
        gsl_matrix_complex * S_l = gsl_matrix_complex_calloc(chi_l*dim, dim*chi_l);
        gsl_matrix_complex * S_r = gsl_matrix_complex_calloc(chi_r*dim, dim*chi_r);

        gsl_matrix_complex * Id_l = gsl_matrix_complex_calloc(chi_l, chi_l);
        gsl_matrix_complex_set_identity(Id_l);
        gsl_matrix_complex * Id_r = gsl_matrix_complex_calloc(chi_r, chi_r);
        gsl_matrix_complex_set_identity(Id_r);

        gsl_matrix_complex * temp = gsl_matrix_complex_calloc(chi_l*dim, dim*chi_r);
        for(size_t i=0; i<3; i++){ // Start from 1 bc we need only Sx, Sy, Sz and not Id in O[]
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

        // Hl phi 1r
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), HLo, in_mat, gsl_complex_rect(1,0), out_mat);
        // 1l phi Hr
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), in_mat, HoR, gsl_complex_rect(1,0), out_mat);  
        
        // Reshape out: gsl_matrix -> vector 
        res_mat_vec(out_mat, out);

    };

    // // std::vector<double> test_in{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};  // TEST
    // std::vector<double> test_in{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};            // TEST
    // // std::vector<double> test_in{0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};            // TEST
    // std::vector<double> test_out{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};           // TEST
    // cout << "BEFORE:" << endl;                                                              // TEST
    // for(size_t i=0; i<size(test_out); i++) {cout << test_out[i] << endl;}                   // TEST
    // cout << endl;                                                                           // TEST
    // mv_mul(test_in, test_out);                                                              // TEST
    // cout << "AFTER: " << endl;                                                              // TEST
    // for(size_t i=0; i<size(test_out); i++) {cout << test_out[i] << endl;}                   // TEST
    // cout << endl;                                                                           // TEST

    // Lanczos
    lambda_lanczos::LambdaLanczos<double> engine(mv_mul, n, false, 1); 
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

/* SVD and resize -> define RL and RR*/

void sys::compute_Rmat()
{
    // Set dimensions
    int chi_l = L->getChi();
    int chi_r = R->getChi();
    int dim = M->getDim();
    // int n = dim*dim*chi_l*chi_r;

    // Definition of temporary U and V, to be truncated after SVD
    gsl_matrix * tempU = gsl_matrix_alloc(chi_l*dim, dim*chi_r);
    conv_comp_real(tempU, GS);
    gsl_matrix * tempV = gsl_matrix_alloc(dim*chi_r, dim*chi_r);
    gsl_vector * tempS = gsl_vector_alloc(dim*chi_r);
    gsl_vector * tempW = gsl_vector_alloc(dim*chi_r);

    gsl_linalg_SV_decomp(tempU, tempV, tempS, tempW);

    // cout << "U = " << endl;     // TEST
    // gsl_matrix_print(tempU);    // TEST
    // cout << endl;               // TEST
    // cout << "V = " << endl;     // TEST
    // gsl_matrix_print(tempV);    // TEST
    // cout << endl;               // TEST
    // cout << "S = " << endl;     // TEST
    // gsl_vector_print(tempS);    // TEST

    // Truncation
    if(tempS->size > chimax)
    {
        //RL = Utruncated
        //RR = Vtruncated
    }
    else
    {
        //RL = Utemp
        //RR = Vtemp
    }
}