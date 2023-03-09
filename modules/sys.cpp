#define SYS_CPP

#include "DMRG.hpp"
#include "service.hpp"

sys::sys(double Jx_, double Jy_, double Jz_, double h_, int dim_)
{
    M = new model(Jx_, Jy_, Jz_, h_, dim_);
    L = new block(*M);
    R = new block(*M);
}



/* Add site to the left block to build HLo from HL */

void sys::add_site_L(gsl_matrix_complex * H)
{
    int chi_l = L->getChi();
    gsl_matrix_complex * Id = gsl_matrix_complex_calloc(chi_l, chi_l)

    gsl_blas_zgetp(gsl_complex_rect(1, 0), L->getH(), M->getO(0), H);
    gsl_blas_zgetp(gsl_complex_rect(1, 0), Id , M->getHs(0), H);
    
    for(size_t i=1; i<4; i++)
    { gsl_blas_zgetp(gsl_complex_rect(1, 0), L->getS(L->getl(), i-1), M->getO(i), H); }
}

/* Add site to the left block to build HoR from HR */

void sys::add_site_R(gsl_matrix_complex * H)
{
}

/* GS computation */

void sys::compute_GS()
{
    // Set dimensions
    int chi_l = L->getChi();
    int chi_r = R->getChi();
    int dim = M->getDim();
    int n = dim*dim*chi_l*chi_r;

    // Define Hlo and Hor starting from HL and HR
    gsl_matrix_complex * HLo = gsl_matrix_complex_calloc(chi_l*dim, chi_l*dim);
    add_site_L(HLo);
    gsl_matrix_complex * HoR = gsl_matrix_complex_calloc(dim*chi_r, dim*chi_r); 
    add_site_R(HoR);

    // Define mv_mul
    auto mv_mul = [&](const std::vector<double>& in, std::vector<double>& out) 
    {
        gsl_matrix_complex * in_mat = gsl_matrix_complex_calloc(chi_l*dim, dim*chi_r);
        gsl_matrix_complex * out_mat = gsl_matrix_complex_calloc(chi_l*dim, dim*chi_r);
        
        // Reshape in: vector -> gsl_matrix 
        res_vec_mat(in, in_mat);
        
        // Sl phi Sr
        gsl_matrix_complex * temp = gsl_matrix_complex_calloc(chi_l*dim, dim*chi_r);
        gsl_matrix_complex * S_l = gsl_matrix_complex_calloc(chi_l*dim, dim*chi_l);
        gsl_matrix_complex * S_r = gsl_matrix_complex_calloc(chi_r*dim, dim*chi_r);

        gsl_matrix_complex * Id_l = gsl_matrix_complex_calloc(chi_l, chi_l);
        gsl_matrix_complex_set_identity(Id_l);
        gsl_matrix_complex * Id_r = gsl_matrix_complex_calloc(chi_r, chi_r);
        gsl_matrix_complex_set_identity(Id_r);

        for(size_t i=1; i<4; i++){ // Start from 1 bc we need only Sx, Sy, Sz and not Id in O[]
            gsl_blas_zgetp(gsl_complex_rect(1,0), Id_l, M->getO(i), S_l);
            gsl_blas_zgetp(gsl_complex_rect(1,0), M->getO(i), Id_r, S_r);
            gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, M->getJ(i-1), S_l, in_mat, gsl_complex_rect(0,0), temp);
            gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), temp, S_r, gsl_complex_rect(0,0), out_mat);
            gsl_matrix_complex_set_zero(temp);
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

    // Lanczos
    lambda_lanczos::LambdaLanczos<double> engine(mv_mul, n, false, 1); 
    std::vector<double> eigenvalues;
    std::vector<std::vector<double>> eigenvectors;
    engine.run(eigenvalues, eigenvectors);
    res_vec_mat(eigenvectors[0], GS);
}



// void sys::Ham_mul(gsl_matrix_complex * in_mat, gsl_matrix_complex * out_mat){
//     gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1, 0), )
// }
