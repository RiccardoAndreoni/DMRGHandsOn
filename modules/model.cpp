#define MODEL_CPP

#include "DMRG.hpp"
#include "service.hpp"

model::model(double Jx_, double Jy_, double Jz_, double h_, int dim_)
{
	// Set Hamiltonian parameters
	J = new gsl_complex[3];
	J[0] = gsl_complex_rect(Jx_, 0);	//read(Jx, parfile)
	J[1] = gsl_complex_rect(Jy_, 0);	//read(Jy, parfile)
	J[2] = gsl_complex_rect(Jz_, 0);	//read(Jz, parfile)
	h  = gsl_complex_rect(h_, 0);	//read(h, parfile)

	// Set local dimension
	dim = dim_;	//read(dim, parfile)

	// Define single site operators
	InitOps();

	// Initialize single site Hamiltonian
	InitH();
}

void model::InitOps()
{
	O = new gsl_matrix_complex* [3]; // O = {Sx, Sy, Sz}
	for(int a=0; a<3; a++){
		O[a] = gsl_matrix_complex_calloc(dim,dim);
	}

	Id = gsl_matrix_complex_alloc(dim,dim);
	// Identity
	gsl_matrix_complex_set_identity(Id);

	// Pauli
	switch(dim) 
	{
		case 3:
			// Sx
			gsl_matrix_complex_set(O[0], 0, 1, gsl_complex_rect(1/sqrt(2), 0));
			gsl_matrix_complex_set(O[0], 1, 0, gsl_complex_rect(1/sqrt(2), 0));
			gsl_matrix_complex_set(O[0], 1, 2, gsl_complex_rect(1/sqrt(2), 0));
			gsl_matrix_complex_set(O[0], 2, 1, gsl_complex_rect(1/sqrt(2), 0));
			// Sy
			gsl_matrix_complex_set(O[1], 0, 1, gsl_complex_rect(0, -1/sqrt(2)));
			gsl_matrix_complex_set(O[1], 1, 0, gsl_complex_rect(0, 1/sqrt(2)));
			gsl_matrix_complex_set(O[1], 1, 2, gsl_complex_rect(0, -1/sqrt(2)));
			gsl_matrix_complex_set(O[1], 2, 1, gsl_complex_rect(0, 1/sqrt(2)));
			// Sz
			gsl_matrix_complex_set(O[2], 0, 0, gsl_complex_rect(1, 0));
			gsl_matrix_complex_set(O[2], 2, 2, gsl_complex_rect(-1, 0));
		break;
		
		case 2:
			// Sx
			gsl_matrix_complex_set(O[0], 0, 1, gsl_complex_rect(1, 0));
			gsl_matrix_complex_set(O[0], 1, 0, gsl_complex_rect(1, 0));
			// Sy
			gsl_matrix_complex_set(O[1], 0, 1, gsl_complex_rect(0, -1));
			gsl_matrix_complex_set(O[1], 1, 0, gsl_complex_rect(0, 1));
			// Sz
			gsl_matrix_complex_set(O[2], 0, 0, gsl_complex_rect(1, 0));
			gsl_matrix_complex_set(O[2], 1, 1, gsl_complex_rect(-1, 0));
		break;
		
		default:
			error_message("spin value not yet implemented (only 1 and 1/2 are accepted)");
	}
}

void model::InitH(){
	// Hs = h*Sz
	Hs = gsl_matrix_complex_calloc(dim,dim);
	gsl_matrix_complex_memcpy(Hs, O[2]);
	gsl_matrix_complex_scale(Hs, h);

	// Hint = Jx*SxSx + Jy*SySy + Jz*SzSz 
	// Hint = gsl_matrix_complex_calloc(dim*dim,dim*dim);
	// gsl_blas_zgetp(J[0], O[0], O[0], Hint);
	// gsl_blas_zgetp(J[1], O[1], O[1], Hint);
	// gsl_blas_zgetp(J[2], O[2], O[2], Hint);
	// gsl_blas_zgetp(h, Id, O[2], Hint);
	// gsl_blas_zgetp(h, O[2], Id, Hint);
}