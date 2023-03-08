#define MODEL_CPP

#include "DMRG.hpp"
#include "service.hpp"

model::model(double Jx_, double Jy_, double Jz_, double h_, int dim_)
{
	// Set Hamiltonian parameters 
	Jx = gsl_complex_rect(Jx_, 0);	//read(Jx, parfile)
	Jy = gsl_complex_rect(Jy_, 0);	//read(Jy, parfile)
	Jz = gsl_complex_rect(Jz_, 0);	//read(Jz, parfile)
	h  = gsl_complex_rect(h_, 0);	//read(h, parfile)

	// Set local dimension
	dim = dim_;	//read(dim, parfile)

	// Define single site operators
	InitOps();
	// TEST: print ss ops
		// cout << "Id = " << endl;
		// gsl_matrix_complex_print(O[0]);
		// cout << "Sx = " << endl;
		// gsl_matrix_complex_print(O[1]);
		// cout << "Sy = " << endl;
		// gsl_matrix_complex_print(O[2]);
		// cout << "Sz = " << endl;
		// gsl_matrix_complex_print(O[3]);
	//

	// Initialize single site Hamiltonian
	InitH();
}

void model::InitOps()
{
	O = new gsl_matrix_complex* [4]; // O = {1, Sx, Sy, Sz}
	for(int a=0; a<4; a++){
		O[a] = gsl_matrix_complex_calloc(dim,dim);
	}

	// Identity
	gsl_matrix_complex_set_identity(O[0]);

	// Pauli
	switch(dim) 
	{
		case 3:
			// Sx
			gsl_matrix_complex_set(O[1], 0, 1, gsl_complex_rect(1/sqrt(2), 0));
			gsl_matrix_complex_set(O[1], 1, 0, gsl_complex_rect(1/sqrt(2), 0));
			gsl_matrix_complex_set(O[1], 1, 2, gsl_complex_rect(1/sqrt(2), 0));
			gsl_matrix_complex_set(O[1], 2, 1, gsl_complex_rect(1/sqrt(2), 0));
			// Sy
			gsl_matrix_complex_set(O[2], 0, 1, gsl_complex_rect(0, -1/sqrt(2)));
			gsl_matrix_complex_set(O[2], 1, 0, gsl_complex_rect(0, 1/sqrt(2)));
			gsl_matrix_complex_set(O[2], 1, 2, gsl_complex_rect(0, -1/sqrt(2)));
			gsl_matrix_complex_set(O[2], 2, 1, gsl_complex_rect(0, 1/sqrt(2)));
			// Sz
			gsl_matrix_complex_set(O[3], 0, 0, gsl_complex_rect(1, 0));
			gsl_matrix_complex_set(O[3], 2, 2, gsl_complex_rect(-1, 0));
		break;
		
		case 2:
			// Sx
			gsl_matrix_complex_set(O[1], 0, 1, gsl_complex_rect(1, 0));
			gsl_matrix_complex_set(O[1], 1, 0, gsl_complex_rect(1, 0));
			// Sy
			gsl_matrix_complex_set(O[2], 0, 1, gsl_complex_rect(0, -1));
			gsl_matrix_complex_set(O[2], 1, 0, gsl_complex_rect(0, 1));
			// Sz
			gsl_matrix_complex_set(O[3], 0, 0, gsl_complex_rect(1, 0));
			gsl_matrix_complex_set(O[3], 1, 1, gsl_complex_rect(-1, 0));
		break;
		
		default:
			error_message("spin value not yet implemented (only 1 and 1/2 are accepted)");
	}
}

void model::InitH(){
	// Hs = h*Sx
	Hs = gsl_matrix_complex_calloc(dim,dim);
	gsl_matrix_complex_memcpy(Hs, O[3]);
	gsl_matrix_complex_scale(Hs, h);

	// Hint = Jx*SxSx + Jy*SySy + Jz*SzSz
	Hint = gsl_matrix_complex_calloc(dim*dim,dim*dim);
	gsl_blas_zgetp(Jx, O[1], O[1], Hint);
	gsl_blas_zgetp(Jy, O[2], O[2], Hint);
	gsl_blas_zgetp(Jz, O[3], O[3], Hint);

	// TEST: print Hs and Hint
		cout << "Hs = " << endl;
		gsl_matrix_complex_print(Hs);
		cout << "Hint = " << endl;
		gsl_matrix_complex_print(Hint);
	//
}

// Get parameters

int model::D(){ return dim;}
gsl_complex model::JX(){ return Jx;}
gsl_complex model::JY(){ return Jy;}
gsl_complex model::JZ(){ return Jz;}