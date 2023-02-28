#include "DMRG.h"

block::block(){
    l = 1; 

    // Set Operators    
    O = new complex<double>* [4];
    O[0] = new complex<double>[dim*dim] {1., 0., 0., 0., 1., 0., 0., 0., 1.};
    O[1] = new complex<double>[dim*dim] {0, 1/sqrt(2), 0, 1/sqrt(2), 0, 1/sqrt(2), 0, 1/sqrt(2), 0};
    O[2] = new complex<double>[dim*dim] {0, -1i/sqrt(2), 0, 1i/sqrt(2), 0, -1i/sqrt(2), 0, 1i/sqrt(2), 0};
    O[3] = new complex<double>[dim*dim] {1, 0, 0, 0, 0, 0, 0, 0, -1};

    // Initialize single site Hamiltonian
    InitHamiltonian();

}

void block::InitHamiltonian(){
    H = new complex<double>[dim*dim];

    for (int j=0; j<dim*dim; j++)
    {
        H[j] = h * O[3][j];
    }
}

