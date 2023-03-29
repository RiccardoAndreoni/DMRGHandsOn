#define MAIN

#include "DMRG.hpp"
#include "service.hpp"

int main(){



    DMRG * dmrg = new DMRG(1,2,3,4,2);
    
    dmrg -> Infinite();
    double En;
    En = dmrg -> getEgs();

    cout << "Groundstate energy: " << En << endl;
    
    gsl_matrix_complex_print(dmrg->getSYS()->getL()->getH());
    cout<<endl;

    gsl_matrix_complex_print(dmrg->getSYS()->getL()->getS(0,1));
    cout<<endl;

    cout << dmrg->getSYS()->getL()->getChi() << endl;
    cout<<endl;


// ------------------------------End of first iteration--------------

    dmrg -> Infinite();
    cout << "here 1" << endl;

    En = dmrg -> getEgs();
    cout << "here 2" << endl;

    cout << "Groundstate energy: " << En << endl;

    return 0;
}