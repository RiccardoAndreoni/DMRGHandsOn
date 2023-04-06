#define MAIN

#include "DMRG.hpp"
#include "service.hpp"

int main(){



    DMRG * dmrg = new DMRG(1,2,3,4,2);
    
    dmrg -> Infinite();
    double En;
    En = dmrg -> getEgs();

    cout << "Groundstate energy 4 sites: " << En << endl;
    
    //gsl_matrix_complex_print(dmrg->getSYS()->getL()->getH());
    //cout<<endl;
    //gsl_matrix_complex_print(dmrg->getSYS()->getL()->getS(0,1));
    //cout<<endl;
    //cout << dmrg->getSYS()->getL()->getChi() << endl;
    //cout<<endl;


// ------------------------------End of first iteration--------------

    dmrg -> Infinite();

    En = dmrg -> getEgs();

    cout << "Groundstate energy 6 sites: " << En << endl;

    dmrg -> Infinite();

    En = dmrg -> getEgs();

    cout << "Groundstate energy 8 sites: " << En << endl;

    dmrg -> Infinite();

    En = dmrg -> getEgs();

    cout << "Groundstate energy 10 sites: " << En << endl;

    return 0;
}