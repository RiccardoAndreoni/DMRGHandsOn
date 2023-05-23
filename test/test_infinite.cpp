#define MAIN

#include "DMRG.hpp"
#include "service.hpp"

int main(){



    DMRG * dmrg = new DMRG(1,2,3,4,2);
    
    double En;

    std::vector<gsl_matrix_complex*> HmatsR;
    std::vector<gsl_matrix_complex*> HmatsL;
    std::vector<gsl_matrix_complex*> RmatsR;
    std::vector<gsl_matrix_complex*> RmatsL;



int it=2;
for(int i=0; i<it; i++){

    dmrg -> Infinite();
    En = dmrg -> getEgs();

    cout << "Groundstate energy "<< 4+2*i<<" sites: " << En << endl;
    
    gsl_matrix_complex_print(dmrg->getSYS()->getL()->getH());
    cout<<endl;
    gsl_matrix_complex_print(dmrg->getSYS()->getL()->getS(0,1));
    cout<<endl;
    cout << dmrg->getSYS()->getL()->getChi() << endl;
    cout<<endl;

    HmatsL= dmrg -> getHL();
    HmatsR= dmrg -> getHR();
    RmatsL= dmrg -> getRL();
    RmatsR= dmrg -> getRR(); 

    cout << "================ Hamiltonian Memory test ===================" << endl;
    gsl_matrix_complex_print(HmatsL[i]);
    cout << endl;
    gsl_matrix_complex_print(HmatsR[i]);
    cout << endl;
    cout << "================ Renormalization Memory test ===================" << endl;
    cout << endl;
    gsl_matrix_complex_print(RmatsL[i]);
    cout<< endl;
    gsl_matrix_complex_print(RmatsR[i]);

}

    // cout << "len = " << HmatsL.size() << endl;
    // cout << endl;

    // cout << "HmatsL[0]" << endl;
    // gsl_matrix_complex_print(HmatsL[0]);
    // cout << "HmatsR[0]" << endl;
    // gsl_matrix_complex_print(HmatsR[0]);

    // for (int i=0; i<it;i++){

        // cout << "================ Hamiltonian Memory test ===================" << endl;
        // gsl_matrix_complex_print(HmatsL[i]);
        // cout << endl;
        // gsl_matrix_complex_print(HmatsR[i]);
        // cout << endl;
        // cout << "================ Renormalization Memory test ===================" << endl;
        // cout << endl;
        // gsl_matrix_complex_print(RmatsL[i]);
        // cout<< endl;
        // gsl_matrix_complex_print(RmatsR[i]);

    // }

    return 0;
}