#define MAIN

#include "DMRG.hpp"
#include "service.hpp"

int main(){

    DMRG * dmrg = new DMRG(1,2,3,4,2);
        
    double En;
    char dir = 'r';


    int it=30;
    for(int i=0; i<it; i++){
        dmrg -> Infinite();
        En = dmrg -> getEgs();

        // cout << "======================================" << endl;
        // cout << "======================================" << endl;

        // cout << "Groundstate energy "<< 4+2*i<<" sites: " << En << endl;

        // cout << "================ L ===================" << endl;
        // gsl_matrix_complex_print(dmrg->getSYS()->getL()->getH());
        // cout<<endl;
        // cout << "ChiL = " << dmrg->getSYS()->getL()->getChi() << endl;
        // cout << "lL = "   << dmrg->getSYS()->getL()->getl() << endl;
        // cout<<endl;
        // cout << "================ R ===================" << endl;
        // gsl_matrix_complex_print(dmrg->getSYS()->getR()->getH());
        // cout<<endl;
        // cout << "ChiR = " << dmrg->getSYS()->getR()->getChi() << endl;
        // cout << "lR = "   << dmrg->getSYS()->getR()->getl() << endl;
        // cout<<endl;
    }

    cout << "==============BEFORE==============" << endl;
    cout << "memRL.size = " << dmrg->getSYS()->getRL().size() << endl;
    cout << "memRR.size = " << dmrg->getSYS()->getRR().size() << endl;
    cout << "memHL.size = " << dmrg->getSYS()->getHL().size() << endl;
    cout << "memHR.size = " << dmrg->getSYS()->getHR().size() << endl;

    cout << "########### Infinite done ###########" << endl;

    cout << "POINTERS " << dmrg->getSYS()->getL()->getH() << " " << dmrg->getSYS()->getR()->getH() << endl;
    cout << "POINTERS " << dmrg->getSYS()->getHL().back() << " " << dmrg->getSYS()->getHR().back() << endl;
    
    cout << "==============AFTER===============" << endl;
    cout << "memRL.size = " << dmrg->getSYS()->getRL().size() << endl;
    cout << "memRR.size = " << dmrg->getSYS()->getRR().size() << endl;
    cout << "memHL.size = " << dmrg->getSYS()->getHL().size() << endl;
    cout << "memHR.size = " << dmrg->getSYS()->getHR().size() << endl;

    ///////////////////////////////////////////////////////

    cout << "RL sizes: ";
    for(size_t i=0; i<dmrg->getSYS()->getRL().size(); i++)
    {
        cout << "(" << dmrg->getSYS()->getRL()[i]->size1 << "x" << dmrg->getSYS()->getRL()[i]->size2 << "), ";
    }
    cout << endl;

    cout << "RR sizes: ";
    for(size_t i=0; i<dmrg->getSYS()->getRR().size(); i++)
    {
        cout << "(" << dmrg->getSYS()->getRR()[i]->size1 << "x" << dmrg->getSYS()->getRR()[i]->size2 << "), ";
    }
    cout << endl;

    ///////////////////////////////////////////////////////

    cout << "HL sizes: ";
    for(size_t i=0; i<dmrg->getSYS()->getHL().size(); i++)
    {
        cout << "(" << dmrg->getSYS()->getHL()[i]->size1 << "x" << dmrg->getSYS()->getHL()[i]->size2 << "), ";
    }
    cout << endl;

    cout << "HR sizes: ";
    for(size_t i=0; i<dmrg->getSYS()->getHR().size(); i++)
    {
        cout << "(" << dmrg->getSYS()->getHR()[i]->size1 << "x" << dmrg->getSYS()->getHR()[i]->size2 << "), ";
    }
    cout << endl;

    ///////////////////////////////////////////////////////

    // cout << "POINTERS " << dmrg->getSYS()->getL()->getH() << " " << dmrg->getSYS()->getR()->getH() << endl;
    // cout << "POINTERS " << dmrg->getSYS()->getHL().back() << " " << dmrg->getSYS()->getHR().back() << endl;


    for(int i=0; i<it; i++){

        cout << "======================================" << endl;
        cout << "======================================" << endl;
        cout << "===== FINITE STEP "<< i+1 << " ===== " << endl;
        cout << "======================================" << endl;
        cout << "======================================" << endl;

        dmrg -> Finite(dir);
        En = dmrg -> getEgs();

        cout << "======================================" << endl;
        cout << "======================================" << endl;

        cout << "Groundstate energy "<< 4+2*i<<" sites: " << En << endl;

        cout << "================ L ===================" << endl;
        gsl_matrix_complex_print(dmrg->getSYS()->getL()->getH());
        cout<<endl;
        cout << "ChiL = " << dmrg->getSYS()->getL()->getChi() << endl;
        cout << "lL = "   << dmrg->getSYS()->getL()->getl() << endl;
        cout<<endl;
        cout << "================ R ===================" << endl;
        gsl_matrix_complex_print(dmrg->getSYS()->getR()->getH());
        cout<<endl;
        cout << "ChiR = " << dmrg->getSYS()->getR()->getChi() << endl;
        cout << "lR = "   << dmrg->getSYS()->getR()->getl() << endl;
        cout<<endl;

        cout << "==============MEMORY===============" << endl;
        cout << "memRL.size = " << dmrg->getSYS()->getRL().size() << endl;
        cout << "memRR.size = " << dmrg->getSYS()->getRR().size() << endl;
        cout << "memHL.size = " << dmrg->getSYS()->getHL().size() << endl;
        cout << "memHR.size = " << dmrg->getSYS()->getHR().size() << endl;

        cout << "RL sizes: ";
        for(size_t i=0; i<dmrg->getSYS()->getRL().size(); i++)
        {
            cout << "(" << dmrg->getSYS()->getRL()[i]->size1 << "x" << dmrg->getSYS()->getRL()[i]->size2 << "), ";
        }
        cout << endl;

        cout << "RR sizes: ";
        for(size_t i=0; i<dmrg->getSYS()->getRR().size(); i++)
        {
            cout << "(" << dmrg->getSYS()->getRR()[i]->size1 << "x" << dmrg->getSYS()->getRR()[i]->size2 << "), ";
        }
        cout << endl;

        cout << "HL sizes: ";
        for(size_t i=0; i<dmrg->getSYS()->getHL().size(); i++)
        {
            cout << "(" << dmrg->getSYS()->getHL()[i]->size1 << "x" << dmrg->getSYS()->getHL()[i]->size2 << "), ";
        }
        cout << endl;

        cout << "HR sizes: ";
        for(size_t i=0; i<dmrg->getSYS()->getHR().size(); i++)
        {
            cout << "(" << dmrg->getSYS()->getHR()[i]->size1 << "x" << dmrg->getSYS()->getHR()[i]->size2 << "), ";
        }
        cout << endl;


    }

        return 0;
}