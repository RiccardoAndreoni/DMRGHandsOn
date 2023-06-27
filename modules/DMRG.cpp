#define DMRG_CPP

#include "DMRG.hpp"
#include "service.hpp"

DMRG::DMRG(double Jx_, double Jy_, double Jz_, double h_, int dim_)
{
    // construct sys
    S = new sys(Jx_, Jy_, Jz_, h_, dim_);
}

// void DMRG::run_DMRG(){
        // Infinite() until L
        
        // Remove last HR from memory
        // HR.popback();

        // Remove last RR from memory
        // RR.popback();

        // Finite until convergence (dE<eps)
        // Measure
// }

/****** Perform one step of the infinite algorithm ******/
void DMRG::Infinite()
{
    /* Renormalize */ // If nothing in memory == first step => no renorm needed
    if(!S->getRL().empty()) S->getL()->Renormalize(S->getRL().back());
    if(!S->getRR().empty()) S->getR()->Renormalize(S->getRR().back());

    /* Add sites */
    gsl_matrix_complex * temp_HLo = S->getL()->AddSite();
    gsl_matrix_complex * temp_HoR = S->getR()->AddSite();
    
    /* Compute GS */
    Egs = S->compute_GS(false); 

    /* Compute renormalization matrices */
    std::pair<gsl_matrix_complex*, gsl_matrix_complex*> R;
    R = S->compute_Rmat();

    /* Store R mats in memory vectors */
    S->saveR('l', R.first);
    S->saveR('r', R.second);

    /* Store H mats in memory vectors */
    S->saveH('l', temp_HLo);
    S->saveH('r', temp_HoR);


    // cout << "POINTERS " << S->getL()->getH() << " " << S->getR()->getH() << endl;
    // cout << "POINTERS " << S->getHL().back() << " " << S->getHR().back() << endl;


}


void DMRG::Finite(char dir)
{
    cout << "Finite.1" << endl; // TEST
    
    /* Delete memory in direction of the sweep*/
    getSYS()->DeleteMem(dir);

    // cout << "POINTERS " << S->getL()->getH() << " " << S->getR()->getH() << endl;
    // cout << "POINTERS " << S->getHL().back() << " " << S->getHR().back() << endl;

    /* Enlarge and break blocks*/
    gsl_matrix_complex * temp_Henl;
    gsl_matrix_complex * tempH;
    switch(dir) 
	{
		case 'r':	
			// Enlarge on left
                // cout << "0. " << S->getL()->getH() << " " << S->getR()->getH() << endl;

            S->getL()->Renormalize(S->getRL().back());  //Renormalize 
                // cout << "01 " << S->getL()->getH() << " " << S->getR()->getH() << endl;

            temp_Henl = S->getL()->AddSite();           // block -> block+site
                cout << "\tHL(" << S->getL()->getH()->size1 << "x" << S->getL()->getH()->size2 << ")" << endl; // TEST

            S->saveH('l', temp_Henl);            //Store left hamiltonians in memory vectors

                // cout << "1. " << S->getL()->getH() << " " << S->getR()->getH() << endl;

            // Shrink and break on right
            tempH = S->loadH('r');
                // cout << "2. " << temp_Henl << " " << tempH << endl;
                // cout << "   " << S->getL()->getH() << " " << S->getR()->getH() << endl;

                cout << "\t\tHl(" << S->getL()->getH()->size1 << "x" << S->getL()->getH()->size2 << ")" << endl;
            S->getR()->BreakBlock(tempH);
                // cout << "3. " << temp_Henl << " " << tempH << endl;
                // cout << "   " << S->getL()->getH() << " " << S->getR()->getH() << endl;
                cout << "\t\tHl(" << S->getL()->getH()->size1 << "x" << S->getL()->getH()->size2 << ")" << endl;
		break;
		case 'l':
			// Enlarge on right
            S->getR()->Renormalize(S->getRR().back());
            temp_Henl = S->getR()->AddSite();
            S->saveH('r', temp_Henl);

            // Shrink and break on left
            // getSYS()->DeleteMem(dir);

            tempH = S->loadH('l');
            S->getL()->BreakBlock(tempH); 
		break;
		default:
			error_message("Direction not allowed in Finite");
	}
    
    cout << "Finite.2" << endl; // TEST
	// cout << "\t\tHl(" << S->getL()->getH()->size1 << "x" << S->getL()->getH()->size2 << ")" << endl; // TEST

	cout << "\t\tRL(" << S->getRL().back()->size1 << "x" << S->getRL().back()->size2 << ")" << endl;    // TEST
	cout << "\t\tRR(" << S->getRR().back()->size1 << "x" << S->getRR().back()->size2 << ")" << endl;    // TEST

    /* Guess for the GS */
    S->GStoGuess(dir);
    // gsl_matrix_complex_print(getSYS()->getGS());                            // TEST
    // cout << "\tNORM = " << gsl_matrix_complex_norm(getSYS()->getGS());      // TEST
    cout << endl;
    // gsl_matrix_complex_print(getSYS()->getGUESS());                         // TEST
    // cout << "\tNORM = " << gsl_matrix_complex_norm(getSYS()->getGUESS());   // TEST


	cout << "\t\tH_l(" << S->getL()->getH()->size1 << "x" << S->getL()->getH()->size2 << ")" << endl;
	cout << "\t\tH_r(" << S->getR()->getH()->size1 << "x" << S->getR()->getH()->size2 << ")" << endl;

    cout << "Finite.3" << endl; // TEST
    
    /* Compute GS */
    cout << "size HL = " << S->getHL().size() << ", size HR = " << S->getHR().size() << endl;
    Egs = S->compute_GS(true);

    cout << "Finite.4" << endl; // TEST
    
    /* Compute renormalization matrices */
    std::pair<gsl_matrix_complex*, gsl_matrix_complex*> R;
    R = S->compute_Rmat();

    cout << "Finite.5" << endl; // TEST
    /* Store R in memory */
    switch(dir)
    {
        case 'r':
            S->saveR('l', R.first);              //Storing the left one, right one already freed in GStoGuess
        break;
        case 'l':
            S->saveR('r', R.second);
        break;
        default:
            error_message("Direction not allowed in Finite");
    }
}

void DMRG::Sweeps()
{
    while(S->getR()->getl()>1) Finite('r');
    while(S->getL()->getl()>1) Finite('l');
    while(S->getR()->getl()!=S->getL()->getl()) Finite('r');
}