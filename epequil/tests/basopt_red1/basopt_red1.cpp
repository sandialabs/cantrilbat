/*
 *  $Author: hkmoffa $
 *  $Date: 2013-01-07 15:54:04 -0700 (Mon, 07 Jan 2013) $
 *  $Revision: 508 $
 *
 */
#ifndef DEBUG_HKM
#define DEBUG_HKM
#endif

#include "cantera/IdealGasMix.h"
#include "cantera/equilibrium.h"

using namespace std;
using namespace Cantera;

int main(int argc, char **argv) {
  try {
    IdealGasMix g("red1.xml", "gri30_mix");

#ifdef DEBUG_BASISOPTIMIZE
   Cantera::BasisOptimize_print_lvl = 1;
#endif
#ifdef DEBUG_CHEMEQUIL
   Cantera::ChemEquil_print_lvl = 1;
#endif

    double pres = 1.0E5;
    g.setState_TPX(2000.0, pres, "C2H2:0.9, CH:0.1");

    MultiPhase mphase;
    mphase.addPhase(&g, 10.0);
    mphase.init();
    int usedZeroedSpecies = 0;
    vector<size_t> orderVectorSpecies;
    vector<size_t> orderVectorElements;

    bool doFormMatrix = true;
    vector_fp formRxnMatrix;

    int nc = BasisOptimize(&usedZeroedSpecies, doFormMatrix,
		           &mphase, orderVectorSpecies,
                           orderVectorElements,
                           formRxnMatrix);

    cout << "number of components = " << nc << endl;

    /*
     * The ChemEquil solver throws an error for this case.
     * The MultiPhaseEquil solver just gets the wrong result.
     */
    equilibrate(g, "TP", -1);
    cout << g;
 
    return 0;
  }
  catch (CanteraError) {
    showErrors(cerr);
    cerr << "program terminating." << endl;
    return -1;
  }
}
