/*
 *  $Author: hkmoffa $
 *  $Date: 2013-01-07 15:54:04 -0700 (Mon, 07 Jan 2013) $
 *  $Revision: 508 $
 *
 */
#ifndef DEBUG_HKM
#define DEBUG_HKM
#endif

#include "zuzax/IdealGasMix.h"
#include "zuzax/equilibrium.h"

using namespace std;
#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif

int main(int argc, char **argv) {
  try {
    IdealGasMix g("red1.xml", "gri30_mix");

#ifdef DEBUG_BASISOPTIMIZE
   Zuzax::BasisOptimize_print_lvl = 1;
#endif
#ifdef DEBUG_CHEMEQUIL
   Zuzax::ChemEquil_print_lvl = 1;
#endif

    double pres = 1.0E5;
    g.setState_TPX(2000.0, pres, "C2H2:0.9, CH:0.1");

    MP_EquilStatic mphase;
    mphase.addPhase(&g, 10.0);
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
  catch (ZuzaxError) {
    showErrors(cerr);
    cerr << "program terminating." << endl;
    return -1;
  }
}
