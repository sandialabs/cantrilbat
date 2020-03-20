/*
 *  $Author: hkmoffa $
 *  $Date: 2013-01-07 15:54:04 -0700 (Mon, 07 Jan 2013) $
 *  $Revision: 508 $
 *
 */

#include "zuzax/IdealGasMix.h"
#include "zuzax/equilibrium.h"

#include <cstdio>

using namespace std;
#ifdef useZuzaxNamespace
using namespace Zuzax;
#define Zuzax Zuzax
#else
using namespace Cantera;
#define Zuzax Cantera
#endif

int main(int argc, char **argv) {
  try {
#ifdef DEBUG_BASISOPTIMIZE
    Zuzax::BasisOptimize_print_lvl = 1;
#endif
    IdealGasMix g("gri30.xml", "gri30_mix");

    double pres = 1.0E5;
    g.setState_TPX(1500.0, pres, "CH4:0.3, O2:0.3, N2:0.4");

    MultiPhase mphase;
    mphase.addPhaseWithMoles(&g, 10.0);
    mphase.init();
    int usedZeroedSpecies = 0;
    vector<size_t> orderVectorSpecies;
    vector<size_t> orderVectorElements;

    bool doFormMatrix = false;
    vector_fp formRxnMatrix;

    int nc = Zuzax::BasisOptimize(&usedZeroedSpecies, doFormMatrix,
			   &mphase, orderVectorSpecies, orderVectorElements,
			   formRxnMatrix);
    cout << "number of components = " << nc << endl;

    vector_fp elementAbundances;
    int nct = Zuzax::ElemRearrange(nc, elementAbundances, &mphase, 
				     orderVectorSpecies, orderVectorElements);
    if (nc != nct) {
      printf("ERROR\n");
      exit(-1);
    }

    doFormMatrix = true;
    nc = BasisOptimize(&usedZeroedSpecies, doFormMatrix,
			   &mphase, orderVectorSpecies, orderVectorElements,
			   formRxnMatrix);

    cout << "number of components = " << nc << endl;

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
