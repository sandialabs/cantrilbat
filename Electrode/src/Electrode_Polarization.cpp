
#include "Electrode_Polarization.h"
#include "cantera/base/ctexceptions.h"

#include <cstdio>

//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
//==================================================================================================================================
PolarizationSurfRxnResults::PolarizationSurfRxnResults(int electrodeDomainNumber, int electrodeCellNumber, 
                                                       size_t surfIndex, size_t rxnIndex) :
    isurf(surfIndex),
    iRxnIndex(rxnIndex),
    electrodeDomainNumber_(electrodeDomainNumber),
    electrodeCellNumber_(electrodeCellNumber)
{

}
//==================================================================================================================================
void PolarizationSurfRxnResults::addSubStep(struct PolarizationSurfRxnResults& sub)
{


}
//==================================================================================================================================

} 
//----------------------------------------------------------------------------------------------------------------------------------

