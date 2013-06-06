/*
 * $Id: Electrode_Jacobian.cpp 496 2013-01-07 21:15:37Z hkmoffa $
 */
#include "cantera/base/mdp_allo.h"


#include "cantera/thermo/FixedChemPotSSTP.h"

#include "Electrode_Jacobian.h"

using namespace Cantera;
using namespace std;

#ifndef SAFE_DELETE
#define SAFE_DELETE(x)  if (x) { delete x;  x = 0;}
#endif
#ifndef MAX
#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#endif

namespace Cantera {

//====================================================================================================================
Electrode_Jacobian::Electrode_Jacobian(Electrode* elect) :
                ee_(elect),
                printLvl_(0)
{
}
//====================================================================================================================
Electrode_Jacobian::Electrode_Jacobian(const Electrode_Jacobian& right) :
                ee_(right.ee_),
                printLvl_(0)
{
    /*
     * Call the assignment operator.
     */
    operator=(right);
}
//====================================================================================================================
// Destructor
Electrode_Jacobian::~Electrode_Jacobian()
{
}
//======================================================================================================================
// Assignment operator
/*
 *  @param right object to be copied
 */
Electrode_Jacobian& Electrode_Jacobian::operator=(const Electrode_Jacobian& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }
    /*
     *  Do a shallow copy of the Electrode pointer. This is all that is necessary
     */
    ee_ = right.ee_;


    printLvl_ = right.printLvl_;

    /*
     * Return the reference to the current object
     */
    return *this;
}

//====================================================================================================================
}// End of namespace Cantera
//======================================================================================================================

