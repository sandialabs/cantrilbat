/*
 * $Id: Electrode_SurfaceRegion.cpp 298 2012-08-08 20:15:48Z hkmoffa $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tok_input_util.h"



#include "cantera/base/mdp_allo.h"

#include "Electrode_SurfaceRegion.h"
#include "Electrode_RadialDiffRegions.h"
#include "cantera/integrators.h"

#include "BlockEntryGlobal.h"

using namespace Cantera;
using namespace std;
using namespace BEInput;
using namespace TKInput;
using namespace mdpUtil;


#ifndef MAX
#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x,y) (( (x) < (y) ) ? (x) : (y))
#endif

namespace Cantera
{


//======================================================================================================================
/*
 *  ELECTRODE_INPUT: constructor
 *
 *  We initialize the arrays in the structure to the appropriate sizes.
 *  And, we initialize all of the elements of the arrays to defaults.
 */
Electrode_SurfaceRegion::Electrode_SurfaceRegion() :
    Electrode_Integrator(),
    rsd_(0),
    indexSurfaceRegion_(-1),
    numSurfaceSpecies_(0),
    numSurfacePhases_(0),
    numKRSpeciesLeft_(0),
    indexRNode_(-1),
    rnodePos_final_(0.0),
    rnodePos_init_(0.0),
    rnodePos_init_init_(0.0),
    rRefPos_final_final_(0.0),
    rRefPos_final_(0.0),
    rRefPos_init_(0.0),
    rRefPos_init_init_(0.0)

{


}

//======================================================================================================================
// Copy Constructor
/*
 * @param right Object to be copied
 */
Electrode_SurfaceRegion::Electrode_SurfaceRegion(const Electrode_SurfaceRegion& right) :
    Electrode_Integrator(),
    rsd_(0),
    indexSurfaceRegion_(-1),
    numSurfaceSpecies_(0),
    numSurfacePhases_(0),
    numKRSpeciesLeft_(0),
    indexRNode_(-1),
    rnodePos_final_(0.0),
    rnodePos_init_(0.0),
    rnodePos_init_init_(0.0),
    rRefPos_final_final_(0.0),
    rRefPos_final_(0.0),
    rRefPos_init_(0.0),
    rRefPos_init_init_(0.0)
{
    /*
     * Call the assignment operator.
     */
    *this = operator=(right);
}
//======================================================================================================================
// Assignment operator
/*
 *  @param right object to be copied
 */
Electrode_SurfaceRegion&
Electrode_SurfaceRegion::operator=(const Electrode_SurfaceRegion& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    Electrode_Integrator::operator=(right);

    rsd_                      = right.rsd_;
    indexSurfaceRegion_       = right.indexSurfaceRegion_;
    numSurfaceSpecies_        = right.numSurfaceSpecies_;
    numSurfacePhases_         = right.numSurfacePhases_;
    numKRSpeciesLeft_         = right.numKRSpeciesLeft_;
    indexRNode_               = right.indexRNode_;
    rnodePos_final_           = right.rnodePos_final_;
    rnodePos_init_            = right.rnodePos_init_;
    rnodePos_init_init_       = right.rnodePos_init_init_;
    rRefPos_final_final_      = right.rRefPos_final_final_;
    rRefPos_final_            = right.rRefPos_final_;
    rRefPos_init_             = right.rRefPos_init_;
    rRefPos_init_init_        = right.rRefPos_init_init_;

    /*
     * Return the reference to the current object
     */
    return *this;
}
//======================================================================================================================
/*
 *
 * Destructor
 *
 * We need to manually free all of the arrays.
 */
Electrode_SurfaceRegion::~Electrode_SurfaceRegion()
{


}
//======================================================================================================================
//    Return the type of electrode
/*
 *  Returns the enum type of the electrode. This is used in the factory routine.
 *
 *  @return Returns an enum type, called   Electrode_Types_Enum
 */
Electrode_Types_Enum Electrode_SurfaceRegion::electrodeType() const
{
    return SIMPLE_DIFF_ET;
}
//======================================================================================================================
int
Electrode_SurfaceRegion::electrode_model_create(ELECTRODE_KEY_INPUT* ei)
{

    /*
     * Number of cells - hard code for now
     */


    /*
     * Initialize the arrays in this object now that we know the number of equations
     */
    init_sizes();




    return 0;
}

//=====================================================================================================================
int
Electrode_SurfaceRegion::setInitialConditions(ELECTRODE_KEY_INPUT* eibase)
{

    /*
     *  Downcast the Key input to make sure we are being fed the correct child object
     */
    ELECTRODE_RadialDiffRegions_KEY_INPUT* ei = dynamic_cast<ELECTRODE_RadialDiffRegions_KEY_INPUT*>(eibase);
    if (!ei) {
        throw CanteraError(" Electrode_SurfaceRegion::electrode_model_create()",
                           " Expecting a child ELECTRODE_KEY_INPUT object and didn't get it");
    }


    return 0;
}
//====================================================================================================================
void
Electrode_SurfaceRegion::init_sizes()
{


}

//====================================================================================================================
//    The internal state of the electrode must be kept for the initial and final times of an integration step.
/*
 *  This function advances the initial state to the final state that was calculated
 *  in the last integration step.
 *
 * @param Tinitial   This is the New initial time. This time is compared against the "old"
 *                   final time, to see if there is any problem.
 */
void  Electrode_SurfaceRegion::resetStartingCondition(double Tinitial, bool doResetAlways)
{
    /*
    * If the initial time is input, then the code doesn't advance
    */
    double tbase = MAX(t_init_init_, 1.0E-50);
    if (fabs(Tinitial - t_init_init_) < (1.0E-9 * tbase) && !doResetAlways) {
        return;
    }

}
//====================================================================================================================
//! Take the state (i.e., the final state) within the Electrode_Model and push it down
//! to the ThermoPhase objects and propogate it to all other aspects of the final state
/*!
 *  (virtual function from Electrode)
 *  This virtual function should be written so that child routines are not called by parent routines.
 *
 *  We take the values of spMoles_final_[] and propagate them down to the ThermoPhase
 *  objects in the electrode.
 *
 *  We also take the state of the electrode as described by the mole numbers and mole fractions
 *  and calculate the geometrical details involved with the model. This includes the radii and
 *  thicknesses of the regions and the existence of the annular regions with their associated boolean flags.
 *
 *  All of these properties are defined for the _final_ state.
 *
 *  Thus, this is the main routine that reconciles all of the state information within the object.
 *  At the end of this routine, all aspects of the final state are consistent with each other.
 *
 *  prerequisites: The object must have been already created.
 *
 *  Fundamental Variables:
 *       concKRSpecies_Cell_final_[]
 *
 */
void Electrode_SurfaceRegion::updateState()
{



}
//================================================================================================
/*
 * There is a small dependence on mf_external and mf_internal exhibited by this function
 */
void  Electrode_SurfaceRegion::extractInfo(std::vector<int>& justBornMultiSpecies)
{


}
//====================================================================================================================
void Electrode_SurfaceRegion::printElectrode(int pSrc, bool subTimeStep)
{

}
//===================================================================================================================

void Electrode_SurfaceRegion::printElectrodePhase(int iph, int pSrc, bool subTimeStep)
{


}


} // End of namespace Cantera
//======================================================================================================================
